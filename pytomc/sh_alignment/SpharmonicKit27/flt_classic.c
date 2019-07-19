/***************************************************************************
  **************************************************************************
  
                Spherical Harmonic Transform Kit 2.7
  
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 1997-2003  Sean Moore, Dennis Healy,
                        Dan Rockmore, Peter Kostelec
  
  
   Copyright 2004  Peter Kostelec, Dan Rockmore


     SpharmonicKit is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     SpharmonicKit is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
  
   Commercial use is absolutely prohibited.
  
   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*********************************************************

  This file contains the code which is actually used in
  the fast Legendre transform (and not for precomputing
  data), as described in Section 3 of the tech report

       FFTs for the 2-Sphere - Improvements and Variations

       by D.M. Healy, Jr., D. Rockmore and Sean S.B. Moore
       Department of Computer Science
       Dartmouth College
       Technical Report PCS-TR96-292

  Since it is the algorithm developed in Section 3, I
  call it "classic". It also has this name because this
  algorithm was the starting point in the development of
  more computationally (as opposed to asymptotically,
  which is what the original algorithm is) efficient
  and stable algorithms. This algorithm is unstable for
  large order m's !

  The notation and terms I use in this code are based on
  those found in the tech report.


  *********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* to declare memcpy */

#include "csecond.h"
#include "newFCT.h"
#include "precomp_flt_classic.h"
#include "primitive.h"
#include "weights.h"

#ifdef FFTPACK
#include "fftpack.h"
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#define mymin(a, b) ((a) < (b) ? (a) : (b))
#define mymax(a, b) ((a) > (b) ? (a) : (b))


/*******************************************************************

  This is the script-L operation defined in the technical report on
  page 13. It is the "critically sampled lowpass operator". To quote
  from the technical report: "The effect of L_p^N is to first compute
  the cosine representation of a vector (s) of length N, then remove
  all frequency components beyond p from s (smoothing) and finally,
  to keep only those samples necessary to represent this smoothed
  version (subsampling)."

  input = input signal, length n
  output = smoothed version of input, length p
  workspace = of size 3 * n

  NOTE: assuming that n >= p

  ******************************************************************/
#ifndef FFTPACK

static void CritSamLow( double *input,
			double *output,
			int n,
			int p,
			double *workspace )
{
  double *scratch, *tmpout;
  int permflag;

  /* assign workspace */
  tmpout = workspace;    /* needs to be of size p */
  scratch = tmpout + p;

  /* need to permute data */
  permflag = 1;
  
  /* do forward transform, return only p many coefficients */
  kFCT(input, tmpout, scratch, n, p, permflag);

  /* now reconstruct smoothed version */
  ExpIFCT(tmpout, output, scratch, p, p, permflag);

  /* that's it */
}

#endif /* FFTPACK */

/******************************************************************

  Just like above EXCEPT that it won't do the inverse cosine
  transform if p = length ... this saves some cpus; will also
  use lowhigh structs to get a tinier bit faster FCT 

  lowpoly, highpoly = arrays of type struct lowhigh, each of
                      length n

  CoswSave, CoswSave2 = the two precomputed data arrays that are needed
         for the forward (CoswSave) and inverse (CoswSave2) dcts. Two
	 are needed because there are two different sized transforms
	 to do

 ***************************************************************/

#ifndef FFTPACK
static void CritSamLowX( int length,
			 double *input,
			 double *output,
			 int n,
			 int p,
			 double *workspace,
			 struct lowhigh *lowpoly,
			 struct lowhigh *highpoly )
#else
static void CritSamLowX( int length,
			 double *input,
			 double *output,
			 int n,
			 int p,
			 double *workspace,
			 double *CoswSave,
			 double *CoswSave2 )
#endif
{
  double *scratch, *tmpout;
  int permflag;

  /* assign workspace */
  tmpout = workspace;    /* needs to be of size p */
  scratch = tmpout + p;

  /* need to permute data */
  permflag = 1;
  
  /* do forward transform, return only p many coefficients */
#ifndef FFTPACK
  kFCTX(input, tmpout, scratch, n, p, 1,
	lowpoly, highpoly);
#else
  memcpy(tmpout, input, sizeof(double) * n);
  DCTf( tmpout, n, p, CoswSave );
#endif

  if(p != length)
    {
      /* now reconstruct smoothed version */
#ifndef FFTPACK
      ExpIFCT(tmpout, output, scratch, p, p, permflag);
#else
      memcpy(output, tmpout, sizeof(double) * p);
      DCTb( output, p, CoswSave2 );
#endif
    }
  else
    memcpy(output, tmpout, sizeof(double) * p);

  /* that's it */
}


/*******************************************************

  Ok, this is it. SplitOp is the function which, given
  epsilon, applies the splitting operator to a vector.
  This splitting operator (to quote the tech report),
  "appropriately weight(s) and then lowpass(es) the
  input to provide the data for the low and high order
  transforms." In other words, this function is the
  "divide" in "divide-and-conquer".

  In matrix notation:

            [ Z1 ]    [CritSamLow   0]                  [ Z1 ]
  SplitOp . [    ]  = [              ] . ShiftDiagMat . [    ]
            [ Z0 ]    [0   CritSamLow]                  [ Z0 ]


	              [CritSamLow   0]   [ W1 ]   [ Z1OUT ]
		    = [              ] . [    ] = [       ]
		      [0   CritSamLow]   [ W0 ]   [ Z0OUT ]

  Z1,Z0 = input vectors, each of length LENGTH_IN

  ShiftDiagMat = 2x2 block matrix with LENGTH_IN x LENGTH_IN
                 diagonal blocks defined by sampled values of
		 shifted Legendre polynomials

  CritSamLow = (LENGTH_IN/2) x LENGTH_IN lowpass operator

  Z1OUT,Z0OUT = weighted and smoothed output vectors, each of
                length LENGTH_OUT x LENGTH_OUT

  Note: LENGTH_OUT = LENGTH_IN / 2

  Now for the function arguments ...


  m = order of the problem

  length = length of the FINAL segments when I switch to
           semi-naive -> I don't want to call ExpIFCT if
	   if have to (I don't have to when length = length_out)

  low_in, high_in = lower and higher degree vector inputs
                    (Z0, Z1 respectively)
  low_out, high_out = lower and higher degree vector outputs
                      (Z0OUT, Z1OUT, respectively)

  here = degree shifting from (these two arguments are NOT      )
  there = degree shifting to  (used; they are here for the      )
                              (historic reason that they helped )
			      (me keep track of things when I   )
			      (was writing FLT_classic          )

  length_in = length of input vectors
  length_out = length of output vectors

  matrix = array of double ptrs, pointing to the shifted Legendres
  workspace = of size 10 * length_in

  lowpoly, highpoly = arrays of type struct lowhigh, each of
                      length length_in; used in the FCT

  CoswSave, CoswSave2 = the two precomputed data arrays that are needed
         for the forward (CoswSave) and inverse (CoswSave2) dcts. Two
	 are needed because there are two different sized transforms
	 to do


*/
#ifndef FFTPACK
static void SplitOp( int m,
		     int length,
		     double *low_in,
		     double *high_in,
		     double *low_out,
		     double *high_out,
		     int here,
		     int there,
		     int length_in,
		     int length_out,
		     double **matrix,
		     double *workspace,
		     struct lowhigh *lowpoly,
		     struct lowhigh *highpoly )
#else
static void SplitOp( int m,
		     int length,
		     double *low_in,
		     double *high_in,
		     double *low_out,
		     double *high_out,
		     int here,
		     int there,
		     int length_in,
		     int length_out,
		     double **matrix,
		     double *workspace,
		     double *CoswSave,
		     double *CoswSave2 )
#endif
{

  int i;
  double *wdata0, *wdata1, *matrixspace, *scratch;
  double *ul, *ur, *ll, *lr;
  /* assign ptrs */
  
  wdata0 = workspace ;
  wdata1 = workspace + length_in ;
  matrixspace = wdata1 + length_in ;
  scratch = matrixspace + (4 * length_in) ;

  /*
    here's the key: ul = upper left block
                    ur = upper right block
		    ll = lower left block
		    lr = lower right block
		    */

  ul = matrix[0];
  ur = matrix[1];
  ll = matrix[2];
  lr = matrix[3];

  /*** now multiply the inputs by the shifted Legendre matrix ***/
  for(i = 0 ; i < length_in ; i++)
    {
      wdata1[i] = high_in[i] * ul[i] + low_in[i] * ur[i];
      wdata0[i] = high_in[i] * ll[i] + low_in[i] * lr[i];
    }

  /*** now apply the critical sample lowpass operator ***/

#ifndef FFTPACK
  /* first to the top half of the weighted data */
  CritSamLowX(length, wdata1, high_out, length_in, length_out, scratch,
	      lowpoly, highpoly);

  /* now to the bottom half of the weighted data */
  CritSamLowX(length, wdata0, low_out, length_in, length_out, scratch,
	      lowpoly, highpoly);
#else
  /* first to the top half of the weighted data */
  CritSamLowX(length, wdata1, high_out, length_in, length_out, scratch,
	      CoswSave, CoswSave2);

  /* now to the bottom half of the weighted data */
  CritSamLowX(length, wdata0, low_out, length_in, length_out, scratch,
	      CoswSave, CoswSave2);
#endif

  /* and that's it */

}

/************************************************************/
/************************************************************/
/*

  Computes the Legendre transform via the recursive splitting.
  The recursive splitting is based on Lemma 5 in the tech report.

  data = input array of size 2 * bw
  result = output array of size bw - m
  m = order of the transform
  bw = bandwidth of the problem
  loops = how many times to loop through, for timing or
          error purposes ?
  Z = array of 6 double ptrs, each pointing to an array
      of size 2 * bw
  numsplits = how many levels down the "divide-and-conquer"
              tree do you want to go? CANNOT BE GREATER THAN 5!

  workspace = of size 20 * bw
  lowpoly, highpoly = arrays of type struct lowhigh, each of
                      length bw; used in the FCT
  errflag = 0: interested in timing, = 1: interested in
          determining accuracy

  */
void FLT_classic( double *data,
		  double *result,
		  int m,
		  int bw,
		  int loops,
		  double **Z,
		  int numsplits,
		  double *workspace,
		  struct lowhigh *lowpoly,
		  struct lowhigh *highpoly,
		  int errorflag)
{
  int i, j, length, which_loop;
  double *z0, *z1, *scratch;
  double **matrix;

  double *lowptr, *highptr;

  double tstart, tstop;
  
  int zindex;
  const double *weights;
  double *tmpweights, *tmpdata, *eval_pts, *pmm;
  double *storeshift, *storeshift2;
  int fudge;
  double halfbw; /* fudge factor for later */
  int pow_of_two;
  double tmp0, tmp1;
  double *shiftstore, *bigmatrixstore;
  int *shiftstore_ctr, *shiftstore_ctr_ptr;
  int *splat, **Zsplat, fudge2, k;
  int save_this_much;
  double *shiftstoreptr, *lowshiftptr, *highshiftptr;
  int level, forlength, ctr;
  int fudge_a, fudge_b;
  double *result_ptr;

  int n;
  double **CoswSavedptr, *CoswSave, *tmp_ptr;

#ifdef FFTPACK

  n = 2 * bw;
  CoswSavedptr = (double **) malloc(sizeof(double *) * (numsplits + 2));
  
  /* ok, now precompute data for the forward dcts */
  for (i = 0 ; i < numsplits + 2 ; i ++)
    {
      CoswSavedptr[i] = precomp_dct( n );
      n /= 2;
    }
  n = 2 * bw;
  
#endif

  /* silly way to get correct power of two */
  pow_of_two = 1;
  for(i = 0; i < numsplits; i++)
    pow_of_two *= 2;

  /* length = after numsplits-many divisions, what's
              the length of the segments I'm left with? */

  length = bw / pow_of_two;

  fprintf(stdout,"FLT_CLASSIC\n");
  fprintf(stderr,"bw = %d\t m = %d\t", bw, m);
  fprintf(stderr,"numsplits = %d\tnumcoef = %d\n",
	  numsplits, length);

  /***
    allocate space for shifted legendre polynomials
    (those used when calculating the coefficients)
    ***/

  shiftstore = (double *) malloc(sizeof(double) *
				 (4 * pow_of_two *
				  (length + (length*(m+1))/2 +
				   (length-1)*length/4)));
  
  shiftstore_ctr = (int *) malloc(sizeof(int) * 4 * bw);
  
  /*** alocate space for the shifted Legendre polynomials
    that are used in the matrices, i.e. those matrices
    used when I'm dividing the problem into smaller
    sizes ***/
  bigmatrixstore = (double *) malloc(sizeof(double) *
				     2 * numsplits * 4 * bw);


  tmpweights = (double *) malloc(sizeof(double) * 2 * bw );

  halfbw = ((double) bw) * 0.5 ;

  /* assign ptrs */
  z0 = workspace;
  z1 = z0 + bw;
  storeshift = z1 + bw;
  storeshift2 = storeshift + bw;
  pmm = storeshift2 + bw;
  eval_pts = pmm + 2*bw;
  tmpdata = eval_pts + 2 * bw;
  scratch = tmpdata + 2*bw;

  result_ptr = result;
  
  splat = (int *) malloc(sizeof(int) * 132);
  Zsplat = (int **) malloc(sizeof(int *) * 6);

  /* assign pointers by hand */
  Zsplat[0] = splat;
  Zsplat[1] = Zsplat[0] + 3;
  Zsplat[2] = Zsplat[1] + 5;
  Zsplat[3] = Zsplat[2] + 9;
  Zsplat[4] = Zsplat[3] + 17;
  Zsplat[5] = Zsplat[4] + 33;

  /*******
    the format of the split locations will be like a pyramid
    
    Example: Let m = 0, bw = 512, numsplits = 3.
    Then
    
    Zsplat[0] = {0, 256, 512}
    Zsplat[1] = {0, 128, 256, 384, 512}
    Zsplat[2] = {0, 64, 128, 192, 256, 320, 384, 448, 512}
    
    *******/

  /* calculate split locations */
  Zsplat[0][0] = m ; Zsplat[0][2] = bw;
  fudge = 2;

  for(i = 1 ; i < numsplits+1; i++)
    {
      Zsplat[i][0] = Zsplat[0][0];
      Zsplat[i][2 * fudge] = Zsplat[0][2];
      fudge *= 2;
    }
  fudge = 2;
  for(i = 0 ; i < numsplits+1; i++)
    {
      for(j = 0 ; j < fudge - 1; j += 2)
	{
	  Zsplat[i][j+1] = (Zsplat[i][j+2]-Zsplat[i][j])/2 +
	    Zsplat[i][j];
	  fudge2 = 2;
	  for(k = i+1 ; k < numsplits+1; k++)
	    {
	      Zsplat[k][fudge2*(j+1)] = Zsplat[i][j+1];
	      fudge2 *= 2;
	    }
	}
      fudge *= 2;
    }

  length = bw / fudge2;

  matrix = (double **) malloc(sizeof(double *) * 4);

  i = 0 ;
  pow_of_two = 1;
  while ( i < numsplits )
    {
      pow_of_two *= 2;
      i ++ ;
    }

  /* precompute shifted legendre diagonal matrices */
  ShiftDiagAll(m,
	       bw,
	       numsplits,
	       Zsplat,
	       matrix, bigmatrixstore,
	       scratch);

  /* precompute the shifted legendres I'll need
     for semi-naiving */
  ShiftSemi(m,
	    bw,
	    Zsplat[numsplits - 1],
	    pow_of_two,
	    length,
	    shiftstore,
	    shiftstore_ctr,
	    scratch);

  /* Look up quadrature weights and compute Pmm */
  weights = get_weights(bw);


  /* precompute pmm in case m > 0 */
  if (m > 0)
    {
      ArcCosEvalPts(2*bw, eval_pts);
      Pmm_L2(m, eval_pts, 2*bw, pmm);
      
      for(i = 0 ; i < 2 * bw; i++)
	tmpweights[i] = weights[i] * pmm[i];
    }
  else
    memcpy ( tmpweights , weights , sizeof(double) * 2 * bw );

  /*
    turn on timer
    */
  tstart = csecond();

  for(which_loop = 0 ; which_loop < loops ; which_loop++)
    {
      /* reset counters, length */
      zindex = 0;
      pow_of_two = 1;
      pow_of_two <<= numsplits;
      length = bw / pow_of_two;

      /*********************/
      /* initialize arrays */
      /*********************/
      
      memset(Z[0], 0, sizeof(double) * bw);

      matrix[0] = bigmatrixstore;
      matrix[1] = matrix[0] + pow_of_two;
      matrix[2] = matrix[1] + pow_of_two;
      matrix[3] = matrix[2] + pow_of_two;


      /* weight the data */
      if(errorflag == 1)
	{
	  for (i = 0 ; i < 2 * bw ; i ++ )
	    tmpdata[i] = *data++ * tmpweights[i];
	}
      else
	{
	  for (i = 0 ; i < 2 * bw ; i ++ )
	    tmpdata[i] = data[i] * tmpweights[i];
	}

#ifndef FFTPACK
      /* first, smooth and subsample the data */
      CritSamLowX( length, tmpdata, Z[0] + bw, 2 * bw, bw, scratch,
		   lowpoly, highpoly );
#else
      /* first, smooth and subsample the data */
      CritSamLowX( length, tmpdata, Z[0] + bw, 2 * bw, bw, scratch,
		   CoswSavedptr[0], CoswSavedptr[1] );
#endif

      /* since I won't always be working with segments
	 whose lengths are powers of 2, I have to find
	 the first power of 2 greater than the length
	 of the two initial segments (there are two
	 segments, depending on where I make the first
	 split). This is because FCTs work only in powers
	 of 2. */

      fudge = mymax(Zsplat[0][1], bw - Zsplat[0][1]);

      /***
	save_this_much = first power of 2 greater than
	the maximum number of coefficients that I need
	at the first split ... should equal bw (always?)
	***/

      save_this_much = 1;
      while (save_this_much <= fudge)
	save_this_much *= 2;

    /* ok, now have to apply the forward and reverse
       shifted Legendre polynomials */
    
    /* Good luck */

    pow_of_two = 1;
    
    for(i = 0 ; i < numsplits; i++)
      {
	
	for(j = 0 ; j < pow_of_two; j ++ )
	  {
	    /* apply the splitting operator to go
	       in the forward direction */

#ifndef	FFTPACK
	    SplitOp(m,
		    length,
		    Z[i] + (save_this_much * 2 * j),
		    Z[i] + (save_this_much * 2 * j) +
		    save_this_much,
		    Z[i+1] + (save_this_much * (2 * j)),
		    Z[i+1] + (save_this_much * (2 * j)) +
		    (save_this_much / 2),
		    Zsplat[i][2*j],
		    Zsplat[i][2*j],
		    save_this_much,
		    save_this_much / 2,
		    matrix,
		    scratch,
		    lowpoly, highpoly);
#else
	    SplitOp(m,
		    length,
		    Z[i] + (save_this_much * 2 * j),
		    Z[i] + (save_this_much * 2 * j) +
		    save_this_much,
		    Z[i+1] + (save_this_much * (2 * j)),
		    Z[i+1] + (save_this_much * (2 * j)) +
		    (save_this_much / 2),
		    Zsplat[i][2*j],
		    Zsplat[i][2*j],
		    save_this_much,
		    save_this_much / 2,
		    matrix,
		    scratch,
		    CoswSavedptr[i+1],
		    CoswSavedptr[i+2]);
#endif	    

	    matrix[0] += 4 * save_this_much;
	    matrix[1] = matrix[0] + save_this_much;
	    matrix[2] = matrix[1] + save_this_much;
	    matrix[3] = matrix[2] + save_this_much;

#ifndef FFTPACK	  
	    SplitOp(m,
		    length,
		    Z[i] + (save_this_much * 2 * j),
		    Z[i] + (save_this_much * 2 * j) +
		    save_this_much,
		    Z[i+1] + (save_this_much * (2 * j + 1)),
		    Z[i+1] + (save_this_much * (2 * j + 1)) +
		    (save_this_much / 2),
		    Zsplat[i][2*j],
		    Zsplat[i][2*j + 1],
		    save_this_much,
		    save_this_much / 2,
		    matrix,
		    scratch,
		    lowpoly, highpoly);
#else
	    SplitOp(m,
		    length,
		    Z[i] + (save_this_much * 2 * j),
		    Z[i] + (save_this_much * 2 * j) +
		    save_this_much,
		    Z[i+1] + (save_this_much * (2 * j + 1)),
		    Z[i+1] + (save_this_much * (2 * j + 1)) +
		    (save_this_much / 2),
		    Zsplat[i][2*j],
		    Zsplat[i][2*j + 1],
		    save_this_much,
		    save_this_much / 2,
		    matrix,
		    scratch,
		    CoswSavedptr[i+1],
		    CoswSavedptr[i+2]);
#endif

	    if(j == (pow_of_two-1))
	      {
		if(i != (numsplits - 1))
		  {
		    matrix[0] += 4 * save_this_much;
		    matrix[1] = matrix[0] + save_this_much/2;
		    matrix[2] = matrix[1] + save_this_much/2;
		    matrix[3] = matrix[2] + save_this_much/2;
		  }
	      }
	    else
	      {
		matrix[0] += 4 * save_this_much;
		matrix[1] = matrix[0] + save_this_much;
		matrix[2] = matrix[1] + save_this_much;
		matrix[3] = matrix[2] + save_this_much;
	      }
	  }
	
	save_this_much /= 2;
	pow_of_two *= 2;
	
      }
    
    /* Now the array Z[numsplits] contains the coefficients
       I need to semi-naive. But before semi-naiving, I need
       to massage (i.e. scale) them a little ...
       */
    for(i = 0 ; i < 2 * bw ; i += length)
      {
	Z[numsplits][i] *= 2.0;
	for(j=0;j<length;j++)
	  Z[numsplits][i+j] *=  halfbw ;
      }


    /*** ok, here goes ... now I have to semi-naive the
      rest of the way ... almost there ! ***/

    /*****
      Ok, to compute a coefficient, here's what I do
      (roughly):

      <f,P_{L+r}> = <A_r^L , f P_L > + < B_r^L , f P_{L-1}>

      lowptr points to the (f P_{L-1}) (corresponding to the
      lower degree Legendre polynomial) and highptr points
      to the (f P_L) (corresponding to the higher degree
      Legendre polynomial).
      *****/
    
    lowptr = Z[numsplits];
    highptr = lowptr + save_this_much;
        
    zindex = 0;

    /*** point to the array of fct'd shifted legendres,
    and to their quantities ***/

    shiftstoreptr = shiftstore;
    shiftstore_ctr_ptr = shiftstore_ctr;


    /*****
      set up ptrs to the shifted Legendre polynomials that
      I've precomputed in order to seminaive the coefficients:
      
      We have
      
      P_{L+r} = A_r^L P_L + B_r^L P_{L-1}
      
      highshiftptr points to A_r^L
      
      lowshiftptr points to B_r^L
      
      *****/

    highshiftptr = shiftstoreptr;
    lowshiftptr = highshiftptr + *shiftstore_ctr_ptr;

    while(zindex <  pow_of_two  )
      {
	/* which degree am I starting from ? */
	level = Zsplat[numsplits-1][zindex];

	/* and how far to the next degree ? */
	forlength = (Zsplat[numsplits-1][zindex+1] -
		     Zsplat[numsplits-1][zindex]);

	/*** do semi-naive in the forward direction ***/
	for(i = 0; i < forlength; i++)
	  {
	    tmp0 = tmp1 = 0.0;
	    ctr = 0;

	    fudge2 = mymin(i + m + 1, save_this_much);
	    fudge_a = (i % 2);
	    fudge_b = 1 - fudge_a;

	    if((fudge2 % 2) == 0)
	      {
		for(j = 1; j < fudge2; j += 2)
		  {
		    tmp0 += highshiftptr[ctr] * highptr[j-fudge_b];
		    tmp1 += lowshiftptr[ctr] * lowptr[j-fudge_a];
		    ctr++;
		  }
	      }
	    else
	      {
		if(fudge2 == 1)
		  tmp0 = highshiftptr[0] * highptr[0];
		else
		  {
		    for(j = 1; j < fudge2 - 1; j += 2)
		      {
			tmp0 += highshiftptr[ctr] * highptr[j-fudge_b];
			tmp1 += lowshiftptr[ctr] * lowptr[j-fudge_a];
			ctr++;
		      }
		    if(fudge_a==0)
		      tmp0 += highshiftptr[fudge2/2] * highptr[fudge2-1];
		    else
		      tmp0 += lowshiftptr[fudge2/2] * lowptr[fudge2-1];
		  }
	      }
	    /* save the result */
	    result_ptr[level+i-m] = tmp0+tmp1;

	    /*** now have to shift the ptrs CAREFULLY:
	    watch the pointer arithmetic !!! ***/

	    highshiftptr += shiftstore_ctr_ptr[0] + shiftstore_ctr_ptr[1];
	    lowshiftptr += shiftstore_ctr_ptr[1] + shiftstore_ctr_ptr[2];
	    shiftstore_ctr_ptr+=2;

	  }


	lowptr += 2 * save_this_much;
	highptr += 2 * save_this_much;

	level = Zsplat[numsplits-1][zindex+1];
	forlength = (Zsplat[numsplits-1][zindex+2] -
		     Zsplat[numsplits-1][zindex+1]);

	/*** do semi-naive in the forward direction ***/
	for(i = 0; i < forlength; i++)
	  {
	    tmp0 = tmp1 = 0.0;
	    ctr = 0;

	    fudge2 = mymin(i + m + 1, save_this_much);
	    fudge_a = (i % 2);
	    fudge_b = 1 - fudge_a;

	    if((fudge2 % 2) == 0)
	      {
		for(j = 1; j < fudge2; j += 2)
		  {
		    tmp0 += highshiftptr[ctr] * highptr[j-fudge_b];
		    tmp1 += lowshiftptr[ctr] * lowptr[j-fudge_a];
		    ctr++;
		  }
	      }
	    else
	      {
		if(fudge2 == 1)
		  tmp0 = highshiftptr[0] * highptr[0];
		else
		  {
		    for(j = 1; j < fudge2 - 1; j += 2)
		      {
			tmp0 += highshiftptr[ctr] * highptr[j-fudge_b];
			tmp1 += lowshiftptr[ctr] * lowptr[j-fudge_a];
			ctr++;
		      }
		    if(fudge_a==0)
		      tmp0 += highshiftptr[fudge2/2] * highptr[fudge2-1];
		    else
		      tmp0 += lowshiftptr[fudge2/2] * lowptr[fudge2-1];
		  }
	      }

	    /* save the result */
	    result_ptr[level+i-m] = tmp0+tmp1;

	    /*** now have to shift the ptrs CAREFULLY:
	    watch the pointer arithmetic !!! ***/

	    if(( i != (forlength - 1)) || ( zindex != (pow_of_two - 2)))
	      {

		highshiftptr += shiftstore_ctr_ptr[0] + shiftstore_ctr_ptr[1];
		lowshiftptr += shiftstore_ctr_ptr[1] + shiftstore_ctr_ptr[2];
		shiftstore_ctr_ptr+=2;
		
	      }
	  }

	lowptr += 2 * save_this_much;
	highptr += 2 * save_this_much;
	
	zindex += 2;
      }
    
    /* normalize the result */
    if (m != 0)
      for(i=0;i<bw-m;i++)
	result_ptr[i] *= 2.0;

    if(errorflag == 1)
      result_ptr += (bw - m);


  } /* end of which_loop loop */

  tstop = csecond();
  
#ifndef WALLCLOCK
  fprintf(stdout,
	  "loops = %d:\ntotal cpu time =\t%.4e secs\n",
	  loops, tstop-tstart);
  fprintf(stdout, "avg cpu time:\t\t%.4e secs\n",
	  (tstop-tstart)/((double) loops));
#else
  fprintf(stdout,
	  "loops = %d:\ntotal wall time =\t%.4e secs\n",
	  loops, tstop-tstart);
  fprintf(stdout, "avg wall time:\t\t%.4e secs\n",
	  (tstop-tstart)/((double) loops));
#endif


  free(Zsplat) ; free(splat);
  free(tmpweights);
  free(bigmatrixstore);
  free(shiftstore_ctr);
  free(shiftstore);
  free(matrix);

#ifdef FFTPACK

  for (i = 0 ; i < numsplits + 2 ; i ++)
    free(CoswSavedptr[i]);
  
  free( CoswSavedptr );

#endif
  
}


