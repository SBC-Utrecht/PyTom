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
  performing the fast Legendre transform (and not for
  precomputing data), as described in the tech report

       FFTs for the 2-Sphere - Improvements and Variations

       by D.M. Healy, Jr., D. Rockmore and Sean S.B. Moore
       Department of Computer Science
       Dartmouth College
       Technical Report PCS-TR96-292

  The algorithm used in this file to compute the Legendre
  transform is the Bounded DH-Mid variation described in
  the tech report. This algorithm is unstable for large
  order m's !

  The notation and terms I use in this code are based on
  those found in the tech report.


  *********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* to declare memcpy */

#include "csecond.h"
#include "newFCT.h"
#include "precomp_flt_dhmid.h"
#include "primitive.h"
#include "weights.h"

#ifdef FFTPACK
#include "fftpack.h"
#endif

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define mymin(a, b) ((a) < (b) ? (a) : (b))
#define PI 3.141592653589793

/***************************************************/
/* look up table for the power of 2 which is >= i */

static int Power2Ceiling( int i )
{

  int pow_2;
  
  pow_2 = 1;

  if (i == 0)
    return 0;
  else
    {
      while ( pow_2 < i )
	pow_2 *= 2;
      return pow_2 ;
    }

}

/************************************************************************/
/***
  This is the script-L operation defined in the technical report on page 13.
  It is the "critically sampled lowpass operator". To quote from the technical
  report: "The effect of L_p^N is to first compute the cosine representation
  of a vector (s) of length N, then remove all frequency components beyond p
  from s (smoothing) and finally, to keep only those samples necessary to
  represent this smoothed version (subsampling)."

  input = input signal, length n
  output = smoothed version of input, length p
  workspace = of size 3 * n

  NOTE: assuming that n >= p

  ***/

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
#endif

/**********************************************************************/
/**********************************************************************/
/*** just like above EXCEPT that it won't do the inverse cosine
  transform if p = length ... this saves some cpus 

  lowpoly, highpoly = arrays of type struct lowhigh, each of
                      length n

  CoswSave, CoswSave2 = the two precomputed data arrays that are needed
         for the forward (CoswSave) and inverse (CoswSave2) dcts. Two
	 are needed because there are two different sized transforms
	 to do


  ***/
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
  int i;

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

/**********************************************************/
/*
  Ok, this is it. SplitOp is the function which, given
  epsilon, applies the splitting operator to a vector


  In matrix notation:

            [ Z1 ]    [CritSamLow   0]                  [ Z1 ]
  SplitOp . [    ]  = [              ] . ShiftDiagMat . [    ]
            [ Z0 ]    [0   CritSamLow]                  [ Z0 ]


	              [CritSamLow   0]   [ W1 ]   [ Z1OUT ]
		    = [              ] . [    ] = [       ]
		      [0   CritSamLow]   [ W0 ]   [ Z0OUT ]

  m = order of the problem

  length = length of the FINAL segments when I switch to
           semi-naive -> I don't want to call ExpIFCT if
	   if have to (I don't have to when length = length_out)

  low_in, high_in = lower and higher degree Legendre inputs
  low_out, high_out = lower and higher degree Legendre outputs

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
		     struct lowhigh *highpoly)
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


/**********************************************************/
/*

  SplitOpRev is the function which, given epsilon, applies
  the splitting operator to a vector. It is used when splitting
  in the REVERSE direction!


  In matrix notation:

            [ Z1 ]    [CritSamLow   0]                  [ Z1 ]
  SplitOp . [    ]  = [              ] . ShiftDiagMat . [    ]
            [ Z0 ]    [0   CritSamLow]                  [ Z0 ]


	              [CritSamLow   0]   [ W1 ]   [ Z1OUT ]
		    = [              ] . [    ] = [       ]
		      [0   CritSamLow]   [ W0 ]   [ Z0OUT ]

  m = order of the problem

  length = length of the FINAL segments when I switch to
           semi-naive -> I don't want to call ExpIFCT if
	   if have to (I don't have to when length = length_out)

  low_in, high_in = lower and higher degree Legendre inputs
  low_out, high_out = lower and higher degree Legendre outputs

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
static void SplitOpRev( int m,
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
			struct lowhigh *highpoly)
#else
static void SplitOpRev( int m,
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
      wdata0[i] = high_in[i] * ul[i] + low_in[i] * ur[i];
      wdata1[i] = high_in[i] * ll[i] + low_in[i] * lr[i];
    }

  /*** now apply the critical sample lowpass operator ***/

#ifndef FFTPACK
  /* first smooth and subsample the lower degree fct */
  CritSamLowX(length, wdata0, low_out, length_in, length_out, scratch,
	      lowpoly, highpoly);

  /* now smooth and subsample the higher degree fct */
  CritSamLowX(length, wdata1, high_out, length_in, length_out, scratch,
	      lowpoly, highpoly);
#else
  /* first smooth and subsample the lower degree fct */
  CritSamLowX(length, wdata0, low_out, length_in, length_out, scratch,
	      CoswSave, CoswSave2 );

  /* now smooth and subsample the higher degree fct */
  CritSamLowX(length, wdata1, high_out, length_in, length_out, scratch,
	      CoswSave, CoswSave2 );
#endif

  /* and that's it */

}


/************************************************

  This is the Bounded DH_Mid algorithm. After
  dividing and conquering so many levels down, it will
  then switch to semi-naive.

  data = input array of size 2 * bw
  result = output array of size bw - m
  m = order of the transform
  bw = bandwidth of the problem
  Z = array of 6 double ptrs, each pointing to an array
      of size 2 * bw

  numsplits = how many levels down the "divide-and-conquer"
              tree do you want to go? CANNOT BE GREATER THAN 5!

	      Example: if bw = 256 and length_of_split = 2,
	      the ShiftDiagAll produces the matrices needed
	      to perform 2 levels of divides, so the problem
	      is reduced from two problems of size 128 to
	      four problems of size 64. (Recall that in the
	      DHMID algorithm, I'm spltting in the middle:
	      I'm using forward AND reverse three-term
	      recurrences.)

  workspace = of size 40 * bw

  lowpoly, highpoly = arrays of type struct lowhigh, each of
                      length bw; used in the FCT

  loops = how many times to loop through, for timing or
          error purposes ?

  errflag = 0: interested in timing, = 1: interested in
          determining accuracy

 **********************************************/

void FLT_DHMID( double *data,
		double *result,
		int m,
		int bw,
		double **Z,
		int numsplits,
		double *workspace,
		struct lowhigh *lowpoly,
		struct lowhigh *highpoly,
		int loops,
		int errorflag)
{
  int i, j, k, length, which_loop, zindex;
  int endpt, revlength, forlength;
  double *z0, *z1, *scratch;
  double **matrix, *lowptr, *highptr;
  const double *weights;
  double *eval_pts, *pmm, *storeshift, *storeshift2;
  int fudge, fudge_a, fudge_b, fudge2;
  double halfbw, tmp0, tmp1, *pml, *pmlp1;
  double *blank, *tmpdata;
  double tstart, tstop;
  int save_this_much, pow_of_two;
  int level, ctr;
  double *result_ptr;
  double *lowshiftptr, *highshiftptr;
  double *bigmatrixstore;
  double *shiftstore;
  int *shiftstore_ctr, *shiftstore_ctr_ptr;
  int *split, *splat;
  int **Zsplit, **Zsplat;
  int n;
  double **CoswSavedptr, *CoswSave, *tmp_ptr;

#ifdef FFTPACK
  n = 2 * bw;

  CoswSavedptr = (double **) malloc(sizeof(double *) * (numsplits + 2));
  
  for (i = 0 ; i < numsplits + 2 ; i ++)
    {
      CoswSavedptr[i] = precomp_dct( n );
      
      n /= 2;

      if ( i == 0 )
	n /= 2;

    }
  n = 2 * bw;
#endif




  splat = (int *) malloc(sizeof(int) * 132);
  Zsplat = (int **) malloc(sizeof(int *) * 6);
  split = (int *) malloc(sizeof(int) * 63);
  Zsplit = (int **) malloc(sizeof(int *) * 6);

  /* assign pointers by hand */
  Zsplat[0] = splat;
  Zsplat[1] = Zsplat[0] + 3;
  Zsplat[2] = Zsplat[1] + 5;
  Zsplat[3] = Zsplat[2] + 9;
  Zsplat[4] = Zsplat[3] + 17;
  Zsplat[5] = Zsplat[4] + 33;

  Zsplit[0] = split;
  Zsplit[1] = Zsplit[0] + 1;
  Zsplit[2] = Zsplit[1] + 2;
  Zsplit[3] = Zsplit[2] + 4;
  Zsplit[4] = Zsplit[3] + 8;
  Zsplit[5] = Zsplit[4] + 16;

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
  for(i = 1 ; i < numsplits + 1; i++)
    {
      Zsplat[i][0] = Zsplat[0][0];
      Zsplat[i][2 * fudge] = Zsplat[0][2];
      fudge *= 2;
    }
  fudge = 2;
  for(i = 0 ; i < numsplits + 1; i++)
    {
      for(j = 0 ; j < fudge - 1; j += 2)
	{
	  Zsplat[i][j+1] = (Zsplat[i][j+2]-Zsplat[i][j])/2 +
	    Zsplat[i][j];
	  Zsplit[i][j/2] = Zsplat[i][j+1];
	  fudge2 = 2;
	  for(k = i+1 ; k < numsplits + 1; k++)
	    {
	      Zsplat[k][fudge2*(j+1)] = Zsplat[i][j+1];
	      fudge2 *= 2;
	    }
	}
      fudge *= 2;
    }


  halfbw = ((double) bw) * 0.5 ;
  matrix = (double **) malloc(sizeof(double *) * 4);

  /* silly way to get correct power of two */
  i = 0;
  pow_of_two = 1;
  while ( i < numsplits )
    {
      pow_of_two *= 2 ;
      i ++ ;
    }

  /* length = after numsplits-many divisions, what's
              the length of the segments I'm left with? */
  length = bw / pow_of_two;
  fprintf(stdout,"FLT_DHMID\n");
  fprintf(stderr,"bw = %d\t m = %d\t", bw, m);
  fprintf(stderr,"numsplits = %d\tnumcoef = %d\n",
	  numsplits, length);

  /* space for precomputing the shifted legendres I'll
     need for semi-naiving; space to save number of
     elements written */

  shiftstore = (double *) malloc(sizeof(double) *
			       (4 * pow_of_two *
				(length + (length*(m+1))/2 +
				 (length-1)*length/4)));
  
  shiftstore_ctr = (int *) malloc(sizeof(int) * 4 * bw);

  /* space for all the matrices */
  bigmatrixstore = (double *) malloc(sizeof(double) *
				     16 * numsplits * bw);


  /* assign ptrs */
  blank = workspace;
  tmpdata = blank + 2 * bw;
  z0 = tmpdata + 2 * bw;;
  z1 = z0 + bw;
  storeshift = z1 + bw;
  storeshift2 = storeshift + bw;
  pmm = storeshift2 + bw;
  pml = pmm + 2 * bw;
  pmlp1 = pml + 2 * bw;
  eval_pts = pmlp1 + 2*bw;
  scratch = eval_pts + 2*bw;
  result_ptr = result;

  /* Look up quadrature weights */
  weights = get_weights(bw);

  /* generate the pmls */
  PmlPair(m, Zsplit[0][0] - 1, bw, pml, pmlp1, scratch);

  /* precompute the shifted legendres I'll need
     for semi-naiving */
  ShiftSemi(m,
	    bw,
	    Zsplit[numsplits - 1],
	    pow_of_two,
	    length,
	    shiftstore,
	    shiftstore_ctr,
	    scratch);

  
  /*** point to the array of quantities of the above coefs ***/
  shiftstore_ctr_ptr = shiftstore_ctr;

  /* precompute shifted legendre diagonal matrices */
  ShiftDiagAll(m,
	       bw,
	       numsplits,
	       Zsplit,
	       matrix, bigmatrixstore,
	       workspace);

  /* turn on the stopwatch */
  tstart = csecond();

  for(which_loop = 0 ; which_loop < loops ; which_loop++) {
    
    shiftstore_ctr = shiftstore_ctr_ptr;
    
    matrix[0] = bigmatrixstore;
    matrix[1] = matrix[0] + bw/2;
    matrix[2] = matrix[1] + bw/2;
    matrix[3] = matrix[2] + bw/2;


    /* now weight the data */
    if(errorflag == 0)
      for (i=0; i<2*bw; i++)
	tmpdata[i] = data[i] *  weights[i];
    else
      for (i = 0 ; i < 2 * bw ; i ++)
	tmpdata[i] = *data++ * weights[i];


    /* now I have to figure out how much to smooth
       and subsample. If the initial split is in the
       middle and it's a power of 2, I just save that
       power of 2 many coefficients. If I'm not in the
       middle, I have to take the larger portion and
       then the next power of 2 greater than that many
       coefficients. It's a pain but I don't know how
       else to keep the flexibility of shifting split
       points. */
    
    fudge = mymin(Zsplit[0][0], bw - Zsplit[0][0]);
    save_this_much = Power2Ceiling(fudge);
    

    /* now multiply the weighted data by pml;
       smooth, subsample and save */
    for(i = 0; i < 2 * bw; i++)
      blank[i] = tmpdata[i] * pml[i];
#ifndef FFTPACK
    CritSamLowX(length, blank, Z[0], 2 * bw, save_this_much,
		scratch, lowpoly, highpoly);
#else
    CritSamLowX(length, blank, Z[0], 2 * bw, save_this_much,
		scratch, CoswSavedptr[0], CoswSavedptr[1]);
#endif


    /* now do the same for pmlp1 */
    /* multiply the weighted data by pml;
       smooth, subsample and save */
    for(i = 0; i < 2 * bw; i++)
      blank[i] = tmpdata[i] * pmlp1[i];

#ifndef FFTPACK
    CritSamLowX(length, blank, Z[0] + save_this_much,
		2 * bw, save_this_much,
		scratch, lowpoly, highpoly);
#else
    CritSamLowX(length, blank, Z[0] + save_this_much,
		2 * bw, save_this_much,
		scratch, CoswSavedptr[0], CoswSavedptr[1]);
#endif

    
    /* ok, now have to apply the forward and reverse
       shifted Legendre polynomials */
    
    pow_of_two = 1;
    
    for(i = 0 ; i < numsplits - 1; i++)
      {
	
	for(j = 0 ; j < pow_of_two; j++)
	  {
	    /* first apply the splitting operator to go
	       in the reverse direction */
#ifndef FFTPACK	    
	    SplitOpRev(m,
		       length,
		       Z[i] + (save_this_much * 2 * j),
		       Z[i] + (save_this_much * 2 * j) +
		       save_this_much,
		       Z[i+1] + (save_this_much * 2 * j),
		       Z[i+1] + (save_this_much * 2 * j) +
		       (save_this_much / 2),
		       Zsplit[i][j],
		       Zsplit[i+1][2*j],
		       save_this_much,
		       save_this_much / 2,
		       matrix,
		       scratch,
		       lowpoly, highpoly);
#else
	    SplitOpRev(m,
		       length,
		       Z[i] + (save_this_much * 2 * j),
		       Z[i] + (save_this_much * 2 * j) +
		       save_this_much,
		       Z[i+1] + (save_this_much * 2 * j),
		       Z[i+1] + (save_this_much * 2 * j) +
		       (save_this_much / 2),
		       Zsplit[i][j],
		       Zsplit[i+1][2*j],
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

	    /* now apply the splitting operator to go
	       in the forward direction */
#ifndef FFTPACK	    
	    SplitOp(m,
		    length,
		    Z[i] + (save_this_much * 2 * j),
		    Z[i] + (save_this_much * 2 * j) +
		    save_this_much,
		    Z[i+1] + (save_this_much * (2 * j + 1)),
		    Z[i+1] + (save_this_much * (2 * j + 1)) +
		    (save_this_much / 2),
		    Zsplit[i][j],
		    Zsplit[i+1][2*j + 1],
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
		    Zsplit[i][j],
		    Zsplit[i+1][2*j + 1],
		    save_this_much,
		    save_this_much / 2,
		    matrix,
		    scratch,
		    CoswSavedptr[i+1],
		    CoswSavedptr[i+2]);
#endif	    
	    
	    if(j == (pow_of_two-1))
	      {
		matrix[0] += 4 * save_this_much;
		matrix[1] = matrix[0] + save_this_much/2;
		matrix[2] = matrix[1] + save_this_much/2;
		matrix[3] = matrix[2] + save_this_much/2;
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
    
    /* Now the array Z[numsplits-1] contains the coefficients
       I need to semi-naive. But before semi-naiving, I need to
       massage them a little ...
       */
    for(i = 0 ; i < 1 * bw ; i += length)
      {
	Z[numsplits-1][i] *= 2.0;
	for(j=0;j<length;j++)
	  Z[numsplits-1][i+j] *=  halfbw ;
      }
 

    /*****
      Ok, to compute a coefficient, here's what I do
      (roughly):

      <f,P_{L+r}> = <A_r^L , f P_L > + < B_r^L , f P_{L-1}>

      lowptr points to the (f P_{L-1}) (corresponding to the
      lower degree Legendre polynomial) and highptr points
      to the (f P_L) (corresponding to the higher degree
      Legendre polynomial).
      *****/
    lowptr = Z[numsplits - 1];
    highptr = lowptr + save_this_much;

    /* point to the array of fct'd shifted legendres */
    highshiftptr = shiftstore;
    lowshiftptr = highshiftptr + *shiftstore_ctr_ptr;        

    /*** point to the array of quantities of the above coefs ***/
    shiftstore_ctr_ptr = shiftstore_ctr;

    
    /* note: pow_of_two = 2^(numsplits - 1) */
    
    endpt = m;
    zindex = 0;

    while(zindex <  pow_of_two )
      {
	/* what degree is the split ? */
	level = Zsplit[numsplits-1][zindex];

	/* revlength, forlength = how many shifts
	   I perform from the current level using the
	   reverse and forward shifted Legendre polynomials
	   */	
	revlength = level - endpt;
	if (zindex != (pow_of_two-1))
	  forlength = (Zsplit[numsplits-1][zindex+1] -
		       Zsplit[numsplits-1][zindex]) / 2;
	else
	  forlength = (bw - Zsplit[numsplits-1][zindex]);
	
	endpt = level + forlength;
	

	/*************
	  ptrs to the REVERSE shifted Legendre polynomials
	  that I've precomputed in order to seminaive the coefficients:
	  
	  We have
      
	  P_{L-r-1} = (revB_{-r-1}^L  * P_L) + (revA_{-r-1}^L * P_{L-1})
      
	  highshiftptr points to revA_{-r-1}^L
      	  lowshiftptr points to revB_{-r-1}^L
      
	  *****/


	/*** do semi-naive in the reverse direction ***/
	for(i = 0; i < revlength ; i++)
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
		    tmp0 += highshiftptr[ctr] * lowptr[j-fudge_b];
		    tmp1 += lowshiftptr[ctr] * highptr[j-fudge_a];
		    ctr++;
		  }
	      }
	    else
	      {
		if(fudge2 == 1)
		  tmp0 = highshiftptr[0] * lowptr[0];
		else
		  {
		    for(j = 1; j < fudge2 - 1; j += 2)
		      {
			tmp0 += highshiftptr[ctr] * lowptr[j-fudge_b];
			tmp1 += lowshiftptr[ctr] * highptr[j-fudge_a];
			ctr++;
			
		      }
		    if(fudge_a==0)
		      tmp0 += highshiftptr[fudge2/2] * lowptr[fudge2-1];
		    else
		      tmp0 += lowshiftptr[fudge2/2] * highptr[fudge2-1];
		  }
	      }
	    /* save result */
	    result_ptr[level-i-1-m] = tmp0+tmp1;

	    /*** now have to shift the ptrs CAREFULLY:
	    watch the pointer arithmetic !!! ***/

	    highshiftptr += shiftstore_ctr[0] + shiftstore_ctr[1];
	    lowshiftptr += shiftstore_ctr[1] + shiftstore_ctr[2];
	    shiftstore_ctr+=2;

	  }

	/*************
	  ptrs to the FORWARD shifted Legendre polynomials
	  that I've precomputed in order to seminaive the coefficients:
	  
	  We have
      
	  P_{L+r} = A_r^L P_L + B_r^L P_{L-1}
      
	  highshiftptr points to A_r^L
	  lowshiftptr points to B_r^L
      
	  *****/


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
	    /* save result */
	    result_ptr[level+i-m] = tmp0+tmp1;

	    /*** now have to shift the ptrs CAREFULLY:
	    watch the pointer arithmetic !!! ***/

	    highshiftptr += shiftstore_ctr[0] + shiftstore_ctr[1];
	    lowshiftptr += shiftstore_ctr[1] + shiftstore_ctr[2];
	    shiftstore_ctr+=2;

	  }

	lowptr += 2 * save_this_much;
	highptr += 2 * save_this_much;
	
	zindex++;
      }
    
    /* normalize the result */
    if (m != 0)
      for(i=0;i<bw-m;i++)
	result_ptr[i] *= 2.0;

    if(errorflag == 1)
      result_ptr += (bw - m);
	

    
  } /* end of which_loop loop */
  
  /* turn off stopwatch */
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
    
  free(bigmatrixstore);
  free(shiftstore_ctr_ptr);
  free(shiftstore);
  free(matrix);
  free(Zsplit);free(split);free(Zsplat);free(splat);

#ifdef FFTPACK

  for (i = 0 ; i < numsplits + 2 ; i ++)
    free(CoswSavedptr[i]);
  
  free( CoswSavedptr );

#endif
  
}
