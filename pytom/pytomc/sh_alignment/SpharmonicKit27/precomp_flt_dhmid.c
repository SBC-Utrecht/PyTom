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


/**********************************************************
  This file contains some (all?) of the functions necessary
  in precomputing the data that the FLT_DHMID() routine
  requires. Details of FLT_DHMID() can be found in
  flt_dhmid.c and the tech report

       FFTs for the 2-Sphere - Improvements and Variations

       by D.M. Healy, Jr., D. Rockmore and Sean S.B. Moore
       Department of Computer Science
       Dartmouth College
       Technical Report PCS-TR96-292

  The notation and terms I use in this code are based on
  those found in the tech report.

  Since the code in this file is used only in computing
  the precomputed data, it's not optimized very much, if
  at all.

*********************************************************/


#include <math.h>
#include <string.h>  /* to declare memcpy */
#include "newFCT.h"
#include "primitive.h"
#include "weights.h"

#ifdef FFTPACK
#include "fftpack.h"
#endif



#ifndef PI
#define PI 3.141592653589793
#endif

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define mymin(a, b) ((a) < (b) ? (a) : (b))




/************************************************************************/
/* SHIFTED LEGENDRE FUNCTION GENERATION CODE */
/************************************************************************/

/***
  In the relation
  
  P_{L+r} = A_r^L P_L + B_r^L P_{L-1}
  
  this routine generates A_r^L . Initial conditions are
  
  A_0^L = 1, A_{-1}^L = 0.
  
  m = order of the problem
  r = how much to shift
  l = where to start shift (i.e. l = L in above *)
  n = number of points desired
  
  result is double array of size n 
  workspace is double array of size 4 * n
  ***/


static void Ashifted( int m,
		      int r,
		      int l,
		      int n,
		      double *result,
		      double *workspace )
{
  double *cospts, *pi, *piminus1, *pcurr;
  double di, dn, aval, cval;
  int i, j;

  /* allocate work space */
  cospts = workspace;
  pi = cospts + n;
  piminus1 = pi + n;
  pcurr = piminus1 + n;

  /* generate sample points */
  dn = (double) n;
  for (i=0; i<n; i++) {
    di = (double) i;
    cospts[i] = cos(
		    (((2.0 * di)+1.0) * PI)/(2.0 * dn));
  }
  
  if (r < 0) {
      for (i=0; i<n; i++)
	result[i] = 0.0;
  }

  /* return (zeros(n)); */

  /* piminus1 starts off as all ones */
  for (i=0; i<n; i++)
    piminus1[i] = 1.0;

  if (r == 0) {
      for (i=0; i<n; i++)
	result[i] = 1.0;
  }

  /* return (piminus1); */

  /* generate initial value of pi */
  for (i=0; i<n; i++)
    pi[i] = cospts[i] * L2_an(m,m+l);
  
  if (r == 1) {
      for (i=0; i<n; i++)
	result[i] = pi[i];
  }

  /* return(pi); */

  /* main loop */
  if (r > 1) {
      for (i=m+1; i<m+r; i++) {
	  aval = L2_an(m,i+l);
	  cval = L2_cn(m,i+l);
	  for (j=0; j<n; j++)
	    pcurr[j] = aval * cospts[j] * pi[j] +
	      cval * piminus1[j];
	  /* need to copy values */
	  memcpy(piminus1, pi, sizeof(double) * n);
	  memcpy(pi, pcurr, sizeof(double) * n);

      }
      memcpy(result, pcurr, sizeof(double) * n);
  }

/*  return (pcurr); */

}


/************************************************************************/
/***
  In the relation

  P_{L+r} = A_r^L P_L + B_r^L P_{L-1}
  
  this routine generates B_r^L . Initial conditions are
  
  B_0^L = 0, B_{-1}^L = 1.
  
  m = order of the problem
  r = how much to shift
  l = where to start shift (i.e. l = L in above *)
  n = number of points desired
  
  result is double array of size n 
  workspace is double array of size 4 * n
  ***/

static void Bshifted( int m,
		      int r,
		      int l,
		      int n,
		      double *result,
		      double *workspace )
{
  double *cospts, *pi, *piminus1, *pcurr;
  double di, dn, aval, cval;
  int i, j;

  /* allocate work space */
  cospts = workspace;
  pi = cospts + n;
  piminus1 = pi + n;
  pcurr = piminus1 + n;

  /* generate sample points */
  dn = (double) n;
  for (i=0; i<n; i++) {
    di = (double) i;
    cospts[i] = cos(
		    (((2.0 * di)+1.0) * PI)/(2.0 * dn));
  }
  
  if (r < 0) {
      for (i=0; i<n; i++)
	result[i] = 1.0;
  }

  /* piminus1 starts off as all zeros */
  for (i=0; i<n; i++)
    piminus1[i] = 0.0;

  if (r == 0) {
      for (i=0; i<n; i++)
	result[i] = 0.0;
  }

  /* return (piminus1); */

  /* generate initial value of pi */
  for (i=0; i<n; i++)
    pi[i] = L2_cn(m,m+l);
  
  if (r == 1) {
      for (i=0; i<n; i++)
	result[i] = pi[i];
  }

  /* return(pi); */

  /* main loop */
  if (r > 1) {
      for (i=m+1; i<m+r; i++) {
	  aval = L2_an(m,i+l);
	  cval = L2_cn(m,i+l);
	  for (j=0; j<n; j++)
	    pcurr[j] = aval * cospts[j] * pi[j] +
	      cval * piminus1[j];
	  /* need to copy values */
	  memcpy(piminus1, pi, sizeof(double) * n);
	  memcpy(pi, pcurr, sizeof(double) * n);

      }
      memcpy(result, pcurr, sizeof(double) * n);
  }

/*  return (pcurr); */

}

/************************************************************************

  now to generate the A's and the B's for the REVERSE recurrence

 ************************************************************************/

/****

  In the relation

  P_{L-r} = A_{-r}^L P_L + B_{-r}^L P_{L+1}
  
  
  this routine generates A_{-r}^L . Initial conditions are
  
  A_0^L = 1, A_{-1}^L = 0.
  
  m = order of the problem
  r = how much to shift FROM L (NOTE: input r will be a positive number!)
  l = where to start shift (l = L above)
  n = number of points desired
  
  result is double array of size n 
  workspace is double array of size 4 * n
***/

static void A_revshifted( int m,
			  int r,
			  int l,
			  int n,
			  double *result,
			  double *workspace )
{
  double *cospts, *pi, *piplus1, *pcurr;
  double di, dn, aval, cval;
  int i, j;

  /* allocate work space */
  cospts = workspace;
  pi = cospts + n;
  piplus1 = pi + n;
  pcurr = piplus1 + n;

  /* generate sample points */
  dn = (double) n;
  for (i=0; i<n; i++) {
    di = (double) i;
    cospts[i] = cos((((2.0 * di)+1.0) * PI)/(2.0 * dn));
  }

  
  if (r < 0) {
      for (i=0; i<n; i++)
	result[i] = 0.0;
  }

  /* return (zeros(n)); */

  /* piplus1 starts off as all ones */
  for (i=0; i<n; i++)
    piplus1[i] = 1.0;

  if (r == 0) {
      for (i=0; i<n; i++)
	result[i] = 1.0;
  }

  /* return (piplus1); */

  /* generate initial value of pi */
  /* aval = -L2_an(m,m+l); */
  /* cval = 1.0 / L2_cn(m,m+l); */
  aval = L2_ancn(m,m+l); 
  for (i=0; i<n; i++)
    pi[i] = cospts[i] * aval;
  
  if (r == 1) {
      for (i=0; i<n; i++)
	result[i] = pi[i];
  }

  /* return(pi); */

  /* main loop */
  if (r > 1) {
    for(i = 1 ; i < r ; i ++){
      /* aval = -L2_an(m, m + l - i); */
      aval = L2_ancn(m, m + l - i);
      /* cval = 1.0/L2_cn(m, m + l - i); */
      cval = L2_cn_inv(m,m+l-i);
      for(j = 0 ; j < n ; j ++)
	pcurr[j] = aval * cospts[j] * pi[j] +
	  cval * piplus1[j];
      /* need to copy values */
      memcpy(piplus1, pi, sizeof(double) * n);
      memcpy(pi, pcurr, sizeof(double) * n);
    }
    /*** now copy the final result ***/
    memcpy(result, pcurr, sizeof(double) * n);
  }

/*  return (pcurr); */

}

/************************************************************************/
/************************************************************************
  In the relation
  
  P_{L-r} = A_{-r}^L P_L + B_{-r}^L P_{L+1}
  
  this routine generates B_{-r}^L . Initial conditions are
  
  B_0^L = 0, B_{-1}^L = 0 (?????)
  
  m = order of the problem
  r = how much to shift FROM L (NOTE: input r will be a positive number!)
  l = where to start shift (i.e. l = L in above *)
  n = number of points desired
  
  result is double array of size n 
  workspace is double array of size 4 * n

*****/

static void B_revshifted( int m,
			  int r,
			  int l,
			  int n,
			  double *result,
			  double *workspace )
{
  double *cospts, *pi, *piplus1, *pcurr;
  double di, dn, aval, cval;
  int i, j;

  /* allocate work space */
  cospts = workspace;
  pi = cospts + n;
  piplus1 = pi + n;
  pcurr = piplus1 + n;

  /* generate sample points */
  dn = (double) n;
  for (i=0; i<n; i++) {
    di = (double) i;
    cospts[i] = cos(
		    (((2.0 * di)+1.0) * PI)/(2.0 * dn));
  }
  
  if (r < 0) {
      for (i=0; i<n; i++)
	result[i] = 0.0;
  }

  /* piplus1 starts off as all zeros */
  for (i=0; i<n; i++)
    piplus1[i] = 0.0;

  if (r == 0) {
      for (i=0; i<n; i++)
	result[i] = 0.0;
  }

  /* return (piplus1); */

  /* generate initial value of pi */
  for (i=0; i<n; i++)
    {
      /* pi[i] = 1.0 / L2_cn(m,m+l); */
      pi[i] = L2_cn_inv(m, m+l);
    }
  
  if (r == 1) {
      for (i=0; i<n; i++)
	result[i] = pi[i];
  }

  /* return(pi); */

  /* main loop */
  if (r > 1) {
    for ( i = 1 ; i < r ; i ++) {
      /* aval = -L2_an(m,m+l-i); */
      aval = L2_ancn(m,m+l-i);
      /* cval = 1.0 / L2_cn(m,m+l-i); */
      cval = L2_cn_inv(m,m+l-i);
      for ( j = 0 ; j < n ; j ++)
	pcurr[j] = aval * cospts[j] * pi[j] +
	  cval * piplus1[j];
      /* need to copy values */
      memcpy(piplus1, pi, sizeof(double) * n);
      memcpy(pi, pcurr, sizeof(double) * n);
    }
    /*** now copy the final result ***/
    memcpy(result, pcurr, sizeof(double) * n);
  }

/*  return (pcurr); */

}


/************************************************************************/
/* This is a special-purpose function.

   Given order m, degree l, bandwidth bw, this function produces
   the associated Legendre functions

   P(m,l)
   P(m,l+1)

   sampled at (2 * bw) points.

   The results are stored in the two arrays (each of size 2 * bw)
   which are pointed to by

   pml
   pmlp1   (read "pml plus 1")

   NOTE: This function assumes that when it's called, l >= m.

   workspace points to a double array of size (16 * bw)

*/

void PmlPair(int m,
	     int l,
	     int bw,
	     double *pml,
	     double *pmlp1,
	     double *workspace)
{
  double *prev, *prevprev;
  double *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
  int i, n;

  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);

  n = 2 * bw;

  /* main loop */

  /* Set the initial number of evaluation points to appropriate
     amount */

  /* now get the evaluation nodes */
  EvalPts(n,x_i);
  ArcCosEvalPts(n,eval_args);
    
  /* set initial values of first two Pmls */
  for (i=0; i<n; i++) 
    prevprev[i] = 0.0;
  if (m == 0) 
    {
      for (i=0; i<n; i++) 
	prev[i] = 1.0;
    }
  else 
    Pmm_L2(m, eval_args, n, prev);

  /* now generate remaining pmls  */

  for (i=0; i<l-m+1; i++) 
    {
      vec_mul(L2_cn(m,m+i),prevprev,temp1,n);
      vec_pt_mul(prev, x_i, temp2, n);
      vec_mul(L2_an(m,m+i), temp2, temp3, n);
      vec_add(temp3, temp1, temp4, n); /* temp4 now contains P(m,m+i+1) */

      /* now update P(m,m+i) and P(m,m+i+1) */
      memcpy(prevprev, prev, sizeof(double) * n);
      memcpy(prev, temp4, sizeof(double) * n);
    }

  /* save pml , pmlp1 */
  memcpy(pml, prevprev, sizeof(double) * n);
  memcpy(pmlp1, prev, sizeof(double) * n);

  /* that's it */

}

/************************************************************************/

/**********************************************************************

  ShiftDiagMat is the routine which generates the matrix M^epsilon
  (as defined in the tech report) found in the matrix S^epsilon_{L,M}
  defined in part 2 of Lemma 5, page 19.

  M^epsilon is a 2x2 block matrix with MxM diagonal blocks. The
  diagonals in the blocks are the appropriately sampled values
  of the appropriate shifted Legendre polynomials.

  m = order of polynomials
  matrix = array of 4 double pointers, each pointing to
           block of the ShiftDiagMatrix
  here = degree of where you are shifting from
  there = degree of where you want to shift to
  k = how many points to evaluate at;
  matrixspace = of size 4 * k
  workspace = of size 4 * k


  The matrix this routine produces will take you from
  degrees (here, here-1) to degrees (there, there-1)

  The blocks will be arranged so that the higher degree
  associated Legendre function will be on top, i.e.

  [ul  ur] [P_here    ]   [P_there    ]
  [      ].[          ] = [           ]
  [ll  lr] [P_{here-1}]   [P_{there-1}]


  ****/

static void ShiftDiagMat( int m,
			  double **matrix,
			  int here,
			  int there,
			  int k,
			  double *matrixspace,
			  double *workspace )
{

  int shift_amount;

  matrix[0] = matrixspace;
  matrix[1] = matrixspace+k;
  matrix[2] = matrixspace+(2*k);
  matrix[3] = matrixspace+(3*k);

  shift_amount = there - here;

  /*** A_r^L -> ul of M^epsilon ***/
  Ashifted(m,shift_amount,here-m,k,
	   matrix[0],workspace);
  
  /*** B_r^L -> ur of M^epsilon ***/
  Bshifted(m,shift_amount,here-m,k,
	   matrix[1],workspace);

  /*** A_(r-1)^L -> ll of M^epsilon ***/
  Ashifted(m,shift_amount-1,here-m,k,
	   matrix[2],workspace);

  /*** B_(r-1)^L -> lr of M^epsilon ***/
  Bshifted(m,shift_amount-1,here-m,k,
	   matrix[3],workspace);

  /* that's it ... trivial but hope it works */

}


/**********************************************************************/
/*
  ShiftDiagMatRev is just like the above routine EXCEPT
  it produces the matrix needed to go go BACKWARDS from a given
  degree. That is, it applies the REVERSE-shifted Legendre
  recurrences, NOT THE FORWARD.

  M^epsilon is a 2x2 block matrix with MxM diagonal blocks. The
  diagonals in the blocks are the appropriately sampled values
  of the appropriate reverse shifted Legendre polynomials.

  m = order of polynomials
  matrix = array of 4 double pointers, each pointing to
           block of the ShiftDiagMatrix_reverse
  here = degree to start shift from
  there = degree that want to shift to
  k = how many points to evaluate at;

  matrixspace = of size 4 * k
  workspace = of size 4 * k

  The matrix this routine produces will take you from
  degrees (here, here-1) to degrees (there-1, there)

  NOTE that since I'm going in reverse, there-1 is
  larger than there

  The blocks will be arranged so that the lower degree
  associated Legendre function will be on top, i.e.

  [ul  ur] [P_here    ]   [P_{there-1}]
  [      ].[          ] = [           ]
  [ll  lr] [P_{here-1}]   [P_there    ]

  ****/

static void ShiftDiagMatRev( int m,
			     double **matrix,
			     int here,
			     int there,
			     int k,
			     double *matrixspace,
			     double *workspace )
{

  int shift_amount;

  matrix[0] = matrixspace;
  matrix[1] = matrixspace+k;
  matrix[2] = matrixspace+(2*k);
  matrix[3] = matrixspace+(3*k);

  shift_amount = here - there;

  /*** B_r^L -> ul of M^epsilon ***/
  B_revshifted(m,shift_amount,here-1-m,k,
	       matrix[0],workspace);

  /*** A_r^L -> ur of M^epsilon ***/
  A_revshifted(m,shift_amount,here-1-m,k,
	       matrix[1],workspace);

  /*** B_(r-1)^L -> ll of M^epsilon ***/
  B_revshifted(m,shift_amount-1,here-1-m,k,
	       matrix[2],workspace);

  /*** A_(r-1)^L -> lr of M^epsilon ***/
  A_revshifted(m,shift_amount-1,here-1-m,k,
	       matrix[3],workspace);

  /* that's it ... trivial but hope it works */

}


/******************************************

  ShiftDiagAll will precompute all the shifted legendre polynomial
  matrices necessary for the divide-part of FLT_DHMID, i.e. these
  matrices will be used in SplitOp and SplitOpRev

  m = order of shifted legendre polynomials
  bw = bandwidth of problem
  length_of_split = how many times to divide in the
                    divide-and-conquer?

		    Example: if bw = 256 and length_of_split = 2,
		    the ShiftDiagAll produces the matrices needed
		    to perform 2 levels of divides, so the problem
		    is reduced from two problems of size 128 to
		    four problems of size 64. (Recall that in the
		    DHMID algorithm, I'm spltting in the middle:
		    I'm using forward AND reverse three-term
		    recurrences.)
		    
  split = array of ptrs to ptrs of ints; these contain the
          "here" and "there" that ShiftDiagMat needs

  matrix = array of 4 double pointers, each will (in the
           function) be pointing to block of a ShiftDiagMatrix

  bigmatrixstore = the double array which will contain all
                   sampled shifted legendre polynomials, of
		   size 8 * numsplits * bw (which should be
		   more than enough)

  workspace = array of size 8 * bw

  ***************************/

void ShiftDiagAll(int m,
		  int bw,
		  int length_of_split,
		  int **split,
		  double **matrix,
		  double *bigmatrixstore,
		  double *workspace)
{
  int i, j, length, pow_of_two;
  double *bigptr, *matrixspace, *scratch;

  length = bw/2;

  matrixspace = workspace;
  scratch = matrixspace + 2 * bw;

  bigptr = bigmatrixstore;

  for(i = 0 ; i < length_of_split - 1; i++)
    {
      /* stupid way to get correct power of two */

      pow_of_two = 1;
      for(j = 0; j < i; j++)
	pow_of_two *= 2;

      for(j = 0 ; j < pow_of_two; j++)
	{
	  /*** first calculate the matrix for the
	    reverse recurrence ***/

	  ShiftDiagMatRev(m, matrix,
			  split[i][j], split[i+1][2*j],
			  length,
			  matrixspace, scratch);
	  
	  /*** now save in the bigspace ***/
	  memcpy(bigptr, matrixspace, sizeof(double) * 4 * length);

	  /*** now advance ptr for the next matrix ***/
	  bigptr += 4 * length;

	  /*** now calculate the matrix for the
	    forward recurrence ***/

	  ShiftDiagMat(m, matrix,
		       split[i][j], split[i+1][2*j + 1],
		       length,
		       matrixspace, scratch);
	  
	  /*** now save in the bigspace ***/
	  memcpy(bigptr, matrixspace, sizeof(double) * 4 * length);
	  
	  /*** now advance ptr for the next matrix ***/
	  bigptr += 4 * length;

	}
      length /= 2;
    }
}

/************************************************************************/
/****
  ShiftSemi precomputes the fcts of shifted legendre polynomials needed
  for the semi-naive portion of FLT_DHMID.

  NOTE that since the fcts of the shifted legendre fcts are
  zero-striped, I'll be saving every-other non-zero coefficient.
  Memory doesn't grow on trees, you know.

  Inputs are:

  m = order of transform

  splits = split locations at FINAL division (i.e. the last
           row of the split[][] array as defined in FLT_DHMID)

  lsplits = number of elements in splits

  length = how many samples of the shifted legendre fcts should
           I take? This number will also be the input and output
	   lengths of the kFCT routine.
	   This number is "save_this_much".

  workspace = workspace area of size 4*length + 2*length + 2*length
              which equals 8 * length

  bigstore = where I'll store things, has size (are you ready?)

             4*lsplits*(length + (length*(m+1))/2 +
                        (length-1)*length/4)

  bigstore_ctr = array of ints of size 2 * splits * length

  I'm hoping this is a little more space than I need. We'll
  see.

  I really do hope there's more space that I need. Since I'm
  having trouble counting things, after I write the coefficients
  of a given shifted legendre, in a separate array I'll write
  the value of ctr -> the number of coefficients I've just
  saved in bigstore; that way I can shift my pointers correctly
  (as opposed to incorrectly) in FLT_DHMID

  ***/

void ShiftSemi(int m,
	       int bw,
	       int *splits,
	       int lsplits,
	       int length,
	       double *bigstore,
	       int *bigstore_ctr,
	       double *workspace)
{

  int zindex, fudge2, fudge_a, fudge_b, level;
  int i, j, ctr;
  double *bigstoreptr, *storeshift, *storeshift2, *scratch;
  int revlength, forlength, endpt;
  double *CoswSave;


  /* assign pointers */
  storeshift = workspace;
  storeshift2 = storeshift + length;
  scratch = storeshift2 + length;
  bigstoreptr = bigstore;

#ifdef FFTPACK
  CoswSave = precomp_dct( length );
#endif


  zindex = 0;
  endpt = m;
  while(zindex < (lsplits / 2))
    {
      /* what degree is the split? */
      level = splits[zindex];

      /* revlength, forlength = how many shifts
	 I perform from the current level using the
	 reverse and forward shifted Legendre polynomials
	 */

      revlength = level - endpt;
      if (zindex != (lsplits / 2 -1))
	forlength = (splits[zindex+1] -
		     splits[zindex]) / 2;
      else
	forlength = (bw - splits[zindex]);
      
      endpt = level + forlength;
      
      /* generate shifted legendres used in REVERSE
	 recurrence, take their cosine transform and
	 save those coefficients */

      for(i = 0; i < revlength ; i++)
	{
	  fudge2 = mymin(i + m + 1, length);
	  fudge_a = (i % 2);
	  fudge_b = 1 - fudge_a;

	  /* reset counter */
          ctr = 0;

	  A_revshifted(m, i, level - 1 - m, length,
		       storeshift, scratch);

#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );
#endif

	  for(j = fudge_a; j < fudge2; j+=2)
	    {
	      bigstoreptr[j/2] = storeshift2[j];
	      ctr ++;
	    }

	  /* save number of coeffs just saved */
	  *bigstore_ctr = ctr;
	  bigstore_ctr++;

	  /* advance pointer to store the next shifted guys */
	  bigstoreptr += ctr;

	  /* reset counter */
	  ctr = 0;

	  B_revshifted(m, i, level - 1 - m, length,
		       storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );
#endif


	  for(j = fudge_b; j < fudge2; j+=2)
	    {
	      bigstoreptr[j/2] = storeshift2[j];
	      ctr ++;
	    }

	  /* save number of coeffs just saved */
	  *bigstore_ctr = ctr;
	  bigstore_ctr++;

	  /* advance pointer to store the next shifted guys */
	  bigstoreptr += ctr;

	}


      /* generate shifted legendres used in FORWARD
	 recurrence, take their cosine transform and
	 save those coefficients */

      for(i = 0; i < forlength; i++)
	{

	  fudge2 = mymin(i + m + 1, length);
	  fudge_a = (i % 2);
	  fudge_b = 1 - fudge_a;

	  /* reset counter */
          ctr = 0;

	  Ashifted(m, i, level - m, length,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );
#endif
	  for(j = fudge_a; j < fudge2; j+=2)
	    {
	      bigstoreptr[j/2] = storeshift2[j];
	      ctr ++;
	    }

	  /* save number of coeffs just saved */
	  *bigstore_ctr = ctr;
	  bigstore_ctr++;

	  /* advance pointer to store the next shifted guys */
	  bigstoreptr += ctr;

	  /* reset counter */
	  ctr = 0;

	  Bshifted(m, i, level - m, length,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );
#endif

	  for(j = fudge_b; j < fudge2; j+=2)
	    {
	      bigstoreptr[j/2] = storeshift2[j];
	      ctr ++;
	    }

	  /* save number of coeffs just saved */
	  *bigstore_ctr = ctr;
	  bigstore_ctr++;

	  /* advance pointer to store the next shifted guys */
	  bigstoreptr += ctr;
	}

      zindex++;
    }

#ifdef FFTPACK
  free(CoswSave);
#endif


}

