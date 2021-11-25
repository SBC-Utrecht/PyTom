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


/**********************************************************/
/***********************************************************

  This file contains some (all?) of the functions necessary
  in precomputing the data that the FLT_classic() routine
  requires. Details of FLT_classic() can be found in
  flt_classic.c and the tech report

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

**********************************************************/

#include <math.h>
#include <string.h> /* to declare memcpy */
#include "newFCT.h"
#include "primitive.h"
#include "weights.h"

#ifndef PI
#define PI 3.141592653589793
#endif

#define mymin(a, b) ((a) < (b) ? (a) : (b))
#define mymax(a, b) ((a) > (b) ? (a) : (b))

/************************************************************************/
/* SHIFTED LEGENDRE FUNCTION GENERATION CODE */
/************************************************************************/

/************************************************************************/
/************************************************************************/
/* In the relation

   P_{L+r} = A_r^L P_L + B_r^L P_{L-1}

   this routine generates A_r^L . Initial conditions are

   A_0^L = 1, A_{-1}^L = 0.

   m = order of the problem
   r = how much to shift
   l = where to start shift (i.e. l = L in above *)
   n = number of points desired
   result = double array of size n 
   workspace = double array of size 4 * n

   */

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
/************************************************************************/
/* In the relation

   P_{L+r} = A_r^L P_L + B_r^L P_{L-1}

   this routine generates B_r^L . Initial conditions are

   B_0^L = 0, B_{-1}^L = 1.

   m = order of the problem
   r = how much to shift
   l = where to start shift (i.e. l = L in above *)
   n = number of points desired

   result = double array of size n 
   workspace =  double array of size 4 * n
   */


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

/**********************************************************************/
/**********************************************************************/
/*
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


  */

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

/******************************************

  ShiftDiagAll will precompute all the shifted legendre polynomial
  matrices necessary for the divide-part of FLT_classic, i.e. this
  matrices will be used in SplitOp().

  m = order of shifted legendre polynomials
  bw = bandwidth of problem
  length_of_split = how many times to divide in the
                    divide-and-conquer?

		    Example: if bw = 256 and length_of_split = 2,
		    then ShiftDiagAll produces the matrices needed
		    to perform 2 levels of divides, so the one
		    problem of size 256 is reduced to 4 problems
		    each of size 64 (one 256 -> two 128 -> four 64)
		    
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

  length = bw;

  matrixspace = workspace;
  scratch = matrixspace + 4 * bw;

  bigptr = bigmatrixstore;

  for(i = 0 ; i < length_of_split; i++)
    {
      /*
	at split level i, need to generate 2^(i+1)
	many shifted matrices
	*/

      pow_of_two = 2;
      pow_of_two <<= i; /* so pow_of_two = 2^(i+1) */

      for(j = 0 ; j < pow_of_two; j += 2 )
	{
	  /*** calculate the matrix for the
	    forward recurrence, epsilon = 0 ***/

	  ShiftDiagMat(m, matrix,
		       split[i][j], split[i][j],
		       length,
		       matrixspace, scratch);

	  /*** now save in the bigspace ***/
	  memcpy(bigptr, matrixspace, sizeof(double) * 4 * length);

	  /*** now advance ptr for the next matrix ***/
	  bigptr += 4 * length;


	  /*** calculate the matrix for the
	    forward recurrence, epsilon = 1 ***/

	  ShiftDiagMat(m, matrix,
		       split[i][j], split[i][j + 1],
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

/****************************************************************/

/****
  ShiftSemi precomputes the fcts of shifted legendre polynomials needed
  for the semi-naive portion of FLT_classic. That is, after
  dividing so many levels, I need to complete the algorithm and
  compute the coefficients by applying the seminaive algorithm
  to remaining subproblems.

  NOTE that since the fcts of the shifted legendre fcts are
  zero-striped, I'll be saving every-other non-zero coefficient.
  Memory doesn't grow on trees, you know.

  Inputs are:

  m = order of transform

  splits = split locations at FINAL division (i.e. the last
           row of the split[][] array as defined in FLT_classic)

  lsplits = number of elements in splits = 2^numsplits (i.e. so
            there are 2^numsplits many subproblems left after
	    dividing)

  length = how many samples of the shifted legendre fcts should
           I take? This number will also be the input and output
	   lengths of the kFCT routine.
	   This number is "save_this_much".

  workspace = workspace area of size 8*length

  bigstore = where I'll store things, has size (are you ready?)

             4*lsplits*(length + (length*(m+1))/2 +
                        (length-1)*length/4)

  bigstore_ctr = array of ints of size 2 * splits * length

  I'm hoping this is a little more space than I need. We'll
  see.

  Since I'm having trouble counting things, after I write
  the coefficients of a given shifted legendre polynomial,
  in a separate array I'll write the value of

        ctr = the number of coefficients I've just
              saved in bigstore;

  that way I can shift my pointers correctly (as opposed to
  incorrectly) in FLT_classic

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

  int zindex, fudge, fudge_a, fudge_b, level;
  int n, i, j, ctr;
  double *bigstoreptr, *storeshift, *storeshift2, *scratch;
  int forlength;
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

  /***
    This is where I do the precomputing. Each time I loop
    through the while-construction, I precompute the shifted
    legendre polynomials needed for TWO subproblems.

    Example: If m = 0, bw = 256 and numsplits = 2, then I have
    FOUR subproblems, each of length 64 to seminaive. So I
    stay in this while-loop for TWO iterations before leaving
    ( for each iteration, zindex += 2)

    Sorry for the use of variables like "fudge". It makes
    the code clunky, but it was the only way I could keep
    track of what degree poly I was computing and which
    coefficients (even or odd-numbered) needed saving.
    Remember: the fcts of the shifted Legendre polynomials
    are zero-striped!

    ***/

  while(zindex < lsplits )
    {

      /***
	level = degree that I am shifting from
	forlength = how many shifts I perform
	***/

      level = splits[zindex];
      forlength = splits[zindex+1] - splits[zindex];

      /* generate shifted legendres used in
	 recurrence, take their cosine transform and
	 save those coefficients */

      /**
	Notation: starting at degree L, to shift forward
	by r degrees:

	P_{L+r} = A_r^L P_L + B_r^L P_{L-1}

	***/

      for(i = 0; i < forlength; i++)
	{

	  fudge = mymin(i + m + 1, length); /* I'll never need more
					       than length-many coeffs */

	  fudge_a = (i % 2);
	  fudge_b = 1 - fudge_a;

	  /* reset counter */
          ctr = 0;

	  /* calculate the A-shifts */
	  Ashifted(m, i, level - m, length,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );
#endif

	  /* since every other coefficient is 0,
	     don't bother saving it */

	  for(j = fudge_a; j < fudge; j+=2)
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

	  /* calculate the B-shifts */

	  Bshifted(m, i, level - m, length,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, storeshift2, scratch,
	       length, length, 1);
#else
	  memcpy(storeshift2, storeshift, sizeof(double) * length);
	  DCTf( storeshift2, length, length, CoswSave );;
#endif

	  /* since every other coefficient is 0,
	     don't bother saving it */

	  for(j = fudge_b; j < fudge; j+=2)
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


      /*** now for the next subproblem ***/

      level = splits[zindex+1];
      forlength = splits[zindex+2] - splits[zindex+1];


      for(i = 0; i < forlength; i++)
	{

	  fudge = mymin(i + m + 1, length);
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
	  DCTf( storeshift2, length, length, CoswSave );;
#endif

	  for(j = fudge_a; j < fudge; j+=2)
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
	  DCTf( storeshift2, length, length, CoswSave );;
#endif
	  for(j = fudge_b; j < fudge; j+=2)
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
      zindex+=2 ;
    }

#ifdef FFTPACK
  free(CoswSave);
#endif

}
