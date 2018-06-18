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
  in precomputing the data that the FLT_HYBRID() routine
  requires. Details of FLT_HYBRID() can be found in
  flt_hybrid.c and the tech report

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


  All of the functions in this file are, of course, important.
  But the following are especially important because they deal
  with how the hybrid algorithm determines switch points (i.e.
  at a given order m, when does hybrid switch from the semi-naive
  portion of its algorithm to the simple-split portion), numsplits
  (how many splits the hybrid algorithm should perform in the
  simple-split portion of the algorithm), where the splits are,
  how much memory to allocate, etc etc. The functions are:

  1) SplitLocales() - determines where the switch point is and where
        the splits are, given order m and numsplits
  2) HybridPts() - how many splits to do at various orders; used
        in the forward spherical transform
  3) HybridPrecomp() - precomputes data used by the hybrid by calling
        the functions CosPmlTableGenLim() (in cospmls.c), SplitPml()
	(in this file), and ShiftSemiY() (in this file).

  The user is advised to read the documentation in this file describing
  how the Big Three functions work!
  

*********************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h> /* for memcpy */
#include "cospmls.h"
#include "newFCT.h"
#include "primitive.h"

#ifdef FFTPACK
#include "fftpack.h"
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define mymin(a, b) ((a) < (b) ? (a) : (b))


/*
  undefine if want to see some diagnostics
  from SplitLocales
*/
#define NOPRINT

/*********************************************************************

  SHIFTED LEGENDRE FUNCTION GENERATION CODE

  *********************************************************************/

/********************************************************************

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

  *******************************************************************/

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


/*********************************************************************
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

*******************************************************************/
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

/**********************************************************************

  now to generate the A's and the B's for the REVERSE recurrence
  
  *********************************************************************/

/*********************************************************************
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

  ******************************************************************/

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

/***********************************************************************

  In the relation

  P_{L-r} = A_{-r}^L P_L + B_{-r}^L P_{L+1}
  
  this routine generates B_{-r}^L . Initial conditions are
  
  B_0^L = 0, B_{-1}^L = 0 (?????)
  
  m = order of the problem
  r = how much to shift FROM L (NOTE: input r will be a positive number!)
  l = where to start shift (i.e. l = L in above *)
  n = number of points desired
  
  result is double array of size n 
  workspace is double array of size 4* n

  ***********************************************************************/


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

static void PmlPair( int m,
		     int l,
		     int bw,
		     double *pml,
		     double *pmlp1,
		     double *workspace)
{
  double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
  int i, n;

  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);



  /* Set the initial number of evaluation points to appropriate
     amount */
  n = 2 * bw;

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

  if ((m % 2) == 1) { /* need to divide out sin x */
    for (i=0; i<2*bw; i++)
      prev[i] /= sin(eval_args[i]);
  }


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

/***************************************************
  This function is exactly like the one above except
  in two aspects:
   
  1. It samples each of the associated Legendre functions
     at bw-many points, NOT (2*bw)-many

  2. It saves pml and pmlp1 in an "interlaced" manor.
     That is, in memory, the array looks like
     [ pml[0],pmlp1[0],pml[1],pmlp1[1],...,pml[bw-1],pmlp1[bw-1] ]
     This is done in an attempt to make more efficient use
     of the cache when I'm doing the Legendre transform.

  As a result, this function has one argument-less than
  PmlPair. Instead of *pml and *pmlp1 each pointing to
  double arrays of length bw, in this function *pml
  points to one array of length 2*bw.

  In its current form, FLT_HYBRID EXPECTS these associated
  Legendre functions to be saved in this interlaced manor.
  Change at your own peril !!!

   */

static void PmlPairX(int m,
		     int l,
		     int bw,
		     double *pml,
		     double *workspace)
{
  double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4;
  double *x_i, *eval_args;
  int i, n;

  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);


  /* Set the initial number of evaluation points to appropriate
     amount */
  n = bw;

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


  if ((m % 2) == 1) { /* need to divide out sin x */
    for (i=0; i<bw; i++)
      prev[i] /= sin(eval_args[i]);
  }
 

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
  for (i = 0 ; i < n ; i ++)
    {
      pml[2*i] = prevprev[i];
      pml[2*i+1] = prev[i];
    }

  /* that's it */

}


/************************************************************************/

/******

  THIS IS AN IMPORTANT FUNCTION !!! READ CAREFULLY !!!

  SplitLocales() determines where to put the split points,
  as well as things like numcoef, etc. It also determines
  the size of the arrays that will hold the precomputed
  Legendres for semi-naive and simple-split portions of
  the hybrid algorithm.

  In order to allocate the proper amount of space, this
  function must be executed first!!!

  A brief description of the hybrid algorithm is in
  order. Fix a bandwidth bw, an order m, a switch-point
  sw, and the number of splits numsplits. The hybrid
  algorithm computes a Legendre transform as follows.
  It applies the standard semi-naive algorithm to compute
  the coefficients f-hat[m] through f-hat[sw-1]. For
  the remaining coefficients, f-hat[sw] through f-hat[bw-1],
  the simple split algorithm is applied: numsplit-many
  split points are evenly distributed between sw and bw-1
  and the remaining coefficients are computed using the
  forward and reverse three-term recurrences.

  Example: bw = 256, m = 0, sw = 128, numsplits = 1
           
           Given this input, a split point is placed
	   at degree 192. The seminaive algorithm is
	   used to compute coefficients f-hat[0] through
	   f-hat[127]. Then the simple-split algorithm
	   computes coefficients f-hat[128] through f-hat[192]
	   using the reverse shifted Legendre polynomials
	   and f-hat[193] through f-hat[255] using the
	   forward shifted Legendre polynomials.

  In general, the switch point sw can be placed anywhere.
  If sw = m, then the hybrid algorithm becomes a pure
  simple-split algorithm. As it's written now, the function
  sets the switch point according to the following formula:

    case 16: case 32: *sw = m + (bw - m) / 2;
      break;
    case 64: case 128:
    case 256: case 512:
      if (m < (bw / 16))
	*sw = m + (bw - m) / 2;
      else
	*sw = m + 3 * (bw - m) / 4;
     break;
    case 1024:
      *sw = m + (bw - m) / 2;
      break;


  The goal you want is to minimize the total runtime of the
  algorithm. There are two parameters that can be adjusted:
  the location of the switch point and the number of splits.
  Our experiments have indicated that these settings are
  COMPUTER ARCHITECTURE-DEPENDENT !!! The optimum settings
  for a given bandwidth and order on a DEC Alpha may not
  be the optimum settings on an SGI Origin. For bw <= 512,
  the above formula seems to give the fastest timings when
  running on a DEC Alpha. When running on an HP Exemplar
  at bw = 1024, the other formula proved better. Other settings
  may give better times.

  The point is: your mileage may vary! To determine what's
  best for your platform, change the formula for the switch-point,
  try semi-naiving more (or fewer) coefficients before
  simple-splitting the rest. Change the number of split-points.

  Unfortunately, at the present time, there is no real guide
  in finding the optimal settings, except that it seems to be
  the case that doing fcts (which are an integral part of the
  simple-split portion of the hybrid algorithm) of lengths
  shorter that 64 (i.e. too many split-points) is probably not
  wise - the overhead begins to slow things down.
  

  The first three arguments are expected from the user:

  bw = bandwidth of problem
  m = order of transform   
  numsplits = number of splits desired


  The remaining arguments are set by SplitLocales(); they
  will be the inputs necessary to do a hybrid transform.

  sw = degree of where switching from semi-naive to
       simple-split

  quant_locs = how many splits there are (MAY NOT EQUAL numsplits
       because of rounding)

  loc = pointer to array of split locations, size numsplits+1;
      will be saved in the following format:

      [location of first split point, location of second split point,
       ...,location of last split point]
       

  cos_size = I have to precompute fcts of legendres for
             semi-naive ... how big should the array be?
	     It should be cos_size .

  simsplit_size = just like above but the fcts of the shifted
                  Legendres polynomials I'll need for the
		  simple-split portion of the algorithm.

  *****/

void SplitLocales( int bw,
		   int m,
		   int numsplits,
		   int *sw,
		   int *loc,
		   int *quant_locs,
		   int *cos_size,
		   int *simsplit_size )
{

  int curloc, adjust, i;
  int spacing, numcoef;

#ifndef NOPRINT    
  fprintf(stderr, "order = %d\t numsplits = %d\n\n", m, numsplits);
#endif

  /* setting switch point (where to switch from
     semi-naive to hybrid); set *sw = m for pure simple-split */ 

  switch ( bw )
    {
    case 16: case 32: *sw = m + (bw - m) / 2;
      break;
    case 64: case 128:
    case 256: case 512:
      if (m < (bw / 16))
	*sw = m + (bw - m) / 2;
      else
	*sw = m + 3 * (bw - m) / 4;
     break;
    case 1024:
      *sw = m + (bw - m) / 2;
      break;
    default:
      fprintf(stderr,"error in SplitLocales: can't do bw = %d\n",bw);
      exit ( 0 ) ;
    }

#ifndef NOPRINT   
  fprintf(stderr,"sw = %d\t", *sw);
#endif
  
  /* calculate the even spacing between split points,
     given the switch point and number of splits;
     and I want the spacing to be even */
  
  spacing = (bw - *sw) / (2 * numsplits);

#ifndef NOPRINT
    fprintf(stderr, "spacing = %d\t", spacing);
#endif

  while( ( spacing % 2) != 0 )
     (spacing)++;

#ifndef NOPRINT
    fprintf(stderr, "spacing = %d\n", spacing);
#endif

  /* set the number of coefficients I'll need
     when doing the sums in the frequency domain */

  numcoef = Power2Ceiling( spacing );
  if ( numcoef > bw )
    numcoef = bw;

  /*
    now calculate exactly where the split points
    are and store them
    */
  *quant_locs = 0;
  curloc = *sw + spacing;
  while (curloc < bw )
    {

#ifndef NOPRINT
      fprintf(stderr,"quant_locs = %d\t curloc = %d\n",
	      *quant_locs, curloc);
#endif

      loc[ *quant_locs ] = curloc;
      curloc += 2 * (spacing);
      (*quant_locs) ++;
    }
  
  /*** 
     quant_locs might not equal numsplits because
     of rounding
  
     In order to prevent aliasing of the very last
     coefficient, I have to make sure that the last potential
     coefficient is no further than numcoef away from bw;
     if it is I'll move sw up by just enough to compensate,
     so semi-naive might compute a little more than you think
     it might
     ***/

  adjust = bw - loc[ *quant_locs - 1 ];

#ifndef NOPRINT
    fprintf(stderr,"numcoef = %d\t adjust = %d\n", numcoef, adjust);
#endif

  if (adjust > numcoef)
    {
      *sw += (adjust - numcoef);

#ifndef NOPRINT
      fprintf(stderr,"new sw = %d\n", *sw);
#endif

      for(i = 0 ; i < *quant_locs; i++)
	{
	  loc[i] += (adjust - numcoef);
#ifndef NOPRINT
	  fprintf(stderr, "new loc[%d] = %d\n", i, loc[i]);
#endif
	}
    }

  *cos_size = ( (*sw + 1) + (*sw) * (*sw + 1) / 4 -
		(m - 1) * (m) / 4 );

  *simsplit_size = *quant_locs * ( (numcoef) * (numcoef) +
				 5 * (numcoef) + 4 );

}

/*******************************************************/
/*******************************************************/
/*******************************************************/

/****

  This is another VERY IMPORTANT function for the hybrid
  forward spherical transform !!!

  definitions (to repeat from above ... it's important!):

  lim = at what order the spherical transform algorithm
        goes from using the hybrid Legendre transform to
	the seminaive Legendre transform.

  switch point = at a given order m, the switch point is
        the degree of the first coefficient computed by
	the simple-split portion of the hybrid algorithm

  split points = at a given order m, the split points are
        the locations (in degrees) of where the splits
	used in the simple-split portion of the hybrid
	algorithm are made

  Ok.

  It is here where lim , the switch points and split points
  are all set. The hybrid spherical transform algorithm
  has a great many switches and knobs to set, in order to
  optimize performance on your platform. This function,
  along with SplitLocales(), is where the switches and
  knobs are set.

  In this function is defined the number of splits to
  perform in the hybrid algorithm for a range of orders
  (and different bandwidths). You may take these values
  as a starting point when optimizing them for your
  platform. The settings here are not guaranteed to
  minimize the hybrid spherical transform's runtime.

  Feel free to change the settings here and in SplitLocales()!

  YOUR MILEAGE MAY VARY. You need to determine in
  advance the number of splits at each order that will
  a) optimize runtime
  b) still be faster than the pure seminaive algorithm

  It may be in your case that for small bandwidths, you
  may want lim to be smaller than is set here, say bw/4.
  Setting lim higher than bw/2 will probably result in a
  slower algorithm. This is directly related to the number
  of split points at each order, and how they are placed.
  (See SplitLocales() above for details.) As the order
  increases, more split points are needed at the higher
  orders than the lower in order to maintain reasonable
  numerical stability. However, introducing more split-points
  introduces more overhead which makes the hybrid transform
  (at high orders) less competitive than the seminaive
  algorithm (which is stable for all orders). This is why
  the hybrid algorithm is not used at all orders when
  doing a forward spherical transform.

  So you must balance these two issues, lim and number
  of splits per order, when determining suitable parameters
  which would be optimal for your platform. In our tests,
  for bw = 64 through 256, a pure seminaive forward spherical
  transform was usually faster than the hybrid forward
  spherical transform; at bw = 512 they were roughly equal.
  At bw = 1024, the hybrid version was approximately 30% faster
  on an HP Exemplar.

  
  Use test_flt_hybrid to determine your optimal lim and
  numsplit settings.


  input: bw = bandwidth

  lim = pointer to ONE int: the final order that will use
        the hybrid Legendre transform algorithm in the forward
	spherical transform; for orders lim < m < bw it will
	use the seminaive algorithm; CAN BE CHANGED BY THE
	USER

  max_split = pointer to ONE int: for all the orders that
        use the hybrid Legendre transform, this is the maximum
	number of splits that are ever taken in a single
	order; CAN BE CHANGED BY THE USER

  The function assigns the following values:

  cos_size = int array containing the sizes of precomputed
        data (the fcts of the associated Legendre functions)
	used in the semi-naive portion of the hybrid
	algorithm (needed when writing to disk); of size
	bw

  shift_size =  int array containing the sizes of precomputed
        data (the fcts of the forward and reverse shifted
	Legendre polynomials) used in the simple-split
	portion of the hybrid algorithm (needed when writing
	to disk); of size bw

  max_cos_size, max_shift_size = each is a pointer to a different
        int: the max size of the precomputed data arrays that
	will be written to disk (max taken over all orders
	m = 0 through m = lim); needed when writing to disk

  This function also returns a pointer to an int array
  which contains all of the switch points, splitpoints
  and quantlocs;


  ******************************************************/

int *HybridPts( int bw, int *lim,
		int *cos_size, int *shift_size,
		int *max_split,
		int *max_cos_size,
		int *max_shift_size )
{
  int i;
  int space_needed;
  int *numsplits;
  int *loc; /* will be returned ! */
  int *loc_ptr;


  /*
    max_split: maximum number of splits at any order;
    CAN BE CHANGED BY THE USER
    */
    
  /*** where lim is set; CAN BE CHANGED BY THE USER ***/
  switch( bw )
    {
    case 16: case 32:
      *lim = bw / 4;
      break;
    case 64: case 128:
    case 256: case 512:
      *lim = bw / 4;
      break;
    case 1024:
      *lim = bw / 2;  /* set = bw / 4 if on puddleby  */
      break;
    default:
      fprintf(stderr,"error in HybridPts: can't do bw = %d\n",bw);
      exit ( 0 );
      break;
    }


  space_needed = 0;


  /***
    the numsplits array will contain the number of
    splits I do at each order I use the hybrid Legendre
    transform algorithm
    ***/

  numsplits = (int *) malloc(sizeof(int) * ((*lim) + 1));

  /***

    this is where I define how many splits per order <= lim
    
    numsplits[i] = the number of splits that will be done in the
          hybrid Legendre transform algorithm for order m = i.
          CAN BE CHANGED BY THE USER

    *****/

  switch( bw )
    {
    case 16:
      *max_split = 1;
      for ( i = 0 ; i <= (*lim) ; i ++)
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      break;
    case 32:
      *max_split = 1;
      for ( i = 0 ; i <= (*lim) ; i ++)
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      break;
    case 64:
      *max_split = 1;
      for ( i = 0 ; i <= (*lim) ; i ++)
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      break;
    case 128:
      *max_split = 1;
      for ( i = 0 ; i <= (*lim) ; i ++)
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      break;
    case 256:
      *max_split = 1;
      for ( i = 0 ; i < (*lim)/2 ; i ++ )
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      for ( i = (*lim)/2 ; i <= (*lim) ; i ++ )
	{
	  numsplits[i] = 1;
	  space_needed += numsplits[i];
	}
      break;
    case 512:
      *max_split = 3;
      for ( i = 0 ; i < 40 ; i ++ )
        {
          numsplits[i] = 1;
          space_needed += numsplits[i];
        }      
      for ( i = 40 ; i < 100 ; i ++ )
        {
          numsplits[i] = 2;
          space_needed += numsplits[i];
        }
      for ( i = 100 ; i <= *lim ; i ++)
        {
          numsplits[i] = 3;
          space_needed += numsplits[i];
        }
      break;
    case 1024:
      *max_split = 8 ;
      for (i = 0 ; i <= 69 ; i ++)
	{
	  numsplits[i] = 2;
	  space_needed += numsplits[i];
	}
      for (i = 70 ; i <= 115 ; i ++)
	{
	  numsplits[i] = 3;
	  space_needed += numsplits[i];
	}
      for (i = 116 ; i <= 187 ; i ++)
	{
	  numsplits[i] = 4;
	  space_needed += numsplits[i];
	}
      for (i = 188 ; i <= 265 ; i ++)
	{
	  numsplits[i] = 5;
	  space_needed += numsplits[i];
	}
      for (i = 266 ; i <= 454 ; i ++)
	{
	  numsplits[i] = 6;
	  space_needed += numsplits[i];
	}
      for (i = 455 ; i <= 512 ; i ++)
	{
	  numsplits[i] = 7;
	  space_needed += numsplits[i];
	}
    break;
    default:
      fprintf(stderr,"error in HybridPts: can't do bw = %d\n",bw);
      exit ( 0 );
      break;
    }


  /**
    loc will be used for sw's, quant_loc's and the
    loc's themselves:

    sw = switch points at the different orders
    quant_locs = how many splits done at the different
         orders
    locs = where the splits are
    ***/
  
  loc = (int *) malloc(sizeof(int) * (space_needed +
				      (2 * ((*lim) + 1))));


  /* ok, now I have to calculate the sizes of the
     arrays shiftstore, cos_pml_table, split_pml_table */

  loc_ptr = loc;
  for(i = 0 ; i <= (*lim) ; i++)
    {
      SplitLocales(bw, i, numsplits[i],
		   loc_ptr, loc_ptr+2, loc_ptr+1,
		   &cos_size[i], &shift_size[i]);
     
      loc_ptr += 2 + *(loc_ptr+1);

    }

  /*** find max cos_size, shift_size ***/
  *max_cos_size = 0;
  *max_shift_size = 0;
  for(i = 0 ; i <= (*lim) ; i ++)
    {
      *max_cos_size = mymax(*max_cos_size, cos_size[i]);
      *max_shift_size = mymax(*max_shift_size, shift_size[i]);
    }


  free(numsplits);


  return ( loc );

}


/**********************************************/

/**********************************************
  
  Hybrid_SetDptrs() is a utility function which assigns
  (double **) pointers to point at the different parts
  of the hybrid algorithm's precomputed data.

  bw = bandwidth

  lim = order where switching from hybrid algorithm to
      seminaive algorithm in the forward spherical transform

  loc = the int array produced by HybridPts(); contains
      switch points, number of splits, where the splits are

  precomp_space = double array which will contain the
      precomputed data for the hybrid algorithm;

  cos_size = int array containing the sizes of precomputed
        data (the fcts of the associated Legendre functions)
	used in the semi-naive portion of the hybrid
	algorithm (needed when writing to disk); of size
	bw

  shift_size =  int array containing the sizes of precomputed
        data (the fcts of the forward and reverse shifted
	Legendre polynomials) used in the simple-split
	portion of the hybrid algorithm (needed when writing
	to disk); of size bw

  cos_dptr = (double **) ptrs which will point to the precomputed
        data required by the seminaive portion of the hybrid
	algorithm; will point to arrays containing the fcts
	of the associated Legendre functions

  split_dptr = (double **) ptrs which will point to the precomputed
        data required by the simple-split portion of the hybrid
	algorithm; will point to arrays containing the
	associated Legendre functions

  shift_dptr = (double **) ptrs which will point to the precomputed
        data required by the simple-split portion of the hybrid
	algorithm; will point to arrays containing the fcts of
	the shifted Legendre polynomials

 ****************************************************/

   


void Hybrid_SetDptrs( int bw, int lim,
		      int *loc, double *precomp_space,
		      int *cos_size, int *shift_size,
		      double **cos_dptr,
		      double **split_dptr,
		      double **shift_dptr )
{
  int i , *loc2; 

  /***
    advance ptr to how many locs at first order
    ***/

  loc++;

  /*** now start assigning ***/
  cos_dptr[0] = precomp_space;
  precomp_space += cos_size[0];
  split_dptr[0] = precomp_space;
  precomp_space += 2 * bw * loc[0];
  shift_dptr[0] = precomp_space;
  precomp_space += shift_size[0];

  loc2 = loc; /* to point at the old quant_locs */
  loc += *loc + 2; /* to point at the next quant_locs */

  for (i = 1 ; i <= lim ; i ++)
    {

      cos_dptr[i] = precomp_space;
      precomp_space += cos_size[i];
      split_dptr[i] = precomp_space;
      precomp_space += 2 * bw * loc[0];
      shift_dptr[i] = precomp_space;
      if(i != lim)
	precomp_space += shift_size[i];

      /** point to the next quant_locs **/
      if (i != lim)
	{
	  loc += 2 + loc[0];
	  loc2 += 2 + loc2[0];
	}
    }
}

/*********************************************/

/***
  SplitPml() precomputes the Pmls at the split points;
  needed for the simple-split algorithm

  bw, m = bandwidth, order of problem
  quant_locs = how many splits?
  loc = an array of length quant_locs which lists
        where the splits are; has length = quant_locs;
	values in this array are set by SplitLocales()
  split_pml_table = where the pmls will be stored,
        an array of size 2 * bw * quant_locs

  workspace = of size 16 * bw

  The pmls will be stored in the following order:

  {p_{loc[0]-1} and p_{loc[0]} interlaced}, followed by 
  {p_{loc[1]-1} and p_{loc[1]} interlaced}, followed by
  ...
  {p_{loc[quant_locs]-1} and p_{loc[quant_locs]} interlaced}

  ****/

void SplitPml( int bw,
	       int m, 
	       int quant_locs,
	       int *loc,
	       double *split_pml_table,
	       double *workspace )
{
  int i;
  double *pml_ptr;

  pml_ptr = split_pml_table;

  for(i = 0 ; i < quant_locs ; i++)
    {
      PmlPairX(m, loc[i]-1, bw,
	       pml_ptr,
	       workspace);

      pml_ptr += (2 * bw);

    }

}

/*********************************************/

/****
  ShiftSemiY precomputes the fcts of shifted legendre
  polynomials needed for semi-naiving the segments of
  the simple-split portion of the hybrid algorithm.

  NOTE that since the fcts of the shifted legendre fcts are
  zero-striped, I'll be saving every-other non-zero coefficient.
  Memory doesn't grow on trees, you know.

  How this code works is not at all clear. (Hmmm, that's
  encouraging.) A lot of work has gone into making the
  hybrid algorithm as efficient as possible. One of the
  things done (as seen already in PmlPairX) is the "interlacing"
  of precomputed data. This was done with the intent of
  trying to use the cache as efficiently as possible. That
  means, I want to limit (if I can) the amount of thrashing
  that's done when I'm computing coefficients. By interlacing
  the precomputed data, I hope to lessen the thrashing.

  The code is precomputes as follow: at a given split point,
  let me first precompute (and take the fct of) all the
  reverse shifted Legendre transforms I'll need to calculate
  coefficients to the "left" of the split point. Not only will
  I interlace the precomputed data, but I'll also save it
  with the following assumption: the hybrid algorithm will
  first calculate the coefficient furthest to the left of
  the split point, then the furthest-1 coefficient, then
  the furthest-2, ..., the coefficient immediately to the left
  of the split point, and then at the split point. (So if
  my split point is at 31 and the furthest coefficient I
  expect to compute from that split-point is 0, the hybrid
  algorithm first computes f-hat[0], then f-hat[1], ...,
  f-hat[30], f-hat[31].) After precomputing all the reverse
  shifted Legendre polynomials needed at this split point,
  I precompute those I need go forward, to the right, interlacing
  them and saving them with the assumption that the hybrid
  algorithm will compute the other coefficients in order.
  (So in the example, after computing f-hat[31] using
  reverse shifted Legendre polynomials, the code will use
  forward shifted Legendres to compute f-hat[32], then
  f-hat[33], etc etc). So things are saved with the expectation
  that I always compute my coefficients from left to right,
  from lowest degree to highest.

  After I do the above, I go to the next split-point and repeat.


  So that is in words what the mysterious looking code below
  does.
  


  Inputs are:

  bw = bandwidth

  m = order of transform

  sw = where switching from semi-naive to simple-spit;
       it's the degree of the first coefficient I compute
       with simple-split

  loc = int array of size quant_locs, has locations
        of the splits

  quant_locs = number of split points

  bigstore = where I'll store things, has size
             simsplit_size (as calculated by SplitLocales)

  workspace = workspace area of size 4*length + 2*length + 2*length
              which equals 8 * length
             
  ***/

void ShiftSemiY(int bw,
		int m,
		int sw,
		int *loc,
		int quant_locs,
		double *bigstore,
		double *workspace)
{
  
  int ctr, endpt;
  int i, j;
  double *big_ptr, *storeshift, *scratch;
  double *AO, *BO, *AE, *BE;   /* O = odd degree, E = even */
  int forlength, revlength;
  int spacing, numcoef;
  double *CoswSave;

  /* compute spacing between split points and
     the number of coefficients I'll need given
     this spacing */

  spacing = (loc[0] - sw);

  numcoef = 2;
  while (numcoef < spacing )
    numcoef *= 2;

  if (numcoef > bw)
    numcoef = bw;

  spacing *= 2;

  /* assign pointers */
  storeshift = workspace;
  AO = storeshift + numcoef;
  BO = AO + numcoef;
  AE = BO + numcoef;
  BE = AE + numcoef;
  scratch = BE + numcoef;
  big_ptr = bigstore;

#ifdef FFTPACK
  CoswSave = precomp_dct( numcoef );
#endif

  ctr = 0;
  endpt = sw;

  while(ctr < quant_locs)
    {
      /* set lengths of the for-loops: how many
	 coefficients to compute in this split ? */
      
      revlength = loc[ctr] - endpt;
      
      if (ctr != (quant_locs - 1) )
	forlength = (spacing/2);
      else
	forlength = bw - loc[ctr];
      
      /* advance endpt for next loop */
      endpt += spacing;
      
      
      /* generate shifted legendres used in REVERSE
	 recurrence, take their cosine transform and
	 save those coefficients */

      /* now save for the case where i = revlength - 1 */
      A_revshifted(m, revlength - 1, loc[ctr] - 1 - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
      kFCT(storeshift, AO, scratch, numcoef, numcoef, 1);
#else
      memcpy(AO, storeshift, sizeof(double) * numcoef);
      DCTf( AO, numcoef, numcoef, CoswSave );
#endif
      
      B_revshifted(m, revlength - 1, loc[ctr] - 1 - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
      kFCT(storeshift, BO, scratch, numcoef, numcoef, 1);
#else
      memcpy(BO, storeshift, sizeof(double) * numcoef);
      DCTf( BO, numcoef, numcoef, CoswSave );
#endif


      for (j = 0 ; j < (revlength - 1) / 2 + 1 ; j ++ )
	{
	  big_ptr[0] = AO[2 * j + 1];
	  big_ptr[1] = BO[2 * j];
	  big_ptr += 2;
	}
 
      for(i = revlength - 1 ; i > 1 ; i -= 2)
	{
	  /* generate the polys and their coefficients */

	  A_revshifted(m, i-2, loc[ctr] - 1 - m, numcoef,
		       storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, AO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(AO, storeshift, sizeof(double) * numcoef);
	  DCTf( AO, numcoef, numcoef, CoswSave );
#endif

	  B_revshifted(m, i-2, loc[ctr] - 1 - m, numcoef,
		       storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, BO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(BO, storeshift, sizeof(double) * numcoef);
	  DCTf( BO, numcoef, numcoef, CoswSave );
#endif

	  A_revshifted(m, i-1, loc[ctr] - 1 - m, numcoef,
		       storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, AE, scratch, numcoef, numcoef, 1);
#else
	  memcpy(AE, storeshift, sizeof(double) * numcoef);
	  DCTf( AE, numcoef, numcoef, CoswSave );
#endif

      
	  B_revshifted(m, i-1, loc[ctr] - 1 - m, numcoef,
		       storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, BE, scratch, numcoef, numcoef, 1);
#else
	  memcpy(BE, storeshift, sizeof(double) * numcoef);
	  DCTf( BE, numcoef, numcoef, CoswSave );
#endif

	  /* now save stuff */
	  for(j = 0 ; j < i / 2 + 1 ; j++)
	    {
	      big_ptr[0] = AO[2 * j + 1];
	      big_ptr[1] = BO[2 * j];
	      big_ptr[2] = AE[2 * j];
	      big_ptr[3] = BE[2 * j + 1];
	      big_ptr += 4;
	    }
	  big_ptr[0] = AE[i + 1];
	  big_ptr ++;

	}
     
      /* store the lone coefficient needed when computing
	 at a split point */

      A_revshifted(m, 0, loc[ctr] - 1 - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
      kFCT(storeshift, AE, scratch, numcoef, numcoef, 1);
#else
      memcpy(AE, storeshift, sizeof(double) * numcoef);
      DCTf( AE, numcoef, numcoef, CoswSave );
#endif

      big_ptr[0] = AE[0];

      /* advance pointer, get it ready to start
	 saving more coefs */
      big_ptr ++;

 


      /*********************************************/


      /* generate shifted legendres used in FORWARD
	 recurrence, take their cosine transform and
	 save those coefficients */

      /* store the lone coefficient needed when computing
	 at a split point */

      Ashifted(m, 0, loc[ctr] - m, numcoef,
	       storeshift, scratch);
#ifndef FFTPACK
      kFCT(storeshift, AE, scratch, numcoef, numcoef, 1);
#else
      memcpy(AE, storeshift, sizeof(double) * numcoef);
      DCTf( AE, numcoef, numcoef, CoswSave );
#endif


      big_ptr[0] = AE[0];

      /* advance pointer, get it ready to start
	 saving more coefs */
      big_ptr ++;

      for(i = 1; i < forlength - 1; i += 2)
	{
	  /* generate the polys and their coefficients */

	  Ashifted(m, i, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, AO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(AO, storeshift, sizeof(double) * numcoef);
	  DCTf( AO, numcoef, numcoef, CoswSave );
#endif


	  Bshifted(m, i, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, BO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(BO, storeshift, sizeof(double) * numcoef);
	  DCTf( BO, numcoef, numcoef, CoswSave );
#endif

	  Ashifted(m, i+1, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, AE, scratch, numcoef, numcoef, 1);
#else
	  memcpy(AE, storeshift, sizeof(double) * numcoef);
	  DCTf( AE, numcoef, numcoef, CoswSave );
#endif
      
	  Bshifted(m, i+1, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, BE, scratch, numcoef, numcoef, 1);
#else
	  memcpy(BE, storeshift, sizeof(double) * numcoef);
	  DCTf( BE, numcoef, numcoef, CoswSave );
#endif

	  /* now save stuff */
	  for(j = 0 ; j < i / 2 + 1 ; j++)
	    {
	      big_ptr[3] = AO[2 * j + 1];
	      big_ptr[2] = BO[2 * j];
	      big_ptr[1] = AE[2 * j];
	      big_ptr[0] = BE[2 * j + 1];
	      big_ptr += 4;
	    }
	  big_ptr[0] = AE[i + 1];
	  big_ptr ++;

	}

      /* now save for the case where forlength is even */
      if( forlength % 2 == 0)
	{
	  Ashifted(m, forlength - 1, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, AO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(AO, storeshift, sizeof(double) * numcoef);
	  DCTf( AO, numcoef, numcoef, CoswSave );
#endif

	  Bshifted(m, forlength - 1, loc[ctr] - m, numcoef,
		   storeshift, scratch);
#ifndef FFTPACK
	  kFCT(storeshift, BO, scratch, numcoef, numcoef, 1);
#else
	  memcpy(BO, storeshift, sizeof(double) * numcoef);
	  DCTf( BO, numcoef, numcoef, CoswSave );
#endif

	  for (j = 0 ; j < (forlength - 1) / 2 + 1 ; j ++ )
	    {
	      big_ptr[0] = AO[2 * j + 1];
	      big_ptr[1] = BO[2 * j];
	      big_ptr += 2;
	    }
	}
      ctr++;
    }

#ifdef FFTPACK
  free( CoswSave );
#endif


}



/************************************************************************/


/***************************************************************

  HybridPrecomp will actually produce the precomputed data
  necessary for the forward hybrid spherical transform

  bw = bandwidth

  lim = order where switching from hybrid algorithm to
      seminaive algorithm in the forward spherical transform

  loc = the int array produced by HybridPts(); contains
      switch points, number of splits, where the splits are

  cos_dptr = (double **) ptrs which will point to the precomputed
        data required by the seminaive portion of the hybrid
	algorithm; will point to arrays containing the fcts
	of the associated Legendre functions

  split_dptr = (double **) ptrs which will point to the precomputed
        data required by the simple-split portion of the hybrid
	algorithm; will point to arrays containing the
	associated Legendre functions

  shift_dptr = (double **) ptrs which will point to the precomputed
        data required by the simple-split portion of the hybrid
	algorithm; will point to arrays containing the fcts of
	the shifted Legendre polynomials

  workspace = double array of size 16 * bw


  ************************************************************/


void HybridPrecomp( int bw, int lim,
		    int *loc,
		    double **cos_dptr,
		    double **split_dptr,
		    double **shift_dptr,
		    double *workspace )
{
  int i;

  fprintf(stdout,"Precomputing for the hybrid ...\n");

  for (i = 0 ; i <= lim ; i ++)
    {
      /*** precompute for the seminaive portion of the hybrid algorithm ***/
      CosPmlTableGenLim(bw, i, loc[0] + 1, cos_dptr[i], workspace);
      
      /*** now the pmls at the split points ***/
      SplitPml(bw, i, loc[1], loc + 2,
	       split_dptr[i],
	       workspace);
      
      /*** now the shifted Legendres to be used in the simple-split
	portion of the algorithm ***/
      ShiftSemiY( bw, i, loc[0], loc+2,
		  loc[1],
		  shift_dptr[i],
		  workspace);
      
      loc += 2 + loc[1];
      
    }
  
}
