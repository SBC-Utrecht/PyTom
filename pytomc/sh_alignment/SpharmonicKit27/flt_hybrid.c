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
  performing the fast hybrid Legendre transform (and NOT for
  precomputing data), both the "stand-alone" version
  and the spherical transform version. The algorithm is
  described in the tech report

       FFTs for the 2-Sphere - Improvements and Variations

       by D.M. Healy, Jr., D. Rockmore and Sean S.B. Moore
       Department of Computer Science
       Dartmouth College
       Technical Report PCS-TR96-292

  The notation and terms I use in this code are based on
  those found in the tech report.


  *********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memcpy */

#include "csecond.h"
#include "cospmls.h"
#include "newFCT.h"
#include "oddweights.h"
#include "seminaive.h"
#include "weights.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define mymin(a, b) ((a) < (b) ? (a) : (b))
#define PI 3.141592653589793


/************************************************************************/


/************************************************

  Fast Legendre Transform, hybrid version: apply semi-naive
  algorithm to compute first group of coefficients; compute
  remainder via simple-split.

  This code has been written to (some extent) limit the
  amount of thrashing that is done in the cache. Arrays
  the can be interlaced (like the precomputed data) probably
  are. The simple-split portion of code is a bit mysterious.
  There were a lot of iterations between the initial,
  understandable version and this final, non-inituitive
  one. Actually, it is intuitive, but not quite clear.
  Here's a description of how the hybrid code does the
  simple-split portion of the algorithm, how it computes
  the coefficients:

  The code is computes as follow: at a given split point,
  let me first compute (and take the fct of) all the
  coefficients to the "left" of the split point. Since
  the precomputed data is interlaced (see ShiftSemiY() in
  precomp_flt_hybrid.c), I have to be very careful with
  my pointers. First I calculate the coefficient furthest
  to the left of the split point, then the furthest-1 coefficient,
  then the furthest-2, ..., the coefficient immediately to the left
  of the split point, and then at the split point. (So if
  my split point is at 31 and the furthest coefficient I
  expect to compute from that split-point is 0, the hybrid
  algorithm first computes f-hat[0], then f-hat[1], ...,
  f-hat[30], f-hat[31].) After computing all the coefficients
  to the left of the split point, I compute all those to
  the right. (So I compute f-hat[32], then f-hat[33], ...)
  Again, I have to be careful with pointers since my
  precomputed data is interlaced.

  After I do the above, I go to the next split-point and repeat.
  Note that I'm always computing the coefficients in order,
  from lowest degree to highest. In one sense this seems
  reasonable, in another sense it does not. But this is
  the way it's done.


  m = order of transform
  bw = bandwidth
  data = input pointer to an array of length 2 * bw
  result = output pointer to an array of length bw - m
  Z = pointer to an array of pointers, of length quant_locs
  workspace = array of size 23 * bw (more than enough)
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT
  loc = pointer to where the split points are, set by SplitLocale()
  shiftstore = pointer to precomputed fcts of shifted Legendre
             polynomials used in the simple-split portion of
	     the hybrid algorithm
  cos_pml_table = pointer to precomputed fcts of the Legendre
             functions used in the seminaive portion of the
	     hybrid algorithm
  split_pml_table = pointer to the precomputed Legendre functions
             used in the hybrid algorithm; used to weight the
	     data at the split points
  loops = how many times to loop through, for timing or
          error purposes ?
  errflag = 0: interested in timing, = 1: interested in
          determining accuracy

  ***/

void FLT_HYBRID( int m,
		 int bw,
		 int sw,
		 double *data,
		 double *result,
		 double **Z,
		 double *workspace,
		 struct lowhigh *poly,
		 int quant_locs,
		 int *loc,
		 double *shiftstore,
		 double *cos_pml_table,
		 double *split_pml_table,
		 int loops,
		 int errorflag)
{
  int n, i, j;
  int which_loop;
  int spacing, numcoef, ctr, endpt;
  int revlength, forlength;
  int my_ctr;

  double halfbw, tmp0, tmp1, tmp2, tmp3;
  double tstart, tstop;
  double *result_ptr;
  double *Zptr;
  double *tmp_result;

  const double *weights;
  double *split_pml_ptr;
  double *lowptr;
  double *blank, *tmpdata, *datafct;
  double *scratch;
  double *shiftptr;

  double times_two ; /* to adjust in case the order m = 0 */

  struct lowhigh *lowpoly ;
  double *CoswSave, *CoswSave2;

  /* assign ptrs for struct */
  lowpoly = poly;

  /* assign ptrs for the double arrays */
  blank = workspace;
  tmpdata = blank + 2 * bw;
  datafct = tmpdata + 2 * bw;
  scratch = datafct + bw;
  result_ptr = result;

  halfbw = ((double) bw) * 0.5 ;
  n = 2 * bw;

#ifdef FFTPACK

  CoswSave = precomp_dct( n );
  CoswSave2 = precomp_dct( bw );

#endif

  /*** get the weights I'll need ***/
  if (m % 2 == 0)
    weights = get_weights(bw);
  else
    weights = get_oddweights(bw);

  /* get spacing of split points, number of coefficients
   to be computed per segment (during the simple-split stage) */
  spacing = loc[0] - sw;
  numcoef = 2;
  while ( numcoef < spacing )
    numcoef *= 2;


#ifndef NOPRINTX
  fprintf(stderr,"FLT_SIMSPLIT\n");
  fprintf(stderr,"bw = %d\t m = %d\t", bw, m);
  fprintf(stderr,"sw = %d\t numsplits = %d\tnumcoef = %d\n",
	  sw, quant_locs, numcoef);
#endif

  /* fudge factor to normalize coeffs - depends on
     order m */
  if(m == 0)
    times_two = (double) 1;
  else
    times_two = (double) 2;


  /* turn on the stopwatch */
  tstart = csecond();

  for(which_loop = 0 ; which_loop < loops ; which_loop++) {

#ifndef FFTPACK    
    /* weight the data */
    if(errorflag == 0)
      for (i=0; i<2*bw; i++)
	tmpdata[i] = data[i] * weights[i];
    else
      for(i = 0; i < 2 * bw ; i ++)
	tmpdata[i] = *data++ * weights[i];

    /* take FCTX of the weighted data */
    kFCTX(tmpdata, datafct, scratch, 2 * bw, bw, 1,
	  lowpoly, lowpoly + bw);
#else
    /* note that I can safely use datafct to store
       something of length 2*bw, at least for a
       little bit */

    /* weight the data */
    if(errorflag == 0)
      for (i=0; i<2*bw; i++)
	datafct[i] = data[i] * weights[i];
    else
      for(i = 0; i < 2 * bw ; i ++)
	datafct[i] = *data++ * weights[i];
    /* take FCTX of the weighted data */
    DCTf( datafct, n, bw, CoswSave );
#endif

    /* need to normalize */
    for (i=0; i<bw; i++)
      datafct[i] *= halfbw;    

    /* now to get a subsampled signal */
#ifndef FFTPACK
    ExpIFCT(datafct, tmpdata, scratch, bw, bw, 1);
#else
    memcpy(tmpdata, datafct, sizeof(double) * bw );
    DCTb( tmpdata, bw, CoswSave2 );
#endif

    /* the subsampled data in the tmpdata array will be
       used in the simple-split portion of the algorithm */

    /**** now I normalize the first coefficient in datafct,
    preparing it for the seminaive portion of the algorithm ****/
    datafct[0] *= 2.0;

    /*****************************************

      Now to compute the first group of coefficients
      via seminaive

      SemiNaiveReducedX assumes that the data is ALREADY
      cosine-transformed. Otherwise, it's just
      like SemiNaiveReduced.

      NOTE: sw-m is HOW MANY coefficients I want
      to compute !!!

      NOTE I can use the blank array now for "scratch"
      area

      ******************************************/

    SemiNaiveReducedX(datafct, bw, m, sw - m,
		      result_ptr,
		      cos_pml_table,
		      blank); 

    /*******************************************

      Having seminaived the first portion of the coefficients,
      now have to get around to computing the remaining
      coefficients by simple-split

      ******************************************/


    /* assign point to precomputed pmls */
    split_pml_ptr = split_pml_table;
  
    /* now to pointwise multiply the smoothed data
       with the pmls that I've save at the split points
       (there are two pmls per split point, a lower
       degree one and a higher ) */
    for ( i = 0 ; i < quant_locs ; i++)
      {
	Zptr = Z[i];
	for ( j = 0 ; j < bw ; j++)
	  {
	    Zptr[j] = tmpdata[j] * split_pml_ptr[0];
	    Zptr[bw + j] = tmpdata[j] *
	      split_pml_ptr[1];
	    split_pml_ptr += 2;
	  }
      }


    /* now that I've point-wise multiplied, I have
       take the cosine transforms of this data.
       I'll just keep using the same space. */

    for (i = 0 ; i < quant_locs ; i++)
      {
#ifndef FFTPACK
	/*** the low degree poly  ***/
	kFCTX(Z[i], blank, scratch, bw, numcoef, 1,
	      lowpoly, lowpoly + bw/2);
	blank[0] *= 2.0;

	/*** the high degree poly  ***/
	kFCTX(Z[i]+bw, blank + numcoef, scratch, bw, numcoef, 1,
	      lowpoly, lowpoly + bw/2);
	blank[numcoef] *= 2.0;
#else
	/*** the low degree poly  ***/
	memcpy(blank, Z[i], sizeof(double) * bw);
	DCTf( blank, bw, numcoef, CoswSave2 );
	blank[0] *= 2.0;
	/*** the high degree poly  ***/
	memcpy(blank + numcoef, Z[i]+bw, sizeof(double) * bw);
	DCTf( blank+numcoef, bw, numcoef, CoswSave2 );
	blank[numcoef] *= 2.0;
#endif

	/* at any given split point, there is associated
	   with it two pmls: note that I'm saving the
	   coefficients associated with each pair in
	   an interlaced fashion. Again, this is to
	   try to minimize cache-thrashing. */

	Zptr = Z[i];
	for(j = 0 ; j < numcoef ; j += 2){
	  Zptr[0] = blank[j+1];
	  Zptr[1] = blank[j+numcoef];
	  Zptr[2] = blank[j];
	  Zptr[3] = blank[j+numcoef+1];
	  Zptr += 4;
	}
      }

    /* point to array of precomputed fct'd shifted Legendres */
    shiftptr = shiftstore;

    /* set pointer to correct location in the results
       array */
    tmp_result = result_ptr + loc[0] - m - spacing;

    /*************************************************

      now to compute the remaining coefficients by
      applying the simple-split algorithm

      *************************************************/

    /* counter for locations-int array */
    ctr = 0;
    endpt = sw;

    while(ctr <  quant_locs )
      {

	/** set pointer to the cosine coeffients of
	  the weighted data (weighted by p_(l-1) and
	  p_l at the split-point location loc[ctr]) **/

	lowptr = Z[ctr];
	
	/* set lengths of the for-loops: how many
	   coefficients to compute in this split ? */
	revlength = spacing;

	if (ctr != (quant_locs - 1) )
	  forlength = spacing;
	else
	  forlength = bw - loc[ctr];
	
	/* advance endpt for next loop */
	endpt += 2 * spacing;

	/*********************************************
	*********************************************/

	/*** do semi-naive in the reverse direction ***/

	/**********************************************
	*********************************************/

	tmp0 = 0.0 ; tmp1 = 0.0;
	my_ctr = 0;
	for( j = 0 ; j < (revlength - 1)/2 + 1 ; j ++) {
	  tmp0 += shiftptr[my_ctr] * lowptr[my_ctr+my_ctr];
	  tmp1 += shiftptr[my_ctr+1] * lowptr[my_ctr+my_ctr+1];
	  my_ctr += 2;
	}
	shiftptr += 2 * j;

	/*** save result ***/
	*tmp_result++ = times_two * ( tmp0 + tmp1 );

	i = revlength - 1;
	do
	  {

	    tmp0 = tmp1 = tmp2 = tmp3 = 0.0;

	    j = 0;
	    do {
	      tmp0 += shiftptr[j+0] * lowptr[j];
	      tmp1 += shiftptr[j+1] * lowptr[j+1];
	      tmp2 += shiftptr[j+2] * lowptr[j+2];
	      tmp3 += shiftptr[j+3] * lowptr[j+3];
	      j += 4;
	    } while (j < 4 * (i/2 + 1) ) ;

	    shiftptr += j;
	    tmp2 += *shiftptr++ * lowptr[2 * (i + 1) + 2];

	    /*** save results ***/
	    *tmp_result++ = times_two * ( tmp2 + tmp3 );
	    *tmp_result++ = times_two * ( tmp0 + tmp1 );

	    i -= 2;
	  } while (i > 1);

	/*** save result ***/
	*tmp_result++ = times_two * *shiftptr++ * lowptr[2];

	/*********************************************
	**********************************************/

	/*** do semi-naive in the forward direction ***/

	/**********************************************
	**********************************************/

	/*** save result ***/
	*tmp_result++ = times_two * shiftptr[0] * lowptr[1];

	shiftptr ++;

	for(i = 1 ; i < forlength - 1 ; i += 2)
	  {
	    tmp0 = tmp1 = tmp2 = tmp3 = 0.0;

	    j = 0;
	    do
	      {
		tmp3 += shiftptr[j] * lowptr[j];
		tmp2 += shiftptr[j+1] * lowptr[j+1];
		tmp1 += shiftptr[j+2] * lowptr[j+2];
		tmp0 += shiftptr[j+3] * lowptr[j+3];
		j += 4;
	      } while (j < 4 * (i/2 + 1) );
	    shiftptr += j;

	    tmp2 += *shiftptr++ * lowptr[2 * (i + 1) + 1];

	    /*** save result ***/
	    *tmp_result++ = times_two * ( tmp0 + tmp1 );
	    *tmp_result++ = times_two * ( tmp2 + tmp3 );

	}
	
	if ( (forlength % 2) == 0 )
	  {
	    tmp0 = tmp1 = 0.0; 

	    my_ctr = 0;
	    for(j = 0 ; j < (forlength - 1)/2 + 1 ; j ++)
	      {
		tmp0 += shiftptr[my_ctr] * lowptr[my_ctr+my_ctr+3];
		tmp1 += shiftptr[my_ctr+1] * lowptr[my_ctr+my_ctr+2];
		my_ctr += 2;
	      }
	    shiftptr += 2*j;

	    /*** save result ***/
	    *tmp_result++ = times_two * tmp0 + times_two * tmp1;
	  }
	
	/** advance counter, to do the next split-point **/	
	ctr++;
      }
    
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

#ifdef FFTPACK
  free(CoswSave);
  free(CoswSave2);
#endif

}



/******************************************

  This function is just like the one above EXCEPT for
  two things:

  1. It's designed to be used in the forward spherical
     transform.

  2. It lacks the following arguments: sw, shiftstore,
     quant_locs, loops and errflag:

     a) loops and errflag are not needed
     b) sw and quant_locs have been saved in the loc
        array (along with the split point locations,
	unlike in above where they are separate)
     c) by the way things are precomputed in
        test_FST_simsplit_timing.c, at each order,
        the shiftstore array ALWAYS follows the split_pml
	array

     These differences were made with the intention of
     passing as few arguments as possible (since so many
     are already being passed)

  NOTE: data input is partially written over

  m = order of transform
  bw = bandwidth
  data = input pointer to an array of length 2 * bw
  result = output pointer to an array of length bw - m
  Z = pointer to an array of pointers, of length quant_locs
  workspace = array of size 23 * bw (more than enough)
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT
  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform
  cos_pml_table = pointer to precomputed fcts of the Legendre
             functions used in the seminaive portion of the
	     hybrid algorithm
  split_pml_table = pointer to the precomputed Legendre functions
             and the precomputed fcts of the shifted Legendre
	     polynomials which are used in the simple-split
	     portion of the hybrid algorithm

  CoswSave, CoswSave2 = the two precomputed data arrays that are needed
         for the forward (CoswSave) and inverse (CoswSave2) dcts. Two
	 are needed because there are two different sized transforms
	 to do: CoswSave for a size 2*bw dct, CoswSave2 for a size bw dct

  *****/
#ifndef FFTPACK
void FLT_HYBRID_FST( int m,
		     int bw,
		     double *data,
		     double *result,
		     int *loc,
		     double **Z,
		     double *workspace,
		     struct lowhigh *poly,
		     double *cos_pml_table,
		     double *split_pml_table )
#else
void FLT_HYBRID_FST( int m,
		     int bw,
		     double *data,
		     double *result,
		     int *loc,
		     double **Z,
		     double *workspace,
		     double *CoswSave,
		     double *CoswSave2,
		     double *cos_pml_table,
		     double *split_pml_table )
#endif
{
  register int i, j;
  int spacing, numcoef, ctr, endpt;
  int revlength, forlength;
  int n, sw, quant_locs;
  double halfbw;
  double tmp0, tmp1, tmp2, tmp3;
  const double *weights;
  double *lowptr;
  double *blank, *datafct;
  double *scratch;
  double *shift_ptr;
  double *Zptr;
  double *tmp_result;
  double times_two ; /* to adjust in case the order m = 0 */

  /* assign ptrs for the double arrays */
  blank = workspace;
  datafct = blank + bw;
  scratch = datafct + bw;
  halfbw = ((double) bw) * 0.5 ;

  n = 2 * bw;

  /* get switch point, the number of split points, and
     set pointer to where the locations of the split points
     begin */
  sw = loc[0];
  quant_locs = loc[1];
  loc += 2;

  /*** get the weights I'll need ***/
  if (m % 2 == 0)
    weights = get_weights(bw);
  else
    weights = get_oddweights(bw);


  /* get spacing of split points, number of coefficients */
  spacing = loc[0] - sw;
  numcoef = 2;
  while ( numcoef < spacing )
    numcoef *= 2;

  /* fudge factor to normalize coeffs - depends on
     order m */
  if(m == 0)
    times_two = (double) 1;
  else
    times_two = (double) 2;
    
#ifndef FFTPACK
  /* now weight the data   */  
  for(i = 0 ; i < 2 * bw ; i ++)
    data[i] *= weights[i];

  /* take FCT of the weighted data, return only bw-many coeffs */
  kFCTX(data, datafct, scratch, 2 * bw, bw, 1,
	poly, poly + bw);
  
  /* a slower alternative to kFCTX? may depend on the platform that
     this code is run 
     kFCT(data, datafct, scratch, 2 * bw, bw, 1);
     */
#else
  /* note that I can safely use datafct to store
     something of length 2*bw, at least for a
     little bit */
  
  /* now weight the data   */  
  for(i = 0 ; i < 2 * bw ; i ++)
    datafct[i] = data[i] * weights[i];

  /* take FCTX of the weighted data */
  DCTf( datafct, n, bw, CoswSave );

#endif
  
  /* normalize */
  for (i=0; i<bw; i++)
    datafct[i] *= halfbw;
  
  /*
    Ok, the semi-naive part is all set-up,
    Now I have to do the simple-split stuff.
    */

  /* get a subsampled signal back from the coefficients */
#ifndef FFTPACK
  ExpIFCT(datafct, data, scratch, bw, bw, 1);
#else
  memcpy(data, datafct, sizeof(double) * bw);
  DCTb( data, bw, CoswSave2 );
#endif


  /**** now I normalize the first coefficient in datafct ****/
  datafct[0] *= 2.0;

    /*****************************************

      Now to compute the first group of coefficients
      via seminaive

      SemiNaiveReducedX assumes that the data is ALREADY
      cosine-transformed. Otherwise, it's just
      like SemiNaiveReduced.

      NOTE: sw-m is HOW MANY coefficients I want
      to compute !!!
      
      ******************************************/

  SemiNaiveReducedX(datafct, bw, m, sw - m,
		    result,
		    cos_pml_table,
		    blank); 
  
  /*******************************************
    
    Having seminaived the first portion of the coefficients,
    now have to get around to computing the remaining
    coefficients by simple-split
    
    ******************************************/
  

  /* now to pointwise multiply the smoothed data
     with the pmls that I've save at the split points
     (there are two pmls per split point, a lower
     degree one and a higher ) */
  
  for ( i = 0 ; i < quant_locs ; i++)
    {
      Zptr = Z[i];
      for ( j = 0 ; j < bw ; j++)
	{
	  Zptr[j] = split_pml_table[0] * data[j];
	  Zptr[j+bw] =
	    split_pml_table[1] * data[j];
	  split_pml_table += 2;
	}
    }

  /***
    NOTE THAT NOW, given the way I stored the precomputed
    guys, SPLIT_PML_TABLE now points to the beginning of
    the fct'd shifted Legendre polynomials array.
    ***/

  /* now that I've point-wise multiplied, I have
     take the cosine transforms of this data.
     I'll just keep using the same space. */
  
  for (i = 0 ; i < quant_locs ; i++)
    {
#ifndef FFTPACK
      /**** *****/
      kFCTX(Z[i], blank, scratch, bw, numcoef, 1,
	    poly, poly + bw/2);
      blank[0] *= 2.0;
      kFCTX(Z[i]+bw, blank + numcoef, scratch, bw, numcoef, 1,
	    poly, poly + bw/2);
      blank[numcoef] *= 2.0;
      /*****  *****/

      /*** OR a slower alternative to kFCTX ? ***/

      /****
      kFCT(Z[i], blank, scratch, bw, numcoef, 1);
      blank[0] *= 2.0;
      kFCT(Z[i]+bw, blank + numcoef, scratch, bw, numcoef, 1);
      blank[numcoef] *= 2.0;
      *****/
#else
      memcpy(blank, Z[i], sizeof(double) * bw);
      DCTf( blank, bw, numcoef, CoswSave2 );
      blank[0] *= 2.0;
      memcpy(blank + numcoef, Z[i]+bw, sizeof(double) * bw);
      DCTf( blank + numcoef, bw, numcoef, CoswSave2 );
      blank[numcoef] *= 2.0;
#endif


	/* at any given split point, there is associated
	   with it two pmls: note that I'm saving the
	   coefficients associated with each pair in
	   an interlaced fashion. Again, this is to
	   try to minimize cache-thrashing. */

      Zptr = Z[i];
      for(j = 0 ; j < numcoef ; j += 2){
	*Zptr = blank[j+1];
	*(Zptr+1) = blank[j+numcoef];
	*(Zptr+2) = blank[j];
	*(Zptr+3) = blank[j+numcoef+1];
	Zptr+=4;
      }
      
    }
  
  /* set the pointer to the beginning of the array of
     precomputed fct'd shifted Legendres */
  shift_ptr = split_pml_table;

  /* set pointer to correct location in the results
     array */
  tmp_result = result + loc[0] - m - spacing;

  /* counter for locations-int array */
  ctr = 0;
  endpt = sw;

  while(ctr <  quant_locs )
    {
      /** set pointer to the cosine coeffients of
	the weighted data (weighted by p_(l-1) and
	p_l at the current split-point location ) **/
      
      lowptr = Z[ctr];
      
      /* set lengths of the for-loops: how many
	 coefficients to compute in this split ? */
      revlength = loc[ctr] - endpt;
      
      if (ctr != (quant_locs - 1) )
	forlength = spacing;
      else
	forlength = bw - loc[ctr];
      
      /* advance endpt for next loop */
      endpt += 2 * spacing;
      
      /**********************************************/
      /**********************************************/
      
      /*** do semi-naive in the reverse direction ***/
      
      /**********************************************/
      /**********************************************/
	tmp0 = 0.0 ; tmp1 = 0.0;
	i = 0;
	for( j = 0 ; j < (revlength - 1)/2 + 1 ; j ++) {
	  tmp0 += shift_ptr[i] * lowptr[i+i];
	  tmp1 += shift_ptr[i+1] * lowptr[i+i+1];
	  i+=2;
	}
	shift_ptr += 2 * j;

	/* save result */
	*tmp_result++ = times_two * ( tmp0 + tmp1 );

	i = revlength - 1;
	do
	  {
	    
	    tmp0 = 0.0 ; tmp1 = 0.0;
	    tmp2 = 0.0 ; tmp3 = 0.0;

	    j = 0;
	    do {
	      tmp0 += shift_ptr[j] * lowptr[j];
	      tmp1 += shift_ptr[j+1] * lowptr[j+1];
	      tmp2 += shift_ptr[j+2] * lowptr[j+2];
	      tmp3 += shift_ptr[j+3] * lowptr[j+3];
	      j += 4;
	    } while (j < 4 * (i/2 + 1) ) ;

	    shift_ptr += j;
	    tmp2 += shift_ptr[0] * lowptr[2 * (i + 1) + 2];
	    shift_ptr ++;

	    /** save results **/
	    *tmp_result++ = times_two * ( tmp2 + tmp3 );
	    *tmp_result++ = times_two * ( tmp0 + tmp1 );

	    i -= 2;
	  } while (i > 1);

	/** save results **/
	*tmp_result = times_two * shift_ptr[0] * lowptr[2];
	
	shift_ptr ++;
	tmp_result ++;
	
	/**********************************************/
	/**********************************************/
	
	/*** do semi-naive in the forward direction ***/
	
	/**********************************************/
	/**********************************************/

	/* save result */
	*tmp_result++ = times_two * shift_ptr[0] * lowptr[1];
	shift_ptr ++;


	for(i = 1 ; i < forlength - 1 ; i += 2)
	  {


	    tmp0 = 0.0 ; tmp1 = 0.0;
	    tmp2 = 0.0 ; tmp3 = 0.0;

	    j = 0;
	    do
	      {
		tmp3 += shift_ptr[j] * lowptr[j];
		tmp2 += shift_ptr[j+1] * lowptr[j+1];
		tmp1 += shift_ptr[j+2] * lowptr[j+2];
		tmp0 += shift_ptr[j+3] * lowptr[j+3];
		j += 4;
	      } while (j < 4 * (i/2 + 1) );
	    shift_ptr += j;

	    tmp2 += shift_ptr[0] * lowptr[2 * (i + 1) + 1];
	    shift_ptr ++;

	    /** save results **/
	    *tmp_result++ = times_two * ( tmp0 + tmp1 );
	    *tmp_result++ = times_two * ( tmp2 + tmp3 );

	}
	
	if ( (forlength % 2) == 0 )
	  {
	    tmp0 = 0.0 ; tmp1 = 0.0;

	    for(j = 0 ; j < (forlength - 1)/2 + 1 ; j ++)
	      {
		tmp0 += shift_ptr[2*j] * lowptr[4*j+3];
		tmp1 += shift_ptr[2*j+1] * lowptr[4*j+2];
	      }
	    shift_ptr += 2*j;

	    /** save result **/
	    *tmp_result++ = times_two * tmp0 + times_two * tmp1;
	  }
	/** advance counter, to do the next split-point **/
	ctr++;
    }
  
}

