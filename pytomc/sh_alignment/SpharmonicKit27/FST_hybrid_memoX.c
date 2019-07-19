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


/**************************************************************************

  FST_hybrid_memoX.c - routines to perform convolutions on the
  2-sphere using a combination of hybrid and semi-naive algorithms.

  ASSUMES THAT ALL PRECOMPUTED-DATA IS IN MEMORY, AND NOT TO BE
  READ FROM THE DISK.

  FUNCTIONS WILL WRITE OVER INPUT DATA

  The primary functions in this package are

  1) FST_hybrid_memo() - computes the spherical harmonic expansion
                    using the seminaive and hybrid algorithms
  2) InvFST_semi_memo() - computes the inverse spherical harmonic transform
                    using the seminaive algorithm
  3) FZT_hybrid_memo() - computes the zonal harmonic transform
                    using the hyrbid algorithm      
  4) TransMult() - Multiplies harmonic coefficients using Driscoll-Healy
                    result.  Dual of convolution in "time" domain.
  5) Conv2Sphere_hyb_memo() - Convolves two functions defined on the 2-sphere,
                    using seminaive and seminaive/hybrid transforms


   For descriptions on calling these functions, see the documentation
   preceding each function.


*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "FFTcode.h"
#include "cospmls.h"
#include "fftpack.h"
#include "flt_hybrid.h"
#include "primitive.h"
#include "primitive_FST.h"
#include "seminaive.h"

#ifdef FFTPACK
#include "fft_grids.h"
#endif

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))


#ifndef FFTPACK

/**********************************************************************
  Computes the fourier transform of each row of the grid.  This
  is NOT the same as a 2-D Fourier transform.

  Used by forward spherical transform routine.

  Since this will be input to an associated legendre transform,
  the lines of longitude, or zones, are loaded into the rows
  in a transposed fashion for easy access by the Legendre
  transform procedure.  The grid is expected to
  be size * size, which is probably (2*bw) * (2*bw). 

  realgrid, imaggrid - (size x size) arrays of real and imag
                       input

  RESULT WILL BE WRITTEN OVER INPUT ARRAYS
  
  workspace - double pointer to array of (6 * size) = (12 * bw)

  ***************************************************************/

static void grid_fourier(double *realgrid,
			 double *imaggrid,
			 int size,
			 double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */

  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */
 
  for (i=0; i<size; i++) 
    {
      FFTInterp(realgrid+(i*size), imaggrid+(i*size),
		rout, iout, size, size, scratchpad, 1);

      /* store the results     */      
      memcpy(realgrid+(i*size),rout,sizeof(double)*size);
      memcpy(imaggrid+(i*size),iout,sizeof(double)*size);
    }
  
  /* now transpose the results */
  transpose(realgrid,size);
  transpose(imaggrid,size);

}

/**********************************************************************

  same as above except for inverse Fourier transform

  used by inverse spherical transform procedure 

  RESULT WILL BE WRITTEN OVER INPUT ARRAYS


  workspace = (12 * bw)

  *****************************************************************/

static void grid_invfourier(double *realgrid,
			    double *imaggrid,
			    int size,
			    double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */

  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */

  for (i=0; i<size; i++) {
    FFTEval(realgrid+(i*size), imaggrid+(i*size),
	    rout, iout, size, size, scratchpad, 1);
    
    /* now load results into output buffer */
    memcpy(realgrid+(i*size),rout,sizeof(double) * size);
    memcpy(imaggrid+(i*size),iout,sizeof(double) * size);
    
    
  }
  
}

#endif /* FFTPACK */



#ifndef FFTPACK

/*********************************************************************

  Performs a forward spherical harmonic transform using
  the hybrid and semi-naive algorithms.

  Size is the dimension of the input array (size x size) and it is
  expected that size=2*bw.  The inputs rdata and idata are expected
  to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
  to be pointers to bw x bw arrays, and will contain the harmonic
  coefficients in a "linearized" form.

  seminaive_table should be a (double **) pointer to the
  arrays of precomputed data, both for the hybrid and
  seminaive portion of the transform. By its construction,
  these (double **) ptrs, along with the (double **) ptrs
  split_dptr and cos_dptr, are pointing to various locations of
  a SINGLE BLOCK OF MEMORY.

  split_dptr and cos_dptr point to the precomputed data
  which is used in the hybrid algorithm for orders m = 0
  through m = lim. The seminaive_table ptrs point to the
  data that will be used in the seminaive algorithm. This
  arrangement is admittedly awkward and a bit hokey, but
  it allows me to allocate less memory than I would otherwise
  have to do if I decided to have separate arrays for the
  hybrid portion of the algorithm and the seminaive portion
  of the algorithm (given the ways the data is precomputed).


  workspace needs to be a double pointer to an array of size
  26 * bw.
 
  dataformat = 1: input data is real, = 0: input data complex

  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used

  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform

  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT

  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform


  
  */

/* 
   Output Ordering of coeffs f(m,l) is
   f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
                 f(bw-2,bw-2), f(bw-2,bw-1)
		               f(bw-1,bw-1)
			       f(-(bw-1),bw-1)
		 f(-(bw-2),bw-2) f(-(bw-2),bw-1)
	  etc.
	          f(-2,2) ... f(-2,bw-1)
	  f(-1,1) f(-1,2) ... f(-1,bw-1)
    
   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.

*/

/*      dataformat =0 -> samples are complex, =1 -> samples real */

void FST_hybrid_memo(double *rdata, double *idata,
		     double *rcoeffs, double *icoeffs,
		     int size, double *workspace, int dataformat,
		     int lim, double **Z, struct lowhigh *poly,
		     int *loc, double **seminaive_table,
		     double **split_dptr,
		     double **cos_dptr)
{
  int bw, m, i, j;
  double *fltres, *scratchpad;
  int tmpindex[2];

  double *cos_even;
  double pow_one;
  double *rres, *ires;
  int *loc_home, *loc_places;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* assign space */
  if(0)
    {
      rres = workspace;   /* needs (size * size) = (4 * bw^2) */
      ires = rres + (size * size); /* needs (size * size) = (4 * bw^2) */
      fltres = ires + size*size; /* needs bw  */
      cos_even = fltres + bw; /* needs bw */
      scratchpad = cos_even + bw; /* needs (24 * bw)  */
      /* total workspace is (8 * bw^2 + 29 * bw) */

    }
  else
    {
      fltres = workspace; /* needs bw  */
      cos_even = fltres + bw; /* needs bw */
      scratchpad = cos_even + bw; /* needs (24 * bw)  */
    }      


  loc_places = (int *) malloc(sizeof(int) * (lim + 1) );
  loc_home = loc;

  /* do the FFTs along phi */
  grid_fourier(rdata, idata, size, scratchpad);

  for(m = 0 ; m < lim + 1 ; m ++) {
    /*** first get location of file ptr, and locs ptr ***/
    loc_places[m] = (loc - loc_home);

    /* do the real part */
    FLT_HYBRID_FST( m, bw, rdata+(m*size), rcoeffs,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_dptr[m],
		    split_dptr[m]);
    rcoeffs += bw-m;
        
    /* now do the imaginary part */
    FLT_HYBRID_FST( m, bw, idata+(m*size), icoeffs,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_dptr[m],
		    split_dptr[m]);
    
    icoeffs += bw-m;

    /***
      advance to the next sw, quant_locs and locs
      ***/

    loc += 2 + loc[1];

  }

  for (m = lim + 1 ; m < bw ; m ++ ) {
    
    /* do the real part */
    SemiNaiveReduced(rdata+(m*size), 
		     bw, 
		     m, 
		     rcoeffs, 
		     seminaive_table[m],
		     scratchpad,
		     cos_even);
    rcoeffs += bw-m;
    
    /* do imaginary part */
    SemiNaiveReduced(idata+(m*size),
		     bw, 
		     m, 
		     icoeffs, 
		     seminaive_table[m],
		     scratchpad,
		     cos_even);

    icoeffs += bw-m;
    
  }

  
  /*** now do upper coefficients ****/
  
  /* now if the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m)),
     so use that to get the rest of the coefficients
     
     dataformat =0 -> samples are complex, =1 -> samples real
     
     */
  
  if( dataformat == 0 ){
    
    /* note that m is greater than bw here, but this is for
       purposes of indexing the input data arrays.  
       The "true" value of m as a parameter for Pml is
       size - m  */

    for ( m = bw + 1 ; m < size - lim  ; m ++ ) {

      /* do the real part */
      SemiNaiveReduced(rdata+(m*size),
		       bw, 
		       size-m, 
		       rcoeffs, 
		       seminaive_table[size-m],
		       scratchpad,
		       cos_even);
      
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  rcoeffs[i] *= -1.0;
      }
      rcoeffs += m-bw;
      
      /* do the imaginary part */
      SemiNaiveReduced(idata+(m*size),
		       bw, 
		       size-m, 
		       icoeffs, 
		       seminaive_table[size-m],
		       scratchpad,
		       cos_even);
      /* now load imag part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  icoeffs[i] *= -1.0;
      }
      icoeffs += m-bw;

    }
    
    for ( m = size - lim ; m < size ; m ++ ) {


      /*** first set location of file ptr ***/
      loc = loc_home + loc_places[size-m];

      /* do the real part */
      FLT_HYBRID_FST( size-m, bw, rdata+(m*size), rcoeffs,
		      loc,
		      Z, scratchpad,
		      poly,
		      cos_dptr[size-m],
		      split_dptr[size-m]);
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  rcoeffs[i] *= -1.0;
      }
      rcoeffs += m-bw;
      
      /* now do the imaginary part */
      FLT_HYBRID_FST( size-m, bw, idata+(m*size), icoeffs,
		      loc,
		      Z, scratchpad,
		      poly,
		      cos_dptr[size-m],
		      split_dptr[size-m]);
      /* now load imag part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  icoeffs[i] *= -1.0;
      }
      icoeffs += m-bw;
    }

  }
  else        /**** if the data is real ****/
    {            
      /*** have to move the pointers back to the beginning
      of the array; don't want to mess up the indexing
      that follows ***/

      rcoeffs -= (bw * (bw + 1)) / 2;
      icoeffs -= (bw * (bw + 1)) / 2;

      pow_one = 1.0;
      for(i = 1; i < bw; i++){
	pow_one *= -1.0;
	for( j = i; j < bw; j++){	
	  seanindex2(i, j, bw, tmpindex);
	  rcoeffs[tmpindex[1]] =
	    pow_one * rcoeffs[tmpindex[0]];
	  icoeffs[tmpindex[1]] =
	    -1.0 * pow_one * icoeffs[tmpindex[0]];
	}
      }
    }

  free(loc_places);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );
  tmpB = sqrt( 2. * PI );

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	    
	  }
      }
  
}


/**********************************************************************
  Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
  are harmonic coefficients stored in (bw * bw) arrays in the order
  spec'ed above.  rdata and idata are (size x size) arrays with
  the transformed result.  size is expected to be 2 * bw.

  transpose_spharmonic_pml_table should be the (double **) result of a call
  to Transpose_Spharmonic_Pml_Table(). This is the precomputed data
  for the inverse transform.
  
  dataformat = 1: input data is real, = 0: input data complex
  
  workspace is (8 * bw^2) + (33 * bw)

  ***************************************************************/

void InvFST_semi_memo(double *rcoeffs, double *icoeffs, 
		      double *rdata, double *idata, 
		      int size, 
		      double **transpose_seminaive_table,
		      double *workspace,
		      int dataformat)
{
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *rfourdata, *ifourdata;
  double *rinvfltres, *iminvfltres, *scratchpad;
  double *sin_values, *eval_pts;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* allocate space */
  if (0)
    {
      rfourdata = workspace; /* needs (size * size) */
      ifourdata = rfourdata + size * size; /* needs (size * size) */
      rinvfltres = ifourdata + size * size;  /* needs (2 * bw) */
      iminvfltres = rinvfltres + (2 * bw);    /* needs (2 * bw) */
      sin_values = iminvfltres + (2 * bw);    /* needs (2 * bw) */
      eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
      scratchpad = eval_pts + (2 * bw);       /* needs (24 * bw) */
      /* total workspace = (8 * bw^2) + (32 * bw) */
    }
  else
    {
      rinvfltres = workspace;  /* needs (2 * bw) */
      iminvfltres = rinvfltres + (2 * bw);    /* needs (2 * bw) */
      sin_values = iminvfltres + (2 * bw);    /* needs (2 * bw) */
      eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
      scratchpad = eval_pts + (2 * bw);       /* needs (24 * bw) */
      /* total workspace = (32 * bw) */
    }

  /* load up the sin_values array */
  n = 2*bw;
  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /**********************

  New in 2.6: Need to massage the coefficients, to get that
              right factor in there (how embarrassing)

  ***********************/

  tmpA = 1./(2. * sqrt(PI) );
  tmpB = 1./sqrt(2. * PI ) ;

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
      }

  /* Now do all of the inverse Legendre transforms */
  rdataptr = rcoeffs;
  idataptr = icoeffs;
  for (m=0; m<bw; m++) {

    /* fprintf(stderr,"m = %d\n",m); */
      
      /* do real part first */ 
      
      InvSemiNaiveReduced(rdataptr,
			  bw,
			  m,
			  rinvfltres,
			  transpose_seminaive_table[m],
			  sin_values,
			  scratchpad);
      
      /* now do imaginary part */
      
      InvSemiNaiveReduced(idataptr,
			  bw,
			  m,
			  iminvfltres,
			  transpose_seminaive_table[m],
			  sin_values,
			  scratchpad);
      
      /* will store normal, then tranpose before doing inverse fft */
      memcpy(rdata+(m*size), rinvfltres, sizeof(double) * size);
      memcpy(idata+(m*size), iminvfltres, sizeof(double) * size);
      
      /* move to next set of coeffs */
      
      rdataptr += bw-m;
      idataptr += bw-m;
      
  }
  
  /* now fill in zero values where m = bw (from problem definition) */
  memset(rdata + (bw * size), 0, sizeof(double) * size);
  memset(idata + (bw * size), 0, sizeof(double) * size);
  
  
  /* now if the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     invf-hat(l,-m) = conjugate(invf-hat(l,m)),
     so use that to get the rest of the real data
     
     dataformat =0 -> samples are complex, =1 -> samples real
     
     */
  if(dataformat == 0){
    
    /* now do negative m values */
    
    for (m=bw+1; m<size; m++) 
      {
	/*
	  fprintf(stderr,"m = %d\n",-(size-m));
	  */
	
	/* do real part first */
	
	InvSemiNaiveReduced(rdataptr,
			    bw,
			    size - m,
			    rinvfltres,
			    transpose_seminaive_table[size - m],
			    sin_values,
			    scratchpad);
	
	/* now do imaginary part */
	
	InvSemiNaiveReduced(idataptr,
			    bw,
			    size - m,
			    iminvfltres,
			    transpose_seminaive_table[size - m],
			    sin_values,
			    scratchpad);
	

	/* will store normal, then tranpose before doing inverse fft    */
	if ((m % 2) != 0)
	  for(i=0; i< size; i++){
	    rinvfltres[i] *= -1.0;
	    iminvfltres[i] *= -1.0;
	  }
	
	memcpy(rdata + (m*size), rinvfltres, sizeof(double) * size);
	memcpy(idata + (m*size), iminvfltres, sizeof(double) * size);
	
	
	/* move to next set of coeffs */
	rdataptr += bw-(size-m);
	idataptr += bw-(size-m);
	
      } /* closes m loop */
  }
  else {
    for(m = bw + 1; m < size; m++){

      memcpy(rdata+(m*size), rdata+((size-m)*size), sizeof(double) * size);
      memcpy(idata+(m*size), idata+((size-m)*size), sizeof(double) * size);
      for(i = 0; i < size; i++)
	idata[(m*size)+i] *= -1.0;
    }

  }

  /** now transpose **/
  transpose(rdata, size);
  transpose(idata, size);
  
  
  /* now do inverse fourier grid computation */
  grid_invfourier(rdata, idata, size, scratchpad);
  
  /* amscray */
  
}
/*********************************************************************

  Zonal Harmonic transform using hybrid algorithm - used in convolutions
  rdata and idata should be pointers to size x size arrays.
  rres and ires should be pointers to double arrays of size bw.
  size = 2 * bw
 
  FZT_hybrid_memo only computes spherical harmonics for m=0.

  dataformat = 1: input data is real, = 0: input data complex
  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used
  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT
  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform
  disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.

  
  workspace needed is (12 * bw) 

*/

void FZT_hybrid_memo(double *rdata, double *idata,
		     double *rres, double *ires,
		     int size,
		     double *workspace, int dataformat,
		     int lim, double **Z, struct lowhigh *poly,
		     int *loc, double **seminaive_table,
		     double **split_dptr,
		     double **cos_dptr)
{
  
  int i, j, bw;
  double *r0, *i0, dsize;
  double tmpreal, tmpimag;
  double *scratchpad;
  double *cos_even;
  int l, dummy ;
  double tmpA ;


  bw = size/2;

  /* assign memory */

  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (8 * bw) */
                              /* total workspace = 13*bw */

  dsize = (double) size;
  /* compute the m = 0 components */

  for (i=0; i<size; i++) {
    tmpreal = 0.0;
    tmpimag = 0.0;
    for(j=0; j<size; j++) {
      tmpreal += rdata[(i*size)+j];
      tmpimag += idata[(i*size)+j];
    }
    /* normalize */
    r0[i] = tmpreal/dsize;
    i0[i] = tmpimag/dsize;
  }

  /* do the real part */
  FLT_HYBRID_FST( 0, bw, r0, rres,
		  loc,
		  Z, scratchpad,
		  poly,
		  cos_dptr[0],
		  split_dptr[0]);
 
  if(dataformat == 0)   /* do imaginary part */
    FLT_HYBRID_FST( 0, bw, i0, ires,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_dptr[0],
		    split_dptr[0]);
  else                 /* otherwise set coefficients = 0 */
    memset(ires, 0, sizeof(double) * size);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );

  for(l=0;l<bw;l++)
    {
      dummy = seanindex(0,l,bw);

      rres[dummy] *= tmpA ;
      ires[dummy] *= tmpA ;
    }

}
/************************************************************************/
/* multiplies harmonic coefficients of a function and a filter */
/* See convolution theorem of Driscoll and Healy */
/* datacoeffs should be output of an FST, filtercoeffs the 
   output of an FZT.  There should be (bw * bw) datacoeffs,
   and bw filtercoeffs.
   rres and ires should point to arrays
   of dimension bw * bw. size parameter is 2*bw  
*/

void TransMult(double *rdatacoeffs, double *idatacoeffs, 
	       double *rfiltercoeffs, double *ifiltercoeffs, 
	       double *rres, double *ires,
	       int size)
{

  int m, l, bw;
  double *rdptr, *idptr, *rrptr, *irptr;

  bw = size/2;

  rdptr = rdatacoeffs;
  idptr = idatacoeffs;
  rrptr = rres;
  irptr = ires;

  for (m=0; m<bw; m++) {
    for (l=m; l<bw; l++) {
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-m], idptr[l-m],
	       rrptr[l-m], irptr[l-m]);

      rrptr[l-m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += bw-m; idptr += bw-m;
    rrptr += bw-m; irptr += bw-m;
  }
  for (m=bw+1; m<size; m++) {
    for (l=size-m; l<bw; l++){
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-size+m], idptr[l-size+m],
	       rrptr[l-size+m], irptr[l-size+m]);

      rrptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += m-bw; idptr += m-bw;
    rrptr += m-bw; irptr += m-bw;
  }

}

/************************************************************************
  Here's the big banana
  Convolves two functions defined on the 2-sphere.
  Uses hybrid and seminaive algorithms for spherical harmonic transforms
  size = 2*bw
  

  Inputs:
  rdata, idata - (size * size) arrays containing real and
               imaginary parts of sampled function.
  rfilter, ifilter - (size * size) arrays containing real and
               imaginary parts of sampled filter function.
  rres, ires - (size * size) arrays containing real and
                  imaginary parts of result function.
  size - should be 2 * bw

  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used
  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform
  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT


  seminaive_table should be a (double **) pointer to the
  arrays of precomputed data, both for the hybrid and
  seminaive portion of the transform. By its construction,
  these (double **) ptrs, along with the (double **) ptrs
  split_dptr and cos_dptr, are pointing to various locations of
  a SINGLE BLOCK OF MEMORY.

  split_dptr and cos_dptr point to the precomputed data
  which is used in the hybrid algorithm for orders m = 0
  through m = lim. The seminaive_table ptrs point to the
  data that will be used in the seminaive algorithm. This
  arrangement is admittedly awkward and a bit hokey, but
  it allows me to allocate less memory than I would otherwise
  have to do if I decided to have separate arrays for the
  hybrid portion of the algorithm and the seminaive portion
  of the algorithm (given the ways the data is precomputed).

  transpose_spharmonic_pml_table should be the (double **) result of a call
  to Transpose_Spharmonic_Pml_Table(). This is the precomputed data
  for the inverse transform.


  workspace needs to be a double pointer to an array of size
	16 * bw * bw + 33 * bw.



   ASSUMPTIONS:
   1. data is strictly REAL

*/

void Conv2Sphere_hyb_memo( double *rdata, double *idata,
			   double *rfilter, double *ifilter,
			   double *rres, double *ires,
			   int size, int lim,
			   int *loc, double **Z,
			   struct lowhigh *poly,
			   double **seminaive_table,
			   double **trans_seminaive_table,
			   double **split_dptr,
			   double **cos_dptr,
			   double *workspace)
{
  int bw;
  double *frres, *fires, *filtrres, *filtires, *trres, *tires;
  double *scratchpad;

  bw = size / 2;

  frres = workspace;   /* needs (bw*bw) */
  fires = frres + (bw*bw);    /* needs (bw*bw) */
  trres = fires + (bw*bw);    /* needs (bw*bw) */
  tires = trres + (bw*bw);    /* needs (bw*bw) */
  filtrres = tires + (bw*bw); /* needs bw */
  filtires = filtrres + bw;   /* needs bw */
  scratchpad = filtires + bw; /* needs (8 * bw^2) + (32 * bw) */


  FST_hybrid_memo( rdata, idata,
		   frres, fires,
		   size,
		   scratchpad, 1,
		   lim, Z, poly,
		   loc, seminaive_table,
		   split_dptr, cos_dptr);
    
  FZT_hybrid_memo( rfilter, ifilter,
		   filtrres, filtires,
		   size,
		   scratchpad, 1,
		   lim, Z, poly,
		   loc, seminaive_table, 
		   split_dptr, cos_dptr);

  TransMult(frres, fires, filtrres, filtires, trres,
	    tires, size);

  InvFST_semi_memo(trres, tires,
		   rres, ires,
		   size, trans_seminaive_table,
		   scratchpad, 1);


}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

#else

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

/*********************************************************************

  Performs a forward spherical harmonic transform using
  the hybrid and semi-naive algorithms.

  Size is the dimension of the input array (size x size) and it is
  expected that size=2*bw.  The inputs rdata and idata are expected
  to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
  to be pointers to bw x bw arrays, and will contain the harmonic
  coefficients in a "linearized" form.

  seminaive_table should be a (double **) pointer to the
  arrays of precomputed data, both for the hybrid and
  seminaive portion of the transform. By its construction,
  these (double **) ptrs, along with the (double **) ptrs
  split_dptr and cos_dptr, are pointing to various locations of
  a SINGLE BLOCK OF MEMORY.

  split_dptr and cos_dptr point to the precomputed data
  which is used in the hybrid algorithm for orders m = 0
  through m = lim. The seminaive_table ptrs point to the
  data that will be used in the seminaive algorithm. This
  arrangement is admittedly awkward and a bit hokey, but
  it allows me to allocate less memory than I would otherwise
  have to do if I decided to have separate arrays for the
  hybrid portion of the algorithm and the seminaive portion
  of the algorithm (given the ways the data is precomputed).


  workspace needs to be a double pointer to an array of size
  26 * bw.
 
  dataformat = 1: input data is real, = 0: input data complex

  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used

  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform

  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT

  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform


  
  */

/* 
   Output Ordering of coeffs f(m,l) is
   f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
                 f(bw-2,bw-2), f(bw-2,bw-1)
		               f(bw-1,bw-1)
			       f(-(bw-1),bw-1)
		 f(-(bw-2),bw-2) f(-(bw-2),bw-1)
	  etc.
	          f(-2,2) ... f(-2,bw-1)
	  f(-1,1) f(-1,2) ... f(-1,bw-1)
    
   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.

*/

/* data assumed to be real */

void FST_hybrid_memo(double *rdata,
		     double *rcoeffs,
		     double *icoeffs,
		     int size,
		     int lim,
		     int *loc,
		     double **Z,
		     double *workspace,
		     double **seminaive_table,
		     double **split_dptr,
		     double **cos_dptr,
		     double *wSave,
		     double *CoswSave,
		     double *CoswSave2 )
{
  int bw, m, i, j;
  double *fltres, *scratchpad;
  int tmpindex[2];
  double *cos_even;
  double pow_one;
  double *rres;
  int *loc_home, *loc_places;
  double *tmprdata_ptr;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* assign space */
  cos_even = workspace; /* needs bw */
  scratchpad = cos_even + bw; /* needs (24 * bw)  */


  loc_places = (int *) malloc(sizeof(int) * (lim + 1) );
  loc_home = loc;

  /* do the FFTs along phi */
  grid_fourier_FFTPACK(rdata, size, wSave);
  tmprdata_ptr = rdata;

  for(m = 0 ; m < lim + 1 ; m ++) {
    /*** first get location of file ptr, and locs ptr ***/
    loc_places[m] = (loc - loc_home);

    /* do the real part */
    FLT_HYBRID_FST( m,
		    bw,
		    tmprdata_ptr,
		    rcoeffs,
		    loc,
		    Z,
		    scratchpad,
		    CoswSave,
		    CoswSave2,
		    cos_dptr[m],
		    split_dptr[m]);
    rcoeffs += bw-m;
    tmprdata_ptr += size;
        
    if ( m != 0 )
      {
	/* now do the imaginary part */
	FLT_HYBRID_FST( m,
			bw,
			tmprdata_ptr,
			icoeffs,
			loc,
			Z,
			scratchpad,
			CoswSave,
			CoswSave2,
			cos_dptr[m],
			split_dptr[m]);
	
	icoeffs += bw-m;
	tmprdata_ptr += size;
      }
    else
      {
	memset(icoeffs, 0, sizeof(double) * bw );
	icoeffs += bw;
      }

    /***
      advance to the next sw, quant_locs and locs
      ***/

    loc += 2 + loc[1];

  }

  for (m = lim + 1 ; m < bw ; m ++ ) {
    
    /* do the real part */
    SemiNaiveReduced(tmprdata_ptr,
		     bw, 
		     m, 
		     rcoeffs, 
		     seminaive_table[m],
		     scratchpad,
		     cos_even,
		     CoswSave);
    rcoeffs += bw-m;
    tmprdata_ptr += size;

    if ( m != 0 )
      {
	/* do imaginary part */
	SemiNaiveReduced(tmprdata_ptr,
			 bw, 
			 m, 
			 icoeffs, 
			 seminaive_table[m],
			 scratchpad,
			 cos_even,
			 CoswSave);
	
	icoeffs += bw-m;
	tmprdata_ptr += size;
      }
    else
      {
	memset(icoeffs, 0, sizeof(double) * bw);
	icoeffs += bw;
      }
    
  }

  
  /*** now do upper coefficients ****/

  /*** have to move the pointers back to the beginning
    of the array; don't want to mess up the indexing
    that follows ***/
  
  rcoeffs -= (bw * (bw + 1)) / 2;
  icoeffs -= (bw * (bw + 1)) / 2;
  
  pow_one = 1.0;
  for(i = 1; i < bw; i++)
    {
      pow_one *= -1.0;
      for( j = i; j < bw; j++)
	{	
	  seanindex2(i, j, bw, tmpindex);
	  rcoeffs[tmpindex[1]] =
	    pow_one * rcoeffs[tmpindex[0]];
	  icoeffs[tmpindex[1]] =
	    -1.0 * pow_one * icoeffs[tmpindex[0]];
	}
    }


  free(loc_places);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );
  tmpB = sqrt( 2. * PI );

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	    
	  }
      }
  
}


/**********************************************************************
  Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
  are harmonic coefficients stored in (bw * bw) arrays in the order
  spec'ed above.  rdata and idata are (size x size) arrays with
  the transformed result.  size is expected to be 2 * bw.

  transpose_spharmonic_pml_table should be the (double **) result of a call
  to Transpose_Spharmonic_Pml_Table(). This is the precomputed data
  for the inverse transform.
  
  dataformat = 1: input data is real, = 0: input data complex
  
  workspace is (8 * bw^2) + (33 * bw)

  ***************************************************************/

void InvFST_semi_memo(double *rcoeffs,
		      double *icoeffs, 
		      double *rdata,
		      int size, 
		      double **transpose_seminaive_table,
		      double *workspace,
		      double *wSave )
{
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *scratchpad;
  double *sin_values, *eval_pts;
  double *tmp_ptr;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* allocate space */
  sin_values = workspace ;    /* needs (2 * bw) */
  eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
  scratchpad = eval_pts + (2 * bw);       /* needs (16 * bw) */



  /* load up the sin_values array */
  n = 2*bw;
  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 1./(2. * sqrt(PI) );
  tmpB = 1./sqrt(2. * PI ) ;

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
      }

  tmp_ptr = rdata;

  /* Now do all of the inverse Legendre transforms */
  rdataptr = rcoeffs;
  idataptr = icoeffs;
  for (m=0; m<bw; m++) {

    /* fprintf(stderr,"m = %d\n",m); */
      
      /* do real part first */ 
      
      InvSemiNaiveReduced(rdataptr,
			  bw,
			  m,
			  tmp_ptr,
			  transpose_seminaive_table[m],
			  sin_values,
			  scratchpad);
      tmp_ptr += size;
      /* move to next set of coeffs */
      rdataptr += bw - m;

      if ( m != 0 )
	{
	  /* now do imaginary part */
	  
	  InvSemiNaiveReduced(idataptr,
			      bw,
			      m,
			      tmp_ptr,
			      transpose_seminaive_table[m],
			      sin_values,
			      scratchpad);
	  tmp_ptr += size;
	  /* move to next set of coeffs */
	  idataptr += bw - m;
	}
      else
	{
	  /* move to next set of coeffs */
	  idataptr += bw - m;
	}	  
  }
  
  /* now fill in zero values where m = bw (from problem definition) */
  memset(tmp_ptr, 0, sizeof(double) * size);


  /** now transpose **/
  transpose(rdata, size);
  
  /* now do inverse fourier grid computation */
  grid_invfourier_FFTPACK(rdata, size, wSave );
  
  /* amscray */
  
}
/*********************************************************************

  Zonal Harmonic transform using hybrid algorithm - used in convolutions
  rdata and idata should be pointers to size x size arrays.
  rres and ires should be pointers to double arrays of size bw.
  size = 2 * bw
 
  FZT_hybrid_memo only computes spherical harmonics for m=0.

  dataformat = 1: input data is real, = 0: input data complex
  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used
  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT
  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform
  disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.

  
  workspace needed is (12 * bw) 

*/

void FZT_hybrid_memo(double *rdata,
		     double *rres,
		     double *ires,
		     int size,
		     int lim,
		     int *loc,
		     double **Z,
		     double *workspace,
		     double **seminaive_table,
		     double **split_dptr,
		     double **cos_dptr,
		     double *CoswSave,
		     double *CoswSave2 )
{
  
  int i, j, bw;
  double *r0, *i0, dsize;
  double tmpreal, tmpimag;
  double *scratchpad;
  double *cos_even;
  int l, dummy ;
  double tmpA ;

  bw = size/2;

  /* assign memory */

  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (8 * bw) */
                              /* total workspace = 13*bw */

  dsize = (double) size;
  /* compute the m = 0 components */

  for (i=0; i<size; i++) {
    tmpreal = 0.0;
    for(j=0; j<size; j++) {
      tmpreal += rdata[(i*size)+j];
    }
    /* normalize */
    r0[i] = tmpreal/dsize;
  }

  /* do the real part */
  FLT_HYBRID_FST( 0,
		  bw,
		  r0,
		  rres,
		  loc,
		  Z,
		  scratchpad,
		  CoswSave,
		  CoswSave2,
		  cos_dptr[0],
		  split_dptr[0]);

  /* set imaginary part = to zero */
  memset(ires, 0, sizeof(double) * size);
 
  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );

  for(l=0;l<bw;l++)
    {
      dummy = seanindex(0,l,bw);

      rres[dummy] *= tmpA ;
    }

}
/************************************************************************/
/* multiplies harmonic coefficients of a function and a filter */
/* See convolution theorem of Driscoll and Healy */
/* datacoeffs should be output of an FST, filtercoeffs the 
   output of an FZT.  There should be (bw * bw) datacoeffs,
   and bw filtercoeffs.
   rres and ires should point to arrays
   of dimension bw * bw. size parameter is 2*bw  
*/

void TransMult(double *rdatacoeffs, double *idatacoeffs, 
	       double *rfiltercoeffs, double *ifiltercoeffs, 
	       double *rres, double *ires,
	       int size)
{

  int m, l, bw;
  double *rdptr, *idptr, *rrptr, *irptr;

  bw = size/2;

  rdptr = rdatacoeffs;
  idptr = idatacoeffs;
  rrptr = rres;
  irptr = ires;

  for (m=0; m<bw; m++)
    {
      for (l=m; l<bw; l++)
	{
	  compmult(rfiltercoeffs[l], ifiltercoeffs[l],
		   rdptr[l-m], idptr[l-m],
		   rrptr[l-m], irptr[l-m]);

	  rrptr[l-m] *= sqrt(4*M_PI/(2*l+1));
	  irptr[l-m] *= sqrt(4*M_PI/(2*l+1));


	}
      rdptr += bw-m; idptr += bw-m;
      rrptr += bw-m; irptr += bw-m;
    }
  for (m=bw+1; m<size; m++)
    {
      for (l=size-m; l<bw; l++)
	{
	  compmult(rfiltercoeffs[l], ifiltercoeffs[l],
		   rdptr[l-size+m], idptr[l-size+m],
		   rrptr[l-size+m], irptr[l-size+m]);

	  rrptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));
	  irptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));


	}
      rdptr += m-bw; idptr += m-bw;
      rrptr += m-bw; irptr += m-bw;
    }

}

/************************************************************************
  Here's the big banana
  Convolves two functions defined on the 2-sphere.
  Uses hybrid and seminaive algorithms for spherical harmonic transforms
  size = 2*bw
  

  Inputs:
  rdata, idata - (size * size) arrays containing real and
               imaginary parts of sampled function.
  rfilter, ifilter - (size * size) arrays containing real and
               imaginary parts of sampled filter function.
  rres, ires - (size * size) arrays containing real and
                  imaginary parts of result function.
  size - should be 2 * bw

  lim = last (i.e. highest) order that the hybrid Legendre
       transform algorithm will be used
  loc = pointer to array containing locations of switch points,
        split-points, and quant_locs for all the orders that
	the hybrid transform will be used when performing the
	forward spherical transform
  Z = pointer to an array of pointers, of length (max_splits + 1);
      these pointers will be pointing to double arrays, each
      of length 2 * bw, that are used for temporary storage in
      the hybrid Legendre transform
  poly = array of type struct lowhigh of length 4 * bw;
         used in the FCT


  seminaive_table should be a (double **) pointer to the
  arrays of precomputed data, both for the hybrid and
  seminaive portion of the transform. By its construction,
  these (double **) ptrs, along with the (double **) ptrs
  split_dptr and cos_dptr, are pointing to various locations of
  a SINGLE BLOCK OF MEMORY.

  split_dptr and cos_dptr point to the precomputed data
  which is used in the hybrid algorithm for orders m = 0
  through m = lim. The seminaive_table ptrs point to the
  data that will be used in the seminaive algorithm. This
  arrangement is admittedly awkward and a bit hokey, but
  it allows me to allocate less memory than I would otherwise
  have to do if I decided to have separate arrays for the
  hybrid portion of the algorithm and the seminaive portion
  of the algorithm (given the ways the data is precomputed).

  transpose_spharmonic_pml_table should be the (double **) result of a call
  to Transpose_Spharmonic_Pml_Table(). This is the precomputed data
  for the inverse transform.


  workspace needs to be a double pointer to an array of size
	16 * bw * bw + 33 * bw.



   ASSUMPTIONS:
   1. data is strictly REAL

*/

void Conv2Sphere_hyb_memo( double *rdata,
			   double *rfilter,
			   double *rres,
			   int size,
			   int lim,
			   int *loc,
			   double **Z,
			   double *workspace,
			   double **seminaive_table,
			   double **trans_seminaive_table,
			   double **split_dptr,
			   double **cos_dptr,
			   double *wSave,
			   double *CoswSave,
			   double *CoswSave2 )
{
  int bw;
  double *frres, *fires, *filtrres, *filtires, *trres, *tires;
  double *scratchpad;

  bw = size / 2;

  frres = workspace;   /* needs (bw*bw) */
  fires = frres + (bw*bw);    /* needs (bw*bw) */
  trres = fires + (bw*bw);    /* needs (bw*bw) */
  tires = trres + (bw*bw);    /* needs (bw*bw) */
  filtrres = tires + (bw*bw); /* needs bw */
  filtires = filtrres + bw;   /* needs bw */
  scratchpad = filtires + bw; /* needs (8 * bw^2) + (32 * bw) */


  FST_hybrid_memo( rdata,
		   frres,
		   fires,
		   size,
		   lim,
		   loc,
		   Z,
		   scratchpad,
		   seminaive_table,
		   split_dptr,
		   cos_dptr,
		   wSave,
		   CoswSave,
		   CoswSave2);
    
  FZT_hybrid_memo( rfilter,
		   filtrres,
		   filtires,
		   size,
		   lim,
		   loc,
		   Z,
		   scratchpad,
		   seminaive_table, 
		   split_dptr,
		   cos_dptr,
		   CoswSave,
		   CoswSave2);

  TransMult(frres, fires, filtrres, filtires, trres,
	    tires, size);

  InvFST_semi_memo(trres, tires,
		   rres, 
		   size,
		   trans_seminaive_table,
		   scratchpad,
		   wSave );


}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

#endif
