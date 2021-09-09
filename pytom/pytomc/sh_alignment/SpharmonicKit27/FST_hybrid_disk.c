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


/************************************************************

  FST_hybrid_disk.c - routines to perform convolutions on the
  2-sphere using a combination of hybrid and semi-naive algorithms.


  ASSUMES THAT ALL PRECOMPUTED-DATA IS TO BE READ FROM THE
  DISK!!!!!

  The primary functions in this package are

  1) FST_hybrid_disk() - computes the spherical harmonic expansion
                    using the seminaive and hybrid algorithms
  2) InvFST_semi_disk() - computes the inverse spherical harmonic transform
                    using the seminaive algorithm
  3) FZT_hybrid_disk() - computes the zonal harmonic transform
                    using the hyrbid algorithm      
  4) TransMult() - Multiplies harmonic coefficients using Driscoll-Healy
                    result.  Dual of convolution in "time" domain.
  5) Conv2Sphere_hyb_disk() - Convolves two functions defined on the 2-sphere,
                    using seminaive and seminaive/hybrid transforms

  For descriptions on calling these functions, see the documentation
  preceding each function.


  ***********************************************************/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "cospmls.h"
#include "fft_grids.h"
#include "flt_hybrid.h"
#include "primitive.h"
#include "primitive_FST.h"
#include "seminaive.h"

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))

#ifndef FFTPACK


/*********************************************************************


  Performs a spherical harmonic transform using the hybrid and
  seminaive algorithms

  Size is the dimension of the input array (size x size) and it is
  expected that size=2*bw.  The inputs rdata and idata are expected
  to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
  to be pointers to bw x bw arrays, and will contain the harmonic
  coefficients in a "linearized" form.
   
  workspace needs to be a double pointer to an array of size
  8 * bw * bw + 33 * bw.
 
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

void FST_hybrid_disk( double *rdata, double *idata,
		      double *rcoeffs, double *icoeffs,
		      int size, double *workspace,
		      int dataformat, int lim,
		      double **Z, struct lowhigh *poly,
		      int *loc, double *disk_array)
{
  FILE *fp;
  int bw, m, i, j;
  double *fltres, *scratchpad;
  int tmpindex[2];
  double *cos_even;
  double pow_one;
  long int *ptr_locs;
  int data_amount, data_amount2;
  double *cos_data, *split_data;
  int *loc_home, *loc_places;
  double *rres, *ires;
  char strm[128], tmpstr[128];
  int l, dummy ;
  double tmpA, tmpB ;

  cos_data = disk_array;
  bw = size/2;

  /* assign space */
  rres = workspace;   /* needs (size * size) = (4 * bw^2) */
  ires = rres + (size * size);
  fltres = ires + (size * size); /* needs bw  */
  cos_even = fltres + bw; /* needs bw */
  scratchpad = cos_even + bw; /* needs (24 * bw)  */
 
  /* total workspace is (29 * bw) */

  ptr_locs = (long int *) malloc(sizeof(long int) * bw);
  loc_places = (int *) malloc(sizeof(int) * (lim + 1) );
  loc_home = loc;

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Hybrid_bw%d.dat", bw);
  strcat(strm, tmpstr);

  fp = fopen(strm, "r");
  if ( fp == NULL )
    {
      fprintf( stderr,"Error in FST_hybrid_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }


  /* do the FFTs along phi */
  grid_fourier(rdata, idata, rres, ires, size, scratchpad);


  /* do the legendre transforms for orders m = 0
     through m = lim via the hybrid algorithm */
  
  for(m = 0 ; m < lim + 1 ; m ++) {

    /*** have to copy data to arrays ***/
    
    /*** first get location of file ptr, and locs ptr;
    I have to do this in case the input data is complex, it
    makes the book-keeping involved with doing the negative
    orders easier***/
    
    ptr_locs[m] = ftell( fp ) ;
    loc_places[m] = (loc - loc_home);

    /* first get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );

    /* assign the split_ptr now */
    split_data = cos_data + data_amount ;
    /* how many elements for the split array ? */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* read in that many */
    fread ( split_data, sizeof(double) , data_amount, fp );

    /* ditto for the fct'd shifted legendres; note where I'm
       saving the shifted data ... FLT_HYBRID_FST
       expects this Legendres to immediately follow the
       splits */
    fread ( &data_amount2, sizeof(int), 1, fp );
    fread ( split_data + data_amount, sizeof(double),
	    data_amount2 , fp );


    /* do the real part */
    FLT_HYBRID_FST( m, bw, rres+(m*size), rcoeffs,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_data,
		    split_data);
    
    rcoeffs += bw-m;
    
    /* now do the imaginary part */
    FLT_HYBRID_FST( m, bw, ires+(m*size), icoeffs,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_data,
		    split_data);

    icoeffs += bw-m;

    /***
      advance to the next switch point, quant_locs and locs
      ***/
    loc += 2 + loc[1];

  }


  /** now do the Legendre transforms with seminaive
    algorithm for orders m = lim + 1 through bw - 1 **/

  for (m = lim + 1 ; m < bw ; m ++ ) {

    /*** have to copy data to arrays ***/

    /*** first get location of file ptr; I have to do
      this in case the input data is complex, it
      makes the book-keeping involved with doing the negative
      orders easier***/
   
    /*** first get location of file pointer ***/
    ptr_locs[m] = ftell( fp );

    /* first get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );

    /* do the real part */
    SemiNaiveReduced(rres+(m*size), 
		     bw, 
		     m, 
		     rcoeffs, 
		     cos_data,
		     scratchpad,
		     cos_even);

    rcoeffs += bw-m;
    
    /* do imaginary part */
    SemiNaiveReduced(ires+(m*size),
		     bw, 
		     m, 
		     icoeffs, 
		     cos_data,
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

      /*** first set location of file pointer ***/
      fseek( fp , ptr_locs[size-m] , 0);
      
      /* now get how many elements in cos_array */
      fread ( &data_amount, sizeof(int), 1, fp );
      /* now read in that amount */
      fread (cos_data, sizeof(double), data_amount, fp );
      
      /* do the real part */
      SemiNaiveReduced(rres+(m*size),
		       bw, 
		       size-m, 
		       rcoeffs, 
		       cos_data,
		       scratchpad,
		       cos_even);
      
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  rcoeffs[i] *= -1.0;
      }
      rcoeffs += m-bw;
      
      /* do the imaginary part */
      SemiNaiveReduced(ires+(m*size),
		       bw, 
		       size-m, 
		       icoeffs, 
		       cos_data,
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
      fseek ( fp , ptr_locs[size-m] , 0);
      loc = loc_home + loc_places[size-m];

      /* first get how many elements in cos_array */
      fread ( &data_amount, sizeof(int), 1, fp );
      /* now read in that amount */
      fread (cos_data, sizeof(double), data_amount, fp );
      
      
      /* assign the split_ptr now */
      split_data = cos_data + data_amount ;
      /* how many elements for the split array ? */
      fread ( &data_amount, sizeof(int), 1, fp );
      /* read in that many */
      fread ( split_data, sizeof(double) , data_amount, fp );
      
      /* ditto for the fct'd shifted legendres; note where I'm
	 saving the shifted data ... FLT_HYBRID_FST
	 expects this Legendres to immediately follow the
	 splits */
      fread ( &data_amount2, sizeof(int), 1, fp );
      fread ( split_data + data_amount, sizeof(double),
	      data_amount2 , fp );
      

      /* do the real part */
      FLT_HYBRID_FST( size-m, bw, rres+(m*size), rcoeffs,
		      loc,
		      Z, scratchpad,
		      poly,
		      cos_data,
		      split_data);
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  rcoeffs[i] *= -1.0;
      }
      rcoeffs += m-bw;
      
      /* now do the imaginary part */
      FLT_HYBRID_FST( size-m, bw, ires+(m*size), icoeffs,
		      loc,
		      Z, scratchpad,
		      poly,
		      cos_data,
		      split_data);
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

  fclose(fp);
  free(loc_places);
  free(ptr_locs);

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


/************************************************************************/
/* Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
   are harmonic coefficients stored in (bw * bw) arrays in the order
   spec'ed above.  rdata and idata are (size x size) arrays with
   the transformed result.  size is expected to be 2 * bw.

   workspace is 8 * bw * bw + 33 * bw

   dataformat = 1: input data is real, = 0: input data complex

   disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.

   

*/


void InvFST_semi_disk(double *rcoeffs, double *icoeffs, 
		      double *rdata, double *idata, 
		      int size, double *workspace,
		      int dataformat, double *disk_array)
{
  FILE *fp;
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *rinvfltres, *iminvfltres, *scratchpad;
  double *sin_values, *eval_pts;
  double *rfourdata, *ifourdata;
  long int *ptr_locs;
  double *cos_data;
  int data_amount;
  char strm[128], tmpstr[128];
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;
  cos_data = disk_array;

  /* allocate space */

  ptr_locs = (long int *) malloc(sizeof(long int) * bw);

  rfourdata = workspace; /* needs (size * size) */
  ifourdata = rfourdata + size * size; /* needs (size * size) */
  rinvfltres = ifourdata + size * size ;  /* needs (2 * bw) */
  iminvfltres = rinvfltres + (2 * bw);    /* needs (2 * bw) */
  sin_values = iminvfltres + (2 * bw);    /* needs (2 * bw) */
  eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
  scratchpad = eval_pts + (2 * bw);       /* needs (24 * bw) */
  
  /* total workspace = (8 * bw^2) + (32 * bw) */

  /* load up the sin_values array */
  n = 2*bw;
  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm,PRECOMP_DIR);
  sprintf(tmpstr,"InvSemi_bw%d.dat", bw);
  strcat(strm, tmpstr);

  fp = fopen(strm, "r");

  if ( fp == NULL )
    {
      fprintf( stderr,"Error in InvFST_semi_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }

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

    /*** have to copy data to arrays ***/
    
    /*** first get location of file pointer ***/
    ptr_locs[m] = ftell( fp );
    /* get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );


      /* do real part first */ 
      
      InvSemiNaiveReduced(rdataptr,
			  bw,
			  m,
			  rinvfltres,
			  cos_data,
			  sin_values,
			  scratchpad);
      
      /* now do imaginary part */
      
      InvSemiNaiveReduced(idataptr,
			  bw,
			  m,
			  iminvfltres,
			  cos_data,
			  sin_values,
			  scratchpad);
      
      /* will store normal, then tranpose before doing inverse fft */
      memcpy(rfourdata+(m*size), rinvfltres, sizeof(double) * size);
      memcpy(ifourdata+(m*size), iminvfltres, sizeof(double) * size);
      
      /* move to next set of coeffs */
      
      rdataptr += bw-m;
      idataptr += bw-m;
      
  }
  
  /* now fill in zero values where m = bw (from problem definition) */
  memset(rfourdata + (bw * size), 0, sizeof(double) * size);
  memset(ifourdata + (bw * size), 0, sizeof(double) * size);
  
  
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
      /*** first set location of file pointer ***/
      fseek( fp , ptr_locs[size-m] , 0);
      /* now get how many elements in cos_array */
      fread ( &data_amount, sizeof(int), 1, fp );
      /* now read in that amount */
      fread (cos_data, sizeof(double), data_amount, fp );


	/* do real part first */
	
	InvSemiNaiveReduced(rdataptr,
			    bw,
			    size - m,
			    rinvfltres,
			    cos_data,
			    sin_values,
			    scratchpad);
	
	/* now do imaginary part */
	
	InvSemiNaiveReduced(idataptr,
			    bw,
			    size - m,
			    iminvfltres,
			    cos_data,
			    sin_values,
			    scratchpad);
	

	/* will store normal, then tranpose before doing inverse fft    */
	if ((m % 2) != 0)
	  for(i=0; i< size; i++){
	    rinvfltres[i] *= -1.0;
	    iminvfltres[i] *= -1.0;
	  }
	
	memcpy(rfourdata + (m*size), rinvfltres, sizeof(double) * size);
	memcpy(ifourdata + (m*size), iminvfltres, sizeof(double) * size);
	
	
	/* move to next set of coeffs */
	rdataptr += bw-(size-m);
	idataptr += bw-(size-m);
	
      } /* closes m loop */
  }
  else {

    for(m = bw + 1; m < size; m++){

      memcpy(rfourdata+(m*size), rfourdata+((size-m)*size),
	     sizeof(double) * size);
      memcpy(ifourdata+(m*size), ifourdata+((size-m)*size),
	     sizeof(double) * size);
      for(i = 0; i < size; i++)
	ifourdata[(m*size)+i] *= -1.0;
    }

  }

  /** now transpose **/
  transpose(rfourdata, size);
  transpose(ifourdata, size);
  
  
  /* now do inverse fourier grid computation */
  grid_invfourier(rfourdata, ifourdata, rdata, idata, size, scratchpad);
  
  free(ptr_locs);
  fclose(fp);

  /* amscray */
  
}
/************************************************************************
  Zonal Harmonic transform using hybrid algorithm - used in convolutions

  rdata and idata should be pointers to size x size arrays.
  rres and ires should be pointers to double arrays of size bw.
  size = 2 * bw

  workspace needed is (12 * bw) 

  FZT_hybrid_disk only computes spherical harmonics for m=0.
  
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
  
  */

void FZT_hybrid_disk( double *rdata, double *idata,
		      double *rres, double *ires,
		      int size, double *workspace,
		      int dataformat, int lim,
		      double **Z, struct lowhigh *poly,
		      int *loc, double *disk_array)
{
  FILE *fp;
  int bw, i, j;
  int data_amount, data_amount2;
  double *scratchpad, *cos_even, dsize;
  double *cos_data, *split_data;
  double tmpreal, tmpimag, *r0, *i0;
  char strm[128], tmpstr[128];
  int l, dummy ;
  double tmpA ;


  cos_data = disk_array;
  bw = size/2;
  dsize = (double) size;

  /* assign space */
  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (23 * bw), more than enough */

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Hybrid_bw%d.dat", bw);
  strcat(strm, tmpstr);
  fp = fopen(strm, "r");

  if ( fp == NULL )
    {
      fprintf( stderr,"Error in FZT_hybrid_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }

  /*** have to copy data to arrays ***/
  /* first get how many elements in cos_array */
  fread ( &data_amount, sizeof(int), 1, fp );
  /* now read in that amount */
  fread (cos_data, sizeof(double), data_amount, fp );
  
  /* assign the split_ptr now */
  split_data = cos_data + data_amount ;
  /* how many elements for the split array ? */
  fread ( &data_amount, sizeof(int), 1, fp );
  /* read in that many */
  fread ( split_data, sizeof(double) , data_amount, fp );
  
  /* ditto for the fct'd shifted legendres; note where I'm
     saving the shifted data ... FLT_HYBRID_FST
     expects this Legendres to immediately follow the
     splits */

  fread ( &data_amount2, sizeof(int), 1, fp );
  fread ( split_data + data_amount, sizeof(double),
	  data_amount2 , fp );

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
		  cos_data,
		  split_data);
  if(dataformat == 0) /* do the imaginary part */
    FLT_HYBRID_FST( 0, bw, i0, ires,
		    loc,
		    Z, scratchpad,
		    poly,
		    cos_data,
		    split_data);
  else                 /* otherwise set coefficients = 0 */
    memset(ires, 0, sizeof(double) * size);

  fclose(fp);

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

/************************************************************************

  multiplies harmonic coefficients of a function and a filter
  
  See convolution theorem of Driscoll and Healy
  
  datacoeffs should be output of FST_hybrid_disk,
  filtercoeffs the output of an FZT. 
  There should be (bw * bw) datacoeffs and bw filtercoeffs.
  rres and ires should point to arrays of dimension bw * bw.
  size parameter is 2*bw  
  
  *****************************/

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
  Uses hydrid and seminaive algorithms for spherical harmonic transforms
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
  disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.
  workspace needs to be a double pointer to an array of size
	16 * bw * bw + 33 * bw.


   ASSUMPTIONS:
   1. data is strictly REAL

*/

void Conv2Sphere_hyb_disk( double *rdata, double *idata,
			   double *rfilter, double *ifilter,
			   double *rres, double *ires,
			   int size, int lim,
			   int *loc, double **Z,
			   struct lowhigh *poly,
			   double *disk_array,
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


  FST_hybrid_disk( rdata, idata,
		   frres, fires,
		   size, scratchpad,
		   1, lim,
		   Z, poly,
		   loc, disk_array);
    

  FZT_hybrid_disk( rfilter, ifilter,
		   filtrres, filtires,
		   size, scratchpad, 1,
		   lim, Z,
		   poly, loc, disk_array) ;

  TransMult(frres, fires, filtrres, filtires, trres,
	    tires, size);

  InvFST_semi_disk(trres, tires,
		   rres, ires,
		   size, scratchpad,
		   1, disk_array);
  

}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

#else /* FFTPACK is defined */

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/




/*********************************************************************


  Performs a spherical harmonic transform using the hybrid and
  seminaive algorithms

  Size is the dimension of the input array (size x size) and it is
  expected that size=2*bw.  The inputs rdata and idata are expected
  to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
  to be pointers to bw x bw arrays, and will contain the harmonic
  coefficients in a "linearized" form.
   
  workspace needs to be a double pointer to an array of size
  8 * bw * bw + 33 * bw.
 
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

   input data assumed real

*/

void FST_hybrid_disk( double *rdata,
		      double *rcoeffs,
		      double *icoeffs,
		      int size,
		      int lim,
		      int *loc,
		      double **Z,
		      double *workspace,
		      double *wSave,
		      double *CoswSave,
		      double *CoswSave2,
		      double *disk_array)
{
  FILE *fp;
  int bw, m, i, j;
  double *scratchpad;
  int tmpindex[2];
  double *cos_even;
  double pow_one;
  long int *ptr_locs;
  int data_amount, data_amount2;
  double *cos_data, *split_data;
  int *loc_home, *loc_places;
  double *rres;
  char strm[128], tmpstr[128];
  double *tmprdata_ptr;
  int l, dummy ;
  double tmpA, tmpB ;

  cos_data = disk_array;
  bw = size/2;

  /* assign space */
  rres = workspace;   /* needs (size * size) = (4 * bw^2) */
  cos_even = rres + (size*size); /* needs bw */
  scratchpad = cos_even + bw; /* needs (24 * bw)  */
 
  /* total workspace is (29 * bw) */

  ptr_locs = (long int *) malloc(sizeof(long int) * bw);
  loc_places = (int *) malloc(sizeof(int) * (lim + 1) );
  loc_home = loc;

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Hybrid_bw%d.dat", bw);
  strcat(strm, tmpstr);

  fp = fopen(strm, "r");
  if ( fp == NULL )
    {
      fprintf( stderr,"Error in FST_hybrid_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }

  /* first save input data: don't want to write over it! */
  memcpy(rres, rdata, sizeof(double) * size * size);
  tmprdata_ptr = rres; /* to be used for advancing */

  /* do the FFTs along phi */
  grid_fourier_FFTPACK(rres, size, wSave);


  /* do the legendre transforms for orders m = 0
     through m = lim via the hybrid algorithm */
  
  for(m = 0 ; m < lim + 1 ; m ++) {

    /*** have to copy data to arrays ***/
    
    /*** first get location of file ptr, and locs ptr;
    I have to do this in case the input data is complex, it
    makes the book-keeping involved with doing the negative
    orders easier***/
    
    ptr_locs[m] = ftell( fp ) ;
    loc_places[m] = (loc - loc_home);

    /* first get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );

    /* assign the split_ptr now */
    split_data = cos_data + data_amount ;
    /* how many elements for the split array ? */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* read in that many */
    fread ( split_data, sizeof(double) , data_amount, fp );

    /* ditto for the fct'd shifted legendres; note where I'm
       saving the shifted data ... FLT_HYBRID_FST
       expects this Legendres to immediately follow the
       splits */
    fread ( &data_amount2, sizeof(int), 1, fp );
    fread ( split_data + data_amount, sizeof(double),
	    data_amount2 , fp );


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
		    cos_data,
		    split_data);
    
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
			cos_data,
			split_data);
	
	icoeffs += bw-m;
	tmprdata_ptr += size;
      }
    else
      {
	  memset(icoeffs, 0, sizeof(double) * bw);
	  icoeffs += bw;
      }


    /***
      advance to the next switch point, quant_locs and locs
      ***/
    loc += 2 + loc[1];

  }


  /** now do the Legendre transforms with seminaive
    algorithm for orders m = lim + 1 through bw - 1 **/

  for (m = lim + 1 ; m < bw ; m ++ ) {

    /*** have to copy data to arrays ***/

    /*** first get location of file ptr; I have to do
      this in case the input data is complex, it
      makes the book-keeping involved with doing the negative
      orders easier***/
   
    /*** first get location of file pointer ***/
    ptr_locs[m] = ftell( fp );

    /* first get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );

    /* do the real part */
    SemiNaiveReduced(tmprdata_ptr,
		     bw, 
		     m, 
		     rcoeffs, 
		     cos_data,
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
			 cos_data,
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
  
  /* now since the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m)),
     so use that to get the rest of the coefficients
     
     */

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


  fclose(fp);
  free(loc_places);
  free(ptr_locs);

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


/************************************************************************/
/* Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
   are harmonic coefficients stored in (bw * bw) arrays in the order
   spec'ed above.  rdata and idata are (size x size) arrays with
   the transformed result.  size is expected to be 2 * bw.

   workspace is 8 * bw * bw + 33 * bw

   dataformat = 1: input data is real, = 0: input data complex

   disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.

   

*/


void InvFST_semi_disk(double *rcoeffs,
		      double *icoeffs, 
		      double *rdata,
		      int size,
		      double *workspace,
		      double *wSave,
		      double *disk_array)
{
  FILE *fp;
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *scratchpad;
  double *sin_values, *eval_pts;
  long int *ptr_locs;
  double *cos_data;
  int data_amount;
  char strm[128], tmpstr[128];
  double *tmp_ptr;
  double tmpA, tmpB ;
  FILE *fp;

  bw = size/2;
  cos_data = disk_array;

  /* allocate space */

  ptr_locs = (long int *) malloc(sizeof(long int) * bw);
  sin_values = workspace;                  /* needs (2 * bw) */
  eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
  scratchpad = eval_pts + (2 * bw);       /* needs (16 * bw) */
  
  /* total workspace = 33 * bw */

  /* load up the sin_values array */
  n = 2*bw;
  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm,PRECOMP_DIR);
  sprintf(tmpstr,"InvSemi_bw%d.dat", bw);
  strcat(strm, tmpstr);

  fp = fopen(strm, "r");

  if ( fp == NULL )
    {
      fprintf( stderr,"Error in InvFST_semi_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }

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

    /*** have to copy data to arrays ***/
    
    /*** first get location of file pointer ***/
    ptr_locs[m] = ftell( fp );
    /* get how many elements in cos_array */
    fread ( &data_amount, sizeof(int), 1, fp );
    /* now read in that amount */
    fread (cos_data, sizeof(double), data_amount, fp );


      /* do real part first */ 
      
      InvSemiNaiveReduced(rdataptr,
			  bw,
			  m,
			  tmp_ptr,
			  cos_data,
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
			      cos_data,
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
  
  
  /* now since the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     invf-hat(l,-m) = conjugate(invf-hat(l,m)),
     so use that to get the rest of the real data
     
     dataformat =0 -> samples are complex, =1 -> samples real
     
     */

  /** now transpose **/
  transpose(rdata, size);
  
  /* now do inverse fourier grid computation */
  grid_invfourier_FFTPACK(rdata, size, wSave );
  
  free(ptr_locs);
  fclose(fp);

  /* amscray */
  
}
/************************************************************************
  Zonal Harmonic transform using hybrid algorithm - used in convolutions

  rdata and idata should be pointers to size x size arrays.
  rres and ires should be pointers to double arrays of size bw.
  size = 2 * bw

  workspace needed is (12 * bw) 

  FZT_hybrid_disk only computes spherical harmonics for m=0.
  
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
  
  */

void FZT_hybrid_disk( double *rdata,
		      double *rres,
		      double *ires,
		      int size,
		      int lim,
		      int *loc,
		      double **Z,
		      double *workspace,
		      double *CoswSave,
		      double *CoswSave2,
		      double *disk_array)
{
  FILE *fp;
  int bw, i, j;
  int data_amount, data_amount2;
  double *scratchpad, *cos_even, dsize;
  double *cos_data, *split_data;
  double tmpreal, tmpimag, *r0, *i0;
  char strm[128], tmpstr[128];
  int l, dummy ;
  double tmpA ;


  cos_data = disk_array;
  bw = size/2;
  dsize = (double) size;

  /* assign space */
  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (23 * bw), more than enough */

  /*********************************/
  /*                               */
  /* OPEN FILE OF PRECOMPUTED DATA */
  /*                               */
  /*********************************/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Hybrid_bw%d.dat", bw);
  strcat(strm, tmpstr);
  fp = fopen(strm, "r");

  if ( fp == NULL )
    {
      fprintf( stderr,"Error in FZT_hybrid_disk -->\n" );
      perror( strm );
      exit( 1 ) ;
    }

  /*** have to copy data to arrays ***/
  /* first get how many elements in cos_array */
  fread ( &data_amount, sizeof(int), 1, fp );
  /* now read in that amount */
  fread (cos_data, sizeof(double), data_amount, fp );
  
  /* assign the split_ptr now */
  split_data = cos_data + data_amount ;
  /* how many elements for the split array ? */
  fread ( &data_amount, sizeof(int), 1, fp );
  /* read in that many */
  fread ( split_data, sizeof(double) , data_amount, fp );
  
  /* ditto for the fct'd shifted legendres; note where I'm
     saving the shifted data ... FLT_HYBRID_FST
     expects this Legendres to immediately follow the
     splits */

  fread ( &data_amount2, sizeof(int), 1, fp );
  fread ( split_data + data_amount, sizeof(double),
	  data_amount2 , fp );

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
		  cos_data,
		  split_data);

  /* set imaginary part = to zero */
  memset(ires, 0, sizeof(double) * size);

  fclose(fp);

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

/************************************************************************

  multiplies harmonic coefficients of a function and a filter
  
  See convolution theorem of Driscoll and Healy
  
  datacoeffs should be output of FST_hybrid_disk,
  filtercoeffs the output of an FZT. 
  There should be (bw * bw) datacoeffs and bw filtercoeffs.
  rres and ires should point to arrays of dimension bw * bw.
  size parameter is 2*bw  
  
  *****************************/

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
  Uses hydrid and seminaive algorithms for spherical harmonic transforms
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
  disk_array = points to a double array of length
        TableSize(0,bw) + TableSize(1,bw); the data that is
	read off the disk will be stored here.
  workspace needs to be a double pointer to an array of size
	16 * bw * bw + 33 * bw.


   ASSUMPTIONS:
   1. data is strictly REAL

*/

void Conv2Sphere_hyb_disk( double *rdata,
			   double *rfilter,
			   double *rres,
			   int size,
			   int lim,
			   int *loc,
			   double **Z,
			   double *disk_array,
			   double *workspace,
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


  FST_hybrid_disk( rdata,
		   frres,
		   fires,
		   size, lim,
		   loc,
		   Z,
		   scratchpad,
		   wSave,
		   CoswSave,
		   CoswSave2,
		   disk_array);
    

  FZT_hybrid_disk( rfilter,
		   filtrres,
		   filtires,
		   size,
		   lim,
		   loc,
		   Z,
		   scratchpad,
		   CoswSave,
		   CoswSave2,
		   disk_array );


  TransMult(frres, fires, filtrres, filtires, trres,
	    tires, size);

  InvFST_semi_disk(trres,
		   tires,
		   rres,
		   size,
		   scratchpad,
		   wSave,
		   disk_array);
}

#endif /* FFTPACK */
