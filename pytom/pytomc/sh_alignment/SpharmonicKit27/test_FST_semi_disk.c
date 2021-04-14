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


/*************************************************************************/

/* 

Source code to test full spherical harmonic transform
using the seminaive algorithm coded up during October, 1995.


WILL READ PRECOMPUTED DATA OFF DISK !!!



Idea is to record the execution time for the spectral->grid->spectral
roundtrip, after precomputations have been done.

The strategy is to generate random coefficients representing
a complete spherical harmonic series of funcions Y(m,l).
The ordering of the these coefficients is assumed to be the
same as that of the output of the FST_seminaive() routine, namely


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
    

This means that there are (bw*bw) - (2*bw) coefficients.

Once the coefficients are generated, the corresponding
function is synthesized using InvFST_seminaive(), then
transformed (analyzed) using FST_seminaive(). Timing data
is printed.

Sample call

% test_FST_semi_disk bw loops [error_file]

NOTE: In this program, the coefficients generated are such that the grid points
(sample values) produced are REAL.

Appropriate timimg data will be printed out.

*/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "FST_semi_disk.h"
#include "cospmls.h"
#include "csecond.h"
#include "primitive_FST.h"

#ifdef FFTPACK
#include "newFCT.h"
#include "fft_grids.h"
#endif

#define max(A, B) ((A) > (B) ? (A) : (B))

/**************************************************************/
/**************************************************************/


void main(int argc, char **argv)
{
  FILE *fp, *errorsfp;
  int i, j, l, m, bw, size, loops, dummy;
  double tstart, tstop, total_time, for_time, inv_time;
  double *rcoeffs, *icoeffs, *rdata, *idata, *rresult, *iresult;
  double *workspace, *disk_array;
  double *relerror, *curmax, granderror, grandrelerror;
  double ave_error, ave_relerror;
  double stddev_error, stddev_relerror;
  double realtmp, imagtmp,origmag, tmpmag;
  double dumx, dumy, fudge;
  time_t seed;
  double *wSave, *CoswSave; /* for fftpack */
  double *rcoeffs2, *icoeffs2;

  granderror = 0.0; grandrelerror = 0.0; 


  if (argc < 3)
    {
      fprintf(stdout,"Usage: test_FST_semi_disk bw loops [error_file]\n");
      exit(0);
    }

  bw = atoi(argv[1]);
  loops = atoi(argv[2]);

  size = 2*bw;

  total_time = 0.0;
  for_time = 0.0;
  inv_time = 0.0;

#ifdef FFTPACK
  fprintf(stdout,"precomputing for fft, dct\n");

  wSave = precomp_fft( size );

  CoswSave = precomp_dct( size );

#endif
  

  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  rdata = (double *) malloc(sizeof(double) * (size * size));
#ifndef FFTPACK
  idata = (double *) malloc(sizeof(double) * (size * size));
#endif
  rresult = (double *) malloc(sizeof(double) * (bw * bw));
  iresult = (double *) malloc(sizeof(double) * (bw * bw));

#ifndef FFTPACK
  workspace = (double *) malloc(sizeof(double) * 
				((8 * (bw*bw)) + 
				 (33 * bw)));
#else
  workspace = (double *) malloc(sizeof(double) * 
				((4 * (bw*bw)) + 
				 (16 * bw)));
#endif


  /* should be PLENTY of space */
  disk_array = (double *) malloc(sizeof(double) *
				 2 * TableSize(0,bw)) ;



  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);


  /****
    At this point, check to see if all the memory has been
    allocated. If it has not, there's no point in going further.
    ****/
#ifndef FFTPACK
  if ( (rcoeffs == NULL) || (icoeffs == NULL) ||
       (rdata == NULL) || (idata == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (workspace == NULL) ||
       (disk_array == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#else
  if ( (rcoeffs == NULL) || (icoeffs == NULL) ||
       (rdata == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (workspace == NULL) ||
       (disk_array == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#endif



  /*** generate a seed, needed to generate random data ***/
  time(&seed);
  srand48(seed );

  /***
      Change for Version 2.6

      Since I want to calculate the error between the original
      coefficients and the calculated ones, I need space to copy the
      originals. Why? Because I'm renormalizing the coefficients
      inside the inverse transform. If I wasn't interested in the
      error, I wouldn't have to do this.
  ***/

  rcoeffs2 = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs2 = (double *) malloc(sizeof(double) * (bw * bw));

  fprintf(stdout,"about to enter loop\n\n");
  for(i=0; i<loops; i++){

    /**** loop to generate spherical harmonic coefficients *****/
    for(m=0;m<bw;m++)
      for(l=m;l<bw;l++){
	dumx = 2.0 * (drand48()-0.5);
	dumy = 2.0 * (drand48()-0.5);
	dummy = seanindex(m,l,bw);
	rcoeffs[dummy] = dumx;
	icoeffs[dummy] = dumy;
	dummy = seanindex(-m,l,bw);
	rcoeffs[dummy] = ((double) pow(-1.0, (double) m)) * dumx;
	icoeffs[dummy] = ((double) pow(-1.0, (double) (m + 1))) * dumy;
      }

    /* have to zero out the m=0 coefficients, since those are real */
    for(m=0;m<bw;m++)
      icoeffs[m] = 0.0;


    /*
      Change for Version 2.6:
      copy for the error checking later
    */
    memcpy( rcoeffs2, rcoeffs, sizeof(double)*bw*bw);
    memcpy( icoeffs2, icoeffs, sizeof(double)*bw*bw);


    /* do the inverse spherical transform */

    tstart = csecond();
#ifndef FFTPACK
    InvFST_semi_disk(rcoeffs2, icoeffs2, 
		     rdata, idata, 
		     size,
		     workspace,
		     1,
		     disk_array);
#else
    InvFST_semi_disk(rcoeffs2, icoeffs2, 
		     rdata,
		     size,
		     workspace,
		     wSave,
		     disk_array);
#endif
    tstop = csecond();
    inv_time += (tstop - tstart);
    
    fprintf(stdout,"inv time \t = %.4e\n", tstop - tstart);

    /* now do the forward spherical transform */

    tstart = csecond();
#ifndef FFTPACK
    FST_semi_disk(rdata, idata,
		  rresult, iresult,
		  size,
		  workspace,
		  1,
		  disk_array);
#else
    FST_semi_disk(rdata,
		  rresult, iresult,
		  size,
		  workspace,
		  wSave,
		  CoswSave,
		  disk_array);
#endif
    tstop = csecond();
    for_time += (tstop - tstart);
    
    fprintf(stdout,"forward time \t = %.4e\n", tstop - tstart);

    /* now to compute the error */
    relerror[i] = 0.0;
    curmax[i] = 0.0;
    for(j=0;j<(bw*bw);j++){
      realtmp = rresult[j]-rcoeffs[j];
      imagtmp = iresult[j]-icoeffs[j];
      origmag = sqrt((rcoeffs[j]*rcoeffs[j]) + (icoeffs[j]*icoeffs[j]));
      tmpmag  = sqrt((realtmp*realtmp) + (imagtmp*imagtmp));
      relerror[i] = max(relerror[i],tmpmag/(origmag + pow(10.0, -50.0)));
      curmax[i]  = max(curmax[i],tmpmag);
    }
    
    fprintf(stdout,"r-o error\t = %.12f\n", curmax[i]);
    fprintf(stdout,"(r-o)/o error\t = %.12f\n\n", relerror[i]);
    
    
    granderror += curmax[i];
    grandrelerror += relerror[i];
  
  }

  total_time = inv_time + for_time;

  ave_error = granderror / ( (double) loops );
  ave_relerror = grandrelerror / ( (double) loops );
  stddev_error = 0.0 ; stddev_relerror = 0.0;
  
  for( i = 0 ; i < loops ; i ++ )
    {
      stddev_error += pow( ave_error - curmax[i] , 2.0 );
      stddev_relerror += pow( ave_relerror - relerror[i] , 2.0 );
    }
  /*** this won't work if loops == 1 ***/
  if( loops != 1 )
    {
      stddev_error = sqrt(stddev_error / ( (double) (loops - 1) ) );
      stddev_relerror = sqrt(stddev_relerror / ( (double) (loops - 1) ) );
    }
  
  fprintf(stderr,"Program: test_FST_semi_disk\n");
  fprintf(stderr,"Bandwidth = %d\n", bw);
  
#ifndef WALLCLOCK
  fprintf(stderr,"Total elapsed cpu time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stderr,"Average cpu forward per iteration:\t %.4e seconds.\n",
	  for_time/((double) loops));  
  fprintf(stderr,"Average cpu inverse per iteration:\t %.4e seconds.\n",
	  inv_time/((double) loops));
#else
  fprintf(stderr,"Total elapsed wall time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stderr,"Average wall forward per iteration:\t %.4e seconds.\n",
	  for_time/((double) loops));  
  fprintf(stderr,"Average wall inverse per iteration:\t %.4e seconds.\n",
	  inv_time/((double) loops));
#endif

  fprintf(stderr,"Average r-o error:\t\t %.4e\t",
	  granderror/((double) loops));
  fprintf(stderr,"std dev: %.4e\n",stddev_error);
  fprintf(stderr,"Average (r-o)/o error:\t\t %.4e\t",
	  grandrelerror/((double) loops));
  fprintf(stderr,"std dev: %.4e\n\n",stddev_relerror);
  
  if (argc == 4)
    {
      fudge = -1.0;
      errorsfp = fopen(argv[3],"w");
      for(m = 0 ; m < bw ; m++ )
	{
	  fudge *= -1.0;
	  for(l = m ; l< bw ; l++ )
	    {
	      dummy = seanindex(m,l,bw);
	      
	      fprintf(errorsfp,
		      "dummy = %d\t m = %d\tl = %d\t%.10f  %.10f\n",
		      dummy, m, l,
		      fabs(rcoeffs[dummy] - rresult[dummy]),
		      fabs(icoeffs[dummy] - iresult[dummy]));
	      
	      dummy = seanindex(-m,l,bw);
	      
	      fprintf(errorsfp,
		      "dummy = %d\t m = %d\tl = %d\t%.10f  %.10f\n",
		      dummy, -m, l,
		      fabs(rcoeffs[dummy] - rresult[dummy]),
		      fabs(icoeffs[dummy] - iresult[dummy]));
	      
	    }
	}
      
      fclose(errorsfp);
      
    }


  free(icoeffs2); free(rcoeffs2);

  free(curmax); free(relerror);
  free(disk_array);
  free(workspace);
  free(iresult); free(rresult);
#ifndef FFTPACK
  free(idata);
#endif
  free(rdata);
  free(icoeffs);free(rcoeffs);

#ifdef FFTPACK
  free(CoswSave);
  free(wSave);
#endif
  

}


