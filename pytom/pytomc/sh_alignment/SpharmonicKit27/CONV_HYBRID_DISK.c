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


/*****************************************************

  Source code to convolve two functions defined on the 2-sphere.
  Uses hybrid/seminaive algorithm in forward transform direction,
  and seminaive in the inverse direction.

  
  ASSUMES ALL PRECOMPUTED DATA IS FOUND ON THE DISK. NONE
  IS GENERATED IN THIS PROGRAM. 
  

  Reads in a function and filter from files specified at
  shell level, and dumps output into a
  specified file.  
  Both function and filter must be
  (size x size) arrays, where
  size = 2*bandwidth.  
  
  User specifies bandwidth as THIRD argument
  
  Sample call:
  
  CONV_HYBRID_DISK sphere_bw64.dat filter_bw64.dat 64

  Output is dumped to a file called convres.dat in Mathematica
  format.
  
  In this example, the sphere function and the filter are stored
  on disk files in Mathematica arrays, and are real-valued, not
  complex-valued.  Some sample files and filters are
  
  sphere_bw64.dat
  sphere_bw128.dat
  
  filter_bw64.dat
  filter_bw128.dat

  sphere_bw* = sample values of a bump on a "noisey" sphere
  filete_bw* = sample values of the filter: a smooth, symmetric
               bump at the north pole

  convres.dat = location of maximum value tells me where the
               bump is

  ********************************************************/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FST_hybrid_disk.h"
#include "MathFace.h"
#include "cospmls.h"
#include "precomp_flt_hybrid.h"

#ifdef FFTPACK
#include "newFCT.h"
#include "fft_grids.h"
#endif

/**************************************************************/
/**************************************************************/

void main(int argc, char **argv)
{
  FILE *signalfp, *filterfp, *resultfp;
  int i, bw, size;
  int lim, max_split, max_cos_size, max_shift_size;
  int *loc, *cos_size, *shift_size,  width, height;
  double *rfilter, *ifilter, *rsignal, *isignal, *rresult, *iresult;
  double **Z, *zspace, *Zptr, *disk_array, *workspace;
  struct lowhigh *poly;
  double *wSave, *CoswSave, *CoswSave2;


  if (argc != 4)
    {
      fprintf(stdout,"Usage: CONV_HYBRID_DISK signal_file filter_file bw\n");
      exit(0);
    }

 
  bw = atoi(argv[3]);

  size = 2*bw;

#ifdef FFTPACK
  fprintf(stdout,"precomputing for fft, dct\n");

  wSave = precomp_fft( size );

  CoswSave = precomp_dct( size );

  CoswSave2 = precomp_dct( bw );

#endif



  /***
    allocate memory for data, coefficients, and results
    ***/
  rsignal = (double *) malloc(sizeof(double) * size * size);
  rfilter = (double *) malloc(sizeof(double) * size * size);
  rresult = (double *) malloc(sizeof(double) * size * size);

#ifndef FFTPACK
  isignal = (double *) malloc(sizeof(double) * size * size);
  ifilter = (double *) malloc(sizeof(double) * size * size);
  iresult = (double *) malloc(sizeof(double) * size * size);
#endif

  /** needed for newFCTX **/
  poly = (struct lowhigh *) malloc(sizeof(struct lowhigh) * 4 * bw );

  /*** temporary workspace ***/
  workspace = (double *) malloc(sizeof(double) * 
				( 16 * bw * bw  + 33 * bw));

  /*** when precomputed data is read off the disk, it will
    be saved in this array; at any given order, the most
    precomputed needed is, at most, slightly greater than
    TableSize(0,bw); let me allocate roughly twice as much
    to avoid segmentation faults ***/

  disk_array = (double *) malloc( sizeof(double) * 2 *
				  TableSize(0,bw) ) ;
    
    
    
  /****
    At this point, check to see if all the memory has been
    allocated. If it has not, there's no point in going further.
    ****/
#ifndef FFTPACK  
  if ( (rsignal == NULL) || (isignal == NULL) ||
       (rfilter == NULL) || (ifilter == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (workspace == NULL) ||
       (disk_array == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#else
  if ( (rsignal == NULL) ||
       (rfilter == NULL) ||
       (rresult == NULL) ||
       (workspace == NULL) ||
       (disk_array == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#endif

  /*** allocate memory for sizes of arrays at the different
    orders; no matter what lim is, lim <= bw; I really don't
    need these arrays in this particular case, since I'm reading
    data off the disk; I need them because I use HybridPts in
    the function which precomputes and writes to disk the
    hybrid and seminaive data - then I DO need these int
    arrays; for similar reasons I really don't need max_cos_size
    or max_shift_size here, either ***/

  cos_size = (int *) malloc(sizeof(int) * bw);
  shift_size = (int *) malloc(sizeof(int) * bw);

  /**
    loc will contain the switch points, the quantity of splits,
    and the locations of the splits, for all the orders that
    the hybrid algorithm is used **/

  loc = HybridPts( bw, &lim,
		   cos_size, shift_size,
		   &max_split,
		   &max_cos_size,
		   &max_shift_size ) ;



  /* This is space for the intermediate results of the Legendre
     transform algorithm */
  zspace = (double *) malloc(sizeof(double) *
			     (2 * (max_split + 1) * bw));
  
  /*** to be used by FLT_HYBRID_FST ***/
  Z = (double **) malloc(sizeof(double *) * (max_split + 1));
  Z[0] = zspace;
  Zptr = zspace;
  for(i = 1 ; i < max_split + 1; i++)
    {
      Zptr += 2 * bw;
      Z[i] = Zptr;
    }


  /*****
    read in samples and filters
    ****/

  signalfp = fopen(argv[1],"r");
  fprintf(stdout,"Reading signal file...\n");
  readMMRealTable(signalfp, rsignal, size*size, &width, &height);
  fclose(signalfp);

  filterfp = fopen(argv[2],"r");
  fprintf(stdout,"Reading filter file...\n");
  readMMRealTable(filterfp, rfilter, size*size, &width, &height);
  fclose(filterfp);

#ifndef FFTPACK
  /* set imaginary coefficients = 0 */
  memset(isignal, 0, sizeof(double) * size * size );
  memset(ifilter, 0, sizeof(double) * size * size);
#endif

  fprintf(stdout,"Calling Conv2Sphere_hyb_disk()\n");

#ifndef FFTPACK
  Conv2Sphere_hyb_disk( rsignal, isignal,
			rfilter, ifilter,
			rresult, iresult,
			size, lim,
			loc, Z,
			poly,
			disk_array,
			workspace );
#else
  Conv2Sphere_hyb_disk( rsignal,
			rfilter,
			rresult,
			size,
			lim,
			loc,
			Z,
			disk_array,
			workspace,
			wSave,
			CoswSave,
			CoswSave2 );
#endif  

  resultfp = fopen("convres.dat","w");
  fprintf(stdout,"Writing output file...\n");
  /* result should be real-valued */
  seanprintMMRealTable(resultfp, rresult, size, size);

  fclose(resultfp);

#ifdef FFTPACK
  free( CoswSave2 );
  free( CoswSave );
  free( wSave );
#endif



}

  

