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


/****************************************************************

  Source code to convolve two functions defined on the 2-sphere.
  Uses seminaive algorithms.

  PRECOMPUTES IN MEMORY ALL NECESSARY DATA PRIOR TO
  TRANSFORMING.

  Reads in a function and filter from files specified at
  shell level, and dumps output into a
  specified file.  
  Both function and filter must be
  (size x size) arrays, where
  size = 2*bandwidth.  
  
  User specifies bandwidth as third argument
  
  Sample call:
  
  CONV_SEMI_MEMO sphere_bw64.dat filter_bw64.dat 64

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

 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "FST_semi_memo.h"
#include "MathFace.h"
#include "cospmls.h"

#ifdef FFTPACK
#include "newFCT.h"
#include "fft_grids.h"
#endif

void main(int argc, char **argv)
{
  FILE *signalfp, *filterfp, *resultfp;
  int bandwidth, bw, size, sizesquared, width, height;
  register int i;
  double *rsignal, *isignal, *rfilter, *ifilter, *rresult, *iresult;
  double *workspace;
  double *wSave, *CoswSave;

  if (argc != 4)
    {
      fprintf(stdout,"Usage: CONV_SEMI_MEMO signal_file filter_file bw\n");
      exit(0);
    }


  bandwidth = atoi(argv[3]);
  bw = bandwidth;
  size = 2*bandwidth;


#ifdef FFTPACK
  fprintf(stdout,"precomputing for fft, dct\n");

  wSave = precomp_fft( size );

  CoswSave = precomp_dct( size );

#endif
 

  rsignal = (double *) malloc(sizeof(double) * size * size);
  rfilter = (double *) malloc(sizeof(double) * size * size);
  rresult = (double *) malloc(sizeof(double) * size * size);

#ifndef FFTPACK
  isignal = (double *) malloc(sizeof(double) * size * size);
  ifilter = (double *) malloc(sizeof(double) * size * size);
  iresult = (double *) malloc(sizeof(double) * size * size);
#endif


  workspace = (double *) malloc(sizeof(double) 
				* ((2 * (Spharmonic_TableSize(bw)))  +
				   (8 * (bw*bw))  +
				   (2 * bw) +
				   (8 * (bw*bw)) + 
				   (32 * bw)));

#ifndef FFTPACK
  /****
    At this point, check to see if all the memory has been
    allocated. If it has not, there's no point in going further.
    ****/

  if ( (rsignal == NULL) || (isignal == NULL) ||
       (rfilter == NULL) || (ifilter == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (workspace == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#else
  /***
    At this point, check to see if all the memory has been
    allocated. If it has not, there's no point in going further.
    ****/

  if ( (rsignal == NULL) ||
       (rfilter == NULL) ||
       (rresult == NULL) ||
       (workspace == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#endif /* FFTPACK */


  signalfp = fopen(argv[1],"r");
  fprintf(stdout,"Reading signal file...\n");
  readMMRealTable(signalfp, rsignal, size*size, &width, &height);
  fclose(signalfp);
  sizesquared = size*size;


  filterfp = fopen(argv[2],"r");
  fprintf(stdout,"Reading filter file...\n");
  readMMRealTable(filterfp, rfilter, size*size, &width, &height);
  fclose(filterfp);

#ifndef FFTPACK
  for (i=0; i<sizesquared; i++)
    isignal[i] = 0.0;

  for (i=0; i<sizesquared; i++)
    ifilter[i] = 0.0;

  fprintf(stdout,"Calling Conv2Sphere_semi()\n");
  Conv2Sphere_semi_memo(rsignal, isignal, 
			rfilter, ifilter, 
			rresult, iresult, 
			size,
			workspace);
#else
  fprintf(stdout,"Calling Conv2Sphere_semi()\n");
  Conv2Sphere_semi_memo(rsignal,
			rfilter,
			rresult,
			size,
			workspace,
			wSave,
			CoswSave );
#endif

  resultfp = fopen("convres.dat","w");
  fprintf(stdout,"Writing output file...\n");
  /* result should be real-valued */
  seanprintMMRealTable(resultfp, rresult, size, size);

  fclose(resultfp);

#ifdef FFTPACK
  free( CoswSave );
  free( wSave );
#endif


}




