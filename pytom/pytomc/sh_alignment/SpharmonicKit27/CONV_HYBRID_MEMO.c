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


  PRECOMPUTES IN MEMORY ALL NECESSARY DATA PRIOR TO
  TRANSFORMING.

  See important information below and in SplitLocales()
  in precomp_flt_hybrid.c, relating to use of the hybrid
  algorithm. There are lots of switches to set.


  Reads in a function and filter from files specified at
  shell level, and dumps output into a
  specified file.  
  Both function and filter must be
  (size x size) arrays, where
  size = 2*bandwidth.  
  
  User specifies bandwidth as THIRD argument
  
  Sample call:
  
  CONV_HYBRID_MEMO sphere_bw64.dat filter_bw64.dat 64

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

#include "FST_hybrid_memo.h"
#include "MathFace.h"
#include "cospmls.h"
#include "precomp_flt_hybrid.h"

#ifdef FFTPACK
#include "fft_grids.h"
#include "newFCT.h"
#endif

/**************************************************************/
/**************************************************************/


void main(int argc, char **argv)
{
  FILE *signalfp, *filterfp, *resultfp;
 
  int i, bw, size, lim, max_split, width, height;
  int max_cos_size, max_shift_size;
  int *loc, *cos_size, *shift_size;
  double *rfilter, *ifilter, *rsignal, *isignal, *rresult, *iresult;
  double *seminaive_tablespace, *trans_seminaive_tablespace, *workspace;
  double **seminaive_table, **trans_seminaive_table;
  double **split_dptr, **cos_dptr, **shift_dptr;
  double **Z, *zspace, *Zptr;
  struct lowhigh *poly;
  double *wSave, *CoswSave, *CoswSave2;

  if (argc != 4)
    {
      fprintf(stdout,"Usage: CONV_HYBRID_MEMO signal_file filter_file bw\n");
      exit(0);
    }


 
  bw = atoi(argv[3]);

  /* bw has to be at least 64 */
  if (bw < 64)
    {
      fprintf(stdout,"Error in main: bw must be at least 64\n");
      exit ( 0 );
    }


  size = 2*bw;
 


#ifdef FFTPACK
  fprintf(stdout,"precomputing for fft, dct\n");

  wSave = precomp_fft( size );

  CoswSave = precomp_dct( size );

  CoswSave2 = precomp_dct( bw );

#endif

  /***
    let's malloc stuff
    ***/

  rsignal = (double *) malloc(sizeof(double) * size * size);
  rfilter = (double *) malloc(sizeof(double) * size * size);
  rresult = (double *) malloc(sizeof(double) * size * size);

#ifndef FFTPACK
  isignal = (double *) malloc(sizeof(double) * size * size);
  ifilter = (double *) malloc(sizeof(double) * size * size);
  iresult = (double *) malloc(sizeof(double) * size * size);
#endif

  seminaive_tablespace =
    (double *) malloc(sizeof(double) * Spharmonic_TableSize(bw));
  
  trans_seminaive_tablespace =
    (double *) malloc(sizeof(double) * Spharmonic_TableSize(bw));

  workspace = (double *) malloc(sizeof(double) * 
				(16 * bw * bw + 33 * bw));

  /** needed for newFCTX as used in FLT_SIMSPL **/
  poly = (struct lowhigh *) malloc(sizeof(struct lowhigh) * 4 * bw );


  /****
    At this point, check to see if all the memory has been
    allocated. If it has not, there's no point in going further.
    ****/
#ifndef FFTPACK
  if ( (rsignal == NULL) || (isignal == NULL) ||
       (rfilter == NULL) || (ifilter == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (seminaive_tablespace == NULL) ||
       (trans_seminaive_tablespace == NULL) ||
       (workspace == NULL) || (poly == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#else
  if ( (rsignal == NULL) ||
       (rfilter == NULL) ||
       (rresult == NULL) ||
       (seminaive_tablespace == NULL) ||
       (trans_seminaive_tablespace == NULL) ||
       (workspace == NULL) || (poly == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }
#endif

  /*** precompute for the seminaive portion of the forward
    spherical transform algorithm (i.e. for those
    orders greater than lim) ***/

  fprintf(stdout,"Generating seminaive tables...\n");
  seminaive_table =
    Spharmonic_Pml_Table( bw, seminaive_tablespace, workspace );
					  

  /*** transpose the above precomputed data...this will
    be the precomputed data for the inverse spherical transform ***/

  fprintf(stdout,"Generating trans_seminaive tables...\n");
  trans_seminaive_table = 
    Transpose_Spharmonic_Pml_Table(seminaive_table,
				   bw,
				   trans_seminaive_tablespace,
				   workspace) ;


  /***
    now to precompute for the hybrid algorithm
    ***/


  /***
    allocate some space needed by HybridPts and
    get lim, locations of switch points, etc
    ***/

  cos_size = (int *) malloc(sizeof(int) * bw );
  shift_size = (int *) malloc(sizeof(int) * bw );

  loc = HybridPts( bw, &lim,
		   cos_size, shift_size,
		   &max_split,
		   &max_cos_size,
		   &max_shift_size ) ;


  /* This is space for the intermediate results
     of the hybrid algorithm */  
  zspace = (double *) malloc(sizeof(double) *
			     (2 * (max_split) * bw));
  Z = (double **) malloc(sizeof(double *) * (max_split));
  
  Z[0] = zspace;
  Zptr = zspace;
  for(i = 1 ; i < max_split + 0; i++)
    {
      Zptr += 2 * bw;
      Z[i] = Zptr;
    }
  

  /**
    I'm pointing it to the seminaive_tablespace.

    THIS IS OK PROVIDED THAT BW >= 64.

    At smaller bandwidths, use pure semi-naive.


    I first needed to precompute the seminaive data
    in order to take the transpose (to derive the precomputed
    data for the inverse spherical transform). Having done that,
    I don't care about the precomputed data at the start of this
    array. I can store the precomputed data for the hybrid algorithm
    here. Why? Because, basically (and with quite reasonable
    assumptions ... I think)

    sizeof(precomputed data for hybrid algorithm for orders m = 0 though lim)

    is less than

    sizeof(precomputed data for seminaive algorithm for orders m = 0 though lim)

    I could allocate separate space, but already at bw = 512 I'm
    allocating ALOT of memory. Might as well try to cut it down
    where I can.

    **/
  
  /*
    (double **)  ptrs which will be pointed to the various
    locations of the hybrid algorithm's precomputed data */

  split_dptr = (double **) malloc(sizeof(double *) * 3 * (lim + 1));
  cos_dptr = split_dptr + lim + 1;
  shift_dptr = cos_dptr + lim + 1;					    

  Hybrid_SetDptrs( bw,
		   lim,
		   loc, seminaive_tablespace,
		   cos_size, shift_size,
		   cos_dptr,
		   split_dptr,
		   shift_dptr ) ;




  /****

    finally can precompute the data !

    *****/

  HybridPrecomp( bw, lim,
		 loc,
		 cos_dptr,
		 split_dptr,
		 shift_dptr,
		 workspace );



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
  memset(isignal, 0, sizeof(double) * size * size);
  memset(ifilter, 0, sizeof(double) * size * size);
#endif

  fprintf(stdout,"Calling Conv2Sphere_hyb_memo()\n");

#ifndef FFTPACK
  Conv2Sphere_hyb_memo( rsignal, isignal,
			rfilter, ifilter,
			rresult, iresult,
			size, lim,
			loc, Z,
			poly,
			seminaive_table,
			trans_seminaive_table,
			split_dptr,
			cos_dptr,
			workspace) ;
#else
  Conv2Sphere_hyb_memo( rsignal,
			rfilter,
			rresult,
			size, lim,
			loc, Z,
			workspace,
			seminaive_table,
			trans_seminaive_table,
			split_dptr,
			cos_dptr,
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


