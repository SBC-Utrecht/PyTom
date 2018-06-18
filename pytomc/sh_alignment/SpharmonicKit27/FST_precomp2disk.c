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


/****

  Source code to generate and save to disk precomputed data
  necessary for the seminaive forward and inverse spherical
  transforms and the hybrid forward spherical transform.

  input: bw = bandwidth of transform

  Sample call
   
  % FST_precomp2disk bw

  IMPORTANT IMPLEMENTATION NOTES:

  1. set the PRECOMP_SEMINAIVE and PRECOMP_HYBRID values
     below to the desired settings PRIOR TO compilation

  2. This routine saves files to disk. Check the setting
     PRECOMP_DIR in config.h to make sure that is the
     directory where the data will be saved (and where
     the data can be read). Currently, the files created
     will live in the same directory as the executable.
     Make changes before compiling.

 
  To understand the hybrid fst settings, read the documentation
  below and also that of SplitLocales() and HybridPts() in
  precomp_flt_hybrid.c.
  
  *****/

/***

  define to be 1 if wish to precompute and save
  to disk for that algorithm:

  PRECOMP_SEMINAIVE: precompute for seminaive fst (forward and reverse)
  PRECOMP_HYBRID: precompute for the hybrid fst (forward ONLY)

  ***/

#define PRECOMP_SEMINAIVE 1
#define PRECOMP_HYBRID 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "cospmls.h"
#include "precomp_flt_hybrid.h"

/***********************************************************************

  This function computes and writes to disk all of the
  (cosine transforms of) Legendre functions necessary
  to do a full spherical harmonic transform, i.e., it calls
  CosPmlTableGen for each value of m less than bw, and writes
  the result (the decimated cosine series of pml (even m)
  or gml (odd m) ) to disk.

  See CosPmlTableGen for further clarification.
  
  Inputs - the bandwidth bw of the problem

  resultspace,invspace - need to allocate
           2 * TableSize(0, bw) amount of space
	   for each array, to store results for writing to disk

  workspace - needs to be (16 * bw)

  WILL BE WRITING IN BINARY
  
*/

void Write_Seminaive_Tables(int bw,
			    double *resultspace,
			    double *invspace,
			    double *workspace)
{
  FILE *arrayfile_forward, *arrayfile_inverse;
  int i, result_size;
  char strm[128], strm2[128], tmpstr[64];

  /****************************************/
  /*                                      */
  /* OPEN FILES TO STORE PRECOMPUTED DATA */
  /*                                      */
  /****************************************/

  /**** where the forward seminaive fst data will be written ****/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Semi_bw%d.dat", bw);
  strcat(strm, tmpstr);

  /**** where the inverse seminaive fst data will be written ****/
  strcpy(strm2,PRECOMP_DIR);
  sprintf(tmpstr,"InvSemi_bw%d.dat", bw);
  strcat(strm2, tmpstr);

  arrayfile_forward = fopen(strm,"w");
  arrayfile_inverse = fopen(strm2,"w");

	    
  for ( i = 0 ; i < bw ; i ++ )
    {
      fprintf(stdout,"m = %d\n", i);
      result_size = TableSize(i, bw);

      /* to be on the safe side, let me initialize = 0 */
      memset(resultspace, 0, sizeof(double) * result_size);
      memset(invspace, 0, sizeof(double) * result_size);

      /* calculate the cosine transforms */
      CosPmlTableGen(bw, i, resultspace, workspace);

      /* first save size of the array to disk */
      fwrite(&result_size, sizeof(int), 1, arrayfile_forward);

      /* now save the array */
      fwrite(resultspace, sizeof(double), result_size, arrayfile_forward);

      /* now transpose it */
      Transpose_CosPmlTableGen(bw, i, resultspace, invspace);

      /* first save size of the array to disk */
      fwrite(&result_size, sizeof(int), 1, arrayfile_inverse);

      /* now save the array */
      fwrite(invspace, sizeof(double), result_size, arrayfile_inverse);

      /* that's it */
    }


  fclose(arrayfile_inverse);
  fclose(arrayfile_forward);

}

/**************************************************************/
/****

  This function computes and writes to disk the precomputed
  data required by the hybrid spherical transform.

  bw = bandwidth

  seminaive_tablespace = double array used as a temp space when
        precomputing for the seminaive algorithm, of size
	2 * TableSize(0,bw) (almost twice as much as needed,
	but why allocated less and worry about segmentation
	faults ?)

  workspace = double array of size 16 * bw;

  The following arguments were set by a call to HybridPts()
  ( which is defined in precomp_flt_hybrid.c )

  lim = order where switching from hybrid to seminaive algorithm
        when doing the forward spherical transform

  cos_size = int array containing the sizes of the data arrays
        to save when precomputing for the semi-naive portion
	of the hybrid algorithm
  
  cos_array = double array used as a tmp space when precomputing
        for the seminaive portion of the hybrid algorithm

  split_array = double array used as a tmp space when precomputing
        for the simple-split portion of the hybrid algorithm

  shift_size = int array containing the sizes of the data arrays
        to save when precomputing for the simple-split portion
	of the hybrid algorithm
  
  split_array = double array used as a tmp space when precomputing
        for the simple-split portion of the hybrid algorithm

  WILL BE WRITING IN BINARY

  *********************************************/



void Write_Hybrid_Tables( int bw, int lim,
			  int *loc, double *seminaive_tablespace,
			  int *cos_size, double *cos_array,
			  double *split_array,
			  int *shift_size, double *shift_array,
			  double *workspace )
{
  FILE *arrayfile;
  char strm[128], tmpstr[128];
  int i, tmp_int;
    
  fprintf(stdout,"\nPrecomputing for the hybrid...\n");
  fprintf(stdout,"\nFirst the simple-split portion\n");
  
  
  /****************************************/
  /*                                      */
  /* OPEN FILES TO STORE PRECOMPUTED DATA */
  /*                                      */
  /****************************************/

  /**** where the forward hybrid fst data will be written ****/
  strcpy(strm, PRECOMP_DIR);
  sprintf(tmpstr,"Hybrid_bw%d.dat", bw);
  strcat(strm, tmpstr);

  arrayfile = fopen(strm,"w");
  
  for (i = 0 ; i <= lim ; i ++)
    {
      
      fprintf(stdout,"m = %d\n", i);
      
      /* to be on the safe side, let me initialize = 0 */
      memset(cos_array, 0, sizeof(double) * cos_size[i]);
      
      /* now calculate array */
      CosPmlTableGenLim(bw, i, loc[0] + 1, cos_array, workspace);
      
      /* first save size of the array to disk */
      fwrite(&cos_size[i], sizeof(int), 1, arrayfile);
      
      /* now save the array */
      fwrite(cos_array, sizeof(double), cos_size[i], arrayfile);
      
      /*** NEXT ***/
      
      tmp_int = 2 * bw * loc[1];
      
      /* to be on the safe side, let me initialize = 0 */
      memset(split_array, 0, sizeof(double) * tmp_int);
      
      /* now calculate array */
      SplitPml(bw, i, loc[1], loc + 2,
	       split_array,
	       workspace);
      
      /* first save size of the array to disk */
      fwrite(&tmp_int, sizeof(int), 1, arrayfile);
      
      /* now save the array */
      fwrite(split_array, sizeof(double), tmp_int, arrayfile);
      
      
      /*** NEXT ***/
      
      /* to be on the safe side, let me initialize = 0 */
      memset(shift_array, 0, sizeof(double) * shift_size[i]);
      
      /* now calculate array */
      ShiftSemiY(bw, i, loc[0], loc+2,
		 loc[1],
		 shift_array,
		 workspace);
      
      /* first save size of the array to disk */
      fwrite(&shift_size[i], sizeof(int), 1, arrayfile);
      
      /* now save the array */
      fwrite(shift_array, sizeof(double), shift_size[i], arrayfile);
      
      loc += 2 + loc[1];
      
    }
  
  /* now do the semi-naive part */
  
  fprintf(stdout,"Now the semi-naive portion\n");
  
  for ( i = lim + 1 ; i < bw ; i ++ )
    {
      fprintf(stdout,"m = %d\n", i);
      tmp_int = TableSize(i, bw);
      
      /* to be on the safe side, let me initialize = 0 */
      memset(seminaive_tablespace, 0, sizeof(double) * tmp_int);
      
      /* calculate the cosine transforms */
      CosPmlTableGen(bw, i, seminaive_tablespace, workspace);
      
      /* first save size of the array to disk */
      fwrite(&tmp_int, sizeof(int), 1, arrayfile);
      
      /* now save the array */
      fwrite(seminaive_tablespace, sizeof(double), tmp_int, arrayfile);
      
      /* that's it */
    }
  
  fclose(arrayfile);
}




/**************************************************************/
/**************************************************************/


void main(int argc, char **argv)
{
  int bw;
  int max_cos_size, max_shift_size;
  int lim, max_split;
  int *cos_size, *shift_size, *loc;
  double *seminaive_tablespace, *trans_seminaive_tablespace, *workspace;
  double *cos_array, *shift_array, *split_array;


  if (argc != 2)
    {
      fprintf(stdout,"Usage: FST_precomp2disk bw\n");
      exit(0);
    }



  bw = atoi(argv[1]);
 
  /***
    allocating memory to hold precomputed data
    ***/

  seminaive_tablespace =
    (double *) malloc(sizeof(double) * 2 * TableSize(0, bw) );

  trans_seminaive_tablespace =
    (double *) malloc(sizeof(double) * 2 * TableSize(0, bw) );


  workspace = (double *) malloc(sizeof(double) *  16 * bw );


  /*****

    First precompute for the seminaive algorithm

    *****/


  if( PRECOMP_SEMINAIVE )
    {
      fprintf(stdout,"\nGenerating seminaive tables...\n");
      Write_Seminaive_Tables(bw,
			     seminaive_tablespace,
			     trans_seminaive_tablespace,
			     workspace);
    }

  /*****

    Now precompute for the hybrid algorithm

    *****/

  if( PRECOMP_HYBRID )
    {
      /*** allocate memory for sizes of arrays at the different
	orders; no matter what lim is, lim <= bw ***/

      cos_size = (int *) malloc(sizeof(int) * bw);
      shift_size = (int *) malloc(sizeof(int) * bw);
      
      loc = HybridPts( bw, &lim, cos_size, shift_size,
		       &max_split, 
		       &max_cos_size,
		       &max_shift_size );
            
      /*** allocate maximum memory needed ***/
      cos_array = (double *) malloc(sizeof(double) * max_cos_size);
      shift_array = (double *) malloc(sizeof(double) * max_shift_size);
      split_array = (double *) malloc(sizeof(double) * 2 * bw * max_split);
      
      Write_Hybrid_Tables( bw, lim,
			   loc, seminaive_tablespace,
			   cos_size, cos_array,
			   split_array,
			   shift_size, shift_array,
			   workspace ) ;
      
      free(split_array); free(shift_array); free(cos_array);
      free(shift_size); free(cos_size);
      free(loc);
    }      
  
  free(workspace);
  free(trans_seminaive_tablespace); free(seminaive_tablespace);

}



