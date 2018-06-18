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

/* test_naive.c - top-level code for computing Legendre transform using the
   naive algorithm 

   m    - order of the problem
   bw   - bandwidth
   timing - 1 print out timing data
   loops - number of loops thru timed portion of code.  Intended
           to reduce noise due to multiprocessing and 
	   discretization errors
   infile - filename for input data
   outfile - filename for result
   
   Sample calls:

   test_naive m bw timing loops infile outfile

   test_naive 0 32 1 32 coeffs_flat.mca outdata.dat

*/

#include <stdio.h>
#include <stdlib.h>
#include "naive_synthesis.h"


int main(argc, argv)
     int argc;
     char **argv;
{
  int i;
  int bw, m, timing, loops;
  FILE *infile, *outfile;
  double *data, *result, *workspace;
  double runtime;

  if (argc < 7)
    {
      fprintf(stdout,"Usage: test_naive m bw timing loops input_file output_file\n");
      return(0);
    }

  m = atoi(argv[1]);
  bw = atoi(argv[2]);
  timing = atoi(argv[3]);
  loops = atoi(argv[4]);

  data = (double *) malloc(sizeof(double) * 2 * bw);
  result = (double *) malloc(sizeof(double) * bw);
  workspace = (double *) malloc(sizeof(double) * 18 * bw);

  /* read in function */
  infile = fopen(argv[5],"r");
  for (i=0; i<(2*bw); i++)
    fscanf(infile,"%lf",data+i);
  fclose(infile);

  /* let's initialize everything to 0 */
  for(i=0;i<bw;i++)
    result[i] = 0.0;

  /* do it */
  /*** time the grinding out of the Legendre functions
  Naive_Analysis_Timing(data, bw, m, result,
	    timing, &runtime, loops, workspace);
  *****/

  /*** precompute the Legendre functions, then turn
    on the stopwatch ***/
  Naive_Analysis_TimingX(data, bw, m, result,
	    timing, &runtime, loops, workspace);


  outfile = fopen(argv[6],"w");
  for (i=0; i<bw; i++)
    fprintf(outfile,"%f\n",result[i]);
  fclose(outfile);

  printf("Done\n");

  free(workspace); free(result); free(data);

  return(0);
}










