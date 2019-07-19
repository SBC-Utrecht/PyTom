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

/* semi.c - top-level code for computing Legendre transform using the
   seminaive algorithm 

   bw   - bandwidth
   m    - order of the problem
   k - number of desired coeffcients  1 <= k <= bw-m 
   infile - filename for input data
   outfile - filename for result
   timing - 1 print out timing data
*/

#include <stdio.h>
#include <stdlib.h>

#include "cospmls.h"
#include "seminaive.h"

main(argc, argv)
     int argc;
     char **argv;
{
  int i;
  int bw, m, k, timing;
  FILE *infile, *outfile;
  double *data, *result, *cos_pml_table, *workspace;

  if (argc < 6)
    {
      fprintf(stdout,"Usage: semi bw m k input_file output_file timing\n");
      return(0);
    }

  bw = atoi(argv[1]);
  m = atoi(argv[2]);
  k = atoi(argv[3]);

  if (argc == 7)
    timing = atoi(argv[6]);
  else
    timing = 0;

  data = (double *) malloc(sizeof(double) * 2 * bw);
  result = (double *) malloc(sizeof(double) * bw);
  workspace = (double *) malloc(sizeof(double) * 16 * bw);
  cos_pml_table = (double *) malloc(sizeof(double) * TableSize(m,bw));


  /* read in function */
  infile = fopen(argv[4],"r");
  for (i=0; i<(2*bw); i++)
    fscanf(infile,"%lf",data+i);
  fclose(infile);

  /* precompute cosine series for Pml (Gml) functions */
  CosPmlTableGen(bw, m, cos_pml_table, workspace);

  /* make sure k is a legal value */
  if ((k < 1) || (k > bw-m))
    k = bw-m;

  /* do it */
  SemiNaive(data, bw, m, k, result, cos_pml_table, timing, workspace);

  outfile = fopen(argv[5],"w");
  for (i=0; i<k; i++)
    fprintf(outfile,"%lf\n",result[i]);
  fclose(outfile);

  printf("Done\n");

  free(cos_pml_table);
  free(workspace);
  free(result);
  free(data);

}

