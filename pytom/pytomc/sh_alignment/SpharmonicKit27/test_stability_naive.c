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

/* test_stability_naive.c - top-level code for computing error data
   for the Legendre transform using the naive algorithm, i.e.,
   the three-term recurrence

   bw   - bandwidth
   m    - order of the problem
   trials - number of trials
   infile - filename for input data
   outfile - filename for result

   Sample Calls:

   test_stability_naive bw m trials infile outfile

   test_stability_naive 128 64 10 coeffs_flat.mca errordata.dat


*/

#include <math.h>
#include <stdio.h>

#include "naive_synthesis.h"

#define max(A, B) ((A) > (B) ? (A) : (B))
   
int main(argc, argv)
     int argc;
     char **argv;
{
  int i, j;
  int bw, m, trials;
  FILE *fpin, *fpout;
  double *coeffs, *computed_coeffs, *errors, *function;
  double *workspace;
  double relerror, curmax, granderror, grandrelerror;
  double mag, *moreerrors;

  if (argc != 6)
    {
      fprintf(stdout,"Usage: test_stability_naive bw m trials infile outfile.\n");
      return(0);
    }

  bw = atoi(argv[1]);
  m = atoi(argv[2]);
  trials = atoi(argv[3]);


  /* allocate memory */
  coeffs = (double *) malloc(sizeof(double) * (bw-m));
  computed_coeffs = (double *) malloc(sizeof(double) * (bw-m));
  errors = (double *) malloc(sizeof(double) * (bw-m));
  moreerrors = (double *) malloc(sizeof(double) * (bw-m));
  function = (double *) malloc(sizeof(double) * bw * 2);
  workspace = (double *) malloc(sizeof(double) * 32 * bw);

  /* open input file of coefficients */
  fpin = fopen(argv[4], "r");

  fprintf(stderr,"\n");

  /* open and set up output file */
  fpout = fopen(argv[5], "w");
  fprintf(stderr,"Naive Error Test \n");
  fprintf(stderr,"bw = %d\n", bw);
  fprintf(stderr,"m = %d\n", m);
  fprintf(stderr,"#trials = %d\n", trials);
  /*  fprintf(stderr,"--------\n"); */

  granderror = 0.0;
  grandrelerror = 0.0;


  /* main trial loop */
  for (i=0; i<trials; i++)
    {
      /* get a new set of coefficients */
      for (j=0; j<(bw-m); j++)
	fscanf(fpin, "%lf", coeffs+j);

      /* synthesize the corresponding function */
      Naive_Synthesize(bw, m, coeffs, function, workspace);

      /* Compute the transform */
      Naive_Analysis(bw, m, function, computed_coeffs, workspace);

      /*** compute absolute errors	***/
      for (j=0; j<(bw-m); j++)
	moreerrors[j] = fabs(computed_coeffs[j]-coeffs[j]);


      relerror = 0.0;
      curmax = 0.0;
      for(j = 0; j < (bw - m); j++)
	{
	  errors[0] = fabs(computed_coeffs[j]-coeffs[j]);
	  mag = fabs(coeffs[j]);
	  relerror = max(relerror,errors[0]/(mag + pow(10.0, -50.0)));
	  curmax = max(curmax,errors[0]);
	}

      granderror += curmax;
      grandrelerror += relerror;


      /*** output (abserror, coeff) pairs       ***/
      for (j=0; j<(bw-m); j++)
	fprintf(fpout,"%d\t %17.15e  %17.15e\n", j, moreerrors[j], coeffs[j]);
      fprintf(fpout,"\n");


    }

  fprintf(stderr,"Average r-o error:\t\t %.2e\n", granderror/((double) trials));
  fprintf(stderr,"Average (r-o)/o error:\t\t %.2e\n\n", grandrelerror/((double) trials));


  fclose(fpin);
  fclose(fpout);

  free(workspace); free(function); free(moreerrors);
  free(errors); free(computed_coeffs);
  free(coeffs);

  return(0);

}

