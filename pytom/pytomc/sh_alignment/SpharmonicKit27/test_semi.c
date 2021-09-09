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


/*
  test_semi.c - top-level code for computing Legendre transform using the
  seminaive algorithm 

  m    - order of the problem
  bw   - bandwidth
  k - number of desired coeffcients  1 <= k <= bw-m 

  loops - number of loops thru timed portion of code.  Intended
          to reduce noise due to multiprocessing and 
          discretization errors
  errflag - = 1: test for error (reads in more of norm.dat). If
                 your testing for error, loops shouldn't exceed
		 10. Depending on the bw and order, if loops > 10
		 then the routine can't allocate enough memory
		 (given how the routine computes the error)
            = 0: don't test for error (reads in just enough
	         from norm.dat); can do as many loops as
		 you like
		
                 
*/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* to declare memcpy */

#include "cospmls.h"
#include "primitive.h"
#include "seminaive.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define DATAFILE "norm.dat"   /* file which contains random numbers */

#ifdef FFTPACK
#include "newFCT.h"
#include "fftpack.h"
#endif

/**************************************************/
/**************************************************/

void main( int argc, char **argv)
{
  int i, j;
  int n, bw, m, k, loops;
  FILE *infile;
  double *data, *result, *cos_pml_table, *workspace;
  double *coeffs, *eval_args;
  double runtime;
  
  double tmp_error, max_error, ave_error, sum_error;
  double tmp_relerror, max_relerror, ave_relerror, sum_relerror;
  double stddev_error, stddev_relerror;
  int *bad_coeff, sum_bad_coeff, errorflag;
  double *error, *relerror, ave_bad_coeff;

  double *wSave;

  if (argc < 6)
    {
      fprintf(stdout,"Usage: test_semi m bw k loops errflag\n");
      exit(0);
    }

  m = atoi(argv[1]);
  bw = atoi(argv[2]);
  k = atoi(argv[3]);
  loops = atoi(argv[4]);
  errorflag = atoi(argv[5]);

#ifdef FFTPACK
  n = 2 * bw;

  wSave = precomp_dct( n );

#endif



  workspace = (double *) malloc(sizeof(double) * 16 * bw);
  cos_pml_table = (double *) malloc(sizeof(double) * TableSize(m,bw));

  /* space for error stuff */
  error = (double *) malloc(sizeof(double) * loops);
  relerror = (double *) malloc(sizeof(double) * loops);
  bad_coeff = (int *) malloc(sizeof(int) * loops);
  
  /* need at most bw coefficients */
  coeffs = (double *) malloc(sizeof(double) * (bw - m) *
			     (loops*errorflag + 1));
  
  /* space for result of Legendre transform */
  result = (double *) malloc(sizeof(double) * (bw - m) *
			     (loops*errorflag + 1));
  
  /* space for data */
  data = (double *) malloc(sizeof(double) * 2 * bw *
			   (loops*errorflag + 1));
  
  
  /* used to synthesize a test function using P_eval */
  eval_args = (double *) malloc(sizeof(double) * 2 * bw);

  /* generate arguments (arccos of sample points) for synthesizing
     a test function */  
  ArcCosEvalPts(2*bw,eval_args);
  
  /**********************************************/
  /* Here I'm just setting all the coefficients */
  /**********************************************/
  /*
    coeffs[0] = 0.0;
    for (i=1; i<bw; i++)
    coeffs[i] = (double) i;
    */
  infile = fopen(DATAFILE,"r");
  if ( infile == NULL )
    {
      fprintf( stderr,"Error in test_flt_classic -->\n" );
      perror( DATAFILE );
      exit( 1 ) ;
    }
  for (i=0; i<(bw-m)*(loops*errorflag+1); i++)
    fscanf(infile,"%lf",coeffs+i);
  fclose(infile);

  /* synthesize the corresponding function */
  for(j = 0 ; j < loops*errorflag+1 ; j++)
    P_eval(m, coeffs+j*(bw-m), eval_args,
	   data+j*2*bw, workspace, bw);

  /* precompute cosine series for Pml (Gml) functions */
  CosPmlTableGen(bw, m, cos_pml_table, workspace);

  /* make sure k is a legal value */
  if ((k < 1) || (k > bw-m))
    k = bw-m;

  /* do it */
  SemiNaiveX(data, bw, m, k, result, cos_pml_table, 
	     &runtime, loops, workspace,
	     errorflag);

  /* now to compute the error     */
  /*** have to reset loops if I'm looping over the same data
    (which I am when timing ***/

  if(errorflag == 0)
    loops = 1;
  
  for(i = 0 ; i < loops ; i ++)
    {
      max_error = 0.0 ; max_relerror = 0.0;
      for(j = 0 ; j < (bw - m) ; j ++)
	{
	  tmp_error = fabs(coeffs[i*(bw-m)+j] - result[i*(bw-m)+j]);
	  tmp_relerror = tmp_error /(fabs(coeffs[i*(bw-m)+j])+
				     pow(10.0, -50.0));
	  max_error = mymax(max_error, tmp_error);
	  max_relerror = mymax(max_relerror, tmp_relerror);

	  /***
	    if you want to keep track of which coefficient
	    has the greatest error (per loop)
	    ***/

	  if (max_error == tmp_error)
	    bad_coeff[i] = j;

	}
      error[i] = max_error;
      relerror[i] = max_relerror;
    }

  sum_error = 0.0 ; sum_relerror = 0.0; sum_bad_coeff = 0;
  for(i = 0 ; i < loops ; i ++)
    {
      sum_error += error[i];
      sum_relerror += relerror[i];
      sum_bad_coeff += bad_coeff[i];
    }
  ave_error = sum_error / ( (double) loops );
  ave_relerror = sum_relerror / ( (double) loops );
  ave_bad_coeff = ( (double) sum_bad_coeff ) / ( (double) loops );

  stddev_error = 0.0 ; stddev_relerror = 0.0;
  if ( (errorflag == 1) && (loops > 1) )
    {
      for( i = 0 ; i < loops ; i ++ )
	{
	  stddev_error += pow( ave_error - error[i] , 2.0 );
	  stddev_relerror += pow( ave_relerror - relerror[i] , 2.0 );
	}
      if (loops != 1)
	{
	  stddev_error = sqrt(stddev_error / ( (double) (loops - 1) ) );
	  stddev_relerror = sqrt(stddev_relerror / ( (double) (loops - 1) ) );
	}
    }

  fprintf(stderr,"error r-o:\t\t %.4e\t",ave_error);
  fprintf(stderr,"std dev: %.4e\n",stddev_error);
  fprintf(stderr,"rel error (r-o)/o:\t %.4e\t",ave_relerror);
  fprintf(stderr,"std dev: %.4e\n\n",stddev_relerror);
  /**
    fprintf(stderr,"bad_coef: %.4f\n\n", ave_bad_coeff);
    **/

  free(eval_args); free(data); free(result); free(coeffs);
  free(bad_coeff); free(relerror); free(error);
  free(cos_pml_table); free(workspace);

#ifdef FFTPACK
  free(wSave);
#endif


}

