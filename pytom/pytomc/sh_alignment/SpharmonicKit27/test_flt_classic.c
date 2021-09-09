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


/******

   Legendre transform example procedure. Basically, this is
   the routine which times and tests for error the Legendre
   transform algorithm FLT_classic .

   The major reference is the tech report

       FFTs for the 2-Sphere - Improvements and Variations

       by D.M. Healy, Jr., D. Rockmore and Sean S.B. Moore
       Department of Computer Science
       Dartmouth College
       Technical Report PCS-TR96-292


   Also, see precomp_flt_classic.c, flt_classic.c for details.

   Sample call:

   % test_flt_classic m bw k numsplits errflag

   m = order, bw = bandwidth, k = # of loops
   numsplits = how many times want to divide in the
               divide-and-conquer, CANNOT EXCEED 5

	       Example: if bw = 256, m = 0 and numsplits = 2,
	       then the problem of size 256 is reduced to four
	       problems, each of size 64
	       ( one 256 -> two 128 -> four 64)

   errflag = 0 if want to time, = 1 if want to determine
             errors (std devs, etc)


****/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flt_classic.h"
#include "primitive.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))
#define DATAFILE "norm.dat"   /* file which contains random numbers */

void main(int argc, char **argv)
{
  int m, bw, i, j, k, numsplits;
  double *coeffs, *data, *eval_args, *workspace, *result, *zspace, **Z;
  FILE *fp;
  struct lowhigh *lowpoly, *highpoly;

  double tmp_error, max_error, ave_error, sum_error;
  double tmp_relerror, max_relerror, ave_relerror, sum_relerror;
  double stddev_error, stddev_relerror;
  int *bad_coeff, sum_bad_coeff, errorflag;
  double *error, *relerror, ave_bad_coeff;

  if (argc != 6 )
    {
      fprintf(stdout,"Usage: test_flt_classic m bw k numsplits errflag\n");
      exit(0);
    }

  m  = atoi(argv[1]);
  bw = atoi(argv[2]);
  k  = atoi(argv[3]);   /* number of loops */
  numsplits = atoi(argv[4]);
  errorflag = atoi(argv[5]);

  /* space for error stuff */
  error = (double *) malloc(sizeof(double) * k);
  relerror = (double *) malloc(sizeof(double) * k);
  bad_coeff = (int *) malloc(sizeof(int) * k);

  /* allocate space for coefficients */
  coeffs = (double *) malloc(sizeof(double) * (bw - m) *
			     (k*errorflag + 1));

  /* space for result of Legendre transform */
  result = (double *) malloc(sizeof(double) * (bw - m) *
			     (k*errorflag + 1));

  /* space for data */
  data = (double *) malloc(sizeof(double) * 2 * bw *
			   (k*errorflag + 1));
    

  /* used to synthesize a test function using P_eval */
  eval_args = (double *) malloc(sizeof(double) * 2 * bw);

  /* forward transform needs scratch workspace
     of size 20 * bw */
  workspace = (double *) malloc(sizeof(double) * 20 * bw);

  /* This is space for the intermediate results of the Legendre
     transform algorithm. That is, this space stores the results
     of the divide-and-conquer portion of the algorithm 
     */

  zspace = (double *) malloc(sizeof(double) * (12 * bw));

  /* set up Z space; need an array of double pointers */
  Z = (double **) malloc(sizeof(double *) * 6 ); 

  Z[0] = zspace;
  for ( i = 1 ; i < 6 ; i ++ )
    Z[i] = Z[i - 1] + (2 * bw);

  /** needed for FCT **/
  lowpoly = (struct lowhigh *) malloc(sizeof(struct lowhigh) * (bw));
  highpoly = (struct lowhigh *) malloc(sizeof(struct lowhigh) * (bw));

  /**********************************************/
  /* Here I'm just setting all the coefficients */
  /**********************************************/
  /*
    coeffs[0] = 0.0;
    for (i=1; i<bw; i++)
    coeffs[i] = (double) i;
    */

  fp = fopen(DATAFILE,"r");
  if ( fp == NULL )
    {
      fprintf( stderr,"Error in test_flt_classic -->\n" );
      perror( DATAFILE );
      exit( 1 ) ;
    }
  for (j=0; j<(bw-m)*(k*errorflag+1); j++)
    fscanf(fp, "%lf", coeffs+j);
  fclose(fp);

  /* generate arguments (arccos of sample points) for synthesizing
     a test function */
  ArcCosEvalPts(2*bw,eval_args);

  /* synthesize the corresponding function */
  for(j = 0 ; j < k*errorflag+1 ; j++)
    P_eval(m, coeffs+j*(bw-m), eval_args,
	   data+j*2*bw, workspace, bw);

  FLT_classic(data, result, m, bw, k,
	      Z, numsplits, workspace,
	      lowpoly, highpoly, errorflag);

  /* now to compute the error     */
  /* k = number of loops */
  /*** have to reset k if I'm looping over the same data
    (which I am when timing) ***/
    
  if(errorflag == 0)
    k = 1;
    
  for(i = 0 ; i < k ; i ++)
    {
      max_error = 0.0 ; max_relerror = 0.0;
      for(j = 0 ; j < (bw - m) ; j ++)
	{
	  tmp_error = fabs(coeffs[i*(bw-m)+j] - result[i*(bw-m)+j]);
	  tmp_relerror = tmp_error /(fabs(coeffs[i*(bw-m)+j])+
				     pow(10.0, -50.0));
	  max_error = mymax(max_error, tmp_error);
	  max_relerror = mymax(max_relerror, tmp_relerror);
	    
	  if (max_error == tmp_error)
	    bad_coeff[i] = j;
	    
	}
      error[i] = max_error;
      relerror[i] = max_relerror;
    }
    
  sum_error = 0.0 ; sum_relerror = 0.0; sum_bad_coeff = 0;
  for(i = 0 ; i < k ; i ++)
    {
      sum_error += error[i];
      sum_relerror += relerror[i];
      sum_bad_coeff += bad_coeff[i];
    }
  ave_error = sum_error / ( (double) k );
  ave_relerror = sum_relerror / ( (double) k );
  ave_bad_coeff = ( (double) sum_bad_coeff ) / ( (double) k );
    
  stddev_error = 0.0 ; stddev_relerror = 0.0;
  if ( errorflag == 1 )
    {
      for( i = 0 ; i < k ; i ++ )
	{
	  stddev_error += pow( ave_error - error[i] , 2.0 );
	  stddev_relerror += pow( ave_relerror - relerror[i] , 2.0 );
	}
      if ( k != 1 )
	{
	  stddev_error = sqrt(stddev_error / ( (double) (k - 1) ) );
	  stddev_relerror = sqrt(stddev_relerror / ( (double) (k - 1) ) );
	}
    }
    
  fprintf(stderr,"error r-o:\t\t %.4e\t",ave_error);
  fprintf(stderr,"std dev: %.4e\n",stddev_error);
  fprintf(stderr,"rel error (r-o)/o:\t %.4e\t",ave_relerror);
  fprintf(stderr,"std dev: %.4e\n\n",stddev_relerror);

  free(highpoly);free(lowpoly);free(Z);
  free(zspace);free(workspace);
  free(eval_args);free(data);free(result);free(coeffs);
  free(bad_coeff);free(relerror);free(error);
    
    
}


