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
/*************************************************************************/

/* This is code to debug the InvSemiNaiveReduced
   routine.

   Start by reading in a file of Legendre coefficients
   Synthesize the function f.
   Apply forward seminaive algorithm - should get Legendre coeffs.
   Apply inverse seminaive algorithm - should get f.


   Sample calls: 

   test_semi_roundtrip m bw coeffs_file output_file

   test_semi_roundtrip 0 128 coeffs_flat.mca output_file.dat

 */

#include <math.h>
#include <stdio.h>

#include "naive_synthesis.h"
#include "cospmls.h"
#include "seminaive.h"

#ifdef FFTPACK
#include "newFCT.h"
#include "fftpack.h"
#endif


int main(argc, argv)
     int argc;
     char **argv;
{
  int n, i, m, bw;
  double *coeffs;
  double *function, *semiresult, *invsemiresult;
  double *cospmltable, *cos_even, *sin_values, *eval_pts;
  double *trans_cospml_table;
  double *workspace;
  FILE *infile, *outfile;
  double *wSave;

  if (argc != 5)
    {
      fprintf(stdout,"Usage: test_semi_roundtrip m bw input_file output_file.\n");
      return(0);
    }


  m = atoi(argv[1]);
  bw = atoi(argv[2]);
  infile = fopen(argv[3], "r");

  
  outfile = fopen(argv[4], "w");
  
  /*
  outfile = stdout;
  */

  fprintf(outfile,"bw = %d\n", bw);
  fprintf(outfile,"m = %d\n", m);

#ifdef FFTPACK
  n = 2 * bw;

  wSave = precomp_dct( n );

#endif

  coeffs = (double *) malloc(sizeof(double) * bw);
  function = (double *) malloc(sizeof(double) * 2 * bw);
  semiresult = (double *) malloc(sizeof(double) * bw);
  invsemiresult = (double *) malloc(sizeof(double) * 2 * bw);
  cospmltable = (double *) malloc(sizeof(double) * TableSize(m,bw));
  trans_cospml_table = (double *) malloc(sizeof(double) * TableSize(m,bw));
  cos_even = (double *) malloc(sizeof(double) * bw);
  sin_values = (double *) malloc(sizeof(double) * 2 * bw);
  eval_pts = (double *) malloc(sizeof(double) * 2 * bw);
  workspace = (double *) malloc(sizeof(double) * 32 * bw);

  /* read in coeffs */

  for (i=0; i<(bw-m); i++)
    fscanf(infile, "%lf", coeffs+i);

  /* print them out */
  fprintf(outfile, "Input Legendre coefficients\n");
  for (i=0; i<bw-m; i++)
    fprintf(outfile, "%f\n", coeffs[i]);

  /* synthesize the corresponding function */

  Naive_Synthesize(bw, m, coeffs, function, workspace);

  /* print out function */
  fprintf(outfile, "Synthesized function\n");
  for (i=0; i<2*bw; i++)
    fprintf(outfile, "%f\n", function[i]);

  /* now do the SemiNaive transform */
  
  /* get the table of cospmls */
  CosPmlTableGen(bw, m, cospmltable, workspace);

  /* print out the cospmltable */
  /*
  fprintf(outfile,"CosPmlTable values\n");
  tablesize_limit = TableSize(m,bw);
  for (i=0; i<tablesize_limit; i++)
    fprintf(outfile,"%lf\n",cospmltable[i]);
  */


  /* compute the table of sin_values */
  
  ArcCosEvalPts(2*bw, eval_pts);
  for (i=0; i<2*bw; i++)
    sin_values[i] = sin(eval_pts[i]);

  /* do it */
#ifndef FFTPACK
  SemiNaiveReduced(function, 
		   bw, 
		   m, 
		   semiresult, 
		   cospmltable, 
		   workspace,
		   cos_even);
#else
  SemiNaiveReduced(function, 
		   bw, 
		   m, 
		   semiresult, 
		   cospmltable, 
		   workspace,
		   cos_even,
		   wSave);
#endif
  /* print out the result */

  fprintf(outfile, "Legendre coeffs computed by SemiNaiveReduced\n");
  for (i=0; i<(bw-m); i++)
    fprintf(outfile,"%f\n",semiresult[i]);

  /* now apply the inverse seminaive transform */

  /* get the transposed table values */

  Transpose_CosPmlTableGen(bw, m, cospmltable, trans_cospml_table);

  /* print out the costrans_pmltable */
  /*
  fprintf(outfile,"Transpose_CosPmlTable values\n");
  for (i=0; i<tablesize_limit; i++)
    fprintf(outfile,"%lf\n",trans_cospml_table[i]);
  */

  InvSemiNaiveReduced(semiresult, 
		      bw, 
		      m, 
		      invsemiresult, 
		      trans_cospml_table, 
		      sin_values,
		      workspace);

  /* print it */
  fprintf(outfile,"Function computed by InvSemiNaiveReduced\n");
  for (i=0; i<2*bw; i++)
    fprintf(outfile,"%f\n", invsemiresult[i]);

  /* amscray */

  fclose(infile);
  fclose(outfile);

  free(workspace); free(eval_pts); free(sin_values);
  free(cos_even); free(trans_cospml_table);
  free(cospmltable); free(invsemiresult);
  free(semiresult); free(function);
  free(coeffs);

#ifdef FFTPACK
  free(wSave);
#endif

  return(0);


}
