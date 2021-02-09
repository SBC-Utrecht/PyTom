/***************************************************************************
  **************************************************************************

                     Python interface of NFFT

   AUTHOR:
   Yuxiang Chen 
   Contact: chenyxkk@hotmail.com 

   Copyright 2011 Yuxiang Chen

   HOW TO USE:
   1. run mkswig.sh, which should generate shared library called _swig_nufft.so
   2. start python, import swig_nufft and enjoy!

  ************************************************************************
  ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

/**
  v: input, the input volume (linearlized)
  n1, n2, n3: input, dimension of the volume
  x: input, the fourier sampling nodes (linearlized)
  res: output, the result fourier coefficients at those nodes (linearlized)
 */
int fourier_interpolate_3d(float *v, int dim1, int n1, int n2, int n3, double *x, int dim2, float *res, int dim3)
{
  int N[3] = {n1, n2, n3};
  int n[3] = {n1*2, n2*2, n3*2}; // oversampling ratio = 2
  int M = dim2/3;
  int i;

  nfft_plan p;

  nfft_init_guru(&p, 3, N, M, n, 6, // window cutoff = 6
     PRE_PHI_HUT| PRE_PSI | MALLOC_F_HAT| MALLOC_X| MALLOC_F |
     FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  
  // assign f_hat
  for (i = 0; i < p.N_total; ++i)
  {
    p.f_hat[i] = v[i];
  }

  // assign nodes, for some stupid reason the x and z should be swtiched!
  for (i = 0; i < p.M_total; ++i)
  {
    p.x[3*i] = x[3*i+2];
    p.x[3*i+1] = x[3*i+1];
    p.x[3*i+2] = x[3*i];
  }
  
  if(p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);
  
  // do the transform
  nfft_trafo(&p);

  // assign output
  for (i = 0; i < p.M_total; ++i)
  {
    res[2*i] = creal(p.f[i]);
    res[2*i+1] = cimag(p.f[i]);
  }

  nfft_finalize(&p);

  return 1;
}

int fourier_rotate_vol(float *v, int dim1, int n1, int n2, int n3, double *rot, int dim2, float *res, int dim3)
{
  int i,j,k;
  int center_x = n1/2;
  int center_y = n2/2;
  int center_z = n3/2;

  // prepare the plan
  int N[3] = {n1, n2, n3};
  int n[3] = {n1*2, n2*2, n3*2}; // oversampling ratio = 2
  int M = n1*n2*(n3/2+1); // number of sampling points in Fourier space
  nfft_plan p;

  nfft_init_guru(&p, 3, N, M, n, 6, // window cutoff = 6
     PRE_PHI_HUT| PRE_PSI | MALLOC_F_HAT| MALLOC_X| MALLOC_F |
     FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  
  // assign f_hat
  for (i = 0; i < p.N_total; ++i)
  {
    p.f_hat[i] = v[i];
  }

  // calculate the nodes
  double xx, yy, zz;
  for (k = 0; k < n3/2+1; ++k) // calculate the half
  {
    for (j = 0; j < n2; ++j)
    {
      for (i = 0; i < n1; ++i)
      {
        // caculate the rotated position in original volume
        xx = (i-center_x)*rot[0]+(j-center_y)*rot[1]+(k-center_z)*rot[2];
        yy = (i-center_x)*rot[3]+(j-center_y)*rot[4]+(k-center_z)*rot[5];
        zz = (i-center_x)*rot[6]+(j-center_y)*rot[7]+(k-center_z)*rot[8];

        p.x[3*(k*n1*n2+j*n1+i)]   = xx/n1;
        p.x[3*(k*n1*n2+j*n1+i)+1] = yy/n2;
        p.x[3*(k*n1*n2+j*n1+i)+2] = zz/n3;
      }
    }
  }

  if(p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);

  // do the transform
  nfft_trafo(&p);

  // assign output
  int x,y,z;
  for (i = 0; i < p.M_total; ++i)
  {
    res[2*i] = creal(p.f[i]);
    res[2*i+1] = cimag(p.f[i]);

    // fill the other half
    z = i / (n1*n2); // get the current point 3d coordinate
    if(z == 0 && (n3%2) == 0)
      continue; // do not need to map this dimension
    y = (i - n1*n2*z) / n1;
    x = i - n1*n2*z - n1*y;
    x = (2*center_x - x) % n1; // get the counter point 3d coordinate
    y = (2*center_y - y) % n2;
    z = 2*center_z - z;
    j = z*n1*n2+y*n1+x; // get the linear index back
    res[2*j] = creal(p.f[i]);
    res[2*j+1] = -cimag(p.f[i]); // set the conjugate
  }

  nfft_finalize(&p);

  return 1;
}

/**
  Fourier iterative reconstruction method.
  real: real part of Fourier coefficients
  imag: imag part of Fourier coefficients
  N: image height
  Z: number of slices
  M: number of knots in Fourier space
  weights: knots weights in Fourier space
  kx: knots position in Fourier space
  ky: knots position in Fourier space
  kz: knots position in Fourier space
  out_real: real part of result
  out_imag: imag part of result
  iteration: number of iterations to run
  */
int fourier_3d_iter_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3, float *kx, int dim4, float *ky, int dim5, float *kz, int dim6, float *out_real, int dim7, float *out_imag, int dim8, int iteration)
{
  int j,k,z,l;                  /* some variables  */
  nfft_plan my_plan;            /* plan for the two dimensional nfft  */
  solver_plan_complex my_iplan;          /* plan for the two dimensional infft */
  int my_N[3],my_n[3];          /* to init the nfft */
  double epsilon=0.0000003;     /* tmp to read the obsolent z from 700.acs
                                   epsilon is a the break criterion for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft */

  /* initialise my_plan, specific.
     we don't precompute psi */
  my_N[0]=Z; my_n[0]=ceil(Z*1.2); // for some reason this oversampling factor suffices
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  my_N[2]=N; my_n[2]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 3, my_N, M, my_n, 6,
                      PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* precompute lin psi */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);

  infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  solver_init_advanced_complex(&my_iplan,(mv_plan_complex*)(&my_plan), infft_flags);

  /* set the weights */
  for(j=0;j<M;j++)
  {
    my_iplan.w[j] = weights[j];
  }

  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        for(z=0;z<N;z++) {
        int j2= j-N/2;
        int k2= k-N/2;
        int z2= z-N/2;
        double r=sqrt(j2*j2+k2*k2+z2*z2);
        if(r>(double) N/2)
          my_iplan.w_hat[z*N*N+j*N+k]=0.0;
        else
          my_iplan.w_hat[z*N*N+j*N+k]=1.0;
        }
      }
    }
  }

  /* set x,y,freal and fimag from the knots */
  for(j=0;j<M;j++)
  {
    // the order should be correct now
    // the only confusing thing is the inconsistency with previous functions?
    my_plan.x[3*j] = kx[j];
    my_plan.x[3*j+1] = ky[j];
    my_plan.x[3*j+2] = kz[j];
    my_iplan.y[j] = real[j] + _Complex_I*imag[j];
  }

  /* precompute psi */
  if(my_plan.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&my_plan);

  /* precompute full psi */
  if(my_plan.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k]=0.0;

  /* inverse trafo */
  solver_before_loop_complex(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
      break;
    // fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),l+1,iteration);
    solver_loop_one_step_complex(&my_iplan);
  }

  for(l=0;l<Z;l++)
  {
    for(k=0;k<N*N;k++)
    {
      /* write every Layer to the output */
      out_real[k+l*N*N] = creal(my_iplan.f_hat_iter[ k+N*N*l ]);
      out_imag[k+l*N*N] = cimag(my_iplan.f_hat_iter[ k+N*N*l ]);
    }
  }

  solver_finalize_complex(&my_iplan);
  nfft_finalize(&my_plan);

  return 1;
}

int fourier_2d1d_gridding_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3, float *kx, int dim4, float *ky, int dim5, float *out_real, int dim6, float *out_imag, int dim7)
{
  fftw_complex *mem;
  fftw_plan plan;

  /* Allocate memory to hold every slice in memory after the 2D-infft */
  mem = (fftw_complex*) nfft_malloc(sizeof(fftw_complex) * N * N * Z);

  /* Create plan for the 1d-ifft */
  plan = fftw_plan_many_dft(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_BACKWARD, FFTW_MEASURE);



  /* execute the 2d-nfft's */
  int i,j,k,z;               /* some variables  */
  nfft_plan my_plan;       /* plan for the two dimensional nfft  */
  int my_N[2],my_n[2];     /* to init the nfft */

  /* initialise my_plan */
  my_N[0]=N; my_n[0]=ceil(N*1.2); // for some reason this oversampling factor suffices
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 2, my_N, M/Z, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                        MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                        FFTW_INIT| FFT_OUT_OF_PLACE,
                        FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* precompute lin psi if set */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);

  for(z=0;z<Z;z++) {
    // set the input value, knots and weights
    for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j+0] = kx[j];
      my_plan.x[2*j+1] = ky[j];
      my_plan.f[j] = (real[M/Z*z+j] + _Complex_I*imag[M/Z*z+j])*weights[j];
    }

    /* precompute psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    /* precompute full psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

    /* compute the adjoint nfft */
    nfft_adjoint(&my_plan);

    for(k=0;k<my_plan.N_total;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)] = my_plan.f_hat[k];
    }
  }

  nfft_finalize(&my_plan);



  /* execute the 1d-fft's */
  fftw_execute(plan);

  /* fill in the results */
  for(i=0;i<Z;i++) {
    for (j=0;j<N*N;j++) {
      out_real[i*N*N+j] = creal(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z;
      out_imag[i*N*N+j] = cimag(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z;
    }
  }

  /* free memory */
  nfft_free(mem);
  fftw_destroy_plan(plan);

  return 1;
}


int fourier_2d1d_iter_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3, float *kx, int dim4, float *ky, int dim5, float *out_real, int dim6, float *out_imag, int dim7, int iteration, double threshold, float *damping, int dim8)
{
  fftw_complex *mem;
  fftw_plan plan;

  /* Allocate memory to hold every slice in memory after the 2D-infft */
  mem = (fftw_complex*) nfft_malloc(sizeof(fftw_complex) * N * N * Z);

  /* Create plan for the 1d-ifft */
  plan = fftw_plan_many_dft(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_BACKWARD, FFTW_MEASURE);



  /* execute the 2d-nfft's */
  int i,j,k,l,z;               /* some variables  */
  nfft_plan my_plan;       /* plan for the two dimensional nfft  */
  solver_plan_complex my_iplan; /* plan for the two dimensional infft */
  int my_N[2],my_n[2];     /* to init the nfft */
  double epsilon=0.0000003;
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT; /* flags for the infft */
  // unsigned infft_flags = CGNR | PRECOMPUTE_WEIGHT; /* flags for the infft */

  if(threshold > 0)
  {
    epsilon = threshold;
  }

  /* initialise my_plan */
  my_N[0]=N; my_n[0]=ceil(N*1.2); // for some reason this oversampling factor suffices
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 2, my_N, M/Z, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                        MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                        FFTW_INIT| FFT_OUT_OF_PLACE,
                        FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* precompute lin psi if set */
  if(my_plan.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);

  /* initialise my_iplan, advanced */
  solver_init_advanced_complex(&my_iplan,(mv_plan_complex*)(&my_plan), infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    for(j=0;j<my_plan.M_total;j++)
    {
        my_iplan.w[j] = weights[j];
    }
  }

  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        my_iplan.w_hat[j*N+k] = damping[j*N+k];
        // int j2= j-N/2;
        // int k2= k-N/2;
        // double r=sqrt(j2*j2+k2*k2);
        // if(r>(double) N/2)
        //   my_iplan.w_hat[j*N+k]=0.0;
        // else
        //   my_iplan.w_hat[j*N+k]=1.0;
      }
    }
  }
  double residual = 0.0;
  for(z=0;z<Z;z++) {
    // fprintf(stderr,"%d\n",z);
    // set the input value, knots
    for(j=0;j<my_plan.M_total;j++)
    {
      my_plan.x[2*j+0] = kx[j];
      my_plan.x[2*j+1] = ky[j];
      my_iplan.y[j] = real[M/Z*z+j] + _Complex_I*imag[M/Z*z+j];
    }

    /* precompute psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.nfft_flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    /* precompute full psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.nfft_flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

    /* init some guess */
    for(k=0;k<my_plan.N_total;k++)
      my_iplan.f_hat_iter[k]=0.0;

    /* inverse trafo */
    solver_before_loop_complex(&my_iplan);
    for(l=0;l<iteration;l++)
    {
      /* break if dot_r_iter is smaller than epsilon*/
      if(my_iplan.dot_r_iter<epsilon)
      break;
      // fprintf(stdout,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),iteration*z+l+1,iteration*Z);
      solver_loop_one_step_complex(&my_iplan);
    }
    residual += my_iplan.dot_r_iter;

    for(k=0;k<my_plan.N_total;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)] = my_iplan.f_hat_iter[k];
    }
  }
  // fprintf(stderr,"%e\n",sqrt(residual));

  /* finalize the infft */
  solver_finalize_complex(&my_iplan);

  /* finalize the nfft */
  nfft_finalize(&my_plan);



  /* execute the 1d-fft's */
  fftw_execute(plan);

  /* fill in the results */
  for(i=0;i<Z;i++) {
    for (j=0;j<N*N;j++) {
      out_real[i*N*N+j] = creal(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z;
      out_imag[i*N*N+j] = cimag(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z;
    }
  }

  /* free memory */
  nfft_free(mem);
  fftw_destroy_plan(plan);

  return 1;
}

