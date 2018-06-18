%module swig_nufft 
%{
#define SWIG_FILE_WITH_INIT
int fourier_interpolate_3d(float *v, int dim1, int n1, int n2, int n3, double *x, int dim2, float *res, int dim3);
int fourier_rotate_vol(float *v, int dim1, int n1, int n2, int n3, double *rot, int dim2, float *res, int dim3);
int fourier_3d_iter_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3,
  float *kx, int dim4, float *ky, int dim5, float *kz, int dim6, float *out_real, int dim7, float *out_imag, int dim8, int iteration);
int fourier_2d1d_gridding_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3,
  float *kx, int dim4, float *ky, int dim5, float *out_real, int dim6, float *out_imag, int dim7);
int fourier_2d1d_iter_reconstruct(float *real, int dim1, float *imag, int dim2, int N, int Z, int M, float *weights, int dim3,
  float *kx, int dim4, float *ky, int dim5, float *out_real, int dim6, float *out_imag, int dim7, int iteration, double threshold, float *damping, int dim8);
%}
%include "numpy.i"
%init %{
    import_array();
%}

int fourier_interpolate_3d(float *IN_ARRAY1, int DIM1, int n1, int n2, int n3, double *IN_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1);

int fourier_rotate_vol(float *IN_ARRAY1, int DIM1, int n1, int n2, int n3, double *IN_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1);

int fourier_3d_iter_reconstruct(float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, int N, int Z, int M, float *IN_ARRAY1, int DIM1,
  float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1, int iteration);

int fourier_2d1d_gridding_reconstruct(float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, int N, int Z, int M, float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1);

int fourier_2d1d_iter_reconstruct(float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, int N, int Z, int M, float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *IN_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1, float *INPLACE_ARRAY1, int DIM1, int iteration, double threshold, float *IN_ARRAY1, int DIM1);