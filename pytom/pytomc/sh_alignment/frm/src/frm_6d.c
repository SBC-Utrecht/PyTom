/**********************************************************************
 *                              FRM6D    2/20/03                      *
 **********************************************************************
 * Program is part of the Situs package URL: http://situs.scripps.edu *
 * (c) Julio Kovacs, Pablo Chacon and Willy Wriggers, 2001-2002       *
 **********************************************************************
 *                                                                    *
 * General 6D (5 angles + rho ) rotational accelerated fitting tool.  * 
 *                                                                    *
 **********************************************************************
 * Due to the use of the FFTW library we are forced to distribute     *
 * this file under the GNU General Public License (GPL) appended to   *
 * the source code below. The GNU GPL license replaces the legal      *
 * statement distributed with the Situs package.                      *
 *********************************************************************/


#include "situs.h"
#include <rfftw.h>
#define max(A,B) ((A)>(B) ? (A):(B))
#define min(A,B) ((A)<(B) ? (A):(B))


/* global variables */
static double g_sqt2;                        /* sqrt(2) */
static double g_target_res;                  /* resolution in A */
static double g_target_ani;                  /* resolution anisotropy factor */
static double g_low_cutoff;                  /* low res. map cutoff */
static int g_corr_mode;                      /* correlation mode */
static int g_rho_mode;                       /* rho scanning mode */
static unsigned g_extx;                      /* map extent */
static unsigned g_exty;                      /* map extent */
static unsigned g_extz;                      /* map extent */
static unsigned long g_nvox;                 /* number of voxels */
static double g_width;                       /* voxel size in Angstroms */
static double g_gridx;                       /* map origin */
static double g_gridy;                       /* map origin */
static double g_gridz;                       /* map origin */
static double *g_phi_lo;                     /* low resolution map */
static double *g_phi_hi;                     /* high resolution map */
static double *g_phi_du;                     /* dummy map */
static double *g_phi_la;                     /* filter kernel */
static double *g_phi_ga;                     /* low-pass (Gaussian) kernel */
static double *g_phi_fx;                     /* filtered low-pass kernel */
static double g_norm_hi, g_norm_lo;          /* normalization factors */
static unsigned g_ext_ga;                    /* low-pass kernel linear extent */
static unsigned g_ext_la;                    /* filter kernel linear extent */
static unsigned g_ext_fx;                    /* filtered low-pass kernel linear extent */
static unsigned long g_nvox_ga;              /* low-pass kernel voxel count */
static unsigned long g_nvox_la;              /* filter kernel voxel count */
static unsigned g_ignored[3];                /* zero margin ignored in fast kernel convolution */  

static FILE *g_outfile;                      /* frequently used output file */ 



static unsigned g_num_atoms;                 /* number of PDB atoms */
static PDB *g_pdb_original;                  /* PDB coordinates */
static PDB *g_pdb_move;                      /* dummy PDB coordinates */
static unsigned g_num_explored;              /* number of explored best fits */



/* external functions */
extern void sgenrand (unsigned long);
extern double genrand ();
extern void powell (int *, double *, double *, int, double(*)(double *), unsigned, FILE *, double, double *);


extern void zero_vect (double *, unsigned long);
extern void do_vect (double **, unsigned long);
extern void cp_vect (double **, double **, unsigned long);
extern void cp_vect_destroy (double **, double **, unsigned long);
extern void create_padded_map (double **, unsigned *, unsigned *, unsigned *, double *, double *, 
			       double *, unsigned long *, double *, unsigned, unsigned, unsigned, 
			       double, double, double, double, double, double, unsigned *);
extern void shrink_margin (double **, unsigned *, unsigned *, unsigned *, double *, double *, 
			   double *, unsigned long *, double *, unsigned, unsigned, unsigned, 
			   double, double, double, double, double, double);
extern void shrink_to_sigma_factor (double **, unsigned *, double *, unsigned, double, double);
extern void interpolate_map (double **, unsigned *, unsigned *, unsigned *, double *, double *, double *, 
			     double, double, double, double *, unsigned, unsigned, unsigned, double, double, double,
			     double, double, double);
extern double calc_total (double *, unsigned long);
extern double calc_norm (double *, unsigned long);
extern void print_map_info (double *, unsigned long);
extern void threshold (double *, unsigned long, double);
extern void normalize (double *, unsigned long, double);
extern void create_gaussian (double **, unsigned long *, unsigned *, double, double);
extern void create_identity (double **, unsigned long *, unsigned *);
extern void create_laplacian (double **, unsigned long *, unsigned *);
extern void convolve_kernel_inside (double **, double *, unsigned, unsigned, unsigned, double *, unsigned);  
extern void convolve_kernel_inside_erode (double **, double *, unsigned, unsigned, unsigned, double *, unsigned);  
extern void convolve_kernel_inside_fast (double **, double *, unsigned, unsigned, 
					 unsigned, double *, unsigned, double, unsigned *);  

extern double calc_sphere (PDB *, unsigned, double, double, double);
extern double calc_mass (PDB *, unsigned); 
extern double calc_center_mass (PDB *, unsigned, double *, double *, double *);
extern void rot_euler (PDB *, PDB *, unsigned, double, double, double);
extern void translate (PDB *, PDB *, unsigned, double, double, double);
extern void project_mass (double **, unsigned long, double, double, double, 
			  unsigned, unsigned, unsigned, PDB *, unsigned, double *);
extern void read_and_weight_mass (char *, unsigned *, PDB **); 
extern void get_rot_matrix (double [3][3], double, double, double);
extern the_time get_the_time (void);
extern the_time time_diff (the_time, the_time);
extern char *smart_sprint_time (double);
extern double time_to_sec (the_time);
extern void readvol (char *, double *, double *, double *, double *, int *, int *, int *, double **);
extern void writevol (char *, double, double, double, double, int, int, int, double *);
extern void appendpdb (char *, int, PDB *);

extern void relax_laplacian (double **,  unsigned , unsigned , unsigned, unsigned*, double);


/* from sphkit */
extern int seanindex(int, int, int);


/* list of functions defined here */
static void read_options (int,char **);
static void get_centered_structure_and_radius (char *,double *);
static void get_lowmap(char *);
static void get_highmap(char *);
static void draw_line ();
static void cen_of_mass(double *, double *, double *, double *);
static double longitude(int, int);
static double latitude(int, int);
static void wigner(int, double, double *);
static int prime(int, int);

static void legendre(int, int, double *);

static void two_cent_int(int, int, int, int, int, int, int,
                         double[], long[], double *,
                         double *[], double *[], double *[], double *[],
                         double *, double *);
static void euler_comp(double, double, double, double, double, double,
                       double *, double *, double *);
static double dist_mot(double, long, double *[], long, double *[]);

static void volumetric6d_sampling(double *, double, double, double, double *[], int);

static void rmaxmin(double *[], int, int[]);
static void threshold_down (double *, unsigned long, double);

main(int argc, char **argv){

  /* general paramenters */
  int log2s;         /* log2 of size */
  int tlb;           /* log2(bw^2)=2*log2(bw) */
  int bw;            /* bandwidth for sampling on the sphere */
  int bx;            /* bandwidth for sampling the d^l_{n0}(beta) */
  int size;          /* no. of points in each angular coordinate = 2*bw */
  int spp;           /* size+2 */

  /* timers */
  the_time itime, etime;
  the_time time00, time01, time02, time03, time04, time05, time06, time07, time08, time09, time10, time11; 

  double struct_max_radius; 
  double sigma1d;                  /* sigma of the Gaussian used for blurrong */
  double zero_shift[3] = {0.0, 0.0, 0.0};
  int zeropad[3] = {0, 0, 0};
  double grid_extent_pdb;
  char out_string[20], corr_string[70];

  int bin, h, hp, m, mp, n, jp, kp, flag1, flag2;
  long i, j, k, l, lp, tl, ind1, ind2, ind3, ind4, r2, size2, size3;
  long size4, size5, tsize5, ssize, ssize2, ssize3, ssize4;
  long npeaks, npeakx, npeakxa, ncorr, ncorrx;
  long offset, offsetp, offset2, offset3, offset4, offset5, factor, factorp;
  unsigned long init, index, indexb, index1, index2, index3, index3b, itemp;

  double comx_hi, comy_hi, comz_hi;  /* center of mass of high res. map */
  double comx_lo, comy_lo, comz_lo;  /* center of mass of low res. map */

 
  /* radial integral variables */
  int ifit, isign, isigm, r, r_a, r_b;
  int rmaxmin_hi[6], rmaxmin_lo[6];     /* max and min r's */
  int rhomin, rhomax, rhobar, rhobar0, rhobarini, rhobarmin, rhobarmax;
  double rho, rhoini;
  double paramf[5];     /* parameters to be passed to two_cent_int */
  long parami[7];
  double integralr, integrali;     /* real & imag parts of the 2-center integrals */

  rfftwnd_plan plan_fftw;          /* plan for the FFTW routine */

  /* variables used in the peak search */
  double temp, tempi, tempj, tempk, tempjp, tempkp, tempd;
  double sqll;                     /* to store sqrt((l+0.5)*(lp+0.5)) */
  double f00000, f10000, f20000, f01000, f02000, f00100, f00200;
  double f00010, f00020, f00001, f00002;
  double maxpeak, avepeak, lpeak;  /* of the correlation function */
  double maxcor[1000];             /* correlation maxima for each rho */
  double maxcor_all;               /* max correlation found so far */
  double stepsize;                 /* for integration (unless <1) */

  double phi, th, psi, phip, thp, psip;      /* Euler angles */
  double sigma, eta, ohm, etap, ohmp;        /* angular parameters */
  double phic, thc, psic;                    /* Euler angles of R^{-1}*R' */
  double vx, vy, vz;                         /* translation vector R^{-1}*(0,0,rho) */

  double debe;               /* delta beta = M_PI/bw */
  double dist_cut;           /* cutoff distance for eliminating nearby peaks */
  int *szs;                  /* number of sampling points for each angular parameter */
  int  max_r;                /* maximum number of r values for the sampling the densities, minus 1 */


  /* variables used in the spherical harmonic analysis with SpharmonicKit */
  int cut_naive;     /* switch from semi-naive to naive at order bw - cut_naive */
  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table;
  double **sph_sampl_lo;  /* array that will contain samplings of the low res. map */
                          /* [0] corresponds to r=0, [1] to r=1 and so on */
  double **sph_sampl_hi;  /* array that will contain samplings of the high res. map */
                          /* [0] corresponds to r=0, [1] to r=1 and so on */
  double *idata;                    /* imaginary part of data */
  double **coeffr_hi, **coeffi_hi;  /* spherical harmonics coefficients for high res. map */
  double **coeffr_lo, **coeffi_lo;  /* spherical harmonics coefficients for low res. map */


  /* d-functions variables */
  double *ddd;               /* piramid array that will contain the d-functions */
  double *dd0;               /* prism array for the d^l_{n,0}(beta_j) */
  double *dtd, *dtdp;        /* products of 2 d-functions */
  long *inits;               /* initial indices in ddd of d^l_{mn}, for l=0,...,bw-1 */


  /* correlation-related variables */
  float *coef_corr;          /* Fourier coefficients of correlation function */
  double *corr[7];           /* correlation peaks at each rho, and corresponding Euler angles */
  double *corrx[7];          /* final correlation maxima (for all rho) and Euler angles */
                             /* the actual # of valid peaks in it is npeakx */
  long *inco;                /* indices of the peaks in corr */
  long *incox;               /* indices of the peaks in corrx */
  int coco[1000];            /* flags indicating whether the correlations for each rho have
                                been computed: 1:they have; 0:they haven't */

 

  /* ===================== INITIALIZATION =====================*/ 

  g_sqt2=sqrt(2.0);

  itime = get_the_time();

  /* set a variety of global variables */
  draw_line();  
  read_options(argc,&(*argv)); 

  /* load data */
  draw_line();
  get_lowmap(argv[1]);         
  draw_line(); 
  get_centered_structure_and_radius (argv[2], &struct_max_radius);
  draw_line();

  /* calculate sigma of kernel (1D) */
  sigma1d = g_target_res / (2.0*g_width*sqrt(3.0)); 

  /* create Gaussian kernel */  
  create_gaussian (&g_phi_ga, &g_nvox_ga, &g_ext_ga, sigma1d, 3.0);
 

  /* create filter (_fi) kernel (e.g. Laplacian) and indicate (_fx) sigma factor */
  switch (g_corr_mode) {
  case 0: 
    strcpy(corr_string,"REMARK    Standard Linear Correlation"); 
    break;
  case 1:
    create_laplacian (&g_phi_la, &g_nvox_la, &g_ext_la);
    strcpy(corr_string,"REMARK    Laplacian Correlation");
    break;
  }


  /* check dimensions of the grid respect pdb max radius */
  grid_extent_pdb = ceil(struct_max_radius/g_width)+1; 
  if (grid_extent_pdb > g_extx)  zeropad[0] = ceil((grid_extent_pdb-g_extx)/2);
  if (grid_extent_pdb > g_exty)  zeropad[1] = ceil((grid_extent_pdb-g_exty)/2);
  if (grid_extent_pdb > g_extz)  zeropad[2] = ceil((grid_extent_pdb-g_extz)/2);

  /* save g_ignored margin width */
  for (m=0;m<3;m++) g_ignored[m] = zeropad[m];

  /* finally, we need to add half of the convolution kernel size to zeropad */
  for (m=0;m<3;m++) zeropad[m] += (g_ext_ga-1)/2;


  /* now add extra padding to g_phi_lo, changing the map size parameters */  
  create_padded_map (&g_phi_du, &g_extx, &g_exty, &g_extz, &g_gridx, &g_gridy, &g_gridz, &g_nvox, 
		     g_phi_lo, g_extx, g_exty, g_extz, g_gridx, g_gridy, g_gridz, 
		     g_width, g_width, g_width*g_target_ani, zeropad); 

  /* get the enlarged g_phi_lo */
  cp_vect_destroy(&g_phi_lo,&g_phi_du,g_nvox);
 
  draw_line();

  printf("frm6d> Projecting probe structure onto lattice...\n");
  do_vect(&g_phi_hi,g_nvox); 
  project_mass (&g_phi_hi, g_nvox, g_width, g_width, g_width*g_target_ani, g_extx, g_exty, g_extz, g_pdb_original, g_num_atoms, zero_shift);

  printf("frm6d> Low pass filtering probe map...\n");
  convolve_kernel_inside (&g_phi_hi, g_phi_hi, g_extx, g_exty, g_extz, g_phi_ga, g_ext_ga);

  /* normalization */
  printf("frm6d> Normalizing target and probe maps...\n");
  g_norm_hi = calc_norm(g_phi_hi, g_nvox);
  g_norm_lo = calc_norm(g_phi_lo, g_nvox);
  normalize(g_phi_hi, g_nvox, g_norm_hi); 
  normalize(g_phi_lo, g_nvox, g_norm_lo); 

  printf("frm6d> Target and probe maps:\n");
  print_map_info(g_phi_lo, g_nvox);
  print_map_info(g_phi_hi, g_nvox);


  /* ==================== FRM  PRECOMPUTATIONS =======================*/ 
 
  draw_line();
  fprintf(stderr, "frm6d> Precomputations for FST and d-matrices...\n");

  /* set parameters */
  log2s=5;                    /* log2(size) */
              /* size=32 is the maximum possible for now due to memory limitations */
              /* size=32 needs about 150MB, while size=64 needs 4.5GB */
  size=ldexp(1,log2s)+0.1;
  tlb=2*(log2s-1);            /* log2(bw^2) */
  bw=size/2;                  /* bandwidth */
  bx=4*bw;                    /* analog to bw for sampling the angular argument
                                 of the Legendre functions */
  debe=M_PI/bw;
  dist_cut=3.0;
  cut_naive=bw;
  size2=size*size;
  size3=size2*size;
  size4=size3*size;
  size5=size4*size;
  tsize5=size5+2*size4;
  spp=size+2;
  ssize=spp*size;
  ssize2=spp*size2;
  ssize3=spp*size3;
  ssize4=spp*size4;
  ncorr=1000;                /* inital size of corr pointer */
  ncorrx=10;                 /* inital size of corrx pointer */

  /* precomputations for the spherical harmonic transform */
  do_vect(&seminaive_naive_tablespace, Reduced_Naive_TableSize(bw,cut_naive) +
	  Reduced_SpharmonicTableSize(bw,cut_naive));
  do_vect(&workspace,(8*bw*bw+16*bw));
  do_vect(&ddd, bw*(4*bw*bw-1)/3);

  seminaive_naive_table = 
    (double **)SemiNaive_Naive_Pml_Table(bw, cut_naive,
					 seminaive_naive_tablespace, workspace); 

  /*  precomputation of the d-matrices */            
  wigner(bw, 0.5*M_PI, ddd); 

  /* precomputation of the d^l_{n,0}(beta_j) */
  do_vect(&dd0, (2*bx+1)*bw*bw);
  legendre(bw, bx, dd0); 
   

  /* precomputations for the FFTW */  
  if((szs=(int *) malloc(5*sizeof(int)))==NULL){
    perror("Error in allocating memory 1\n");
    exit(1);
  }
  for(i=0;i<5;i++)
    *(szs+i)=size;

  plan_fftw=rfftwnd_create_plan(5, szs, FFTW_COMPLEX_TO_REAL,
				FFTW_ESTIMATE | FFTW_IN_PLACE);    

 
  time00 = get_the_time();
  printf("frm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time00, itime))));

  
  /* ==================== MAPS SAMPLING  =======================*/
 
  draw_line();
  fprintf(stderr, "frm6d> Sampling  maps \nfrm6d>\n");

  /* determine COM of high-res map */
  fprintf(stderr, "frm6d> Computing COM of high-res map...\n");
  cen_of_mass(g_phi_hi, &comx_hi, &comy_hi, &comz_hi);  

  /* determine COM of low-res map */
  fprintf(stderr, "frm6d> Computing COM of low-res  map...\n");
  cen_of_mass(g_phi_lo, &comx_lo, &comy_lo, &comz_lo);  


  /* calculate max radius & allocate memory */
  max_r=ceil(sqrt((double)(g_extx*g_extx+g_exty*g_exty+g_extz*g_extz)))+10;
  do_mat(&sph_sampl_lo,max_r,size*size);
  do_mat(&sph_sampl_hi,max_r,size*size);


  /* sample high res. map and get radii */
  /* (if Laplacian, this sampling is used just to obtain rmaxmin,
     and a second sampling will be done after the filtering) */
  fprintf(stderr, "frm6d> Sampling high resolution map...\n");
  volumetric6d_sampling(g_phi_hi, comx_hi, comy_hi, comz_hi, sph_sampl_hi, bw);
  rmaxmin(sph_sampl_hi, bw, rmaxmin_hi);

  /* sample low res. map and get radii */
  /* (if Laplacian, this sampling is used just to obtain rmaxmin,
     and a second sampling will be done after the filtering) */
  fprintf(stderr, "frm6d> Sampling low resolution map...\n");
  volumetric6d_sampling(g_phi_lo, comx_lo, comy_lo, comz_lo, sph_sampl_lo, bw);
  rmaxmin(sph_sampl_lo, bw, rmaxmin_lo);


  /* if Laplacian correlation requested, filter maps with Laplacian kernel,
     and do the sampling again */
  printf("frm6d> Applying filters to target and probe maps...\n");
  
  if(g_corr_mode==1) {

    relax_laplacian(&g_phi_lo, g_extx, g_exty, g_extz, g_ignored, 5.0); 
    convolve_kernel_inside_erode (&g_phi_lo, g_phi_lo, g_extx, g_exty, g_extz, g_phi_la, g_ext_la);
    
    convolve_kernel_inside (&g_phi_hi, g_phi_hi, g_extx, g_exty, g_extz, g_phi_la, g_ext_la);
   

    /* sample filtered high res. map */
    fprintf(stderr, "frm6d> Sampling filtered high resolution map...\n");
    volumetric6d_sampling(g_phi_hi, comx_hi, comy_hi, comz_hi, sph_sampl_hi, bw);

    /* sample filtered low res. map */
    fprintf(stderr, "frm6d> Sampling filtered low resolution map...\n");
    volumetric6d_sampling(g_phi_lo, comx_lo, comy_lo, comz_lo, sph_sampl_lo, bw);

    /* increase or decrease the rmaxmin's by 1 */
    rmaxmin_lo[0] += 1; rmaxmin_lo[1] += 1;
    rmaxmin_lo[2] -= 1; if(rmaxmin_lo[2]<0) rmaxmin_lo[2]=0;
    rmaxmin_lo[3] -= 1; if(rmaxmin_lo[3]<0) rmaxmin_lo[3]=0;
    rmaxmin_lo[4] += 1; rmaxmin_lo[5] += 1;
    rmaxmin_hi[0] += 1; rmaxmin_hi[1] += 1;
    rmaxmin_hi[2] -= 1; if(rmaxmin_hi[2]<0) rmaxmin_hi[2]=0;
    rmaxmin_hi[3] -= 1; if(rmaxmin_hi[3]<0) rmaxmin_hi[3]=0;
    rmaxmin_hi[4] += 1; rmaxmin_hi[5] += 1;
  }

  free(g_phi_hi);
  free(g_phi_lo);

  time01=get_the_time();
  printf("frm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time01, time00))));

  draw_line();

  stepsize=min(rmaxmin_hi[1], rmaxmin_lo[1])/((double) bw);

  /* step size for the numerical computation of the 2-center integrals: */
  if(stepsize<1.0)
    paramf[0]=1.0;
  else
    paramf[0]=stepsize;

  printf("frm6d> Step size for integration= %.3f pixels\n", paramf[0]);

  paramf[1]=(rmaxmin_lo[1]/paramf[0])*(rmaxmin_lo[1]/paramf[0]);
  paramf[2]=(rmaxmin_hi[1]/paramf[0])*(rmaxmin_hi[1]/paramf[0]);
  parami[0]=tlb;


  /* calculate max and min rho to be scanned */
  rhomax=rmaxmin_lo[1]-rmaxmin_hi[0];

  if(rmaxmin_lo[2]>0){
    if(rmaxmin_hi[2]>0){
      rhomin=max(rmaxmin_lo[2]-rmaxmin_hi[3], 0);
      rhoini=0.0;
    }
    else{
      rhomin=rmaxmin_lo[2]+rmaxmin_hi[4];
      rhoini=0.5*(rhomin+rhomax);
    }
  }
  else{
    if(g_corr_mode==1){
      rhomin=max(rmaxmin_lo[4]-rmaxmin_hi[1], 0);
      rhoini=max(0.5*(rmaxmin_lo[1]+rmaxmin_lo[0]-rmaxmin_hi[1]-rmaxmin_hi[0]), 0);
    }
    if(g_corr_mode==0){
      rhomin=0;
      rhoini=0.0;
    }
  }

  /* the 'bar' variables are in units of the step size */
  rhobarini=floor(rhoini/paramf[0]);
  rhobarmin=max(ceil((rhomin-2*stepsize)/paramf[0]),0); 
  rhobarmax=floor((rhomax+2*stepsize)/paramf[0]);

  printf("frm6d> rmaxmin_hi=%d %d %d %d %d %d\n",rmaxmin_hi[0], rmaxmin_hi[1], rmaxmin_hi[2], rmaxmin_hi[3], rmaxmin_hi[4], rmaxmin_hi[5]);
  printf("frm6d> rmaxmin_lo=%d %d %d %d %d %d\n",rmaxmin_lo[0], rmaxmin_lo[1], rmaxmin_lo[2], rmaxmin_lo[3], rmaxmin_lo[4], rmaxmin_lo[5]);
  printf("frm6d> rhomin = %.3f pixels,  rhomax = %.3f pixels\n", rhobarmin*paramf[0], rhobarmax*paramf[0]);

  if(g_rho_mode==0){
    if((rmaxmin_lo[1]>2*rmaxmin_hi[1] || rmaxmin_lo[1]>2*rmaxmin_lo[0]) &&
       rmaxmin_hi[0]>0.8*rmaxmin_lo[0]){
      g_rho_mode=2;       /* scan the whole range, from rhomin to rhomax */
    }
    else{
      g_rho_mode=1;       /* search around rhoini */
    }
  }

  if(g_rho_mode==2)
    printf("frm6d>\nfrm6d> Scan whole rho range\n");
  else
    printf("frm6d>\nfrm6d> Search around rho_initial\n");


  /* ==================== FRM COMPUTATIONS =======================*/ 
  draw_line();
  printf("frm6d> FRM 6D computations\nfrm6d>\n");
  do_vect(&idata,size*size);
  do_mat(&coeffr_lo,max_r,bw*bw);
  do_mat(&coeffi_lo,max_r,bw*bw); 
  do_mat(&coeffr_hi,max_r,bw*bw);
  do_mat(&coeffi_hi,max_r,bw*bw); 
  

  fprintf(stderr, "frm6d> Computing spherical coefficients...\n");
 
  /* get coefficients for high res. map */ 
  for(r=0;r<=rmaxmin_hi[1];r++){                      
    FST_semi_memo(sph_sampl_hi[r], idata,
		  coeffr_hi[r], coeffi_hi[r],
		  size, seminaive_naive_table, workspace,
		  1, cut_naive);
    /* correction factor for m=0 coefficients */
    for(l=0;l<bw;l++)
      *(coeffr_hi[r]+l) *= g_sqt2;
  }
 
  /* get coefficients for low res. map */ 
  for(r=1;r<=rmaxmin_lo[1];r++){                      
    FST_semi_memo(sph_sampl_lo[r], idata,
		  coeffr_lo[r], coeffi_lo[r],
		  size, seminaive_naive_table, workspace,
		  1, cut_naive);
    /* correction factor for m=0 coefficients */
    for(l=0;l<bw;l++)
      *(coeffr_lo[r]+l) *= g_sqt2;
  }


  /* allocate memory */ 
  ncorr=10000;ncorrx=1000;
  for(i=0;i<7;i++){
    corr[i]  = (double *) malloc((ncorr+1) *sizeof(double));
    corrx[i] = (double *) malloc((ncorrx+1)*sizeof(double));
    if( (corr[i] == NULL) || (corrx[i] == NULL)){
      perror("Error in allocating memory corr-x\n");
      exit(1);
    }
  }
 
  inco = (long *) malloc((ncorr+1)*sizeof(long));
  incox = (long *) malloc((ncorrx+1)*sizeof(long));
  if((inco == NULL)||(incox == NULL)){
    perror("Error in allocating memory inco-x\n");
    exit(1);
  }
  
  if ((coef_corr = (float *) malloc(tsize5*sizeof(float)))==NULL){
    perror("Error in allocating memory coef_corr\n");
    exit(1);
  }
 

  do_vect(&dtd,bw);
  do_vect(&dtdp,bw);

  inits = (long *) malloc(bw*sizeof(long));
  for(l=0;l<bw;l++)
    *(inits+l)=l*(4*l*l-1)/3;

  time02=get_the_time();
  printf("frm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time02, time01))));

  draw_line();


  /****** MAIN LOOP ******/
  fprintf(stderr, "frm6d> Entering main loop...\n");
  npeakx=0; maxcor_all=-1.0E99;
  
  for(rhobar=0;rhobar<1000;rhobar++)
    coco[rhobar]=0;

  if(g_rho_mode==1)
    rhobar=rhobarini;
  else
    rhobar=rhobarmin;

  do{
    rho=rhobar*paramf[0];
    parami[1]= ceil(max(-rmaxmin_lo[1], rho-rmaxmin_hi[1])/paramf[0]);  /* lower limit for z-sum  */
    parami[2]=floor(min( rmaxmin_lo[1], rho+rmaxmin_hi[1])/paramf[0]);  /* upper limit for z-sum  */
    paramf[3]=rhobar;
    paramf[4]=rhobar*rhobar;
    
    fprintf(stderr, "frm6d>\nfrm6d> rho=%.3f pixels\n", rho);
    fprintf(stderr, "frm6d> Computing Fourier coefficients of the correlation function...\n");

    for(index=0;index<tsize5;index++)
      *(coef_corr+index)=0.0;

    for(l=0;l<bw;l++){   /* compute T(n,h,m,hp,mp;rho) for h>=0 and hp>=0 */
      offset=l*(l+1);
      factor=(l<<1)+1;
      for(lp=0;lp<bw;lp++){
        offsetp=lp*(lp+1);
        factorp=(lp<<1)+1;
        sqll=sqrt((l+0.5)*(lp+0.5));
        isign=-1;
        for(n=0;n<=l && n<=lp;n++){
          isign=-isign;                /* (-1)^n */
          offset2=(n+l)*factor+l;
          offset4=(-n+lp)*factorp+lp;
          parami[3]=offset+n;
          parami[4]=offsetp+n;
          for(m=-l;m<=l;m++){
            index1=(n<<1)+prime(m,size)*ssize2;
            parami[5]=seanindex(m,l,bw);
            offset3=m+l;
            for(h=0;h<=l;h++)
              *(dtd+h) = ddd[*(inits+l)+offset2+h] * ddd[*(inits+l)+(h+l)*factor+offset3];
            for(mp=-lp;mp<=lp;mp++){
              index2=index1+prime(mp,size)*spp;
              parami[6]=seanindex(mp,lp,bw);
              two_cent_int(bw, bx, l, lp, m, n, mp, paramf, parami, dd0,
                           coeffr_hi, coeffi_hi, coeffr_lo, coeffi_lo,
                           &integralr, &integrali);     /* 2-center integral */
              integralr *= sqll;       /* multiply returned values by sqrt((l+0.5)*(lp+0.5)) */
              integrali *= sqll;       /*(to get actual integral (except for a factor h^3 (disregarded))) */
              if(isign<0){             /* multiply integral by (-1)^n */
                integralr = -integralr;
                integrali = -integrali;
              }
              offset5=mp+lp;
              for(hp=0;hp<=lp;hp++)
                *(dtdp+hp) = ddd[*(inits+lp)+offset4+hp] * ddd[*(inits+lp)+(hp+lp)*factorp+offset5];
              for(h=0;h<=l;h++){
                index3=index2+h*ssize3;
                for(hp=0;hp<=lp;hp++){
                  index=index3+hp*ssize;
                  tempd = *(dtd+h) * *(dtdp+hp);
                  *(coef_corr+index)   += tempd * integralr;
                  *(coef_corr+index+1) += tempd * integrali;
                }   /*  hp */
              }     /*  h  */
            }       /*  mp */
          }         /*  m  */
        }           /*  n  */
      }             /*  lp */
    }               /*  l  */

    isign=-1;
    for(n=0;n<bw;n++){          /* fill in the h<0, hp>=0 part of T(n,h,m,hp,mp;rho) */
      isign=-isign;             /* (-1)^n */
      isigm=1;
      for(m=-bw+1;m<bw;m++){
        isigm=-isigm;           /* (-1)^m */
        index1=2*n+prime(m,size)*ssize2;
        for(mp=-bw+1;mp<bw;mp++){
          index2=index1+prime(mp,size)*spp;
          for(h=-bw+1;h<0;h++){
            itemp=h*ssize3;
            index3 =index2+itemp+ssize4;
            index3b=index2-itemp;
            for(hp=0;hp<bw;hp++){
              itemp=hp*ssize;
              index =index3 +itemp;
              indexb=index3b+itemp;
              if(isign == isigm){           /* multiply by (-1)^{n+m} */
                *(coef_corr+index)   = *(coef_corr+indexb);
                *(coef_corr+index+1) = *(coef_corr+indexb+1);
              }
              else{
                *(coef_corr+index)   = -*(coef_corr+indexb);
                *(coef_corr+index+1) = -*(coef_corr+indexb+1);
              }
            }   /*  hp */
          }     /*  h  */
        }       /*  mp */
      }         /*  m  */
    }           /*  n  */

   
    isign=-1;
    for(n=0;n<bw;n++){            /* fill in the hp<0 part of T(n,h,m,hp,mp;rho) */
      isign=-isign;               /* (-1)^n */
      for(m=-bw+1;m<bw;m++){
        index1=2*n+prime(m,size)*ssize2;
        isigm=1;
        for(mp=-bw+1;mp<bw;mp++){
          isigm=-isigm;           /* (-1)^mp */
          index2=index1+prime(mp,size)*spp;
          for(h=-bw+1;h<bw;h++){
            index3=index2+prime(h,size)*ssize3;
            for(hp=-bw+1;hp<0;hp++){
              itemp=hp*ssize;
              index =index3+itemp+ssize2;
              indexb=index3-itemp;
              if(isign == isigm){           /* multiply by (-1)^{n+mp} */
                *(coef_corr+index)   = *(coef_corr+indexb);
                *(coef_corr+index+1) = *(coef_corr+indexb+1);
              }
              else{
                *(coef_corr+index)   = -*(coef_corr+indexb);
                *(coef_corr+index+1) = -*(coef_corr+indexb+1);
              }
            }   /*  hp */
          }     /*  h  */
        }       /*  mp */
      }         /*  m  */
    }           /*  n  */

    /* compute inverse FFT using FFTW */ 
    fprintf(stderr, "frm6d> Computing inverse FFT...\n");
    rfftwnd_one_complex_to_real(plan_fftw,  (fftw_complex *) coef_corr, NULL);


    /*=========================== DETECT AND SORT PEAKS  ============================*/ 

    /* if peak at rhobar+-1 is smaller than new
       correlation value at same angles, delete that peak, i.e. corrx[][i] */
    i=0;
    while(i<npeakx){
      if(coef_corr[incox[i]]>corrx[0][i] && fabs(corrx[6][i]-rhobar)<1.5){
        for(j=i+1;j<npeakx;j++){
          corrx[0][j-1]=corrx[0][j]; corrx[1][j-1]=corrx[1][j];
          corrx[2][j-1]=corrx[2][j]; corrx[3][j-1]=corrx[3][j];
          corrx[4][j-1]=corrx[4][j]; corrx[5][j-1]=corrx[5][j];
          corrx[6][j-1]=corrx[6][j]; incox[j-1]=incox[j];
        }
        npeakx--;
      }
      else
        i++;
    }

    coco[rhobar]=1;       /* this indicates that the correlation function for this rhobar
                             has been computed */

    fprintf(stderr,"frm6d> Finding peaks of the correlation function...\n");

    /* begin finding peaks of the correlation function for the current value of rho */
    /* store them in the array corr */

    /* first pass: compute max and average of the peaks */
    /* j ranges up to bw because eta<=M_PI */
    /* jp ranges up to bw because etap<=M_PI */

    npeaks=0; maxpeak=coef_corr[0]; avepeak=0.0;

    for(i=0;i<size;i++)             
      for(j=0;j<=bw;j++){           
        ind1=i+j*ssize3;
        for(k=0;k<size;k++){
          ind2=ind1+k*ssize2;
          for(jp=0;jp<=bw;jp++){    
            ind3=ind2+jp*ssize;
            for(kp=0;kp<size;kp++){
              ind4=ind3+kp*spp;

              if((f02000=coef_corr[ind4+ssize3])<=(f00000=coef_corr[ind4]))
                if ((f00020=coef_corr[ind4+ssize])<=f00000) {
                  if(i==0) f10000=coef_corr[ind4+size-1];
                  else f10000=coef_corr[ind4-1];
                  if (f10000<=f00000) {
                    if (i==size-1) f20000=coef_corr[ind4-size+1];
                    else  f20000=coef_corr[ind4+1];
                    if (f20000<=f00000) {
                      if(j==0) f01000=coef_corr[ind4+ssize4-ssize3];
                      else f01000=coef_corr[ind4-ssize3];
                      if (f01000<=f00000) {
                        if(k==0) f00100=coef_corr[ind4+ssize3-ssize2];
                        else f00100=coef_corr[ind4-ssize2];
                        if (f00100<=f00000) {
                          if(k==size-1) f00200=coef_corr[ind4-ssize3+ssize2];
                          else f00200=coef_corr[ind4+ssize2];
                          if (f00200<=f00000) {
                            if(jp==0) f00010=coef_corr[ind4+ssize2-ssize];
                            else f00010=coef_corr[ind4-ssize];
                            if (f00010<=f00000) {
                              if(kp==0) f00001=coef_corr[ind4+ssize-spp];
                              else f00001=coef_corr[ind4-spp];
                              if (f00001<=f00000) {
                                if(kp==size-1) f00002=coef_corr[ind4-ssize+spp];
                                else f00002=coef_corr[ind4+spp];
                                if (f00002<=f00000) {             
                                  npeaks++; avepeak+=f00000;
                                  if(f00000>maxpeak) maxpeak=f00000;
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
  
            }
          }
        }
      }

    maxcor[rhobar]=maxpeak;                      /* maximum for current rho */
    if(npeaks>0) avepeak /= npeaks;              /* average of all peaks at this rho */
    if(maxpeak>maxcor_all){
      maxcor_all=maxpeak;                        /* max peak so far */
      rhobar0=rhobar;
    }
    lpeak = 0.8*maxcor_all;                      /* cutoff level for the peaks */
    printf("frm6d> max peak for this rho: %f,  cutoff level: %f\n", maxpeak, lpeak);

    /* second pass: store peaks >lpeak in corr */
    /* j ranges up to bw because eta<=M_PI */
    /* jp ranges up to bw because etap<=M_PI */

    npeaks=0;    

    for(i=0;i<size;i++)            
      for(j=0;j<=bw;j++){
        ind1=i+j*ssize3;
        for(k=0;k<size;k++){
          ind2=ind1+k*ssize2;
          for(jp=0;jp<=bw;jp++){
            ind3=ind2+jp*ssize;
            for(kp=0;kp<size;kp++){
              ind4=ind3+kp*spp;
              
              if ((f00000=coef_corr[ind4])>lpeak)
                if ((f02000=coef_corr[ind4+ssize3])<=f00000)
                  if ((f00020=coef_corr[ind4+ssize])<=f00000) {
                    if(i==0) f10000=coef_corr[ind4+size-1];
                    else f10000=coef_corr[ind4-1];
                    if (f10000<=f00000) {
                      if (i==size-1) f20000=coef_corr[ind4-size+1];
                      else  f20000=coef_corr[ind4+1];
                      if (f20000<=f00000) {
                        if(j==0) f01000=coef_corr[ind4+ssize4-ssize3];
                        else f01000=coef_corr[ind4-ssize3];
                        if (f01000<=f00000) {
                          if(k==0) f00100=coef_corr[ind4+ssize3-ssize2];
                          else f00100=coef_corr[ind4-ssize2];
                          if (f00100<=f00000) {
                            if(k==size-1) f00200=coef_corr[ind4-ssize3+ssize2];
                            else f00200=coef_corr[ind4+ssize2];
                            if (f00200<=f00000) {
                              if(jp==0) f00010=coef_corr[ind4+ssize2-ssize];
                              else f00010=coef_corr[ind4-ssize];
                              if (f00010<=f00000) {
                                if(kp==0) f00001=coef_corr[ind4+ssize-spp];
                                else f00001=coef_corr[ind4-spp];
                                if (f00001<=f00000) {
                                  if(kp==size-1) f00002=coef_corr[ind4-ssize+spp];
                                  else f00002=coef_corr[ind4+spp];
                                  if (f00002<=f00000) {

                                    if(npeaks>ncorr){
                                      ncorr++;
                                      for(l=0;l<7;l++)
					corr[l]=(double *) realloc(corr[l],(ncorr+1)*sizeof(double));
                                      inco=(long *) realloc(inco,(ncorr+1)*sizeof(long));
                                    }

                                    sigma= i*M_PI/bw;
                                    eta  = j*M_PI/bw;
                                    ohm  = k*M_PI/bw;
                                    etap =jp*M_PI/bw;
                                    ohmp =kp*M_PI/bw;

                                    temp=f00000;

                                    psi =M_PI/2.0-sigma; psi=psi-2*M_PI*floor(psi/(2*M_PI));
                                    th  =M_PI-eta;
                                    phi =M_PI-ohm; phi=phi-2*M_PI*floor(phi/(2*M_PI));

                                    psip=M_PI/2.0;
                                    thp =M_PI-etap;
                                    phip=M_PI-ohmp; phip=phip-2*M_PI*floor(phip/(2*M_PI));
  
                                    corr[0][npeaks]=temp;
                                    corr[1][npeaks]=psi;
                                    corr[2][npeaks]=th;
                                    corr[3][npeaks]=phi;
                                    corr[4][npeaks]=thp;
                                    corr[5][npeaks]=phip;
                                    corr[6][npeaks]=rhobar;
                                    inco[npeaks]=ind4;
                                    npeaks++;

                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }

            }
          }
        }
      }
    fprintf(stderr, "frm6d> %d peaks detected\n", npeaks );

    /* sort corr in decreasing order */
    for(i=0;i<npeaks;i++)         
      for(j=i+1;j<npeaks;j++)
        if(corr[0][i]<corr[0][j]){
          for(k=0;k<7;k++)  
            SWAPPING(corr[k][j], corr[k][i], double);
          SWAPPING(inco[j], inco[i], long);   
        }

    /* eliminate close-by solutions */
    i=0; 
    while(i<npeaks){
      j=i+1;
      while(j<npeaks){
        if(dist_mot(debe, i, corr, j, corr)<dist_cut){
          for(k=j+1;k<npeaks;k++){
            for(l=0;l<7;l++)
              corr[l][k-1]=corr[l][k];
            inco[k-1]=inco[k];
	  }
	  npeaks--;
        }
        else
          j++;
      }
      i++;
    }

    fprintf(stderr,"frm6d> %d distance-filtered peaks\n", npeaks);

    /* save value of npeakx because it will be updated in this loop */
    npeakxa=npeakx;

    /* append new peaks whose distance to all bigger ones at rhobar+-1 is large */
    for(i=0;i<npeaks;i++){
      /* if flag1 stays =0 in the j-loop, append corr[][i] to corrx */
      flag1=0;
      for(j=0;j<npeakxa;j++){
        if(fabs(corrx[6][j]-rhobar)<1.5 && corrx[0][j]>corr[0][i] && dist_mot(debe, i, corr, j, corrx)<dist_cut){
          /* in this case corr[][i] is not appended to corrx */
          flag1=1;           
          break;
        }
      }
      if(flag1==0){
        if(npeakx>ncorrx){
          ncorrx++;
          for(l=0;l<7;l++){
            corrx[l]=(double *) realloc(corrx[l],(ncorrx+1)*sizeof(double));
          }
          incox=(long *) realloc(incox,(ncorrx+1)*sizeof(long));
        }
        for(k=0;k<7;k++)
          corrx[k][npeakx]=corr[k][i];
        incox[npeakx]=inco[i];
        npeakx++;
      }
    }

    fprintf(stderr,"frm6d> %d cumulative peaks\n", npeakx);

    /* test for terminating rho loop: */

    if(g_rho_mode==1){
      if(rhobar0>rhobarmin && coco[rhobar0-1]==0){
        rhobar=rhobar0-1;          /* compute for this rhobar */
        flag2=0;
      }
      else if(rhobar0<rhobarmax && coco[rhobar0+1]==0){
        rhobar=rhobar0+1;          /* compute for this rhobar */
        flag2=0;
      }
      else
        flag2=1;                   /* rhobar0-1 and rhobar0+1 are done... */
      /* ...so check to the right and left of rhobar0 for 90%: */
      if(flag2==1){                /* check to the right */
        if(rhobar0<rhobarmax) flag2=0;
        for(rhobar=rhobar0+1;rhobar<=rhobarmax;rhobar++){
          if(coco[rhobar]==1){
            if(maxcor[rhobar]<0.9*maxcor_all || rhobar==rhobarmax){
              flag2=1;             /* OK to the right of rhobar0 */
              break;
            }
          }
          else
            break;
        }
      }
      if(flag2==1){                /* check to the left */
        if(rhobar0>rhobarmin) flag2=0;
        for(rhobar=rhobar0-1;rhobar>=rhobarmin;rhobar--){
          if(coco[rhobar]==1){
            if(maxcor[rhobar]<0.9*maxcor_all || rhobar==rhobarmin){
              flag2=1;             /* OK to the left of rhobar0 */
              break;
            }
          }
          else
            break;
        }
      }
    }
    else{
      if(rhobar<rhobarmax){
        rhobar++;
        flag2=0;
      }
      else
        flag2=1;
    }
  } while(flag2==0);             /****** end of rho loop ******/

  fprintf(stderr, "frm6d>\nfrm6d> Exiting main loop...\n");

  rfftwnd_destroy_plan(plan_fftw);

  time03=get_the_time();
  printf("frm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time03, time02))));

  draw_line();

  /* free space */
  free(workspace); free(ddd); free(dd0);
  free(seminaive_naive_tablespace);
  for(i=0;i<7;i++) free(corr[i]);

  for(i=0;i<max_r;i++){ 
    free(sph_sampl_hi[i]); 
    free(sph_sampl_lo[i]); 
    free(coeffr_hi[i]); free(coeffi_hi[i]); 
    free(coeffr_lo[i]); free(coeffi_lo[i]); 
  } 
  free(idata);
  free(coef_corr);

  fprintf(stderr, "frm6d> Generating list of peaks...\n");

  /* delete peaks smaller than lpeak */ 
  i=0;
  while(i<npeakx-1){
    if(corrx[0][i]<lpeak){        
      for(j=i+1;j<npeakx;j++){
        for(l=0;l<7;l++)
          corrx[l][j-1]=corrx[l][j];
        incox[j-1]=incox[j];
      }
      npeakx--;
    }
    else
      i++;
  }

  /* sort corrx in increasing order (to do the next step properly) */
  for(i=0;i<npeakx;i++)
    for(j=i+1;j<npeakx;j++)
      if(corrx[0][i]>corrx[0][j]){
        for(k=0;k<7;k++)
          SWAPPING(corrx[k][i], corrx[k][j], double);
        SWAPPING(incox[i], incox[j], double);
      }

  /* check for relative maxima in the rho direction */
  /* and delete those points that are not */
  i=0;                   
  while(i<npeakx-1){
    j=i+1;
    while(j<npeakx){
      if(fabs(corrx[1][j]-corrx[1][i])+fabs(corrx[2][j]-corrx[2][i])+
         fabs(corrx[3][j]-corrx[3][i])+fabs(corrx[4][j]-corrx[4][i])+
         fabs(corrx[5][j]-corrx[5][i])<0.5*debe &&
         fabs(corrx[6][j]-corrx[6][i])<1.2){    /* in this case delete corrx[][i] */
        for(k=i+1;k<npeakx;k++){
          for(l=0;l<7;l++)
            corrx[l][k-1]=corrx[l][k];
          incox[k-1]=incox[k];
        }
        npeakx--;
        j=i+1;
      }
      else
        j++;
    }
    i++;
  }

  /* sort corrx in decreasing order */
  for(i=0;i<npeakx;i++)
    for(j=i+1;j<npeakx;j++)
      if(corrx[0][i]<corrx[0][j]){
        for(k=0;k<7;k++)
          SWAPPING(corrx[k][j], corrx[k][i], double);
        SWAPPING(incox[j], incox[i], long);
      }

  /* eliminate close-by solutions */
  i=0; 
  while(i<npeakx){
    j=i+1;
    while(j<npeakx){
      if(dist_mot(debe, i, corrx, j, corrx)<dist_cut){
	for(k=j+1;k<npeakx;k++){
          for(l=0;l<7;l++)
            corrx[l][k-1]=corrx[l][k];
          incox[k-1]=incox[k];
	}
	npeakx--;
      }
      else
        j++;
    }
    i++;
  }


  /* show results */
  printf("frm6d>\nfrm6d>        Psi    Theta    Phi       X       Y       Z      Correlation\nfrm6d>\n");
  for(i=0; i<npeakx && i<20; i++){
    euler_comp(M_PI-corrx[3][i], corrx[2][i], M_PI-corrx[1][i], M_PI/2.0, corrx[4][i], corrx[5][i],
               &psic, &thc, &phic);
    vx= (corrx[6][i]*paramf[0])*sin(corrx[2][i])*sin(corrx[3][i]);
    vy=-(corrx[6][i]*paramf[0])*sin(corrx[2][i])*cos(corrx[3][i]);
    vz= (corrx[6][i]*paramf[0])*cos(corrx[2][i]);
    printf("frm6d> %3d  %6.2f  %6.2f  %6.2f   %6.2f  %6.2f  %6.2f  %14.6E\n",
           i+1, psic*180/M_PI, thc*180/M_PI, phic*180/M_PI, vx, vy, vz, corrx[0][i]);
  }

  time04=get_the_time();
  printf("frm6d>\nfrm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time04, time03))));

  draw_line();


  /*=========================== SAVE BEST RESULTS ============================*/ 
 
  /* create dummy structure g_pdb_move */
  g_pdb_move = (PDB *) malloc(g_num_atoms*sizeof(PDB));
  for(i=0;i<g_num_atoms;i++)
    g_pdb_move[i]=g_pdb_original[i];

  if (npeakx < g_num_explored) g_num_explored=npeakx;

  fprintf(stderr, "frm6d> Saving the %d best results of %d peaks found...\n", g_num_explored, npeakx);


  for(j=0;j<g_num_explored;j++){

    sprintf(out_string, "frm6_best%ld.pdb",j+1);

    g_outfile = fopen(out_string, "w");
    if (g_outfile == NULL) {
      fprintf(stderr, "frm6d> Error: can't open file! %s [e.c. 80320]\n",out_string); 
      exit(80320);
    }
    fprintf(g_outfile, "REMARK\nREMARK    File name %s\n",out_string);
    fprintf(g_outfile, "REMARK    FRM6D low-resolution fit of structure %s into map %s\n",argv[2],argv[1]);
    fprintf(g_outfile, "%s\n",corr_string);
    fprintf(g_outfile, "REMARK    Resolution anisotropy factor (Z vs. XY): %f\n",g_target_ani);
    fprintf(g_outfile, "REMARK    Resolution: %f Angstrom; density cutoff: %f\n",g_target_res,g_low_cutoff);
    fprintf(g_outfile, "REMARK    Normalized correlation coefficient: %f\n",
	    corrx[0][j]);
    fprintf(g_outfile, "REMARK    \n");
    fclose(g_outfile);

    euler_comp(M_PI-corrx[3][j], corrx[2][j], M_PI-corrx[1][j], M_PI/2.0, corrx[4][j], corrx[5][j],
               &psic, &thc, &phic);
    vx= (corrx[6][j]*paramf[0])*sin(corrx[2][j])*sin(corrx[3][j]);
    vy=-(corrx[6][j]*paramf[0])*sin(corrx[2][j])*cos(corrx[3][j]);
    vz= (corrx[6][j]*paramf[0])*cos(corrx[2][j]);


    fprintf(g_outfile, "REMARK    Euler angles (Psi, Theta, Phi):  %7.3f %7.3f %7.3f deg.\n",psic, thc, phic);
    fprintf(g_outfile, "REMARK    New COM of PDB (X,Y,Z):          %7.3f %7.3f %7.3f A\n",  (vx+comx_lo)*g_width+g_gridx,
	    (vy+comy_lo)*g_width+g_gridy,
	    (vz+comz_lo)*g_width*g_target_ani+g_gridz);
    fprintf(g_outfile, "REMARK    \n");


    /* shift to com */
    translate(g_pdb_original, g_pdb_move, g_num_atoms,
              -(comx_hi-0.5*g_extx)*g_width,
              -(comy_hi-0.5*g_exty)*g_width,
              -(comz_hi-0.5*g_extz)*g_width*g_target_ani);

    /* rotate pdb wrt its com=(0,0,0) */
    rot_euler(g_pdb_move, g_pdb_move, g_num_atoms, psic, thc, phic);

    /* translate pdb to its final position */
    translate(g_pdb_move, g_pdb_move, g_num_atoms,
              (vx+comx_lo)*g_width+g_gridx,
              (vy+comy_lo)*g_width+g_gridy,
              (vz+comz_lo)*g_width*g_target_ani+g_gridz);

    /* save pdb */   
    appendpdb(out_string, g_num_atoms, g_pdb_move);
  }

  time05=get_the_time();
  printf("frm6d>   Lap time: %s\n", smart_sprint_time(time_to_sec(time_diff(time05, time04))));

  free(g_pdb_original); free(g_pdb_move);

  etime = get_the_time();
  printf("frm6d>\nfrm6d> Total time:  %s\nfrm6d>\n", 
         smart_sprint_time(time_to_sec(time_diff(etime, itime))));

  return 0;
}


/*====================================================================*/
void volumetric6d_sampling(double *phi, double cx, double cy, double cz, double *sample[], int bw){
/* performs the sampling in spherical coordinates */

  int i, j, r, size;
  double lambda, theta, initlambda;          /* longitude and latitude */
  double stsl, stcl, st, ct;
  double x, y, z;
  double a, b, c;
  int x1, y1, z1;
  double maxext;
  double f111, f112, f121, f122, f211, f212, f221, f222;
  unsigned long index;

  size=(bw<<1);
  maxext=ceil(sqrt((double)(g_extx*g_extx+g_exty*g_exty+g_extz*g_extz)));

  initlambda=M_PI/bw;
  /* i represents latitude values  */
  /* j represents longitude values */
  for(i=0;i<size;i++){           
    theta=latitude(i, bw);
    st=sin(theta); ct=cos(theta);
    index=i*size;
    for(j=0, lambda=0;j<size;j++, index++, lambda+=initlambda) {        
      stsl=st*sin(lambda); stcl=st*cos(lambda);
      for(r=0, x = cx, y=cy, z=cz;r<=maxext;r++, x+=stcl, y+=stsl, z+=ct){
        if(x<0.0 || x>=g_extx-1.0 || y<0.0 || y>=g_exty-1.0 || z<0.0 || z>=g_extz-1.0){
          for(;r<=maxext;r++)
            *(sample[r]+index)=0.0;
          break;
        }
        x1=floor(x); y1=floor(y); z1=floor(z);
        a=x-x1; b=y-y1; c=z-z1;
        f111=*(phi+gidz_general(z1, y1, x1, g_exty, g_extx));
        f112=*(phi+gidz_general(z1+1, y1, x1,g_exty,g_extx));
        f121=*(phi+gidz_general(z1, y1+1, x1,g_exty,g_extx));
        f122=*(phi+gidz_general(z1+1, y1+1, x1,g_exty,g_extx));
        f211=*(phi+gidz_general(z1, y1, x1+1,g_exty,g_extx));
        f212=*(phi+gidz_general(z1+1, y1, x1+1,g_exty,g_extx));
        f221=*(phi+gidz_general(z1, y1+1, x1+1,g_exty,g_extx));
        f222=*(phi+gidz_general(z1+1, y1+1, x1+1,g_exty,g_extx));
        *(sample[r]+index)=(1-a)*((1-b)*((1-c)*f111+c*f112)+b*((1-c)*f121+c*f122))+
	  a*((1-b)*((1-c)*f211+c*f212)+b*((1-c)*f221+c*f222));
      }
    }
  }
}


void rmaxmin(double *sample[], int bw, int rmm[6]){
/* computes various radii of the sampled map */

  int i, j, r, size;
  double maxext;
  double f111, f112, f121, f122, f211, f212, f221, f222;
  unsigned long index;

  size=(bw<<1);
  maxext=ceil(sqrt((double)(g_extx*g_extx+g_exty*g_exty+g_extz*g_extz)));

  rmm[0]=maxext;      /* rmin,e: minimum distance from COM to starred hull of map */
  rmm[1]=0;           /* rmax,e: maximum distance from COM to starred hull of map */

  for(i=0;i<size;i++){
    index=i*size;
    for(j=0;j<size;j++,index++)
      for(r=maxext;r>=0;r--)
        if (*(sample[r]+index) != 0.0) {
          if(r<rmm[0]) rmm[0]=r;
          if(r>rmm[1]) rmm[1]=r;
          break;
        }
  }

  rmm[2]=maxext;      /* r'min,i: minimum distance from COM to points of the map with nonzero density */
  rmm[3]=0;           /* r'max,i: maximum distance from COM to points of the map with nonzero density */

  for(i=0;i<size;i++){
    index=i*size;
    for(j=0;j<size;j++,index++)
      for(r=0;r<=maxext;r++)
        if(*(sample[r]+index) != 0.0){
          if(r<rmm[2]) rmm[2]=r;
          if(r>rmm[3]) rmm[3]=r;
          break;
        }
  }

  rmm[4]=maxext;        /* rmin,i: minimum distance from COM to the first point of zero density */
  rmm[5]=0;             /* rmax,i: maximum distance from COM to the first point of zero density */

  for(i=0;i<size;i++){
    index=i*size;
    for(j=0;j<size;j++,index++)
      for(r=0;r<=maxext;r++)
        if(*(sample[r]+index)==0.0){
          if(r-1<rmm[4]) rmm[4]=r-1;
          if(r-1>rmm[5]) rmm[5]=r-1;
          break;
        }
  }
  if(rmm[4]<0) rmm[4]=0;
}


/*====================================================================*/
void legendre(int bw, int bx, double *dd0){
  /* computes the d-functions d^l_{n,0}(theta_j), 
     for l=0,...,bw-1; abs(n)<=l; j=0,...,2*bx. 
     theta_j is defined by: cos(theta_j)=j/bx-1 */
  /* they are returned in the (3D) prism array dd0 */
  /* Reference: T. Risbo, Journal of Geodesy (1996) 70:383-396 */

  double p, q, pc, qc, theta, temp, *d;
  int size, index, i, j, k, l, l2;
  unsigned long init1, init2, bw2;
  double fact1, fact2;

  size=2*bw; bw2=bw*bw;

  d = (double *) malloc(3*size*sizeof(double));
 
  if( d == NULL ){
    perror("Error in allocating memory");
    exit(1);
  }

  for(j=0;j<=2*bx;j++){                 /* loop over argument theta */
    init1=j*bw2;
    theta=acos(((double) j)/bx-1.0);
    p=sin(theta/2.); q=cos(theta/2.);   /* Cayley-Klein parameters */
    pc=p; qc=q;

    *(dd0+init1)=1.0;                   /* d^0_{0,0}(theta) */
  
    /* d-functions of degree 1 */
    *(d+0  +1) =  g_sqt2*q*p;
    *(d+1*3+1) = -p*pc+q*qc;
    *(d+2*3+1) = -g_sqt2*pc*qc;

    /* copy the d^1 to the proper location in dd0 */
    *(dd0+init1+1) = *(d+1);
    *(dd0+init1+2) = *(d+4);
    *(dd0+init1+3) = *(d+7);
    
    for(l=2;l<bw;l++){                /* l = degree of d-functions to be computed next */
      l2=(l<<1);
      for(i=0;i<l2;i++){
        *(d+i*3+0)=0.0;
        *(d+i*3+2)=0.0;
      }
      for(i=0;i<l2-1;i++){            /* balanced Edmonds scheme */
        index=i*3;
        temp=*(d+index+1);
        fact1=sqrt((l2-1-i)/(double)l)*temp;
        fact2=sqrt((i+1)   /(double)l)*temp;
        *(d+index  ) +=  fact1*q;
        *(d+index+3) += -fact2*pc;
        *(d+index+2) +=  fact1*p;
        *(d+index+5) +=  fact2*qc;
      }
      for(i=0;i<=l2;i++)
        *(d+i*3+1)=0.0;
      for(i=0;i<l2;i++){             /* binomial recursion */
        index=i*3;
        fact1=sqrt((l2-i)*l)/l2;
        fact2=sqrt((i+1) *l)/l2;
        *(d+index+1) +=  fact1* *(d+index+2)*q;
        *(d+index+4) += -fact2* *(d+index+2)*pc;
        *(d+index+1) +=  fact1* *(d+index  )*p;
        *(d+index+4) +=  fact2* *(d+index  )*qc;
      }
      init2=l*l+init1;
      for(i=0;i<=l2;i++)           /* copy d to the proper location in dd0 */
        *(dd0+init2+i) = *(d+i*3+1);
    }
  }
  free(d);
}


/*====================================================================*/
void two_cent_int(int bw, int bx, int l, int lp, int m, int n, int mp,
                  double paramf[], long parami[], double *dd0,
                  double *coeffr_hi[], double *coeffi_hi[],
                  double *coeffr_lo[], double *coeffi_lo[],
                  double *tcir, double *tcii){
  /* computes the 2-center integral I^{l,lp}_{m,n,mp}(rho), except for the factor
     h^3 * sqrt((l+0.5)*(lp+0.5)) (h^3 is never considered because it is constant throughout) */
  /* Note: rho is rhobar*h, where rhobar comes in paramf[3] and h in paramf[0] */
  /* the real and imag parts of the integral are returned in tcir and tcii */

  int r, rp, i, j, sbar, zbar, f1;
  long sbar2, rbar2, rpbar2;
  double Sr, Si;            /* real & imag parts of zbar-sum */
  double fact;
  double rbar, rpbar;
  double fr, fi, gr, gi;

  *tcir=0.0; *tcii=0.0;
  for(sbar=1;sbar<=bw;sbar++){
    sbar2=sbar*sbar;
    Sr=0.0; Si=0.0; f1=0;
    for(zbar=parami[1];zbar<=parami[2];zbar++){
      rbar2=sbar2+zbar*zbar;
      if(rbar2>paramf[1])
        if(f1==0) continue; else break;
      rpbar2=rbar2-paramf[3]*(zbar<<1)+paramf[4];
      if(rpbar2>paramf[2])
        if(f1==0) continue; else break;
      f1=1;  /* current point is in the intersection of both circles */
      rbar=sqrt(rbar2); rpbar=sqrt(rpbar2);
      r =floor(paramf[0]*rbar +0.5);
      rp=floor(paramf[0]*rpbar+0.5);
      i=floor((zbar/rbar+1)*bx+0.5);
      j=floor(((zbar-paramf[3])/rpbar+1)*bx+0.5);
      fr = *(coeffr_lo[r] +parami[5]); fi = *(coeffi_lo[r] +parami[5]);
      gr = *(coeffr_hi[rp]+parami[6]); gi = *(coeffi_hi[rp]+parami[6]);
      fact = *(dd0+parami[3]+(i<<parami[0])) * *(dd0+parami[4]+(j<<parami[0]));
      Sr +=  (fr*gr - fi*gi)*fact;
      Si += -(fi*gr + fr*gi)*fact;
    }
    if(f1==0) break;
    *tcir += Sr*sbar;
    *tcii += Si*sbar;
  }
}


/*====================================================================*/
void euler_comp(double psi1, double th1, double phi1, double psi2, double th2, double phi2,
                double *psi, double *th, double *phi) {
  /* computes the Euler angles of the composition R=R1*R2 of the rotations given by 
     R1=(phi1,th1,psi1) and R2=(phi2,th2,psi2)  */
 
  int i, j, k;
  double matrix1[3][3], matrix2[3][3], rot[3][3], temp;  

  get_rot_matrix(matrix1,psi1,th1,phi1);
  get_rot_matrix(matrix2,psi2,th2,phi2);
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      temp=0.0;
      for(k=0;k<3;k++)
        temp += matrix1[i][k]*matrix2[k][j];
      rot[i][j]=temp;
    }

  if(rot[2][2]>=1.0)
    *th=0.0;
  else if(rot[2][2]<=-1.0)
    *th=M_PI;
  else
    *th=acos(rot[2][2]);

  if(fabs(rot[0][2])<0.001 && fabs(rot[1][2])<0.001){    
    temp=atan2(rot[1][0], rot[1][1]);
    *psi=M_PI/2.0;
    *phi=3.0*M_PI/2.0-temp; *phi -= 2*M_PI*floor(*phi/(2*M_PI)); 
  }
  else{
    temp=atan2(rot[1][2], rot[0][2]);
    *psi=M_PI/2.0-temp; *psi -= 2*M_PI*floor(*psi/(2*M_PI));
    temp=atan2(rot[2][1], -rot[2][0]);
    *phi=3.0*M_PI/2.0-temp; *phi -= 2*M_PI*floor(*phi/(2*M_PI));
  }
}


/*====================================================================*/
double dist_mot(double debe, long i, double *corr[], long j, double *corrx[]){
  /* returns the distance between the motions given by the Euler angles and
     rhos in corr[1..6][i] and corrx[1..6][j] */

  int k, l, t;
  double product[3][3], rot[3][3], tm[3][3], temp, trace;
  double sin_phi, cos_phi, sin_th, cos_th, sin_psi, cos_psi;
  double angle, length, c_th1p, s_th1p, c_th2p, s_th2p;

  get_rot_matrix(product,corrx[1][j],corrx[2][j],corrx[3][j]); /* R2 */
  get_rot_matrix(rot,corr[1][i],corr[2][i],corr[3][i]);        /* R1 */
  
  for(k=0;k<3;k++)             /* tm = R2 * transpose(R1) */
    for(l=0;l<3;l++){
      tm[k][l]=0.0;
      for(t=0;t<3;t++)
        tm[k][l] += product[k][t]*rot[l][t];
    }
  for(k=0;k<3;k++)             /* move tm to product */
    for(l=0;l<3;l++)
      product[k][l]=tm[k][l]; 

  get_rot_matrix(rot, M_PI/2.0, corr[4][i], corr[5][i]);  /* R1p */

  for(k=0;k<3;k++)              /* tm = R2 * transpose(R1) * R1p */
    for(l=0;l<3;l++){
      tm[k][l]=0.0;
      for(t=0;t<3;t++)
        tm[k][l] += product[k][t]*rot[t][l];
    }
  for(k=0;k<3;k++)              /* move tm to product */
    for(l=0;l<3;l++)
      product[k][l]=tm[k][l];
  
  get_rot_matrix(rot, M_PI/2.0, corrx[4][j], corrx[5][j]);  /* R2p */
 
  
  trace=0.0;     /* trace(R2 * transpose(R1) * R1p * transpose(R2p)) */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      trace += product[k][l]*rot[k][l];

  temp=0.5*(trace-1.0);
  if(temp>=1.0)
    angle=0.0;
  else if(temp<=-1.0)
    angle=M_PI;
  else
    angle=acos(temp);

  length=corr[6][i]*corr[6][i]+corrx[6][j]*corrx[6][j]-2.*corr[6][i]*corrx[6][j]*
    (cos(corr[4][i])*cos(corrx[4][i])+sin(corr[4][i])*sin(corrx[4][i])*cos(corrx[5][i]-corr[5][j]));

  if(length<0.0) length=0.0;

  return sqrt(length+angle*angle/(debe*debe));
}


/*====================================================================*/
static void cen_of_mass(double *phi, double *cx, double *cy, double *cz){
  /* Computes the center of mass of the density phi */

  unsigned long i, j, k, ind;
  double mass, temp;
  double cxap, cyap, czap;

  cxap=0.0; cyap=0.0; czap=0.0; mass=0.0;
  for (k=0;k<g_extz;k++) 
    for (j=0;j<g_exty;j++)
      for (i=0;i<g_extx;i++){
        ind=gidz_general(k, j, i, g_exty, g_extx);
        temp=fabs(*(phi+ind));
        cxap += temp*i;
        cyap += temp*j;
        czap += temp*k;
        mass += temp;
      }

  cxap /= mass; cyap /= mass; czap /= mass;

  *cx=cxap; *cy=cyap; *cz=czap;
}


/*====================================================================*/
static double longitude(int j, int bw){
  return M_PI*j/bw;
}


/*====================================================================*/
static double latitude(int i, int bw){
  return 0.5*M_PI*(i+0.5)/bw;
}


/*====================================================================*/
static void wigner(int bw, double theta, double *ddd){
  /* computes the d-functions for argument theta, for degrees 0 through bw-1 */
  /* the matrices are returned in the (3D) piramid ddd */
  /* Reference: T. Risbo, Journal of Geodesy (1996) 70:383-396 */

  double p, q, pc, qc, temp, *d, *dd;
  int size, index, i, j2, k, l, hdeg;
  unsigned long init;
  double fact1, fact2, fact3, fact4;

  size=2*bw;
  
  d  = (double *) malloc(size*size*sizeof(double));
  dd = (double *) malloc(size*size*sizeof(double));

  if( (d == NULL) || (dd == NULL) ){
    perror("Error in allocating memory");
    exit(1);
  }

  p=sin(theta/2); q=cos(theta/2);   /* Cayley-Klein parameters */
  pc=p; qc=q;

  *(ddd+0)=1.0;                     /* d-matrix of degree 0 */
  
  /* d-matrix of degree 1/2 */
  *(d+0)=q;
  *(d+1)=p;
  *(d+1*size+0)=-pc;
  *(d+1*size+1)=qc;
  
  for(l=1;l<bw;l++){       /* l is the degree of the next d-matrix to be saved in ddd */
    j2=(l<<1)-1;           /* j2 will be twice the degree of the d-matrix about to be computed */
    for(hdeg=0;hdeg<2;hdeg++){
      fact1=q/ ++j2;
      fact2=pc/j2;
      fact3=p/j2;
      fact4=qc/j2;
      for(i=0;i<=j2+1;i++)
        for(k=0;k<=j2+1;k++)
          *(dd+i*size+k)=0.0;
      for(i=0;i<j2;i++)
        for(k=0;k<j2;k++){
          index=i*size+k;
          temp=*(d+index);
          *(dd+index)        += sqrt((j2-i)*(j2-k))*temp*fact1;
	  *(dd+index+size)   += -sqrt((i+1)*(j2-k))*temp*fact2;
	  *(dd+index+1)      += sqrt((j2-i)*(k+1)) *temp*fact3;
	  *(dd+index+size+1) += sqrt((i+1) *(k+1)) *temp*fact4;
        }
      for(i=0;i<=j2;i++)        /* move dd to d */
        for(k=0;k<=j2;k++)
          *(d+i*size+k)=*(dd+i*size+k);
      if(hdeg==0){              /* if degree is integer, copy d to the proper location in ddd */
        init=l*(4*l*l-1)/3;
        for(i=0;i<=j2;i++)
          for(k=0;k<=j2;k++)
            *(ddd+init+i*(j2+1)+k) = *(d+i*size+k);
      }
      if(l==bw-1) break;     /* for l=bw-1 do hdeg=0 only */
    }
  }
  free(d); free(dd);
}



/*====================================================================*/
static int prime(int n, int s){
  /* returns n if it is >=0, or n+s otherwise */

  if(n>=0) return n;
  return n+s;
}


/*====================================================================*/
static void read_options(int argc, char **argv) {
  /* print usage info and read options from input arguments */
  
  char *pos;
  int i;
  char option_string[199];
  
  if(argc<2) { /* print usage and options info */

    printf("frm6d> USAGE:  frm6 <Situs density map> <PDB structure> <options>\n"); 
    printf("frm6d>\n");
    printf("frm6d> OPTIONS:\n"); 
    printf("frm6d>\n");
    printf("frm6d>\t-res <float>\t   Target resolution in A [default: -res 15]\n");  
    printf("frm6d>\t-ani <float>\t   Resolution anisotropy factor [default: -ani 1]\n"); 


    printf("frm6d>\t-cutoff <float>\t   Density map cutoff value [default: -cutoff 0.0]\n");
    printf("frm6d>\t-corr <int>\t   Correlation method [default: -corr 1]\n"); 
    printf("frm6d>\t\t\t\t-corr 0  -->  Volumetric correlation\n");
    printf("frm6d>\t\t\t\t-corr 1  -->  Laplacian correlation\n");
    printf("frm6d>\t-mode <int>\t   rho scanning mode [default: -mode 0]\n");
    printf("frm6d>\t\t\t\t-mode 0  -->  Automatically select scanning mode\n");
    printf("frm6d>\t\t\t\t-mode 1  -->  Search around rho_initial\n");
    printf("frm6d>\t\t\t\t-mode 2  -->  Scan the whole rho range\n");
    printf("frm6d>\t-explor <int>\t   Number of local maxima explored [default: -explor 5]\n");
    printf("frm6d>\n");
    draw_line();
    exit(0);
  }

  /* now read options from arguments and assign global variables */
  
  sprintf(option_string,"");
  printf("frm6d> Options read:\n");
  for(i=1;i<argc;i++) sprintf(option_string,"%s %s",option_string,argv[i]);
  
  /* target resolution */
  g_target_res=15.0;
  if((pos=(char *)strstr(option_string," -res"))!=NULL)  
    sscanf(pos+5,"%lf",&g_target_res);
  printf("frm6d> Target resolution %3.3f\n",g_target_res);

  /* resolution anisotropy */
  g_target_ani=1.0;
  if((pos=(char *)strstr(option_string," -ani"))!=NULL)  
    sscanf(pos+5,"%lf",&g_target_ani);
  if (g_target_ani<0.001 || g_target_ani>1000) {
    fprintf(stderr, "frm6d> Error: Anisotropy out of range [e.c. 80410]\n"); 
    exit(80410);   
  }
  printf("frm6d> Resolution anisotropy %3.3f\n",g_target_ani);


  /* cutoff of the low-resolution map */
  g_low_cutoff=0.0;
  if((pos=(char *)strstr(option_string," -cutoff"))!=NULL)  
    sscanf(pos+8,"%lf",&g_low_cutoff);
  printf("frm6d> Low-resolution map cutoff %3.3f\n",g_low_cutoff);

  /* correlation mode */
  g_corr_mode=1;
  if((pos=(char *)strstr(option_string," -corr"))!=NULL)  
    sscanf(pos+6,"%d",&g_corr_mode);
  switch (g_corr_mode) {
  case 0: printf("frm6d> Volumetric correlation\n"); break;
  case 1: printf("frm6d> Laplacian correlation\n"); break;
  default:
    fprintf(stderr, "frm6d> Error: Could not identify -corr option [e.c. 80410]\n"); 
    exit(80410);   
  }

  /* rho scanning mode */
  g_rho_mode=0;
  if((pos=(char *)strstr(option_string," -mode"))!=NULL)  
    sscanf(pos+6,"%d",&g_rho_mode);
  switch (g_rho_mode) {
  case 0: printf("frm6d> Automatic selection of rho scanning mode\n"); break;
  case 1: printf("frm6d> Search around rho_initial\n"); break;
  case 2: printf("frm6d> Scan the whole rho range\n"); break;
  default:
    fprintf(stderr, "frm6d> Error: Could not identify -mode option [e.c. 80411]\n"); 
    exit(80411);   
  }

  /* number of best fits */
  g_num_explored=6;
  if((pos=(char *)strstr(option_string," -explor"))!=NULL)
    sscanf(pos+8,"%d",&g_num_explored);
  printf("frm6d> Number of best fits explored %2d\n",g_num_explored);

}


/*====================================================================*/
static void get_lowmap(char *file_name) {
  /* reads low resolution map from file_name */

  printf("colores> Processing low-resolution map.\n");

  readvol(file_name, &g_width, &g_gridx, &g_gridy, &g_gridz, &g_extx, &g_exty, &g_extz, &g_phi_lo);
  g_nvox = g_extx*g_exty*g_extz;
  
  if (g_target_ani!=1.0) {
    interpolate_map (&g_phi_du, &g_extx, &g_exty, &g_extz, &g_gridx, &g_gridy, &g_gridz, 
		     g_width, g_width, g_width*g_target_ani, g_phi_lo, g_extx, g_exty, g_extz, g_gridx, 
		     g_gridy, g_gridz, g_width, g_width, g_width);  
    g_nvox = g_extx*g_exty*g_extz;
    cp_vect_destroy(&g_phi_lo,&g_phi_du,g_nvox); 
  }
  
  /* set density values below g_low_cutoff to zero */
  threshold(g_phi_lo, g_nvox, g_low_cutoff);

  /* shrink map about non-zero density and resize to odd intervals */
  shrink_margin(&g_phi_du, &g_extx, &g_exty, &g_extz, &g_gridx, &g_gridy, &g_gridz, &g_nvox, 
		g_phi_lo, g_extx, g_exty, g_extz, g_gridx, g_gridy, g_gridz, 
		g_width, g_width, g_width*g_target_ani);
  cp_vect_destroy(&g_phi_lo,&g_phi_du,g_nvox); 

  print_map_info(g_phi_lo, g_nvox);
}


/*====================================================================*/
static void get_centered_structure_and_radius (char *file_name, double *max_radius) {
  /* reads atomic structure and centers it */

  double com_x; 
  double com_y;
  double com_z;

  printf("frm6d> Processing atomic structure.\n");  

  /* read PDB file */
  read_and_weight_mass (file_name, &g_num_atoms, &g_pdb_original);
   
  calc_center_mass (g_pdb_original, g_num_atoms, &com_x, &com_y, &com_z);
  *max_radius = calc_sphere(g_pdb_original,g_num_atoms,com_x,com_y,com_z);
  printf ("frm6d> COM: %6.3f %6.3f %6.3f, radius: %6.3f Angstrom\n",com_x,com_y,com_z,*max_radius);
  
  /* center structure */ 
  translate (g_pdb_original, g_pdb_original, g_num_atoms, -com_x, -com_y, -com_z);
  

}


/*====================================================================*/
static void draw_line() {
  printf("_______________________________________________________________________________________\n\n");
}
