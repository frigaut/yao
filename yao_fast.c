/*
 * yao_fast.c
 * 
 * yao wavefront sensors engines
 * 
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2012, Francois Rigaut
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include "ydata.h"
#include "yapi.h"

#include "play.h"  // for p_wall_secs() and p_cpu_secs()

#define FFTWOPTMODE FFTW_EXHAUSTIVE
// use FFTW_PATIENT for thread optimization (see below):
//#define FFTWOPTMODE FFTW_PATIENT

// static int n_threads = 1;

int use_sincos_approx_flag = 0;

void _eclat_float(float *ar, int nx, int ny);
void _poidev(float *xmv, long n);
void _gaussdev(float *xmv, long n);

// int Y__fftw_init_threads(void)
// {
  // return fftwf_init_threads();
// }
// 
// 
// void Y_fftw_set_n_threads(int nargs)
// {
  // n_threads = ygets_i(0);
  // If you plan with FFTW_PATIENT, it will automatically
  // disable threads for sizes that don't benefit from parallelization.
  // fftwf_plan_with_nthreads(n_threads);
// }
// 
// void Y_fftw_get_n_threads(int nargs)
// {
  // ypush_int(n_threads);
// }


void _set_sincos_approx(int flag)
{
  use_sincos_approx_flag = flag;
}

void Y__get_sincos_approx(int nargs)
{
  ypush_int(use_sincos_approx_flag);
}

float sine(float x)
{
  const float pi = 3.141592653589793f;
  const float B = 4.0f/pi;
  const float C = -4.0f/(pi*pi);
  const float D = 2*pi;
  const float P = 0.225;

  x = x - roundf(x/D)*D;
  
  float y = B * x + C * x * fabsf(x);
  
  y = P * (y * fabsf(y) - y) + y;
  
  return y;
}

float cosine(float x)
{
  const float hpi = 3.141592653589793f/2.0f;

  float y = sine(x+hpi);
  
  return y;
}

void _sinecosinef(float x, float *s, float *c)
{
  float sc;
  *s = sine(x);
  *c = cosine(x);
  sc = (*c)*(*c) + (*s)*(*s);

  sc = 1.0f/sqrtf(sc);
  *s = *s * sc;
  *c = *c * sc;
}

int _import_wisdom(char *wisdom_file)
{
  FILE *fp;
  int  status;

  if((fp = fopen(wisdom_file, "r"))==NULL) return (1);

  status = 1-fftwf_import_wisdom_from_file(fp);

  fclose(fp);

  return status;
}

int _export_wisdom(char *wisdom_file)
{
  FILE *fp;

  if((fp = fopen(wisdom_file, "w"))==NULL) return (1);

  fftwf_export_wisdom_to_file(fp);

  fclose(fp);

  return (0);
}

int _init_fftw_plans(int nlimit)
{
  int n, size;
  fftwf_complex *inf,*outf;
  fftwf_plan p1, p2;
  int plan_mode;

  plan_mode = FFTWOPTMODE;

  size=1;
  for (n=0;n<nlimit+1;n++) {
    printf("Optimizing 2D FFT - size = %d\n",size);
    inf  = fftwf_malloc(size*size*sizeof(fftwf_complex));
    outf = fftwf_malloc(size*size*sizeof(fftwf_complex));
    p1 = fftwf_plan_dft_2d(size, size, inf, outf, FFTW_FORWARD, plan_mode);
    p2 = fftwf_plan_dft_2d(size, size, inf, outf, FFTW_BACKWARD, plan_mode);
    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_free(inf);
    fftwf_free(outf);

    size*=2;
  }

  size = 1;
  for (n=0;n<nlimit+1;n++) {
    printf("Optimizing 1D FFT - size = %d\n",size);
    inf  = fftwf_malloc(size*sizeof(fftwf_complex));
    outf = fftwf_malloc(size*sizeof(fftwf_complex));
    p1 = fftwf_plan_dft_1d(size, inf, outf, FFTW_FORWARD, plan_mode);
    p2 = fftwf_plan_dft_1d(size, inf, outf, FFTW_BACKWARD, plan_mode);
    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_free(inf);
    fftwf_free(outf);

    size*=2;
  }

  return(0);
}

int _init_fftw_plan(int size)
{
  fftwf_complex *inf;
  fftwf_complex *outf;
  fftwf_plan p1,p2;
  float *ptr;
  int i;
  int plan_mode;

  plan_mode = FFTWOPTMODE;

  printf("Optimizing 2D FFT - size = %d\n",size);
  inf  = fftwf_malloc(size*size*sizeof(fftwf_complex));
  outf = fftwf_malloc(size*size*sizeof(fftwf_complex));

  ptr = (void *)inf;
  for (i=0;i<2*size*size;i++) *(ptr++) = 0.0f;

  p1 = fftwf_plan_dft_2d(size, size, inf, outf, FFTW_FORWARD, plan_mode);
  p2 = fftwf_plan_dft_2d(size, size, inf, outf, FFTW_BACKWARD, plan_mode);

  fftwf_destroy_plan(p1);
  fftwf_destroy_plan(p2);
  fftwf_free(inf);
  fftwf_free(outf);

  return(0);
}

/**************************************************************
 * The following function computes the PSF, given a pupil and *
 * a phase, both float. phase can be a 3 dimensional array,   *
 * in which case the routine compute and returns an image for *
 * each plan. This routine uses the FFTW routines.            *
 * Called by calc_psf_fast(pupil,phase,scale=) in yorickVE.i      *
 * Author: F.Rigaut                                           *
 * Date: Dec 12 2003.                                         *
 * modified 2004apr3 to adapt FFTW from apple veclib          *
 **************************************************************/

int _calc_psf_fast(float *pupil, /* pupil image, dim [ 2^n2 , 2^n2 ] */
                   float *phase, /* phase cube,  dim [ 2^n2 , 2^n2 , nplans ] */
                   float *image, /* output images, dim [ 2^n2 , 2^n2 , nplans ] */
                   int n2,       /* log2 of side linear dimension (e.g. 6 for 64x64 arrays) */
                   int nplans,   /* number of plans (min 1, mostly to spped up calculations */
                                 /* by avoiding multiple setups and mallocs */
                   float scal)   /* phase scaling factor */
{
  /* Declarations */

  fftwf_complex *in, *out;
  float         *ptr;
  long          i,k,koff,n;
  fftwf_plan     p;
  float         ppsin,ppcos;
      
  /* Set the size of 2d FFT. */
  n = 1 << n2;

  // fftwf_plan_with_nthreads(n_threads);

  /* Allocate memory for the input operands and check its availability. */
  in  = fftwf_malloc(sizeof(fftwf_complex) * n * n);
  out = fftwf_malloc(sizeof(fftwf_complex) * n * n );

  if ( in == NULL || out == NULL ) { return (-1); }
  
  /* Set up the required memory for the FFT routines */
  p = fftwf_plan_dft_2d(n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  /* Main loop on plan #, in case several phases are input to the routine */
  for ( k=0; k<nplans; k++ ) {

    koff = k*n*n;
    //    ptr  = &(in[0]);
    ptr  = (void *)in;

    for ( i=0; i<n*n; i++ ) {
      if ( pupil[i] != 0.f ) {
        //~ *(ptr)   = pupil[i] * cosf( phase[koff+i] * scal );
        //~ *(ptr+1) = pupil[i] * sinf( phase[koff+i] * scal );
        //~ *(ptr)   = pupil[i] * cosine( phase[koff+i] * scal );
        //~ *(ptr+1) = pupil[i] * sine( phase[koff+i] * scal );
        if (use_sincos_approx_flag) _sinecosinef(phase[koff+i] * scal, &ppsin, &ppcos);
        else sincosf(phase[koff+i] * scal, &ppsin, &ppcos);
        *(ptr)   = pupil[i] * ppcos;
        *(ptr+1) = pupil[i] * ppsin;
      } else {
        *(ptr)   = 0.0f;
        *(ptr+1) = 0.0f;
      }
      ptr +=2;
    }

    /* Carry out a Forward 2d FFT transform, check for errors. */

    fftwf_execute(p); /* repeat as needed */

    //    ptr = &(out[0]);
    ptr  = (void *)out;
    for ( i = 0; i < n*n; i++ ) {
      image[koff+i] = ( *(ptr) * *(ptr) +
                        *(ptr+1) * *(ptr+1) );
      ptr +=2;
    }
    /* swap quadrants. */
    _eclat_float(&(image[koff]),n,n);
  }

  /* Free allocated memory. */
  fftwf_destroy_plan(p);
  fftwf_free(in);
  fftwf_free(out);
  return (0);
}




int _fftVE(float *rp,
           float *ip,
           int   n,
           int   dir)      /* forward (1) or reverse (-1) */
{
  /* Declarations */

  fftwf_complex *in, *out;
  float         *ptr;
  long          i;
  fftwf_plan    p;

  /* Allocate memory for the input operands and check its availability. */
  in  = fftwf_malloc(sizeof(fftwf_complex) * n * n);
  out = fftwf_malloc(sizeof(fftwf_complex) * n * n );

  if ( in == NULL || out == NULL ) { return (-1); }
  
  /* Set up the required memory for the FFT routines */
  if (dir == 1) {
    p = fftwf_plan_dft_2d(n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    p = fftwf_plan_dft_2d(n, n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  /* fill input */
  ptr  = (void *)in;

  for ( i=0; i<n*n; i++ ) {
    *(ptr)   = rp[i];
    *(ptr+1) = ip[i];
    ptr +=2;
  }

  /* Carry out a Forward 2d FFT transform, check for errors. */

  fftwf_execute(p); /* repeat as needed */

  ptr = (void *)out;
  for ( i = 0; i < n*n; i++ ) {
    rp[i] = *(ptr);
    ip[i] = *(ptr+1);
    ptr +=2;
  }

  /* Free allocated memory. */
  fftwf_destroy_plan(p);
  fftwf_free(in);
  fftwf_free(out);
  return (0);
}



void _fftVE2(fftwf_complex *in,
             fftwf_complex *out,
             int n,       /* log2 of side linear dimension (e.g. 6 for 64x64 arrays) */
             int dir)      /* forward (1) or reverse (-1) */
/* do not use this one, it was an experiment, but _fftVE is actually faster */
{
  /* Declarations */
  fftwf_plan     p;
      
  /* Set up the required memory for the FFT routines */
  if (dir == 1) {
    p = fftwf_plan_dft_2d(n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  } else {
    p = fftwf_plan_dft_2d(n, n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  
  /* Carry out a Forward 2d FFT transform, check for errors. */
  fftwf_execute(p); /* repeat as needed */

  /* Free allocated memory. */
  fftwf_destroy_plan(p);

}

int embed_image(float *inim, // Input (origin) image
   int indx,      // X dim of origin image
   int indy,      // Y dim of origin image
   float *outim,  // Output (destination) image
   int outdx,     // X dim of destination image
   int outdy,     // Y dim of destination image
   int destx,     // 0-based X indice of origin image pixel (0,0) in dest. image
   int desty,     // 0-based Y indice of origin image pixel (0,0) in dest. image
   int roll)      // roll the input image before embedding?
{
  int i,j,ioff,joff;

  if (roll==0) {
    for ( j=0 ; j<indy ; j++ ) {
      joff = j+desty;
      if (joff<0) continue;
      if (joff>=outdy) break;
      joff *= outdx;
      for ( i=0 ; i<indx ; i++ ) {
        ioff = i+destx;
        if (ioff<0) continue;
        if (ioff>=outdx) break;
        outim[ioff+joff] += inim[i+j*indx];
      }
    }
    return 0;
  }
  
  int hx,hy;
  hx = indx/2;
  hy = indy/2;
  
  for ( j=0 ; j<hy ; j++ ) {
    joff = j+desty;
    if (joff<0) continue;
    if (joff>=outdy) break;
    joff *= outdx;
    for ( i=0 ; i<hx ; i++ ) {
      ioff = i+destx;
      if (ioff<0) continue;
      if (ioff>=outdx) break;
      outim[ioff+joff] += inim[i+hx+(j+hy)*indx];
    }
  }

  for ( j=hy ; j<indy ; j++ ) {
    joff = j+desty;
    if (joff<0) continue;
    if (joff>=outdy) break;
    joff *= outdx;
    for ( i=0 ; i<hx ; i++ ) {
      ioff = i+destx;
      if (ioff<0) continue;
      if (ioff>=outdx) break;
      outim[ioff+joff] += inim[i+hx+(j-hy)*indx];
    }
  }

  for ( j=0 ; j<hy ; j++ ) {
    joff = j+desty;
    if (joff<0) continue;
    if (joff>=outdy) break;
    joff *= outdx;
    for ( i=hx ; i<indx ; i++ ) {
      ioff = i+destx;
      if (ioff<0) continue;
      if (ioff>=outdx) break;
      outim[ioff+joff] += inim[i-hx+(j+hy)*indx];
    }
  }

  for ( j=hy ; j<indy ; j++ ) {
    joff = j+desty;
    if (joff<0) continue;
    if (joff>=outdy) break;
    joff *= outdx;
    for ( i=hx ; i<indx ; i++ ) {
      ioff = i+destx;
      if (ioff<0) continue;
      if (ioff>=outdx) break;
      outim[ioff+joff] += inim[i-hx+(j-hy)*indx];
    }
  }
  return 0;
  
}



/* Shack- Hartmann coded in C
   pass one phase array and a set of indices (start and end of each subapertures)
   then this routine puts the phase sections in one larger phase and does a serie
   of FFTs to compute the images. Then each image is binned according to the map
   "binindices". It then put these images back into a large image
   that contains all subaperture image, suitable for visualization.
   flow:
   - set up of a couple of constants and memory allocation
   - loop on subaperture number:
   - put pupil and phase in complex A[ns,ns]
   - do FFT
   - compute image as square(FFT result complex) (no eclat)
   - compute binned/resampled subaperture image using binindices (has 
   to made up for the missing eclat at previous step).
   - compute X and Y centroids
   - optionaly stuff this image in a larger image (fimage) for display
   - free memory

   The array used in here are:

   pupil (dim*dim): pupil image, float
   phase (dim*dim): phase image, float

   virtual(nx*ny) : wavefront-let, i.e. part of the pupil/phase contained
                    within the subaperture. 

   A      (ns*ns) : Complex array into which the amplitude and phase 
                    corresponding to  one subaperture wavefront-let are placed. 
                    Dimension n2*n2 e.g. 16x16 if 5x5 pixels/subaperture
   simage (ns*ns) : Focal plane image obtained after FFT of A.

   ximage (nx*nx) : Extended focal plane array, into which simage is embedded 
                    (+ shifted) to allow extended field of view.

   bimage (nb*nb) : Image obtained after binning of ximage to take into account 
                    bigger pixels than the one delivered after the FFT. 

   fimage (nf*nf) : Final image into which all the subaperture images have been 
                    placed at pre-defined positions (see imistart)
*/
int _shwfs_phase2spots(
   float *pupil,        // input pupil
   float *phase,        // input phase, in microns
   float phasescale,    // phase scaling factor: microns -> radians @ wfs.lambda
   float *phaseoffset,  // input phase offset
   int   dim,           // X or Y dim of phase. Used as stride for extraction
   
   int   *istart,       // vector of i starts of each subaperture
   int   *jstart,       // vector of j starts of each subaperture
   int   nsx,           // subaperture i size
   int   nsy,           // subaperture j size
   int   nsubs,         // # of subapertures
   
   int   sdimpow2,      // dimension of small array for FFTs, in power of 2
   
   long  domask,        // go through amplitude mask loop (0/1).
   float *submask,      // subaperture mask. Corresponds/applied to simage.
   
   float *kernel,       // to convolve the (s)image with. dim 2^sdimpow2
                        // compute dynamically at each call, i.e. can change 
   float *kernels,      // to convolve the (s)image with, one per subaperture
                        // dimension: 2^sdimpow2 * 2 * nsubs. FFTs precomputed
                        // at init call.
   float *kerfftr,      // real part of kernels FFT. dim: same as kernels
   float *kerffti,      // imaginary part of kernels FFT. same dim as kerfftr
   int   initkernels,   // init kernels: pre-compute FFTs
   int   kernconv,      // convolve with kernel?
   
   int   *binindices,   // int 2d array containing the indices of the binned 
                        // subaperture image, i.e. to which pixel in the binned
                        // image should one pixel in the FFT'd image be added.
   int   nb,            // dimension of binned subaperture image (bimage)
   int   rebinfactor,   // rebin factor from small to big pixels
   int   nx,            // dimension of extended subaperture image
   float *unittip,      // nsx*nsy array with unit tip (1,2,3..)
   float *unittilt,     // nsx*nsy array with unit tilt
   float *lgs_prof_amp, // vector of lgs profile amplitudes
   float *lgs_defocuses,// vector of lgs profile altitudes
   int   n_in_profile,  // dimension of lgs_prof_amp and lgs_prof_alt
   float *unit_defocus,  // as it says, same dimsof as phase
   float *fimage,       // final image with spots
   int   *svipc_subok,  // to skip (0) subap for svipc partial spot comput.
   int   *imistart,     // vector of i starts of each image
   int   *imjstart,     // vector of j starts of each image
   int   fimnx,         // final image X dimension
   int   fimny,         // final image Y dimension
   
   float *flux,         // vector of flux (input), dim nsubs
   float *rayleighflux, // vector of flux for rayleigh (input), dim nsubs
   float *skyflux,      // vector of flux for sky (input), dim nsubs
   float darkcurrent,   // dark current, e-/pix/frame
   
   int   rayleighflag,  // enable rayleigh processing
   float *rayleigh,     // rayleigh background, nx*nx*nsubs
                        // here I separated background and rayleigh. the background includes
                        // not only the rayleigh, but also the sky and the dark current.
   int   bckgrdinit,    // init background processing. fill bckgrdcalib
   
   int   counter,       // current counter (in number of cycles)
   int   niter)         // total # of cycles over which to integrate
           
{
  /* Declarations */

  fftwf_complex *A, *result, *Ax, *Kx, *resultx;
  fftwf_plan    fftps,fftpx,fftpxi;
  float         *ptr,*ptr1,*ptr2;
  float         *phase_scaled;
  float         *simage;
  float         *bimage;
  float         *ximage; // extended image for one subap, un-rebinned
  float         *bsubmask;
  float         *brayleigh;
  float         tot, totrayleigh, krp, kip, sky;
  float         corfact;
  float         dx,dxp,dy,dyp;
  float         lgsdef,lgsamp,pp,ppsin,ppcos;
  long          log2nr, log2nc, n, ns;
  int           i,j,k,l,koff,kk,nalt;
  int           idxp,idyp,ndx,ndy,dynrange;
  int           debug=0;
  int           nxdiff;
  double        sys,cpu0,cpu1,cpu2,cpu3,cpu4,cpu5,cpu6;
  double        cpu10,cpu21,cpu32,cpu43,cpu54,cpu65;
  const float   pi = 3.141592653589793f;
  const float   twopi = 2*pi;
  
  cpu10=0.0;cpu21=0.0;cpu32=0.0;cpu43=0.0;cpu54=0.0;cpu65=0.0;
  
  //======================
  // Global setup for FFTs:
  //======================

  // fftwf_plan_with_nthreads((int)1);

  // Set the size of 2d FFT.
  log2nr = sdimpow2; 
  log2nc = sdimpow2;
  n      = 1 << ( log2nr + log2nc ); // total number of pixels in small array
  ns     = 1 << log2nr;              // side of small array, pixels, n=ns*ns

  // enlarge dynamical range?
  if (nx==ns) dynrange=0; else dynrange=1;

  // integrate = 1; // force to pass by the end (subap overlap upgrade)
  // if (niter > 1) {integrate = 1;}  // we are in "integrating mode"

  if (debug>1) printf("here1\n");

  // Allocate memory for the input operands and check its availability.
  A            = fftwf_malloc ( n * sizeof ( fftwf_complex ) );
  result       = fftwf_malloc ( n * sizeof ( fftwf_complex ) );
  Ax           = fftwf_malloc ( nx *nx * sizeof ( fftwf_complex ) );
  Kx           = fftwf_malloc ( nx *nx * sizeof ( fftwf_complex ) );
  resultx      = fftwf_malloc ( nx *nx * sizeof ( fftwf_complex ) );
  phase_scaled = ( float* ) malloc ( dim * dim * sizeof ( float ) );
  simage       = ( float* ) malloc ( ns * ns * sizeof ( float ) );
  ximage       = ( float* ) malloc ( nx * nx * sizeof ( float ) );
  bimage       = ( float* ) malloc ( nb * nb * sizeof ( float ) );
  brayleigh    = ( float* ) malloc ( nb * nb * sizeof ( float ) );
  bsubmask     = ( float* ) malloc ( nb * nb * sizeof ( float ) );
  
  
  if ( A == NULL || result == NULL || Ax == NULL || \
       Kx == NULL || resultx == NULL || phase_scaled == NULL || 
       bimage == NULL || simage == NULL || ximage == NULL || 
       brayleigh == NULL || bsubmask == NULL ) { return (1); }

  //Zero out final image if first iteration
  if (counter == 1) {
    // in fact, let's put the dark current value in there:
    //~ for (i=0;i<(fimnx*fimny);i++) { fimage[i] = 0.0f; }
    // these 2 lines temporarily disabled for svipc shm_var of ffimage:
    // FIXME (remove comments):FIXME FIXME FIXME
    //~ totdark = (float) niter * darkcurrent;
    //~ for (i=0;i<(fimnx*fimny);i++) { fimage[i] = totdark; }
    // 2012sep17: the problem above it seems is that the loop is over the
    // entire array. If several forks are accessing it, -> problem.
    // the simplest might be to add it at the yorick level before
    // calling _shwfs_spots2slopes().
  }

  // compute scaled phase ( we do this operation several times below):
  for ( i=0 ; i<dim*dim ; i++ ) \
     phase_scaled[i] = (phase[i]+ phaseoffset[i]) * phasescale;

  // compute bsubmask from submask and binindices:
  for ( i=0 ; i<nb*nb ; i++ ) { bsubmask[i] = (float)(1-domask); }
  if (domask) {
    for ( i=0 ; i<nx*nx ; i++ ) {
      if (binindices[i] >= 0) {
        bsubmask[binindices[i]] += (float)submask[i];
      }
    }
    corfact = 1.0f / (float)rebinfactor / (float)rebinfactor;
    for ( i=0 ; i<nb*nb ; i++ ) { bsubmask[i] *= corfact; }
  }
  
  // in the following, we'll need to shift slightly where we embed simage
  // into ximage.
  // this is linked to the number of -1 pixels at the end of binindices.
  nxdiff = 0;
  // do that for first row only:
  for (i=0;i<nx;i++) if (binindices[i]==-1) nxdiff++;
  
  // Set up the required memory for the FFT routines and 
  // check its availability.
  fftpx = fftwf_plan_dft_2d(nx, nx, Ax, Kx, FFTW_FORWARD, FFTW_ESTIMATE);

  if (initkernels == 1) {
    // Transform kernels, store and return for future use
    for ( l=0 ; l<nsubs ; l++ ) {
      ptr = (void *)Ax;
      for ( i=0 ; i<nx*nx ; i++ ) {
        *(ptr) = kernels[i+l*nx*nx];
        *(ptr+1) = 0.0f;
        ptr += 2;
      }
      fftwf_execute(fftpx);

      ptr = (void *)Kx;
      for ( i=0 ; i<nx*nx ; i++ ) {
        kerfftr[i+l*nx*nx] = *(ptr);
        kerffti[i+l*nx*nx] = *(ptr+1);
        ptr += 2;
      }
    }
  }

  if (debug>1) printf("here2\n");

  if (kernconv == 1) {
    // Transform kernel
    ptr = (void *)Ax;
    for ( i = 0; i < nx*nx; i++ ) {
      *(ptr)   = kernel[i];
      *(ptr+1) = 0.0f;
      ptr +=2;
    }
    fftwf_execute(fftpx);
  }
  fftwf_destroy_plan(fftpx);

  fftps  = fftwf_plan_dft_2d(ns, ns, A, result, FFTW_FORWARD, FFTW_ESTIMATE);
  fftpx  = fftwf_plan_dft_2d(nx, nx, Ax, resultx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftpxi = fftwf_plan_dft_2d(nx, nx, Ax, resultx, FFTW_BACKWARD, FFTW_ESTIMATE);

  //=====================
  // LOOP ON SUBAPERTURES
  //=====================
  for ( l=0 ; l<nsubs ; l++ ) {

    // zero out ximage:
    for ( i=0 ; i<nx*nx ; i++ ) ximage[i] = 0.0f;
    
    //====================
    // LOOP ON LGS PROFILE    
    //====================
    for ( nalt=0; nalt<n_in_profile ; nalt++ ) {
      
      lgsdef = lgs_defocuses[nalt];
      lgsamp = lgs_prof_amp[nalt];
      if ( n_in_profile==1 ) lgsamp = 1.0f;
  
      if ( svipc_subok[l]==0 ) continue;
      
      cpu0  = p_cpu_secs(&sys);

      // reset A and result
      ptr = (void *)A;
      for ( i=0; i<n ; i++ ) { *(ptr)   = 0.0f; *(ptr+1) = 0.0f; ptr +=2; }
      
      ptr = (void *)result;
      for ( i=0; i<n ; i++ ) { *(ptr)   = 0.0f; *(ptr+1) = 0.0f; ptr +=2; }

      // indice offset of phaselet in phase/pupil array
      koff = istart[l] + jstart[l]*dim;
  
      // START section to allow larger dynamical range by
      // subtracting a tilt to the phase and moving later on
      // the image in the big image
      // (declarations on top of function)
      if (dynrange) {
        // compute approximate average slope over the subaperture
        // it doesn't matter if this is not the exact value
        // as it only serves to determine if we should offset,
        // but the end result should be the same. Because of edge subapertures,
        // we'll have to do the whole average gradient calculation
        dx = 0.0f; dy = 0.0f;
        ndx = 0; ndy = 0;
        for ( j=0; j<(nsy-1); j++ ) {
          for ( i=0; i<(nsx-1) ; i++ ) {
            k = koff + i + j*dim;
            if ( pupil[k] && pupil[k+1] ) {
              dx +=            phase_scaled[k+1] - phase_scaled[k] + \
                      lgsdef*( unit_defocus[k+1] - unit_defocus[k] );
              ndx++;
            }
            if ( pupil[k] && pupil[k+dim] ) {
              dy +=            phase_scaled[k+dim] - phase_scaled[k] + \
                      lgsdef*( unit_defocus[k+dim] - unit_defocus[k] );
              ndy++;
            }
          }
        }
        if (ndx) dx = dx / (float)(ndx); // in radian/pixel
        if (ndy) dy = dy / (float)(ndy);
        // now we need to transform this average phase gradient into
        // pixel motion.
        // So how many small pixels we expect the spot to move?
        // if dx = 2*pi, the spot will wrap all the way back to center,
        // i.e. will move by ns pixels:
        dxp = dx / twopi * (float)(ns);
        dyp = dy / twopi * (float)(ns);
        // round to the nearest (small) pixel: 
        idxp = lroundf(dxp);
        idyp = lroundf(dyp);
        if (idxp!=0) dx = dx * (float)(idxp) / dxp; else dx = 0.0f;
        if (idyp!=0) dy = dy * (float)(idyp) / dyp; else dy = 0.0f;
         if ((debug>1)&&((idxp!=0)||(idyp!=0))) printf("idxp = %d, idyp = %d, dx=%f, dy=%f\n",idxp,idyp,dx,dy);
      } else {
        dx = 0.0f; dy = 0.0f; 
        idxp = 0; idyp = 0;
      }
      // END section to allow larger dynamical range
      // (more below to add to phase and to shift imagelets)
  
      cpu1  = p_cpu_secs(&sys);
      cpu10 += cpu1-cpu0;

      // fill in the complex wavefront array for this subaperture:
      // cos & sin are very costly, so we use an aproximation:
      ptr = (void *)A;
      if (dynrange) {
        
        if (debug>1) if (l==10) printf("here, dynrange enabled\n");
        for ( j=0; j<nsy ; j++ ) {
          for ( i=0; i<nsx ; i++ ) {
            k = koff + i + j*dim;
            kk = i + j*nsx;
            pp = phase_scaled[k] + lgsdef * unit_defocus[k] \
                     - dx * unittip[kk] - dy * unittilt[kk];
            if (use_sincos_approx_flag) _sinecosinef(pp,&ppsin,&ppcos);
            else sincosf(pp,&ppsin,&ppcos);
            *(ptr + 2*(i+j*ns))   = pupil[k] * ppcos;
            *(ptr + 2*(i+j*ns)+1) = pupil[k] * ppsin;
          }
        }

      } else {

        if (debug>1) if (l==10) printf("here, dynrange disabled\n");
        for ( j=0; j<nsy ; j++ ) {
          for ( i=0; i<nsx ; i++ ) {
            k = koff + i + j*dim;
            pp = phase_scaled[k];
            if (use_sincos_approx_flag) _sinecosinef(pp,&ppsin,&ppcos);
            else sincosf(pp,&ppsin,&ppcos);
            *(ptr + 2*(i+j*ns))   = pupil[k] * ppcos;
            *(ptr + 2*(i+j*ns)+1) = pupil[k] * ppsin;
          }
        }

      }

      if (debug>1) printf("here3 ");

      cpu2  = p_cpu_secs(&sys);
      cpu21 += cpu2-cpu1;

      // Carry out a Forward 2d FFT transform, check for errors.
      fftwf_execute(fftps);
  
      // compute image from complex image object:
      ptr = (void *)result;
      for ( i=0; i<n; i++ ) {
        simage[i] = (*(ptr) * *(ptr) + *(ptr+1) * *(ptr+1) );
        simage[i] *= lgsamp;
        if ( (debug>10) && ( l==10 ) ) printf("%f ",simage[i]);
        ptr +=2;
      }
      
      if ( (debug>1) && ( l==10 ) ) printf(" ");
      
      cpu3  = p_cpu_secs(&sys);
      cpu32 += cpu3-cpu2;

      // Embed (and add to) this simage into ximage, the extended field 
      // of view image for this subaperture (with shifts computed above):
      if ( (debug>1) && (l==10) ) printf("\nns=%d nx=%d\n",(int)ns,(int)nx);
      embed_image(simage,ns,ns,ximage,nx,nx,(nx-ns)/2+idxp-nxdiff,(nx-ns)/2+idyp-nxdiff,1);
    
    } // END LOOP ON LGS PROFILE
    
    if ( (debug>1) && (l==10) ) printf("here4\n");

    if ((debug>10)&&(l==10)) {
      FILE *fp;
      fp=fopen("xim-pre.dat", "w");
      fprintf(fp, "%d\n",(int)nx);
      for ( i=0 ; i<nx*nx ; i++ ) fprintf(fp, "%f\n",ximage[i]);
      fclose(fp);
    }
    
    // Carry out convolution by kernel if required
    if (kernconv == 1) {
      // Transform ximage
      ptr = (void *)Ax;
      for ( i=0 ; i<nx*nx ; i++ ) {
        *(ptr)   = ximage[i]; 
        *(ptr+1) = 0.0f;
        ptr +=2;
      }
      fftwf_execute(fftpx);
  
      // multiply by kernel transform:
      ptr  = (void *)Kx;
      ptr1 = (void *)Ax;
      ptr2 = (void *)resultx;
      for ( i=0 ; i<nx*nx ; i++ ) {
        // this is FFT(kernel) * FFT(kernelS)
        krp = *(ptr)*kerfftr[i+l*nx*nx]- *(ptr+1)*kerffti[i+l*nx*nx];
        kip = *(ptr)*kerffti[i+l*nx*nx]+ *(ptr+1)*kerfftr[i+l*nx*nx];
        // and next we multiply by FFT(image):
        *(ptr1)   = *(ptr2)*krp   - *(ptr2+1)*kip;
        *(ptr1+1) = *(ptr2+1)*krp + *(ptr2)*kip;
        ptr +=2; ptr1 +=2; ptr2 +=2;
      }
      // Transform back:
      fftwf_execute(fftpxi);
  
      ptr = (void *)resultx;
      for ( i=0 ; i<nx*nx ; i++ ) {
        ximage[i] = ( *(ptr) * *(ptr) + *(ptr+1) * *(ptr+1) );
        ptr +=2;
      }
    }
  
    cpu4  = p_cpu_secs(&sys);
    cpu43 += cpu4-cpu3;

    // FLUX NORMALIZATION FOR STAR. Has to be done *before* applying fieldstop
    // will be used a bit below
    // LGS FIXME FIXME FIXME: flux totally screwed up w/ new lgs_prof_amp!
    tot = 0.0f;
    for ( i=0 ; i<nx*nx ; i++ ) tot += ximage[i];
      
    // APPLY FIELD STOP / AMPLITUDE MASK
    // For instance to take into account the central dark spot of STRAP,
    // or more generally a field stop
    if (domask == 1) {
      for ( i=0 ; i<nx*nx ; i++ ) {
        ximage[i] = ximage[i] * submask[i];
      }
    }
  
    // IF BACKGROUND CALIBRATION, NULL STAR SIGNAL
    if (bckgrdinit) {
      for ( i=0 ; i<nx*nx ; i++ ) ximage[i] = 0.0f;
    }
  
    // PUT THIS SUBAPERTURE'S XIMAGE INTO BIMAGE (binned image)
    for ( i=0 ; i<nb*nb ; i++ ) { bimage[i] = 0.0f; }
  
    for ( i=0 ; i<nx*nx ; i++ ) {
      if (binindices[i]<0) continue;
      //~ if (binindices[i]>nb*nb) {
        //~ printf("binindices[%d] = %d > nb*nb (%d). It shouldn't be so! (ns=%d, nx=%d, nb=%d)\n",i,binindices[i],nb*nb,ns,nx,nb);
        //~ return 1;
      //~ }
      bimage[binindices[i]] += ximage[i];
      //~ if (ximage[i]==0) printf("ximage[%d]==0 ",i);
    }
  
    if ((debug>10)&&(l==10)) {
      FILE *fp;
      fp=fopen("sim.dat", "w");
      fprintf(fp, "%d\n",(int)ns);
      for ( i=0 ; i<ns*ns ; i++ ) fprintf(fp, "%f\n",simage[i]);
      fclose(fp);
      fp=fopen("xim.dat", "w");
      fprintf(fp, "%d\n",nx);
      for ( i=0 ; i<nx*nx ; i++ ) fprintf(fp, "%f\n",ximage[i]);
      fclose(fp);
      fp=fopen("bim.dat", "w");
      fprintf(fp, "%d\n",nb);
      for ( i=0 ; i<nb*nb ; i++ ) fprintf(fp, "%f\n",bimage[i]);
      fclose(fp);
      fp=fopen("kre.dat", "w");
      fprintf(fp, "%d\n",nx);
      ptr  = (void *)Kx;
      for ( i=0 ; i<nx*nx ; i++ ) { fprintf(fp, "%f\n",*ptr); ptr+=2; }
      fclose(fp);
      fp=fopen("kim.dat", "w");
      fprintf(fp, "%d\n",nx);
      ptr  = (void *)Kx; ptr++;
      for ( i=0 ; i<nx*nx ; i++ ) { fprintf(fp, "%f\n",*ptr); ptr+=2; }
      fclose(fp);
    }
  
    // NORMALIZE FLUX FOR STAR
    if (tot>0.0f) {
      tot = flux[l]/tot;
      for ( i=0 ; i<nb*nb ; i++ ) { bimage[i] = bimage[i]*tot; }
    }
    
    // COMPUTE RAYLEIGH BACKGROUND
    // I have to do this after the convolution because otherwise there is a lot
    // of wrapping/ringing. I could do it once the image is binned. Saved for future
    // upgrade (will save a bit of time).
    for ( i=0 ; i<nb*nb ; i++ ) { brayleigh[i] = 0.0f; }
    if (rayleighflag==1) {
      for ( i=0 ; i<nx*nx ; i++ ) {
        if (binindices[i]>=0) {
          brayleigh[binindices[i]] += rayleigh[i+l*nx*nx];
        }
      }
      // NORMALIZE FLUX FOR RAYLEIGH
      if (debug) printf("here4-1\n");
      totrayleigh = 0.0f;
      for ( i=0 ; i<nb*nb ; i++ ) { totrayleigh += brayleigh[i]; }
      if (debug) printf("l=%d, totrayleigh=%f\n",l,totrayleigh);
  
      if (totrayleigh > 0.0f) {
        totrayleigh = rayleighflux[l]/totrayleigh;
        for ( i = 0; i < nb*nb; i++ ) brayleigh[i] = brayleigh[i]*totrayleigh;
      }
      if (debug) printf("here4-2\n");
    }
  
    // NORMALIZE FLUX FOR SKY
    //    sky = skyflux[l] / (float)(nb);  // sky per rebinned pixel, e-/frame
    sky = skyflux[l];  // sky per rebinned pixel, e-/frame
  
    if (debug) printf("here4-3\n");
    for ( i=0 ; i<nb*nb ; i++ ) { 
      // bimage[i] += darkcurrent; // nope. has to be added only once/pixel!
      bimage[i] += ( sky + brayleigh[i] ) * bsubmask[i]; 
    }
    if (debug) printf("here4-4\n");
    
    cpu5  = p_cpu_secs(&sys);
    cpu54 += cpu5-cpu4;
  
    // put image where it belongs in large image
    koff = imistart[l] + (imjstart[l])*fimnx;
  
    for ( j=0 ; j<nb ;j++) {
      for ( i=0 ; i<nb ;i++) {
        k = koff + i + j*fimnx;
        *(fimage+k) += bimage[i+j*nb];
        //~ if (bimage[i+j*nb]==0) printf("bimage[%d]==0 ",i+j*nb);
      }
    }
    
    cpu6  = p_cpu_secs(&sys);
    cpu65 += cpu6-cpu5;
    
  }  // END LOOP ON SUBAPERTURES
  
  if (debug) {
    printf("\n1-0 %.3f  2-1(%d,%d) %.3f  3-2 %.3f  4-3 %.3f  5-4 %.3f  6-5 %.3f\n",\
     cpu10*1000.,dynrange,n_in_profile,cpu21*1000.,cpu32*1000.,cpu43*1000.,cpu54*1000.,cpu65*1000.);
  }
  //============================
  // END OF LOOP ON SUBAPERTURES
  //============================


  fftwf_destroy_plan(fftps);
  fftwf_destroy_plan(fftpx);
  fftwf_destroy_plan(fftpxi);

  fftwf_free ( A );
  fftwf_free ( result );
  fftwf_free ( Ax );
  fftwf_free ( Kx );
  fftwf_free ( resultx );
  free ( phase_scaled );
  free ( simage );
  if (debug>1) printf("here5\n");
  free ( ximage );
  free ( bimage );
  free ( brayleigh ); /* fixed memleak 2008nov18 */
  free ( bsubmask );

  if (debug>1) printf("here6\n");


/*
all time used by loop, good:
1-0 0.000  2-1 0.000  3-2 73.329  4-3 0.000
*/  


  return (0);
}



int _shwfs_spots2slopes(
    float    *fimage,       // final image with spots
    int      *imistart2,    // vector of i starts of each image
    int      *imjstart2,    // vector of j starts of each image
    int      nsubs,         // # of subapertures
    int      binxy2,        // side size of the image extracted from the
                            // final image which only corresponds to the
                            // subaperture (= wfs._npixels)
    int      fimnx,         // final image X dimension
    int      fimny,         // final image Y dimension
    int      yoff,          // y offset (to process only part of the image, when using svipc)
    float    *centroidw,    // centroid weight vector for centroid computations, X & Y
    long     shthmethod,    // threshold method (1: yao default, 2: podium, 3: brightest pixels)
    float    *threshold,    // vector of threshold (input), dim nsubs
    float    *bias,         // bias array, dim = nsubs*binxy*binxy, in e-
    float    *flat,         // flat array, dim = nsubs*binxy*binxy, normalized @ 1.
    float    ron,           // read-out noise
    long     noise,         // compute noise
    float    *bckgrdcalib,  // background calib array, binxy*binxy*nsubs
    int      bckgrdinit,    // init background processing. fill bckgrdcalib
    int      bckgrdsub,     // subtract background
    int      *validsubs,    // valid subaps within the set of subap 
                            // for which an image is computed (0/1)
    int      *svipc_subok,  // to skip (0) subap for svipc partial spot comput.
    int      niter,         // total # of cycles over which to integrate
    float    *mesvec)       // final measurement vector
           
{
  /* Declarations */
  const         int sizeThArray = binxy2*binxy2;

  float         *bimage2;
  float         *fnoise;
  float         centx, centy, centi, tmpf, thback;
  float         val[sizeThArray];
  float         minval;
  long          nb2;
  int           i, j, k, l, koff, xyoff;
  int           nvalidsubs, vsc; // vsc = valid sub counter
  int           ind[sizeThArray];
  int           ii, kk, minind;
  int           met3threshold;
  
  // Set the size of 2d FFT.
  nb2 = binxy2 * binxy2;
  vsc = -1;
  xyoff = yoff * fimnx;
  //  printf("shthemthod = %ld\n",shthmethod);
  
  // Allocate memory for the input operands and check its availability.
  bimage2      = ( float* ) malloc ( nb2 * sizeof ( float ) );
  fnoise       = ( float* ) malloc ( fimnx * fimny * sizeof ( float ) );
  
  if ( bimage2 == NULL || fnoise == NULL ) { return (1); }


  // how many valid subs?
  nvalidsubs = 0;
  for( i=0; i<nsubs; i++) nvalidsubs += validsubs[i];

  // zero out mesurement vector ( now done in sh_wfs() )

  // ok, so at this point, we have, for niter frames, done the following:
  // added darkcurrent to entire fimage
  // cumulated star images in each subapertures
  // cumulated sky
  // cumulated rayleight, if needed.
  // applied field stop to star, sky and rayleight.
  // we haven't done anything yet about bias and flat
  // 
  // first, if this is a call to establish the calibration background, 
  // let's do it (before we apply bias, flat and noise):
  if (bckgrdinit) {
    for ( i=xyoff; i<(xyoff+fimnx*fimny); i++ ) bckgrdcalib[i] = fimage[i];
  }
    
  
  // first we have to multiply by the gain (flat) before applying
  // the noise:
  for ( i=xyoff; i<(xyoff+fimnx*fimny); i++) fimage[i] *= flat[i];
                
  // Then apply photon and gaussian (read-out) noise
  if ( noise==1 ) {
    // poisson noise
    _poidev(fimage+xyoff,fimnx*fimny);
    // add gaussian noise
    _gaussdev(fnoise,fimnx*fimny);
    for ( i=xyoff; i<(xyoff+fimnx*fimny); i++ ) fimage[i] += ron*fnoise[i-xyoff];
  }
  
  // Then add the bias (counted as many time as there were frames)
  // No noise on the bias, as this is a bias *error* (see below)
  for ( i=xyoff; i<(xyoff+fimnx*fimny); i++) fimage[i] += (float) niter * bias[i];

  // let's clarify something, as I'm not sure it was as I understood it
  // initially: the bias and flat entries in the parfile are for
  // *errors* in the bias and flat. Those are supposed to be calibrated
  // out, hence here we don't correct for the given value, as it is
  // supposed to result from miscalibrations
    
  // if requested, we should now subtract the background:
  if (bckgrdsub) {
    for ( i=xyoff; i<(xyoff+fimnx*fimny); i++ ) fimage[i] -= bckgrdcalib[i];
  }

  //thresholding in case of default yao thresholding or "podium" thresholding
  if (shthmethod == 1 || shthmethod == 2) {
    // 2010apr16: see note on outside pixels not being thresholded
    // in yao_wfs.i (func sh_wfs). not a big issue, but we solve it here.
    // we need to threshold the fimage pixels outside the valid
    // subap. Let's remove this value here and then adjust by the
    // difference below:
    thback = (float) niter * threshold[nsubs];
    for ( i=xyoff; i<(xyoff+fimnx*fimny); i++ ) fimage[i] -= thback;

  
    // OK, now let's finally compute the slopes over all subaps.
    for ( l=0 ; l<nsubs ; l++ ) {

      if (validsubs[l]) vsc++;
    
      if (svipc_subok[l]==0) continue;

      // retrieve the final/integrated subaperture spot image from fimage:
      koff = imistart2[l] + imjstart2[l]*fimnx;  
      for ( j=0 ; j<binxy2 ;j++) {
        for ( i=0 ; i<binxy2 ;i++) {
          k = koff + i + j*fimnx;
          bimage2[i+j*binxy2] = fimage[k];
        }
      }

      // Apply threshold
       tmpf = (float) niter * threshold[l];
       if (shthmethod ==1) {
         for ( i=0; i<nb2; i++ ) { 
           // see above for thback meaning
           bimage2[i] = bimage2[i] + thback - tmpf; 
           if (bimage2[i] < 0.0f) bimage2[i] = 0.0f;
         }
       } else { 
         for ( i=0; i<nb2; i++ ) { 
           // see above for thback meaning
           bimage2[i] = bimage2[i] + thback; 
           if (bimage2[i] < tmpf) bimage2[i] = 0.0f;
         }
       }
    
      // put back image where it belongs in large image for displays
      // a big stupid: the only operation we did was the threshold...
      koff = imistart2[l] + imjstart2[l]*fimnx;
    
      for ( j=0 ; j<binxy2 ;j++) {
        for ( i=0 ; i<binxy2 ;i++) {
          k = koff + i + j*fimnx;
          fimage[k] = bimage2[i+j*binxy2];
        }
      }
    
      // is that a valid subaperture?
      if (validsubs[l]) {
        // Compute centroids
        centx = centy = centi = 0.0f;
        for ( i=0 ; i<binxy2 ; i++ ) {
          for ( j=0 ; j<binxy2 ; j++ ) {
            centx += centroidw[i]*bimage2[i+j*binxy2];
            centy += centroidw[j]*bimage2[i+j*binxy2];
            centi += bimage2[i+j*binxy2];
          }
        }
      
        // Compute measurements
        if (centi > 0.0f) {
          mesvec[vsc]   = centx/centi;
          mesvec[nvalidsubs+vsc] = centy/centi;
        } // otherwise stay at zero (init value)
        //      vsc++; // increment valid subaperture counter
      }

    }
    // thresholding in case of "bright pixels" thresholding
  } else if (shthmethod == 3) {
    // OK, now let's finally compute the slopes over all subaps.
    for ( l=0 ; l<nsubs ; l++ ) {
      
      met3threshold = (int) threshold[l];
      // make sure we don't select a number of brightest pixel
      // <1 or >total number of pixels in subap:
      if (met3threshold < 1)   met3threshold = 1;
      if (met3threshold > nb2) met3threshold = nb2;
      
      if (validsubs[l]) vsc++;
    
      if (svipc_subok[l]==0) continue;
      // retrieve the final/integrated subaperture spot image from fimage
      // search for the brightest pixels
      // start putting back image where it belongs in large image
      // for displays by "zeroing" fimage
      koff = imistart2[l] + imjstart2[l]*fimnx;    
      minval = fimage[koff];
      minind = 0;
      val[0] = minval;
      ind[0] = minind;
      for ( j=0 ; j<binxy2 ;j++) {
        for ( i=0 ; i<binxy2 ;i++) {
          k = koff + i + j*fimnx;
          ii = i+j*binxy2;
          bimage2[ii] = fimage[k];

          // put the first "threshold" values of bimage2 in val and the
          // subscript values in ind look for the minimum value of val
          // and its subscript minind
          if (ii < met3threshold && ii > 0) {
            val[ii] = bimage2[ii];
            ind[ii] = ii;
            if (val[ii] < minval) {
              minval = val[ii];
              minind = ind[ii];
            }

            // compare each following value of bimage2 to the min value of
            // val and if it is larger, replace it
          } else if (ii > 0 && bimage2[ii] > minval) {
            val[minind] = bimage2[ii];
            ind[minind] = ii;

            // look for the minimum value of val and its subscript minind
            minval = val[0];
            minind = 0;
            for ( kk=1 ; kk < met3threshold ; kk++) {
              if (val[kk] < minval) {
                minval = val[kk];
                minind = kk;
              }
            }   
          }
          // zero fimage before filling it by the brightest pixels
          fimage[k] = 0;
        }
      }
      centx = centy = centi = 0.0f;
       
      for (ii=0 ; ii < met3threshold ; ii++) {
        i = ind[ii]%binxy2;
        j = ind[ii]/binxy2;
        k = koff + i + j*fimnx;
        // fill fimage with the brightest pixels
        fimage[k] = val[ii];
        // is that a valid subaperture?
        if (validsubs[l]) {
          // Compute centroids
          centx += centroidw[i]*val[ii];
          centy += centroidw[j]*val[ii];
          centi += val[ii];
        }
      }
      // Compute measurements
      if (centi > 0.0f) {
        mesvec[vsc]   = centx/centi;
        mesvec[nvalidsubs+vsc] = centy/centi;
      } // otherwise stay at zero (init value)
 
    }
  }
  // with respect to the noise in the outskirt of the fimage:
  // the threshold has been applied everywhere, but the condition
  // if (pixel<0) pixel=0 only has been applied in the valid
  // subapertures. Let's make sure it is applied everywhere here:
  for ( i = xyoff; i < (xyoff+fimnx*fimny); i++ ) {
    if (fimage[i] < 0.0f) fimage[i] = 0.0f;
  }

  free ( bimage2 );
  free ( fnoise );

  return (0);
}


int _shwfs_simple(float *pupil,      // input pupil
                  float *phase,      // input phase
                  float phasescale,  // phase scaling factor
                  float *phaseoffset,// input phase offset
                  int dimx,          // X dim of phase. Use as stride for extraction
                  int dimy,          // Y dim of phase. Use as stride for extraction
     
                  int *istart,       // vector of i starts of each subaperture
                  int *jstart,       // vector of j starts of each subaperture
                  int nx,            // subaperture i size
                  int ny,            // subaperture j size
                  int nsubs,         // # of subapertures
     
                  float toarcsec,    // conversion factor to arcsec
     
                  float *mesvec)     // final measurement vector

{
  /* Declarations */
  float         avgx, avgy, avgi;
  int           i,j,k,l,koff;

  // loop on subapertures:
  for ( l=0 ; l<nsubs ; l++ ) {

    koff = istart[l] + jstart[l]*dimx;
    
    avgx = 0.0f;
    avgy = 0.0f;
    avgi = 0.0f;
    
    for ( j=0; j<ny ; j++ ) {
      for ( i=0; i<nx ; i++ ) {
  
        k = koff + i + j*dimx;
        // the term between parenthesis is the 2 point estimate 
        // of the phase X derivative
        // take care of outliers:
        if ( (istart[l]==0) & (i==0) ) { // start of a row
          avgx += pupil[k] * phasescale * \
            (phase[k+1]-phase[k]+phaseoffset[k+1]-phaseoffset[k]);
        } else if ( ((istart[l]+nx)>=dimx) & (i==(nx-1)) ) { // end of a row
          avgx += pupil[k] * phasescale * \
            (phase[k]-phase[k-1]+phaseoffset[k]-phaseoffset[k-1]);
        } else if (pupil[k-1]==0) { // edge of the pupil 
          avgx += pupil[k] * phasescale * \
            (phase[k+1]-phase[k]+phaseoffset[k+1]-phaseoffset[k]);
        } else if (pupil[k+1]==0) { // edge of the pupil 
          avgx += pupil[k] * phasescale * \
            (phase[k]-phase[k-1]+phaseoffset[k]-phaseoffset[k-1]);
        } else { // then 2 neightbors derivative estimate
          avgx += pupil[k] * phasescale * \
            (phase[k+1]-phase[k-1]+phaseoffset[k+1]-phaseoffset[k-1])/2.;
        }
        // same for y:
        if ( (jstart[l]==0) & (j==0) ) { // start of a column
          avgy += pupil[k] * phasescale * \
            (phase[k+dimx]- phase[k]+phaseoffset[k+dimx]-phaseoffset[k]);
        } else if ( ((jstart[l]+ny)>=dimy) & (j==(ny-1)) ) { // end of a column
          avgy += pupil[k] * phasescale * \
            (phase[k]- phase[k-dimx]+phaseoffset[k]-phaseoffset[k-dimx]);
        } else if (pupil[k-dimx]==0) { // edge of the pupil
          avgy += pupil[k] * phasescale * \
            (phase[k+dimx]- phase[k]+phaseoffset[k+dimx]-phaseoffset[k]);
        } else if (pupil[k+dimx]==0) { // edge of the pupil
          avgy += pupil[k] * phasescale * \
            (phase[k]- phase[k-dimx]+phaseoffset[k]-phaseoffset[k-dimx]);
        } else {
          avgy += pupil[k] * phasescale * \
            (phase[k+dimx]- phase[k-dimx]+phaseoffset[k+dimx]-phaseoffset[k-dimx])/2.;
        }
        avgi += pupil[k];
      }
    }
    
    if (avgi > 0.0f) {
      mesvec[l]   = avgx/avgi*toarcsec;
      mesvec[nsubs+l] = avgy/avgi*toarcsec;
    } else {
      mesvec[l]   = 0.0f;
      mesvec[nsubs+l] = 0.0f;
    }
  }

  return (0);
}


int _cwfs (float *pupil,      // input pupil
           float *phase,      // input phase
           float phasescale,  // phase scaling factor
           float *phaseoffset,// input phase offset
           float *cxdef,      // cos(defoc)
           float *sxdef,      // sin(defoc)
           int dimpow2,       // dim of phase in powers of 2
           int *sind,         // array containing the indices of subaperture#i
           int *nsind,        // nsind(i) = number of valid pixels in sind(i)
           int nsubs,         // number of subapertures
           float *fimage1,    // final image1 with spots
           float *fimage2,    // final image2 with spots
           float nphotons,    // total number of photons, for noise
           float skynphotons, // total number of photons from sky, for noise
           float ron,         // read-out noise, in case.
           float darkcurrent, // dark current, per detector/full sample
                              // note it will be 1/2 of given value per image
                              // ron and darkcurrent are only added if noise = 1
           int noise,         // enable noise ?
           float *mesvec)     // final measurement vector

{
  fftwf_complex *A, *B, *result;
  float         *ptr,*ptr1;
  fftwf_plan     p,p1;
  float         *x1,*x2,*gnoise;
  float         tot;
  long          log2nr, log2nc, n, ns;
  int           i,k,sindstride,koff;
  float         pp,ppsin,ppcos;

  /*
    Global setup for FFTs:
  */

  // fftwf_plan_with_nthreads(n_threads);

  sindstride = 0;
  for ( i=0 ; i<nsubs ; i++ ) {
    if (nsind[i] > sindstride) { sindstride = nsind[i];}
  }

  // Set the size of 2d FFT.
  log2nr = dimpow2; 
  log2nc = dimpow2;
  n  = 1 << ( log2nr + log2nc ); // total number of pixels in small array
  ns = 1 << log2nr; // total number of pixels in small array

  // Allocate memory for the input operands and check its availability.
  x1     = ( float* ) malloc ( nsubs * sizeof ( float ) );
  x2     = ( float* ) malloc ( nsubs * sizeof ( float ) );
  gnoise = ( float* ) malloc ( nsubs * sizeof ( float ) );
  A      = fftwf_malloc ( n * sizeof ( fftwf_complex ) );
  B      = fftwf_malloc ( n * sizeof ( fftwf_complex ) );
  result = fftwf_malloc ( n * sizeof ( fftwf_complex ) );

  if ( A == NULL || B == NULL || result == NULL ) { return (1); }

  for (i=0;i<nsubs;i++) {
    x1[i] = 0.0f;
    x2[i] = 0.0f;
  }
      
  // Set up the required memory for the FFT routines and 
  // check its availability.
  p  = fftwf_plan_dft_2d(ns, ns, A, B, FFTW_FORWARD, FFTW_ESTIMATE);
  p1 = fftwf_plan_dft_2d(ns, ns, A, result, FFTW_BACKWARD, FFTW_ESTIMATE);


  // intermediate FFT, to find intermediate complex amplitude
  // fill A
  ptr = (void *)A;
  for ( i=0; i<n ; i++ ) {
    if (pupil[i] != 0.0f) {
      pp = phasescale*(phase[i]+phaseoffset[i]);
      if (use_sincos_approx_flag) _sinecosinef(pp,&ppsin,&ppcos);
      else sincosf(pp,&ppsin,&ppcos);
      *(ptr)   = pupil[i]*ppcos;
      *(ptr+1) = pupil[i]*ppsin;
    } else {
      *(ptr)   = 0.0f;
      *(ptr+1) = 0.0f;
    }
    ptr += 2;
  }


  // Carry out a Forward 2d FFT transform, check for errors.
  fftwf_execute(p); /* repeat as needed */

  // image #1:
  ptr = (void *)A;
  ptr1 = (void *)B;
  for ( i=0; i<n ; i++ ) {
    *(ptr)   = *(ptr1)*cxdef[i] - *(ptr1+1)*sxdef[i];
    *(ptr+1) = *(ptr1)*sxdef[i] + *(ptr1+1)*cxdef[i];
    ptr += 2; ptr1 += 2;
  }

  fftwf_execute(p1); /* repeat as needed */

  ptr = (void *)result;
  for ( i=0; i<n; i++ ) {
    fimage1[i] = (*(ptr) * *(ptr) + *(ptr+1) * *(ptr+1) );
    ptr += 2;
  }

  // image #2:
  ptr = (void *)A;
  ptr1 = (void *)B;
  for ( i=0; i<n ; i++ ) {
    *(ptr)   = *(ptr1)*cxdef[i]   + *(ptr1+1)*sxdef[i];
    *(ptr+1) = *(ptr1+1)*cxdef[i] - *(ptr1)*sxdef[i];
    ptr += 2; ptr1 += 2;
  }

  fftwf_execute(p1); /* repeat as needed */

  ptr = (void *)result;
  for ( i=0; i<n; i++ ) {
    fimage2[i] = (*(ptr) * *(ptr) + *(ptr+1) * *(ptr+1) );
    ptr += 2;
  }

  // now we got to sum the relevant pixels:
  for ( k=0 ; k<nsubs ; k++ ) {
    koff = k*sindstride;
    //printf("nsind[%d]=%d\n",k,nsind[k]);
    for ( i=0 ; i<nsind[k] ; i++ ) {
      x1[k] += fimage1[sind[koff+i]];
      x2[k] += fimage2[sind[koff+i]];
    }
  }

  // NOISE:
  if (noise == 1) {
    // x1:
    tot = 0.0f;
    // compute total of current flux vector
    for ( i=0 ; i<nsubs ; i++ ) { tot += x1[i]; }
    if (tot > 0.0f) { 
      tot = nphotons/2.0f/tot;
      // normalize so that sum(x1) = nphotons
      for ( i=0 ; i<nsubs ; i++ ) { x1[i] = x1[i]*tot ; }
      // add darkcurrent:
      for ( i=0 ; i<nsubs ; i++ ) { x1[i] += darkcurrent/2.0f ; }
      // add sky:
      for ( i=0 ; i<nsubs ; i++ ) { x1[i] += skynphotons/2.0f ; }
      // apply poisson noise
      _poidev(x1,nsubs);
    }
    if (ron > 0.0f) {
      // set up gaussian noise vector
      for ( i=0 ; i<nsubs ; i++ ) { gnoise[i] = ron; }
      // draw gaussian noise
      _gaussdev(gnoise,nsubs);
      // add to x1
      for ( i=0 ; i<nsubs ; i++ ) { x1[i] += gnoise[i]; }
    }

    // x2, same comment as above for x1:
    tot = 0.0f;
    for ( i=0 ; i<nsubs ; i++ ) { tot += x2[i]; }
    if (tot > 0.0f) { 
      tot = nphotons/2.0f/tot;
      for ( i=0 ; i<nsubs ; i++ ) { x2[i] = x2[i]*tot ; }
      for ( i=0 ; i<nsubs ; i++ ) { x2[i] += darkcurrent/2.0f ; }
      for ( i=0 ; i<nsubs ; i++ ) { x2[i] += skynphotons/2.0f ; }
      _poidev(x2,nsubs);
    }
    if (ron > 0.0f) {
      for ( i=0 ; i<nsubs ; i++ ) { gnoise[i] = ron; }
      _gaussdev(gnoise,nsubs);
      for ( i=0 ; i<nsubs ; i++ ) { x2[i] += gnoise[i]; }
    }    
  }

  for ( i=0 ; i<nsubs ; i++ ) {
    if ( (x1[i]+x2[i]) == 0.0f ) {
      mesvec[i] = 0.0f;
    } else {
      mesvec[i] = (x1[i]-x2[i])/(x1[i]+x2[i]);
    }
  }

  fftwf_destroy_plan(p);
  fftwf_destroy_plan(p1);
  free ( x1 );
  free ( x2 );
  free ( gnoise );
  fftwf_free ( A );
  fftwf_free ( B );
  fftwf_free ( result );

  return (0);
}
