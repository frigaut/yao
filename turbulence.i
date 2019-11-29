/*
 * turbulence.i
 *
 * A collection of routines about turbulence.
 *
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2017, Francois Rigaut
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

// define function names. Will be redefined by GUI routine if needed:
func null (arg,..) { return 0; }
pyk_error=pyk_info=pyk_warning=gui_message=gui_progressbar_frac=gui_progressbar_text=null;

//+++++++++++++++++++++++++++

func phase_struct_func(phase,npt,step,&xsep,plot=,color=)
  /* DOCUMENT phase_struct_func(phase,npt,step,plot=)
     Compute (and plot) structure function along the first dimension
     of array. Transpose if you need to plot structure function along
     another dimension. The structure function is defined as:
     SF(j) = avg( ( phase - shift_along_dim1_by_j (phase) )^2. )
     F.Rigaut, 2001/11/10.
     SEE ALSO: -
  */

{
  if (is_void(step)) {step=1;}
  xm = min([npt,((dimsof(phase))(2))/step-1]);
  psfx = array(float,xm);
  for (i=step;i<=xm*step;i=i+step) {
    psfx(i/step) = avg((phase(1:-i,..)-phase(1+i:,..))^2.);
  }
  psfx = grow(0.,psfx);
  xsep = step*(indgen(xm+1)-1.);
  if (is_set(plot)) { plg,psfx,xsep,color=color; }
  return psfx;
}
phaseStructFunc = phase_struct_func;

//+++++++++++++++++++++++++++

func create_phase_screens(dimx,dimy,l0=,prefix=,nalias=,no_ipart=,nogsl=,silent=)
  /* DOCUMENT create_phase_screens(dimx,dimy,prefix=)
     Create phase screens and save them in fits files.
     The saved phase screens have a dimension dimx*dimy.
     Number of phase screens = dimx/dimy.
     The phase screens are normalized so that 1 pixel = 1 r0,
     i.e. the variance of the squared difference of the screen
     with itself at one pixel interval is 6.88 (rd^2).

     dimx = long dimension of result screens
     dimy = short dimension of result screens
     prefix = Prefix to filename. if prefix is not set, the
       screens are returned by not saved.
     no_ipart = set to loose the imaginary part. Gain some RAM
       to reach larger dimension.

     Example:
     create_phase_screens,2048,256,prefix="screen256"

     F.Rigaut, 2001/11/10.
     modify 2003 Feb 24 to add dimy (before dimy=256) and prefix
     SEE ALSO: generate_phase, phase_struct_func.
  */

{
  if (is_void(l0)) l0 = 0.;
  if (yaopy) gui_progressbar_text,"Generating the screen power spectrum";
  nps = (no_ipart?1:2);
  nscreen = dimx/dimy*nps;
  pscreen = generate_phase_with_L0(dimx,l0,nalias=nalias,no_ipart=no_ipart,silent=silent);
  if (yaopy) gui_progressbar_frac,0.25;

  if (!silent) write,"Normalizing phase screens";
  if (yaopy) gui_progressbar_text,"Normalizing phase screens";
  off = [1,5]; // spatial offset for structure function normalization
  psfunc = array(float,off(2));
  for (i=off(1);i<=off(2);i++ ) {
    // need to split over phase screens to reduce load on RAM
    fsx = fsy = 0.;
    for (n=1;n<=nps;n++) {
      fsx += avg((pscreen(1:-i,,n)-pscreen(i+1:,,n))^2.);
      fsy += avg((pscreen(,1:-i,n)-pscreen(,i+1:,n))^2.);
    }
    fsx /= nps; fsy /= nps;
    psfunc(i) = sqrt((fsx+fsy)/2.);
    if (yaopy) gui_progressbar_frac,0.25+0.6*(i-off(1))/(off(2)-off(1));
  }

  c = (24./5.*gamma(6/5.))^(5/6.);
  r = float(indgen(off(2)));
  if (l0 == 0){
    theo = sqrt(2*c*r^(5./3.));
  } else {
    f0 = 1./l0;
    c = (24./5.*gamma(6/5.))^(5/6.);
    r = float(indgen(off(2)));

    include, "gsl.i", 3; // loads the Bessel functions, but does not crash if function does not exist
    if ((gsl_sf_bessel_Knu == [])||(nogsl)||(dont_use_gsl)) {
      // function does not exist, use an approximation
      // flags for comparison tests
      write, "*** WARNING: gsl_sf_bessel_Knu is not defined ***";
      write, "Normalization might be wrong with finite outer scale";
      write, "Please install ygsl (https://github.com/emmt/ygsl)";
      write, "*************************************************";
      theo = psfunc; // do not renormalize
    } else { // bessel function is defined
      theo = sqrt(2*c*gamma(11./6.)/(2^(5./6)*pi^(8./3))*(f0)^(-5./3)*(gamma(5./6.)/2^(1./6.) - (2*pi*r*f0)^(5./6.)*gsl_sf_bessel_Knu(5./6., 2*pi*r*f0)));
    }
  }
  nfact = avg(psfunc(off(1):off(2))/theo(off(1):off(2)));
  if (!silent) write,format="normalization factor (actual/theo)= %f\n",
    avg(psfunc(off(1):off(2))/theo(off(1):off(2)));
  if (!silent) write,psfunc(off(1):off(2))/theo(off(1):off(2));

  pscreen(*) = pscreen(*)/float(nfact);

  if (!silent) write,"Splitting and saving phase screens";
  if (yaopy) gui_progressbar_text,"Splitting and saving phase screens";
  pscreen = reform(pscreen,[3,dimx,dimy,nscreen]);
  if (!is_void(prefix)) {
    for (i=1;i<=nscreen;i++) {
      fname = prefix+((nscreen==1) ? "":swrite(i,format="%i"))+".fits";
      fits_write,fname,pscreen(,,i),overwrite=1;
      if (yaopy) gui_progressbar_text,swrite(format="Saving %s",fname);
      if (yaopy) gui_progressbar_frac,0.85+0.15*i/nscreen;
    }
  }

  if (yaopy) after,4,clean_progressbar;
  return pscreen;
}
createPhaseScreens = create_phase_screens;

func clean_progressbar(void)
{
  gui_progressbar_text,"";
  gui_progressbar_frac,0.;
}

//+++++++++++++++++++++++++++

func generate_phase(dim,silent=)
  /* DOCUMENT generate_phase(size)
     Generate by Fourier an un-normalized 2D phase screen from the
     -11/3 amplitude and a randomn phase component. Only returns the
     real part. Beware that these screens have a effective outer scale
     of about half the length of the screen.
     F.Rigaut, 2001/11/10.
     SEE ALSO: create_phase_screens, phase_struct_func.
  */

{
  if (!do_not_randomize_pscreens) randomize;
  if (!silent) write,"Creating arrays";
  if (!silent) write,"Creating amplitude";
  tmp   = clip(dist(dim),1e-8,);
  amp = eclat(tmp^(-11.f/6.f));
  p = array(complex,dim,dim);
  p.re  = amp;
  p.im  = amp;
  amp = [];
  if (!silent) write,"Creating phase";
  pha = float(random(dim,dim)*2.*pi);
  p.re  = p.re*cos(pha);
  p.im  = p.im*sin(pha);
  pha = [];
  p.re(1,1)= 0.;
  p.im(1,1)= 0.;
  if (!silent) write,"Doing FFT...";
  phaout  = float(fft(p,1));
  p = [];
  return phaout;
}

//+++++++++++++++++++++++++++

func generate_von_karman_spectrum(dim,k0,nalias=,silent=)
/* DOCUMENT func generate_von_karman_spectrum(sdim,bdim,k0)
   generate correct von Karman spectrum including aliasing.

   SEE ALSO:
 */
{
  if (is_void(nalias)) nalias = 0;
  if (!silent) for (i=1;i<=2*nalias+1;i++) write,format="%s","#";

  res = array(float,[2,dim,dim]);
  for (i=-nalias;i<=nalias;i++) {
    for (j=-nalias;j<=nalias;j++) {
      if ((i==0) && (j==0)) {
        // bug ((1,1) pixel shift fixed 2009oct13)
        tmp = sqrt(dist(dim)^2.f+k0^2.);
      } else {
        tmp = sqrt(dist(dim,xc=i*dim+dim/2,yc=j*dim+dim/2)^2.f+k0^2.);
      }
      if ((i==0) && (j==0)) tmp = clip(tmp,1.,);
      amp = (6.88*0.00969)*tmp^(-11.f/6.f);
      res += amp;
    }
    if (!silent) write,format="%s","#";
  }
  if (!silent) write,format="%s\n"," > Done";
  roll,res;
  res = float(res);
  return res;
}
//+++++++++++++++++++++++++++


func generate_phase_with_L0(dim,l0,nalias=,silent=,no_ipart=)
  /* DOCUMENT generate_phase(size,l0)
     Generate by Fourier an un-normalized 2D phase screen from the
     -11/3 amplitude and a randomn phase component. Returns the real and
     complex parts.
     Uses fftVE and cosf/sinf to keep floats for RAM use consideration
     (the previous version of this routine was using the yorick fft,
     thus double complex, which limits things on my machine to 4096 screens).

     dim: desired dimension of the result phase screens
     l0: outer scale IN PIXELS
     no_ipart = loose the im part to gain RAM

     F.Rigaut, 2001/11/10.
     SEE ALSO: create_phase_screens, phase_struct_func.
  */

{
  if (l0 == 0.) { k0=0.0f; } else { k0 = float(dim)/l0; }
  if (!do_not_randomize_pscreens) randomize;
  gui_progressbar_frac,0.01;
  if (!silent) {
    write,"Creating arrays";
    write,"Creating amplitude";
  }
  //  tmp = clip(float(sqrt(dist(dim)^2.f+k0^2.)),1e-8,);
  //  amp = 6.88*0.00969*dim*eclat(tmp^(-11.f/6.f));
  //  amp = eclat(generate_von_karman_spectrum(dim,clip(2*dim,,2048),k0));
  //  amp = dim*eclat(generate_von_karman_spectrum(dim,1024,k0));
  amp = dim*generate_von_karman_spectrum(dim,k0,nalias=nalias,silent=silent);
  gui_progressbar_frac,0.10;
  amp = float(amp);
  // normalized so that the structure function is correct for r0 = 1 pixel
  tmp = [];
  if (!silent) write,"Creating phase";
  pha = reform(float(2*pi*random(dim*dim)),dimsof(amp));
  re  = amp*cosf(pha);
  im  = amp*sinf(pha);
  amp = pha = [];
  re(1,1) = 0.0f;
  im(1,1) = 0.0f;
  if (!silent) write,"Doing FFT...";
  phaout = fftVE(re,im,1);
  if (no_ipart) phaout = phaout(,,1);
  gui_progressbar_frac,0.20;
  re = im = [];
  return phaout;
}
