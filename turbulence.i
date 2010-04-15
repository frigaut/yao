/*
 * TURBULENCE.I
 * A collection of routines about turbulence.
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: turbulence.i,v 1.3 2010-04-15 02:36:53 frigaut Exp $
 *
 * Copyright (c) 2002-2007, Francois Rigaut
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
 * $Log: turbulence.i,v $
 * Revision 1.3  2010-04-15 02:36:53  frigaut
 *
 *
 * final commit to upgrade this repo to yao 4.5.1
 *
 * Revision 1.2  2007/12/19 13:18:59  frigaut
 * - explicit message when screens are not present/found
 * - more messages in statusbar
 * - added statusbar1 (that can hide/show) for strehl status header
 *
 * Revision 1.1.1.1  2007/12/12 23:29:13  frigaut
 * Initial Import - yorick-yao
 *
 *
 */

// define function names. Will be redefined by GUI routine if needed:
func null (arg,..) { return 0; }
pyk_error=pyk_info=pyk_warning=gui_message=gui_progressbar_frac=gui_progressbar_text=null; 

//+++++++++++++++++++++++++++

func phase_struct_func(phase,npt,step,plot=)
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
  if (is_set(plot)) { fma; plg,psfx,step*(indgen(xm+1)-1.); }
  return psfx;
}
phaseStructFunc = phase_struct_func;

//+++++++++++++++++++++++++++

func create_phase_screens(dimx,dimy,l0=,prefix=,nalias=)
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

     Example:
     create_phase_screens,2048,256,prefix="screen256"
     
     F.Rigaut, 2001/11/10.
     modify 2003 Feb 24 to add dimy (before dimy=256) and prefix
     SEE ALSO: generate_phase, phase_struct_func.
  */ 

{
  if (is_void(l0)) l0 = 0.;
  gui_progressbar_text,"Generating the screen power spectrum";
  nscreen = dimx/dimy*2;
  pscreen = generate_phase_with_L0(dimx,l0,nalias=nalias);
  gui_progressbar_frac,0.25;

  print,"Normalizing phase screens";
  gui_progressbar_text,"Normalizing phase screens";
  off = [1,5]; // spatial offset for structure function normalization
  psfunc = array(float,off(2));
  for (i=off(1);i<=off(2);i++ ) {
    fsx     = avg((pscreen(1:-i,,)-pscreen(i+1:,,))^2.);
    fsy     = avg((pscreen(,1:-i,)-pscreen(,i+1:,))^2.);
    psfunc(i)      = sqrt((fsx+fsy)/2.);
    gui_progressbar_frac,0.25+0.6*(i-off(1))/(off(2)-off(1));
  }
  theo = sqrt(6.88*indgen(off(2))^1.66);

  //  write,format="normalization factor = %f\n",sqrt(6.88*off^1.666/fs);
  write,format="normalization factor (actual/theo)= %f\n",
    (nfact=avg(psfunc(off(1):off(2))/theo(off(1):off(2))));
  write,psfunc(off(1):off(2))/theo(off(1):off(2));

  pscreen = pscreen/nfact;
  

  print,"Sectionning and saving phase screens";
  gui_progressbar_text,"Sectionning and saving phase screens";
  pscreen = reform(pscreen,[3,dimx,dimy,nscreen]);
  if (!is_void(prefix)) {
    for (i=1;i<=nscreen;i++) {
      fname = prefix+((nscreen==1) ? "":swrite(i,format="%i"))+".fits";
      fitsWrite,fname,pscreen(,,i);
      gui_progressbar_text,swrite(format="Saving %s",fname);
      gui_progressbar_frac,0.85+0.15*i/nscreen;
    }
  }

  after,4,clean_progressbar;
  return pscreen;
}
createPhaseScreens = create_phase_screens;

func clean_progressbar(void)
{
  gui_progressbar_text,"";
  gui_progressbar_frac,0.;
}

//+++++++++++++++++++++++++++

func generate_phase(dim)
  /* DOCUMENT generate_phase(size)
     Generate by Fourier an un-normalized 2D phase screen from the 
     -11/3 amplitude and a randomn phase component. Only returns the 
     real part. Beware that these screens have a effective outer scale 
     of about half the length of the screen. 
     F.Rigaut, 2001/11/10.
     SEE ALSO: create_phase_screens, phase_struct_func.
  */

{
  randomize;
  print,"Creating arrays";
  print,"Creating amplitude";
  tmp   = clip(dist(dim),1e-8,);
  amp = eclat(tmp^(-11.f/6.f));
  p = array(complex,dim,dim);
  p.re  = amp;
  p.im  = amp;
  amp = [];
  print,"Creating phase";
  pha = float(random(dim,dim)*2.*pi);
  p.re  = p.re*cos(pha);
  p.im  = p.im*sin(pha);
  pha = [];
  p.re(1,1)= 0.;
  p.im(1,1)= 0.;
  print,"Doing FFT...";
  phaout  = float(fft(p,1));
  p = [];
  return phaout;
}

//+++++++++++++++++++++++++++

func generate_von_karman_spectrum(dim,k0,nalias=,silent=)
/* DOCUMENT func generate_von_karman_spectrum(sdim,bdim,k0)
   generate correct VoKarman spectrum including aliasing.
     
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


func generate_phase_with_L0(dim,l0,nalias=,silent=)
  /* DOCUMENT generate_phase(size,l0)
     Generate by Fourier an un-normalized 2D phase screen from the 
     -11/3 amplitude and a randomn phase component. Returns the real and 
     complex parts.
     Uses fftVE and cosf/sinf to keep floats for RAM use consideration
     (the previous version of this routine was using the yorick fft,
     thus double complex, which limits things on my machine to 4096 screens).

     dim: desired dimension of the result phase screens
     l0: outer scale IN PIXELS
     
     F.Rigaut, 2001/11/10.
     SEE ALSO: create_phase_screens, phase_struct_func.
  */

{
  if (l0 == 0.) { k0=0.0f; } else { k0 = float(dim)/l0; }
  randomize;
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
  gui_progressbar_frac,0.20;
  re = im = [];
  return phaout;
}

