/*
 * TURBULENCE.I
 * A collection of routines about turbulence.
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: turbulence.i,v 1.1 2007-12-12 23:29:13 frigaut Exp $
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
 * Revision 1.1  2007-12-12 23:29:13  frigaut
 * Initial revision
 *
 *
 */

//+++++++++++++++++++++++++++

func PhaseStructFunc(phase,npt,step,plot=)
  /* DOCUMENT PhaseStructFunc(phase,npt,step,plot=)
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
phaseStructFunc = PhaseStructFunc;

//+++++++++++++++++++++++++++

func CreatePhaseScreens(dimx,dimy,l0=,prefix=,nalias=)
  /* DOCUMENT CreatePhaseScreens(dimx,dimy,prefix=)
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
     CreatePhaseScreens,2048,256,prefix="screen256"
     
     F.Rigaut, 2001/11/10.
     modify 2003 Feb 24 to add dimy (before dimy=256) and prefix
     SEE ALSO: generate_phase, PhaseStructFunc.
  */ 

{
  if (is_void(l0)) l0 = 0.;
  
  nscreen = dimx/dimy*2;
  pscreen = generatePhaseWithL0(dimx,l0,nalias=nalias);

  print,"Normalizing phase screens";

  off = [1,5];
  psfunc = array(float,off(2));
  for (i=off(1);i<=off(2);i++ ) {
    fsx     = avg((pscreen(1:-i,,)-pscreen(i+1:,,))^2.);
    fsy     = avg((pscreen(,1:-i,)-pscreen(,i+1:,))^2.);
    psfunc(i)      = sqrt((fsx+fsy)/2.);
  }
  theo = sqrt(6.88*indgen(off(2))^1.66);

  //  write,format="normalization factor = %f\n",sqrt(6.88*off^1.666/fs);
  write,format="normalization factor (actual/theo)= %f\n",
    (nfact=avg(psfunc(off(1):off(2))/theo(off(1):off(2))));
  write,psfunc(off(1):off(2))/theo(off(1):off(2));

  pscreen = pscreen/nfact;
  
  print,"Sectionning and saving phase screens";
  pscreen = reform(pscreen,[3,dimx,dimy,nscreen]);
  if (!is_void(prefix)) {
    for (i=1;i<=nscreen;i++) {
      fname = prefix+((nscreen==1) ? "":swrite(i,format="%i"))+".fits";
      fitsWrite,fname,pscreen(,,i);
    }
  }

  return pscreen;
}
createPhaseScreens = CreatePhaseScreens;

//+++++++++++++++++++++++++++

func generate_phase(dim)
  /* DOCUMENT generate_phase(size)
     Generate by Fourier an un-normalized 2D phase screen from the 
     -11/3 amplitude and a randomn phase component. Only returns the 
     real part. Beware that these screens have a effective outer scale 
     of about half the length of the screen. 
     F.Rigaut, 2001/11/10.
     SEE ALSO: CreatePhaseScreens, PhaseStructFunc.
  */

{
  randomize;
  print,"Creating arrays";
  print,"Creating amplitude";
  tmp   = clip(dist(dim),1e-8,);
  amp	= eclat(tmp^(-11.f/6.f));
  p	= array(complex,dim,dim);
  p.re	= amp;
  p.im	= amp;
  amp	= [];
  print,"Creating phase";
  pha	= float(random(dim,dim)*2.*pi);
  p.re	= p.re*cos(pha);
  p.im	= p.im*sin(pha);
  pha	= [];
  p.re(1,1)= 0.;
  p.im(1,1)= 0.;
  print,"Doing FFT...";
  phaout	= float(fft(p,1));
  p	= [];
  return phaout;
}

//+++++++++++++++++++++++++++

func generateVKspectrum(dim,k0,nalias=)
/* DOCUMENT func generateVKspectrum(sdim,bdim,k0)
   generate correct VoKarman spectrum including aliasing.
     
   SEE ALSO:
 */
{
  if (is_void(nalias)) nalias = 0;
  for (i=1;i<=2*nalias+1;i++) {write,format="%s","#";}
  write,format="%s\n"," < - number or row to do";
  
  res = array(float,[2,dim,dim]);
  for (i=-nalias;i<=nalias;i++) {
    for (j=-nalias;j<=nalias;j++) {
      tmp = sqrt(dist(dim,xc=i*dim+dim/2,yc=j*dim+dim/2)^2.f+k0^2.);
      if ((i==0) && (j==0)) tmp = clip(tmp,1.,);
      amp = (6.88*0.00969)*tmp^(-11.f/6.f);
      res += amp;
      //print,i,j;tv,res; pause,500;
    }
    write,format="%s","#";
  }
  write,format="%s\n"," > Done";
  roll,res;
  res = float(res);
  return res;
}
//+++++++++++++++++++++++++++

func generateVKspectrum_old(sdim,bdim,k0)
/* DOCUMENT func generateVKspectrum(sdim,bdim,k0)
   generate correct VoKarman spectrum including aliasing.
     
   SEE ALSO:
 */
{
  tmp = clip(float(sqrt(dist(bdim)^2.f+k0^2.)),1.f,);
  amp = 6.88*0.00969*tmp^(-11.f/6.f);
  amp = roll(amp,[(sdim-bdim)/2,(sdim-bdim)/2]);
  tmp = amp;
  for (i=1;i<=bdim/sdim-1;i++) {
    tmp += roll(amp,[i*sdim,0]);
    //tv,tmp; //pause,500;
  }
  amp = tmp;
  for (j=1;j<=bdim/sdim-1;j++) {
    tmp += roll(amp,[0,j*sdim]);
    //tv,tmp; //pause,500;
  }
  return tmp(1:sdim,1:sdim);
}
  
func generatePhaseWithL0(dim,l0,nalias=)
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
     SEE ALSO: CreatePhaseScreens, PhaseStructFunc.
  */

{
  if (l0 == 0.) { k0=0.0f; } else { k0 = float(dim)/l0; }
  randomize;
  print,"Creating arrays";
  print,"Creating amplitude";
  //  tmp = clip(float(sqrt(dist(dim)^2.f+k0^2.)),1e-8,);
  //  amp = 6.88*0.00969*dim*eclat(tmp^(-11.f/6.f));
  //  amp = eclat(generateVKspectrum(dim,clip(2*dim,,2048),k0));
  //  amp = dim*eclat(generateVKspectrum(dim,1024,k0));
  amp = dim*generateVKspectrum(dim,k0,nalias=nalias);
  amp = float(amp); sum(amp);
  // normalized so that the structure function is correct for r0 = 1 pixel
  tmp = [];
  print,"Creating phase";
  pha = reform(float(2*pi)*gaussdev(dim*dim),dimsof(amp));
  re  = amp*cosf(pha);
  im  = amp*sinf(pha);
  amp = pha = [];
  re(1,1) = 0.0f;
  im(1,1) = 0.0f;
  print,"Doing FFT...";
  phaout = fftVE(re,im,1);
  re = im = [];
  return phaout;
}

