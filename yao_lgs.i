/*
 * yao_lgs.i
 *
 * Compilation of functions related to LGS and LGS profile operations
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
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

func fit_lgs_profile(amp,alt,npt,&fitamp,&fitalt,&fitdepth)
/* DOCUMENT fit_lgs_profile(amp,alt,npt,&fitamp,&fitalt,&fitdepth)
  amp = [0.00125154,0.00143618,0.00161906,0.0481773,0.043927,0.0533849,0.0932741,
      0.0816419,0.0849489,0.155098,0.146013,0.130336,0.096709,0.022861,0.0130669,
      0.00831713,0.00523775,0.0045743,0.0048842,0.00324208];
  alt = [90,91.3158,92.6316,93.9474,95.2632,96.5789,97.8947,99.2105,100.526,101.842,
        103.158,104.474,105.789,107.105,108.421,109.737,111.053,112.368,113.684,115]*1e3;
  fit_lgs_profile,amp,alt,7,fitamp,fitalt,fitdepth;
  sigma = 0.987607
  [0.0541135,-0.0128751,0.0859853,0.0686155,0.168526,0.130719,0.0148146]
  [94.5293,145.704,97.6732,99.8664,102.383,105.023,108.472]
  fit_lgs_profile,amp,alt,5,fitamp,fitalt,fitdepth;
  sigma = 1.080501
  [0.0511654,0.0832791,0.0511242,0.159607,0.12324]
  [94.5549,97.813,99.847,102.332,105.047]
  for (i=1;i<=6;i++) {
    wfs(i).lgs_prof_amp= &float(fitamp);
    wfs(i).lgs_prof_alt= &float(fitalt);
    wfs(i).gsdepth     = 
    wfs(i).lgs_focus_alt = 0.; // will be set by shwfs_comp_lgs_defocuses()
  }
  // and this to compute the defocus vector associated with the altitudes:
  shwfs_comp_lgs_defocuses,indgen(6);
  // re-compute kernels as gsdepth is likely to have changed after the fit:
  for (i=1;i<=6;i++) shwfs_init_lgs_kernels,i;
  // possibly sync children threads:
  sync_wfs;
 */
{
  alt = alt/1e3;
  delta = max(alt)-min(alt);
  a = _(1,array(1.0,npt),span(min(alt)+0.3*delta,max(alt)-0.3*delta,npt));
  r= lmfit(foo_lgs_profile,alt,a,amp,tol=1e-12);
  if (sim.debug) {
    plot,amp,alt;
    plg,foo_lgs_profile(alt,a),alt,color="red";
    pltitle_vp,"LGS profile (fg) and fit (red)";
    xytitles_vp,"Altitude [km]","Na density [fraction total]",[0.010,0.015];
  }
  fitamp = abs(a(2:npt+1));
  fitalt = clip(a(npt+2:),min(alt),max(alt))*1e3;
  fitdepth = a(1)*2.35*1e3;
  if (sim.verbose>1) {
    write,"Na profile fit results:";
    write,format="Depth = %f\n",a(1)*2.35;
    write,format="%s: ","altitudes"; fitalt;
    write,format="%s: ","amplitudes"; fitamp;
  }
}

//----------------------------------------------------
func foo_lgs_profile(x,a)
{
  // x vector of altitude
  npt = (numberof(a)-1)/2;
  sig = a(1);
  amp = abs(a(2:npt+1));
  pos = clip(a(npt+2:),min(x),max(x));
  res = x*0;
  for (i=1;i<=npt;i++) {
    res += amp(i)*exp(-((x-pos(i))/sqrt(2.)/sig)^2.);
  }
  return res;
}

//----------------------------------------------------
func shwfs_comp_lgs_defocuses(ns)
/* DOCUMENT shwfs_comp_lgs_defocuses,ns
 * This routine computes the defocus coefficients wfs._lgs_defocuses,
 * used by sh_wfs(), from the user-defined wfs.lgs_prof_alt vector.
 * Example:
 * wfs(1).lgs_prof_alt = &float([90,92,95,100,102]*1e3]);
 * shwfs_comp_lgs_defocuses,1;
 * add sync_wfs if using svipc.
 * Note that ns can be a vector:
 * shwfs_comp_lgs_defocuses,indgen(6);
 */
/* tests on 2012oct09:
 * diam=25, off-axis = 11.7m, delta_alt = 10km, sep should be 11.7/4.=2.92
 * lambda=0.65, pixsize=0.2145, sep = 2.82 < yeah.
 * lambda=1.0, pixsize=0.1980, sep = 2.76 < yeah
 * lambda=2.0, pixsize=0.2640, sep =2.71
 * I think we can say it works. let's try for larger pixsize:
 * lambda=2.0, pixsize=0.3960, sep=2.76, good.
 */
{
  extern wfs;
  
  for (i=1;i<=numberof(ns);i++) {
    // ok, so according to my calculations, we have:
    // a4[m] = D^2/(16*sqrt(3)) * 1/RoC
    // a4[rd] = a4[m]*2*pi/lambda = D^2/(16*sqrt(3)) * 2*pi/lambda * 1/RoC
    // and f = RoC/2
    if ((*wfs(ns(i)).lgs_prof_alt==[])||(allof(*wfs(ns(i)).lgs_prof_alt==0.))) {
      if (sim.verbose>0) \
        write,format="shwfs_comp_lgs_defocuses: wfs(%d).lgs_prof_alt undefined\n",ns(i);
      wfs(ns(i))._lgs_defocuses = &float([0.]);
      return;
    }
    tmp = tel.diam^2./(16*sqrt(3.)) * 2*pi/wfs(ns(i)).lambda/1e-6;
    a4s = tmp * 1./ *wfs(ns(i)).lgs_prof_alt;
    if (wfs(ns(i)).lgs_focus_alt==0) {
      wfs(ns(i)).lgs_focus_alt = sum(*wfs(ns(i)).lgs_prof_alt * (*wfs(ns(i)).lgs_prof_amp))/sum(*wfs(ns(i)).lgs_prof_amp);
    }
    a4cur = tmp * 1./ wfs(ns(i)).lgs_focus_alt;
    a4s = a4s-a4cur;
    wfs(ns(i))._lgs_defocuses = &float(a4s);
  }
}  

//----------------------------------------------------
func comp_turb_lgs_kernel(ns,init=)
/* DOCUMENT comp_turb_lgs_kernel(ns,init=)
   Compute the LGS object kernel correspoding to uplink seeing.
   takes 120 microsecond like this, probably not worth xferring to C.
 */
{
  extern wfs;

  llt_pscreen_dim = 256;

  if (init) {

    // first look up if this is not already an existing LLT:
    if (ns>1) {
      for (i=1;i<ns;i++) {
        // if (i==ns) continue;
        if (allof(wfs(ns).LLTxy==wfs(i).LLTxy)&&(wfs(ns).LLTr0==wfs(i).LLTr0)) {
          wfs(ns)._LLT_use = i;
          return;
        }
      }
    } 
    wfs(ns)._LLT_use = ns;
    // we didn't find a matching existing LLT. Let's init this one

    // Find out how to scale pixel size etc for the LLT:
    // first, the created turbulent kernel will be used in 
    // _shwfs_phase2spots to convolve an ximage array of size _nx4fft.
    // however, the pixel size in this array is set by:
    // quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/wfs(ns)._sdim;
    // now if f1 is size in physical units [m] of input phase array,
    // p1 the pixels size [m], NxN the size of the array,
    // and f2 and p2 the corresponding field size and pixel size in the
    // image plane, we have: f2 = lambda/p1 and p2=f2/N, thus p2 = lambda/(p1*N)
    // thus p1 = lambda/(p2*N)
    // ps = wfs(ns).lambda/wfs(ns)._nx4fft/quantumPixelSize/4.848;
    // which means, by the way:
    // ps = wfs(ns).lambda/wfs(ns)._nx4fft/quantumPixelSize/4.848;
    //    = wfs(ns).lambda/wfs(ns)._nx4fft/4.848/wfs(ns).lambda*(tel.diam/sim.pupildiam)*4.848*sdim
    //    = (tel.diam/sim.pupildiam)*(wfs(ns)._sdim/wfs(ns)._nx4fft)
    ps = (tel.diam/sim.pupildiam)*(float(wfs(ns)._sdim)/wfs(ns)._nx4fft);
    // now build up the "pupil":
    xy = dist(wfs(ns)._nx4fft)*ps;
    pup = exp(-8*(xy/wfs(ns).LLT1overe2diam)^2.);
    pup *= (xy<(wfs(ns).LLTdiam/2.));
    wfs(ns)._LLT_pupil = &float(pup);

    wfs(ns)._LLT_phase = &float(pup*0.0f); // just allocate space

    // do we have a phase screen?
    if (!wfs(1)._LLT_pscreen_name) { // we don't, let's look for it
      llt_pscreen = YAO_SAVEPATH+parprefix+"_llt_screen";
      if (noneof(findfiles(llt_pscreen+".fits"))) {
        // then we have to create it:
        if (sim.verbose) write,format="Creating %s.fits\n",llt_pscreen;
        create_phase_screens,llt_pscreen_dim,llt_pscreen_dim, \
                       prefix=llt_pscreen,no_ipart=1,silent=1;
      }
      wfs(1)._LLT_pscreen_name = llt_pscreen+".fits";
    }
    // now we can read it:
    pscreen = yao_fitsread(wfs(1)._LLT_pscreen_name);
    // scale it (done below)
    // pscreen *= (ps/wfs(ns).LLTr0)^(5./6.);
    // pad it:
    pscreen = _(pscreen,pscreen(,1:wfs(ns)._nx4fft));
    pscreen = transpose(pscreen);
    pscreen = _(pscreen,pscreen(,1:wfs(ns)._nx4fft));      
    // write,"TODO: Check normalization of llt phase screen !";
    wfs(1)._LLT_pscreen = &float(pscreen);

    if (!wfs(ns)._nx4fft) {
      // then the WFS has likely not been initialize.
      // we can't do anything until it is:
      error,swrite(format="WFS#%d has not been initialized",ns);
    }

    // this is the starting position for lower left phase to be extracted
    wfs(ns)._LLT_pos = float([1.,1 + (ns-1)*llt_pscreen_dim/nwfs]);

    return;
  } // end init
  // Compute turbulent kernel:
  // again, check that we have to:
  // below, assumes that the other one has already been computed, but
  // it should be if the user is following the ns order.
  if (wfs(ns)._LLT_use!=ns) return *wfs(wfs(ns)._LLT_use)._LLT_kernel;

  // we'll do the scaling here, so that the user can change the r0 on the fly
  ps = (tel.diam/sim.pupildiam)*(float(wfs(ns)._sdim)/wfs(ns)._nx4fft);
  scal = (ps/wfs(ns).LLTr0)^(5./6.);

  // Else we need to compute it:
  // v bar:
  wspeed = sum((*atm.layerfrac)*(*atm.layerspeed));
  // displacements w.r.t. last time we were called:
  wfs(ns)._LLT_pos += wspeed*loop.ittime*[1,float(wfs(ns)._nx4fft)/llt_pscreen_dim];

  // wrap if needed:
  wfs(ns)._LLT_pos = wfs(ns)._LLT_pos % llt_pscreen_dim;

  xshifts = float(wfs(ns)._LLT_pos(1)+indgen(wfs(ns)._nx4fft)-1);
  ishifts = int(xshifts); xshifts -= ishifts;
  yshifts = float(wfs(ns)._LLT_pos(2)+indgen(wfs(ns)._nx4fft)-1);
  jshifts = int(yshifts); yshifts -= jshifts;
  d = int(dimsof(*wfs(1)._LLT_pscreen));

  *wfs(ns)._LLT_phase *= 0.0f;

  err = _get2dPhase(wfs(1)._LLT_pscreen,d(2),d(3),1,
                    wfs(ns)._LLT_phase,wfs(ns)._nx4fft,wfs(ns)._nx4fft,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  psf = calc_psf_fast(*wfs(ns)._LLT_pupil,*wfs(ns)._LLT_phase,scale=scal,noswap=1);

  wfs(ns)._LLT_kernel = &psf;

  return psf;
}
