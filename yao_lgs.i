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
  plot,amp,alt;
  plg,foo_lgs_profile(alt,a),alt,color="red";
  write,format="depth = %f\n",a(1)*2.35;
  fitamp = abs(a(2:npt+1));
  fitalt = clip(a(npt+2:),min(alt),max(alt))*1e3;
  fitdepth = a(1)*2.35*1e3;
  fitamp;
  fitalt;
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
