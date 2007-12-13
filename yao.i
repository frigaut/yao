/*
 * yao.i
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: yao.i,v 1.3 2007-12-13 16:06:58 frigaut Exp $
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
 * Initial release F.Rigaut, June 2002.
 * see Release notes in README
 * all documentation at http://www.maumae.net/yao/aosimul.html
 *
 * version 4.1.1: 2007dec13:
 *  - upgrade/update of Makefile
 *
 * version 4.1: 2007dec12:
 *  - cleaned up all things for the debian release and a fresh pkg_mngr
 *    release
 *  - cleaned up the LICENSE issue (GPLv2)
 *  - suppressed warning message if fftw_wisdom.dat not found. In my
 *    experience, on modern hardware and for fftw3, this does not improve
 *    speed. One can still generate and use one if wanted (run
 *    init_fftw_wisdom)
 *  - produced man page (doc/yao.1)
 *
 * version 4.0: 2007jun03: 
 *  - Added all GTK GUI. No further modifs in yao,
 *    hence yao can still be used safely without a gui (scripts,
 *    scripts).
 *  - migrated to a flat directory structure for yao. There is no
 *    more yao_fftw subdirectory
 *  - This is now meant to be installed in i0, so the path to the yao
 *    required functions has been changed from yao/file.i to file.i
 * 
 * version 3.81, 2007jun02: added gui_message
 * 
 * version 3.80, 2007apr19: now 64 bits safe !
 * 
 * version 3.75, 2007feb22: added angle offset for curvature wfs and dms
 * 
 * 2006mar01: corrected bug. reset put command to 0,
 * while it needed to zero dm._command 
 *
 * Old CVS Log: aosimul.i
 * Revision 1.26  2004/10/19 00:18:04  frigaut
 * replaced strJoin -> strjoin and the like for compat with new string_ops.i
 *
 * Revision 1.25  2004/10/18 21:57:02  frigaut
 * fixed the mess that I had the swap screen update committed from lono
 * and I just committed the coupling update with a version that did
 * not contain the swap screen update.
 * Now I believe this version contains both.
 *
 * Revision 1.23  2004/09/29 03:59:37  frigaut
 * aosimul modified to accomodate screen swaps
 * turbulence.i modified for outer scale
 * aoutil.i modified to not impose the gain value, allowing 0.
 * it prints out a simple warning if loop.gain=0 (allow open
 * loop runs).
 *
 * Revision 1.22  2004/09/14 04:32:56  frigaut
 * several modifs to do with the creation of turbulent phase screens
 * - implemented cosf and sinf which take and return floats to save space
 * - started coding generatePhaseWithL0 in turbulence.i, not finished
 * - modif YORICKEXE in makefiles, just "yorick" did not do it on osx
 * - modifs ok for both veclib and fftw implementations
 *
 * Revision 1.21  2004/09/14 00:27:34  frigaut
 * refresh/update from web site files
 *
 * Revision 1.20  2004/09/11 01:28:17  frigaut
 * - added a printout of the yao version used in yao_fast.i's
 * - added a make check-plug in Makefiles
 * - updated version and date in aosimul.i
 *
 * Revision 1.19  2004/09/11 00:21:37  frigaut
 * modified to include yao and yao_fast.i at parse
 * This is the first version compatible with plugins.
 *
 * Revision 1.18  2004/08/02 22:58:02  frigaut
 * modified tic and tac to accomodate multiple counters
 * modified aoloop to use it.
 * now the "Time Left" functionality works again.
 *
 * Revision 1.17  2004/08/02 08:19:18  frigaut
 * global update to be sure I did not forget a file
 *
 * Revision 1.16  2004/08/02 07:59:05  frigaut
 * Updated version number and date in file header.
 *
 * Revision 1.15  2004/08/02 07:54:11  frigaut
 * Fixed a bug in ShWfsInit noted by Miska when running odd number of
 * subaperture shack-hartmann systems. There was a division by zero
 * in the call to atan to determine the angle of the kernel. fixed
 * by calling atan(y,x) instead of atan(y/x)
 *
 * Revision 1.14  2004/08/02 07:21:25  frigaut
 * Corrected a nasty bug in the final adjustment of yposvec in
 * getTurbPhaseInit (used xmargin instead of ymargin)
 *
 * Revision 1.13  2004/08/02 07:11:13  frigaut
 * Added call to getTurbPhaseInitCheckOverflow in getTurbPhaseInit
 *
 * Revision 1.12  2004/07/29 04:06:50  frigaut
 * added cvs dollar Log in header
 *
*/

extern aoSimulVersion, aoSimulVersionDate;
aoSimulVersion = yaoVersion = aoYaoVersion = "4.1.1";
aoSimulVersionDate = yaoVersionDate = aoYaoVersionDate = "2007dec13";

write,format=" Yao version %s, Last modified %s\n",yaoVersion,yaoVersionDate;


plug_in,"yao";

require,"yao_utils.i";
require,"yao_fast.i";
require,"aoutil.i";
require,"yao_gui.i";
require,"utils.i";
require,"newfits.i";
require,"yao_util.i";
require,"turbulence.i";
require,"plot.i";  // in yorick-yutils

func null (arg,..) { return 0; }

pyk_error = pyk_info = pyk_warning = null;
gui_message = gui_progressbar_frac = gui_progressbar_text = null;  // by default. can be redefined by gui routine

//----------------------------------------------------
func compDmShape(nm,command,extrap=)
/* DOCUMENT func compDmShape(nm,command)
   compute the DM shape from the command vector
   branch over _dmsum or _dmsumelt according to case.
   SEE ALSO:
*/
{
  n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = int(n2-n1+1);
  sphase = array(float,[2,nxy,nxy]);

  if (!is_set(extrap)) {

    if (dm(nm).elt == 1) { //use fast dm shape computation
      _dmsumelt, dm(nm)._def, dm(nm)._eltdefsize, dm(nm)._eltdefsize, int(dm(nm)._nact),
        dm(nm)._i1, dm(nm)._j1, command, &sphase,nxy,nxy;
    } else { // use standard
      _dmsum, dm(nm)._def, nxy, nxy, dm(nm)._nact, command, &sphase;
    }
    
  } else { // extrapolated actuators

    if (dm(nm).elt == 1) { //use fast dm shape computation
      _dmsumelt, dm(nm)._edef, dm(nm)._eltdefsize, dm(nm)._eltdefsize, int(dm(nm)._enact),
        dm(nm)._ei1, dm(nm)._ej1, command, &sphase,nxy,nxy;
    } else { // use standard
      _dmsum, dm(nm)._edef, nxy, nxy, dm(nm)._enact, command, &sphase;
    }
    
  }
  
  return sphase;
}
//----------------------------------------------------
func controlScreen(i,init=)
{
  local x0,y0,x,y;
  y0 = 0.84;
  x0 = 0.165;
  ygstep = 0.03;
  mygstep = 0.022;
  sygstep = 0.017;
  csize = 10;

  /*  y0 = 0.86;
      x0 = 0.05;
      ygstep = 0.037;
      mygstep = 0.030;
      sygstep = 0.025;
      csize = 18;
  */
  
  if (is_set(init)) {
    //    dheight = (nwfs+target._nlambda+1)*sygstep + 10*ygstep;
    //    dheight = long(dheight/0.77*570);
    winkill,2;  // make sure
    window,2,width=390,height=350,style="letter.gs",wait=1,dpi=70;
    //    window,2,width=570,height=dheight,style="letter.gs",wait=1,dpi=70;
    limits,0.,1.,0.,1.;
    progressBar,0,pbid1,init=[x0,y0-ygstep,x0+0.25,y0-ygstep+0.015];
    window,0;
    return;
  }

  window,2;
  fma;
  y = y0;
  if (!remainingTimestring) remainingTimestring="N/A";
  plt,sim.name+"  /  Time left: "+remainingTimestring,x0,y,tosys=1,height=csize,justify="LA";
  y -= ygstep;
  it = swrite(format="%d/%d iterations",i,loop.niter);
  plt,it,x0+0.27,y,tosys=1,height=csize,justify="LA";
  y -= ygstep;
  progressBar,i*100./loop.niter,pbid1;
  
  //  s = "WFS     rms     Tip     Tilt    ";
  s = "WFS        Tip    Tilt   ";
  if (anyof(wfs.correctUpTT)) {
    s = s+"UpTip  UpTilt   ";
    if (anyof(wfs.centGainOpt)) {
      s = s+"CentG   ";
    }
  }
  plt,s,x0,y,tosys=1,height=csize,justify="LA",font="courier";
  y -= mygstep;
  for (i=1;i<=nwfs;i++) {
    s = swrite(format="%3i     % 6.3f  % 6.3f  ",i,
               wfs(i)._lastvalidtt(1),wfs(i)._lastvalidtt(2));
    if (anyof(wfs.correctUpTT)) {
      if (wfs(i).correctUpTT) {
        s = s+swrite(format="% 6.3f  % 6.3f  ",
                     wfs(i)._upttcommand(1),wfs(i)._upttcommand(2));
      } else {
        s = s+"   -       -       ";
      }
      if (anyof(wfs.centGainOpt)) {
        if (wfs(i).centGainOpt) {
          s = s+swrite(format="% 6.3f  ",wfs(i)._centroidgain);
        } else {
          s = s+"-       ";
        }
      }
    }
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font="courier";
    y -= sygstep;
  }

  // print out TTM command
  wtt = where(dm.type == "tiptilt");
  if (numberof(wtt) != 0) {
    wtt = wtt(0);
    s = swrite(format="TT Mirror: % 6.3f  % 6.3f",(*dm(wtt)._command)(1),
               (*dm(wtt)._command)(2));
    y -= sygstep;
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font="courier";
    y -= ygstep;
  }

  // print out Strehls:
  plt,"Strehl Ratio  Lambda     Avg     rms     Min     Max",x0,y,tosys=1,height=csize,
    justify="LA",font="courier";
  y -= sygstep;
  
  for (jl=1;jl<=target._nlambda;jl++) {
    strehlle = imav(max,max,,jl)/sairy/(niterok+1e-5);
    s = swrite(format="Since Start   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
               (*target.lambda)(jl),strehlle(avg),strehlle(rms),
               min(strehlle),max(strehlle));
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font="courier";
    y -= sygstep;
    if (jl == target._nlambda) {
      strehlse = im(max,max,)/sairy;
      s = swrite(format="Instantaneous %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
                 (*target.lambda)(jl),strehlse(avg),strehlse(rms),
                 min(strehlse),max(strehlse));
      plt,s,x0,y,tosys=1,height=csize,justify="LA",font="courier";
      y -= sygstep;
    }
  }    
      

  window,0;
}

//----------------------------------------------------

func CurvWfs(pupil,phase,ns,init=,disp=,silent=)
/* DOCUMENT  CurvWfs(pupil,phase,ns,disp=)
   This function computes the signal from a Curvature WFS for a 
   given phase and pupil input and WFS config.
*/
{
  if (is_void(ns)) {ns=1;} // default sensor#1 for one WFS work
  
  size	     = sim._size;
  dimpow2   = int(log(size)/log(2));

  if (init == 1) {
    if ( (sim.verbose>=1) && (!is_set(silent)) ) {write,"> Initializing CurvWfs\n";}
    fratio= 60.;
    defoc = (pi*wfs(ns).lambda*1e-6/(sim._size^2.*(tel.diam/sim.pupildiam)^2.))*
      eclat(dist(sim._size)^2.);
    x	= fratio*tel.diam*(fratio*tel.diam-wfs(ns).l)/wfs(ns).l;
    defoc= defoc*x;
    wfs(ns)._cxdef= &(float(cos(defoc))); 
    wfs(ns)._sxdef= &(float(sin(defoc)));
    wfs(ns)._tiltsh = &(float(defoc*0.));
    wfs(ns)._fimage = &(float(defoc*0.));
    wfs(ns)._fimage2 = &(float(defoc*0.));

    // Work out the total NUMBER OF PHOTONS per sample
    // from star
    if (wfs(ns).gsalt == 0) {

      wfs(ns)._nphotons = gs.zeropoint*10^(-0.4*wfs(ns).gsmag)*
        wfs(ns).optthroughput*                 // include throughput to WFS
        loop.ittime;                           // per iteration time

    } else { // we are dealing with a LGS

      telsurface = pi/4.*tel.diam^2.(1-tel.cobs^2.)*1e4; // in cm^2

      wfs(ns)._nphotons = gs.lgsreturnperwatt*wfs(ns).laserpower*
        telsurface*loop.ittime;

    }
    // from sky, over field stop
    wfs(ns)._skynphotons = gs.zeropoint*10^(-0.4*wfs(ns).skymag)*
      wfs(ns).optthroughput*                 // include throughput to WFS
      loop.ittime*pi/4*wfs(ns).fieldstopdiam^2.;

    if ( (sim.verbose>=1) && (!is_set(silent)) ) {
      write,format="NPhotons/iter from star = %f\n",wfs(ns)._nphotons;
      write,format="NPhotons/iter from sky  = %f\n",wfs(ns)._skynphotons;

    }
    return defoc;
  }

  mesvec = array(float,wfs(ns)._nsub);
  if (typeof(phase) != "float") {
    print,"Phase was not float";
    phase = float(phase);
  }

  phasescale = float(2*pi/wfs(ns).lambda);   // wfs.lambda in microns

  err = _cwfs( &pupil, &phase, phasescale, wfs(ns)._tiltsh, wfs(ns)._cxdef,
               wfs(ns)._sxdef, dimpow2, wfs(ns)._sind, wfs(ns)._nsind,
               wfs(ns)._nsub, wfs(ns)._fimage, wfs(ns)._fimage2,
               wfs(ns)._nphotons, wfs(ns)._skynphotons, float(wfs(ns).ron),
               float(wfs(ns).darkcurrent*loop.ittime), int(wfs(ns).noise), &mesvec);

  return mesvec;
}

//----------------------------------------------------

func mcaoRayleigh(nwfs,xsubap,ysubap,zenith=,fov=,aspp=)
{
  as2rd = dtor/3600.;
  cobs  = 0;
  
  // position of GS/WFS in arcsec:
  w =  where(wfs._gsalt > 0);
  //  xwfs = [0,-30.0,30.0,30.0,-30.0];
  //  ywfs = [0,30.0,30.0,-30.0,-30.0];
  xwfs = wfs(w).gspos(1,);
  ywfs = wfs(w).gspos(2,);
  nbeams = numberof(xwfs);
  
  if (!is_set(zenith)) {zenith = 0.;}
  zenith = zenith*dtor;

  if (!is_set(fov))  { fov  = 2;}  // FoV (arcsec)
  if (!is_set(aspp)) { aspp = 1;} // arcsec per pixel

  // I have not implemented the 4 following parameters in the parfile.
  // one has to fill it b hand in the code for now!!!!!!!
  beamdiameter = 0.3; // fwhm of gaussian laser beam in meter

  laserlambda  = wfs(nwfs).lambda*1e-6; //589e-9;
  diamsubap    = tel.diam/wfs(nwfs).shnxsub; // diameter of a subaperture [m]
  r0           = (tel.diam/atm.dr0at05mic)*cos(zenith)^0.6;
  seeing       = laserlambda/r0/4.848e-6;
  spotsize     = 1.0; // irrelevant for rayleigh in this code.
  d            = 1e-2; // linear size of aperture for flux (/cm^2)

  altsod       = wfs(nwfs)._gsalt;
  fwhmsod      = wfs(nwfs)._gsdepth;

  // definitions of the image array to return to caller
  dim = long(ceil(fov/aspp));
  imrayl = array(float,[2,dim,dim]);
  imstar = array(float,[2,dim,dim]);

  // looking at GS #n with WFS #n -> some angle
  phin   = sqrt(xwfs(nwfs)^2.+ywfs(nwfs)^2.)*as2rd;
  thetan = atan(xwfs(nwfs),ywfs(nwfs));
  // range vector to sample h (for Rayleigh only)
  //  rvec = spanl(2000,altsod,150);
  rvec = spanl(1000,altsod,50);

  // definitions for sodium GS:
  // range vector for sodium:
  rvecsod = span(altsod-fwhmsod,altsod+fwhmsod,20);

  // loop on beams
  for (beam=1;beam<=nbeams;beam++) {

    // define "beam" position angles
    phib = sqrt(xwfs(beam)^2.+ywfs(beam)^2.)*as2rd;
    thetab = atan(xwfs(beam),ywfs(beam));

    // loop on altitude for RAYLEIGH
    for (i=1;i<=numberof(rvec)-1;i++) {

      // range
      r = (rvec(i+1)+rvec(i))/2.;
      // altitude
      z = r*cos(zenith);
      // delta range
      deltar = rvec(i+1)-rvec(i);

      // fit to the lidar equation:
      sbnr = 16.12*exp(-z*0.14177e-3)*1e-6;
      // total number of photons received per cm^2 and period (800Hz)
      rayleigh = wfs(beam).laserpower/(6.62e-34*3e8/589e-9)*sbnr*d^2/(4*pi*r^2.)*
        deltar*loop.ittime*wfs(nwfs).optthroughput;

      // compute position angles of the beam center for this altitude
      alpha = sin(phib)*cos(thetab)-sin(phin)*cos(thetan) - xsubap/r + xsubap/altsod;
      beta  = sin(phib)*sin(thetab)-sin(phin)*sin(thetan) - ysubap/r + ysubap/altsod;
      alpha /= 4.848e-6; // in arcsec
      beta  /= 4.848e-6;
      if (sim.verbose > 1) {
        write,format="range = %f, alpha = %f, beta = %f\n",r,alpha,beta;
      }
      
      // generate gaussian:
      fwhm = beamdiameter/r/as2rd; // in arcsec
      // blur due to finite range of this layer of rayleigh scatering
      blur = (altsod - r)/altsod*diamsubap/r/4.848e-6;
      fwhm = sqrt(fwhm^2. + seeing^2.+blur^2.);

      g = makegaussian(dim,fwhm/aspp,xc=0.5+dim/2+alpha/aspp,yc=0.5+dim/2+beta/aspp,norm=1);
      if (cobs) {
        fwhmobs = fwhm/2;
        fwhmobs = sqrt(fwhmobs^2. + seeing^2.);
        gobs = makegaussian(dim,fwhmobs/aspp,xc=0.5+dim/2+alpha/aspp,yc=0.5+dim/2+beta/aspp);
        gobs = gobs*max(g); g = g-gobs;
      }

      // add contribution from this slab. Correctly normalized
      // in photons/cm^2/exptime/laser_power/telescope+system_throughput
      imrayl += g*rayleigh;
    }

    // loop on altitude for SODIUM

    // Na return detected on WFSCCD in ph/cm^2/exptime/laser_power/telescope+system_throughput:
    nPhotonFromSodStar = wfs(beam).laserpower*gs.lgsreturnperwatt*loop.ittime*
      cos(zenith)*wfs(nwfs).optthroughput;
    sodprofile = exp(-((rvecsod-altsod)/(fwhmsod/2.))^4.);
    sodprofile = sodprofile/sum(sodprofile)*nPhotonFromSodStar;

    for (i=1;i<=numberof(rvecsod)-1;i++) {

      // range
      r = rvecsod(i);
      sodium  = sodprofile(i);

      // compute position angles of the beam center for this altitude
      alpha = sin(phib)*cos(thetab)-sin(phin)*cos(thetan) - xsubap/r + xsubap/altsod;
      beta  = sin(phib)*sin(thetab)-sin(phin)*sin(thetan) - ysubap/r + ysubap/altsod;
      alpha /= 4.848e-6; // in arcsec
      beta  /= 4.848e-6;
      if (sim.verbose > 1) {
        write,format="range = %f, alpha = %f, beta = %f\n",r,alpha,beta;
      }
      
      // generate gaussian:
      fwhm = beamdiameter/r/as2rd; // in arcsec
      // here I'll have to take the defocus of range r into account
      fwhm = sqrt(fwhm^2. + seeing^2.);

      g = makegaussian(dim,fwhm/aspp,xc=0.5+dim/2+alpha/aspp,yc=0.5+dim/2+beta/aspp,norm=1);

      // add contribution from this slab. Correctly normalized
      // in photons/cm^2/exptime/laser_power/telescope+system_throughput
      imstar += g*sodium;
    }
    if (sim.verbose > 1) {
      fma;
      pli,imrayl+imstar,-dim/2*aspp,-dim/2*aspp,+dim/2*aspp,+dim/2*aspp;
      pause,20;
    }
  }
  return [imrayl,imstar];
}  

//----------------------------------------------------

func ShWfsInit(pupsh,ns,silent=,imat=)
{
  extern initkernels;
  
  initkernels = array(1n,nwfs);

  if (is_void(ns)) {ns=1;} // default to wfs#1 for one WFS work.

  if (typeof(pupsh) != "float") {error,"pupsh was not float !";}

  pupd	     = sim.pupildiam;
  size	     = sim._size;
  nxsub	     = wfs(ns).shnxsub(0);
  subsize    = pupd/nxsub;
  fracsub    = wfs(ns).fracIllum;
  sdim       = long(2^ceil(log(subsize)/log(2)+1));
  sdimpow2   = int(log(sdim)/log(2));

  wfs(ns)._centroidgain = 1.f;

  //====================================================================
  // WORK OUT THE NUMBER OF PHOTONS COLLECTED PER SUBAPERTURE AND SAMPLE
  //====================================================================

  telSurf  = pi/4.*tel.diam^2.;
  // from the guide star (computed here as used in wfsCheckPixelSize):
  if (wfs(ns).gsalt == 0) {

    wfs(ns)._nphotons = gs.zeropoint*10^(-0.4*wfs(ns).gsmag)*
      wfs(ns).optthroughput*                 // include throughput to WFS
      (tel.diam/wfs(ns).shnxsub)^2./telSurf* // for unobstructed subaperture
      loop.ittime;                           // per iteration

  } else {  // we are dealing with a LGS

    wfs(ns)._nphotons = gs.lgsreturnperwatt*  // detected by WFS
      wfs(ns).laserpower*                     // ... for given power
      wfs(ns).optthroughput*                  // include throughput to WFS
      (tel.diam/wfs(ns).shnxsub)^2.*1e4*      // for unobstructed subaperture
      loop.ittime;                            // per iteration

  }
  // see below for # of photons from sky

  //=============================================================
  // COMPUTE KERNEL TO CONVOLVE _SHWFS IMAGE FOR IMAT CALIBRATION
  //=============================================================

  // if there is no explicit request for an extended kernel
  // and we are not using LGS (or the depth = 0) then we disable
  // the kernel convolution is _shwfs to gain time (and accuracy)
  // this behavior is overriden in MultWfsIntMat anyway for the purpose
  // of computing the iMat with the correct kernel (to simulate for
  // the extended "seeing" spot
  if ((wfs(ns).kernel == 0) && (wfs(ns)._gsdepth == 0) && (!is_set(imat))) {
    wfs(ns)._kernelconv = 0n;
  } else {
    wfs(ns)._kernelconv = 1n;
  }

  quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;

  write,format="Dimension for optional amplitude mask: %d\n",2^sdimpow2;
  
  // reads out the amplitude mask for the subaperture:
  if (wfs(ns).submask) {
    // read the amplitude image
    tmp = fitsRead(wfs(ns).submask);
    // check that dims are OK:
    if (anyof(dimsof(tmp)!=[2,2^sdimpow2,2^sdimpow2])) {
      error,swrite(format="Bad dimensions for %s. Should be %d, found %d\n",
                   wfs(ns).submask,2^sdimpow2,dimsof(tmp)(2));
    }
    // check of mask is centered (can be a common cause of mistake):
    f1x = sum(tmp(1:2^(sdimpow2-1),));
    f2x = sum(tmp(2^(sdimpow2-1)+1:,));
    f1y = sum(tmp(,1:2^(sdimpow2-1)));
    f2y = sum(tmp(,2^(sdimpow2-1)+1:));
    if ((f1x!=f2x)||(f1y!=f2y)) {
      write,format="%s\n","\nWARNING!";
      write,format="%s\n","The SHWFS amplitude mask is not centered. This can create";
      write,format="%s\n","a bias in the slope calculation. A centered mask should be";
      write,format="%s\n","centered on the 4 central pixels of the mask, not on";
      write,format="%s\n","the (0,0) of the FFT transform. If you did not do that on ";
      write,format="%s\n\n","purpose, you should correct your mask.";
}
    // modify as required:
    wfs(ns)._submask = &(float(roll(tmp)));
    wfs(ns)._domask = 1l;
  } else wfs(ns)._domask = 0l;

  // factor 1.5 to crudely compensate for the fact that the spot is
  // tilt compensated at the focus of the lenslet

  if (is_set(imat)) { // we are in iMat computation. we want
    // the kernel FWHM = seeing + requested kernel fwhm (quadratically)
    dr0 = atm.dr0at05mic*(0.5/wfs(ns).lambda)^1.2/cos(gs.zenithangle*dtor)^0.6;
    fwhmseeing = wfs(ns).lambda/
      (tel.diam/sqrt(wfs(ns).shnxsub^2.+(dr0/1.5)^2.))/4.848;
    kernelfwhm = sqrt(fwhmseeing^2.+wfs(ns).kernel^2.);
  } else { // we are in regular aoloop. no further convolution to account
           // for seeing. However, we want to avoid dividing by zero
           // in makegaussian, so we floor fwhm:
    kernelfwhm = clip(wfs(ns).kernel,1e-8,);
  }

  if ( (sim.verbose >= 1) && (!is_set(silent)) ) {
    write,format="Kernel FWHM for the iMat calibration = %f\n",kernelfwhm;
  }

  tmp = eclat(makegaussian(sdim,kernelfwhm/quantumPixelSize,
                           xc=sdim/2+1,yc=sdim/2+1));
  tmp(1,1) = 1.; // this insures that even with fwhm=0, the kernel is a dirac
  tmp = tmp/sum(tmp);
  wfs(ns)._kernel = &(float(tmp));

  //==================================
  // COMPUTE ISTART AND JSTART VECTORS
  //==================================
  // they are vectors containing the start
  // indices (X and Y) of the subapertures

  is	= size/2+1-pupd/2;
  nsubs = nxsub*nxsub;
  istart = (indgen(nsubs)-1)/nxsub*subsize+is;
  jstart = ((indgen(nsubs)-1)%nxsub)*subsize+is;
  xsub = ((indgen(nsubs)-1)/nxsub+0.5)*tel.diam/nxsub-tel.diam/2.;
  ysub = (((indgen(nsubs)-1)%nxsub)+0.5)*tel.diam/nxsub-tel.diam/2.;

  //==========================================================
  // COMPUTE WHICH SUBAPERTURES ARE ENABLED (FLUX > THRESHOLD)
  //==========================================================

  fluxPerSub = array(float,nsubs);
  for (i=1;i<=nsubs;i++) {
    fluxPerSub(i) = sum(pupsh(istart(i):istart(i)+subsize-1,jstart(i):jstart(i)+subsize-1));
  }
  fluxPerSub = fluxPerSub/subsize^2.;
  // indices of the enabled subapertures: gind
  gind = where(fluxPerSub > fracsub);
  istart = istart(gind);
  jstart = jstart(gind);
  xsub = xsub(gind);
  ysub = ysub(gind);
  fluxPerSub  = fluxPerSub(gind);

  // stuff some of wfs structure for WFS "ns":
  wfs(ns)._istart = &(int(istart-1)); // -1n 'cause C is 0 based
  wfs(ns)._jstart = &(int(jstart-1));
  wfs(ns)._x = &(xsub);
  wfs(ns)._y = &(ysub);
  wfs(ns)._nsub = int(numberof(gind));

  // compute other setup variables for _shwfs:

  //==========================================================================
  // COMPUTE WFS IMAGE KERNELS: SUBAPERTURE DEPENDANT. USED FOR LGS ELONGATION
  //==========================================================================

  if (wfs(ns)._kernelconv != 0n) {
    // if kernelconv is 0, then the _shwfs routine does not use wfs._kernels,
    // so no need to compute it.
    quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
    xy = (indices(sdim)-sdim/2.-1)*quantumPixelSize;  // coordinate array in arcsec

    if (sim.verbose >= 1) {write,"Pre-computing Kernels for the SH WFS";}

    kall = [];

    for (l=1; l<=wfs(ns)._nsub; l++) {
      // for each subaperture, we have to find the elongation and the angle.
      xsub = (*wfs(ns)._x)(l); ysub = (*wfs(ns)._y)(l);
      xllt = wfs(ns).LLTxy(1); yllt = wfs(ns).LLTxy(2);
      d = sqrt((xsub-xllt)^2. +(ysub-yllt)^2.);
      if (d == 0) {
        elong = ang = 0.;
      } else {
        elong = atan((wfs(ns)._gsalt+wfs(ns)._gsdepth)/d)-atan(wfs(ns)._gsalt/d);
        elong /= 4.848e-6; // now in arcsec
        ang = atan(ysub-yllt,xsub-xllt); // fixed 2004aug02
      }
      angp90 = ang+pi/2.;

      expo = 4.;
      alpha = elong/(2.*log(2)^(1/expo));
      beta = 2*quantumPixelSize/(2.*log(2)^(1/expo));
      alpha = clip(alpha,0.5*quantumPixelSize,);
      beta = clip(beta,0.5*quantumPixelSize,alpha);
      
      tmp  = exp(-((cos(ang)*xy(,,1)+sin(ang)*xy(,,2))/alpha)^expo);
      tmp *= exp(-((cos(angp90)*xy(,,1)+sin(angp90)*xy(,,2))/beta)^expo);
      tmp = tmp/sum(tmp);
      grow,kall,(eclat(tmp))(*);
      
      s2 = tel.diam/nxsub*0.9/2.;
      //if (sim.debug >= 2) {fma; pli,tmp,xsub-s2,ysub-s2,xsub+s2,ysub+s2;}
    }
    //if (sim.debug >=2) {hitReturn;}

    wfs(ns)._kernels = &(float(kall));
    wfs(ns)._kerfftr = &(float(kall*0.));
    wfs(ns)._kerffti = &(float(kall*0.));
    kall = [];
  }
  
  //==========================================================================
  // COMPUTE WFS IMAGE KERNELS: SUBAPERTURE DEPENDANT. USED FOR LGS ELONGATION
  //==========================================================================

  rayleighflux = array(0.0f,wfs(ns)._nsub);
  sodiumflux   = array(1.0f,wfs(ns)._nsub);

  if (wfs(ns).rayleighflag == 1n) {
    quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
    xy = (indices(sdim)-sdim/2.-1)*quantumPixelSize;  // coordinate array in arcsec

    if (sim.verbose > 0) {write,"Pre-computing Rayleigh for the SH WFS";}

    kall = [];

    rayfname = parprefix+"-rayleigh-wfs"+swrite(format="%d",ns)+"-zen"+
      swrite(format="%d",long(gs.zenithangle))+".fits";
    isthere = fileExist(rayfname);
    fov    = quantumPixelSize*sdim;
    aspp   = quantumPixelSize;


    if (isthere) {
      
      kall         = fitsRead(rayfname)*1e6;
      rayleighflux = fitsRead(rayfname,hdu=1);
      sodiumflux   = fitsRead(rayfname,hdu=2);

    } else {

      for (l=1; l<=wfs(ns)._nsub; l++) {
        xsub = (*wfs(ns)._x)(l); ysub = (*wfs(ns)._y)(l);
        tmp = mcaoRayleigh(ns,ysub,xsub,fov=fov,aspp=aspp,zenith=gs.zenithangle);
        rayleighflux(l) = sum(tmp(,,1));
        sodiumflux(l)   = sum(tmp(,,2));
        tmp = transpose(tmp(,,1));
        // the switch of xsub <-> ysub and transpose are to accomodate the
        // C vs yorick 0 vs 1 start array index.
        //        tmp = tmp/sum(tmp);
        grow,kall,(eclat(tmp))(*);
      }
      fitsWrite,rayfname,kall;
      fitsWrite,rayfname,rayleighflux,append=1,exttype="image";
      fitsWrite,rayfname,sodiumflux,append=1,exttype="image";

    }

    if (sim.verbose > 0) {
      write,"From Rayleigh calculations:";
      write,format="Rayleigh / Sodium ratio in worse/best subapertures: %f, %f\n",
        max(rayleighflux/sodiumflux),min(rayleighflux/sodiumflux);
      write,format="Rayleigh flux varies between %f and %f /cm2/looptime\n",
        min(rayleighflux),max(rayleighflux);
      write,format="Sodium   flux varies between %f and %f /cm2/looptime\n",
        min(sodiumflux),max(sodiumflux);
    }
    
    wfs(ns)._rayleigh = &(float(kall));
    kall = [];
  }
  
  //================================
  // SUBAPERTURE SIZE AND PIXEL SIZE
  //================================
  wfsCheckPixelSize,ns,binindices,centroidw,printheader=(ns == 1),silent=silent;
  
  binxy = int(int(wfs(ns).npixels));

  // stuff some more of wfs structure for WFS "ns":
  wfs(ns)._binindices = &(int(binindices));
  wfs(ns)._binxy = binxy;
  wfs(ns)._centroidw = &centroidw;
  wfs(ns)._fimage = &(array(float,[2,nxsub*(binxy+1)+1,nxsub*(binxy+1)+1]));

  // 2004mar22: added a guard pixel for each subaperture for the display
  imistart = (istart-min(istart))/subsize*(binxy+1)+1; //same: starts @ 0 'cause C 0 based
  imjstart = (jstart-min(jstart))/subsize*(binxy+1)+1; //same
  
  wfs(ns)._imistart = &(int(imistart));
  wfs(ns)._imjstart = &(int(imjstart));
  wfs(ns)._fimnx = int(nxsub*(binxy+1)+1);
  wfs(ns)._fimny = int(nxsub*(binxy+1)+1);

  // This is the tilt to add to the input phase so that
  // the individual subaperture images fall in between
  // the pixels of the quadcell
  xy     = indices(sim._size);

  wfs(ns)._tiltsh = &(float(-64.*0.098174773*(xy(,,1)+xy(,,2))* \
                            0.5/sdim*wfs(ns).lambda/(2*pi)*(wfs(ns).shmethod == 2)));

  // This tilt array is intended to bring the spot back inbetween 4 pixels
  // instead of centered on dim/2+1 as a result of the regular FFT.
  // this is an achromatic factor of course that just depends on the
  // dimension of the array.
  // the lambda/2pi factor is thus to compensate the x by 2pi/lambda
  // in _shwfs, making tiltsh achromatic.


  // compute # of photons from the sky as above for guide star.
  // skymag is per arcsec, so we have to convert for
  // the subaperture size computed in wfsCheckPixelSize.

  wfs(ns)._skynphotons = gs.zeropoint*10^(-0.4*wfs(ns).skymag)*
    (tel.diam/wfs(ns).shnxsub)^2./telSurf*   // per subaperture
    loop.ittime*                             // per iteration
    wfs(ns).optthroughput*                   // include throughput to WFS
    (wfs(ns).pixsize*wfs(ns).npixels)^2.;    // for the actual subap. size

  if ( wfs(ns).skymag == 0) { wfs(ns)._skynphotons = 0.; } // if skymag not set

  if (sim.verbose == 2) {
    write,format="wfs(%d)._skynphotons = %f\n\n",ns,wfs(ns)._skynphotons;
  }
  
  // for guide star, total "useful" signal per subap.
  wfs(ns)._fluxpersub  = &(float(fluxPerSub*wfs(ns)._nphotons));
  if (sim.verbose > 0) {
    write,"From flux calculations:";
    write,format="Sodium flux varies between %f and %f /subaperture/looptime\n",
      min(*wfs(ns)._fluxpersub),max(*wfs(ns)._fluxpersub);
  }
  
  // for rayleigh, if any:
  wfs(ns)._raylfluxpersub  = &(*wfs(ns)._fluxpersub*float(rayleighflux/sodiumflux));
    
  // for sky:
  wfs(ns)._skyfluxpersub  = &(float(fluxPerSub*wfs(ns)._skynphotons));

  wfs(ns)._bias = &(wfs(ns).biasrmserror *
                    gaussdev(wfs(ns)._nsub*wfs(ns)._binxy*wfs(ns)._binxy));

  wfs(ns)._flat = &(1.0f + wfs(ns).flatrmserror *
                    gaussdev(wfs(ns)._nsub*wfs(ns)._binxy*wfs(ns)._binxy));

  wfs(ns)._bckgrdcalib = &(array(float,wfs(ns)._nsub*wfs(ns)._binxy*wfs(ns)._binxy));
  
  if (sim.verbose == 2) {
    write,format="Dark current wfs#%d / iter / pixel=%f\n",ns,
      float(wfs(ns).darkcurrent*loop.ittime);
  }

  // calibrate background images (have to run shwfs for that):
  wfs(ns)._bckgrdsub  = 1;
  wfs(ns)._bckgrdinit = 1;

  // call ShWfs for calibration of the background
  ShWfs,pupsh,pupsh*0.0f,ns;

  wfs(ns)._bckgrdinit = 0;
  
  return 1;
}
//----------------------------------------------------

func ShWfs(pupsh,phase,ns)
{
  if (is_void(ns)) {ns=1;} // default to wfs#1 for one WFS work.

  if (typeof(pupsh) != "float") {error,"pupsh was not float !";}
  if (typeof(phase) != "float") {error,"Phase was not float !";}
  
  pupd	     = sim.pupildiam;
  size	     = sim._size;
  nxsub	     = wfs(ns).shnxsub(0);
  subsize    = int(pupd/nxsub);

  mesvec = array(float,2*wfs(ns)._nsub);

  // the phase is in microns. this scaling factor restore it in radian
  // at the wfs lambda
  phasescale = float(2*pi/wfs(ns).lambda);   // wfs.lambda in microns

  if (wfs(ns).shmethod == 1) {
    
    toarcsec = float(wfs(ns).lambda/2.0/pi/(tel.diam/sim.pupildiam)/4.848);

    err = _shwfsSimple(&pupsh, &phase, phasescale, wfs(ns)._tiltsh, size, size,
                       wfs(ns)._istart, wfs(ns)._jstart, subsize, subsize,
                       wfs(ns)._nsub, toarcsec, &mesvec);

  } else {

    sdim       = long(2^ceil(log(subsize)/log(2)+1));
    sdimpow2   = int(log(sdim)/log(2));
    threshold = array(float,wfs(ns)._nsub)+wfs(ns).shthreshold;
    
    err = _shwfs(&pupsh, &phase, phasescale, wfs(ns)._tiltsh, int(size), wfs(ns)._istart,
                 wfs(ns)._jstart, subsize, subsize, wfs(ns)._nsub, sdimpow2, wfs(ns)._domask,
                 wfs(ns)._submask, wfs(ns)._kernel, wfs(ns)._kernels, wfs(ns)._kerfftr,
                 wfs(ns)._kerffti, initkernels(ns), wfs(ns)._kernelconv,
                 wfs(ns)._binindices, wfs(ns)._binxy, wfs(ns)._centroidw,
                 wfs(ns)._fimage, wfs(ns)._imistart, wfs(ns)._imjstart,
                 wfs(ns)._fimnx , wfs(ns)._fimny, wfs(ns)._fluxpersub, wfs(ns)._raylfluxpersub,
                 wfs(ns)._skyfluxpersub, &threshold, wfs(ns)._bias, wfs(ns)._flat,
                 float(wfs(ns).ron), float(wfs(ns).darkcurrent*loop.ittime),
                 int(wfs(ns).noise), int(wfs(ns).rayleighflag), wfs(ns)._rayleigh,
                 wfs(ns)._bckgrdcalib, wfs(ns)._bckgrdinit, wfs(ns)._bckgrdsub,
                 &mesvec, wfs(ns)._cyclecounter, wfs(ns).nintegcycles);

    initkernels(ns) = 0n;

    wfs(ns)._cyclecounter += 1;
    if (wfs(ns)._cyclecounter > wfs(ns).nintegcycles) {wfs(ns)._cyclecounter = 1;}
  }
  
  if (err != 0) {error,"problem in _shwfs";}

  mesvec *= wfs(ns)._centroidgain;

  // return measurement vector in arcsec (thanks to centroiw):
  return mesvec;
}

//----------------------------------------------------
func doInter(disp=,sleep=)
/* DOCUMENT doInter(disp=)
   
Measure the interaction matrix.
Each actuator are set and a measurement vector is taken.
The reference (for phase=0) is subtracted.

Keyword disp: set to display stuff as it goes.

This routine uses:
- dm._nact, _n1, _n2, _def (extern)
- mircube (extern)

This routine calls:
- multWfsIntMat

This routine sets:
- iMat (extern)
SEE ALSO: prepSVD, buildComMat
*/

{
  indexDm       = array(long,2,ndm);
  indexDm(,1)   = [1,dm(1)._nact];
  for (nm=2;nm<=ndm;nm++) {
    indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
  }

  // Loop on each mirror:
  for (nm=1;nm<=ndm;nm++) {
    n1 = dm(nm)._n1;
    n2 = dm(nm)._n2;
    if (sim.verbose==2) {write,format="Doing DM# %d, actuator %s",nm," ";}
    // Loop on each actuator:
    command = array(float,dm(nm)._nact);
    for (i=1;i<=dm(nm)._nact;i++) {
      if (sim.verbose==2) {write,format="%d ",i;}
      mircube  *= 0.0f; command *= 0.0f;
      command(i) = dm(nm).push4imat;
      mircube(n1:n2,n1:n2,nm) = compDmShape(nm,&command);
      // Fill iMat (reference vector subtracted in multWfsIntMat):
      iMat(,i+indexDm(1,nm)-1) = multWfsIntMat(disp=disp)/dm(nm).push4imat;
      // display, if requested:
      // WFS spots:
      if (is_set(disp)) {
        fma;
        plt,sim.name,0.01,0.01,tosys=0;
        if (!allof(wfs.shmethod ==1)) {
          disp2D,wfs._fimage,wfs.gspos(1,),wfs.gspos(2,),2;
          mypltitle,"WFSs spots",[0.,-0.005],height=12;
        }
        // mirror surface
        plsys,3;
        limits,square=1;
        for (j=1;j<=ndm;j++) {pli,mircube(,,j),j,0.,j+1,1.;}
        mypltitle,"DM(s)",[0.,0.008],height=12;
        if ((dm(nm).type == "aniso") && (sim.debug == 2)) hitReturn;
      }
      if (sleep) usleep,sleep;
    }
    if (sim.verbose==2) {write," ";}
  }


  // Display if needed:
  if ((sim.debug>=1) || (disp == 1)) {
    tv,-iMat,square=1;
    pltitle,"Interaction Matrix";
    if (sim.debug >= 1) typeReturn;
  }
}
//----------------------------------------------------
func prepSVD(imat,disp=,subs=,svd=)
/* DOCUMENT func prepSVD(imat,disp=)
   Does the SVD decomposition and fill out modToAct
   and mesToMod matrices for further use by buildComMat()

   Keyword disp: set if display is required.

   This routine uses:
   - imat (input)

   This routine calls:
   - SVdec

   This routine sets:
   - eigenvalues (extern)
   - modToAct (extern)
   - mesToMod (extern)

   SEE ALSO: doInter, buildComMat
*/
{
  // Define some extern variables:
  extern modToAct,mesToMod,eigenvalues;

  // Decompose to prepare inversion:
  if (sim.verbose==2) {write,"Doing SVD\n";}
  eigenvalues 	= SVdec(imat,u,vt);

  // Some debug output if needed:
  if (sim.verbose>=2) {
    write,"Normalized eigenvalues:";
    write,eigenvalues/max(eigenvalues);
    do {
      plot,eigenvalues/max(eigenvalues);
      plg,array(1/(*mat.condition)(subs),numberof(eigenvalues)),color="red",type=2;
      limits,square=0;
      if (is_set(svd)) {
        ths = "";
        read,prompt=swrite(format="Threshold [%f] (return to continue): ",\
                           1/(*mat.condition)(subs)),ths;
        if ( strlen(ths) == 0 ) {
          th =  1/(*mat.condition)(subs);
          change = 0;
        } else {
          th=0.; sread,ths,th;
          (*mat.condition)(subs) = 1./th;
          change = 1;
        }
      }
    } while (change == 1);
    typeReturn;
  } else if (sim.verbose == 1) {
    write,format="Smallests 2 normalized eigenvalues = %f",
      eigenvalues(-1)/max(eigenvalues),eigenvalues(0)/max(eigenvalues);
  }
  
  // Mode-to-Actuator and Actuator-to_mode matrices:
  // to be used as follow:
  // modes-coef    = actToMod(,+) * command-coef(+)
  // actuator-coef = modToAct(,+) * modes-coef(+)
  modToAct    = transpose(vt);
  //  actToMod    = LUsolve(modToAct);
  mesToMod    = transpose(u);  // used to be called ut

  // Some debug display if needed:
  if (sim.debug==2) {
    tv,modToAct; pltitle,"modToAct Matrix";
    typeReturn;
    tv,modToAct(,+)*modToAct(,+);
    pltitle,"modToAct(,+)*modToAct(,+)";
    typeReturn;

    tv,mesToMod; pltitle,"mesToMod Matrix";
    typeReturn;
    tv,mesToMod(,+)*mesToMod(,+);
    pltitle,"mesToMod(,+)*mesToMod(,+)";
    typeReturn;
  }
}

//----------------------------------------------------

func buildComMat(condition,modalgain,subsystem=,all=,nomodalgain=,disp=)
/* DOCUMENT buildComMat(all=,nomodalgain=,disp=)
   Build the command matrix from V,UT, the eigenvalues and the modal gains
   F.Rigaut, June 17,2002
   Input keywords:
   all: if not set, one (1) mode is forced to be discarded (useful if regular
   AO system using a DM with a piston component and condition number too
   large). Normally, set all=1.
   nomodalgain: if set, the modal gain are not taken into account.
   disp: set to display stuff.

   This routine uses:
   - dm._def, _nact, _n1, _n2 (extern)
   - ipupil (extern)
   - mat.condition (extern)
   - modalgain (extern)
   - eigenvalues (extern)
   - modToAct (extern)
   - mesToMod (extern)
     
   This routine calls:
   - Nothing
     
   This routine sets:
   - NModesControlled (extern)

   This routine returns:
   - cMat

   SEE ALSO: doInter, prepSVD
*/
{
  extern NModesControlled;  

  neigen = numberof(eigenvalues);
  
  mev 	= array(float,neigen,neigen);

  mask = ((eigenvalues/max(eigenvalues)) > (1./condition));
  
  if (is_set(nomodalgain)) {
    ev = eigenvalues;
  } else {
    ev = eigenvalues/modalgain;
  }
  
  // Including the mode gains and eigenvalues here:
  for (i=1;i<=neigen;i++) {
    if (mask(i) == 1) {mev(i,i)=1./ev(i);}
  }

  // the last eigenvalue is filtered except if all is set.
  if (!is_set(all)) {mev(0,0) = 0.;}

  NModesControlled = sum(mev != 0.);

  // Compute the Command matrix:
  cmat = (modToAct(,+)*mev(+,))(,+) * mesToMod(+,);

  if (sim.verbose>=1) {
    write,long(clip(sum(mask == 0),1-long(is_set(all)),)),
      format="%i modes discarded in the inversion\n";
  }

  //=========================================
  // DISPLAY OF THE FILTERED MODES:
  // BEGINNING OF DISPLAY, NOTHING IMPORTANT
  // UNTIL "END OF DISPLAY"
  //=========================================
  
  if ((sim.debug>=1) && (disp >= 1) ) {

    n1 = (sim._size-sim.pupildiam)/2-2;
    n2 = (sim._size+sim.pupildiam)/2+2;
    nxy = n2-n1+1;
    sswfs = where(wfs.subsystem == subsystem);
    cubphase = array(float,[3,nxy,nxy,numberof(sswfs)]);
    subpos = wfs(sswfs).gspos;
    disp2D,cubphase,subpos(1,),subpos(2,),1,zoom=0.9,init=1;
    
    // find out to where is DM N in this subsystem "modToAct" matrix:
    ssdm = where(dm.subsystem == subsystem);
    //    indexDm = array(long,[2,2,numberof(ssdm)]);
    indexDm = array(long,[2,2,ndm]);
    index   = 0;
    for (nm=1;nm<=ndm;nm++) {
      if (dm(nm).subsystem == subsystem) {
        indexDm(,nm) = [index+1,index+dm(nm)._nact];
        index = index+dm(nm)._nact;
      }
    }
    
    write,"Displaying filtered modes, type 'q' to exit display";
    // loop on modes
    mircube *= 0.0f;
    for (i=neigen;i>=NModesControlled+1;i--) {
      // loop over dm to build mircube
      for (nm=1; nm<=ndm; nm++) {
        if (dm(nm).subsystem == subsystem) {
          n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = n2-n1+1;
          scommand = float(modToAct(indexDm(1,nm):indexDm(2,nm),i));

          mircube(n1:n2,n1:n2,nm) = compDmShape(nm,&scommand);
        }
      }

      // Loop over WFS to get integrated phase
      n1 = (sim._size-sim.pupildiam)/2-2;
      n2 = (sim._size+sim.pupildiam)/2+2;
      nxy = n2-n1+1;
      for (ns=1;ns<=nwfs;ns++) {
        if (wfs(ns).subsystem == subsystem) {
          phase = getPhase2dFromDms(ns,"wfs")*ipupil;
          // fill cubphase
          cubphase(,,ns) = phase(n1:n2,n1:n2);
        }
      }

      // display using disp2D of integrated phases
      fma;
      plt,sim.name,0.01,0.01,tosys=0;
      disp2D,cubphase,subpos(1,),subpos(2,),1;
      pltitle,swrite(format="Mode %d, Normalized eigenvalue = %f",
                     i,eigenvalues(i)/max(eigenvalues));

      // display of DM shapes
      plsys,3;
      limits,square=1;
      tmp = [];
      for (nm=1;nm<=ndm;nm++) {
        if (dm(nm).subsystem == subsystem) {
          if (dm(nm).alt == 0.) {mircube(,,nm) *= ipupil;}
          grow,tmp,transpose(mircube(,,nm));
        }
      }
      if (tmp != []) {pli,transpose(tmp);}
      mypltitle,"DM(s)",[0.,0.008],height=12;

      rep = typeReturn();
      if (rep == "q") {break;}
    }
  }
  //=========================================
  //             END OF DISPLAY
  //=========================================


  return cmat;
}

//----------------------------------------------------
func swapScreens
/* DOCUMENT func swapScreens
   Swap the phase screens. This is to get better statistics
   out of fewer phase screens and iterations.
   The 2nd phase screen becomes the 1rst one, 3->2, etc...
   This routine uses the phase screens and normalization
   factor stored in extern by getTurbPhaseInit
   SEE ALSO:
*/
{
  weight = currentScreenNorm;
  // avoid division per zero
  // it's legitimate to have one screen == 0
  weight += (weight == 0);
  // get rid of layer strength normalization factor:
  pscreens = pscreens/weight(-,-,);

  // swap screens:
  pscreens = pscreens(,,((indgen(nscreens)) % nscreens)+1);

  // re-apply normalization:
  weight = currentScreenNorm;
  pscreens = pscreens*weight(-,-,);   
}
//----------------------------------------------------

func getTurbPhaseInit(skipReadPhaseScreens=)
/* DOCUMENT getTurbPhaseInit(void)
   Belongs to the yao package suite.
   Initializes everything for getTurbPhase (see below), which
   returns the interpolated, integrated phase to the loop function
   for iteration number "iter". Returns the "ok" parameter which is
   set to 0 if this step is not to be used to compute statistics.
   AUTHOR: F.Rigaut, June 11, 2002.
   SEE ALSO: aoinit, aoloop, getTurbPhase.
*/
{
  extern pscreens,// contains screens, dim [dim Screen X,dim Screen Y, nscreens]
    nscreens,     // number of screens
    noptics,      // number of optics
    optphasemaps, // phase maps for each optics
    xposvec,      // contains X position vs iteration# of center of beam in
                  // phase screen N, in pixel coords. dim [N iteracurwtions,nscreens]
    yposvec,      // same for Y. This should be cst in the current implementation
                  // as the screens have been cut in Y, therefore are not periodic
                  // and therefore we can not wrap.
    wfsxposcub,   // contains the position of the ray at which a given X pixel
                  // intersects the Nth phase screen, for each WFS GS.
                  // dimension [dimX outphase, nscreens, #GS] 
    wfsyposcub,   // Same for Y
    gsxposcub,    // contains the position of the ray at which a given X pixel
                  // intersects the Nth phase screen, for each star at which the perf
                  // is evaluated dimension [dimX outphase, nscreens, #Star] 
    gsyposcub,    // Same for Y
    dmwfsxposcub, // contains the position of the ray at which a given X pixel
                  // intersects the Nth mirror, for each WFS GS.
                  // dimension [dimX outphase, ndm, #GS] 
    dmwfsyposcub, // Same for Y
    dmgsxposcub,  // contains the position of the ray at which a given X pixel
                  // intersects the Nth mirror, for each star at which the perf
                  // is evaluated dimension [dimX outphase, ndm, #Star] 
    dmgsyposcub,  // Same for Y

    optwfsxposcub,// contains the position of the ray at which a given X pixel
                  // intersects the Nth optics, for each WFS GS.
                  // dimension [dimX outphase, noptics, #GS] 
    optwfsyposcub,// Same for Y
    optgsxposcub, // contains the position of the ray at which a given X pixel
                  // intersects the Nth optics, for each star at which the perf
                  // is evaluated dimension [dimX outphase, noptics, #Star] 
    optgsyposcub, // Same for Y

    xmargins,     // value of xposvec so that no pixel indice < 1 or > dimx(screen)
                  // in any off-axis beam and altitude (low margin, up margin)
    ymargins,     // same for Y. We use that in getTurbPhase to determine when to
                  // wrap.
    statsokvec,   // 1 if it is ok to collect stats at this iteration.
                  // 0 if not, e.g. we just did a jump. dim [iteration]
    inithistory,  // 1 if init has been done.
    screendim;    // [phase screen X dim, Y dim] before they are extended for safe wrapping

  extern currentScreenNorm; // current screen normalization. Used when swaping screen.

  // Define a few variables:

  nscreens = numberof(*atm.screen);
  if (opt!=[]) noptics = numberof(opt.phasemaps); else noptics=0;

  wfsxposcub = wfsyposcub = array(float,[3,_n,nscreens,nwfs]);
  gsxposcub = gsyposcub = array(float,[3,_n,nscreens,target._ntarget]);

  dmwfsxposcub = dmwfsyposcub = array(float,[3,_n,ndm,nwfs]);
  dmgsxposcub = dmgsyposcub = array(float,[3,_n,ndm,target._ntarget]);

  if (noptics) {
    optwfsxposcub = optwfsyposcub = array(float,[3,_n,noptics,nwfs]);
    optgsxposcub = optgsyposcub = array(float,[3,_n,noptics,target._ntarget]);
  }

  //=======================================================
  // READS AND NORMALIZE THE PHASE SCREENS
  // the phase screens now (v2.4) are normalized in microns 
  //=======================================================

  if (!is_set(skipReadPhaseScreens)) {
    // Compute normalization factor:
    (*atm.layerfrac) = (*atm.layerfrac)/sum(*atm.layerfrac);
    weight = float(sqrt(*atm.layerfrac)*(atm.dr0at05mic/
                                         cos(gs.zenithangle*dtor)^0.6/sim.pupildiam)^(5./6.));
    // above: in radian at 0.5 microns
    weight = weight * float(0.5/(2*pi));
    // ... and now in microns.


    /*=====================================================
      How to relate r0(layer) and atm.layerfrac ?
      r0(i)   = r0 of layer i
      r0tot   = total r0
      f(i)    = "fraction" in layer i ( = (*atm.layerfrac)(i) )
      we have:
      weight(i) = sqrt(f(i)) * (D/r0tot)^(5/6.) = (D/r0(i))^(5/6.)
      thus
      f(i) = (r0tot/r0(i))^(5./3)

      inversely, we have:
      r0(i) = r0tot / f(i)^(3./5)
      =====================================================*/
      
    if (sim.verbose>=1) {
      write,format="Reading phase screen \"%s\"\n",(*atm.screen)(1);
    }

    // read the first one, determine dimensions, stuff it:
    tmp      = fitsRead((*atm.screen)(1));
    dimx     = dimsof(tmp)(2);
    dimy     = dimsof(tmp)(3);
    screendim = [dimx,dimy];
    // Extend dimension in X for wrapping issues
    pscreens = array(float,[3,dimx+2*sim._size,dimy,nscreens]);
    // Stuff it
    pscreens(1:dimx,,1) = tmp;
    // free RAM
    tmp      = [];

    // Now read all the other screens and put in pscreens
    for (i=2;i<=nscreens;i++) {
      if (sim.verbose>=1) {
        write,format="Reading phase screen \"%s\"\n",(*atm.screen)(i);
      }
      pscreens(1:dimx,,i) = fitsRead((*atm.screen)(i));
    }

    // Extend the phase screen length for safe wrapping:
    pscreens(dimx+1:,,) = pscreens(1:2*sim._size,,);
    // Can't do in Y as the phase screens are not periodic (they have been cutted)
    //  pscreens(,dimy+1:,) = pscreens(,1:sim._size,);
    
    // apply weights to each phase screens (normalize):
    // the screens are expressed in microns
    pscreens = pscreens*weight(-,-,);
    currentScreenNorm = weight;

    //=============================================
    // READ THE OPTICS PHASE MAPS
    // the dimension of all the optical phase maps
    // should be the same. These phase maps should
    // be big enough so that the interpolation will
    // not go outside of the provided map. That is,
    // the indices in opt**xposcub (and y), with **
    // being wfs and gs, should not define points
    // outside of the provided phase maps.
    // THIS IS LEFT TO THE RESPONSIBILITY OF THE
    // USER. I do however a simple check below.
    // Normaly, this should not be a problem:
    // if some indices are outside the maps, it means
    // either that the maps are not big enough, in
    // which case we should not use it, or that
    // there was an error in the definition of the
    // altitude.
    //=============================================
    if (noptics) {
      if (sim.verbose) {
        psize  = tel.diam/sim.pupildiam;
        write,"Reading optics. I am expecting phase maps with a pixel";
        write,format=" size of %fm/pixel, as projected in the entrance\n",psize;
        write,"pupil plane. It is also your responsibility to provide phase";
        write,"maps of adequate dimension, i.e. large enough";
      }
      if (sim.verbose) {
        write,format="Reading phase map for optics \"%s\"\n",opt(1).phasemaps;
      }
      tmp = fitsRead(opt(1).phasemaps);
      optdims      = dimsof(tmp);
      optdimx      = optdims(2);
      optdimy      = optdims(3);

      if (optdimx != optdimy)
        error,"Optics phase maps should be square arrays";

      optphasemaps = array(float,[3,optdimx,optdimy,noptics]);
      // Stuff it
      optphasemaps(,,1) = float(tmp);
  
      for (i=2;i<=noptics;i++) {
        if (sim.verbose>=1) {
          write,format="Reading phase map for optics \"%s\"\n",opt(i).phasemaps;
        }
        tmp = fitsRead(opt(i).phasemaps);
        if (anyof(dimsof(tmp) != optdims))
          error,"All optics phase maps should have the same dimensions";
        optphasemaps(,,i) = float(tmp);
      }
      opt._cent = optdimx/2.+(sim._cent-sim._size/2.);
      // (sim._cent-sim._size/2.) =  0. or 0.5
    }
  }
  
  //====================================
  // PREPARE GEOMETRY FOR INTERPOLATION
  //====================================

  // Build a vector of position (integer position in equivalent
  // iterations:
  iposvec = indgen(loop.niter);
  iposvec = iposvec+((iposvec-1)/loop.skipevery)*loop.skipby;

  // build a vector of the iterations at which statistics
  // should be accumulated. 1 if ok to accumulate statistics, 0 if not.
  statsokvec   = (indgen(loop.niter)-1) % loop.skipevery;
  statsokvec   = (statsokvec >= loop.startskip);

  // Build the position vector vs iteration by phase screen
  deltax  = sim.pupildiam/tel.diam*(*atm.layerspeed)*loop.ittime;
  // has to start away from edge as this is the coordinate of the beam center
  /// will modified later on when we know the beams geometry. see few lines below.
    xposvec = (1.+iposvec*deltax(-,)*cos(dtor*(*atm.winddir)(-,)));
    yposvec = (1.+iposvec*deltax(-,)*sin(dtor*(*atm.winddir)(-,)));

    psize  = tel.diam/sim.pupildiam;  // pixel in meter

    //===========================================================
    // PRE-COMPUTATION OF THE INTERSECT POSITIONS FOR EACH WFS GS
    //===========================================================

    // zero-centered vector of position (our reference)
    xref = indgen(_n)-(_n+1)/2.;
    yref = indgen(_n)-(_n+1)/2.;

    // loop on WFS
    for (n=1;n<=nwfs;n++) {
    
      // loop on screen
      for (ns=1;ns<=nscreens;ns++) {
      
        // offsets of the center of beam on screen NS
        xoff = (wfs(n).gspos)(1)*4.848e-6*(*atm._layeralt)(ns)/psize;
        yoff = (wfs(n).gspos)(2)*4.848e-6*(*atm._layeralt)(ns)/psize;
      
        // if we are dealing with LGS, there is a geometric
        // factor on the beam size
        if (wfs(n)._gsalt != 0) {
          fact = (wfs(n)._gsalt-(*atm._layeralt)(ns))/wfs(n)._gsalt;
        } else {
          fact = 1.;
        }
      
        // Compute and stuff the array
        wfsxposcub(,ns,n) = xref*fact + xoff;
        wfsyposcub(,ns,n) = yref*fact + yoff;
      }
    
      // loop on DMs
      for (ns=1;ns<=ndm;ns++) {
      
        // offsets of the center of beam on DM NS
        xoff = (wfs(n).gspos)(1)*4.848e-6*(dm.alt)(ns)/psize;
        yoff = (wfs(n).gspos)(2)*4.848e-6*(dm.alt)(ns)/psize;
      
        // if we are dealing with LGS, there is a geometric
        // factor on the beam size
        if (wfs(n)._gsalt != 0) {
          fact = (wfs(n)._gsalt-(dm.alt)(ns))/wfs(n)._gsalt;
        } else {
          fact = 1.;
        }
      
        // Compute and stuff the array
        dmwfsxposcub(,ns,n) = xref*fact + xoff;
        dmwfsyposcub(,ns,n) = yref*fact + yoff;
      }

      // loop on optics
      for (ns=1;ns<=noptics;ns++) {
      
        // offsets of the center of beam on optics #NS
        xoff = (wfs(n).gspos)(1)*4.848e-6*opt(ns).alt/psize;
        yoff = (wfs(n).gspos)(2)*4.848e-6*opt(ns).alt/psize;
      
        // if we are dealing with LGS, there is a geometric
        // factor on the beam size
        if (wfs(n)._gsalt != 0) {
          fact = (wfs(n)._gsalt-opt(ns).alt)/wfs(n)._gsalt;
        } else {
          fact = 1.;
        }
      
        // Compute and stuff the array
        optwfsxposcub(,ns,n) = xref*fact + xoff;
        optwfsyposcub(,ns,n) = yref*fact + yoff;
      }

    }

    // type conversion:
    wfsxposcub = float(wfsxposcub);
    wfsyposcub = float(wfsyposcub);
    dmwfsxposcub = float(dmwfsxposcub);
    dmwfsyposcub = float(dmwfsyposcub);
    optwfsxposcub = float(optwfsxposcub);
    optwfsyposcub = float(optwfsyposcub);
  
  
    //=============================================================
    // PRE-COMPUTATION OF THE INTERSECT POSITION FOR EACH PERF STAR
    //=============================================================

    // loop on target
    for (n=1;n<=target._ntarget;n++) {
    
      // loop on screen
      for (ns=1;ns<=nscreens;ns++) {
      
        // offsets of the center of beam on screen NS
        xoff = (*target.xposition)(n)*4.848e-6*(*atm._layeralt)(ns)/psize;
        yoff = (*target.yposition)(n)*4.848e-6*(*atm._layeralt)(ns)/psize;

        // note: that's target. we can't be dealing with LGS here (no fact)
        // Compute and stuff the array
        gsxposcub(,ns,n) = xref + xoff;
        gsyposcub(,ns,n) = yref + yoff;
      }

      // loop on DMs
      for (ns=1;ns<=ndm;ns++) {
      
        // offsets of the center of beam on DM NS
        xoff = (*target.xposition)(n)*4.848e-6*(dm.alt)(ns)/psize;
        yoff = (*target.yposition)(n)*4.848e-6*(dm.alt)(ns)/psize;

        // note: that's target. we can't be dealing with LGS here (no fact)
        // Compute and stuff the array
        dmgsxposcub(,ns,n) = xref + xoff;
        dmgsyposcub(,ns,n) = yref + yoff;
      }
    
      // loop on optics
      for (ns=1;ns<=noptics;ns++) {
      
        // offsets of the center of beam on DM NS
        xoff = (*target.xposition)(n)*4.848e-6*opt(ns).alt/psize;
        yoff = (*target.yposition)(n)*4.848e-6*opt(ns).alt/psize;

        // note: that's target. we can't be dealing with LGS here (no fact)
        // Compute and stuff the array
        optgsxposcub(,ns,n) = xref + xoff;
        optgsyposcub(,ns,n) = yref + yoff;
      }
    
    }
    // type conversion:
    gsxposcub = float(gsxposcub);
    gsyposcub = float(gsyposcub);
    dmgsxposcub = float(dmgsxposcub);
    dmgsyposcub = float(dmgsyposcub);
    optgsxposcub = float(optgsxposcub);
    optgsyposcub = float(optgsyposcub);


    //======================================
    // SOME CHECKS TO AVOID INDICES OVERFLOW
    //======================================

    // Now we can modify xposvec and yposvec to make sure we are not going
    // out of the phase screen arrays
    xmargins = abs([min(_(wfsxposcub(*),gsxposcub(*))),max(_(wfsxposcub(*),gsxposcub(*)))]);
    ymargins = abs([min(_(wfsyposcub(*),gsyposcub(*))),max(_(wfsyposcub(*),gsyposcub(*)))]);
    xposvec = xposvec - min(xposvec) + xmargins(1) +1;
    yposvec = yposvec - min(yposvec) + ymargins(1) +1;

    // wrap so that it never goes out of bound (including off axis stuff)
    // have to do that for each screens as they are moving at different speeds
    for (ns=1;ns<=nscreens;ns++) {
      // we know that xposvec is now purely positive (lines above),
      // so that it is enough to insure that it will never go above upper limit:
      // If pixel position of beam center (xposvec) is larger than the 
      // initial dimension of the screen + the xmargin necessary so that
      // off-axis beams don't hit negative indices, then it is
      // time to wrap by subtracting the original screen Xdim to
      // xposvec.
      xposvec(,ns) = xmargins(1)+ ( (xposvec(,ns)-xmargins(1)) % screendim(1));
      // don't do for Y as we're not allowed to move along Y
      // (screen not periodic)
    }
  
    // type conversion:
    xposvec = float(xposvec);
    yposvec = float(yposvec);

    // so now we have everything initiliazed, and we will just have to
    // interpolate the phase screen at points
    // xposvec(iteration,screen#) + wfsxposcub(,screen#,wfs#)
    // and corresponding for Y, and integrate on screen#
    // to get the phase for a given WFS
    // this integration is done by the C routine _get2dPhase
    // which take, in addition to the screens and output phase parameters,
    // only a set of X and Y positions as input [the one we just talked
    // about, xposvec(iteration) + wfsxposcub(,,wfs#) ].

    // check that y index does not overflow:
    getTurbPhaseInitCheckOverflow;
  
    inithistory = 1;
    return 1;
}
//----------------------------------------------------

func getTurbPhase(iter,nn,type)

/* DOCUMENT getTurbPhase(iter,n,type)
   Belongs to the yao package suite.
   Returns the interpolated, integrated phase to the loop function
   for iteration number "iter". 
   You have to call getTurbPhaseInit to initialize prior using.
   AUTHOR: F.Rigaut, June 11, 2002.
   Modified 2003 feb 24 to v1.2.
   now takes additional parameters offsets (offsets from field center
   in arcsec) and altgs (altitude of GS in meters).
   Modified 2003 dec 9 to implement call to C routine _get2dPhase
   SEE ALSO: aoinit, aoloop, getTurbPhaseInit.
*/
{
  if (!inithistory) {error,"getTurbPhase has not been initialized !";}
  
  sphase = array(float,_n,_n);
  bphase = array(float,sim._size,sim._size);

  // Now we can call the C interpolation routine and get the integrated
  // phase for this star
  // there are a few things to do to get ready
  psnx = dimsof(pscreens)(2);
  psny = dimsof(pscreens)(3);
  nscreens = dimsof(pscreens)(4);

  // here we have a branch to be able to process wfs and targets with the same
  // subroutine, this one.
  if (type == "wfs") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = wfsxposcub(,,nn)+xposvec(iter,)(-,);
    yshifts = wfsyposcub(,,nn)+yposvec(iter,)(-,);
  } else if ( type == "target") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = gsxposcub(,,nn)+xposvec(iter,)(-,);
    yshifts = gsyposcub(,,nn)+yposvec(iter,)(-,);
  }

  ishifts = int(xshifts);  xshifts = xshifts - ishifts;
  jshifts = int(yshifts);  yshifts = yshifts - jshifts;
 
  err = _get2dPhase(&pscreens,psnx,psny,nscreens,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  if (err != 0) {error,"Error in getTurbPhase";}
  
  bphase(_n1:_n2,_n1:_n2) = sphase;
  
  return bphase;
}

//----------------------------------------------------
func getPhase2dFromDms(nn,type)
/* DOCUMENT func getPhase2dFromDms(mircube, mir_sh)
   adapted from the IDL function of the same name in tomoclose.pro
   mircube(extern)	= cube or image of phase screen (mirrors)
   first plan = first mirror shape
   second plan = second mirror shape
   etc...
   should be of size _n*_n
   mir_sh     = 2d vector of shift in number of pixels per screen
   *_sh([x,y],screen#)
   SEE ALSO: 
*/
{
  
  sphase = array(float,_n,_n);
  bphase = array(float,sim._size,sim._size);
  
  // Now we can call the C interpolation routine and get the integrated
  // phase for this star
  // there are a few things to do to get ready
  psnx = dimsof(mircube)(2);
  psny = dimsof(mircube)(3);
  nmirrors = dimsof(mircube)(4);

  // here we have a branch to be able to process wfs and targets with the same
  // subroutine, this one.
  if (type == "wfs") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = dmwfsxposcub(,,nn)+(sim._cent+dm.misreg(1,)-1)(-,);
    yshifts = dmwfsyposcub(,,nn)+(sim._cent+dm.misreg(2,)-1)(-,);
  } else if ( type == "target") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = dmgsxposcub(,,nn)+(sim._cent+dm.misreg(1,)-1)(-,);
    yshifts = dmgsyposcub(,,nn)+(sim._cent+dm.misreg(2,)-1)(-,);
  }

  ishifts = int(xshifts); xshifts = xshifts - ishifts;
  jshifts = int(yshifts); yshifts = yshifts - jshifts;

  err = _get2dPhase(&mircube,psnx,psny,nmirrors,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  if (err != 0) {error,"Error in getPhase2dFromDms";}
  
  bphase(_n1:_n2,_n1:_n2) = sphase;
  
  return bphase;
}

//----------------------------------------------------
func getPhase2dFromOptics(nn,type)
/* DOCUMENT func getPhase2dFromOptics(nn,type)
   adapted from the IDL function of the same name in tomoclose.pro
   nn = wfs or GS #
   type = "wfs" or "target"
   SEE ALSO: 
*/
{
  if (opt==[]) return 0.0f;
  
  sphase = array(float,_n,_n);
  bphase = array(float,sim._size,sim._size);
  
  // Now we can call the C interpolation routine and get the integrated
  // phase for this star
  // there are a few things to do to get ready
  psnx  = dimsof(optphasemaps)(2);
  psny  = dimsof(optphasemaps)(3);
  nopts = dimsof(optphasemaps)(4);

  // here we have a branch to be able to process wfs and targets with the same
  // subroutine, this one.
  if (type == "wfs") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = optwfsxposcub(,,nn)+(opt._cent+opt.misreg(1,)-1)(-,);
    yshifts = optwfsyposcub(,,nn)+(opt._cent+opt.misreg(2,)-1)(-,);
  } else if ( type == "target") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = optgsxposcub(,,nn)+(opt._cent+opt.misreg(1,)-1)(-,);
    yshifts = optgsyposcub(,,nn)+(opt._cent+opt.misreg(2,)-1)(-,);
  }

  ishifts = int(xshifts); xshifts = float(xshifts - ishifts);
  jshifts = int(yshifts); yshifts = float(yshifts - jshifts);

  err = _get2dPhase(&optphasemaps,psnx,psny,nopts,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  //  if (err != 0) {error,"Error in getPhase2dFromOptics";}
  
  bphase(_n1:_n2,_n1:_n2) = sphase;
  
  return bphase;
}

//----------------------------------------------------
func correctUpLinkTT(phase,ns)
/* DOCUMENT func correctUpLinkTT(phase)
     
SEE ALSO:
*/
{
  wfs(ns)._upttcommand += wfs(ns).uplinkgain * wfs(ns)._tt;

  phase -= wfs(ns)._upttcommand(1) * tip1arcsec;
  phase -= wfs(ns)._upttcommand(2) * tilt1arcsec;

  return phase;
}
//----------------------------------------------------
func splitWfsVector(v)
/* DOCUMENT func splitWfsVector(v)
   splits the single vector (out of multWfs or multWfsIntMat)
   into as many individual wfs vector as there are sensors.
   Return a pointer vector to the individual wfs vectors.
   SEE ALSO: splitDMCommandVector(v)
*/
{
  vp = [];
  vs = v(1:wfs(1)._nmes);
  grow,vp,&vs;
  iend = wfs(1)._nmes;

  for (ns=2;ns<=nwfs;ns++) {
    istart = iend+1;
    iend = istart+wfs(ns)._nmes-1;
    vs = v(istart:iend);
    grow,vp,&vs;
  }

  return vp;
}
//----------------------------------------------------
func splitDMCommandVector(v)
/* DOCUMENT func splitDMCommandVector(v)
   splits the single vector (out of cMat matrix multiply in aoloop)
   into as many individual command vector as there are DMs.
   Return a pointer vector to the individual command vectors.
   SEE ALSO: splitWfsVector
*/
{
  vp = [];
  vs = v(1:dm(1)._nact);
  grow,vp,&vs;
  iend = dm(1)._nact;

  for (ns=2;ns<=ndm;ns++) {
    istart = iend+1;
    iend = istart+dm(ns)._nact-1;
    vs = v(istart:iend);
    grow,vp,&vs;
  }

  return vp;
}
//----------------------------------------------------
func multWfsIntMat(disp=)
/* intmat : set if computing intMat
   as multWfs but special for IntMat acquisition
   for speed in aoloop
*/

{
  mes = [];
  for (ns=1;ns<=nwfs;ns++) {

    // Impose noise = rmsbias = rmsflat = 0 for interaction matrix measurements
    noiseOrig = wfs(ns).noise; wfs(ns).noise = 0n;
    cycleOrig = wfs(ns).nintegcycles; wfs(ns).nintegcycles = 1;
    if (wfs(ns).type == "hartmann" ) {
      kconvOrig = wfs(ns)._kernelconv; wfs(ns)._kernelconv = 1n;
      bias  = *wfs(ns)._bias; *wfs(ns)._bias = *wfs(ns)._bias*0.0f;
      flat  = *wfs(ns)._flat; *wfs(ns)._flat = *wfs(ns)._flat*0.0f+1.0f;
    }

    offsets = wfs(ns).gspos;
    phase   = getPhase2dFromDms(ns,"wfs");
    // uncomment if needed:
    //    phase  += getPhase2dFromOptics(ns,"wfs");

    if (wfs(ns).type == "curvature") {smes = CurvWfs(pupil,phase,ns);}
    if (wfs(ns).type == "hartmann" ) {smes = ShWfs(ipupil,phase,ns);}
    if (wfs(ns).type == "pyramid")   {smes = PyramidWfs(pupil,phase);}
    
    // subtract the reference vector for this sensor:
    smes = smes - *wfs(ns)._refmes;

    // compute the TT and subtract if required:
    if (wfs(ns).filtertilt) {
      wfs(ns)._tt(1) = sum( smes * (*wfs(ns)._tiprefvn) );
      wfs(ns)._tt(2) = sum( smes * (*wfs(ns)._tiltrefvn) );
      smes = smes - wfs(ns)._tt(1) * (*wfs(ns)._tiprefv) \
        - wfs(ns)._tt(2) * (*wfs(ns)._tiltrefv);
    }

    grow,mes,smes;
    
    // restore whatever value was in bias and flat
    if (wfs(ns).type == "hartmann" ) {
      wfs(ns)._bias = &bias; wfs(ns)._flat = &flat;
      wfs(ns)._kernelconv = kconvOrig;
    }
    wfs(ns).noise = noiseOrig;
    wfs(ns).nintegcycles = cycleOrig;

  }
  return mes;
}
//----------------------------------------------------
func multWfs(iter,disp=)
/* 
 */
{
  mes = [];
  for (ns=1;ns<=nwfs;ns++) {

    offsets = wfs(ns).gspos;
    phase   = getPhase2dFromDms(ns,"wfs");
    phase  += getPhase2dFromOptics(ns,"wfs");
    phase  += getTurbPhase(iter,ns,"wfs");

    if (wfs(ns).correctUpTT) {
      phase = correctUpLinkTT(phase,ns);
    }

    // get the measurements:
    if (wfs(ns).type == "curvature") {smes = CurvWfs(pupil,phase,ns);}
    if (wfs(ns).type == "hartmann" ) {smes = ShWfs(ipupil,phase,ns);}
    if (wfs(ns).type == "pyramid")   {smes = PyramidWfs(pupil,phase);}

    // subtract the reference vector for this sensor:
    if (wfs(ns)._cyclecounter == 1) {
      smes = smes - *wfs(ns)._refmes;
    }
    
    // compute the TT and subtract if required:
    wfs(ns)._tt(1) = sum( smes * (*wfs(ns)._tiprefvn) );
    wfs(ns)._tt(2) = sum( smes * (*wfs(ns)._tiltrefvn) );
    if (wfs(ns).filtertilt) {
      smes = smes - wfs(ns)._tt(1) * (*wfs(ns)._tiprefv) \
        - wfs(ns)._tt(2) * (*wfs(ns)._tiltrefv);
    }
    if (wfs(ns)._cyclecounter == 1) {
      wfs(ns)._lastvalidtt = wfs(ns)._tt;
    }

    grow,mes,smes;
  }
  return mes;
}
//----------------------------------------------------
func aoall(parfile,disp=,dpi=,clean=,controlscreen=)
{
  aoread,parfile;
  aoinit,disp=disp,dpi=dpi,clean=clean;
  aoloop,disp=disp,dpi=dpi,controlscreen=controlscreen;
}
//----------------------------------------------------

func aoread(parfile)
/* DOCUMENT func aoread(parfile)
   Define the relevant structure/variable, and reads out the
   AO simulation parameter file (e.g. "sh6.par"),
   Does a check of the WFS pixel/subaperture sizes.
   This is the first function to call in the ao serie (normally,
   aoinit and aoloop follow).
   This routine was kept separate from aoinit to keep the possibility
   to change variables inbetween the aoread and the aoinit calls.
   example:
   > aoread,"sh6.par"
   > loop.niter = 1000
   > aoinit

   This routine uses:
   - the parameter file (e.g. "sh6.par")
     
   This routine calls:
   - wfsCheckPixelSize
     
   This routine sets:
   - the atm,sim,wfs,dm,mat,tel,target,gs,loop structure contents
   - parprefix

   SEE ALSO: aoinit, aoloop
*/
{
  extern atm,opt,sim,wfs,dm,mat,tel,target,gs,loop,parprefix,oparfile;

  write,format="Yao, Version %s, %s\n",aoSimulVersion, aoSimulVersionDate;

  if (Y_VERSION == "1.5.12") {
    error,"SVD is broken in yorick 1.5.12 ! Get a later version.";
  }
  
  if (strmatch(parfile,".par")) {
    tmp = strpart(parfile,strword(parfile,"/",50));
    tmp = tmp(where(tmp));
    tmp = tmp(0);
    tmp = strpart(tmp,strword(tmp,".",5));
    tmp = tmp(where(tmp));
    parprefix = tmp(:-1)(sum);
  } else {
    parprefix = parfile;
  }
  oparfile = parfile;

  // flush any prior assignments to ao members
  atm=opt=sim=wfs=dm=mat=tel=target=gs=loop=cwfs=[];  

  // INIT STRUCTURES:
  atm  = atm_struct();
  opts = opt_struct();
  sim  = sim_struct();
  wfss = wfs_struct();
  dms  = dm_struct();
  mat  = mat_struct();
  tel  = tel_struct();
  target = target_struct();
  gs   = gs_struct();
  loop = loop_struct();
  paramfile = parfile;
  
  if (!fileExist(parfile)) {
    exit,swrite(format="Can not find parameter file %s !",parfile);}

  // read out the parfile. This stuffs values into the structures:
  require,parfile;

  //=====================================
  // PARAMETER CHECKS. SETS SOME DEFAULTS
  //=====================================
  checkParameters;

}

//----------------------------------------------------

func aoinit(disp=,clean=,forcemat=,svd=,dpi=,keepdmconfig=)
/* DOCUMENT func aoinit(disp=,forcemat=)
   Second function of the ao serie.
   Initialize everything in preparation for the loop (aoloop).
   Follows a call to aoread,parfile.
   Keywords:
   disp: set to display stuff
   forcemat: set to force measuring new iMat and computing new cMat.
   keepdmconfig: when forcemat=1, the default behavior is to agregate the
   previously determined active and extrapolated actuators and do the
   new interacton matrix with all of them. The selection of valid/extrapolated
   is then reset (and re-executed when the interaction matrix is done).
   Setting keepdmconfig=1 impose that the valid/extrapolated setting remains
   as it was.

   This routine uses:
   - the ao structures (atm, wfs, ...)

   This routine calls:
   - a bunch of initialization routines

   This routine sets:
   - a bunch of things.

   SEE ALSO: aoread, aoloop
*/
{
  extern pupil,ipupil;
  extern iMat,cMat;
  extern modalgain;
  extern _n,_n1,_n2,_p1,_p2,_p;
  extern def,mircube;
  extern tip1arcsec,tilt1arcsec;

  if ((disp==[])&&(aoinit_disp!=[])) disp=aoinit_disp;
  if ((clean==[])&&(aoinit_clean!=[])) clean=aoinit_clean;
  if ((forcemat==[])&&(aoinit_forcemat!=[])) forcemat=aoinit_forcemat;
  if ((svd==[])&&(aoinit_svd!=[])) svd=aoinit_svd;
  if ((keepdmconfig==[])&&(aoinit_keepdmconfig!=[])) keepdmconfig=aoinit_keepdmconfig;
  
  sphase = bphase = mircube = [];

  hcp_file,parprefix+"init.ps",ps=1;
  
  //=====================================
  // PARAMETER CHECKS. SETS SOME DEFAULTS
  //=====================================
  checkParameters;
  if (!is_set(disp)) {disp = 0;}
  if (!is_set(dpi)) {dpi = 60;}
  if (is_set(clean)) {forcemat=1;}


  // Sets other parameters:
  sim._size = int(2^ceil(log(sim.pupildiam)/log(2)+1));
  size	    = sim._size;
  // mircube will receive the shape of each DM (one plan per DM):
  mircube   = array(float,sim._size,sim._size,ndm);

  // INITIALIZE OUTPUT RESULT FILE parprefix.RES

  if (!fileExist(parprefix+".res")) {
    if (sim.verbose>=1) {
      write,">> File "+parprefix+".res not found, creating one...\n";
    }
    f = open(parprefix+".res","w");
    close,f;
  }
  f	= open(parprefix+".res","a+");
  write,f,format="=============================\n%s\n",sim.name;
  close,f;


  //===================================================================
  // INITIALIZE SOME STUFF FOR OFF ZENITH CONFIGURATIONS:
  wfs._gsalt   = wfs.gsalt / cos(gs.zenithangle*dtor);
  wfs._gsdepth = wfs.gsdepth / cos(gs.zenithangle*dtor);
  // dm.alt is not modified by the zenith angle in the system I simulate.
  atm._layeralt= &(*atm.layeralt / cos(gs.zenithangle*dtor));

  //===================================================================
  // INITIALIZE PUPIL
  // set pupil center: If SH wfs, should be between four central pixels
  // except if nbsub is odd and npixpersub is odd.
  // 02/20/03: FIX ME: Should rework this part to allow for mixture
  //                   of curvature and SH.
  //===================================================================

  cent = [];

  for (i=1;i<=nwfs;i++) {
    wfs(i)._cyclecounter = 1;
    if (wfs(i).type == "hartmann") {
      npixpersub	= float(sim.pupildiam)/wfs(i).shnxsub;
      if (npixpersub != int(npixpersub)) {
        write,"sim.pupildiam should be a multiple of wfs.shnxsub !";
        return -1;
      }
      if (odd(wfs(i).shnxsub) && odd(npixpersub)) {
        grow,cent,sim._size/2+1;
      } else {
        grow,cent,sim._size/2+0.5;
      }
    } else if (wfs(i).type == "curvature") {
      grow,cent,sim._size/2+1;
    } else if (wfs(i).type == "pyramid") {
      grow,cent,sim._size/2+0.5;
    }
  }
  if (anyof(cent != cent(1))) {
    write,"Wrong mix of hartmann/curvature or subaperture odd/even !";
    write,"I can't handle that as some of your selected sensor require";
    write,"the pupil to be centered on a pixel and others in between 2 pixels";
    write,"Sorry :-(";
    exit;
  }
  sim._cent = cent(1);

  // Initialize pupil array:

  pupil	 = float(MakePupil(sim._size,sim.pupildiam,xc=sim._cent,yc=sim._cent,real=1,cobs=tel.cobs));
  ipupil = float(MakePupil(sim._size,sim.pupildiam,xc=sim._cent,yc=sim._cent,cobs=tel.cobs));

  //==================================
  // DEFINE INDICES FOR SUBARRAY WORK:
  //==================================

  _p        = sim.pupildiam;
  _p1       = long(ceil(sim._cent-_p/2.));
  _p2       = long(floor(sim._cent+_p/2.));
  _p        = _p2-_p1+1;
  
  _n        = _p+4;
  _n1       = _p1-2;
  _n2       = _p2+2;

  //==================================
  // INITIALIZE DISPLAYS
  //==================================

  if (is_set(disp) || (sim.debug>0)) {
    if (!yaopy) { // set if using pygtk GUI, which prevents remapping a new window
      // not a problem if not using it.
      winkill,0;
      winkill,2;
      window,0,style="aosimul3.gs",dpi=dpi,width=long(550*(dpi/50.)),height=long(425*(dpi/50.)),wait=1;
    }
    disp2D,wfs._fimage,wfs.gspos(1,),wfs.gspos(2,),2,zoom=wfs.dispzoom,init=1;
  }

  //==================================
  // INITIALIZE PHASE SCREENS
  //==================================

  if (sim.verbose>=1) {write,"\n> INITIALIZING PHASE SCREENS";}
  getTurbPhaseInit;

  //==================================
  // INITIALIZE SENSOR
  //==================================

  if (sim.verbose>=1) {write,"\n> INITIALIZING SENSOR GEOMETRY";}

  for (n=1;n<=nwfs;n++) {

    if (wfs(n).type == "curvature") {

      // build subaperture geometry:
      MakeCurvWfsSubs,n,size,sim.pupildiam;
      // init WFS
      CurvWfs,pupil,pupil*0.0f,n,init=1;
      wfs(n)._nsub = int(sum(*(wfs(n).nsubperring)));
      wfs(n)._nmes = wfs(n)._nsub;

    } else if (wfs(n).type == "hartmann") {

      // init WFS:
      ShWfsInit,ipupil,n,imat=1;
      wfs(n)._nmes = 2*wfs(n)._nsub;

    } else if (wfs(n).type == "pyramid") {

      // init WFS (does not work for version 2)
      error,"Not Upgraded to version 2";
      v = PyramidWfs(pupil,pupil*0.,disp=disp,init=1);
      wfs(n).nsub = numberof(v)/2;
      wfs(n)._nmes = 2*wfs(n)._nsub;

    }
  }

  // set up array for uplink TT compensation:
  xy = indices(sim._size);
  fact = 4.848 * (tel.diam/sim.pupildiam);
  tip1arcsec = float(xy(,,1)*fact);
  tilt1arcsec = float(xy(,,2)*fact);
  
  //===============================
  // GET WFS REFERENCE MEASUREMENTS
  //===============================

  mircube *= 0.0f;
  // disable filtertilt is any
  mem = wfs.filtertilt; wfs.filtertilt *= 0n;

  for (ns=1;ns<=nwfs;ns++) {
    wfs(ns)._refmes = &(array(0.0f,wfs(ns)._nmes));
  }
  refmes = multWfsIntMat(disp=disp);
  wfs._refmes = splitWfsVector(refmes);

  wfs.filtertilt = mem;

  //============================================
  // GET WFS TIP AND TILT REFERENCE MEASUREMENTS
  //============================================

  // disable filtertilt is any
  mem = wfs.filtertilt; wfs.filtertilt *= 0n;

  // step per pixel to have a x" tilt:
  push = 0.05; // in arcsec

  // tip:
  mircube *= 0.0f;
  mircube(,,1) = float(tip1arcsec*push);

  mes = multWfsIntMat(disp=disp)/push;  // mes in arcsec

  wfs._tiprefv = splitWfsVector(mes);

  for (ns=1;ns<=nwfs;ns++) {
    wfs(ns)._tiprefvn = &(*wfs(ns)._tiprefv/sum(*wfs(ns)._tiprefv^2.));
  }

  // now the tip ref vector are normalized, so to compute the tip content
  // one has just to do sum(vector * tiprefv)

  // tilt:
  mircube(,,1) = float(tilt1arcsec*push);

  mes = multWfsIntMat(disp=disp)/push;  // mes in arcsec

  wfs._tiltrefv = splitWfsVector(mes);

  for (ns=1;ns<=nwfs;ns++) {
    wfs(ns)._tiltrefvn = &(*wfs(ns)._tiltrefv/sum(*wfs(ns)._tiltrefv^2.));
  }

  // restore pre-operation filtertilt:
  wfs.filtertilt = mem;

  //==================================
  // INITIALIZE DM INFLUENCE FUNCTIONS
  //==================================

  if (sim.verbose>=1) {write,"\n> Initializing DM influence functions";}

  // loop over DMs:
  for (n=1;n<=ndm;n++) {

    // Set _n1 and _n2, the limit indices
    if (dm(n).type == "stackarray") {
      // find out the support dimension for the given mirror.
      extent = dm(n).pitch*(dm(n).nxact+2.); // + 1.5 pitch each side      
      dm(n)._n1 = long(clip(floor(sim._cent-extent/2.),1,));
      dm(n)._n2 = long(clip(ceil(sim._cent+extent/2.),,sim._size));
    } else {  // we are dealing with a curvature mirror, TT, zernike or aniso:
      dm(n)._n1 = 1;
      dm(n)._n2 = sim._size;
    }
    //special case = only 2 pixels margin each side:
    if (dm(n).alt == 0) {
      extent=sim.pupildiam+4;
      dm(n)._n1 = long(clip(floor(sim._cent-extent/2.),1,));
      dm(n)._n2 = long(clip(ceil(sim._cent+extent/2.),,sim._size));
    }
    
    // compute influence functions:
    // If file exist, read it out:
    if ( (fileExist(dm(n).iffile)) && (!is_set(clean)) ) { 

      if (sim.verbose>=1) {
        write,format="  >> Reading file %s\n",dm(n).iffile;
      }
      dm(n)._def = &(float(fitsRead(dm(n).iffile)));
      dm(n)._nact = dimsof(*(dm(n)._def))(4);
      if ( dm(n).type == "stackarray" ) {
        dm(n)._x = &(fitsRead(dm(n).iffile,hdu=1));
        dm(n)._y = &(fitsRead(dm(n).iffile,hdu=2));
        if (dm(n).elt == 1) {
          dm(n)._eltdefsize = dimsof(*(dm(n)._def))(2);
          dm(n)._i1 = &(int(fitsRead(dm(n).iffile,hdu=3)));
          dm(n)._j1 = &(int(fitsRead(dm(n).iffile,hdu=4)));
        }
      }

      if ( (fileExist(dm(n)._eiffile)) && (!is_set(clean)) ) { 
        if (sim.verbose>=1) {
          write,format="  >> Reading extrapolated actuators file %s\n",dm(n)._eiffile;
        }
        dm(n)._edef = &(float(fitsRead(dm(n)._eiffile)));
        dm(n)._enact = dimsof(*(dm(n)._edef))(4);
        if ( dm(n).type == "stackarray" ) {
          dm(n)._ex = &(fitsRead(dm(n)._eiffile,hdu=1));
          dm(n)._ey = &(fitsRead(dm(n)._eiffile,hdu=2));
          if (dm(n).elt == 1) {
            dm(n)._ei1 = &(int(fitsRead(dm(n)._eiffile,hdu=3)));
            dm(n)._ej1 = &(int(fitsRead(dm(n)._eiffile,hdu=4)));
          }
        }
      }
      
    } else { // else compute the influence functions:

      if (sim.verbose>=1) {
        write,format="  >> Computing Influence functions for mirror # %d\n",n;
      }
      if (dm(n).type == "bimorph") {MakeBimorphIF, n, disp=disp,cobs=tel.cobs;}
      if (dm(n).type == "stackarray") {
        if (dm(n).elt == 1) {
          MakeEltPztIF, n, disp=disp;
        } else {
          MakePztIF, n, disp=disp;
        }
      }
      if (dm(n).type == "zernike") {MakeZernikeIF, n, disp=disp;}
      if (dm(n).type == "tiptilt") {MakeTipTiltIF, n, disp=disp;}
      if (dm(n).type == "aniso") {MakeAnisoIF, n, disp=disp;}
      // the IF are in microns/volt
      if (sim.verbose>=1) {
        write,format="\n  >> I.F. stored in %s\n",dm(n).iffile;
      }
      fitsWrite,dm(n).iffile,*(dm(n)._def);
      if ( dm(n).type == "stackarray" ) {
        fitsWrite,dm(n).iffile,*(dm(n)._x),exttype="IMAGE",append=1;
        fitsWrite,dm(n).iffile,*(dm(n)._y),exttype="IMAGE",append=1;
        if (dm(n).elt == 1) {
          fitsWrite,dm(n).iffile,long(*(dm(n)._i1)),exttype="IMAGE",append=1;
          fitsWrite,dm(n).iffile,long(*(dm(n)._j1)),exttype="IMAGE",append=1;
        }
      }
    }
  }

  //=========================================
  // DO INTERACTION MATRIX WITH ALL ACTUATORS
  //=========================================

  if (!(allof(fileExist(dm.iffile)) && fileExist(mat.file) && (forcemat != 1))) {


    if (!is_set(keepdmconfig)) { // concatenate dm._def and dm._edef
      for (nm=1;nm<=ndm;nm++) {  // loop on DMs
        if ((*dm(nm)._edef) != []) {  // if _edef == [], there is no extrap. act. defined
          dm(nm)._def = &(_(*dm(nm)._def,*dm(nm)._edef));
          dm(nm)._x = &(_(*dm(nm)._x,*dm(nm)._ex));
          dm(nm)._y = &(_(*dm(nm)._y,*dm(nm)._ey));
          if (dm(nm).elt == 1) {
            dm(nm)._i1 = &(int(_(*dm(nm)._i1,*dm(nm)._ei1)));
            dm(nm)._j1 = &(int(_(*dm(nm)._j1,*dm(nm)._ej1)));
          }
          dm(nm)._nact = dimsof(*(dm(nm)._def))(4);
          dm(nm)._edef = dm(nm)._ex = dm(nm)._ey = &([]);
          dm(nm)._enact = 0;
          if (dm(nm).elt == 1) { dm(nm)._ei1 = dm(nm)._ej1 = &([]); }
        }
      }
    }
    
    if (sim.verbose >= 1) {write,"\n> DOING INTERACTION MATRIX";}
    iMat = array(double,sum(wfs._nmes),sum(dm._nact));
    cMat = transpose(iMat);
    
    // measure interaction matrix:
    doInter,disp=disp;

    // select valid actuators by their response, only if
    // 1. at least one of the DM is a stackarray (otherwise it does not make sense)
    // 2. user has not selected "keepdmconfig"
    // otherwise we proceed without filtering any actuators
    if (anyof(dm.type == "stackarray") && !is_set(keepdmconfig)) {

      indexDm       = array(long,[2,2,ndm]);
      indexDm(,1)   = [1,dm(1)._nact];
      for (nm=2;nm<=ndm;nm++) {
        indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
      }

      // response from actuators:
      resp = sqrt((iMat^2.)(sum,));
      resp = (abs(iMat))(max,);

      if (sim.debug >= 2) {fma; plh,resp; limits,square=0; typeReturn;}

      for (nm=1;nm<=ndm;nm++) {

        // filter actuators only in stackarray mirrors:
        if (dm(nm).type == "stackarray") {

          dmx = *dm(nm)._x;
          dmy = *dm(nm)._y;
          if (dm(nm).elt == 1) {
            dmi1 = *dm(nm)._i1;
            dmj1 = *dm(nm)._j1;
          }

          tmp = resp(indexDm(1,nm):indexDm(2,nm));

          if (dm(nm).thresholdresp < 0.) {
            do {
              // criteria for an actuator to be retained as valid:
              ok = where(tmp >  dm(nm).thresholdresp*max(tmp));
              nok= where(tmp <= dm(nm).thresholdresp*max(tmp));
              dm(nm)._x = &(dmx(ok));
              dm(nm)._y = &(dmy(ok));
              if (numberof(nok) != 0) {
                dm(nm)._ex = &(dmx(nok));
                dm(nm)._ey = &(dmy(nok));
              }
              graphicConfig,dm(nm).subsystem,nm;
              write,format="DM #%d: # of valid actuators: %d ",nm,numberof(ok);
              again = "y"; read,prompt="Again ? ",again;
              if (again != "n") {
                read,prompt=swrite(format="Threshold (old = %f) : ",\
                                   dm(nm).thresholdresp),dm(nm).thresholdresp;
              }
            } while (again != "n");
          }

          ok = where(tmp >  dm(nm).thresholdresp*max(tmp));
          nok= where(tmp <= dm(nm).thresholdresp*max(tmp));
          dm(nm)._x = &(dmx(ok));
          dm(nm)._y = &(dmy(ok));
          if (dm(nm).elt == 1) {
            dm(nm)._i1 = &(int(dmi1(ok)));
            dm(nm)._j1 = &(int(dmj1(ok)));
          }
          
          if (numberof(nok) == 0) continue; //yes, it should be here.

          dm(nm)._ex = &(dmx(nok));
          dm(nm)._ey = &(dmy(nok));
          if (dm(nm).elt == 1) {
            dm(nm)._ei1 = &(int(dmi1(nok)));
            dm(nm)._ej1 = &(int(dmj1(nok)));
          }

          dm(nm)._edef = &((*(dm(nm)._def))(,,nok));
          dm(nm)._def = &((*(dm(nm)._def))(,,ok));

          if (sim.verbose>=1) {
            write,format="  DM #%d: # of valid actuators: %d. "+
              "(I got rid of %d actuators after iMat)\n",
              nm,numberof(ok),numberof(nok);
            write,format="  >> valid I.F. stored in %s\n",dm(nm).iffile;
          }

          // rewrite influence function file:
          fitsWrite,dm(nm).iffile,*(dm(nm)._def);
          fitsWrite,dm(nm).iffile,*(dm(nm)._x),exttype="IMAGE",append=1;
          fitsWrite,dm(nm).iffile,*(dm(nm)._y),exttype="IMAGE",append=1;
          if (dm(nm).elt == 1) {
            fitsWrite,dm(nm).iffile,long(*(dm(nm)._i1)),exttype="IMAGE",append=1;
            fitsWrite,dm(nm).iffile,long(*(dm(nm)._j1)),exttype="IMAGE",append=1;            
          }
          dm(nm)._nact = (dimsof(*(dm(nm)._def)))(4);

          // write extrapolated actuator influence functions file:
          fitsWrite,dm(nm)._eiffile,*(dm(nm)._edef);
          fitsWrite,dm(nm)._eiffile,*(dm(nm)._ex),exttype="IMAGE",append=1;
          fitsWrite,dm(nm)._eiffile,*(dm(nm)._ey),exttype="IMAGE",append=1;
          if (dm(nm).elt == 1) {
            fitsWrite,dm(nm)._eiffile,long(*(dm(nm)._ei1)),exttype="IMAGE",append=1;
            fitsWrite,dm(nm)._eiffile,long(*(dm(nm)._ej1)),exttype="IMAGE",append=1;            
          }
          dm(nm)._enact = (dimsof(*(dm(nm)._edef)))(4);

          if (sim.debug >= 1) {
            tv,compDmShape(nm,&(array(1.0f,dm(nm)._nact)));
            pltitle,swrite(format="Pushing on all actuator of DM#%d",nm);
            typeReturn;
          }

          iMat(,indexDm(1,nm)-1+nok) *=0.;

        }
      }
      //      stop;
      iMat = iMat(,where(iMat(rms,) != 0.));
      cMat = transpose(iMat)*0.;
    }
  }

  //=========================================
  // LOAD OR COMPUTE THE EXTRAPOLATION MATRIX
  //=========================================
  // Normally, I don't load the valid to extrapolation actuator file.
  // It is NOT saved anyway.
  // it's fast enough to compute it here from the existing data (position
  // and number of extrapolation actuators).
  // However, if the user want to use its own extrapolation method (the one
  // here, although it works fine, is pretty dumb), I provide a way to do it
  // by providing a file name (dm(nm).ecmatfile). If this file exist, I will
  // use it.

  for (nm=1;nm<=ndm;nm++) {

    if (dm(nm).type != "stackarray") break;  // possible but not implemented.

    if (dm(nm).noextrap == 1) { // we have to get rid of any extrapolated actuator now
      dm(nm)._edef = dm(nm)._ex = dm(nm)._ey = &([]);
      if (dm(nm).elt == 1) {
        dm(nm)._ei1 = dm(nm)._ej1 = &([]);
      }
      dm(nm)._enact = 0;
      break;
    }

    if (dm(nm)._enact == 0) break;  // no extrapolated actuator

    if (fileExist(dm(nm).ecmatfile)) {

      // if defined and exist, we read it and proceed
      if (sim.verbose >= 1) {
        write,format="Reading valid to extrap. matrix %s\n",dm(nm).ecmatfile;
      }
      dm(nm)._extrapcmat = fitsRead(dm(nm).ecmatfile);

    } else {

      // we compute the valid to extrap matrix:
      if (sim.verbose >= 1) {
        write,format="Computing valid to extrap. matrix for DM#%d\n",nm;
      }
      tx = *dm(nm)._x;
      ty = *dm(nm)._y;
      tex = *dm(nm)._ex;
      tey = *dm(nm)._ey;
      disact = sqrt( (tx-tex(-,))^2. + (ty-tey(-,))^2. ); // distance from extrap to valid
      disact = disact/dm(nm).pitch; // now in pitch unit (1,2,..)
      disact = exp(-(disact/0.7)^2.);  // "influence function" for extrap. command
      disact = disact / disact(sum,)(-,);
      dm(nm)._extrapcmat = &(transpose(disact));
    }
  }
    
  
  /* STILL TO BE IMPLEMENTED IN V1.2
     I'll keep it for later.
     // INITIALIZE NON-COMMON PATH REFERENCE MEASUREMENTS:
     if (wfs.type == "hartmann") {
     if (abs(ao.refZernikes)(sum) != 0) {
     prepzernike,ao._size,ao.PupilDiam,ao._cent,ao._cent;
     phase = ipupil*0.;
     for (i=1;i<=numberof(ao.refZernikes);i++) {
     phase += ao.refZernikes(i)*zernike_ext(i+1);}
     noise = ao.WfsNoise; ao.WfsNoise = 0.;
     WfsRefMes = ShWfs(ipupil,phase);
     ao.WfsNoise = noise;
     if (ao.verbose>=1) {write,format="WfsRefMes rms = %f, Min = %f and Max = %f\n",
     WfsRefMes(rms),min(WfsRefMes),max(WfsRefMes);}
     } else {
     WfsRefMes = array(float,ao._WfsNMes);
     if (ao.verbose>=1) {write,"No WFS reference used\n";}
     }
     }
  */ 

  // INITIALIZE MODAL GAINS:

  if (sim.verbose>=1) {write,"\n> INITIALIZING MODAL GAINS";}

  modalgain = [];

  if (fileExist(loop.modalgainfile)) {

    if (sim.verbose>=1) {write,format=" >> Reading file %s\n\n",loop.modalgainfile;}
    modalgain = fitsRead(loop.modalgainfile);

  }
  
  if (numberof(modalgain) != sum(dm._nact)) {

    if (sim.verbose>=1) {
      write,format="I did not find %s or it did not have the right\n ",
        loop.modalgainfile;
      write,format="  number of elements (should be %i), so I have generated\n",
        sum(dm._nact);
      write,"  a modal gain vector with ones in it.\n";
    }

    modalgain = array(1.,sum(dm._nact));

  }

  // INITIALIZE COMMAND MATRIX:

  if (sim.verbose>=1) {write,"\n> INTERACTION AND COMMAND MATRICES";}

  if ((fileExist(mat.file)) && (forcemat != 1)) {

    if (sim.verbose>=1) {write,format="  >> Reading file %s\n",mat.file;}
    // read out mat file and stuff iMat and cMat:
    tmp = fitsRead(mat.file);
    iMat = tmp(,,1);
    cMat = transpose(tmp(,,2));
    tmp = [];

  }

  // in the opposite case, plus if svd=1 (request re-do SVD):
  if ((!fileExist(mat.file)) || (forcemat == 1) || (svd == 1)) {
    
    if (sim.verbose>=1) {
      write,">> Preparing SVD and computing command matrice";
    }
    // do the SVD and build the command matrix:
    sswfs = ssdm = [];
    for (ns=1;ns<=nwfs;ns++) {grow,sswfs,array(wfs(ns).subsystem,wfs(ns)._nmes);}
    for (nm=1;nm<=ndm;nm++)  {grow,ssdm,array(dm(nm).subsystem,dm(nm)._nact);}
    for (ss=1;ss<=max(dm.subsystem);ss++) {

      wsswfs = where(sswfs == ss);
      wssdm  = where(ssdm  == ss);
      imat = iMat(wsswfs,wssdm);
      mg = modalgain(wssdm);

      prepSVD,imat,disp=disp,subs=ss,svd=svd;

      tmpdisp = noneof(dm(where(dm.subsystem == ss)).type == "aniso");
      cmat = buildComMat( (*mat.condition)(ss), mg, subsystem=ss,
                          all=1,disp=tmpdisp);
      cMat(wssdm,wsswfs) = cmat;
      
    }

    // More debug display
    if (sim.debug>=3) {
      tv,cMat(,+)*iMat(+,);
      pltitle,"cMat(,+)*iMat(+,)";
      typeReturn;
      tv,iMat(,+)*cMat(+,);
      pltitle,"iMat(,+)*cMat(+,)";
      typeReturn;
    }
      
    // save the results:
    fitsWrite,mat.file,[iMat,transpose(cMat)];

  }

  //===================================================================
  // COMPUTE THE COMMAND VECTOR FOR OFFLOADING THE ANISOPLANATISM MODES
  //===================================================================

  if (anyof(dm.type == "aniso")) {

    // find which DM is the aniso DM:
    nmaniso = where(dm.type == "aniso");
    if (numberof(nmaniso) != 1) {
      error,"there can be only one aniso DM !";
    }

    // finds which DM is at altitude 0
    w0 = where(dm.alt == 0);
    nmlow = where( (dm(w0).type == "stackarray") | (dm(w0).type == "bimorph") |
                   (dm(w0).type == "zernike") );
    if (numberof(nmlow) == 0) {
      error,"I can not find a DM at altitude 0 to produce the lower "+
        "part of the anisoplanatism modes !";
    }
    if (numberof(nmlow) > 1) {
      write,format="Weird. There are %d high-order DMs at altitude=0:",numberof(nmlow);
      for (i=1;i<=numberof(nmlow);i++) {
        write,format="DM #%d:",nmlow(i);
        print,dm(i);
      }
      rep = "";
      question = "Which one ? ["+strcompress(strjoin(swrite(nmlow),","),all=1)+"] ";
      read,prompt=question,rep;
      tmp=0; sread,rep,tmp;
      if (noneof(nmlow == tmp)) { error,"Invalid selection"; }
      nmlow = tmp;
    }

    // finds which DM is at altitude specified by dm(nmaniso).alt
    wn0 = where(dm.alt == dm(nmaniso).alt);
    nmhigh = where( (dm(wn0).type == "stackarray") | (dm(wn0).type == "bimorph") |
                    (dm(wn0).type == "zernike") );
    if (numberof(nmhigh) == 0) {
      error,"I can not find a DM at the requested altitude to produce the higher "+
        "part of the anisoplanatism modes !";
    }
    if (numberof(nmhigh) > 1) {
      write,format="Weird. There are %d high-order DMs at altitude %.0f:",
        numberof(nmhigh),dm(nmaniso).alt;
      for (i=1;i<=numberof(nmhigh);i++) {
        write,format="DM #%d:",nmhigh(i);
        print,dm(i);
      }
      rep = "";
      question = "Which one ? ["+strcompress(strjoin(swrite(nmhigh),","),all=1)+"] ";
      read,prompt=question,rep;
      tmp=0; sread,rep,tmp;
      if (noneof(nmhigh == tmp)) { error,"Invalid selection"; }
      nmhigh = tmp;
    }
    
    projectAnisoIF,nmaniso(1),nmlow(1),nmhigh(1),disp=0;
  }
  
  //==========================
  // OUTPUT GRAPHIC FOR CONFIG
  //==========================

  if ((disp == 1)&&(!yaopy)) { graphicConfig; }
  hcp_finish;
  
  //===================================
  // PRINT OUT SUMMARY FOR WFSs AND DMs
  //===================================

  if ( sim.verbose > 0) {
    write,"";
    write,format="%s:\n","Summary";
    for (nm=1;nm<=ndm;nm++) {
      write,format="Mirror #%1d, %s, %d actuators, conjugated @ %.0f m\n",
        nm,dm(nm).type,dm(nm)._nact,dm(nm).alt;
    }
    for (ns=1;ns<=nwfs;ns++) {
      write,format="WFS #%2d, %s (meth. %d), %2d subap., offaxis "+
        "(%+4.1f\",%+4.1f\"), noise %s\n",
        ns,wfs(ns).type,wfs(ns).shmethod,wfs(ns)._nsub,wfs(ns).gspos(1),
        wfs(ns).gspos(2),
        ( (wfs(ns).noise && (wfs(ns).shmethod == 2)) ? "enabled" : "disabled" );
    }
    write,format="D/r0 (500nm) = %.1f; %d iterations\n",atm.dr0at05mic/
      cos(gs.zenithangle*dtor)^0.6,loop.niter;
  }
  
  // same in result file:
  f	= open(parprefix+".res","a+");
  write,f,"";
  write,f,format="%s:\n","Summary";
  for (nm=1;nm<=ndm;nm++) {
    write,f,format="Mirror #%1d, %s, %d actuators, conjugated @ %.0f m\n",
      nm,dm(nm).type,dm(nm)._nact,dm(nm).alt;
  }
  for (ns=1;ns<=nwfs;ns++) {
    write,f,format="WFS #%2d, %s (meth. %d), %2d subap., offaxis "+
      "(%+4.1f\",%+4.1f\"), noise %s\n",
      ns,wfs(ns).type,wfs(ns).shmethod,wfs(ns)._nsub,wfs(ns).gspos(1),
      wfs(ns).gspos(2),
      ( (wfs(ns).noise && (wfs(ns).shmethod == 2)) ? "enabled" : "disabled" );
  }
  write,f,format="D/r0 (500nm) = %.1f; %d iterations\n",atm.dr0at05mic/
    cos(gs.zenithangle*dtor)^0.6,loop.niter;
  close,f;
  
}

//---------------------------------------------------------------
func aoloop(disp=,savecb=,dpi=,controlscreen=,nographinit=,gui=)

/* DOCUMENT func aoloop(disp=)
   Belongs to the yao package suite.
   The parameters are entered in a parameter file, called for instance
   "sh12.par". The sequence for running an AO simulation is -to date-
   as follow:
     
   > aoread,"sh12.par"
   > aoinit,disp=1,forcemat=1  (or forcemat=0 if IF and matrices already
   computed)
   > psf = aoloop(disp=1)      (this routine)
     
   This simulates the ao loop. After the initialization, the loop
   steps include:
   - get the turbulent phase from getTurbPhase()
   - subtract the mirror figure (previous iteration)
   - get the WFS measurements (SH or curvature)
   - computes the command vector from command matrix
   with a pure integrator with gain
   - computes the mirror shape
   - computes the PSF
   - accumulate the results if the statistics flag allows it.
   - displays and writes out the results

   I have changed the structure of this routine on June 11, 2002 to
   include the getTurbPhase function. One of the goal of this new
   structure is to allow for better statistics in doing short time
   sequences separated by large gaps. In other words, the time sequence
   is as follow:
   - Start the close loop
   - do ao.LoopStartSkip iteration without accumulating statistics,
   just to allow convergence
   - Start accumulating statistics
   - after ao.LoopSkipEvery steps, jump in time/space by
   ao.LoopSkipBy steps
   - Continue from there, starting by skiping ao.LoopStartSkip steps,
   and then accumulate before the next jump
   This has the great advantage that it provides a much faster convergence
   of the results, as you flip through completely uncorrelated phase
   samples quite fast, instead of having to let them pass at loop speed
   (which is necessarily quite slow as one loop step is typically 1ms).
   It works really nicely and output PSF are much smoother annd
   statistically significant. Also of course, this technique allow to
   account for servolag error (if not, one could get independant phase
   realization for each step and the convergence would be even faster).
   However, this means that you "loose" ao.LoopStartSkip steps each
   time you jump. Also, ao.LoopStartSkip has to be large enough that
   you don't get the tail of the convergence curve, else you will
   get a biased estimate (lower performance than with the regular
   continuous loop scheme)
*/
{
  //  extern mircube,sphase,bphase;
  //  extern fwhm, strehl, e50;
  //  extern ditherMes,sairy,niterok;
  //  extern imav,im;
  extern looptime, mircube, command, wfsMesHistory, cubphase;
  extern im, imav, imtmp, airy, sairy, fwhm, e50;
  extern niterok;
  extern pp, sphase, bphase, imphase, dimpow2, cMat, pupil, ipupil;
  extern time, strehllp, strehlsp, itv, ok, njumpsinceswap;
  extern remainingTimestring  ;
  extern cbmes, cbcom, cberr;
  extern indexDm, waniso, wdmaniso;
  extern ditherPeriod, ditherAmp, ditherGain, cggain, ditherMes;
  extern ditherCosLast, ditherSinLast, ditherMesCos, ditherMesSin;
  extern dispImImav;
  extern loopCounter, nshots;
  extern dispFlag, savecbFlag, dpiFlag, controlscreenFlag, nographinitFlag;

  dispFlag = disp;
  savecbFlag = savecb;
  dpiFlag = dpi;
  controlscreenFlag = controlscreen;
  nographinitFlag = nographinit;


  // Initialize displays:
  if (!is_set(disp)) {disp = 0;}
  if (!is_set(controlscreen)) {controlscreen = 0;}
  if (is_set(disp) && !is_set(nographinit)) {
    if (!yaopy) {
      if (!is_set(dpi)) {
        window,0,style="aosimul3.gs",wait=1;
      } else {
        winkill,0;
        window,0,style="aosimul3.gs",dpi=dpi,width=long(550*(dpi/50.)),
          height=long(425*(dpi/50.)),wait=1;
      }
    }
  }
  if (is_set(controlscreen) && !is_set(nographinit)) {
    controlScreen,0,init=1;
  }
  
  size		= sim._size;

  
  // Some arrays initialization:
  looptime      = 0.;
  mircube	= array(float,[3,size,size,ndm]);
  command       = array(float,sum(dm._nact));
  wfsMesHistory = array(float,[2,sum(wfs._nmes),loop.framedelay+1]);
  cubphase      = array(float,[3,size,size,target._ntarget]);
  im	        = array(float,[3,size,size,target._ntarget]);
  imav	        = array(float,[4,size,size,target._ntarget,target._nlambda]);
  imtmp	        = array(float,[2,size,size]);
  airy	        = calcPSFVE(pupil,pupil*0.);
  sairy         = max(airy);
  fwhm = e50 = array(float,[2,target._ntarget,target._nlambda]);
  pp            = array(complex,[2,size,size]);
  sphase        = array(float,[2,_n,_n]);
  bphase        = array(float,[2,size,size]);
  imphase       = array(float,[2,size,size]);
  dimpow2       = int(log(size)/log(2));
  cMat          = float(cMat);
  pupil         = float(pupil);
  ipupil        = float(ipupil);
  time          = array(float,10);
  strehllp = strehlsp = itv = [];
  ok            = 0;
  niterok       = 0;
  remainingTimestring = "";
  njumpsinceswap = 0; // number of jumps since last screen swap
  if (is_set(savecb)) {
    cbmes         = array(float,[2,sum(wfs._nmes),loop.niter]);
    cbcom         = array(float,[2,sum(dm._nact),loop.niter]);
    cberr         = array(float,[2,sum(dm._nact),loop.niter]);
  }
  

  // Special: re-init the WFS, since it includes the computation of
  // the # of photons/WFS and the WFS bias anf flats, in case something
  // has changed since the aoinit. The # of photons in particular is one
  // parameter we often want to vary without having to go through the
  // whole aoinit again. For a curvature wfs, it re-inits the # of
  // photons and the extra-focal distance, so it is possible to loop
  // on "aoloop" with various l without the need for a aoinit (although
  // you would still be using the cMat determined for the original l,
  // but that might be on purpose.
  for (ns=1;ns<=nwfs;ns++) {
    // define _dispimage
    wfs(ns)._dispimage = &(*wfs(ns)._fimage*0.0f);
    if (wfs(ns).type == "hartmann") {
      ShWfsInit,ipupil,ns,silent=1;
    } else if (wfs(ns).type == "curvature") {
      CurvWfs,pupil,pupil*0.0f,ns,init=1,silent=1;
    }
  }
  // reset cyclecounters
  wfs._cyclecounter = wfs._cyclecounter*0 +1;
  
  indexDm       = array(long,2,ndm);
  indexDm(,1)   = [1,dm(1)._nact];
  for (nm=2;nm<=ndm;nm++) {
    indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
  }
  // initialize hysteresis if needed:
  for (nm=1;nm<=ndm;nm++) {
    dm(nm)._command = &(array(float,dm(nm)._nact));
    if (dm(nm).hyst > 0.) {
      *dm(nm)._command = hysteresis(*dm(nm)._command,nm,first=*dm(nm)._command);
    }
  }

  // Set up for anisoplatism modes
  aniso = anyof(dm.type == "aniso");
  if (aniso) {
    waniso = indexDm(1,where(dm.type == "aniso"))+[0,1,2];
    wdmaniso = where(dm.type == "aniso")(0);
    if (sim.verbose >=2) {
      print,"index of anisoplatism modes in command vector:",waniso;
    }
  }

  // dithering set up for upTT compensation (LGS only)
  wfs._centroidgain = array(1.f,nwfs);
  //  if (anyof(wfs.centGainOpt)) {
  // we dither if there is an uplink, filtertilt is set and centGainOpt is requested
  ditherPeriod = 8;
  ditherAmp = 0.05;  // supposed to be in arcsec
  ditherGain = 0.1;
  cggain = 0.1; // 0.01;
  ditherCosLast = ditherSinLast = 0.0;
  ditherMesCos = ditherMesSin = ditherMes = array(0.,nwfs);
  //}
  
  // verbose
  if (sim.verbose>=1) {
    write,format="\n> Starting loop with %i iterations\n",loop.niter;
  }

  // display
  if (is_set(disp)) {
    if (target._ntarget == 1) {
      // if only one target, display correct plate scale for image display
      *target.dispzoom = [(*target.lambda)(0)*1e-6/4.848e-6/
                          tel.diam*sim.pupildiam/sim._size/2*sim._size];
      // last 2 terms to accomodate internal dispzoom scaling
    }
    disp2D,im,*target.xposition,*target.yposition,1, 
      zoom=*target.dispzoom,init=1;
    disp2D,wfs._fimage,wfs.gspos(1,),wfs.gspos(2,),2,zoom=wfs.dispzoom,
      init=1;
  }

  dispImImav = 1; // 1 is to display im, 2 to display imav

  olddisp = disp;
  oldcontrolscreen = controlscreen;
  
  // Main loop:


  loopCounter=0;
  nshots = -1;
  if (gui) {
    require,"tyk.i";
    tyk_source,"/Users/frigaut/yorick-2.1/contrib/yao/yao.tcl";
  }
}

func go(nshot)
{
  extern niterok, starttime, now, nshots;
  extern dispFlag, savecbFlag, dpiFlag, controlscreenFlag, nographinitFlag;
  extern cbmes, cbcom, cberr;

  if (nshot!=[]) nshots = nshot;
  disp = dispFlag;
  savecb = savecbFlag;
  dpi = dpiFlag;
  controlscreen = controlscreenFlag;
  nographinit = nographinitFlag;

  if (loopCounter==0) {
    // initialize timers
    tic,2; starttime = _nowtime(2); 
  }
  now = tac(2);

  loopCounter++;
  nshots--;

  gui_progressbar_frac,float(loopCounter)/loop.niter;
  gui_progressbar_text,swrite(format="%d out of %d iterations",loopCounter,loop.niter);
  
  if (savecb) comvec=[];
  
  if (loopCounter>loop.niter) error,"LoopCounter > loop.niter !";
  i=loopCounter;
  tic; time(1) = tac();

    
  prevOK  = ok;

  //==========================
  // DITHERING
  // act on wfs._upttcommand
  // sense on wfs._tt
  if (anyof(wfs.centGainOpt)) {
    wdither = where((wfs.correctUpTT == 1) & (wfs.centGainOpt ==1) & 
              (wfs.filtertilt == 1));
    if (numberof(wdither) != 0) { nditherarray = array(1.,numberof(wdither)); }
    ditherCosCurr = ditherAmp*cos(2*pi*i/ditherPeriod);
    ditherSinCurr = ditherAmp*sin(2*pi*i/ditherPeriod);
    // subtract last command and add new. the circular motion is known. radius=ditherAmp.
    wfs(wdither)._upttcommand -= [ditherCosLast,ditherSinLast]*nditherarray(-,);
    wfs(wdither)._upttcommand += [ditherCosCurr,ditherSinCurr]*nditherarray(-,);
    ditherCosLast = ditherCosCurr;
    ditherSinLast = ditherSinCurr;
    // do the temporal filtering on cos and sin measurements (see requirements documents
    // I wrote for tOSC regarding the process)
    ditherMesCos = (1-ditherGain)*ditherMesCos + ditherGain*wfs._tt(1,)*
      cos(2*pi*i/ditherPeriod);
    ditherMesSin = (1-ditherGain)*ditherMesSin + ditherGain*wfs._tt(2,)*
      sin(2*pi*i/ditherPeriod);
    ditherMes = 2.*sqrt(ditherMesCos^2.+ditherMesSin^2.);
    //      clip,ditherMes,ditherAmp/4,ditherAmp*4;
    //      ditherMes(wdither);
    //      wfs._tt;
    if (i > 30) {
      wfs(wdither)._centroidgain += cggain*(ditherAmp-ditherMes(wdither));
    }
  }
  // END OF DITHERING PART
  //===========================
    
  // get the turbulent phase:   
  // Do the wavefront sensing:
  WfsMes = multWfs(i,disp=disp);

  time(2) +=  tac();

  ok  = statsokvec(i);
    
  if ((prevOK-ok)==1) { // After a jump
    njumpsinceswap++;
    if (njumpsinceswap==loop.jumps2swapscreen) {
      if (sim.verbose) write,"Reset and screens swap";
      swapScreens;        // swap screens
      njumpsinceswap = 0; // reset "jumps since last swap" counter
    } else { write,"Reset"; }
    command *=0.0f; // bug noticed on 2006mar01:
    for (nm=1;nm<=ndm;nm++) *dm(nm)._command *=0.0f; // < this should work instead.
    mircube *=0.0f; wfsMesHistory *=0.0f;
  }

  // Handling frame delay (0 -> no frame delay at all, 1 -> regular
  // one frame delay of integrator, 2 -> integrator + one frame
  // computation, etc...
  wfsMesHistory  = roll(wfsMesHistory,[0,-1]);
  wfsMesHistory(,loop.framedelay+1) = WfsMes;
  usedMes        = float(wfsMesHistory(,1));

  time(3) += tac();
    
  // computes the actuator error vector from the measurements:
  err      = cMat(,+)*usedMes(+);

  // get the anisoplanatism mode coefficients and project it
  // in the actuator space
  if (aniso) {
    err -= comaniso(,+)*err(waniso)(+);
    err(waniso) = 0.;
    //      dm(wdmaniso)._command = invcomaniso(+,)*
  }
    
  time(4) += tac();

  // Computes the mirror shape using influence functions:
  for (nm=1; nm<=ndm; nm++) {

    if (dm(nm).type == "aniso") {continue;}

    n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = n2-n1+1;

    // this DM error:
    dmerr = err(indexDm(1,nm):indexDm(2,nm));

    // update the command vector for this DM:
    *dm(nm)._command -= loop.gain*dm(nm).gain*dmerr;
    
    if (dm(nm).hyst > 0.) {*dm(nm)._command = hysteresis(*dm(nm)._command,nm);}
      
    if (dm(nm).maxvolt != 0) {
      dm(nm)._command=&(clip(*dm(nm)._command,-dm(nm).maxvolt,dm(nm).maxvolt));
    }

    mircube(n1:n2,n1:n2,nm) = compDmShape(nm,dm(nm)._command);

    if (savecb) grow,comvec,*dm(nm)._command;

    // extrapolated actuators:
    if ((dm(nm)._enact != 0) && (dm(nm).noextrap == 0)) {
      ecom = float(((*dm(nm)._extrapcmat)(,+))*(*dm(nm)._command)(+));
      mircube(n1:n2,n1:n2,nm) += compDmShape(nm,&ecom,extrap=1);
    }
  }

  time(5) += tac();
    
  ok = ok*((i % loop.stats_every) == 0); // accumulate stats only every 4 iter.
  
  okdisp = ( is_set(disp) && (((i-1) % disp) == 0) );
  if (nshots>0) okdisp=is_set(disp);
  okcscreen = ( is_set(controlscreen) && (((i-1) % controlscreen) == 0) );
  if (is_set(controlscreen) && (i == loop.niter)) okcscreen=1;  // display at last iteration

  // Computes the instantaneous PSF:
  if (ok) {
    // compute integrated phases and fill phase cube
    for (jl=1;jl<=target._nlambda;jl++) {
      for (jt=1;jt<=target._ntarget;jt++) {
        cubphase(,,jt)  = getPhase2dFromDms(jt,"target") + \
          getPhase2dFromOptics(jt,"target") + \
          getTurbPhase(i,jt,"target");
      }
      // compute image cube from phase cube
      status = _calcPSFVE(&pupil,&cubphase,&im,dimpow2,
                       target._ntarget,float(2*pi/(*target.lambda)(jl)));
        
      // Accumulate statistics:
      imav(,,,jl) = imav(,,,jl) + im;
    }
    niterok += 1;
    grow,itv,i;
    grow,strehlsp,im(max,max,1)/sairy;
    grow,strehllp,imav(max,max,1,0)/sairy/(niterok+1e-5);
  }

  time(6) += tac();

  // Displays
  if (disp) {
    for (ns=1;ns<=nwfs;ns++ ) {
      // if cyclecounter = 1, a final image just got computed
      if (wfs(ns)._cyclecounter == 1) {*wfs(ns)._dispimage = *wfs(ns)._fimage;}
    }
  }
    
  if (okdisp) {
    fma;
    // PSF Images
    plt,sim.name,0.01,0.01,tosys=0;
    if (dispImImav == 1) {
      disp2D,im,*target.xposition,*target.yposition,1,power=0.5;
    } else {
      disp2D,imav,*target.xposition,*target.yposition,1,power=0.5;
    }
    for (j=1;j<=nwfs;j++) {
      plg,wfs(j).gspos(2),wfs(j).gspos(1),marker='\2',
        type="none",marks=1,color="red";
    }
    mypltitle,"Instantaneous PSFs",[0.,-0.005],height=12;
    myxytitles,"","arcsec",[0.02,0.02],height=12;

    // WFS spots
    if (!allof(wfs.shmethod ==1)) {
      disp2D,wfs._dispimage,wfs.gspos(1,),wfs.gspos(2,),2;
      mypltitle,"WFSs",[0.,-0.005],height=12;
    }
      
    // mirror surface
    plsys,3;
    limits,square=1;
    pli,mircube(,,1)*ipupil,1,0.,2,1.;
    for (nm=2;nm<=ndm;nm++) {pli,mircube(,,nm),nm,0.,nm+1,1.;}
    mypltitle,"DM(s)",[0.,0.008],height=12;

    // Strehl plots
    if (anyof(itv)) {
      plsys,4;
      plg,strehlsp,itv;
      plg,strehllp,itv,color="red";
      myxytitles,"Iterations",swrite(format="Strehl @ %.2f mic",        \
                                     (*target.lambda)(0)),[0.025,-0.01],height=9;
      range,0;
    }
    if (userPlot != []) userPlot,i,init=(i==1);  // execute user's plot routine if it exists.
  }

  if (okcscreen) { controlScreen,i; }

  time(7) += tac();

  // Fill the circular buffer
  if (is_set(savecb)) {
    cbmes(,i) = WfsMes;
    cbcom(,i) = comvec;
    cberr(,i) = err;  // err was overwritten by calcPSFVE, corrected 2007mar12
  }

  glt      = 0.02;
  if (i==1) {glt=1.;}
    
  looptime = (tac(2)-now)*glt + looptime*(1.-glt);
  remainingTime = float(loop.niter-i)*looptime;
  remainingTimestring = secToHMS(remainingTime);

  // Prints out some results:

  if (((looptime <= 2) && ((i % 500) == 1)) ||
      ((looptime > 2) && ((i % 20) == 1)) || (nshots>=0)) {
    write,"       Short expo. image  Long expos. image ";
    write,"Iter#  Max.S/Min.S/Avg.S  Max.S/Min.S/Avg.S  Time Left";
  }

  if ((looptime > 2) || ((i % 50) == 1) || (nshots>=0)) {
    msg = swrite(format="%5i  %5.3f %5.3f %5.3f  %5.3f %5.3f %5.3f  %s",
                 i,im(max,max,max)/sairy,min(im(max,max,))/sairy,
                 avg(im(max,max,))/sairy,imav(max,max,max,0)/sairy/(niterok+1e-5),
                 min(imav(max,max,,0))/sairy/(niterok+1e-5),
                 avg(imav(max,max,,0))/sairy/(niterok+1e-5),remainingTimestring);
    write,msg;
    gui_message,msg;
      
  }
  time(8) += tac();

  if (nshots==0) return;

//  if (loopCounter<loop.niter) after, 0.001, go;
  if (loopCounter<loop.niter) set_idler, go;
  else after_loop;
}

func stop {
  write,format="Stopping @ iter=%d\n",loopCounter;
  set_idler;
}

func cont {
  write,"Proceeding...";
  set_idler,go;
}

func restart {
  extern imav, ditherMesCos, ditherMesSin;
  extern niterok, loopCounter;
  extern command, mircube, wfsMesHistory;
  extern strehllp,strehsp, itv;
  extern startime;
  loopCounter = 0;
  imav *= 0.0f;
  wfs._upttcommand *= 0.0f;
  wfs._centroidgain *= 0.0f+1.0f;
  ditherMesCos = ditherMesSin = 0.;
  command *=0.0f; mircube *=0.0f; wfsMesHistory *=0.0f;
  for (nm=1;nm<=ndm;nm++) {*dm(nm)._command *= 0.0f;}
  strehllp = strehlsp = itv = [];
  tic,2; starttime = _nowtime(2);
  niterok = 0;
  //          set_idler,go;
}
func disptoggle { dispFlag=1-dispFlag; }
func whereat { write,format="@ iter=%d\n",loopCounter; }
func do_prompt { write,format="%s ",">"; }

func after_loop(void)
{
  extern strehllp,strehsp, itv;
  extern cbmes, cbcom, cberr;
  extern strehl,e50,fwhm;

  savecb = savecbFlag;
  
  time = (time - roll(time,1))/loopCounter;
  timeComments = ["WF sensing","Reset and measurement history handling","cMat multiplication",\
                  "DM shape computation","Target PSFs estimation","Displays",\
                  "Circular buffers, end-of-loop printouts"];
  for (i=2;i<=8;i++) {
    write,format="time(%d-%d) = %5.2f ms  (%s)\n",i-1,i,time(i)*1e3,timeComments(i-1);}
  
  write,format="Finished on %s\n",gettime();
  tic,2; tottime = (_nowtime(2) - starttime);
  write,format="%f iterations/second in average\n",loop.niter/tottime;
  
  // Save the circular buffers:
  if (is_set(savecb)) {
    fitsWrite,"cbmes.fits",cbmes;
    fitsWrite,"cbcom.fits",cbcom;
    fitsWrite,"cberr.fits",cberr;
    write,"cbmes, cbcom and cberr are saved.";
    write,"You can run modalGainOptimization() to optimize and update the gains";
  }
  
  // End of loop calculations (fwhm, EE, Strehl):
  fairy	 = findfwhm(airy);
  e50airy= encircled_energy(airy,ee50);
  strehl = imav(max,max,,)/sairy/(niterok+1e-5);

  // number of corrected modes from strehl
  dr0l   = atm.dr0at05mic/cos(gs.zenithangle*dtor)^0.6*(0.5/(*target.lambda))^1.2;
  if (dr0l(0) > 0.) {
    nmodes = (-log(clip(strehl(,0),0.,0.99))/0.2944/dr0l(0)^1.6666)^(-1./0.85);
  } else {
    nmodes = strehl(,0)*0.;
  }
  psize  = (float(sim.pupildiam)/sim._size)*(*target.lambda)/tel.diam/4.848e-3;

  write,format="\n         lambda   XPos   YPos  FWHM[mas]  Strehl  E50d[mas]  #modes comp.%s\n","";
  for (jl=1;jl<=target._nlambda;jl++) {
    for (jt=1;jt<=target._ntarget;jt++) {
      fwhm(jt,jl) = findfwhm(imav(,,jt,jl),psize(jl));
      encircled_energy,imav(,,jt,jl),tmp;
      e50(jt,jl)  = tmp*psize(jl);

      write,format=
        "Star#%2d   % 5.2f % 6.1f % 6.1f     %6.1f   %.3f     %6.1f        %5.1f\n",
        jt,(*target.lambda)(jl),(*target.xposition)(jt),(*target.yposition)(jt),
        fwhm(jt,jl),strehl(jt,jl),e50(jt,jl),nmodes(jt);
    }
    write,format="Field Avg % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(avg,jl),strehl(avg,jl),e50(avg,jl);
    write,format="Field rms                         %6.1f   %.3f     %6.1f\n",
      fwhm(rms,jl),strehl(rms,jl),e50(rms,jl);
  }
  
  // Some logging of the results in file parprefix+".res":
  f	= open(parprefix+".res","a+");
  write,f,format=
    "\n         lambda   XPos   YPos  FWHM[mas]  Strehl  E50d[mas]  #modes comp.%s\n","";
  for (jl=1;jl<=target._nlambda;jl++) {
    for (jt=1;jt<=target._ntarget;jt++) {
      write,f,format=
        "Star#%2d   % 6.2f % 6.1f % 5.1f     %6.1f   %.3f     %6.1f        %5.1f\n",
        jt,(*target.lambda)(jl),(*target.xposition)(jt),(*target.yposition)(jt),
        fwhm(jt,jl),strehl(jt,jl),e50(jt,jl),nmodes(jt);
    }
    write,f,format="Field Avg % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(avg,jl),strehl(avg,jl),e50(avg,jl);
    write,f,format="Field rms                         %6.1f   %.3f     %6.1f\n",
      fwhm(rms,jl),strehl(rms,jl),e50(rms,jl);
  }
  write,f,format="\nAverage images written in %s\n",parprefix+"-imav.fits";
  write,f,format="Some other graphics in %s\n",parprefix+".ps";
  write,f,format="\nOriginal parameter file: %s:\n",oparfile;
  fin = open(oparfile,"r");
  while (1) { l = rdline(fin); if (!l) {break;}; write,f,l;}
  close,fin;
  write,f,"\n==== dump of structures ====";
  write,f,print(sim);
  write,f,print(atm);
  write,f,print(wfs);
  write,f,print(dm);
  write,f,print(mat);
  write,f,print(tel);
  write,f,print(target);
  write,f,print(gs);
  write,f,print(loop);
  close,f;
  
  fitsWrite,parprefix+"-imav.fits",imav;

  // saved graphics
  window,7,display="",hcp=parprefix+".ps",wait=1,style="work.gs";
  fma;

  disp2D,im,*target.xposition,*target.yposition,1,zoom=*target.dispzoom,init=1;
  for (jl=1;jl<=target._nlambda;jl++) {
    disp2D,imav(,,,jl),*target.xposition,*target.yposition,1,power=0.5;
  }

  for (j=1;j<=nwfs;j++) {
    plg,wfs(j).gspos(2),wfs(j).gspos(1),marker='\2',
      type="none",marks=1,color="red";
  }
  axisLegend,"arcsec","arcsec";
  pltitle=parprefix+"/ Average PSF";
  hcp;
  if (strehlsp != []) {
    fma;
    limits;
    limits,square=0;
    plg,strehlsp,itv;
    if (strehllp!=[]) plg,strehllp,itv,color="red";
    myxytitles,"Iterations","Strehl";
    hcp;
  }
  plt,sim.name,0.01,0.01,tosys=0;
  hcp_finish;
  
  if (curw != -1) {window,curw;}
  //  if (is_set(disp)) {window,style="boxed.gs";}
  return imav;
}

//----------------------------------------------------

/* DOCUMENT wfs, dm, atmospheric, etc.. structures:
   Main structures for the AO simul package parameters
   If additional parameters are needed, they should be entered in these
   structures definition and changes reflected in the parameter file
   (e.g. sh12.par)

   There are several type of entries:
   - long, float, string scalars: e.g.
   > atm.dr0at05mic = 33.;
   - pointers: These are pointers to variables that can have arbitrary
   number of elements. You have to define them in the following way:
   > wfs.nsubperring = &( [6,12,18] );

   Structure members can be accessed in the following way:
   var = dm.type;
   var = *atm.layerfrac;  // "*" dereference a pointer

   If there are several instance of a given structure: for instance, it is
   common to have a system with several dm:
   dm.type is thus a vector of the types of all instance of structure "dm"
   dm(1).type is a scalar
   *wfs(1).nsubperring is a vector.
   (*wfs(1).nsubperring)(1) is the first element of this vector.
   
   The variables with a "_" in frnt of them are internal variables. they are set
   and used while the program is running. You can still access them, and possibly
   modify them to reach a particular purpose.

   The following generic structures are instanced into structures of the same
   name, without the "_struct" appended, when the parameter file is read. For
   instance, "atm" will be the structure containing the atmospheric parameters
   as defined in the generic type atm_struct below. 


   SYNTAX OF THE COMMENTS BELOW:
   For each entries, we give the type (scalar, vectorptr -vector pointer-, etc),
   what the parameter is, possible comments, whether the parameter is
   required or optional, and the default between [].
   As a general comment: when the structures are instanciated, all their elements
   get a default value. This is zero (0) for a float or long (scalar or vector),
   empty string for a string, and 0x0 for a pointer.

   modified 2004oct18 for aosimulv3.3 to 3.5
   modified 2004july for aosimulv3.0 to 3.3
   modified 2004jan-mar for aosimulv2.4 to 3.0
   modified 2003dec20-24 for aosimul-v2.3
   modified 2003dec19 for aosimul-v2.2
   modified 2003feb19 for aosimul-v1.2
   AUTHOR: F.Rigaut, beginning 2002
   SEE ALSO: All ao simul functions (aoread, aoinit, aoloop).
*/

struct sim_struct
{
  string  name;           // A name for this simulation run. Optional [none]
  long    pupildiam;      // Pupil diameter in pixels. Required [none]
  long    debug;          // set the debug level (0:no debug, 1:some, 2: most). Optional [0]
  long    verbose;        // set the verbose level (0:none, 1: some, 2: most). Optional [0]
  // Internal keywords:
  long    _size;          // Internal. Size of the arrays [pixels]
  float   _cent;          // Internal. Pupil is centered on (_cent,_cent)
};

struct atm_struct
{
  float   dr0at05mic;     // Dr0 at sensing wavelength, at zenith. Required [none].
  pointer screen;         // string vectorptr. Screen file names. Required [none].
  pointer layerfrac;      // float vectorptr. Layer fraction. Sum to one is insured
  // in aoinit. Required [none]
  pointer layerspeed;     // float vectorptr. Layer speed. Required [none]
  pointer layeralt;       // float vectorptr. Layer altitude (m). Specified at Zenith.
  // Required [none]
  pointer winddir;        // Wind dir (use 0 for now)
  // Internal variables
  pointer _layeralt;      // float vectorptr. Actual layer altitude (m), from atm.alt & zen.angle
};

struct opt_struct
{
  string   phasemaps;      // filename of phasemap. Z scale should be in microns
  float    alt;            // float. equivalent altitude in m.
  float    misreg(2);      // float vector. misreg. (similar to DM, see below)
  float    _cent;          // center of the phase maps arrays (similar to sim._cent)
};

struct wfs_struct
{
  string  type;           // WFS type: "curvature", "hartmann" or "pyramid". Required [none]
  long    subsystem;      // Subsystem this WFS belongs to. Optional [1]
  float   lambda;         // WFS wavelength in microns. Required [none]
  long    noise;          // Enable noise (photon noise/read out noise). Optional [0=no].
  float   ron;            // Read out noise in electrons. Optional [0]
  float   darkcurrent;    // dark current in e-/sec/pixel or APD. Optional [0]
  float   gspos(2);       // This WFS guide star position (x<y) in arcsec. Optional [0,0]
  float   gsalt;          // This WFS guide star altitude in meter. 0 for infinity.
  // Specified at zenith. Optional [0]
  float   gsdepth;        // This WFS GS depth in meter (e.g. Na layer thickness).
  // Specified at zenith. Optional [0]
  float   laserpower;     // this wfs laser power (Na laser only), in Watts projected on sky.
  // Required when using lasers. Exclusive with gsmag; i.e. define one
  // OR the other. gs.
  float   gsmag;          // This WFS guide star magnitude. Optional [0]. For LGSs, see above.
  float   skymag;         // This WFS sky mag. Optional [no sky]
  long    filtertilt;     // Filter TT on this sensor? Optional [0=no]
  long    correctUpTT;    // Correct up link tp-tilt ? Optional [0=no]
  float   uplinkgain;     // Up link TT loop gain. Optional [0]
  float   dispzoom;       // Zoom factor for the display (typically around 1). Optional [1]
  float   optthroughput;  // optical throughput to WFS. Optional [1.0]
  // Curvature WFS only keywords:
  pointer nsubperring;    // Long vectorptr. # subapertures per ring. Required [none]
  pointer angleoffset;    // float vectorptr. offset angle for first subaperture of ring.
  float   l;              // Extra focal distance in a F/60 beam. Required [none]
  pointer rint;           // float vectorptr. if set, specify the inner radius for each ring
  pointer rout;           // float vectorptr. if set, specify the outer radius for each ring
  float   fieldstopdiam;  // diameter of field stop in arcsec. Optinal [1]. used only
  // to compute sky contribution (with skymag).
  // Shack-Hartmann WFS only keywords:
  long    shmethod;       // 1 = simple gradient average, 2=full propagation. Required [none]
  long    shnxsub;        // # of subaperture in telescope diameter. Required [none]
  float   pixsize;        // Subaperture pixel size in arcsec. Required [none]
  long    npixels;        // Final # of pixels in subaperture. Required [none]
  float   shthreshold;    // Threshold in computation of subaperture signal, >=0. Optional [0]
  float   biasrmserror;   // rms error on WFS bias in electron. Optional [0]
  float   flatrmserror;   // rms error on WFS flat, referenced to 1. Optional [0]
  string  submask;        // fits file with subaperture amplitude mask. It should have
                          // dimension 2^sdimpow2 square. can be float or long.
  // Typical value can be 0.01
  float   kernel;         // FWHM in arcsec of WFS gaussian kernel. Optional.
  // Default is computed as a function of D/r0
  int     nintegcycles;   // # of cycles/iterations over which to integrate. Optional [1]
  float   fracIllum;      // fraction illuminated to be valid. Optional [0.5]
  float   LLTxy(2);       // 2 element vector with (x,y) of laser projector [m]
  long    centGainOpt;    // Centroid Gain optimization flag. only for LGS (correctupTT and 
  // filtertilt must be set). Optional [0]
  int     rayleighflag;   // set to one to take rayleigh into account
  // Internal keywords:
  float   _gsalt;         // This WFS guide star altitude in meter. 0 for infinity.
  float   _gsdepth;       // This WFS GS depth in meter (e.g. Na layer thickness).
  int     _nsub;          // Internal. Tot # of subs.
  long    _nmes;          // internal. Tot # of measurements.
  pointer _sind;          // Internal: see CurvWFS
  pointer _nsind;         // Internal: see CurvWFS
  pointer _cxdef;         // Internal: see CurvWFS
  pointer _sxdef;         // Internal: see CurvWFS
  pointer _tiltsh;        // Internal: see ShWfs
  pointer _masks;         // Internal: see ShWfs
  pointer _fluxpersub;    // Internal: see ShWfs
  pointer _raylfluxpersub;// Internal: see ShWfs
  pointer _skyfluxpersub; // Internal: see ShWfs
  float   _nphotons;      // Internal: see WFS routines
  float   _skynphotons;   // Internal: see WFS routines
  float   _tt(2);         // Internal: WFS measured Tip and tilt
  float   _lastvalidtt(2);// Internal: WFS measured Tip and tilt
  float   _upttcommand(2);// Internal:
  pointer _refmes;        // internal: reference measurement vector
  pointer _tiprefv;       // internal: tip reference meas. vector
  pointer _tiltrefv;      // internal: tilt reference meas. vector
  pointer _tiprefvn;      // internal: tip reference meas. vector. normalized (norm=1)
  pointer _tiltrefvn;     // internal: tilt reference meas. vector. normalized.
  pointer _istart;        //
  pointer _jstart;        //
  pointer _binindices;    //
  int     _binxy;         //
  pointer _centroidw;     //
  pointer _fimage;        //
  pointer _fimage2;       //
  pointer _imistart;      //
  pointer _imjstart;      //
  int     _fimnx;         //
  int     _fimny;         //
  pointer _bias;          // internal: array of bias error
  pointer _flat;          // internal: array of flat error
  long    _domask;        // internal. flag to do submask calculations in _shwfs
  pointer _submask;       // internal: array. subaperture amplitude mask.
  pointer _kernel;        // internal: kernel for _shwfs. use: dointer or LGS uplink im.
  pointer _kernels;       // internal: subaperture dependant image kernel
  pointer _kerfftr;       // internal: storage of FFTs of kernels
  pointer _kerffti;       // internal: storage of FFTs of kernels
  int     _kernelconv;    // interal: convolve with kernel in _shwfs?
  int     _cyclecounter;  // counter in integration sequence (see nintegcycles above)
  pointer _dispimage;     // image to display (same as fimage except if nintegcycles!=1)
  pointer _x;             // shwfs: X positions of subaperture centers
  pointer _y;             // shwfs: Y positions of subaperture centers
  float   _centroidgain;  // internal: centroid gain if dithering on
  pointer _rayleigh;      // pointer to rayleigh images array for this sensor
  pointer _bckgrdcalib;   // pointer to background array calibration
  int     _bckgrdinit;    // set to one to fill calibration array
  int     _bckgrdsub;     // set to one to subtract background (default)
};

struct dm_struct
{
  string  type;           // "bimorph", "stackarray" "tiptilt", "zernike", "aniso". Required [none]
  long    subsystem;      // Subsystem this DM belongs to. Optional [1]
  string  iffile;         // Influence function file name. Leave it alone.
  float   alt;            // Conjugation altitude in meter. Specified @ zenith! Optional [0]
  float   hyst;           // DM actuator hysteresis (0. to 1.). Optional [0]
  float   push4imat;      // Voltage to apply for imat calibration. Optional [20].
  // Note: the default is not OK for many configs. Change at will.
  float   thresholdresp;  // Normalized response threshold for an act. to be valid. Optional [0.3]
  float   unitpervolt;    // Influence function sensitivity in unit/volt. Optional [0.01]
  // Stackarray: mic/volt, Tip-tilt: arcsec/volt.
  float   maxvolt;        // Saturation voltage (- and +) in volt. Optional [none if not set]
  float   gain;           // loop gain for this DM (total = this x loop.gain). Optional [1]
  float   misreg(2);      // dm misregistration (pixels). optional [0,0]
  // Bimorph-only keywords:
  pointer nelperring;     // long vectorptr. # of elec. per ring, e.g &([6,12,18]). Required [none]
  pointer angleoffset;    // float vectorptr. offset angle for first electrode in ring.
  pointer rint;           // float vectorptr. if set, specify the inner radius for each ring
  pointer rout;           // float vectorptr. if set, specify the outer radius for each ring
  float   supportRadius;  // radius of support points, normalized in pupil radius
  float   supportOffset;  // angle offset of first support points, in degree (default=90)
  // Stackarray-only keywords:
  long    nxact;          // # of act. in pupil diameter. Required [none]
  long    pitch;          // Actuator pitch IN PIXEL. Required [none]
  float   pitchMargin;    // margin to include more corner actuators when creating inf.functions
  // optional [1.44]
  float   coupling;       // coupling coef in influence function. optional [0.2].
  // valid values from 0.04 to 0.30.
  string  ecmatfile;      // valid to extrap. actuator matrix (extrap_com). Optional.
  long    noextrap;       // set to 1 if no extrapolated actuators whatsoever are to be used [0]
  long    elt;            // set to 1 if fast dmsum to be used
  // Zernike-only keywords:
  long    nzer;           // Number of modes, including piston. Required [none]
  // Internal keywords:
  float   _alt;           // Actual conjugation altitude in meter, from dm.alt and zen.
  long    _nact;          // Internal. Tot # of actuators.
  pointer _def;           // Internal: Pointer to IF data
  pointer _x;             // Internal: x position of actuator in pixels
  pointer _y;             // Internal: x position of actuator in pixels
  pointer _i1;            // 
  pointer _j1;            // 
  pointer _ei1;           // 
  pointer _ej1;           // 
  string  _eiffile;       // Influence function file name for extrap. actuators
  pointer _edef;          // Internal: Pointer to IF data for extrap. actuators
  pointer _ex;            // Internal: x position of extrap. actuator in pixels
  pointer _ey;            // Internal: x position of extrap. actuator in pixels
  long    _enact;         // Internal. Tot # of extrap. actuators.
  long    _n1;            // Internal: position of leftmost pixel in ao._size^2 array
  long    _n2;            // Internal: position of leftmost pixel in ao._size^2 array
  pointer _vold;          // internal: hysteresis
  pointer _posold;        // internal: hysteresis
  pointer _chpos;         // internal: hysteresis
  pointer _chv;           // internal: hysteresis
  pointer _dir;           // internal: hysteresis
  pointer _delta;         // internal: hysteresis
  pointer _command;       // pointer to command vector
  pointer _extrapcmat;    // extrapolation matrix: extrap_com = extrapmat(,+)*valid_com(+)
  int     _eltdefsize;    // size of def in case elt=1
};

struct mat_struct
{
  pointer condition;      // float vecorptr. Condition numbers for SVD, per subsystem. Required [none]
  string  file;           // iMat and cMat filename. Leave it alone.
};

struct tel_struct
{
  float   diam;           // Telescope diameter in meters. Required [none]
  float   cobs;           // Central obstruction / telescope diameter ratio. Optional [0]
};

struct target_struct
{
  pointer lambda;         // float vectorptr. Image wavelengths in micron. Required [none]
  pointer xposition;      // float vectorptr. X positions in arcsec. Required [none]
  pointer yposition;      // float vectorptr. Y positions in arcsec. Required [none]
  pointer dispzoom;       // float vectorptr. Display zoom (typically around 1.). Optional [1.]
  // Internal keywords
  long    _ntarget;       // Internal: # of target
  long    _nlambda;       // Internal: # of lambda
};

struct gs_struct
{
  float   zeropoint;      // Photometric zero point (#photons@pupil/s/full_aper, mag0 star).
  // Required [none]
  float   zenithangle;    // zenith angle. Optional [0.]. The zenith angle is used to compute:
  // - r0 off-zenith
  // - atmopheric turbulence layer altitude
  // - LGS altitude and thickness of Na Layer
  // - LGS brighness
  // note that dm altitude is unchanged.
  float  lgsreturnperwatt;// Sodium LGS return in photons/cm2/s at entrance pupil.
  // Specified at zenith. Modified by gs.zenithangle. Optional [22.]
  // basically, you have to fold in this the sodium density
  // and your model of return.
};

struct loop_struct
{
  float   gain;            // Loop gain. Optional, but important! [0]
  long    framedelay;      // # of frame delay ON TOP OF normal 1 frame delay. Optional [0]
  long    niter;           // # of total iteration. Required [none]
  float   ittime;          // Iteration time in seconds. Required [none]
  long    startskip;       // # iter to skip before collecting statistics. Optional [10]
  long    skipevery;       // skip by "skipby" every "skipevery" iterations. Optional [0=none]
  long    skipby;          // see above. this is to get better statistical
  long    stats_every;     // compute stats every so many iteration (default 4)
  // coverage. Optional [10000]
  long    jumps2swapscreen;//number of jumps (i.e. niter/skipevery) after which screens
  //will be swapped (rotation, 2->1, 3->2... 1->last
  string  modalgainfile;   // Name of file with mode gains. Optional.
  //  float   dithering;      // TT dithering for centroid gain (volts).
};

dtor = pi/180.;
radeg = 1./dtor;

