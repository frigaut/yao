/*
  YAO WFS functions
*/
  
func shwfs_init(pupsh,ns,silent=,imat=,clean=)
/* DOCUMENT func shwfs_init(pupsh,ns,silent=,imat=)
   pupsh: pupil image. Should be 1/0 for SHWFS. dimension sim._size
   ns: WFS # (for multi-wfs systems, if not, put 1)
   silent: be silent (except for error messages of course)
   imat: will be used for imat computation (will use kernel)
   clean: force recomputing everything (in particular rayleight maps)
   SEE ALSO: sh_wfs
 */

{
  if (silent==[]) silent = (sim.verbose==0);

  wfs._initkernels = array(1n,nwfs);

  if (is_void(ns)) {ns=1;} // default to wfs#1 for one WFS work.

  if (typeof(pupsh) != "float") {error,"pupsh was not float !";}

  pupd       = sim.pupildiam;
  size       = sim._size;
  nxsub      = wfs(ns).shnxsub(0);
  subsize    = pupd/nxsub;
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;
  fracsub    = wfs(ns).fracIllum;
  sdim       = long(2^ceil(log(subsize)/log(2)+1));
  sdimpow2   = int(log(sdim)/log(2));

  wfs(ns)._centroidgain = 1.f;

  //  if (anyof(wfs.svipc>1)) status = quit_wfs_forks();

  //====================================================================
  // WORK OUT THE NUMBER OF PHOTONS COLLECTED PER SUBAPERTURE AND SAMPLE
  //====================================================================

  telSurf  = pi/4.*tel.diam^2.;
  
  // from the guide star (computed here as used in wfs_check_pixel_size):
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
  if (wfs(ns).shmethod ==1 ) wfs(ns)._kernelconv = 0n;

  quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;

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

  is  = size/2+1-pupd/2;
  if (wfs(ns).npixpersub) {
    is  = size/2+1-(subsize*nxsub)/2;
  }

  if (wfs(ns).pupoffset!=[]) {
    
    puppixoffset = long(round(wfs(ns).pupoffset/tel.diam*sim.pupildiam));

  } else puppixoffset = [0,0];

  nsubs = nxsub*nxsub;
//  istart = (indgen(nsubs)-1)/nxsub*subsize  +is+puppixoffset(1);
//  jstart = ((indgen(nsubs)-1)%nxsub)*subsize+is+puppixoffset(2);

  istart = ((indgen(nsubs)-1)%nxsub)*subsize + is + puppixoffset(1);
  jstart = ((indgen(nsubs)-1)/nxsub)*subsize + is + puppixoffset(2);
  xsub = (((indgen(nsubs)-1)%nxsub)+0.5)*tel.diam/nxsub-tel.diam/2.;
  ysub = ((indgen(nsubs)-1)/nxsub+0.5)*tel.diam/nxsub-tel.diam/2.;

  //==========================================================
  // COMPUTE WHICH SUBAPERTURES ARE ENABLED (FLUX > THRESHOLD)
  //==========================================================

  fluxPerSub = array(float,nsubs);
  
  for (i=1;i<=nsubs;i++) {

    fluxPerSub(i) = sum(pupsh(istart(i):istart(i)+subsize-1,
                              jstart(i):jstart(i)+subsize-1));
  }
  
  fluxPerSub = fluxPerSub/subsize^2.;

  // indices of the enabled subapertures: gind
  //  gind     = where(fluxPerSub > fracsub);
  // changed 2009oct07: we display *all* subaps for which there is flux
  // this is to have a more realistic display (wfs._fimage)
  // however, now we have 2 sets of valid subapertures:
  // - the one for which we want to compute an image in _shwfs() -> _nsub4disp
  // - the really valid ones, for which we will compute the slope information -> _nsub

  if (wfs(ns).shmethod==2) {

    gind       = where(fluxPerSub > 0);

  } else { 

    // shmethod=1 -> geometrical SH. 
    // for these, it makes no sense to differentiate display and slope/valid
    // subapertures. We'll make them the same
    gind       = where(fluxPerSub > fracsub);

  }

  // then out of these, we will only compute mesvec for the "valid":
  tmp = fluxPerSub;
  tmp = tmp(gind);
  wfs(ns)._validsubs = &(int(tmp > fracsub));
  
  istart     = istart(gind);
  jstart     = jstart(gind);
  xsub       = xsub(gind);
  ysub       = ysub(gind);
  fluxPerSub = fluxPerSub(gind);

  // stuff some of wfs structure for WFS "ns":
  wfs(ns)._istart = &(int(istart-1)); // -1n 'cause C is 0 based
  wfs(ns)._jstart = &(int(jstart-1));
  wfs(ns)._x = &(xsub);
  wfs(ns)._y = &(ysub);
  wfs(ns)._nsub4disp = int(numberof(gind));
  wfs(ns)._nsub = int(sum(*wfs(ns)._validsubs));

  // compute other setup variables for _shwfs:

  //==========================================================================
  // COMPUTE WFS IMAGE KERNELS: SUBAPERTURE DEPENDANT. USED FOR LGS ELONGATION
  //==========================================================================

  if (wfs(ns)._kernelconv != 0n) {
    // if kernelconv is 0, then the _shwfs routine does not use wfs._kernels,
    // so no need to compute it.

    quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
    xy = (indices(sdim)-sdim/2.-1)*quantumPixelSize;  // coordinate array in arcsec

    if (sim.verbose >= 1) {
      write,format="%s\n","Pre-computing Kernels for the SH WFS";
    }

    kall = [];

    gui_progressbar_frac,0.;
    gui_progressbar_text,swrite(format=\
      "Precomputing subaperture kernels, WFS#%d",ns);
    
    for (l=1; l<=wfs(ns)._nsub4disp; l++) {
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
      // below: this was too many characters to go thru the python pipe.
      // with fast processor, we compute fast and send too fast. 2009oct22.
      //gui_progressbar_text,swrite(format=\
      //  "Precomputing subaperture kernels for LGS elongation, WFS#%d, sub#%d/%d",\
      //  ns,l,wfs(ns)._nsub4disp);
      gui_progressbar_frac,float(l)/wfs(ns)._nsub4disp;

    }

    clean_progressbar;
    //if (sim.debug >=2) {hitReturn;}

    wfs(ns)._kernels = &(float(kall));
    wfs(ns)._kerfftr = &(float(kall*0.));
    wfs(ns)._kerffti = &(float(kall*0.));
    kall = [];

  }
  
  //==========================================================================
  // COMPUTE WFS IMAGE KERNELS: SUBAPERTURE DEPENDANT. USED FOR LGS ELONGATION
  //==========================================================================

  rayleighflux = array(0.0f,wfs(ns)._nsub4disp);
  sodiumflux   = array(1.0f,wfs(ns)._nsub4disp);

  if (wfs(ns).rayleighflag == 1n) {
    
    quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
    xy = (indices(sdim)-sdim/2.-1)*quantumPixelSize;  // coordinate array in arcsec

    if (sim.verbose > 0) {write,"Pre-computing Rayleigh for the SH WFS";}

    kall = [];

    rayfname = parprefix+"-rayleigh-wfs"+swrite(format="%d",ns)+"-zen"+
      swrite(format="%d",long(gs.zenithangle))+".fits";
    isthere = fileExist(YAO_SAVEPATH+rayfname);
    fov    = quantumPixelSize*sdim;
    aspp   = quantumPixelSize;


    if ( (isthere) && (!clean) ) {
      
      kall         = fitsRead(YAO_SAVEPATH+rayfname)*1e6;
      rayleighflux = fitsRead(YAO_SAVEPATH+rayfname,hdu=1);
      sodiumflux   = fitsRead(YAO_SAVEPATH+rayfname,hdu=2);

    } else {

      for (l=1; l<=wfs(ns)._nsub4disp; l++) {
        
        xsub = (*wfs(ns)._x)(l); ysub = (*wfs(ns)._y)(l);
        tmp = mcao_rayleigh(ns,ysub,xsub,fov=fov,aspp=aspp,zenith=gs.zenithangle);
        rayleighflux(l) = sum(tmp(,,1));
        sodiumflux(l)   = sum(tmp(,,2));
        tmp = transpose(tmp(,,1));
        // the switch of xsub <-> ysub and transpose are to accomodate the
        // C vs yorick 0 vs 1 start array index.
        //        tmp = tmp/sum(tmp);
        grow,kall,(eclat(tmp))(*);

      }

      fitsWrite,YAO_SAVEPATH+rayfname,kall;
      fitsWrite,YAO_SAVEPATH+rayfname,rayleighflux,append=1,exttype="image";
      fitsWrite,YAO_SAVEPATH+rayfname,sodiumflux,append=1,exttype="image";

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
  } else wfs(ns)._rayleigh = &([0.0f]); // to avoid type conversion error in _shwfs call
    

  //================================
  // SUBAPERTURE SIZE AND PIXEL SIZE
  //================================
  wfs_check_pixel_size,ns,sdim,rebinFactor,actualPixelSize,\
    printheader=(ns==1),silent=silent;

  // now compute _shwfs C routine internal array size:
  // for bimage (trimmed and rebinned simage):
  wfs(ns)._rebinfactor = rebinFactor;
  nbigpixels = long(sdim/rebinFactor);
  // check eveness of nbigpixels is same as wfs.npixels:
  if (even(wfs(ns).npixels)!=even(nbigpixels)) nbigpixels--;
  
  // subsize of sdim that we will use in bimage:
  rdim = long(rebinFactor*nbigpixels);
  xy = long(indices(rdim)-1.);
  tmp = xy(,,1)/rebinFactor + xy(,,2)/rebinFactor*nbigpixels;
  
  // what's the rebinned pixels size (in small fft pixels)?
  // answer -> rebinFactor
  // how many integer big pixels can we fit in sdim (which ought to be
  // a power of 2)?
  // answer = long ( sdim/rebinFactor )
  // the number of rebinned pixels in the final subaperture is wfs.npixels
  // (1) if wfs.npixels is even, the spot has to be in the center of the subap 
  //   image
  // (2) if wfs.npixel is odd, the spot has to be in the center of the center
  //   rebinned pixel.
  // remember the spot is shifted by one half fft pixels left and down,
  // so it can fall in between original pixels hence in between rebinned
  // pixels. thus by playing here we can only accomodate:
  // (1) wfs.npixel even, both rebinfactor odd and even
  // (2) wfs.npixel odd, only rebinfactor even.
  // (3) with wfs.npixel odd and rebinfactor odd, the spot would have to
  // be centered on the central rebinned pixel, which has and odd number
  // of fft pixels, thus would have to be centered on a fft pixel, which is
  // not originally the case. will have to modify tiltsh for that.
  
  // for now let's take care of how to dimension binindices for all those
  // cases. after the fft, and 1/2 pixel shift, the spot is centered on 
  // the 2^n x 2^n fft array.
  
  
  binindices = array(-1l,[2,sdim,sdim]);
  binindices(1:rdim,1:rdim) = tmp;
  ss = ceil((sdim-rdim)/2.);
  binindices = roll(binindices,[ss,ss]);
  binindices = int(eclat(binindices));
  // stuff some more of wfs structure for WFS "ns":
  wfs(ns)._binindices = &(int(binindices));
  // checked, it seems to work.

  wfs(ns)._binxy = nbigpixels;

  // centroid reference vector, after final extraction of subimage:
  centroidw = indgen(wfs(ns).npixels)-1.-(wfs(ns).npixels/2.-0.5);
  // we might as well express it in arcsec:
  centroidw = float(centroidw*actualPixelSize);
  wfs(ns)._centroidw = &centroidw;

  // 2004mar22: added a guard pixel for each subaperture for the display
  // 2009oct06: removed it, in the process of implementing the optical 
  // coupling between subapertures
  wfsnpix = wfs(ns).npixels;

  imistart = (istart-min(istart))/subsize*(wfsnpix);
  imjstart = (jstart-min(jstart))/subsize*(wfsnpix);  
  wfs(ns)._imistart = &(int(imistart));
  wfs(ns)._imjstart = &(int(imjstart));

  wfs(ns)._imistart2 = &(int(imistart+(nbigpixels-wfsnpix)/2));
  wfs(ns)._imjstart2 = &(int(imjstart+(nbigpixels-wfsnpix)/2));

  fimdim = long(nxsub*wfsnpix+(nbigpixels-wfsnpix));
  wfs(ns)._fimage = &(array(float,[2,fimdim,fimdim]));
  wfs(ns)._fimnx = int(fimdim);
  wfs(ns)._fimny = int(fimdim);

  // This is the tilt to add to the input phase so that
  // the individual subaperture images fall in between
  // the pixels of the quadcell
  xy     = indices(sim._size);

  wfs(ns)._tiltsh = &(float((-64.)*0.098174773*(xy(,,1)+xy(,,2))* \
               0.5/sdim*wfs(ns).lambda/(2*pi)*(wfs(ns).shmethod == 2)));

  // This tilt array is intended to bring the spot back inbetween 4 pixels
  // instead of centered on dim/2+1 as a result of the regular FFT.
  // In other words, added to the subaperture phase, it will shift the
  // image 1/2 FFT pixel left and down.
  // this is an achromatic factor of course that just depends on the
  // dimension of the array.
  // the lambda/2pi factor is thus to compensate the x by 2pi/lambda
  // in _shwfs, making tiltsh achromatic.
  // If the number of (rebinned) pixels (wfs.npixels) is odd,
  // and if the rebinFactor (number of small FFT pixels in a rebinned one)
  // is also odd, then the spot should be centered on a small FFT pixel,
  // not in between 4 of them. 
  if (odd(rebinFactor) && odd(wfs(ns).npixels) && (wfs(ns).shmethod==2)) {
    *wfs(ns)._tiltsh *= 0;
  }


  //================================
  // FIELD STOP / AMPLITUDE MASK
  //================================

  if (sim.debug>1) \
    write,format="Dimension for optional amplitude mask: %d\n",2^sdimpow2;  

  // reads out the amplitude mask for the subaperture:
  if (wfs(ns).fsname) {

    // read the amplitude image
    tmp = fitsRead(YAO_SAVEPATH+wfs(ns).fsname);

    // check that dims are OK:
    if (anyof(dimsof(tmp)!=[2,2^sdimpow2,2^sdimpow2])) {
      error,swrite(format="Bad dimensions for %s. Should be %d, found %d\n",
                   wfs(ns).fsname,2^sdimpow2,dimsof(tmp)(2));
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

  } else if (strlen(wfs(ns).fstop)>0) {

    // make the field stop with wfs.fstop, wfs.fssize and wfs.fsoffset
    make_fieldstop,ns;
    wfs(ns)._domask = 1l;

  } else wfs(ns)._domask = 0l;


  // compute # of photons from the sky as above for guide star.
  // skymag is per arcsec, so we have to convert for
  // the subaperture size computed in wfs_check_pixel_size.;
  // changed on 2009oct12: now _skynphotons is normalized
  // *per rebinned pixel*.
  subsize  = sim.pupildiam/wfs(ns).shnxsub(0);
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;
  // subsize in meters:
  subsize_m = subsize * tel.diam/sim.pupildiam;


  wfs(ns)._skynphotons = gs.zeropoint*10^(-0.4*wfs(ns).skymag)* // #photons/tel/sec/arcsec^2
    subsize_m^2./telSurf*                    // -> per full subaperture
    loop.ittime*                             // -> per iteration
    wfs(ns).optthroughput*                   // -> include throughput to WFS
    wfs(ns).pixsize^2.;                      // -> per rebinned pixel

  if ( wfs(ns).skymag == 0) { wfs(ns)._skynphotons = 0.; } // if skymag not set

  if (sim.verbose == 2) {
    write,format="wfs(%d)._skynphotons = %f\n\n",ns,wfs(ns)._skynphotons;
  }
  
  // for guide star, total "useful" signal per subap.
  wfs(ns)._fluxpersub  = &(float(fluxPerSub*wfs(ns)._nphotons));
  
  if (sim.verbose > 0) {

    gstype = ( (wfs(ns).gsalt>0)?"LGS":"NGS" );
    tmp = fluxPerSub(where(*wfs(ns)._validsubs))*wfs(ns)._nphotons;
    if (min(tmp)>10) fmt="%.0f"; else fmt="%.1f";
    write,format="%s#%d flux varies between "+fmt+" and "+fmt+
      " photon/subap/it\n",gstype,ns,min(tmp),max(tmp);
  }
  
  // for rayleigh, if any:
  wfs(ns)._raylfluxpersub  = &(*wfs(ns)._fluxpersub*float(rayleighflux/sodiumflux));
    
  // for sky (see above, in photon/subap/it/rebinned pixel)
  wfs(ns)._skyfluxpersub  = &(float(fluxPerSub*wfs(ns)._skynphotons));

  // changed 2009oct13: now compute bias and flat for entire _fimage:
  wfs(ns)._bias = &(wfs(ns).biasrmserror *
                    gaussdev([2,wfs(ns)._fimnx,wfs(ns)._fimny]));

  wfs(ns)._flat = &(1.0f + wfs(ns).flatrmserror *
                    gaussdev([2,wfs(ns)._fimnx,wfs(ns)._fimny]));

  //same here, take background calib image for entire _fimage
  wfs(ns)._bckgrdcalib = &(array(float,[2,wfs(ns)._fimnx,wfs(ns)._fimny]));
  
  if (sim.verbose == 2) {
    write,format="Dark current wfs#%d / iter / pixel=%f\n",ns,
      float(wfs(ns).darkcurrent*loop.ittime);
  }

  // CALIBRATE BACKGROUND IMAGES (have to run sh_wfs for that):
  wfs(ns)._bckgrdsub  = 0; // it doesn't matter
  wfs(ns)._bckgrdinit = 1;

  // call sh_wfs for calibration of the background
  // first sync if needed (svipc)
  if (wfs(ns).svipc>1) status = sync_wfs_forks();
  
  sh_wfs,pupsh,pupsh*0.0f,ns;

  wfs(ns)._bckgrdinit = 0;
  wfs(ns)._bckgrdsub  = 1; // now yes, enable it (by default)
  
  // here we need to re-sync to restore bckgrd properties to forks:
  if (wfs(ns).svipc>1) status = sync_wfs_forks();
  
  if (show_background) {
    fma;
    plsys,1;
    pli,*wfs(ns)._bckgrdcalib;
    pltitle,swrite(format="Calibrated background for WFS#%d",ns);
    hitReturn;
  }

  // DISPLAY OF WFS CONFIG (put that somewhere else so it can be 
  // called independently)
  if ( (sim.debug>=1) && (!silent) && (wfs(ns).shmethod==2) ){
    sh_wfs,pupsh,pupsh*0.0f,ns;
    fma;
    plsys,2;
    pli,*wfs(ns)._fimage;
    limits;
    for (i=1;i<=wfs(ns)._nsub4disp;i++) {
      if ((*wfs(ns)._validsubs)(i)==0) continue;
      x1 = (*wfs(ns)._imistart2)(i);
      y1 = (*wfs(ns)._imjstart2)(i);
      l1 = wfs(ns).npixels;
      plg,_(y1,y1,y1+l1,y1+l1,y1),_(x1,x1+l1,x1+l1,x1,x1),color=[80,80,80],marks=0;
      plt,swrite(format="%d",i),x1,y1,tosys=1,color=[80,80,80],height=pltitle_height/2;
    }
    // display the full extent of this example subap FoV
    // to show the user how it overlaps with neightbors
    // first let's highlight the subaperture we're talking about here:
    i = where(*wfs(ns)._validsubs)(1);
    x1 = (*wfs(ns)._imistart2)(i);
    y1 = (*wfs(ns)._imjstart2)(i);
    l1 = wfs(ns).npixels;
    plg,_(y1,y1,y1+l1,y1+l1,y1),_(x1,x1+l1,x1+l1,x1,x1),color="green",width=3;
    // then fisplay the whole subaperture field of view
    x1 = (*wfs(ns)._imistart)(i);
    y1 = (*wfs(ns)._imjstart)(i);
    l1 = wfs(ns)._binxy; //*wfs(ns)._rebinfactor;
    plg,_(y1,y1,y1+l1,y1+l1,y1),_(x1,x1+l1,x1+l1,x1,x1),color="red";
    // field stop:
    if (*wfs(ns)._submask!=[]) {
      fs = roll(*wfs(ns)._submask);
      xyc = indices(dimsof(fs));
      xyc = (xyc-1.)/(dimsof(fs)(2)-1)*(wfs(ns)._binxy)+0.;
      xyc(,,1) += x1;
      xyc(,,2) += y1;
      plc,fs,xyc(,,2),xyc(,,1),levs=[0.001],color="magenta",width=3;

    }
    require,"plvp.i"; // for plmargin
    plmargin;
    xytitles_vp,"pixels","pixels",[0.015,0.015];
    pltitle_vp,escapechar(swrite(format="wfs(%d)._fimage",ns)),0.005;
    plsys,0;
    ybase = 0.72;
    deltay = 0.03;
    deltayt = 0.005;
    x1 = 0.16;
    x2 = x1 + 0.03;
    x3 = x2 + 0.01;
    plg,_(ybase,ybase),_(x1,x2),color=[80,80,80],width=3;
    plt,"Valid subapertures",x3,ybase-deltayt,tosys=0
    ybase -= deltay;
    plg,_(ybase,ybase),_(x1,x2),color="green",width=3;
    plt,"Highlighted subaperture",x3,ybase-deltayt,tosys=0
    ybase -= deltay;
    plg,_(ybase,ybase),_(x1,x2),color="red",width=3;
    plt,"Highlighted subaperture total FoV (overlap)",x3,ybase-deltayt,tosys=0
    ybase -= deltay;
    if (*wfs(ns)._submask!=[]) {
      plg,_(ybase,ybase),_(x1,x2),color="magenta",width=3;
      plt,"Highlighted subaperture field stop",x3,ybase-deltayt,tosys=0;
    } else {
      plt,"NO field stop defined",x3,ybase-deltayt,tosys=0
    }
    hitReturn;
    plsys,2;
    limits;
  }

  // and let's just re-sync for good measure:
  if (anyof(wfs.svipc>1)) status = sync_wfs_forks();
  
  return 1;
}


//----------------------------------------------------
func wfs_check_pixel_size(ns,&sdim,&rebinFactor,&actualPixelSize,printheader=,silent=)
/* DOCUMENT wfs_check_pixel_size(ns,&sdim,&rebinFactor,&actualPixelSize,printheader=,silent=)

   Finds the pixel size for the requested WFS configuration.

   There are constraints:
     - First, the pixel size can not be arbitrary.
       It is defined by lambda_wfs/pixel_size_in_pupil_space/sdim.
     - Second, the max subaperture size is lambda_wfs/pixel_size_in_pupil_space

   I could have preserved wfs.npixels as a strong constraint, but it lead to
   complicated algorithms. Instead, I am finding the closest pixel size to
   wfs.pixsize and then derive the number of pixel and subaperture size.

   In addition, I only allow an even number of pixels in the subaperture
   (otherwise I have to check many more things, as for instance an odd
   number of pixels in the subaperture combined with a odd rebinFactor
   cause problems. < now solved (2009oct07)
   This routine fills and returns binindices and centroidw (now located in 
   shwfs_init, 2009oct07)
   SEE ALSO: sh_wfs
 */
{
  extern wfs;
  
  
  pupd       = sim.pupildiam;
  nxsub      = wfs(ns).shnxsub(0);
  subsize    = pupd/nxsub;
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;
  // sdim is the dimension of "simage" in _shwfs(),
  // i.e. if 2^l is the smallest array size that can contains the subaperture
  // sdim is 2^(l+1) 
  sdim       = long(2^ceil(log(subsize)/log(2)+1));
  err        = 0;

  desiredPixelSize = wfs(ns)._origpixsize;
  desiredNpixels = wfs(ns).npixels;
  quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
  rebinFactor = long(round(desiredPixelSize/quantumPixelSize));
  if ( rebinFactor == 0 ) {rebinFactor = 1l;}
  actualPixelSize  = rebinFactor*quantumPixelSize;

  wfs(ns).pixsize = actualPixelSize;

  while ( (rebinFactor*wfs(ns).npixels) > sdim ) {
    wfs(ns).npixels -= 2;
    err = 1;
  }

  if (!is_set(silent)) {
    f = open(YAO_SAVEPATH+parprefix+".res","a+");
    if (is_set(printheader)) {
      write,"WFS# |       Pixel sizes         | Subap. size | Number of pixels | #photons";
      write,"     | Desired  Quantum  Actual  | Max  Actual | Desired   Actual | /sub/iter";
      write,f,"\nWFS# |       Pixel sizes         | Subap. size | Number of pixels | #photons";
      write,f,"     | Desired  Quantum  Actual  | Max  Actual | Desired   Actual | /sub/iter";
    }
    write,format="%2d      %.5f  %.5f  %.5f   %4.2f  %4.2f    %2dx%2d      %2dx%2d   %.1f\n",
      ns,desiredPixelSize,quantumPixelSize,actualPixelSize,quantumPixelSize*sdim,
      actualPixelSize*wfs(ns).npixels,desiredNpixels,desiredNpixels,
      wfs(ns).npixels,wfs(ns).npixels,wfs(ns)._nphotons;
    write,f,format="%2d      %.5f  %.5f  %.5f   %4.2f  %4.2f    %2dx%2d      %2dx%2d   %.1f\n",
      ns,desiredPixelSize,quantumPixelSize,actualPixelSize,quantumPixelSize*sdim,
      actualPixelSize*wfs(ns).npixels,desiredNpixels,desiredNpixels,
      wfs(ns).npixels,wfs(ns).npixels,wfs(ns)._nphotons;
    close,f;
  }
    
  //  if (err == 1) {
  //    write,format="  WARNING: #pixel/subaperture reduced to %d\n",wfs(ns).npixels;
  //  }
  
  if (wfs(ns).npixels == 0) {
    write,format="\nWFS#%2d: The desired pixel size is too large.\n",ns;
    write,"       I can not even fit 2x2 pixels.";
    write,"       Reduce pixel size in parfile or use a larger sim.pupildiam.\n";
    exit;
  }

  //  write,format="  Final config: %dx%d pixels, pixsize = %f\","+
  //    " updated in RAM but not in parfile.\n",
  //    wfs(ns).npixels,wfs(ns).npixels,wfs(ns).pixsize;


  return err;
}


//----------------------------------------------------

func sh_wfs(pupsh,phase,ns)
/* DOCUMENT func sh_wfs(pupsh,phase,ns)
   pupsh: pupil (dimension sim._size, 0 or 1). Must be of type "float"
   phase: phase. same dim as pupil. Must be of type "float"
   ns: WFS# (for multi WFS systems, default to 1)
   SEE ALSO: shwfs_init
 */

{
  //  extern mesvec;
  
  if (is_void(ns)) {ns=1;} // default to wfs#1 for one WFS work.

  // bail out if bad type (otherwise error in C function call)
  if (typeof(pupsh) != "float") {error,"pupsh was not float !";}
  if (typeof(phase) != "float") {error,"Phase was not float !";}

  // shorthand
  pupd       = sim.pupildiam;
  size       = sim._size;
  nxsub      = wfs(ns).shnxsub(0);
  subsize    = int(pupd/nxsub);
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;

  // The phase is in microns. this scaling factor restore it in radian
  // at the wfs lambda
  phasescale = float(2*pi/wfs(ns).lambda);   // wfs.lambda in microns

  // fast method (geometrical SHWFS)
  if (wfs(ns).shmethod == 1) {
    
    toarcsec = float(wfs(ns).lambda/2.0/pi/(tel.diam/sim.pupildiam)/4.848);

    // define mesvec (alloc space for C function)
    if (mesvec==[]) mesvec = array(float,2*wfs(ns)._nsub);

    err = _shwfs_simple(pupsh, phase, phasescale, *wfs(ns)._tiltsh, size, size,
                        *wfs(ns)._istart, *wfs(ns)._jstart, int(subsize),
                        int(subsize), wfs(ns)._nsub, toarcsec, mesvec);
                       
  // Full diffraction SHWFS
  } else {

    // size of array in which to embed phase image of 1 subaperture
    sdim       = long(2^ceil(log(subsize)/log(2)+1));
    // sdimpow2 such that sdim = 2^sdimpow2
    sdimpow2   = int(log(sdim)/log(2));
    threshold = array(float,wfs(ns)._nsub4disp+1)+wfs(ns).shthreshold;
    // "feature" noticed on 2010apr16 (thanks to Yann Clenet):
    // the threshold was correctly applied to the valid subap,
    // but not to the background surrounding the valid subap (fimage
    // unused pixels). This has no effect on the simulation, but
    // 1) may be surprising/misleading for the user
    // 2) thresholding these pixels may be useful for the user to
    //    decide on a threshold value.
    // hence, on v4.5.2, I have upgraded yao_fast.c to also threshold
    // these pixels. But as we are passing now a array of thresholds
    // (one per subap, in case in the future we want to upgrade the
    //  code to allow this as input), which value to chose? solution:
    // have a vector of threshold of dim # of subap + 1 (the last one
    // value is the one to apply to these outside pixels/area).
    // by default for now, we apply of course the same value as the
    // given scalar wfs.shthreshold. Hence the +1 in the formula above.

    // simple check to see if wfs(ns).svipc has been changed:
    //    if ((nforks_per_wfs!=[]) && (nforks_per_wfs(ns)!=wfs(ns).svipc)) \
    //      status = quit_wfs_forks();
    
    // init in case we use svipc:
    if (wfs(ns).svipc>1) {
      if (!wfs(ns)._svipc_init_done)  {
        require,"yao_svipc.i";
        status = svipc_wfs_init(phase,ns);
        extern shmphase,ffimage;
        wfs(ns)._svipc_subok = &((*wfs(ns)._fork_subs)(,1));
      }
      shm_var,shmkey,swrite(format="wfs%d_phase",ns),shmphase;
      shm_var,shmkey,swrite(format="wfs%d_mesvec",ns),mesvec;
      shm_var,shmkey,swrite(format="wfs%d_fimage",ns),ffimage;
      shmphase(,) = phase;
      ffimage(,)  = 0;
      mesvec()    = 0;

      yoffset = (*wfs(ns)._yoffset)(1);
      fimny   = (*wfs(ns)._fimny2)(1);
      subok2  = (*wfs(ns)._fork_subs2)(,1);
    
      // give trigger
      if (sim.debug>20) write,format="main: Giving trigger on sem %d\n",2*ns;
      sem_give,semkey,20+4*(ns-1),count=wfs(ns).svipc-1;

      
    } else {
      mesvec = array(float,2*wfs(ns)._nsub);
      wfs(ns)._svipc_subok = &array(1n,wfs(ns)._nsub4disp);
      subok2 = array(1n,wfs(ns)._nsub4disp);
      eq_nocopy,ffimage,*wfs(ns)._fimage;
      ffimage(,) = 0;
      yoffset    = 0n;
      fimny      = wfs(ns)._fimny;
    }
    
    // C function calls

    // phase to spot image (returned in ffimage)
    err = _shwfs_phase2spots( pupsh, phase, phasescale,
                 *wfs(ns)._tiltsh, int(size), *wfs(ns)._istart,
                 *wfs(ns)._jstart, int(subsize), int(subsize), 
                 wfs(ns)._nsub4disp, sdimpow2, wfs(ns)._domask, *wfs(ns)._submask, 
                 *wfs(ns)._kernel, *wfs(ns)._kernels, *wfs(ns)._kerfftr,
                 *wfs(ns)._kerffti, wfs(ns)._initkernels, wfs(ns)._kernelconv,
                 *wfs(ns)._binindices, wfs(ns)._binxy, 
                 wfs(ns)._rebinfactor, ffimage, *wfs(ns)._svipc_subok,
                 *wfs(ns)._imistart, *wfs(ns)._imjstart, 
                 wfs(ns)._fimnx , wfs(ns)._fimny, 
                 *wfs(ns)._fluxpersub, *wfs(ns)._raylfluxpersub,
                 *wfs(ns)._skyfluxpersub, float(wfs(ns).darkcurrent*loop.ittime),
                 int(wfs(ns).rayleighflag), 
                 *wfs(ns)._rayleigh, wfs(ns)._bckgrdinit,
                 wfs(ns)._cyclecounter, wfs(ns).nintegcycles);

    if ( wfs(ns).svipc>1 ) {
      if (sim.debug>20) write,format="main: waiting fork ready sem %d\n",2*ns+1;
      sem_take,semkey,20+4*(ns-1)+1,count=wfs(ns).svipc-1;
      sem_give,semkey,20+4*(ns-1)+2,count=wfs(ns).svipc-1;
    }    

    // spot image to slopes:
    
    err = _shwfs_spots2slopes( ffimage,
                  *wfs(ns)._imistart2, *wfs(ns)._imjstart2,
                  wfs(ns)._nsub4disp, wfs(ns).npixels,
                  wfs(ns)._fimnx , fimny, yoffset, 
                  *wfs(ns)._centroidw, threshold, *wfs(ns)._bias,
                  *wfs(ns)._flat, wfs(ns).ron, wfs(ns).noise, 
                  *wfs(ns)._bckgrdcalib, wfs(ns)._bckgrdinit, wfs(ns)._bckgrdsub,
                  *wfs(ns)._validsubs, subok2, wfs(ns).nintegcycles,
                  mesvec);
                  
    if ( wfs(ns).svipc>1 ) {
    sem_take,semkey,20+4*(ns-1)+3,count=wfs(ns).svipc-1;
    }
    
    // new, cause of svipc to keep *all* forks results in this var:              
    if (wfs(ns)._bckgrdinit) *wfs(ns)._bckgrdcalib = ffimage;
    
    //    if (wfs(ns)._refmes) write,mesvec-*wfs(ns)._refmes;
    //    else write,mesvec;
    
    wfs(ns)._fimage = &ffimage;
    
    wfs(ns)._initkernels = 0n;

    // FIXME: pass cyclecounter to child.
    wfs(ns)._cyclecounter += 1;
    if (wfs(ns)._cyclecounter > wfs(ns).nintegcycles) {wfs(ns)._cyclecounter = 1;}
  }
  
  if (err != 0) {error,"problem in _shwfs";}

  mesvec *= wfs(ns)._centroidgain;

  // return measurement vector in arcsec (thanks to centroiw):
  return mesvec;
}


//----------------------------------------------------

func curv_wfs(pupil,phase,ns,init=,disp=,silent=)
/* DOCUMENT curv_wfs(pupil,phase,ns,init=,disp=,silent=)
   This function computes the signal from a Curvature WFS for a 
   given phase and pupil input and WFS config (WFS #ns)
*/
{
  if (is_void(ns)) {ns=1;} // default sensor#1 for one WFS work
  
  size       = sim._size;
  dimpow2   = int(log(size)/log(2));

  if (init == 1) {
    if ( (sim.verbose>=1) && (!is_set(silent)) ) {write,"> Initializing curv_wfs\n";}
    fratio= 60.;
    defoc = (pi*wfs(ns).lambda*1e-6/(sim._size^2.*(tel.diam/sim.pupildiam)^2.))*
      eclat(dist(sim._size)^2.);
    x = fratio*tel.diam*(fratio*tel.diam-wfs(ns).l)/wfs(ns).l;
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
    if ( wfs(ns).skymag == 0) { wfs(ns)._skynphotons = 0.; } // if skymag not set
    
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

  err = _cwfs( ipupil, phase, phasescale, *wfs(ns)._tiltsh, *wfs(ns)._cxdef,
               *wfs(ns)._sxdef, dimpow2, *wfs(ns)._sind, *wfs(ns)._nsind,
               wfs(ns)._nsub, *wfs(ns)._fimage, *wfs(ns)._fimage2,
               float(wfs(ns)._nphotons), float(wfs(ns)._skynphotons),
               float(wfs(ns).ron), float(wfs(ns).darkcurrent*loop.ittime),
               int(wfs(ns).noise), mesvec);

  return mesvec;
}

//----------------------------------------------------

func zernike_wfs(pupsh,phase,ns,init=)
/* DOCUMENT
   Zernike WFS. Returns Zernike coefficients, expansion of the
   input phase onto Zernike modes.
   pupsh  = pupil
   phase  = phase
   ns     = WFS yao #
   init   = set to init the WFS. has to be called at least once.
   SEE ALSO:
 */
{
  // the phase at the input (call from multwfs) is in microns.
  // I have chosen to return coefficients of zernikes in nm (rms)
  
  extern wfs_zer,wfs_wzer,zn12;
  
  if (init) {
    pupd  = sim.pupildiam;
    size  = sim._size;
    nzer  = wfs(ns).nzer(1);
    wfs(ns)._nmes  = wfs(ns).nzer;
    cent  = sim._cent;
    prepzernike,size,pupd,sim._cent,sim._cent;
    wfs_wzer = where(zernike(1)*ipupil);
    wfs_zer = array(float,[2,numberof(wfs_wzer),nzer]);
    for (i=1;i<=nzer;i++) wfs_zer(,i) = zernike_ext(i)(*)(wfs_wzer);
    wfs_zer = LUsolve(wfs_zer(+,)*wfs_zer(+,))(+,)*wfs_zer(,+);
    // wfs_zer(nzer,npt in pupil)
    tmp = where(zernike(1)(avg,));
    zn12 = minmax(tmp);
    if (sim.verbose>=1) write,"Zernike wfs initialized";
    return;
  }

  wfs(ns)._fimage = wfs(ns)._dispimage = \
          &((phase*pupsh)(zn12(1):zn12(2),zn12(1):zn12(2)));
  mes = wfs_zer(,+)*phase(*)(wfs_wzer)(+);

  // returns microns rms (checked 2008apr10) ??? see above comment
  
  return mes;
}


//----------------------------------------------------

func mult_wfs_int_mat(disp=)
/* DOCUMENT func mult_wfs_int_mat(disp=)
   as mult_wfs but special for IntMat acquisition
   for speed in aoloop
   SEE ALSO:
 */
{
  extern wfs;
  mes = [];
  
  for (ns=1;ns<=nwfs;ns++) {

    // Impose noise = rmsbias = rmsflat = 0 for interaction matrix measurements
    // now done within do_imat (no need to do that at each call, and
    // this use to put a large overhead when using wfs.svipc (sync for every
    // actuators).
    //    noiseOrig = wfs(ns).noise; wfs(ns).noise = 0n;
    //    cycleOrig = wfs(ns).nintegcycles; wfs(ns).nintegcycles = 1;
    //    if (*wfs(ns)._skyfluxpersub!=[]) {
    //      skyfluxpersubOrig = *wfs(ns)._skyfluxpersub; *wfs(ns)._skyfluxpersub *= 0;
    //    }
    //    if (wfs(ns).type == "hartmann" ) {
    //      kconv_orig = wfs(ns)._kernelconv; wfs(ns)._kernelconv = 1n;
    //      bias  = *wfs(ns)._bias; *wfs(ns)._bias = *wfs(ns)._bias*0.0f;
    //      flat  = *wfs(ns)._flat; *wfs(ns)._flat = *wfs(ns)._flat*0.0f+1.0f;
    //    }

    //    offsets = wfs(ns).gspos;
    phase   = get_phase2d_from_dms(ns,"wfs");
    // uncomment if needed:
    //    phase  += get_phase2d_from_optics(ns,"wfs");

    // do the wavefront sensing:
    if (wfs(ns).type == "hartmann" ) {
      if (wfs(ns).disjointpup) {
        smes = sh_wfs(disjointpup(,,ns),phase,ns);
      } else {
        smes = sh_wfs(ipupil,phase,ns);
      }
    } else if (wfs(ns).type == "curvature") {
      smes = curv_wfs(pupil,phase,ns);
    } else if (wfs(ns).type == "pyramid") {
      smes = pyramid_wfs(pupil,phase);
    } else if (wfs(ns).type == "zernike") {
      smes = zernike_wfs(ipupil,phase,ns);
    } else {
      // assign user_wfs to requested function/type:
      cmd = swrite(format="user_wfs = %s",wfs(ns).type);
      include,[cmd],1;
      smes = user_wfs(ipupil,phase,ns);
    }

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
    // again, this was now moved to do_imat()
    //    if (wfs(ns).type == "hartmann" ) {
    //      wfs(ns)._bias = &bias; wfs(ns)._flat = &flat;
    //      wfs(ns)._kernelconv = kconv_orig;
    //    }
    //    wfs(ns).noise = noiseOrig;
    //    wfs(ns).nintegcycles = cycleOrig;
    //    if (*wfs(ns)._skyfluxpersub!=[]) *wfs(ns)._skyfluxpersub = skyfluxpersubOrig;
    //  }
    //  if (anyof(wfs.svipc>1)) status = sync_wfs_forks();
  }
  return mes;
}

//----------------------------------------------------

func mult_wfs(iter,disp=)
/* DOCUMENT func mult_wfs(iter,disp=)
   Goes through all WFS and concatenate the resulting measurement vectors.
   SEE ALSO:
 */
{
  mes = [];
  for (ns=1;ns<=nwfs;ns++) {

    //    offsets = wfs(ns).gspos;
    phase   = get_phase2d_from_dms(ns,"wfs");
    phase  += get_phase2d_from_optics(ns,"wfs");
    phase  += get_turb_phase(iter,ns,"wfs");

    if (wfs(ns).correctUpTT) {
      phase = correct_uplink_tt(phase,ns);
    }

    // get the measurements:
    if (wfs(ns).type == "hartmann" ) {
      if (wfs(ns).disjointpup) {
        smes = sh_wfs(disjointpup(,,ns),phase,ns);
      } else {
        smes = sh_wfs(ipupil,phase,ns);
      }
    } else if (wfs(ns).type == "curvature") {
      smes = curv_wfs(pupil,phase,ns);
    } else if (wfs(ns).type == "pyramid") {
      smes = pyramid_wfs(pupil,phase);
    } else if (wfs(ns).type == "zernike") {
      smes = zernike_wfs(ipupil,phase,ns);
    } else {
      // assign user_wfs to requested function/type:
      cmd = swrite(format="user_wfs = %s",wfs(ns).type);
      include,[cmd],1;
      smes = user_wfs(ipupil,phase,ns);
    }

    if (sim.svipc) {
      // hack, this has nothing to do here. FIXME
      shm_write,shmkey,swrite(format="wfs%d_image",ns),wfs(ns)._fimage;
    }
    
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

func shwfs_tests(void, clean=, wfs_svipc=, debug=, verbose=, batch=)
/* DOCUMENT shwfs_tests(void,clean=,wfs_svipc=,debug=,verbose=,batch=)

   Use to check that sh_wfs looks allright (NGS and LGS). Test a variety
   of features (sky, bias, flats, etc...) + background subtraction for
   each case.

   Use wfs_svipc= to test wfs.svipc. You might have to quit the yorick
   session between each test as somehow with these tests the svipc
   sometimes ends up in a non-correct state.
   
   Notes to myself (FR): For NGS, everything seems to be OK.  For LGS,
   the Rayleight looks allright, except that apparently the shadow by
   the central obstruction is not taken into account.
   SEE ALSO:
 */

{
  thispdiam = 120; //180 //240
  pltitle_height=12;
  //===========================
  aoread,"shwfs-tests.par";
  sim.pupildiam = thispdiam;
  sim.debug     = (debug?debug:0);
  sim.verbose   = (verbose?verbose:0);
  wfs.svipc     = (wfs_svipc?wfs_svipc:0)
  wfs(1).gsmag  = 8;
  wfs(1).skymag = 0;
  aoinit,dpi=90,clean=clean;
  winkill;
  window,0,wait=1,dpi=90,width=0,height=0;

  wfs(1).noise=0;
  shwfs_tests_plots,"NGS with no sky, no noise",batch=batch;
  
  wfs(1).noise=1;
  wfs(1).ron=0;
  shwfs_tests_plots,"NGS with no sky, w/ noise but no RON",batch=batch;
  
  wfs(1).ron=4;
  shwfs_tests_plots,"NGS with no sky, w/ noise",batch=batch;
  
  "ADDING SKY";
  wfs(1).skymag = 10;
  sim.debug = 0;
  aoinit;
  
  wfs(1).noise=0;
  shwfs_tests_plots,"NGS with w/ sky, no noise",batch=batch;
  
  wfs(1).noise=1;
  wfs(1).ron=0;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise but no RON",batch=batch;
  
  wfs(1).ron=4;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise",batch=batch;
  
  "TURNING OFF STAR";
  wfs(1).gsmag = 12;
  wfs(1).skymag = 10;
  sim.debug = 0;
  aoinit;
  
  wfs(1).noise=0;
  shwfs_tests_plots,"NO NGS with w/ sky, no noise",batch=batch;
  
  wfs(1).noise=1;
  wfs(1).ron=0;
  shwfs_tests_plots,"NO NGS with w/ sky, w/ noise but no RON",batch=batch;
  
  wfs(1).ron=4;
  shwfs_tests_plots,"NO NGS with w/ sky, w/ noise",batch=batch;
  
  "W/ BIAS";
  wfs(1).gsmag = 8;
  wfs(1).skymag = 10;
  wfs(1).biasrmserror = 50.;
  sim.debug = 0;
  aoinit;
  
  wfs(1).noise=0;
  shwfs_tests_plots,"NGS with w/ sky, no noise, BIAS error",batch=batch;
  
  wfs(1).noise=1;
  wfs(1).ron=0;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise but no RON, BIAS error",batch=batch;
  
  wfs(1).ron=4;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise, BIAS error",batch=batch;
  
  "W/ FLAT";
  wfs(1).gsmag = 8;
  wfs(1).skymag = 10;
  wfs(1).biasrmserror = 0.;
  wfs(1).flatrmserror = 0.3;
  sim.debug = 0;
  aoinit;
  
  wfs(1).noise=0;
  shwfs_tests_plots,"NGS with w/ sky, no noise, FLAT error",batch=batch;
  
  wfs(1).noise=1;
  wfs(1).ron=0;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise but no RON, FLAT error",batch=batch;
  
  wfs(1).ron=4;
  shwfs_tests_plots,"NGS with w/ sky, w/ noise, FLAT error",batch=batch;

  //===========================
  aoread,"shwfs-tests.par";
  loop.niter = 200;
  sim.debug     = (debug?debug:0);
  sim.verbose   = (verbose?verbose:0);
  wfs.svipc     = (wfs_svipc?wfs_svipc:0)
  sim.pupildiam = thispdiam;
  wfs.gsalt     = 90000;
  wfs.gsdepth   = 10000;
  wfs.laserpower = 10.;
  wfs.rayleighflag = 1;
  aoinit,disp=1,clean=clean;
  // aoloop,disp=1;
  
  shwfs_tests_plots,"2 LGSs with Rayleigh",batch=batch;
  
  wfs.noise=0;
  shwfs_tests_plots,"2 LGSs with Rayleigh no noise",batch=batch;
}

func shwfs_tests_plots(name,batch=)
{
  if (batch && (batch>1)) ptime=batch; else ptime=1000;
  
  wfs(1)._bckgrdsub  = 0;
  
  if (anyof(wfs.svipc>1)) status=sync_wfs_forks();
  
  if (wfs(1).disjointpup) {
    sh_wfs,disjointpup(,,1),pupsh*0.0f,1;
  } else sh_wfs,ipupil,ipupil*0.0f,1;
  
  tv,*wfs(1)._fimage;
  pltitle,name;
  stat,*wfs(1)._fimage;
  if (!batch) {
    r = strcase(0,kinput("Proceed/Spydr/show Bckgrdinit/Exit","P"));
    if (r=="s") {
      spydr,*wfs(1)._fimage;
      exit;
    }
    //~ if (r=="p") return;
    if (r=="e") exit;
  } else pause,ptime;
  
  wfs(1)._bckgrdsub  = 1;

  if (anyof(wfs.svipc>1)) status=sync_wfs_forks();
  
  if (wfs(1).disjointpup) {
    sh_wfs,disjointpup(,,1),pupsh*0.0f,1;
  } else sh_wfs,ipupil,ipupil*0.0f,1;
  tv,*wfs(1)._fimage;
  pltitle,name+" w/ bcksub=1";
  stat,*wfs(1)._fimage;
  if (batch) pause,ptime;
  else hitReturn;


  // 
  /*
  im = *wfs(1)._fimage;
  
  wfs(1)._bckgrdsub  = 1;
  wfs(1)._bckgrdinit = 1;
  oldnoise = wfs(1).noise;
  wfs(1).noise = 0;
  // call sh_wfs for calibration of the background
  sh_wfs,ipupil,ipupil*0.0f,1;
  wfs(1)._bckgrdinit = 0;   
  wfs(1).noise = oldnoise;
  tv,*wfs(1)._fimage;
  pltitle,name+" / bkgr calib";
  stat,*wfs(1)._fimage;
  hitReturn;

  tv,im-*wfs(1)._fimage;
  pltitle,name+" / image - bkgr calib";
  stat,im-*wfs(1)._fimage;
  hitReturn;
*/  
  
  
}

func svipc_time_sh_wfs(void,nit=,svipc=,doplot=)
{
  extern wfs;
  
  if (!nit) nit=100;
  if (svipc!=[]) wfs.svipc=svipc;

  t=array(0.,nit);
  prepzernike,sim._size, sim.pupildiam+4;
  phase = float(zernike(5));
  
  for(i=1;i<100;i++) {
    tic; r = sh_wfs(ipupil,phase,1); t(i)=tac()*1000.;
  }
  msg = swrite(format="%-20s %2d threads, sh_wfs() avg=%.2fms, median=%.2fms\n", \
               parprefix+":",wfs.svipc,avg(t(10:)),median(t(10:)));
  write,format="%s",msg;
  if (doplot) { plot,r; pause,20; }
  return msg;
}

func svipc_shwfs_tests(parfile,nfork_max=,doplot=)
/* DOCUMENT svipc_shwfs_tests(parfile,nfork_max=,doplot=)
   nfork_max = max number of fork (will do 1 to this number).
      If not provided, this function tries to get the number
      of processor (works only on linux systems), or failing so,
      will use a default of 4.
   doplot = enable to see the measurement vector after each test
            to see if everything looks right.
   General test routine for wfs.svipc
   use: grep -e "^processor" /proc/cpuinfo | wc
   to get the number of CPU on your machine (generally you'll want
   nfork_max equal to this number or slightly more for tests).
   SEE ALSO:
 */
{
  if (nfork_max==[]) {
    if (catch(0x10)) {
      nfork_max=4;
      aoread,parfile;
      sim.verbose = sim.debug =0;
      aoinit;
      for (nf=1;nf<=nfork_max;nf++) {
        grow,msg,svipc_time_sh_wfs(svipc=nf,doplot=doplot);
        if (nf>1) status = quit_wfs_forks();
      }
      write,msg;
      return;
    }
    nproc     = get_processor_number();
    nfork_max = nproc+2;
    write,format="\n>>> Found %d processor, using nfork_max = %d\n\n",nproc,nfork_max;
  }
  aoread,parfile;
  sim.verbose = sim.debug =0;
  aoinit;
  for (nf=1;nf<=nfork_max;nf++) {
    grow,msg,svipc_time_sh_wfs(svipc=nf,doplot=doplot);
    if (nf>1) status = quit_wfs_forks();
  }
  //  write,msg;
}

func get_processor_number(void)
// on linux
{
  if (f=open("/proc/cpuinfo","r",1)) {
    r = rdfile(f);
    close,f;
    nproc = sum(strgrep("^processor",r)(2,)>0);
    write,format="This machine has %d processors\n",nproc;
    return nproc;
  } else {
    error,"Can't determine number of processors";
  }
}
