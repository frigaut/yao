/*
 * yao.i
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

extern aoSimulVersion, aoSimulVersionDate;
aoSimulVersion = yaoVersion = aoYaoVersion = yao_version = "5.10.8";
aoSimulVersionDate = yaoVersionDate = aoYaoVersionDate = "2019nov29";

write,format=" Yao version %s, Last modified %s\n",yaoVersion,yaoVersionDate;

// find and include yao.conf now if any:
Y_CONF   = get_env("Y_CONF");
y_user   = streplace(Y_USER,strfind("~",Y_USER),get_env("HOME"))
if (noneof(Y_CONF)) \
  Y_CONF = "./:"+y_user+":"+pathform(_(y_user,Y_SITES,Y_SITE)+"conf/");

path2conf = find_in_path("yao.conf",takefirst=1,path=Y_CONF);
if (!is_void(path2conf)) {
  path2conf = dirname(path2conf)+"/";
  write,format=" Found and included yao.conf in %s\n",path2conf;
  require,path2conf+"yao.conf";
 }

// before we do anything else (to avoid forking a large process),
// let's try to determine what OS we are in:
//f = popen("uname",0);
//rep = ""; read,f,format="%s",rep;
//close,f;
//rep = strcase(0,rep);
//if ( (rep=="darwin") || (rep=="linux") ) os_env = rep; else os_env="unknown";
os_env="unknown";

plug_in,"yao";

require,"yao_utils.i";
require,"yao_fast.i";
require,"yao_wfs.i";
require,"yao_dm.i";
require,"aoutil.i";
require,"yao_gui.i";
require,"utils.i";
require,"yao_newfits.i";
require,"yao_util.i";
require,"yao_lgs.i";
require,"turbulence.i";
require,"plot.i";  // in yorick-yutils
require,"yao_structures.i";
require,"yaodh.i";

func LUsolve2(A,b){
/* DOCUMENT LUsolve2(a, b)
         or a_inverse= LUsolve2(a)

   returns the solution x of the matrix equation:
   A(,+)*x(+) = B
   using lpk_gesv in lapack for the matrix inversion

   If A is an n-by-n matrix then B must have length n, and the returned
   x will also have length n.

   SEE ALSO: LUsolve, lpk_gesv
*/
  local b;
  // b can be a matrix or a vector. If it does not exist, make it the identity matrix
  if (b == []){
    nr = (dimsof(A))(2);
    b = unit(nr);
    if (typeof(A) != "double"){
      b = float(unit(nr));
    }
  }

  return lpk_gesv(A,b);
}

// can't call use_sincos_approx directly in yao.conf, as not defined, work around:
if (use_sincos_approx_default!=[]) use_sincos_approx,use_sincos_approx_default;

// compatibility with GUI (yaopy.i)
func null (arg,..) { return 0; }


// set up for notify. Disabled by default. Set "use_notify" to enable
// we have to do that as soon as possible to avoid forking large a process
// this spawned bash session can actually be used for other purposes if
// needed.
func on_bash_out(msg) { write,msg; }
func notify(msg) {
  if (!use_notify) return; // allows to turn off/on notify within the session
  bash,swrite(format="notify-send -i %s \"%s\"\n",icon,msg);
}
icon = "/home/frigaut/Pictures/logos/yao_64_inv.png";
if (use_notify) bash = spawn("/bin/bash",on_bash_out); else bash=null;


func parse_yao_version(ver)
{
  extern yao_major_version,yao_minor_version;
  tmp = strtok(ver,".",3);
  yao_major_version = tonum(tmp(1));
  yao_minor_version = tonum(tmp(2));
}
parse_yao_version,yao_version;

// import fftw wisdom file:
if (fftw_wisdom_file==[]) fftw_wisdom_file=Y_USER+"fftw_wisdom_file.dat";
if (_import_wisdom(expand_path(fftw_wisdom_file))){
  write,format=" Warning: Can't read FFTW wisdom from %s\n",expand_path(fftw_wisdom_file);
  write, "Run init_fftw_wisdom to create a new file";
 }

// All below is designed to be overwritten by appropriate values
// in yaopy.i when using yao through the GUI
pyk_error = pyk_info = pyk_warning = null;
gui_message = gui_message1 = gui_progressbar_frac = gui_progressbar_text = null;
clean_progressbar = gui_show_statusbar1 = gui_hide_statusbar1 = pyk_flush = null;
if (YAO_SAVEPATH==[]) YAO_SAVEPATH = get_cwd();

//----------------------------------------------------
func comp_dm_shape_init(nm)
{
  extern wpupil;
  if (numberof(wpupil)!=ndm) wpupil = array(pointer,ndm);
  wpupil(nm) = &long(where((abs(*dm(nm)._def))(,,sum)));
}

func comp_dm_shape(nm,command,extrap=)
/* DOCUMENT comp_dm_shape(nm,command,extrap=)
   Fast compute of DM #nm shape from a command vector.
   Branch over _dmsum or _dmsumelt according to case.
   nm: DM yao #
   command: POINTER to a float vector containing the commands
            length of vector = numberof(*dm(nm)._command) = dm(nm)._nact
   extrap : Compute for extrapolated only (otherwise compute for valid only)
            so that to compute valid + extrap, you have to make 2 calls,
            one with and one without the extrap keyword set.
   SEE ALSO:
*/
{
  if (dmsum_use_new&&(wpupil==[])) for (i=1;i<=ndm;i++) comp_dm_shape_init,i;

  if (typeof(*command)!="float") error,"command is not float";

  n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = int(n2-n1+1);
  sphase = array(float,[2,nxy,nxy]);

  if (!is_set(extrap)) {
    // pegged valid actuators (right now, set to position 0):
    if (dm(nm).pegged) {
      com = *command;
      com(*dm(nm).pegged) = 0.0f;
      command = &com; // this way this does not go up in integrated commands
    }

    if (dm(nm)._flat_command) {
      com = *command;
      com -= float(*dm(nm)._flat_command);
      command = &com;
    }

    if (dm(nm).elt == 1) { //use fast dm shape computation
      _dmsumelt, dm(nm)._def, dm(nm)._eltdefsize, dm(nm)._eltdefsize, int(dm(nm)._nact),
        dm(nm)._i1, dm(nm)._j1, command, &sphase,nxy,nxy;
    } else { // use standard
      if (dmsum_use_new) _dmsum2, dm(nm)._def, wpupil(nm),              \
                           numberof(*wpupil(nm)), dm(nm)._nact, command, &sphase, numberof(sphase);
      else _dmsum, dm(nm)._def, nxy, nxy, dm(nm)._nact, command, &sphase;
      // }
    }

  } else { // extrapolated actuators

    // pegged extrapolated actuators (right now, set to position 0):
    if (dm(nm).epegged) {
      com = *command;
      com(*dm(nm).epegged) = 0.0f;
      command = &com; // this way this does not go up in integrated commands
    }

    if (dm(nm).elt == 1) { //use fast dm shape computation
      _dmsumelt, dm(nm)._edef, dm(nm)._eltdefsize, dm(nm)._eltdefsize, int(dm(nm)._enact),
        dm(nm)._ei1, dm(nm)._ej1, command, &sphase,nxy,nxy;
    } else { // use standard
      if (dmsum_use_new) _dmsum2, dm(nm)._edef, wpupil(nm), \
        numberof(*wpupil(nm)), dm(nm)._enact, command, &sphase, numberof(sphase);
      else _dmsum, dm(nm)._edef, nxy, nxy, dm(nm)._enact, command, &sphase;
    }

  }

  // temporary method for shifting DMs/influence functions
  // one cannot anymore use the straight (*dm._def)(,,sum),
  // but instead should always use comp_dm_shape.
  if (dm(nm)._puppixoffset!=[]) {
    sphase = roll(sphase,dm(nm)._puppixoffset);
  }

  if (dm(nm).disjointpup) { // we have to filter by disjointpup(,,nm)
    sphase *= disjointpup(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2,nm);
  }

  return sphase;
}
//----------------------------------------------------
func control_screen(i,init=)
/* DOCUMENT control_screen(i,init=)
   While the loop is running, brings up another graphic window in which
   some loop parameters are displayed. I rarely use this feature anymore.
   i is the loop counter. Not intended for use by the end user.
   SEE ALSO:
 */
{
  local x0,y0,x,y;
  y0 = 0.80;
  x0 = 0.150;
  ygstep = 0.03;
  mygstep = 0.022;
  sygstep = 0.017;
  csize = 12;
  if ( (xft!=[]) && (xft()==1) ) {
    csfont = "Courier:antialias=0";
  } else csfont="courier";

  /*  y0 = 0.86;
      x0 = 0.05;
      ygstep = 0.037;
      mygstep = 0.030;
      sygstep = 0.025;
      csize = 18;
  */

  if (is_set(init) || !window_exists(2)) {
    //    dheight = (nwfs+target._nlambda+1)*sygstep + 10*ygstep;
    //    dheight = long(dheight/0.77*570);
    if (window_exists(2)) winkill,2;  // make sure
    window,2,width=410,height=320,style="letter.gs",wait=1,dpi=70;
    //    window,2,width=570,height=dheight,style="letter.gs",wait=1,dpi=70;
    limits,0.,1.,0.,1.;
    progress_bar,0,pbid1,init=[x0,y0-ygstep,x0+0.25,y0-ygstep+0.015];
    if (is_set(init)) {
      window_select,0;
      return;
    }
  }

  window,2;
  fma;
  y = y0;
  if (!remainingTimestring) remainingTimestring="N/A";
  plt,sim.name+"  /  Time left: "+remainingTimestring,x0,y,tosys=1,\
    height=csize,justify="LA";
  y -= ygstep;
  it = swrite(format="%d/%d iterations",i,loop.niter);
  plt,it,x0+0.27,y,tosys=1,height=csize,justify="LA",font=csfont;
  y -= ygstep;
  progress_bar,i*100./loop.niter,pbid1;

  //  s = "WFS     rms     Tip     Tilt    ";
  s = "WFS        Tip    Tilt   ";
  if (anyof(wfs.correctUpTT)) {
    s = s+"UpTip  UpTilt   ";
    if (anyof(wfs.centGainOpt)) {
      s = s+"CentG   ";
    }
  }
  plt,s,x0,y,tosys=1,height=csize,justify="LA",font=csfont;
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
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font=csfont;
    y -= sygstep;
  }

  // print out TTM command
  wtt = where(dm.type == "tiptilt");
  if (numberof(wtt) != 0) {
    wtt = wtt(0);
    s = swrite(format="TT Mirror: % 6.3f  % 6.3f",(*dm(wtt)._command)(1),
               (*dm(wtt)._command)(2));
    y -= sygstep;
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font=csfont;
    y -= ygstep;
  }

  // print out Strehls:
  plt,"Strehl Ratio  Lambda     Avg     rms     Min     Max",x0,y,tosys=1,height=csize,
    justify="LA",font=csfont;
  y -= sygstep;

  for (jl=1;jl<=target._nlambda;jl++) {
    strehlle = imav(max,max,,jl)/sairy/(niterok+1e-5);
    s = swrite(format="Since Start   %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
               (*target.lambda)(jl),strehlle(avg),strehlle(rms),
               min(strehlle),max(strehlle));
    plt,s,x0,y,tosys=1,height=csize,justify="LA",font=csfont;
    y -= sygstep;
    if (jl == target._nlambda) {
      strehlse = im(max,max,)/sairy;
      s = swrite(format="Instantaneous %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
                 (*target.lambda)(jl),strehlse(avg),strehlse(rms),
                 min(strehlse),max(strehlse));
      plt,s,x0,y,tosys=1,height=csize,justify="LA",font=csfont;
      y -= sygstep;
    }
  }


  window_select,0;
}


//----------------------------------------------------

func mcao_rayleigh(nwfs,xsubap,ysubap,fov,aspp,zenith)
/* DOCUMENT mcao_rayleigh(nwfs,xsubap,ysubap,fov,aspp,zenith)
   Computes the rayleigh backscattering from other LGS beams (fratricide)
   Called from shwfs_init()
   SEE ALSO: shwfs_init
 */
{
  extern wfs;

  as2rd = dtor/3600.;
  cobs  = 0;

  zenith *= dtor;

  // position of GS/WFS in arcsec:
  w =  where(wfs._gsalt > 0);
  ns = where(w == nwfs)(1);

  xwfs = wfs(w).gspos(1,);
  ywfs = wfs(w).gspos(2,);
  nbeams = numberof(xwfs);

  // note that the incoming xsubap and ysubap have been switched!
  // hence
  lltx = wfs(w).LLTxy(2,);
  llty = wfs(w).LLTxy(1,);

  if (wfs(w)(ns).LLT1overe2diam==0) {wfs(w)(ns).LLT1overe2diam=0.3;}
  beamdiameter = wfs(w)(ns).LLT1overe2diam; // Gaussian laser beam FWHM [m]

  laserlambda  = wfs(w)(ns).lambda*1e-6; //589e-9;
  diamsubap    = tel.diam/wfs(w)(ns).shnxsub; // side of a subaperture [m]
  r0           = (tel.diam/atm.dr0at05mic)*cos(zenith)^0.6;
  seeing       = laserlambda/r0/4.848e-6;
  d            = 1e-2; // linear size of aperture for flux (/cm^2) <???

  altsod       = wfs(w)(ns)._gsalt;
  fwhmsod      = wfs(w)(ns)._gsdepth;

  // definitions of the image array to return to caller
  dim = long(ceil(fov/aspp));
  imrayl = array(float,[2,dim,dim]);
  imstar = array(float,[2,dim,dim]);

  // looking at GS #n with WFS #n -> some angle
  phin   = sqrt(xwfs(ns)^2.+ywfs(ns)^2.)*as2rd;
  thetan = atan(xwfs(ns),ywfs(ns));
  // range vector to sample h (for Rayleigh only)
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
      rayleigh = wfs(w)(beam).laserpower/(6.62e-34*3e8/589e-9)*sbnr*d^2/(4*pi*r^2.)*deltar*loop.ittime*wfs(w)(ns).optthroughput;

      if (rayleigh_fudge) rayleigh *= rayleigh_fudge;

      // compute position angles of the beam center for this altitude
      alpha = sin(phib)*cos(thetab)-sin(phin)*cos(thetan) - (xsubap-lltx(beam))/r + (xsubap-lltx(beam))/altsod;
      beta  = sin(phib)*sin(thetab)-sin(phin)*sin(thetan) - (ysubap-llty(beam))/r + (ysubap-llty(beam))/altsod;
      alpha /= 4.848e-6; // in arcsec
      beta  /= 4.848e-6;
      if (sim.debug > 3) {
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
    nPhotonFromSodStar = wfs(w)(beam).laserpower*gs.lgsreturnperwatt*loop.ittime*
      cos(zenith)*wfs(w)(ns).optthroughput;
    sodprofile = exp(-((rvecsod-altsod)/(fwhmsod/2.))^4.);
    sodprofile = sodprofile/sum(sodprofile)*nPhotonFromSodStar;

    for (i=1;i<=numberof(rvecsod)-1;i++) {

      // range
      r = rvecsod(i);
      sodium  = sodprofile(i);

      // compute position angles of the beam center for this altitude
      alpha = sin(phib)*cos(thetab)-sin(phin)*cos(thetan) - (xsubap-lltx(beam))/r + (xsubap-lltx(beam))/altsod;
      beta  = sin(phib)*sin(thetab)-sin(phin)*sin(thetan) - (ysubap-llty(beam))/r + (ysubap-llty(beam))/altsod;
      alpha /= 4.848e-6; // in arcsec
      beta  /= 4.848e-6;
      if (sim.debug > 3) {
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
    if (sim.debug > 3) {
      fma;
      pli,imrayl+imstar,-dim/2*aspp,-dim/2*aspp,+dim/2*aspp,+dim/2*aspp;
      pause,20;
    }
  }

  return [imrayl,imstar];
}

//----------------------------------------------------
func do_imat(disp=)
/* DOCUMENT do_imat(disp=)
   Measure the interaction matrix.
   Each actuator are actuated in a row, a measurement vector is taken
   and placed into the iMat. The reference (for phase=0) is subtracted.

   disp       = set to display stuff as it goes.

   This routine uses:
   - dm._nact, _n1, _n2, _def (extern)
   - mircube (extern)

   This routine calls:
   - mult_wfs_int_mat

   This routine sets:
   - iMat (extern)
   SEE ALSO: prep_svd, build_cmat
   */
{
  extern wfs;

  gui_progressbar_frac,0.;
  gui_progressbar_text,"Doing interaction matrix";

  if (mat.method == "mmse-sparse"){
    extern MR, MN;
    MR = mat.sparse_MR;
    MN = mat.sparse_MN;
  }

  indexDm       = array(long,2,ndm);
  indexDm(,1)   = [1,dm(1)._nact];
  for (nm=2;nm<=ndm;nm++) {
    indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
  }

  ntot = sum(dm._nact);
  ncur=1;

  if (disp) { plsys,1; animate,1; }
  // save state of noise/nintegcycle/etc: everything that is not desired
  // when doing the iMat:
  store_noise_etc_for_imat,noise_orig, cycle_orig, kconv_orig, \
              skyfluxpersub_orig, bckgrdcalib_orig, bias_orig, \
              flat_orig,darkcurrent_orig,use_sincos_approx_orig,\
              rayleigh_orig;

  // sync forks if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  // for imat ETA:
  tic; ndone=0;

  // Loop on each mirror:
  for (nm=1;nm<=ndm;nm++) {
    n1 = dm(nm)._n1;
    n2 = dm(nm)._n2;
    //    if (sim.verbose>1) {write,format="\rDoing DM# %d, actuator %s",nm," ";}
    if (sim.verbose) write,"";
    subsys = dm(nm).subsystem;
    // Loop on each actuator:
    command = array(float,dm(nm)._nact);

    for (i=1;i<=dm(nm)._nact;i++) {

      ncur++;
      gui_progressbar_frac,float(ncur)/ntot;
      gui_progressbar_text,\
        swrite(format="Doing interaction matrix, DM#%d, actuator#%d",nm,i);
      if (sim.verbose) {
        eta = (ndone?(ntot-ndone)*tac()/ndone:0.0f);
        write,format="\rDoing DM# %d, actuator %d/%d, ETA %.1fs",\
          nm,i,dm(nm)._nact,eta;
      }
      //      if (sim.verbose>1) {write,format="%d ",i;}

      mircube  *= 0.0f; command *= 0.0f;
      command(i) = float(dm(nm).push4imat);
      mircube(n1:n2,n1:n2,nm) = comp_dm_shape(nm,&command);

      if (mat.method != "mmse-sparse"){
        if (!dm(nm).ncp){
          // Fill iMat (reference vector subtracted in mult_wfs_int_mat):
          iMat(,i+indexDm(1,nm)-1) = mult_wfs_int_mat(disp=disp,subsys=subsys)/dm(nm).push4imat;
        }
      } else {
        if (!dm(nm).ncp){
          rcobuild, iMatSP, float(mult_wfs_int_mat(disp=disp,subsys=subsys)/dm(nm).push4imat), mat.sparse_thresh;
        } else {
          rcobuild, iMatSP, float(mult_wfs_int_mat(disp=disp,subsys=subsys)/dm(nm).push4imat*0.), mat.sparse_thresh;
        }
      }

      // display, if requested:
      // WFS spots:
      if (is_set(disp)) {
        fma;
        plt,escapechar(swrite(format="YAO %s | %s | %s",aoSimulVersion,parprefix,sim.name)),0.01,0.227,tosys=0;
        if (!allof(wfs.shmethod ==1)) {
          if (wfs_display_mode=="spatial") {
            disp2d,wfs._dispimage,wfs.pupoffset(1,),wfs.pupoffset(2,),2;
            mypltitle,"WFSs spots (spatial mode)",[0.,-0.005],height=12;
          } else {
            disp2d,wfs._dispimage,wfs.gspos(1,),wfs.gspos(2,),2;
            mypltitle,"WFSs spots",[0.,-0.005],height=12;
          }
        }
        if ((nm==1) && (i==1)) { plsys,2; limits,square=1; limits; }
        // mirror surface
        plsys,3;
        limits,square=1;
        if (mergedms4disp) { //e.g. GMT case
          tmp = mircube(,,sum);
          pli,(tmp-min(tmp))*ipupil,1,0.,2,1.;
          range,0.25,0.75;
        } else {
          pli,(mircube(,,1)-min(mircube(,,1)))*ipupil,1,0.,2,1.;
          if (ndm==1) range,0.25,0.75;
          for (jj=2;jj<=ndm;jj++) {
            //pli,mircube(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2,nm),nm,0.,nm+1,1.;
            if (dm(jj).alt==0) {
              pli,(mircube(,,jj)-min(mircube(,,jj)))*ipupil,jj,0.,jj+1,1.;
            } else {
              pli,(mircube(,,jj)-min(mircube(,,jj))),jj,0.,jj+1,1.;
            }
          }
        }
        // mypltitle,"DM(s)",[0.,0.008],height=12;
        myxytitles,"","DM(s)",[0.010,0.],height=12;

        if ((dm(nm).type == "aniso") && (sim.debug >= 2)) hitReturn;
      } // end of display section

      if (sleep) usleep,sleep;
      if ((sim.debug>=3) && (i==(dm(nm)._nact/2))) hitReturn;
      ndone++;
    }
    if (sim.verbose) {write," ";}
  }

  // restore original values to WFS structure:
  restore_noise_etc_for_imat,noise_orig, cycle_orig, kconv_orig, \
                 skyfluxpersub_orig,bckgrdcalib_orig, bias_orig, \
                 flat_orig,darkcurrent_orig,use_sincos_approx_orig,\
                 rayleigh_orig;

  // sync forks if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  // Display if needed:
  if ((mat.method != "mmse-sparse") && ((sim.debug>=1) || (disp == 1))) {
    tv,-iMat,square=1;
    mypltitle,"Interaction Matrix",[0.,-0.005],height=12;
    if (sim.debug >= 1) typeReturn;
  }

  if (disp) {plsys,1; animate,0; }

  clean_progressbar;
}


func store_noise_etc_for_imat(&noise_orig, &cycle_orig, &kconv_orig,
                              &skyfluxpersub_orig, &bckgrdcalib_orig,
                              &bias_orig, &flat_orig, &darkcurrent_orig,
                              &use_sincos_approx_orig,&rayleigh_orig)
{
  extern wfs;

  noise_orig = cycle_orig = kconv_orig = array(0n,nwfs);
  darkcurrent_orig = rayleigh_orig = array(float,nwfs);
  skyfluxpersub_orig = bckgrdcalib_orig = bias_orig = flat_orig = array(pointer,nwfs);
  use_sincos_approx_orig = [use_sincos_approx()]; // & need a vector
  use_sincos_approx,0;

  for (ns=1;ns<=nwfs;ns++) {

    // Impose noise = rmsbias = rmsflat = 0 for interaction matrix measurements
    noise_orig(ns) = wfs(ns).noise;
    wfs(ns).noise = 0n;

    cycle_orig(ns) = wfs(ns).nintegcycles;
    wfs(ns).nintegcycles = 1n;

    darkcurrent_orig(ns) = wfs(ns).darkcurrent;
    wfs(ns).darkcurrent = float(0.);

    rayleigh_orig(ns) = wfs(ns).rayleighflag;
    wfs(ns).rayleighflag = 0n;

    if (*wfs(ns)._skyfluxpersub!=[]) {
      skyfluxpersub_orig(ns) = &(*wfs(ns)._skyfluxpersub);
      *wfs(ns)._skyfluxpersub *= 0;
    }

    if (*wfs(ns)._bckgrdcalib!=[]) {
      bckgrdcalib_orig(ns) = &(*wfs(ns)._bckgrdcalib);
      *wfs(ns)._bckgrdcalib *= 0;
    }

    if (wfs(ns).type == "hartmann" ) {
      kconv_orig(ns) = wfs(ns)._kernelconv;
      wfs(ns)._kernelconv = 1n;

      bias_orig(ns)  = &(*wfs(ns)._bias);
      *wfs(ns)._bias = *wfs(ns)._bias*0.0f;

      flat_orig(ns)  = &(*wfs(ns)._flat);
      *wfs(ns)._flat = *wfs(ns)._flat*0.0f+1.0f;
    }
  }
}

func restore_noise_etc_for_imat(noise_orig, cycle_orig, kconv_orig,
                                skyfluxpersub_orig,bckgrdcalib_orig,
                                bias_orig, flat_orig, darkcurrent_orig,
                                use_sincos_approx_orig,rayleigh_orig)
{
  extern wfs;
  use_sincos_approx,use_sincos_approx_orig(1);

  for (ns=1;ns<=nwfs;ns++) {

    wfs(ns).noise = noise_orig(ns);
    wfs(ns).darkcurrent = darkcurrent_orig(ns);
    wfs(ns).rayleighflag = rayleigh_orig(ns);
    wfs(ns).nintegcycles = cycle_orig(ns);

    if (*wfs(ns)._skyfluxpersub!=[])                    \
      wfs(ns)._skyfluxpersub = skyfluxpersub_orig(ns);
    if (*wfs(ns)._bckgrdcalib!=[])                    \
      wfs(ns)._bckgrdcalib = bckgrdcalib_orig(ns);

    if (wfs(ns).type == "hartmann" ) {
      wfs(ns)._kernelconv = kconv_orig(ns);
      wfs(ns)._bias = bias_orig(ns);
      wfs(ns)._flat = flat_orig(ns);
    }
  }
}

//----------------------------------------------------
func prep_svd(imat,subsystem,svd=,disp=)
/* DOCUMENT prep_svd(imat,subsystem,svd=,disp=)
   Does the SVD decomposition and fill out the modToAct (V)
   and mesToMod (UT) matrices for further use by build_cmat()

   imat = the imat to inverse
   disp = set if display is required.

   This routine uses:
   - imat (input)

   This routine calls:
   - SVdec

   This routine sets:
   - eigenvalues (extern)
   - modToAct (extern)
   - mesToMod (extern)

  Note: The Mode-to-Actuator (modToAct) matrix has to be used as follow:
  modes-coef    = actToMod(,+) * command-coef(+)

  SEE ALSO: do_imat, build_cmat
*/
{
  // Define some extern variables:
  extern modToAct,mesToMod,eigenvalues;

  // Decompose to prepare inversion:
  if (sim.verbose>1) {write,"Doing SVD\n";}
  eigenvalues   = SVdec(imat,u,vt);

  // Some debug output if needed:
  if (sim.verbose) {
    write,format="Smallests 2 normalized eigenvalues = %f",
      eigenvalues(-1)/max(eigenvalues),eigenvalues(0)/max(eigenvalues);
  }

  if (sim.verbose>1) {
    write,"Normalized eigenvalues:";
    write,eigenvalues/max(eigenvalues);
    th = 1.0f/(*mat.condition)(subsystem);
    do {
      plot,eigenvalues/max(eigenvalues);
      plg,array(1/(*mat.condition)(subsystem),numberof(eigenvalues)),color="red",type=2;
      limits,square=0;
      mypltitle,"Normalized eigenvalues",[0.,-0.005],height=12;
      if (yaopy) break;
      if (is_set(svd)) {
        th = kinput("New Threshold ? (return to continue)",th);
        change = ( (th==1/(*mat.condition)(subsystem))? 0:1);
        if (change) (*mat.condition)(subsystem) = 1/th;
      }
    } while (change == 1);
    if (!yaopy) typeReturn;
  }

  modToAct    = transpose(vt);
  //  actToMod    = LUsolve(modToAct);
  mesToMod    = transpose(u);  // used to be called ut

  // Some debug display if needed:
  if (sim.debug>=2) {
    tv,modToAct;
    mypltitle,"modToAct Matrix",[0.,-0.005],height=12;
    if (!yaopy) typeReturn;
    if (sim.debug>=3) {
      tv,modToAct(,+)*modToAct(,+);
      mypltitle,"modToAct(,+)*modToAct(,+)",[0.,-0.005],height=12;
      if (!yaopy) typeReturn;
    }

    tv,mesToMod;
    mypltitle,"mesToMod Matrix",[0.,-0.005],height=12;
    if (!yaopy) typeReturn;
    if (sim.debug>=3) {
      tv,mesToMod(,+)*mesToMod(,+);
      mypltitle,"mesToMod(,+)*mesToMod(,+)",[0.,-0.005],height=12;
      if (!yaopy) typeReturn;
    }
  }
}

//----------------------------------------------------
// func svdmodes_variance(void)
// {
//   prepzernike,sim._size,sim.pupildiam,sim._size/2.+0.5,sim._size/2.+0.5;
//   fma; pli,zernike(1)+ipupil
// }
//----------------------------------------------------

func build_cmat(condition,modalgain,subsystem=,all=,nomodalgain=,disp=)
/* DOCUMENT build_cmat(condition,modalgain,subsystem=,all=,nomodalgain=,disp=)
   Build the command matrix from V,UT, the eigenvalues and the modal gains
   F.Rigaut, June 17,2002

   condition = Filter eigenvalues ev (and modes) for max(ev)/ev > condition
   modalgain = vector of system mode gains to pre-multiply the control matrix
               with. rarely used as these modes are not generally natural
               w.r.t. the turbulence.
   all = if not set, one (1) mode is forced to be discarded (useful if regular
         AO system using a DM with a piston component and condition number
         too large). Normally, set all=1.
   nomodalgain = if set, the modal gain are not taken into account.
   disp = set to display stuff.

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

   SEE ALSO: do_imat, prep_svd
*/
{
  extern NModesControlled;

  neigen = numberof(eigenvalues);

  mev   = array(float,neigen,neigen);

  mask = ((eigenvalues/max(eigenvalues)) > (1./condition));

  if (numberof(ev_user_mask)==numberof(mask)) mask = (mask|ev_user_mask);

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

  if (sim.verbose) {
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
    disp2d,cubphase,subpos(1,),subpos(2,),1,zoom=0.9,init=1;

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

          mircube(n1:n2,n1:n2,nm) = comp_dm_shape(nm,&scommand);
        }
      }

      // Loop over WFS to get integrated phase
      n1 = (sim._size-sim.pupildiam)/2-2;
      n2 = (sim._size+sim.pupildiam)/2+2;
      nxy = n2-n1+1;
      nnss = 1;
      for (ns=1;ns<=nwfs;ns++) {
        if (wfs(ns).subsystem == subsystem) {
          phase = get_phase2d_from_dms(ns,"wfs")*ipupil;
          // fill cubphase
          cubphase(,,nnss) = phase(n1:n2,n1:n2);
          nnss++;
        }
      }

      // display using disp2d of integrated phases
      fma;
      plt,escapechar(swrite(format="YAO %s | %s | %s",aoSimulVersion,parprefix,sim.name)),0.01,0.227,tosys=0;
      disp2d,cubphase,subpos(1,),subpos(2,),1;
      mypltitle,swrite(format="Mode %d, Normalized eigenvalue = %f",
                     i,eigenvalues(i)/max(eigenvalues)),[0.,-0.005],height=12;

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
      if (rep == "q") break;
    }
  }
  //=========================================
  //             END OF DISPLAY
  //=========================================


  return cmat;
}

//----------------------------------------------------
func swap_screens(void)
/* DOCUMENT swap_screens
   Swap the phase screens. This is to get better statistics
   out of fewer phase screens and iterations.
   The 2nd phase screen becomes the 1rst one, 3->2, etc...
   This routine uses the phase screens and normalization
   factor stored in extern by get_turb_phase_init
   SEE ALSO:
*/
{
  extern pscreens;
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

func get_turb_phase_init(skipReadPhaseScreens=)
/* DOCUMENT get_turb_phase_init(skipReadPhaseScreens=)
   Initializes everything for get_turb_phase (see below), which
   returns the interpolated, integrated phase to the loop function
   for iteration number "iter". Returns the "ok" parameter which is
   set to 0 if this step is not to be used to compute statistics.
   AUTHOR: F.Rigaut, June 11, 2002.
   SEE ALSO: aoinit, aoloop, get_turb_phase.
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
                  // and therefore we cannot wrap.
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
    ymargins,     // same for Y. We use that in get_turb_phase to determine when to
                  // wrap.
    //~ statsokvec,   // 1 if it is ok to collect stats at this iteration.
                  //~ // 0 if not, e.g. we just did a jump. dim [iteration]
    inithistory,  // 1 if init has been done.
    screendim;    // [phase screen X dim, Y dim] before they are extended for safe wrapping

  extern currentScreenNorm; // current screen normalization. Used when swaping screen.

  // Define a few variables:

  nscreens = numberof(*atm.screen);
  if (opt!=[]) noptics = numberof(opt); else noptics=0;

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

    if (sim.verbose) {
      write,format="Reading phase screen \"%s\"\n",(*atm.screen)(1);
    }

    // read the first one, determine dimensions, stuff it:
    tmp      = yao_fitsread((*atm.screen)(1));
    dimx     = dimsof(tmp)(2);
    dimy     = dimsof(tmp)(3);
    screendim = [dimx,dimy];

    if (dimx == dimy){ // can wrap both x and y
      // Extend dimension in X and Y for wrapping issues
      // Made larger to accommodate GLAO, Marcos van Dam, May 2012
      pscreens = array(float,[3,dimx+4*sim._size,dimy+4*sim._size,nscreens]);
      // Stuff it
      pscreens(1:dimx,1:dimy,1) = tmp;
      // free RAM
      tmp      = [];

      // Now read all the other screens and put in pscreens
      for (i=2;i<=nscreens;i++) {
        if (sim.verbose) {
          write,format="Reading phase screen \"%s\"\n",(*atm.screen)(i);
        }
        pscreens(1:dimx,1:dimy,i) = yao_fitsread((*atm.screen)(i));
      }

      // Extend the phase screen length for safe wrapping:
      pscreens(dimx+1:,,) = pscreens(1:4*sim._size,,);
      pscreens(,dimy+1:,) = pscreens(,1:4*sim._size,);

    } else {
      // Extend dimension in X for wrapping issues
      pscreens = array(float,[3,dimx+2*sim._size,dimy,nscreens]);
      // Stuff it
      pscreens(1:dimx,,1) = tmp;
      // free RAM
      tmp      = [];

      // Now read all the other screens and put in pscreens
      for (i=2;i<=nscreens;i++) {
        if (sim.verbose) {
          write,format="Reading phase screen \"%s\"\n",(*atm.screen)(i);
        }
        pscreens(1:dimx,,i) = yao_fitsread((*atm.screen)(i));
      }

      // Extend the phase screen length for safe wrapping:
      pscreens(dimx+1:,,) = pscreens(1:2*sim._size,,);
      // Can't do in Y as the phase screens are not periodic (they have been cutted)
      //  pscreens(,dimy+1:,) = pscreens(,1:sim._size,);

    }

    // apply weights to each phase screens (normalize):
    // the screens are expressed in microns
    if (*atm.screen_norm!=[]) {
      weight = float(*atm.screen_norm);
      write,format="%s","Using special normalisation: "; weight;
    }
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
      psize  = tel.diam/sim.pupildiam;
      if (sim.verbose) {
        write,"Reading optics. I am expecting phase maps with a pixel";
        write,format=" size of %fm/pixel, as projected in the entrance\n",psize;
        write,"pupil plane. It is also your responsibility to provide phase";
        write,"maps of adequate dimension, i.e. large enough";
        write,format="Reading phase map for optics \"%s\"\n",opt(1).phasemaps;
      }
      tmp = yao_fitsread(YAO_SAVEPATH+opt(1).phasemaps);
      optdims      = dimsof(tmp);
      optdimx      = optdims(2);
      optdimy      = optdims(3);

      if (optdimx != optdimy)
        error,"Optics phase maps should be square arrays";

      optphasemaps = array(float,[3,optdimx,optdimy,noptics]);
      // Stuff it
      optphasemaps(,,1) = float(tmp);

      for (i=2;i<=noptics;i++) {
        if (sim.verbose) {
          write,format="Reading phase map for optics \"%s\"\n",opt(i).phasemaps;
        }
        tmp = yao_fitsread(YAO_SAVEPATH+opt(i).phasemaps);
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
  if (loop.skipevery > 0){
    iposvec = iposvec+((iposvec-1)/loop.skipevery)*loop.skipby;
  }

  // Build the position vector vs iteration by phase screen
  deltax  = sim.pupildiam/tel.diam*(*atm.layerspeed)*cos(dtor*gs.zenithangle)*loop.ittime;
  // has to start away from edge as this is the coordinate of the beam center
  // will modified later on when we know the beams geometry. see few lines below.
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
    if (dimx == dimy){yposvec(,ns) = ymargins(1)+ ( (yposvec(,ns)-ymargins(1)) % screendim(2));}
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
  // and a int vector speficying which layer/dm to skip,
  // only a set of X and Y positions as input [the one we just talked
  // about, xposvec(iteration) + wfsxposcub(,,wfs#) ].

  // check that y index does not overflow:
  if (dimx != dimy){get_turb_phase_initCheckOverflow;}

  inithistory = 1;
  return 1;
}
//----------------------------------------------------

func get_turb_phase(iter,nn,type)
/* DOCUMENT get_turb_phase(iter,n,type)
   Returns the interpolated, integrated phase along the turbulent phase
   screen data cube to the loop function for iteration number "iter".
   You have to call get_turb_phase_init to initialize prior using.
   SEE ALSO: aoinit, aoloop, get_turb_phase_init.
*/
{
  if (!inithistory) {error,"get_turb_phase has not been initialized !";}

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
    // mod 2011jan19 w/ DG to get rid of screens above LGS
    if (wfs(nn).gsalt>0) nscreens = long(sum(*atm.layeralt < wfs(nn).gsalt));
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = wfsxposcub(,,nn)+xposvec(iter,)(-,);
    yshifts = wfsyposcub(,,nn)+yposvec(iter,)(-,);
  } else if ( type == "target") {
    // stuff xshifts with fractionnal offsets, add xposvec for each screen
    xshifts = gsxposcub(,,nn)+xposvec(iter,)(-,);
    yshifts = gsyposcub(,,nn)+yposvec(iter,)(-,);
    ppm = sim.pupildiam/tel.diam;
    if (target.xspeed&&loopCounter) {
      xss = (*target.xspeed)(nn)*4.848e-6*(*atm.layeralt)*loop.ittime*loopCounter*ppm;
      xshifts += float(xss)(-,);
    }
    if (target.yspeed&&loopCounter) {
      yss = (*target.yspeed)(nn)*4.848e-6*(*atm.layeralt)*loop.ittime*loopCounter*ppm;
      yshifts += float(yss)(-,);
    }
  }
  skip = array(0n,nscreens);

  xshifts = xshifts - screendim(1)*(xshifts>=(psnx-1));
  yshifts = yshifts - screendim(2)*(yshifts>=(psny-1));

  ishifts = int(xshifts);  xshifts = xshifts - ishifts;
  jshifts = int(yshifts);  yshifts = yshifts - jshifts;

  err = _get2dPhase(&pscreens,psnx,psny,nscreens,&skip,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  if (err != 0) {error,"Error in get_turb_phase";}

  bphase(_n1:_n2,_n1:_n2) = sphase;

  return bphase;
}

//----------------------------------------------------
func get_phase2d_from_dms(nn,type)
/* DOCUMENT get_phase2d_from_dms(nn, type)
   Similar than get_turb_phase excepts it works on the mirror shape
   data cube. Returns the integrated phase along one given direction.

   nn = yao object # of type as below
   type = "wfs" or "target".
   so get_phase2d_from_dms(2,"wfs") will return the phase integrated along
   the third dimension of the cube (different altitude) in the direction
   of wfs #2.
   uses mircube (extern) = cube or image of phase screen (mirrors)
   first plan = first mirror shape
   second plan = second mirror shape
   etc...
   should be of size sim._size * sim._size
   SEE ALSO:
*/
{
  sphase = array(float,_n,_n);
  bphase = array(float,sim._size,sim._size);

  if (anyof(dm.ncp)||(target.ncpdm)||((type=="wfs")&&(wfs(nn).ncpdm))) {
    mirrorcube = mircube; // create a copy
    dmidx = where(dm.ncp == 0);

    if (type == "target"){
      if (*target.ncpdm != []){
        ncpdm = (*target.ncpdm)(nn);
        if (ncpdm != 0 && noneof(dmidx == ncpdm)){
          grow, dmidx, ncpdm;
        }
      }
    } else { // type == "wfs"
      if (wfs(nn).ncpdm != 0 && noneof(dmidx == wfs(nn).ncpdm)){
        grow, dmidx, wfs(nn).ncpdm;
      }
    }

    mirrorcube = mirrorcube(,,dmidx);
  } else eq_nocopy,mirrorcube,mircube;

  // Now we can call the C interpolation routine and get the integrated
  // phase for this star
  // there are a few things to do to get ready
  psnx = dimsof(mirrorcube)(2);
  psny = dimsof(mirrorcube)(3);
  nmirrors = dimsof(mirrorcube)(4);

  // here we have a branch to be able to process wfs and targets with the same
  // subroutine, this one.
  if (type == "wfs") {
    // stuff xshifts with fractional offsets, add xposvec for each screen
    xshifts = dmwfsxposcub(,,nn)+(sim._cent+dm.misreg(1,)-1)(-,);
    yshifts = dmwfsyposcub(,,nn)+(sim._cent+dm.misreg(2,)-1)(-,);
    skip = int(*wfs(nn)._dmnotinpath);
  } else if ( type == "target") {
    // stuff xshifts with fractional offsets, add xposvec for each screen
    xshifts = dmgsxposcub(,,nn)+(sim._cent+dm.misreg(1,)-1)(-,);
    yshifts = dmgsyposcub(,,nn)+(sim._cent+dm.misreg(2,)-1)(-,);
    ppm = sim.pupildiam/tel.diam;
    if (target.xspeed&&loopCounter) {
      xss = (*target.xspeed)(nn)*4.848e-6*dm.alt*loop.ittime*loopCounter*ppm;
      xshifts += float(xss)(-,);
    }
    if (target.yspeed&&loopCounter) {
      yss = (*target.yspeed)(nn)*4.848e-6*dm.alt*loop.ittime*loopCounter*ppm;
      yshifts += float(yss)(-,);
    }
    skip = array(0n,nmirrors);
  }

  ishifts = int(xshifts); xshifts = xshifts - ishifts;
  jshifts = int(yshifts); yshifts = yshifts - jshifts;

  err = _get2dPhase(&mirrorcube,psnx,psny,nmirrors,&skip,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  if (err != 0) {error,"Error in get_phase2d_from_dms";}

  bphase(_n1:_n2,_n1:_n2) = sphase;

  return bphase;
}

//----------------------------------------------------
func get_phase2d_from_optics(nn,type)
/* DOCUMENT get_phase2d_from_optics(nn,type)
   Same as get_phase2d_from_dms, but for static optical elements defined
   in the opt structure. See get_phase2d_from_dms.
   nn = wfs or GS #
   type = "wfs" or "target"
   SEE ALSO:
*/
{
  if (opt==[]) return 0.0f;

  // select optics in path that apply to type:
  // there are noptics = numberof(opt)
  if (type=="wfs") {
    w = where(opt.path_type != "target");
    if (numberof(w)==0) return 0.0f; // no common or wfs optics, return

    // now for each optics, check this wfs has not been excluded
    for (no=1;no<=numberof(w);no++) {
      if (opt(w(no)).path_type=="common"){continue;}
      if (noneof(*opt(w(no)).path_which==nn)){
        w(no)=-1;
        continue;
      }
      if (opt(w(no)).scale==0){
        w(no)=-1;
        continue;
      }
    }
    w = w(where(w>=0));
    if (numberof(w)==0) return 0.0f;
  }

  if (type=="target") {
    w = where(opt.path_type!="wfs");
    if (numberof(w)==0) return 0.0f; // no common or target optics, return

    // now for each optics, check this target has not been excluded
    for (no=1;no<=numberof(w);no++) {
      if (opt(w(no)).path_type=="common"){continue;}
      if (noneof(*opt(w(no)).path_which==nn)){
        w(no)=-1;
        continue;
      }
      if (opt(w(no)).scale==0){
        w(no)=-1;
        continue;
      }
    }
    w = w(where(w>=0));
    if (numberof(w)==0) return 0.0f;
  }

  sphase = array(float,_n,_n);
  bphase = array(float,sim._size,sim._size);

  if (typeof(optphasemaps) != "float"){optphasemaps = float(optphasemaps);}
  thisphasemaps = optphasemaps(,,w)*opt(w).scale(-,-,);

  // Now we can call the C interpolation routine and get the integrated
  // phase for this star
  // there are a few things to do to get ready
  psnx  = dimsof(thisphasemaps)(2);
  psny  = dimsof(thisphasemaps)(3);
  nopts = dimsof(thisphasemaps)(4);

  // here we have a branch to be able to process wfs and targets with the same
  // subroutine, this one.
  if (type == "wfs") {
    // stuff xshifts with fractional offsets, add xposvec for each screen
    xshifts = optwfsxposcub(,w,nn)+(opt._cent+opt.misreg(1,)-1)(-,w);
    yshifts = optwfsyposcub(,w,nn)+(opt._cent+opt.misreg(2,)-1)(-,w);
  } else if ( type == "target") {
    // stuff xshifts with fractional offsets, add xposvec for each screen
    xshifts = optgsxposcub(,w,nn)+(opt._cent+opt.misreg(1,)-1)(-,w);
    yshifts = optgsyposcub(,w,nn)+(opt._cent+opt.misreg(2,)-1)(-,w);
    ppm = sim.pupildiam/tel.diam;
    if (target.xspeed&&loopCounter) {
      xss = (*target.xspeed)(nn)*4.848e-6*opt.alt*loop.ittime*loopCounter*ppm;
      xshifts += float(xss)(-,);
    }
    if (target.yspeed&&loopCounter) {
      yss = (*target.yspeed)(nn)*4.848e-6*opt.alt*loop.ittime*loopCounter*ppm;
      yshifts += float(yss)(-,);
    }
  }
  skip = array(0n,nopts);

  ishifts = int(xshifts); xshifts = float(xshifts - ishifts);
  jshifts = int(yshifts); yshifts = float(yshifts - jshifts);

  err = _get2dPhase(&thisphasemaps,psnx,psny,nopts,&skip,
                    &sphase,_n,_n,
                    &ishifts,&xshifts,
                    &jshifts,&yshifts);

  //  if (err != 0) {error,"Error in get_phase2d_from_optics";}

  bphase(_n1:_n2,_n1:_n2) = sphase;

  return bphase;
}

//----------------------------------------------------
func correct_uplink_tt(phase, ns)
/* DOCUMENT correct_uplink_tt(phase, ns)
   Apply tiptilt phase correction to a phase term, given
   wfs._upttcommand determined previously.
   SEE ALSO:
*/
{
  wfs(ns)._upttcommand += wfs(ns).uplinkgain * wfs(ns)._tt;

  phase -= wfs(ns)._upttcommand(1) * tip1arcsec;
  phase -= wfs(ns)._upttcommand(2) * tilt1arcsec;

  return phase;
}
//----------------------------------------------------
func split_wfs_vector(v)
/* DOCUMENT split_wfs_vector(v)
   Splits a single vector (out of mult_wfs or mult_wfs_int_mat)
   into as many individual wfs vectors as there are sensors.
   Return a pointer vector to the individual wfs vectors.
   SEE ALSO: split_dm_vector(v)
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
func split_dm_vector(v)
/* DOCUMENT split_dm_vector(v)
   Splits the single vector (out of cMat matrix multiply in aoloop)
   into as many individual command vector as there are DMs.
   Return a pointer vector to the individual command vectors.
   SEE ALSO: split_wfs_vector
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
func check_control_parameters(void,verbose=){
/* DOCUMENT check_control_parameters(void)
   Check the control parameters and update them
   SEE ALSO:
 */

// first, convert the controllers into what is implemented
  gainho = *loop.gainho;
  leakho = *loop.leakho;
  for (nm=1;nm<=ndm;nm++){
    ctrlnum = *dm(nm).ctrlnum;
    ctrlden = *dm(nm).ctrlden;

    if (numberof(ctrlnum) > 10){
      write, format = "dm(%d).ctrlnum must have 10 or fewer values \n",nm;
      exit;
    }

    if (numberof(ctrlden) > 10){
      write, format = "dm(%d).ctrlden must have 10 or fewer values \n",nm;
      exit;
    }

    if ((ctrlnum != []) && (ctrlden != [])){

      // ctrlnum and ctrlden are defined
     if (verbose){write, format = "Using dm.ctrlden and dm.ctrlnum for DM %d \n",nm;}
    } else {
      ctrlnum = [loop.gain*dm(nm).gain];
      if (gainho != []){
        grow, ctrlnum, dm(nm).gain*gainho;
      }

      ctrlden = [1,-1+loop.leak];
      if (leakho != []){
        grow, ctrlden, -1+leakho;
      }
    }
    dm(nm)._ctrlnum = &(ctrlnum);
    dm(nm)._ctrlden = &(ctrlden);
  }
}

//----------------------------------------------------
func aoall(parfile,disp=,dpi=,clean=,controlscreen=)
/* DOCUMENT aoall(parfile,disp=,dpi=,clean=,controlscreen=)
   Shorthand for aoread, aoinit, aoloop and go
   SEE ALSO: aoread, aoinit, aoloop, go
 */
{
  aoread,parfile;
  aoinit,disp=disp,dpi=dpi,clean=clean;
  aoloop,disp=disp,dpi=dpi,controlscreen=controlscreen;
  after,0.1,go;
}
//----------------------------------------------------

func aoread(parfile)
/* DOCUMENT aoread(parfile)
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

   SEE ALSO: aoinit, aoloop, go
*/
{
  extern atm,opt,sim,wfs,dm,mat,tel,target,gs,loop,parprefix,oparfile;
  extern pupil,ipupil,loopCounter,iMat,cMat;

  write,format="Yao, Version %s, %s\n",aoSimulVersion, aoSimulVersionDate;

  if (Y_VERSION == "1.5.12") {
    error,"SVD is broken in yorick 1.5.12 ! Get a later version.";
  }

  // if we are re-reading and had forks, get rid of them:
  if ((wfs!=[])&&(anyof(wfs.svipc>1))) status=quit_wfs_forks();

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
  // reset a number of other array that may have been set previously
  pupil = ipupil = loopCounter = iMat = cMat = [];

  // INIT STRUCTURES:
  atm    = atm_struct();
  opts   = opt_struct(path_type="common",scale=1.0);
  sim    = sim_struct();
  wfss   = wfs_struct(dispzoom=1,_bckgrdsub=1,shcalibseeing=0.667,subsystem=1,excessnoise=1.,framedelay=-1,pad_simage=1);
  dms    = dm_struct(gain=1.,coupling=0.2,thresholdresp=0.3);
  mat    = mat_struct(file="",fit_type="target",fit_which=1);
  tel    = tel_struct();
  target = target_struct();
  gs     = gs_struct();
  loop   = loop_struct(modalgainfile="",stats_every=4,startskip=10);

  if (!fileExist(parfile)) {
    exit,swrite(format="Cannot find parameter file %s !",parfile);}

  // read out the parfile. This stuffs values into the structures:
  include,parfile,1;

  //=====================================
  // PARAMETER CHECKS. SETS SOME DEFAULTS
  //=====================================
  check_parameters;
  wfs._origpixsize = wfs.pixsize;
  gui_message,"Next step: click on aoinit";
}

//----------------------------------------------------

func aoinit(disp=,clean=,forcemat=,svd=,dpi=,keepdmconfig=)
/* DOCUMENT aoinit(disp=,clean=,forcemat=,svd=,dpi=,keepdmconfig=)
   Second function of the ao serie.
   Initialize everything in preparation for the loop (aoloop).
   Follows a call to aoread, parfile.

   disp         = set to display stuff
   clean        = if set, aoinit will start fresh. *Nothing* is kept from
                  previous runs.
   forcemat     = set to force measuring new iMat and computing new cMat.
   svd          = set to enter SVD threshold condition by hand
   dpi          = dpi of graphic window
   keepdmconfig = when forcemat=1, the default behavior is to agregate the
                  previously determined active and extrapolated actuators
                  and do the new interaction matrix with all of them. The
                  selection of valid/extrapolated is then reset (and
                  re-executed when the interaction matrix is done).
                  Setting keepdmconfig=1 impose that the valid/extrapolated
                  setting remains as it was.
   SEE ALSO: aoread, aoloop
*/
{
  extern pupil,ipupil;
  extern iMat,cMat,fMat,dMat;
  extern iMatSP,AtAregSP,fMatSP,GxSP, polcMatSP;
  extern modalgain;
  extern _n,_n1,_n2,_p1,_p2,_p;
  extern def,mircube;
  extern tip1arcsec,tilt1arcsec;
  extern aoinit_disp,aoinit_clean,aoinit_forcemat;
  extern aoinit_svd,aoinit_keepdmconfig;
  extern tipvib, tiltvib;
  extern segmenttiptiltvib,segmentpistonvib;
  extern default_dpi;

  if (lpk_gesv == []){
    // if ylapack not in use, use LUsolve from matrix.i
    autoload, "matrix.i", LUsolve;
  } else {
    // if the ylapack package is installed, use it!
    // Run it by typing:
    // > include, "lapack"
    LUsolve = LUsolve2;
    if (sim.verbose){write, "Using ylapack to do matrix inversions: LUsolve=LUsolve2";}
  }

  disp = ( (disp==[])? (aoinit_disp==[]? 0:aoinit_disp):disp );
  clean = ( (clean==[])? (aoinit_clean==[]? 0:aoinit_clean):clean );
  forcemat = ( (forcemat==[])? (aoinit_forcemat==[]? 0:aoinit_forcemat):forcemat );
  svd = ( (svd==[])? (aoinit_svd==[]? 0:aoinit_svd):svd );
  keepdmconfig = ( (keepdmconfig==[])? (aoinit_keepdmconfig==[]? 0:aoinit_keepdmconfig):keepdmconfig );

  if (anyof(strlen(dm.actlocfile) > 0)){
    if (sim.verbose){
      write, "Using the specified actuator location file";
      write, "WARNING: running aoinit with keepdmconfig = 1";
    }
    keepdmconfig  = 1;
  }

  // initialize DM hysteresis parameters
  for (nm=1;nm<=ndm;nm++){
    dm(nm)._alpha = [0.01,0.2,4];
    dm(nm)._beta = [0.4,0.63,0.89];
    dm(nm)._w = [0.2,0.35,0.45];
  }

  if (mat.method == "mmse-sparse"){

    require,"soy.i";
    extern MR,MN;
    MR = mat.sparse_MR;
    MN = mat.sparse_MN;
    mat.file = parprefix+"-iMat.rco";

  } else {

    mat.file = parprefix+"-mat.fits";

  }

  if (sim.verbose>1) write,format="Starting aoinit with flags disp=%d,clean=%d,"+ \
                       "forcemat=%d,svd=%d,keepdmconfig=%d\n", \
                       disp,clean,forcemat,svd,keepdmconfig;

  sphase = bphase = mircube = [];

  hcp_file,YAO_SAVEPATH+parprefix+"init.ps",ps=1;

  //=====================================
  // PARAMETER CHECKS. SETS SOME DEFAULTS
  //=====================================
  check_parameters;

  wfs._origpixsize = wfs.pixsize;
  wfs._orignpixels = wfs.npixels;

  if (!is_set(disp)) {disp = 0;}
  if (!is_set(dpi)) {dpi = 70;}
  if (is_set(clean)) {forcemat=1;}

  if (!default_dpi) default_dpi=dpi;

  if (anyof(wfs.nintegcycles != 1) && (loop.method == "open-loop")) {
    exit, ">> nintegcycles > 1 not implemented for open-loop, exiting";
    }

  // Sets other parameters:
  if (!sim._size) sim._size = int(2^ceil(log(sim.pupildiam)/log(2)+1));
  size      = sim._size;
  // mircube will receive the shape of each DM (one plan per DM):
  mircube   = array(float,sim._size,sim._size,ndm);

  // INITIALIZE OUTPUT RESULT FILE YAO_SAVEPATH+parprefix.RES

  if (!fileExist(YAO_SAVEPATH+parprefix+".res")) {
    if (sim.verbose) {
      write,">> File "+parprefix+".res not found, creating one...\n";
    }
    f = open(YAO_SAVEPATH+parprefix+".res","w");
    close,f;
  }
  f = open(YAO_SAVEPATH+parprefix+".res","a+");
  write,f,format="=============================\n%s\n",sim.name;
  close,f;


  //===================================================================
  // INIT SVIPC IF NEEDED:
  //===================================================================
  if ( anyof(wfs.svipc>1) || (sim.svipc) ) {
    require,"yao_svipc.i";
    status = svipc_init();
  }


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
    if ( (wfs(i).type == "hartmann") || (wfs(i).type == "pyramid") ) {
      if (wfs(i).npixpersub) {
        npixpersub = wfs(i).npixpersub;
      } else {
        npixpersub  = float(sim.pupildiam)/wfs(i).shnxsub;
        if (npixpersub != int(npixpersub)) {
          write,"sim.pupildiam should be a multiple of wfs.shnxsub !";
          return -1;
        }
      }
      if (odd(wfs(i).shnxsub) && odd(npixpersub)) {
        grow,cent,sim._size/2+1;
      } else {
        grow,cent,sim._size/2+0.5;
      }
    } else if (wfs(i).type == "curvature") {
      grow,cent,sim._size/2+1;
    } else if (wfs(i).type == "zernike") {
      grow,cent,sim._size/2+1;
    } else if (wfs(i).type == "dh") {
      // for dh we dont mind, as they can be produced for
      // any cent.
    } else {
      if (cent == []){
        grow,cent,sim._size/2+0.5;
        write,format ="FIXME: user wfs function: assuming cent is at %.1f\n", float(sim._size/2+1);
      } else {
        grow,cent,cent(1);
        write,format = "FIXME: user wfs function: assuming cent is at %.1f\n",float(cent(1));
      }
    }
  }
  if (nallof(wfs.type=="dh") && anyof(cent != cent(1))) {
    write,"Wrong mix of hartmann/curvature or subaperture odd/even !";
    write,"I can't handle that as some of your selected sensor require";
    write,"the pupil to be centered on a pixel and others in between 2 pixels";
    write,"Sorry :-(";
    exit;
  }
  if (allof(wfs.type=="dh")) cent=sim._size/2+0.5;
  sim._cent = cent(1);

  // Initialize pupil array:

  pupil  = float(make_pupil(sim._size,sim.pupildiam,xc=sim._cent,yc=sim._cent,\
                            real=sim.pupilapod,cobs=tel.cobs)); // changed from real=1 by Marcos.
  ipupil = float(make_pupil(sim._size,sim.pupildiam,xc=sim._cent,yc=sim._cent,\
                          cobs=tel.cobs));

  if (user_pupil) user_pupil;

  pupil = float(pupil);
  ipupil = float(ipupil);

  // if the pupil for the WFSs has not already been defined by the user_pupil function, then set it to be the same as the telescope pupil
  for (ns=1;ns<=numberof(wfs);ns++){
    if (*wfs(ns)._pupil == []){wfs(ns)._pupil = &pupil;}
  }

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
    if (!yaopy) {
      // set if using pygtk GUI, which prevents remapping a new window
      status = create_yao_window();
    }
    if (wfs_display_mode=="spatial") {
      disp2d,wfs._dispimage,wfs.pupoffset(1,),wfs.pupoffset(2,),2,\
                    zoom=wfs.dispzoom,init=1;
    } else {
      disp2d,wfs._dispimage,wfs.gspos(1,),wfs.gspos(2,),2,zoom=wfs.dispzoom,init=1;
    }
  }

  //==================================
  // INITIALIZE PHASE SCREENS
  //==================================

  if (sim.verbose) {write,"\n> INITIALIZING PHASE SCREENS";}
  gui_message,"Initializing phase screens";
  get_turb_phase_init;

  //==================================
  // INITIALIZE SENSOR
  //==================================

  if (sim.verbose) {write,"\n> INITIALIZING SENSOR GEOMETRY";}
  gui_message,"Initializing sensor geometry";

  for (ns=1;ns<=nwfs;ns++) {

    if (wfs(ns).type == "curvature") {

      // build subaperture geometry:
      make_curv_wfs_subs,ns,size,sim.pupildiam;
      // init WFS
      curv_wfs,*wfs(ns)._pupil,*wfs(ns)._pupil*0.0f,ns,init=1;
      wfs(ns)._nsub = int(sum(*(wfs(ns).nsubperring)));
      wfs(ns)._nmes = wfs(ns)._nsub;

    } else if (wfs(ns).type == "hartmann") {

      // init WFS:
      if (wfs(ns).disjointpup) {
        shwfs_init,disjointpup(,,ns),ns,imat=1,clean=clean;
      } else shwfs_init,*wfs(ns)._pupil,ns,imat=1,clean=clean;
      wfs(ns)._nmes = 2*wfs(ns)._nsub;

    } else if (wfs(ns).type == "pyramid") {

      // init WFS
      v = pyramid_wfs(*wfs(ns)._pupil,*wfs(ns)._pupil*0.,ns,disp=disp,init=1);
      wfs(ns)._nsub = numberof(v)/2;
      wfs(ns)._nmes = 2*wfs(ns)._nsub;

    } else if (wfs(ns).type == "zernike") {
      zernike_wfs,*wfs(ns)._pupil,*wfs(ns)._pupil*0.,ns,init=1;

    } else if (wfs(ns).type == "dh") {
      dh_wfs,*wfs(ns)._pupil,*wfs(ns)._pupil*0.,ns,init=1;

    } else {
      // assign user_wfs to requested function/type:
      cmd = swrite(format="user_wfs = %s",wfs(ns).type);
      include,[cmd],1;
      user_wfs,*wfs(ns)._pupil,*wfs(ns)._pupil*0.,ns,init=1;
    }

    if ( (wfs(ns).disjointpup) && (disjointpup==[]) ) \
      error,swrite(format="wfs(%d).disjointpup set but disjointpup does not exist\n",ns);
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

  // save state of noise/nintegcycle/etc: everything that is not desired
  // when doing the iMat:
  store_noise_etc_for_imat,noise_orig, cycle_orig, kconv_orig, \
               skyfluxpersub_orig, bckgrdcalib_orig, bias_orig, \
               flat_orig, darkcurrent_orig, use_sincos_approx_orig,\
               rayleigh_orig;

  // sync forks if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  refmes = mult_wfs_int_mat(disp=disp);
  wfs._refmes = split_wfs_vector(refmes);

  wfs.filtertilt = mem;

  //============================================
  // GET WFS TIP AND TILT REFERENCE MEASUREMENTS
  //============================================

  // disable filtertilt is any
  mem = wfs.filtertilt; wfs.filtertilt *= 0n;

  // step per pixel to have a x" tilt:
  if (!push4ttref) push4ttref=0.1; // in arcsec

  // tip:
  mircube *= 0.0f;
  mircube(,,1) = float(tip1arcsec*push4ttref);

  // we've changed some wfs value. sync svipc if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  mes = mult_wfs_int_mat(disp=disp)/push4ttref;  // mes in arcsec

  wfs._tiprefv = split_wfs_vector(mes);

  for (ns=1;ns<=nwfs;ns++) {
    sumtip2 = sum(*wfs(ns)._tiprefv^2.);
    if (sumtip2 > 0.){
      wfs(ns)._tiprefvn = &(*wfs(ns)._tiprefv/sumtip2);
    } else {
      wfs(ns)._tiprefvn = &(*wfs(ns)._tiprefv);
    }
  }

  // now the tip ref vector are normalized, so to compute the tip content
  // one has just to do sum(vector * tiprefv)

  // tilt:
  mircube(,,1) = float(tilt1arcsec*push4ttref);

  // we've changed some wfs value. sync svipc if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  mes = mult_wfs_int_mat(disp=disp)/push4ttref;  // mes in arcsec

  wfs._tiltrefv = split_wfs_vector(mes);

  for (ns=1;ns<=nwfs;ns++) {
    sumtilt2 = sum(*wfs(ns)._tiltrefv^2.);
    if (sumtilt2 > 0.){
      wfs(ns)._tiltrefvn = &(*wfs(ns)._tiltrefv/sumtilt2);
    } else {
      wfs(ns)._tiltrefvn = &(*wfs(ns)._tiltrefv);
    }
  }

  // restore pre-operation filtertilt:
  wfs.filtertilt = mem;

  // restore original values to WFS structure:
  restore_noise_etc_for_imat,noise_orig, cycle_orig, kconv_orig, \
                 skyfluxpersub_orig,bckgrdcalib_orig, bias_orig, \
                 flat_orig, darkcurrent_orig, use_sincos_approx_orig,\
                 rayleigh_orig;

  // sync forks if needed:
  if ( (anyof(wfs.type=="hartmann"))&&(anyof(wfs.svipc>1))) s = sync_wfs_forks();

  //==================================
  // INITIALIZE DM INFLUENCE FUNCTIONS
  //==================================
  if (sim.verbose) {write,"\n> Initializing DM influence functions";}
  gui_message,"Initializing DMs";

  // loop over DMs:
  for (n=1;n<=ndm;n++) {
    if ( (dm(n).disjointpup) && (disjointpup==[]) ) \
    error,swrite(format="dm(%d).disjointpup set but disjointpup does not exist\n",n);

    if (dm(n).pupoffset!=[]) \
       dm(n)._puppixoffset = long(dm(n).pupoffset/tel.diam*sim.pupildiam);

    if (clean) dm(n)._def = dm(n)._edef = &([]);
    // Set _n1 and _n2, the limit indices
    if (dm(n).type == "stackarray"){
      if (dm(n).alt == 0) {
        // special case = (only) 6 pixels margin each side:
        // note: "upgraded" from 2 to 8 to allow more margin when misregistering
        extent = sim.pupildiam+16;
      } else {
        extent = dm(n).pitch*(dm(n).nxact+2.); // + 1.5 pitch each side
      }
      dm(n)._n1 = long(clip(floor(sim._cent-extent/2.),1,));
      dm(n)._n2 = long(clip(ceil(sim._cent+extent/2.),,sim._size));
    } else {
      // define for all other programmed DMs, but not for user defined DMs, which could have different values
      if (anyof(dm(n).type == ["bimorph","zernike","dh","kl","tiptilt","segmented","aniso"])){
        dm(n)._n1 = 1;
        dm(n)._n2 = sim._size;
      }
    }

    // compute influence functions:
    // If file exist, read it out:
    if ( (fileExist(YAO_SAVEPATH+dm(n).iffile)) && (!is_set(clean)) ) {

      if (sim.verbose) {
        write,format=">> Reading file %s\n",dm(n).iffile;
      }
      if (dm(n).use_def_of) {
        write,format="Replicating influence functions of DM%d\n",dm(n).use_def_of;
        dm(n)._def = dm(dm(n).use_def_of)._def;
      } else {
        dm(n)._def = &(float(yao_fitsread(YAO_SAVEPATH+dm(n).iffile)));
      }
      dm(n)._nact = dimsof(*(dm(n)._def))(4);
      if ( dm(n).type == "stackarray" ) {
        dm(n)._x = &(yao_fitsread(YAO_SAVEPATH+dm(n).iffile,hdu=1));
        dm(n)._y = &(yao_fitsread(YAO_SAVEPATH+dm(n).iffile,hdu=2));
        if (dm(n).elt == 1) {
          dm(n)._eltdefsize = dimsof(*(dm(n)._def))(2);
          dm(n)._i1 = &(int(yao_fitsread(YAO_SAVEPATH+dm(n).iffile,hdu=3)));
          dm(n)._j1 = &(int(yao_fitsread(YAO_SAVEPATH+dm(n).iffile,hdu=4)));
        }
      }

      if ( (fileExist(YAO_SAVEPATH+dm(n)._eiffile)) && (!is_set(clean)) ) {
        if (sim.verbose) {
          write,format=">> Reading extrapolated actuators file %s\n",dm(n)._eiffile;
        }
        dm(n)._edef = &(float(yao_fitsread(YAO_SAVEPATH+dm(n)._eiffile)));
        dm(n)._enact = dimsof(*(dm(n)._edef))(4);
        if ( dm(n).type == "stackarray" ) {
          dm(n)._ex = &(yao_fitsread(YAO_SAVEPATH+dm(n)._eiffile,hdu=1));
          dm(n)._ey = &(yao_fitsread(YAO_SAVEPATH+dm(n)._eiffile,hdu=2));
          if (dm(n).elt == 1) {
            dm(n)._ei1 = &(int(yao_fitsread(YAO_SAVEPATH+dm(n)._eiffile,hdu=3)));
            dm(n)._ej1 = &(int(yao_fitsread(YAO_SAVEPATH+dm(n)._eiffile,hdu=4)));
          }
        }
      }

    } else { // else compute the influence functions:

      if (sim.verbose) {
        write,format=">> Computing Influence functions for mirror # %d\n",n;
      }

      if (fileExist(YAO_SAVEPATH+dm(n)._eiffile)) {// delete the extrapolated influence functions
        remove, YAO_SAVEPATH+dm(n)._eiffile;
      }

      dispif = (disp&&(sim.debug));
      if (dispif) { plsys,1; animate,1; }

      if (dm(n).use_def_of) {

        write,format="Replicating influence functions of DM%d\n",dm(n).use_def_of;
        dm(n)._def = dm(dm(n).use_def_of)._def;
        dm(n)._nact = dm(dm(n).use_def_of)._nact;

      } else {

        if (dm(n).type == "bimorph") {
          make_curvature_dm, n, disp=dispif,cobs=tel.cobs;
        } else if (dm(n).type == "stackarray") {
          if (dm(n).elt == 1) {
            make_pzt_dm_elt, n, disp=dispif;
          } else {
            make_pzt_dm, n, disp=dispif;
          }
        } else if (dm(n).type == "zernike") {
          make_zernike_dm, n, disp=dispif;
        } else if (dm(n).type == "dh") {
          make_dh_dm, n, disp=dispif;
        } else if (dm(n).type == "kl") {
          make_kl_dm, n, disp=dispif;
        } else if (dm(n).type == "tiptilt") {
          make_tiptilt_dm, n, disp=dispif;
        } else if (dm(n).type == "segmented") {
          make_segmented_dm, n, disp=dispif;
        } else if (dm(n).type == "aniso") {
          make_aniso_dm, n, disp=dispif;
        } else {
          // we're dealing with a user defined DM function:
          // assign user_wfs to requested function/type:
          cmd = swrite(format="user_dm = %s",dm(n).type);
          include,[cmd],1;
          user_dm, n;
        }
      }

      if (dm(n).ifunrot) {
        hxy = dimsof(*dm(n)._def)(2)/2.+0.5;
        for (i=1;i<=dm(n)._nact;i++) {
          (*dm(n)._def)(,,i) = rotate2((*dm(n)._def)(,,i),dm(n).ifunrot,xc=hxy,yc=hxy);
        }
        xy = (*dm(n)._x)(,-);
        grow,xy,(*dm(n)._y)(,-);
        xy -= sim._cent;
        xy = transpose(xy);
        xy = mrot(dm(n).ifunrot)(+,) * xy(+,);
        xy += sim._cent;
        (*dm(n)._x) = xy(1,);
        (*dm(n)._y) = xy(2,);
      }

      if (dm(n).xscale) {
        dd = dimsof(*dm(n)._def)(2);
        xx = yy = indgen(dd);
        xx = (xx-dd/2.)*(1.+dm(n).xscale)+dd/2.;
        for (i=1;i<=dm(n)._nact;i++) (*dm(n)._def)(,,i) = bilinear((*dm(n)._def)(,,i),xx,yy,grid=1);
        *dm(n)._x = (*dm(n)._x-sim._cent)*(1-dm(n).xscale)+sim._cent;
      }

      if (dispif) { plsys,1; animate,0; }

      // the IF are in microns/volt
      if (sim.verbose) {
        write,format="\n>> Storing influence functions in %s...",dm(n).iffile;
      }
      if (dm(n).use_def_of) yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,[0.];
      else yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,*(dm(n)._def);
      if (sim.verbose) write,"Done";
      if ( dm(n).type == "stackarray" ) {
        yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,*(dm(n)._x),exttype="IMAGE",append=1;
        yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,*(dm(n)._y),exttype="IMAGE",append=1;
        if (dm(n).elt == 1) {
          yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,long(*(dm(n)._i1)),exttype="IMAGE",append=1;
          yao_fitswrite,YAO_SAVEPATH+dm(n).iffile,long(*(dm(n)._j1)),exttype="IMAGE",append=1;
        }
      }
    }
    // need to set ._n1 and ._n2 for user defined DMs that have already been defined and saved
    if (dm(n)._n1 == 0){
      dmdims = dimsof(*dm(n)._def);
      dm(n)._n1 = long(sim._cent+0.5-dmdims(2)/2);
      dm(n)._n2 = long(sim._cent-0.5+dmdims(2)/2);
    }
  }

  //==============================
  // INITIALIZE VIBRATION SPECTRUM
  //==============================

  if ((tel.tipvib_white_rms) || (tel.tipvib_1overf_rms) || ((*tel.tipvib_peaks)!=[])) {
    // user has defined at least one of the vibration param.
    tipvib = generate_vib_time_serie(loop.ittime,loop.niter,tel.tipvib_white_rms,\
      tel.tipvib_1overf_rms, *tel.tipvib_peaks, *tel.tipvib_peaks_rms, \
      peak_width = *tel.tipvib_peaks_width);
  } else {
    tipvib = [];
  }
  if ((tel.tiltvib_white_rms) || (tel.tiltvib_1overf_rms) || ((*tel.tiltvib_peaks)!=[])) {
    // user has defined at least one of the vibration param.
    tiltvib = generate_vib_time_serie(loop.ittime,loop.niter,tel.tiltvib_white_rms,\
      tel.tiltvib_1overf_rms, *tel.tiltvib_peaks, *tel.tiltvib_peaks_rms, \
      peak_width = *tel.tiltvib_peaks_width);
  } else {
    tiltvib = [];
  }


  //=========================================
  // DO INTERACTION MATRIX WITH ALL ACTUATORS
  //=========================================

  // determine whether a new interaction matrix is needed

  need_new_iMat = (forcemat || anyof(!fileExist(YAO_SAVEPATH+dm.iffile)));

  if (need_new_iMat == 0){
    if (!fileExist(YAO_SAVEPATH+mat.file)){
      if (mat.method == "mmse-sparse" && fileExist(YAO_SAVEPATH + parprefix+"-mat.fits")){ // convert the full matrix into a sparse matrix
        write, "Saving " + parprefix + "-mat.fits" + " as a sparse matrix";
        tmp = yao_fitsread(YAO_SAVEPATH+ parprefix + "-mat.fits");
        iMat = tmp(,,1);
        iMatSP = sprco(float(iMat));
        save_rco,iMatSP,YAO_SAVEPATH+mat.file;
        svd = 1; // need to generate a new reconstructor
      } else if (mat.method != "mmse-sparse" && fileExist(YAO_SAVEPATH + parprefix+"-iMat.rco")){
        write, "Saving " + parprefix + "-iMat.rco" + " as a full matrix";
        iMatSP = restore_rco(YAO_SAVEPATH + parprefix+"-iMat.rco");
        iMat = rcoinf(iMatSP);
        yao_fitswrite, YAO_SAVEPATH + mat.file, [iMat,iMat];
        svd = 1; // need to generate a new reconstructor
      } else {
        need_new_iMat = 1;
      }
    }
  }

  if (need_new_iMat == 1){
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

    if (sim.verbose) {write,"\n> DOING INTERACTION MATRIX";}
    gui_message,"Doing interaction matrix";
    if (mat.method != "mmse-sparse"){
      iMat = array(double,sum(wfs._nmes),sum(dm._nact));
      cMat = array(double,sum(dm._nact),sum(wfs._nmes));

      if (sim.verbose>1) write,format="sizeof(iMat) = %dMB\n",long(sizeof(iMat)/1.e6);
      pause,100;
    } else {
      iMatSP = rco();
    }

    // measure interaction matrix:
    estDMs = where(!dm.dmfit_which); // DMs used to estimate wavefront
    if ( sum(dm(estDMs)._nact) > sum(wfs._nmes) ) {
      write,format="\n\nWarning: Underconstrained problem: nact (%d) > nmes (%d)\n",sum(dm._nact),sum(wfs._nmes);
      if (mat.method == "svd"){
        write,format="%s\n\n","I will not be able to invert the iMat using"+\
                     " the simple yao SVD";
        typeReturn;
      }
    }

    if (use_user_imat) iMat = user_imat;
    else do_imat,disp=disp;

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
      if (mat.method == "mmse-sparse"){
        AtA = rcoata(iMatSP);
        resp = sqrt((*AtA.xd)(1:AtA.r)); // take the diagonal
        actuators_to_remove = [];
      } else {
        resp = sqrt((iMat^2.)(sum,));
        //resp = (abs(iMat))(max,); // marcos commented this. problem?
      }

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
              graphic_config,dm(nm).subsystem,nm;
              write,format="DM #%d: # of valid actuators: %d ",nm,numberof(ok);
              again = "y"; read,prompt="Again ? ",again;
              if (again != "n") {
                read,prompt=swrite(format="Threshold (old = %f) : ",\
                                   dm(nm).thresholdresp),dm(nm).thresholdresp;
              }
            } while (again != "n");
          }

          if (dm(nm).ncp){// if DM used for WFS
            ok = indgen(numberof(tmp));
            nok = [];
          } else {
            ok = where(tmp >  dm(nm).thresholdresp*max(tmp));
            nok= where(tmp <= dm(nm).thresholdresp*max(tmp));
            if (numberof(nok)==0) nok=[];
          }

          dm(nm)._indval = &ok;
          dm(nm)._indext = &nok;

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

          if (sim.verbose) {
            write,format="DM #%d: # of valid actuators: %d. "+
              "(I got rid of %d actuators after iMat)\n",
              nm,numberof(ok),numberof(nok);
            write,format=">> valid I.F. stored in %s\n",dm(nm).iffile;
          }

          // rewrite influence function file:
          yao_fitswrite,YAO_SAVEPATH+dm(nm).iffile,*(dm(nm)._def);
          yao_fitswrite,YAO_SAVEPATH+dm(nm).iffile,*(dm(nm)._x),exttype="IMAGE",append=1;
          yao_fitswrite,YAO_SAVEPATH+dm(nm).iffile,*(dm(nm)._y),exttype="IMAGE",append=1;
          if (dm(nm).elt == 1) {
            yao_fitswrite,YAO_SAVEPATH+dm(nm).iffile,long(*(dm(nm)._i1)),exttype="IMAGE",append=1;
            yao_fitswrite,YAO_SAVEPATH+dm(nm).iffile,long(*(dm(nm)._j1)),exttype="IMAGE",append=1;
          }
          dm(nm)._nact = (dimsof(*(dm(nm)._def)))(4);

          // write extrapolated actuator influence functions file:
          yao_fitswrite,YAO_SAVEPATH+dm(nm)._eiffile,*(dm(nm)._edef);
          yao_fitswrite,YAO_SAVEPATH+dm(nm)._eiffile,*(dm(nm)._ex),exttype="IMAGE",append=1;
          yao_fitswrite,YAO_SAVEPATH+dm(nm)._eiffile,*(dm(nm)._ey),exttype="IMAGE",append=1;
          if (dm(nm).elt == 1) {
            yao_fitswrite,YAO_SAVEPATH+dm(nm)._eiffile,long(*(dm(nm)._ei1)),exttype="IMAGE",append=1;
            yao_fitswrite,YAO_SAVEPATH+dm(nm)._eiffile,long(*(dm(nm)._ej1)),exttype="IMAGE",append=1;
          }
          dm(nm)._enact = (dimsof(*(dm(nm)._edef)))(4);

          if (sim.debug >= 1) {
            tv,comp_dm_shape(nm,&(array(1.0f,dm(nm)._nact)));
            mypltitle,swrite(format="Pushing on all actuator of DM#%d",nm),[0.,-0.005],height=12;
            typeReturn;
          }

          if (mat.method == "mmse-sparse"){
            grow, actuators_to_remove, indexDm(1,nm)-1+nok;
          } else {
            iMat(,indexDm(1,nm)-1+nok) *=0.;
          }

        }
      }
      if (mat.method == "mmse-sparse"){
        // remove all the rows with a response that is too low
        // need to remove them from the end so as not to disturb numbering
        if (numberof(actuators_to_remove) > 0){
          for (act_ii = numberof(actuators_to_remove); act_ii >=1; act_ii--){
            rcodr,iMatSP,actuators_to_remove(act_ii);
          }
        }
      } else {
        iMat = iMat(,where(iMat(rms,) != 0.));
        cMat = transpose(iMat)*0.;
      }
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

    if (dm(nm).type != "stackarray") continue;  // possible but not implemented.

    if (dm(nm).noextrap == 1) { // we have to get rid of any extrapolated actuator now
      dm(nm)._edef = dm(nm)._ex = dm(nm)._ey = &([]);
      if (dm(nm).elt == 1) {
        dm(nm)._ei1 = dm(nm)._ej1 = &([]);
      }
      dm(nm)._enact = 0;
      continue;
    }

    if (dm(nm)._enact == 0) continue;  // no extrapolated actuator

    if (fileExist(YAO_SAVEPATH+dm(nm).ecmatfile)) {

      // if defined and exist, we read it and proceed
      if (sim.verbose) {
        write,format="Reading valid to extrap. matrix %s\n",dm(nm).ecmatfile;
      }
      dm(nm)._extrapcmat = yao_fitsread(YAO_SAVEPATH+dm(nm).ecmatfile);

    } else {

      // we compute the valid to extrap matrix:
      if (sim.verbose) {
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
     WfsRefMes = sh_wfs(ipupil,phase);
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

  if (mat.method == "svd"){
    if (sim.verbose>1) {write,"\n> INITIALIZING MODAL GAINS";}
    gui_message,"Initializing modal gains";

    modalgain = [];

    if (strlen(loop.modalgainfile)&&(fileExist(YAO_SAVEPATH+loop.modalgainfile))) {

      if (sim.verbose) {write,format=">> Reading file %s\n\n",loop.modalgainfile;}
      modalgain = yao_fitsread(YAO_SAVEPATH+loop.modalgainfile);

    }

    if (numberof(modalgain) != sum(dm._nact)) {

      if (sim.verbose>1) {
        write,format="I did not find %s or it did not have the right\n ",
          loop.modalgainfile;
        write,format="  number of elements (should be %i), so I have generated\n",
          sum(dm._nact);
        write,"  a modal gain vector with ones in it.";
      }

      modalgain = array(1.,sum(dm._nact));

    }
  }

  // INITIALIZE COMMAND MATRIX:

  if (sim.verbose) {write,"\n> INTERACTION AND COMMAND MATRICES";}
  gui_message,"Initializing control matrix";

  if ((fileExist(YAO_SAVEPATH+mat.file)) && (forcemat != 1)) {

    if (sim.verbose) {write,format=">> Reading file %s\n",mat.file;}
    if (mat.method != "mmse-sparse") {
      // read out mat file and stuff iMat and cMat:
      tmp = yao_fitsread(YAO_SAVEPATH+mat.file);
      iMat = tmp(,,1);
      cMat = transpose(tmp(,,2));
      tmp = [];

      if (anyof(dm.dmfit_which)){
        if (fileExist(YAO_SAVEPATH+parprefix+"-dMat.fits")){
          dMat = yao_fitsread(YAO_SAVEPATH + parprefix + "-dMat.fits");
        } else {
          svd = 1;
        }
        if (fileExist(YAO_SAVEPATH+parprefix+"-cMat.fits")){
          cMat = yao_fitsread(YAO_SAVEPATH + parprefix + "-cMat.fits");
        } else {
          svd = 1;
        }
        for (nm=1;nm<=numberof(dm);nm++){
          if (dm(nm).dmfit_which){
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".fits");
            if (fileExist(filename)){
              dm(nm)._fMat = &yao_fitsread(filename);
            } else {
              svd = 1;
            }
          }
        }
      }
    } else {
      dMat = [];
      iMatSP = restore_rco(YAO_SAVEPATH+mat.file);
      if (fileExist(YAO_SAVEPATH+parprefix+"-AtAreg.ruo")){
        AtAregSP = restore_ruo(YAO_SAVEPATH+parprefix+"-AtAreg.ruo");
      } else {
        svd = 1;
      }
      if (fileExist(YAO_SAVEPATH+parprefix+"-GxSP.rco")){
        GxSP = restore_rco(YAO_SAVEPATH+parprefix+"-GxSP.rco");
      } else {
        svd = 1;  // need to recreate reconstructors
      }
      if (anyof(dm.dmfit_which)) {
        if (fileExist(YAO_SAVEPATH+parprefix+"-polcMat.rco")){
          polcMatSP = restore_rco(YAO_SAVEPATH+parprefix+"-polcMat.rco");
        } else {
          svd = 1;  // need to recreate reconstructors
        }
        for (nm=1;nm<=numberof(dm);nm++){
          if (dm(nm).dmfit_which){
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".rco");
            if (fileExist(filename)){
              dm(nm)._fMat = &restore_rco(filename);
            } else {
              svd = 1;  // need to recreate reconstructors
            }
          }
        }
      }
    }
  }

  if (svd || forcemat || need_new_iMat){
    // create the regularization matrices
    nDMs = numberof(dm);
    for (nm=1;nm<=nDMs;nm++){
      if (*dm(nm).regmatrix != []){
        dm(nm)._regmatrix = dm(nm).regmatrix;
      } else if (mat.method == "mmse"){
        if ((dm(nm).type == "stackarray") && (dm(nm).regtype == "laplacian")){
          xpos = *dm(nm)._x;
          ypos = *dm(nm)._y;
          pitch2 = dm(nm).pitch^2;

          L = array(float,dm(nm)._nact, dm(nm)._nact);

          for (ii=1;ii<=dm(nm)._nact;ii++){

            xii = xpos(ii);
            yii = ypos(ii);

            dx = xpos - xii;
            dy = ypos - yii;

            dist2= dx^2+dy^2;

            L(ii,ii) = -1.;

            ww = where(dist2 == pitch2);
            if (numberof(ww)){L(ii,ww) = 0.25;}
          }
          dm(nm)._regmatrix = &(L(+,)*L(+,));
        } else { // identity matrix
          dm(nm)._regmatrix = &unit(dm(nm)._nact);
        }
      } else if (mat.method == "mmse-sparse") {
        if ((dm(nm).type == "stackarray") && (dm(nm).regtype == "laplacian")){
          xpos = *dm(nm)._x;
          ypos = *dm(nm)._y;
          pitch2 = dm(nm).pitch^2;

          laplacian_mat = rco();

          for (ii=1;ii<=dm(nm)._nact;ii++) {
            xii = xpos(ii);
            yii = ypos(ii);

            dx = xpos - xii;
            dy = ypos - yii;

            dist2= dx^2+dy^2;

            lvec = array(float,dm(nm)._nact);
            lvec(ii) = -1.;
            ww = where(dist2 == pitch2);
            if (numberof(ww)) lvec(ww) = 0.25;

            rcobuild,laplacian_mat,lvec,mat.sparse_thresh;
          }
          dm(nm)._regmatrix = &rcoata(laplacian_mat);
        } else { // identity matrix
          dm(nm)._regmatrix = &spunit(dm(nm)._nact);
      }
    }
  }

  // create the fitting matrices for tomography
    if (mat.fit_simple == 1){
      for (nm=1;nm<=numberof(dm);nm++){
        if (dm(nm).dmfit_which){

          virtualDMs = int(*dm(nm).dmfit_which);
          nVirtualDMs = numberof(virtualDMs);

          if (dm(nm).type == "stackarray"){
            xloct = *dm(nm)._x; // location of the tomography actuators
            yloct = *dm(nm)._y;

            xlocv = ylocv = []; // location of the virtual actuators

            for (c1=1;c1<=nVirtualDMs;c1++){
              grow, xlocv, *dm(virtualDMs(c1))._x;
              grow, ylocv, *dm(virtualDMs(c1))._y;
            }
          } else if (dm(nm).type == "zernike"){
            xloct = indgen(dm(nm).minzer:dm(nm).nzer);
            yloct = indgen(dm(nm).minzer:dm(nm).nzer);

            xlocv = ylocv = []; // location of the virtual actuators

            for (c1=1;c1<=nVirtualDMs;c1++){
              grow, xlocv, indgen(dm(virtualDMs(c1)).minzer:dm(virtualDMs(c1)).nzer);
              grow, ylocv, indgen(dm(virtualDMs(c1)).minzer:dm(virtualDMs(c1)).nzer);
            }
          }

          if (mat.method == "mmse"){
            dm(nm)._fMat = &array(float,[2,dm(nm)._nact,sum(dm(virtualDMs)._nact)]);
            for (c1=1;c1<=numberof(xlocv);c1++){
              v1 = ((xloct == xlocv(c1)) + (yloct == ylocv(c1)) == 2);
              (*dm(nm)._fMat)(,c1) = float(v1);
            }
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".fits");
            yao_fitswrite,filename,*dm(nm)._fMat;
          } else {
            temp = rco_d();
            for (c1=1;c1<=numberof(xlocv);c1++){
              v1 = ((xloct == xlocv(c1)) + (yloct == ylocv(c1)) == 2);
              rcobuild, temp, float(v1), mat.sparse_thresh;
            }
            dm(nm)._fMat = &rcotr(temp);
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".rco");
            save_rco, *dm(nm)._fMat, filename;
          }
        }
      }
    } else { // mat.fit_simple = 0

      // define a grid on which to sample the phase
      dims = (dimsof(pupil))(2);
      idx =  indgen(int(ceil(mat.fit_subsamp/2.)):dims:mat.fit_subsamp);

      pupil_sub = pupil(idx,idx);
      wpupil_sub = where(pupil_sub);
      npix = numberof(wpupil_sub);

      for (nm=1;nm<=numberof(dm);nm++){
        if (dm(nm).dmfit_which){

          virtualDMs = int(*dm(nm).dmfit_which);
          nVirtualDMs = numberof(virtualDMs);

          // calculate the matrix that shows how the tomographic DM affects the phase
          n1 = dm(nm)._n1;
          n2 = dm(nm)._n2;

          if (mat.method == "mmse"){
            tomoMat =  array(float,[2,npix,dm(nm)._nact]);
          } else {
            tomoMatSP = rco();
          }

          command = array(float,dm(nm)._nact);
          for (i=1;i<=dm(nm)._nact;i++) {
            if (sim.verbose) {
              write,format="\rFitting DM# %d, actuator %d/%d",nm,i,dm(nm)._nact;
            }
            mircube *= 0.0f;
            command *= 0.0f;
            command(i) = float(1.);
            mircube(n1:n2,n1:n2,nm) = comp_dm_shape(nm,&command);
            if (dm(nm).ncp){ // DM on WFS path only
              phase = get_phase2d_from_dms(dm(nm).ncpfit_which,dm(nm).ncpfit_type)-get_phase2d_from_dms(mat.fit_which,mat.fit_type);
            } else {
              phase = get_phase2d_from_dms(mat.fit_which,mat.fit_type);
            }

            phase = phase(idx,idx);
            active_phase =  phase(wpupil_sub);
            if (mat.method == "mmse"){
              tomoMat(:,i) = active_phase;
            } else {
              rcobuild, tomoMatSP, float(active_phase), mat.sparse_thresh;
            }
          }

          // calculate the matrix that shows how the virtual DMs affect the phase
          if (mat.method == "mmse"){
            virtMat = array(float,[2,npix,sum(dm(virtualDMs)._nact)]);
          } else {
            regmatrix = *dm(nm)._regmatrix;
            ruox, regmatrix, dm(nm).regparam;
            tomoFit = ruoadd(rcoata(tomoMatSP),regmatrix);
            regmatrix = [];
            dm_fMatSP = rco();
          }
          row = 0;
          zeros_act = array(float,dm(nm)._nact);

          for (k=1;k<=nVirtualDMs;k++) {
            nv = virtualDMs(k);

            n1 = dm(nv)._n1;
            n2 = dm(nv)._n2;

            command = array(float,dm(nv)._nact);

            for (i=1;i<=dm(nv)._nact;i++) {
              if (sim.verbose) {
                write,format="\rFitting DM# %d, actuator %d/%d",nv,i,dm(nv)._nact;
              }
              row += 1;
              mircube  *= 0.0f;
              command *= 0.0f;
              command(i) = 1.;
              mircube(n1:n2,n1:n2,nv) = comp_dm_shape(nv,&command);
              if (dm(nm).ncp){ // DM on WFS path only
                phase = get_phase2d_from_dms(dm(nm).ncpfit_which,dm(nm).ncpfit_type)-get_phase2d_from_dms(mat.fit_which,mat.fit_type);
              } else {
                phase = get_phase2d_from_dms(mat.fit_which,mat.fit_type);
              }

              phase = phase(idx,idx);
              active_phase =  phase(wpupil_sub);
              if (mat.method == "mmse"){
                virtMat(:,row) = active_phase;
              } else {
                temp = rcoxv(tomoMatSP,active_phase);
                if (sum(temp != 0) == 0){ // sparse CG method does not work if all zeros
                  temp = zeros_act;
                } else {
                  temp =  ruopcg(tomoFit,temp,zeros_act,tol=mat.sparse_pcgtol);
                }
                rcobuild, dm_fMatSP, temp, mat.fit_minval;
              }
            }
          }

          if (mat.method == "mmse"){
            dm(nm)._fMat = &(LUsolve(tomoMat(+,)*tomoMat(+,) + dm(nm).regparam* (*dm(nm)._regmatrix),tomoMat(+,)*virtMat(+,)));
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".fits");
            yao_fitswrite,filename,*dm(nm)._fMat;
          } else {
            tomoMatSP = tomoFit = virtMatSP = [];
            dm(nm)._fMat = &rcotr(dm_fMatSP);
            filename = YAO_SAVEPATH+parprefix+"-fMat"+swrite(nm, format="%i"+".rco");
            save_rco, *dm(nm)._fMat, filename;
            dm_fMatSP = [];
          }
        }
      }
    }
  }

  // in the opposite case, plus if svd=1 (request re-do SVD):
  if ((!fileExist(YAO_SAVEPATH+mat.file)) || (forcemat == 1) || (svd == 1)) {
    if (mat.method == "svd") {
      if (sim.verbose) {
        write,">> Preparing SVD and computing command matrices";
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

        prep_svd,imat,ss,svd=svd,disp=disp;

        tmpdisp = noneof(dm(where(dm.subsystem == ss)).type == "aniso");
        cmat = build_cmat( (*mat.condition)(ss), mg, subsystem=ss,
                          all=1,disp=tmpdisp);
        cMat(wssdm,wsswfs) = cmat;
      }

    } else if (mat.method == "mmse") {
      if (anyof(dm.dmfit_which)){
        estAct = []; // actuators used to estimate wavefront
        realAct = []; // real (non virtual) actuators used to compensate the wavefront
        actno = 0;

        for (nm=1;nm<=ndm;nm++){
          if (!dm(nm).dmfit_which){
            grow, estAct, indgen(1+actno:actno+dm(nm)._nact);
          }
          if (!dm(nm).virtual & !dm(nm).ncp){
            grow, realAct, indgen(1+actno:actno+dm(nm)._nact);
          }
          actno += dm(nm)._nact;
        }

        Ga = iMat(:,realAct);
        Gx = iMat(:,estAct);
        AtA = Gx(+,)*Gx(+,);

        nEstAct = numberof(estAct);
        nRealAct = numberof(realAct);

        Cphi = array(float,[2,nEstAct,nEstAct]);

        mc = 0;
        for (nm=1;nm<=nDMs;nm++){
          if (!dm(nm).dmfit_which){
            Cphi(mc+1:mc+dm(nm)._nact,mc+1:mc+dm(nm)._nact) = \
                          dm(nm).regparam*(*dm(nm)._regmatrix);
            mc += dm(nm)._nact;
          }
        }

        cMat = LUsolve(AtA+Cphi,transpose(Gx));
        yao_fitswrite, YAO_SAVEPATH + parprefix + "-cMat.fits", cMat;

        realDMs = where(!dm.virtual & !dm.ncp); // DMs used to compensate wavefront
        estDMs = where(!dm.dmfit_which); // DMs used to estimate wavefront

        nRealAct = sum(dm(realDMs)._nact);
        nEstAct = sum(dm(estDMs)._nact);

        fMat = array(float,[2,nRealAct, nEstAct]);

        indexDm       = array(long,2,ndm);
        indexDm(,1)   = [1,dm(1)._nact];
        for (nm=2;nm<=ndm;nm++) {
          indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
        }

        startIdx = 1;
        for (nm=1;nm<=nDMs;nm++){
          if (!dm(nm).virtual & !dm(nm).ncp){
            idx = indgen(startIdx:startIdx + dm(nm)._nact - 1);
            startIdx +=  dm(nm)._nact;

            if (!dm(nm).dmfit_which) {
              // real (ordinary) DM
              vidx = indgen(indexDm(1,nm):indexDm(2,nm));
              fMat(idx,vidx) = float(unit(dm(nm)._nact));
            } else {
              // tomographic DM
              vidx = [];
              virtualDMs = int(*dm(nm).dmfit_which);
              for (c1=1;c1<= numberof(virtualDMs);c1++){
                grow, vidx, indgen(indexDm(1,virtualDMs(c1)):indexDm(2,virtualDMs(c1)));
              }
              fMat(idx,vidx) =  *dm(nm)._fMat;
            }
          }
        }

        Dterm = Ga(,+)*fMat(+,)-Gx;
        polcMat = Gx(+,)*Dterm(+,)-Cphi;
        dMat = LUsolve(AtA+Cphi,polcMat);
        yao_fitswrite, YAO_SAVEPATH + parprefix + "-dMat.fits", dMat;
      } else {
        nAct = (dimsof(iMat))(3);
        Cphi = array(float,[2,nAct,nAct]);

        mc = 0; // matrix counter
        nDMs = numberof(dm);
        for (nm=1;nm<=nDMs;nm++){
          if (!dm(nm).dmfit_which){
            Cphi((mc+1):(mc+dm(nm)._nact),(mc+1):(mc+dm(nm)._nact)) = \
                                  (*dm(nm)._regmatrix)*dm(nm).regparam;
            mc += dm(nm)._nact;
          }
        }
        cMat = LUsolve(iMat(+,)*iMat(+,)+Cphi,transpose(iMat));
      }

    } else if (mat.method == "mmse-sparse") {

      // create the regularization matrices for each DM
      nAct = iMatSP.r;
      nDMs = numberof(dm);

      CphiSP = [];

      for (nm=1;nm<=nDMs;nm++) {
        regmatrix = *dm(nm)._regmatrix;
        ruox, regmatrix, dm(nm).regparam; // multiply the matrix by a constant

        if (!dm(nm).dmfit_which){
        if (CphiSP == []) {
          CphiSP = regmatrix;
          } else {
          spcon, CphiSP,regmatrix, diag=1, ruo=1;
          }
        }
      }

      if (anyof(dm.dmfit_which)) {
        estAct = [];  // actuators used to estimate wavefront
        realAct = []; // real (non virtual) actuators used to compensate
                      // the wavefront
        actno = 0;

        for (nm=1;nm<=ndm;nm++){
          if (!dm(nm).dmfit_which){
            grow, estAct, indgen(1+actno:actno+dm(nm)._nact);
          }
          if (!dm(nm).virtual & !dm(nm).ncp){
            grow, realAct, indgen(1+actno:actno+dm(nm)._nact);
          }
          actno += dm(nm)._nact;
        }

        // create copies of iMatSP and remove the rows corresponding to
        // virtual actuators or tomographic actuators. Need to copy the
        // data in the pointers to avoid modifying the data in iMatSP.

        GaSP = iMatSP;
        GaSP.ix = &(*iMatSP.ix);
        GaSP.jx = &(*iMatSP.jx);
        GaSP.xn = &(*iMatSP.xn);

        GxSP = iMatSP;
        GxSP.ix = &(*iMatSP.ix);
        GxSP.jx = &(*iMatSP.jx);
        GxSP.xn = &(*iMatSP.xn);

        for (a1=iMatSP.r;a1>=1;a1--){
          if (noneof(estAct == a1)){
            rcodr, GxSP, a1;
          }
          if (noneof(realAct == a1)){
            rcodr, GaSP, a1;
          }
        }

        save_rco,GxSP,YAO_SAVEPATH+parprefix+"-GxSP.rco";

        // create a global fitting matrix

        fMatSP = [];

        indexDm       = array(long,2,ndm);
        indexDm(,1)   = [1,dm(1)._nact];
        for (nm=2;nm<=ndm;nm++) {
          indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
        }

        nEstAct =sum(dm(where(!dm.dmfit_which))._nact);
        nEstDMs = numberof(where(!dm.dmfit_which));
        nRealAct = sum(dm(where(!dm.virtual & !dm.ncp))._nact);

        for (nm=1;nm<=numberof(dm);nm++) {
          if (dm(nm).virtual){ continue; } // no virtual DMs in the rows
          if (dm(nm).ncp){ continue; } // no non-common path DMs in the rows

          if (dm(nm).dmfit_which) {
            dmfMat = rcotr(*dm(nm)._fMat);
            vdmidx = [];
            for (mm=1;mm<=nEstDMs;mm++){
              if (anyof(mm == *dm(nm).dmfit_which)){
                grow, vdmidx, indgen(indexDm(1,mm):indexDm(2,mm));
              }
            }

            for (c1=1;c1<=dm(nm)._nact;c1++){
              row = array(float,[2,nEstAct,1]);
              vec = array(float,nRealAct);
              vec(c1) = 1.;
              row(vdmidx) = rcoxv(dmfMat, vec); // extract row c1

              if (fMatSP == []){
                fMatSP = sprco(row);
              } else {
                rcobuild, fMatSP, row, mat.sparse_thresh;
              }
            }

          } else { // an ordinary DM
            for (c1=1;c1<=dm(nm)._nact;c1++){
              row = array(float,[2,nEstAct,1]);
              row(indexDm(1,nm)+c1-1)=1;
              if (fMatSP == []){
                fMatSP = sprco(row);
              } else {
                rcobuild, fMatSP, row, mat.sparse_thresh;
              }
            }
          }
        }

        t1 = rcoatb(rcotr(fMatSP),rcotr(GaSP));
        fMatSP = GaSP = [];

        AtA = rcoata(GxSP);
        AtAregSP = ruoadd(AtA,CphiSP);

        AtA = [];

        (*GxSP.xn) *= -1;
        DtermSP = rcoadd(t1,GxSP);
        t1 = [];
        (*GxSP.xn) *= -1;
        t2 = rcoatb(DtermSP,GxSP);

        DtermSP = [];
        CphiSPrco = ruo2rco(CphiSP);
        CphiSP = [];
        *CphiSPrco.xn *= -1;
        polcMatSP = rcotr(rcoadd(t2,CphiSPrco));
        t2 = CphiSPrco = [];
        save_rco,polcMatSP,YAO_SAVEPATH+parprefix+"-polcMat.rco";
      } else {
        GxSP = iMatSP;
        AtA = rcoata(iMatSP);
        AtAregSP = ruoadd(AtA,CphiSP);
        AtA = CphiSP =  [];
        save_rco,GxSP,YAO_SAVEPATH+parprefix+"-GxSP.rco";
      }

      save_rco,iMatSP,YAO_SAVEPATH+mat.file;
      iMatSP = [];
      save_ruo,AtAregSP,YAO_SAVEPATH+parprefix+"-AtAreg.ruo";
    }

    if (mat.method != "mmse-sparse") {

      // More debug display
      if (sim.debug>=3){
        tv,cMat(,+)*iMat(+,);
        mypltitle,"cMat(,+)*iMat(+,)",[0.,-0.005],height=12;
        typeReturn;
        tv,iMat(,+)*cMat(+,);
        mypltitle,"iMat(,+)*cMat(+,)",[0.,-0.005],height=12;
        typeReturn;
      }

      // save the results:
      if (noneof(dm.dmfit_which)){
      yao_fitswrite,YAO_SAVEPATH+mat.file,[iMat,transpose(cMat)];
      } else { // just save the iMat and force the recreation of cMat on reload
        yao_fitswrite,YAO_SAVEPATH+mat.file,[iMat,iMat];
      }
    }
  }

  if (mat.use_cmat_file) {
    if (sim.verbose) write,format=">> Reading requested cmat %s\n",mat.use_cmat_file;
    cMat = yao_fitsread(mat.use_cmat_file);
  }


  //===================================================================
  // COMPUTE THE COMMAND VECTOR FOR OFFLOADING THE ANISOPLANATISM MODES
  //===================================================================

  if (anyof(dm.type == "aniso")) {

    // find which DM is the aniso DM:
    nmaniso = where(dm.type == "aniso");
    if (numberof(nmaniso) != 1) {
      pyk_error,"there can be only one aniso DM !";
      error,"there can be only one aniso DM !";
    }

    // finds which DM is at altitude 0
    w0 = where((dm.alt == 0)*(dm.virtual == 0));
    nmlow = where( (dm(w0).type == "stackarray") | (dm(w0).type == "bimorph") |
                   (dm(w0).type == "zernike") | (dm(w0).type == "dh")  );
    nmlow = w0(nmlow);
    if (numberof(nmlow) == 0) {
      pyk_error,"I cannot find a DM at altitude 0 to produce the lower "+
        "part of the anisoplanatism modes !";
      error,"I cannot find a DM at altitude 0 to produce the lower "+
        "part of the anisoplanatism modes !";
    }
    if (numberof(nmlow) > 1) {
      pyk_warning,swrite(format=\
        "Weird. There are %d high-order DMs at altitude=0 (look at console)",\
        numberof(nmlow));
      write,format="Weird. There are %d high-order DMs at altitude=0:", \
            numberof(nmlow);
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
                    (dm(wn0).type == "zernike") | (dm(wn0).type == "dh")  );
    nmhigh = wn0(nmhigh);
    if (numberof(nmhigh) == 0) {
      pyk_error,"I cannot find a DM at the requested altitude to produce "+ \
                "the higher part of the anisoplanatism modes !";
      error,"I cannot find a DM at the requested altitude to produce "+ \
            "the higher part of the anisoplanatism modes !";
    }
    if (numberof(nmhigh) > 1) {
      pyk_warning,swrite(format=\
        "Weird. There are %d high-order DMs at altitude %.0f (look at console)",
                         numberof(nmhigh),dm(nmaniso).alt);
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
    project_aniso_dm,nmaniso(1),nmlow(1),nmhigh(1),disp=0;
  }

  //==========================
  // OUTPUT GRAPHIC FOR CONFIG
  //==========================

  if ((disp == 1)&&(!yaopy)) graphic_config;
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
        ( (wfs(ns).noise) ? "enabled" : "disabled" );
    }
    write,format="D/r0 (500nm) = %.1f; %d iterations\n",atm.dr0at05mic/
      cos(gs.zenithangle*dtor)^0.6,loop.niter;
  }

  // same in result file:
  f = open(YAO_SAVEPATH+parprefix+".res","a+");
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
      ( (wfs(ns).noise) ? "enabled" : "disabled" );
  }
  write,f,format="D/r0 (500nm) = %.1f; %d iterations\n",atm.dr0at05mic/
    cos(gs.zenithangle*dtor)^0.6,loop.niter;
  close,f;

  // make sure kernelconv is good after the imat:
  for (ns=1;ns<=nwfs;ns++) {
    wfs(ns)._kernelconv = 1n;
    if ((wfs(ns).kernel==0)&&(wfs(ns)._gsdepth==0)) wfs(ns)._kernelconv = 0n;
    if (wfs(ns).shmethod==1) wfs(ns)._kernelconv = 0n;
  }

  // basic initialization in case of svipc use:
  if (sim.svipc) require,"yao_svipc.i";

  if (aoinit_user_func!=[]) status = aoinit_user_func();

  gui_message,"aoinit done: click on aoloop";
  notify,swrite(format="%s: aoinit done",parprefix);
}

//---------------------------------------------------------------
func aoloop(disp=,savecb=,dpi=,controlscreen=,nographinit=,anim=,savephase=,no_reinit_wfs=)
/* DOCUMENT aoloop(disp=,savecb=,dpi=,controlscreen=,nographinit=,anim=,savephase=)
   Prepare all arrays for the execution of the loop (function go()).

   disp          = set to display stuff as the loop goes
   savecb        = set to save "circular buffers"
   dpi           = dpi of graphic window
   controlscreen = set to display an additional graphic window with main
                   loop parameters updating as the loop goes.
   nographinit   = not sure. this must be for the GUI
   anim          = set to 1 to use double buffering. Avoids flicker but
                   might confuse the user if not turned off in normal operation
                   see animate() yorick function. ON by default.
   savephase     = set to save residual wavefront on first target as fits file
   no_reinit_wfs = don't re-init wfs
   SEE ALSO: aoread, aoinit, go, restart
 */
{
  extern looptime, mircube, command, wfsMesHistory, cubphase;
  extern im, imav, airy, sairy, fwhm, e50;
  extern niterok;
  extern pp, sphase, bphase, imphase, dimpow2, cMat, pupil, ipupil;
  extern time, strehllp, strehlsp, rpv, itv, ok, njumpsinceswap;
  extern remainingTimestring  ;
  extern cbmes, cbcom, cberr;
  extern indexDm, aniso, waniso, wdmaniso;
  extern ditherPeriod, ditherAmp, ditherGain, cggain, ditherMes;
  extern ditherCosLast, ditherSinLast, ditherMesCos, ditherMesSin;
  extern dispImImav;
  extern loopCounter, nshots, statsokvec;
  extern savecbFlag, dpiFlag, controlscreenFlag;
  extern dispFlag, nographinitFlag, animFlag;
  extern aoloop_disp,aoloop_savecb,aoloop_no_reinit_wfs;
  extern commb,errmb; // minibuffers
  extern default_dpi;
  extern savephaseFlag;

  if ((disp==[])&&(aoloop_disp!=[])) disp=aoloop_disp;
  if ((savecb==[])&&(aoloop_savecb!=[])) savecb=aoloop_savecb;
  if ((no_reinit_wfs==[])&&(aoloop_no_reinit_wfs!=[])) no_reinit_wfs=aoloop_no_reinit_wfs;
  if (anim==[]) anim=1; // let's make it the default

  dispFlag = disp;
  savecbFlag = savecb;
  dpiFlag = dpi;
  controlscreenFlag = controlscreen;
  nographinitFlag = nographinit;
  animFlag = anim;
  savephaseFlag = savephase;

  gui_message,"Initializing loop";

  check_control_parameters,verbose=1; // update the control laws

  for (ns=1;ns<=nwfs;ns++) {
    if (wfs(ns).framedelay == -1){
      wfs(ns)._framedelay = loop.framedelay;
    } else {
      wfs(ns)._framedelay = wfs(ns).framedelay;
      if ((sim.verbose) && (loop.framedelay != wfs(ns).framedelay)){
        write, format="Using a loop delay of %d frames for wfs(%d)\n",wfs(ns).framedelay, ns;
      }
    }
  }

  // Initialize hysteresis parameters
  for (nm=1;nm<=ndm;nm++){
    dm(nm)._x0 = &array(float,dm(nm)._nact);
    dm(nm)._y0 = &array(float,[2,dm(nm)._nact,3]);
    dm(nm)._xlast = &array(float,dm(nm)._nact);
    dm(nm)._ylast =  &array(float,[2,dm(nm)._nact,3]);
    dm(nm)._signus = &array(long,dm(nm)._nact);
  }

  if (aoloop_user_func!=[]) status = aoloop_user_func();

  // Initialize displays:
  if (!is_set(disp)) {disp = 0;}
  if (!is_set(controlscreen)) {controlscreen = 0;}
  if (is_set(disp) && !is_set(nographinit)) {
    if (!yaopy) status = create_yao_window();
  }
  if (is_set(controlscreen) && !is_set(nographinit)) {
    control_screen,0,init=1;
  }

  size    = sim._size;

  // Some arrays initialization:
  looptime      = 0.;
  mircube       = array(float,[3,size,size,ndm]);
  command       = array(float,sum(dm._nact));
  wfsMesHistory = array(float,[2,sum(wfs._nmes),max(wfs._framedelay)+1]);
  cubphase      = array(float,[3,size,size,target._ntarget]);
  im            = array(float,[3,size,size,target._ntarget]);
  imav          = array(float,[4,size,size,target._ntarget,target._nlambda]);
  if (*target._pupil == []){
    airy          = calc_psf_fast(pupil,pupil*0.);
    sairy         = max(airy);
  } else {
    sairy = array(float,target._ntarget);
    for (nt=1;nt<=target._ntarget;nt++){
      airy = calc_psf_fast((*target._pupil)(,,nt),pupil*0.);
      sairy(nt) = max(airy);
    }
  }

  fwhm   = e50  = array(float,[2,target._ntarget,target._nlambda]);
  pp            = array(complex,[2,size,size]);
  sphase        = array(float,[2,_n,_n]);
  bphase        = array(float,[2,size,size]);
  imphase       = array(float,[2,size,size]);
  dimpow2       = int(log(size)/log(2));
  cMat          = float(cMat);
  pupil         = float(pupil);
  ipupil        = float(ipupil);
  time          = array(float,10);
  strehllp = strehlsp = itv = rpv = rp_tip1d = rp_tilt1d = [];
  ok            = 0;
  niterok       = 0;
  remainingTimestring = "";
  njumpsinceswap = 0; // number of jumps since last screen swap

  // minibuffers for up to 10th order filters (control laws):
  // determine number of actuators and commands
  nComm = 0;
  nAct = 0;

  for (nm=1;nm<=ndm;nm++){
    if (dm(nm).virtual == 0) nComm += dm(nm)._nact;
    if (*dm(nm).dmfit_which == []) nAct += dm(nm)._nact;
  }

  commb          = array(float,[2,nComm,10]);
  errmb          = array(float,[2,nAct,10]);

  if (is_set(savecb)) {
    cbmes         = array(float,[2,sum(wfs._nmes),loop.niter]);
    cbcom         = array(float,[2,sum(nComm),loop.niter]);
    cberr         = array(float,[2,sum(nAct),loop.niter]);
  }

  // Special: re-init the WFS, since it includes the computation of
  // the # of photons/WFS and the WFS bias anf flats, in case something
  // has changed since the aoinit. The # of photons in particular is one
  // parameter we often want to vary without having to go through the
  // whole aoinit again. For a curvature wfs, it re-inits the # of
  // photons and the extra-focal distance, so it is possible to loop
  // on "aoloop" with various l without the need for a aoinit (although
  // you would still be using the cMat determined for the original l,
  // but that might be on purpose).
  // Sometimes, a re-init might take a long time (e.g. for ELT systems),
  // and might be unnecessary, so the user might want to skip it. Use
  // no_reinit_wfs=1 to do so.
  for (ns=1;ns<=nwfs;ns++) {
    // zero fsm! got burn by this one!
    wfs(ns)._tt *= 0;
    wfs(ns)._upttcommand *= 0;
    // define _dispimage
    if (wfs(ns).type=="zernike") continue;
    if (wfs(ns).type=="dh") continue;
    wfs(ns)._dispimage = &(*wfs(ns)._fimage*0.0f);
    if (!no_reinit_wfs) {
      if (wfs(ns).type == "hartmann") {
        if (wfs(ns).disjointpup) {
          shwfs_init,disjointpup(,,ns),ns,silent=1;
        } else shwfs_init,*wfs(ns)._pupil,ns,silent=1;
      } else if (wfs(ns).type == "curvature") {
        curv_wfs,*wfs(ns)._pupil,*wfs(ns)._pupil*0.0f,ns,init=1,silent=1;
      }
    } else { // even if no_reinit asked, we need to do at least
      // the following as it may have been left in an imat state.
      if (wfs(ns).type == "hartmann") \
        shwfs_init_common_kernel,ns,imat=0;
    }
  }

  statsokvec = (indgen(loop.niter)-1);
  if (loop.skipevery > 0){
    statsokvec = statsokvec % loop.skipevery;
  }
  statsokvec = (statsokvec >= loop.startskip);

  // special: re-init phase screens, as often we want to change the
  // number of iterations without re-doing the full aoinit
  // note: ok, but have to quit_forks() explicitely if using svipc.
  if (loop.niter > dimsof(xposvec)(2)) get_turb_phase_init;

  // reset cyclecounters
  wfs._cyclecounter = wfs._cyclecounter*0 +1;

  indexDm       = array(long,2,ndm);
  indexDm(,1)   = [1,dm(1)._nact];
  for (nm=2;nm<=ndm;nm++) {
    indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
  }

  for (nm=1;nm<=ndm;nm++) {
    dm(nm)._command = &(array(float,dm(nm)._nact));
  }

  // Set up for anisoplatism modes
  aniso = anyof(dm.type == "aniso");
  if (aniso) {
    waniso = indexDm(1,where(dm.type == "aniso"))+[0,1,2];
    wdmaniso = where(dm.type == "aniso")(0);
    if (sim.verbose>1) {
      print,"index of anisoplatism modes in command vector:",waniso;
    }
  }

  // dithering set up for upTT compensation (LGS only)
  wfs._centroidgain = array(1.f,nwfs);
  //  if (anyof(wfs.centGainOpt)) {
  // we dither if there is an uplink, filtertilt is set and
  // centGainOpt is requested
  ditherPeriod = 8;
  ditherAmp = 0.05;  // supposed to be in arcsec
  ditherGain = 0.1;
  cggain = 0.1; // 0.01;
  ditherCosLast = ditherSinLast = 0.0;
  ditherMesCos = ditherMesSin = ditherMes = array(0.,nwfs);
  //}

  // verbose
  if (sim.verbose) {
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

    txpos = *target.xposition;
    typos = *target.yposition;
    disp2d,im,txpos,typos,1,zoom=*target.dispzoom,init=1;
      if (wfs_display_mode=="spatial") {
        disp2d,wfs._dispimage,wfs.pupoffset(1,),wfs.pupoffset(2,),2,\
           zoom=wfs.dispzoom,init=1;
      } else {
        disp2d,wfs._dispimage,wfs.gspos(1,),wfs.gspos(2,),2,zoom=wfs.dispzoom,\
           init=1;
      }
  }

  dispImImav = 0; // 0 is to display im, 1 to display imav

  olddisp = disp;
  oldcontrolscreen = controlscreen;

  // Main loop:

  loopCounter=0;
  nshots = -1;

  if (sim.svipc) {
    status = svipc_init();

    status = svipc_start_forks();

    status = init_sync();
    status = reset_strehl();
  }

  if (animFlag && dispFlag) {
    plsys,1; animate,1;
  }

  gui_message,"Initializing loop...done. \"go\" to start, \"pause\" to pause";
}


func go(nshot,all=)
/* DOCUMENT
   go will start or resume the AO loop
   go will go once through this function, and then call itself, until
     loop.niter has been reached or nshot has been done since last user-
     entered go. Because go() comes back to the main yorick prompt after each
     iteration, the user can enter commands, or ask to pause the loop, at any
     time.
   example:
   ... aoloop,disp=1
   go
   ... the loop goes, one iteration after another and prints out thing on screen
   dispFlag=0      // turns off display. the loop keeps going
   loop.gain = 0.2 // change the loop integrator gain
   reset           // flatten the DMs
   stop            // pause the loop
   cont            // continues the loop

   nshot = # of iterations to go through. if not set, will go until
           loop.niter (total # of iterations) has been done
   all   = do all remaining iteration without returning to the yorick prompt.
     This is useful in scripts to avoid messy side effects of set_idler.
   SEE ALSO: aoloop, stop, cont, reset, restart, whereat, after_loop
 */

{
  extern niterok, starttime, endtime, endtime_str, tottime, now2, nshots;
  extern dispFlag, savecbFlag, dpiFlag, controlscreenFlag;
  extern nographinitFlag,savephaseFlag;
  extern cbmes, cbcom, cberr;
  extern iMatSP, AtAregSP
  extern commb,errmb; // minibuffers (last 10 iterations)
  extern im,imav;
  extern iter_per_sec;

  gui_show_statusbar1;

  if (go_user_func!=[]) status = go_user_func();

  if (nshot!=[]) {
    if (animFlag&&dispFlag) {
      plsys,1;
      animate,1;
    }
    nshots = nshot;
  }
  disp = dispFlag;
  savecb = savecbFlag;
  dpi = dpiFlag;
  controlscreen = controlscreenFlag;
  nographinit = nographinitFlag;
  savephase = savephaseFlag;

  if (loopCounter==0) {
    // initialize timers
    tic,2;
    starttime = _nowtime(2);
    tottime = 0.;
  }

  go_start:  now = tac(2);

  if (go2_user_func!=[]) status = go2_user_func();

  loopCounter++;
  nshots--;

  // update gui progress bar, but don't do it too often.
  if (max_progressbar_update_freq==[]) max_progressbar_update_freq=5.; // in Hz
  if (tac(3)>(1./max_progressbar_update_freq)) {
    gui_progressbar_frac,float(loopCounter)/loop.niter;
    gui_progressbar_text,swrite(format="%d out of %d iterations",loopCounter, \
      loop.niter);
    tic,3;
  }

  comvec=[];

  if (loopCounter>loop.niter) {
    exit,"Can't continue: loopCounter > loop.niter !";
  }

  i = loopCounter;
  tic; time(1) = tac();


  prevOK  = ok;

  check_control_parameters; // update the control laws

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
    // subtract last command and add new. the circular motion is
    // known. radius=ditherAmp.
    wfs(wdither)._upttcommand -= [ditherCosLast,ditherSinLast]*nditherarray(-,);
    wfs(wdither)._upttcommand += [ditherCosCurr,ditherSinCurr]*nditherarray(-,);
    ditherCosLast = ditherCosCurr;
    ditherSinLast = ditherSinCurr;
    // do the temporal filtering on cos and sin measurements
    // (see requirements documents I wrote for tOSC regarding the process)
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
  if (anyof([sim.svipc>>0,sim.svipc>>2]&1)) { // parallel computing.
    svipc_wfsmes = topwfs_svipc();
  } else WfsMes = mult_wfs(i,disp=disp);
  time(2) +=  tac();

  ok  = statsokvec(i);

  if ((prevOK-ok)==1) { // After a jump
    njumpsinceswap++;
    if (njumpsinceswap==loop.jumps2swapscreen) {
      if (sim.verbose) write,"Reset and screens swap";
      swap_screens;        // swap screens
      njumpsinceswap = 0; // reset "jumps since last swap" counter
    } else { write,"Reset"; }
    status = reset();
  }

  // Handling frame delay (0 -> no frame delay at all, 1 -> regular
  // one frame delay of integrator, 2 -> integrator + one frame
  // computation, etc...
  wfsMesHistory  = roll(wfsMesHistory,[0,-1]);

  mc = 0; // measurement counter
  for (ns=1;ns<=nwfs;ns++){
    if (anyof([sim.svipc>>0,sim.svipc>>2]&1)) { // parallel computing.
      //wfsMesHistory(,loop.framedelay) = svipc_wfsmes;
      wfsMesHistory(mc+1:mc+wfs(ns)._nmes,wfs(ns)._framedelay) = svipc_wfsmes(mc+1:mc+wfs(ns)._nmes);
    } else {
      //wfsMesHistory(,loop.framedelay+1) = WfsMes;
      wfsMesHistory(mc+1:mc+wfs(ns)._nmes,wfs(ns)._framedelay+1) = WfsMes(mc+1:mc+wfs(ns)._nmes);
    }
    mc += wfs(ns)._nmes;
  }

  usedMes        = float(wfsMesHistory(,1));

  time(3) += tac();

  // RECONSTRUCTION:
  // computes the actuator error vector from the measurements:
  if (mat.method == "mmse-sparse") {
    if (i == 1){
      MR = mat.sparse_MR;
      MN = mat.sparse_MN;
    }

    if (sum(usedMes != 0) == 0){ // sparse CG method does not work if all zeros
      err = array(float,AtAregSP.r);
    } else {
      Ats=rcoxv(GxSP,usedMes)
      if ((loop.method == "pseudo open-loop") && (i > 1) && (polcMatSP != [])){
        // apply polc correction
        polccorr = rcoxv(polcMatSP,estdmcommand);
        if (anyof(wfs.nintegcycles > 1)){
          // if some DMs are not updating because the WFSs have not finished
          // integrating, then the POLC correction should not be applied to
          // those DMs
          // Check to see which DMs have updated; only non tomographic DMs
          mc = 0; // counter
          for (nm=1;nm<=numberof(dm);nm++){
            if (!dm(nm).dmfit_which){
              comms = Ats(mc+1:mc+dm(nm)._nact);
              if (comms(rms) == 0){ // has not updated
                polccorr(mc+1:mc+dm(nm)._nact) = 0;
              }
              mc += dm(nm)._nact;
            }
          }
        }
        Ats -= polccorr;
      }
      err = float(ruopcg(AtAregSP,Ats, array(float,AtAregSP.r), tol=mat.sparse_pcgtol));
    }
  } else {
    // if loop.method = "none", everything is left to the user,
    // thus we skip this matrix multiply.
    if (loop.method != "none") err = cMat(,+) * usedMes(+);
    if ((loop.method == "pseudo open-loop") && (i > 1)){
      // make sure that dMat exists
      if (dMat == []){
        // if not, load it
        if (fileExist(YAO_SAVEPATH+parprefix+"-dMat.fits")){
          dMat = yao_fitsread(YAO_SAVEPATH + parprefix + "-dMat.fits");
        } else {
          error, "loop.method set to pseudo open-loop, but dMat is not defined";
        }
      }
      polccorr = dMat(,+)*estdmcommand(+); //pseudo open loop correction
      if (anyof(wfs.nintegcycles > 1)){
        // if some DMs are not updating because the WFSs have not finished
        // integrating, then the POLC correction should not be applied to
        // those DMs.
        // Check to see which DMs have updated; only non tomographic DMs
        mc = 0; // counter
        for (nm=1;nm<=numberof(dm);nm++){
          if (!dm(nm).dmfit_which){
            comms = err(mc+1:mc+dm(nm)._nact);
            if (comms(rms) == 0){ // has not updated
              polccorr(mc+1:mc+dm(nm)._nact) = 0;
            }
            mc += dm(nm)._nact;
          }
        }
      }
      err -= polccorr;
    }
  }


  if (user_loop_err!=[]) user_loop_err;

  // get the anisoplanatism mode coefficients and project it
  // in the actuator space
  if (aniso) {
    err += dm(wdmaniso).gain*(comaniso(,+) * err(waniso)(+));
    err(waniso) = 0.;
    //      dm(wdmaniso)._command = invcomaniso(+,)*
  }

  time(4) += tac();

  // Computes the mirror shape using influence functions:
  for (nm=1; nm<=ndm; nm++) {

    if (loop.method!="none") {
      if (dm(nm).type == "aniso") {
        grow,comvec,*dm(nm)._command;
        continue;
      }

      if (*dm(nm).dmfit_which == []) { // not a tomographic DM
        // update the command vector for this DM:
        // this DM error:
        ctrlden = *dm(nm)._ctrlden;
        ctrlnum = *dm(nm)._ctrlnum;

        command = *dm(nm)._command*(-ctrlden(2));
        for (order=3;order<=numberof(ctrlden);order++) {
          imb = (i-order+1)%10;
          command -= commb(indexDm(1,nm):indexDm(2,nm),imb)*ctrlden(order);
        }

        dmerr = err(indexDm(1,nm):indexDm(2,nm));
        command -= dmerr*ctrlnum(1);
        for (order=2;order<=numberof(ctrlnum);order++) {
          imb = (i-order+1)%10;
          command -=errmb(indexDm(1,nm):indexDm(2,nm),imb)*ctrlnum(order);
        }

        *dm(nm)._command = command;

      } else { // tomographic DM; DM commands from virtual DMs
        virtualDMs = int(*dm(nm).dmfit_which);
        virtualdmcommand = [];
        for (idx=1;idx<= numberof(virtualDMs);idx++) {
          grow, virtualdmcommand, *dm(virtualDMs(idx))._command;
        }
        if (mat.method == "mmse-sparse") {
          *dm(nm)._command = rcoxv(*dm(nm)._fMat,virtualdmcommand);
        } else {
          *dm(nm)._command = (*dm(nm)._fMat)(,+)*virtualdmcommand(+);
        }
      }

      if (dm(nm).filterpiston) { // filter piston
        if (dm(nm).type == "stackarray") {
          *dm(nm)._command -= avg(*dm(nm)._command);
        }
      }

      if (dm(nm).filtertilt) { // filter tip and tilt
        if (dm(nm).type == "stackarray") {
          // todo: have a variable to store xv and yv to avoid recomputing
          xv = *dm(nm)._x - avg(*dm(nm)._x);
          xv = xv / sqrt(sum(xv*xv));

          yv = *dm(nm)._y - avg(*dm(nm)._y);
          yv = yv / sqrt(sum(yv*yv));

          *dm(nm)._command -= avg(*dm(nm)._command);
          *dm(nm)._command -= (xv(+)*(*dm(nm)._command)(+))*xv;
          *dm(nm)._command -= (yv(+)*(*dm(nm)._command)(+))*yv;
        }

        if (dm(nm).type == "zernike") {
          if (dm(nm).minzer <= 3) {
            (*dm(nm)._command)(1:4-dm(nm).minzer) = 0;
          }
        }
      }

      estdmcommand = [];
      for (idx=1;idx<=ndm;idx++) {
        if (!dm(idx).dmfit_which) {
          grow, estdmcommand, *dm(idx)._command;
        }
      }

      if (dm(nm).maxvolt != 0) {
        dm(nm)._command=&(float(clip(*dm(nm)._command,-dm(nm).maxvolt,dm(nm).maxvolt)));
      }
    }

    n1 = dm(nm)._n1; n2 = dm(nm)._n2; nxy = n2-n1+1;

    if (user_loop_command!=[]) user_loop_command,nm;

    if (dm(nm).virtual == 0) {
      grow,comvec,*dm(nm)._command;
      if (dm(nm).hyst > 0) {
        mircube(n1:n2,n1:n2,nm) = comp_dm_shape(nm,&(hysteresis(*dm(nm)._command,nm)));
      } else {
        mircube(n1:n2,n1:n2,nm) = comp_dm_shape(nm,dm(nm)._command);
      }
      if ( (nm==1) && (add_dm0_shape!=[]) ) mircube(,,1) += add_dm0_shape;
      // extrapolated actuators:
      if ( (dm(nm)._enact != 0) && (dm(nm).noextrap == 0) ) {
        ecom = float(((*dm(nm)._extrapcmat)(,+))*(*dm(nm)._command)(+));
        mircube(n1:n2,n1:n2,nm) += comp_dm_shape(nm,&ecom,extrap=1);
      }
    } else {
      mircube(n1:n2,n1:n2,nm) = 0.; //virtual DM, does not produce a shape
    }
  }

  // fill minibuffers
  if (err!=[]) errmb(,(i%10)) = err; // 10 will go in 0, which is 10.
  commb(,(i%10)) = comvec; // 10 will go in 0, which is 10.

  if (tipvib!=[]) { // add tip vibrations
    // if there is a TTM, add it there
    if (anyof(dm.type=="tiptilt")) nmvib=where(dm.type=="tiptilt")(1);
    else nmvib=1; // otherwise put it to 1.
    mircube(,,nmvib) += tipvib(i)*tip1arcsec;
  }

  if (tiltvib!=[]) { // add tilt vibrations
    // if there is a TTM, add it there
    if (anyof(dm.type=="tiptilt")) nmvib=where(dm.type=="tiptilt")(1);
    else nmvib=1; // otherwise put it to 1.
    mircube(,,nmvib) += tiltvib(i)*tilt1arcsec;
  }

  if (segmenttiptiltvib!=[]) { // add tiptilt vibrations on a segment by segment basis, designed primarily for GMT but could be used for any telescope
    // if there is a TTM, add it there
    if (anyof(dm.type=="tiptilt")) nmvib=where(dm.type=="tiptilt")(1);
    else nmvib=1; // otherwise put it to 1.
    mircube(,,nmvib) += make_segment_tiptilt(segmenttiptiltvib(,i));
  }

  if (segmentpistonvib!=[]) { // add piston vibrations on a segment by segment basis, designed primarily for GMT but could be used for any telescope
    nmvib=1;
    mircube(,,nmvib) += make_segment_piston(segmentpistonvib(,i));
  }

  time(5) += tac();

  ok = ok*((i % loop.stats_every) == 0); // accumulate stats only every 4 iter.
  if (user_go_ok) ok = user_go_ok();

  okdisp = ( is_set(disp) && (((i-1) % disp) == 0) );
  if (nshots>0) okdisp=is_set(disp);
  okcscreen = ( is_set(controlscreen) && (((i-1) % controlscreen) == 0) );
  if (is_set(controlscreen) && (i == loop.niter)) okcscreen=1;  // display at last iteration

  if (savephase||ok||okdisp) {
    // get the residual phase; initially for the first target
    // display and save the residual wavefront
    if (residual_phase_what==[]) residual_phase_what="target";
    if (residual_phase_which==[]) residual_phase_which=1;

    residual_phase=get_phase2d_from_dms(residual_phase_which,residual_phase_what) +
                   get_phase2d_from_optics(residual_phase_which,residual_phase_what) +
                   get_turb_phase(i,residual_phase_which,residual_phase_what);

    residual_phase1d = residual_phase(where(pupil > 0));
    residual_phase1d -= min(residual_phase1d);
    residual_phase(where(pupil > 0)) = residual_phase1d;

    if (savephase){
      fname = swrite(format="%s/%s_rwf%d.fits",YAO_SAVEPATH,parprefix,loopCounter);
      yao_fitswrite,fname,float(residual_phase*pupil);
    }

    if (residual_phase_zcontent) {
      // project residual phase on zernike
      // init zernikes:
      extern rp_zerns,rp_nzerns,rp_w,rp_rms,rp_z;
      if (rp_zerns==[]) {
        if (!rp_nzerns) rp_nzerns = 45;
        rp_zerns = array(float,[3,_n,_n,rp_nzerns]);
        // I've done the zernike thingy in a very stupid way.
        // it might have been used/initialized elsewhere and I can
        // mess up this initialisation here. To avoid this, temporarily:
        if (zrmod!=[]) error,"zernike have already been initialized, can't do!";
        prepzernike,_n,sim.pupildiam;
        for (nz=1;nz<=rp_nzerns;nz++) rp_zerns(,,nz) = zernike(nz);
        rp_zerns = rp_zerns(*,)(where(pupil(_n1:_n2,_n1:_n2)),);
        rp_zerns = QRsolve(transpose(rp_zerns),unit(rp_nzerns));
        rp_w = where(pupil);
      }
      rp1d = residual_phase(rp_w);
      rp_rms = rp1d(rms)*1000.; // in nm
      rp_z = rp1d(+)*rp_zerns(+,)*1000.; // in nm
    } else if (residual_phase_rms_nott) {
      extern rp_tip1d,rp_tilt1d;
      if (rp_tip1d==[]) {
        tmp = indices(dimsof(pupil)(0));
        rp_tip1d  = tmp(,,1)(where(pupil > 0));
        rp_tilt1d = tmp(,,2)(where(pupil > 0));
        rp_tip1d -= avg(rp_tip1d);
        rp_tilt1d -= avg(rp_tilt1d);
        rp_tip1d = rp_tip1d/sqrt(sum(rp_tip1d^2.));
        rp_tilt1d = rp_tilt1d/sqrt(sum(rp_tilt1d^2.));
      }
      rp1d_nott = residual_phase1d;
      rp1d_nott -= sum(rp1d_nott*rp_tip1d)*rp_tip1d;
      rp1d_nott -= sum(rp1d_nott*rp_tilt1d)*rp_tilt1d;
      if (ok) grow,rpv,rp1d_nott(rms)*1000.; // in nm
    }
  }

  // Computes the instantaneous PSF:
  if (ok) {
    if ((sim.svipc>>1)&1) { // use child
      if (psf_child_started) {
        // read previous results:
        if (smdebug) write,"main: waiting for ready from PSF fork";
        sem_take,semkey,4;
        if (smdebug) write,"main: received ready from PSF fork";
        // results are ready.
        im   = shm_read(shmkey,"imsp");
        imav = shm_read(shmkey,"imlp");
        niterok += 1;
        grow,itv,i;
        if (disp_strehl_indice) sind=disp_strehl_indice; else sind=1;
        grow,strehlsp,im(max,max,sind)/sairy(sind);
        grow,strehllp,imav(max,max,sind,0)/sairy(sind)/(niterok+1e-5);
      }
      extern psf_child_started;
      // give the go for next batch:
      shm_write,shmkey,"loop_counter",&([loopCounter]);
      shm_write,shmkey,"mircube",&mircube;
      if (smdebug) write,"main: giving trigger to PSF child";
      sem_give,semkey,3;
      psf_child_started = 1;
    } else {
      // compute integrated phases and fill phase cube
      for (jl=1;jl<=target._nlambda;jl++) {
        for (jt=1;jt<=target._ntarget;jt++) {
          cubphase(,,jt)  = get_phase2d_from_dms(jt,"target") +    \
                            get_phase2d_from_optics(jt,"target") + \
                            get_turb_phase(i,jt,"target");
          // vibration already added to dm1
        }
        if (*target._pupil == []){
          // compute image cube from phase cube
          status = _calc_psf_fast(&pupil,&cubphase,&im,2^dimpow2,
                                  target._ntarget,float(2*pi/(*target.lambda)(jl)),1n);
        } else {
          for (nt=1;nt<=target._ntarget;nt++){
            imtarg = im(,,nt);
            status = _calc_psf_fast(&(*target._pupil)(,,nt),&cubphase(,,nt),&imtarg,2^dimpow2,1,float(2*pi/(*target.lambda)(jl)),1n);
            im(,,nt) = imtarg;
          }
        }
        // Accumulate statistics:
        imav(,,,jl) = imav(,,,jl) + im;
      }
      niterok += 1;
      grow,itv,i;
      if (disp_strehl_indice) sind=disp_strehl_indice; else sind=1;
      grow,strehlsp,im(max,max,sind)/sairy(sind);
      grow,strehllp,imav(max,max,sind,0)/sairy(sind)/(niterok+1e-5);
    }
  }

  time(6) += tac();

  // Displays
  if (disp) {
    for (ns=1;ns<=nwfs;ns++ ) {
      // if cyclecounter = 1, a final image just got computed
      if (wfs(ns).type=="zernike") continue;
      if (wfs(ns)._cyclecounter == wfs(ns).nintegcycles) {
        if (wfs(ns).type!="hartmann") *wfs(ns)._dispimage = *wfs(ns)._fimage;
        // for hartmann, done within shwfs()
      }
    }
  }

  // testing display by external yorick process. Proof of concept.
  if ((display_offline)&&(shm_init_done)) {
    // nothing to do right now as we use variables shm_var'ed from
    // sh_wfs + loopCounter already in shm
  }

  if (okdisp) {
    if (!animFlag) fma;
    // PSF Images
    plt,escapechar(swrite(format="YAO %s | %s | %s",aoSimulVersion,parprefix,sim.name)),0.01,0.227,tosys=0;

    txpos = *target.xposition;
    if (target.xspeed) txpos += *target.xspeed*loopCounter*loop.ittime;
    typos = *target.yposition;
    if (target.yspeed) typos += *target.yspeed*loopCounter*loop.ittime;

    if (dispImImav) {
      disp2d,imav,txpos,typos,1,power=0.5;
      mypltitle,"Average PSFs",[0.,-0.00],height=12;
    } else {
      disp2d,im,txpos,typos,1,power=0.5;
      mypltitle,"Instantaneous PSFs",[0.,-0.00],height=12;
    }
    for (j=1;j<=nwfs;j++) {
      plg,wfs(j).gspos(2),wfs(j).gspos(1),marker='\2',
        type="none",marks=1,color="red";
    }
    myxytitles,"","arcsec",[0.005,0.01],height=12;

    // WFS spots
    if (!allof(wfs.shmethod ==1)) {
      if (wfs_display_mode=="spatial") {
        disp2d,wfs._dispimage,wfs.pupoffset(1,),wfs.pupoffset(2,),2;
        mypltitle,"WFSs (spatial mode)",[0.,-0.005],height=12;
      } else {
        disp2d,wfs._dispimage,wfs.gspos(1,),wfs.gspos(2,),2;
        mypltitle,"WFSs",[0.,-0.00],height=12;
      }
    }

    // mirror surface
    plsys,3;
    limits,square=1;
    if (mergedms4disp) { //e.g. GMT case
      tmp = mircube(,,sum);
      pli,(tmp-min(tmp))*ipupil,1,0.,2,1.;
      range,0.25,0.75;
    } else {
      pli,(mircube(,,1)-min(mircube(,,1)))*ipupil,1,0.,2,1.;
      for (nm=2;nm<=ndm;nm++) {
        if (dm(nm).alt==0) {
          pli,(mircube(,,nm)-min(mircube(,,nm)))*ipupil,nm,0.,nm+1,1.;
        } else {
          pli,(mircube(,,nm)-min(mircube(,,nm))),nm,0.,nm+1,1.;
        }
      }
    }
    myxytitles,"","DM(s)",[0.005,0.],height=12;

    // Strehl plots
    if (anyof(itv)) {
      plsys,4;
      if (plot_rps) {
        extern plot_rps_limit_done;
        if (rp_z!=[]) {
          plh,rp_z(2:),indgen(2:rp_nzerns);
          myxytitles,"Z#","Z value [nm]",[0.005,0.008],height=12;
          if (!plot_rps_limit_done) {
            range,-600,600;
            plot_rps_limit_done=1;
          }
        }
      } else {
        if (rolling_strehl_disp_length) {
          ii = clip(numberof(itv)-rolling_strehl_disp_length/clip(loop.stats_every,1,),1,);
        } else ii=1;
        plg,strehlsp(ii:),itv(ii:),marks=0;
        plg,strehllp(ii:),itv(ii:),marks=0,color="red";
        myxytitles,"",swrite(format="Strehl @ %.2f mic",\
                      (*target.lambda)(0)),[0.005,-0.005],height=12;
        range,0.,1.;
        if (disp_strehl_indice) sind=disp_strehl_indice; else sind=1;
        plt,"Over target indices "+print(sind)(1),0.070,0.6075,tosys=0,justify="LT",height=8;
      }
    }
    // Display residual wavefront
    plsys,5;
    pli, (ipupil*residual_phase)(_n1:_n2,_n1:_n2);
    if (rp_rms!=[]) \
      plt,swrite(format="rms=%dnm",long(rp_rms)),0.450,0.579,\
      tosys=0,height=10,justify="LA";
    mypltitle,swrite(format="Residual phase on %s #%d",residual_phase_what,\
                     residual_phase_which),[0.,-0.0],height=12;

    if (user_plot != []) user_plot,i,init=(i==1);  // execute user's plot routine if it exists.
    if (animFlag && (nshots!=0) && (loopCounter<loop.niter-disp) ) fma;
  }

  if (okcscreen) { control_screen,i; }

  time(7) += tac();

  // Fill the circular buffer
  if (is_set(savecb)) {
    mc = 0; // measurement counter
    for (ns=1;ns<=nwfs;ns++){
      cbmes(mc+1:mc+wfs(ns)._nmes,i) = wfsMesHistory(mc+1:mc+wfs(ns)._nmes,wfs(ns)._framedelay);
      mc += wfs(ns)._nmes;
    }
    //cbmes(,i) = wfsMesHistory(,loop.framedelay);
    cbcom(,i) = comvec;
    cberr(,i) = err;
  }

  glt      = 0.02;
  if (i==1) {glt=1.;}

  nowend = tac(2);
  looptime = (nowend-now)*glt + looptime*(1.-glt);
  tottime += (nowend-now);
  remainingTime = float(loop.niter-i)*looptime;
  remainingTimestring = secToHMS(remainingTime);
  iter_per_sec = loopCounter/tottime;

  // Prints out some results:

  if (!go_quiet) {

  if (((looptime <= 2) && ((i % 500) == 1)) ||
      ((looptime > 2) && ((i % 20) == 1)) || (nshots>=0)) {
    if (numberof(*target.xposition)==1) {
      write,"Iter#  Inst.Strehl  Long expo.Strehl  Time Left  it/s";
    } else {
      write,"       Short expo. image  Long expos. image ";
      write,"Iter#  Max.S/Min.S/Avg.S  Max.S/Min.S/Avg.S  Time Left  it/s";
    }
  }

  if ((looptime > 2) || ((i % 50) == 1) || (nshots>=0)) {
    if (numberof(*target.xposition)==1) {
      msg = swrite(format="%5i  %5.3f        %5.3f             %s",
                   i,(im(max,max,)/sairy)(max),
                   (imav(max,max,,0)/sairy)(max)/(niterok+1e-5),
                   remainingTimestring);
      msg1 = swrite(format=" %.1f",iter_per_sec);
      write,msg+msg1;
      header = "Iter#  Inst.Strehl  Long expo.Strehl  Time Left";
      msg = swrite(format="%5i  %5.3f          %5.3f             %s",
                   i,(im(max,max,)/sairy)(max),
                   (imav(max,max,,0)/sairy)(max)/(niterok+1e-5),
                   remainingTimestring);
      gui_message1,header;
      gui_message,msg;
    } else {
      msg = swrite(format="%5i  %5.3f %5.3f %5.3f  %5.3f %5.3f %5.3f  %s",
                   i,(im(max,max,)/sairy)(max),(im(max,max,)/sairy)(min),
                   avg(im(max,max,)/sairy),(imav(max,max,,0)/sairy)(max)/(niterok+1e-5),
                   (imav(max,max,,0)/sairy)(min)/(niterok+1e-5),
                   (imav(max,max,,0)/sairy)(avg)/(niterok+1e-5),remainingTimestring);
      msg1 = swrite(format=" %.1f",iter_per_sec);
      write,msg+msg1;
      header = "Iter#  Inst:Max.S/Min.S/Avg.S  Avg:Max.S/Min.S/Avg.S  Time Left";
      msg = swrite(format="%5i         %5.3f %5.3f %5.3f        %5.3f %5.3f %5.3f  %s",
                   i,(im(max,max,)/sairy)(max),(im(max,max,)/sairy)(min),
                   avg(im(max,max,)/sairy),(imav(max,max,,0)/sairy)(max)/(niterok+1e-5),
                   (imav(max,max,,0)/sairy)(min)/(niterok+1e-5),
                   (imav(max,max,,0)/sairy)(avg)/(niterok+1e-5),remainingTimestring);
      gui_message1,header;
      gui_message,msg;
    }
  }
  }
  time(8) += tac();

  if (nshots==0) {
    if (animFlag && dispFlag) {
      plsys,1;
      animate,0;
    }
    return;
  }

  if (user_end_go != []) status = user_end_go();  // execute user's routine if it exists.

  tic,2; endtime=_nowtime(2); endtime_str = gettime();
//  if (loopCounter<loop.niter) after, 0.001, go;
  if (loopCounter<loop.niter) {
    if (all) goto go_start;
    else set_idler, go;
  } else {
    if (yaopy) pyk,"aoloop_to_end()";
    gui_hide_statusbar1;
    if (animFlag && dispFlag) {
      plsys,1;
      animate,0;
    }
    if ((sim.svipc>>0)&1) sem_take,semkey,1;
    status = after_loop();
    if (user_wrapup) status = user_wrapup()
    notify,swrite(format="%s: %d iterations completed",parprefix,loopCounter)
  }
  maybe_prompt;
}

func reset(void,nmv=)
/* DOCUMENT reset(void,nmv=)
   Flattens/reset all the DMs (or a subset if nmv not void)
   SEE ALSO:
 */
{
  extern command,mircube,wfsMesHistory,dm,ndm;

  if (nmv==[]) nmv=indgen(ndm);

  if (command!=[]) command *=0.0f;

  for (i=1;i<=numberof(nmv);i++) {
    nm = nmv(i);
    *dm(nm)._command *=0.0f;
  }
  mircube *=0.0f; wfsMesHistory *=0.0f;
}


func stop(void)
/* DOCUMENT stop(void)
   Pause the loop
   SEE ALSO: cont, go, restart
 */
{
  write,format="Stopping @ iter=%d\n",loopCounter;
  gui_hide_statusbar1;
  if (animFlag&&dispFlag) {
    plsys,1;
    animate,0;
  }
  set_idler;
}


func cont(void)
/* DOCUMENT cont(void)
   Resume the loop
   SEE ALSO: stop, restart, go
 */
{
  write,"Proceeding...";
  if (animFlag&&dispFlag) {
    plsys,1;
    animate,1;
  }
  set_idler,go;
}


func restart(void)
/* DOCUMENT restart(void)
   Re-initialize all variable and reset loop counter to zero.
   SEE ALSO: go
 */
{
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


func disptoggle(void) { dispFlag=1-dispFlag; }
func whereat(void) { write,format="@ iter=%d\n",loopCounter; }
func do_prompt(void) { write,format="%s ",">"; }

func after_loop(void)
/* DOCUMENT after_loop(void)
   Wraps up when all iterations have been done. Compute various performance
   criteria (Strehl, FWHM, various performance timers, ...)
   SEE ALSO: aoloop, go
 */
{
  extern strehllp,strehsp, itv;
  extern cbmes, cbcom, cberr;
  extern strehl,e50,fwhm;
  extern iter_per_sec;

  savecb = savecbFlag;

  gui_message,swrite(format="Saving results in %s.res (ps,imav.fits)...",YAO_SAVEPATH+parprefix);
  if (sim.verbose>0) \
    write,format="Saving results in %s.res (ps,imav.fits)...\n",\
      YAO_SAVEPATH+parprefix;

  time2 = (time - roll(time,1))/loopCounter;
  timeComments = ["WF sensing","Reset and measurement history handling","cMat multiplication",\
                  "DM shape computation","Target PSFs estimation","Displays",\
                  "Circular buffers, end-of-loop printouts"];
  if (!go_quiet) {
    for (i=2;i<=8;i++) {
      write,format="time(%d-%d) = %5.2f ms  (%s)\n",i-1,i,time2(i)*1e3,timeComments(i-1);}

    write,format="Finished on %s\n",endtime_str;
  }
  // tottime = (endtime - starttime);
  iter_per_sec = loopCounter/tottime;
  if (!go_quiet) write,format="%f iterations/second on average\n",iter_per_sec;

  // Save the circular buffers:
  if (is_set(savecb)) {
    yao_fitswrite,YAO_SAVEPATH+"cbmes.fits",cbmes;
    yao_fitswrite,YAO_SAVEPATH+"cbcom.fits",cbcom;
    yao_fitswrite,YAO_SAVEPATH+"cberr.fits",cberr;
    write,"cbmes, cbcom and cberr are saved.";
    write,"You can run modal_gain_optimization() to optimize and update the gains";
  }

  // End of loop calculations (fwhm, EE, Strehl):
  fairy  = findfwhm(airy,saveram=1);
  if (!no_ee) e50airy= encircled_energy(airy,ee50); else e50airy=0.;
  strehl = imav(max,max,,)/sairy/(niterok+1e-5);

  psize  = (float(sim.pupildiam)/sim._size)*(*target.lambda)/tel.diam/4.848e-3;
  tmp = "Pixel size in images (e.g., imav):";
  for (nl=1;nl<=numberof(psize);nl++){
    grow, tmp, swrite(format="%.3f",psize(nl));
  }
  grow, tmp, "mas";
  if (!go_quiet) write, tmp;

  write,format="\n         lambda   XPos   YPos  FWHM[mas]  Strehl  E50d[mas]%s\n","";
  for (jl=1;jl<=target._nlambda;jl++) {
    for (jt=1;jt<=target._ntarget;jt++) {
      fwhm(jt,jl) = findfwhm(imav(,,jt,jl),psize(jl),saveram=1);
      if (!no_ee) encircled_energy,imav(,,jt,jl),tmp; else tmp=0.;
      e50(jt,jl)  = tmp*psize(jl);

      write,format=
        "Star#%2d   % 5.2f % 6.1f % 6.1f     %6.1f   %.3f     %6.1f\n",
        jt,(*target.lambda)(jl),(*target.xposition)(jt),(*target.yposition)(jt),
        fwhm(jt,jl),strehl(jt,jl),e50(jt,jl);
    }
    if (target._ntarget > 1){
      write,format="Field Avg % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(avg,jl),strehl(avg,jl),e50(avg,jl);
      write,format="Field rms % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(rms,jl),strehl(rms,jl),e50(rms,jl);
    }
  }

  // Some logging of the results in file parprefix+".res":
  f = open(YAO_SAVEPATH+parprefix+".res","a+");
  write,f,format=
    "\n         lambda   XPos   YPos  FWHM[mas]  Strehl  E50d[mas]%s\n","";
  for (jl=1;jl<=target._nlambda;jl++) {
    for (jt=1;jt<=target._ntarget;jt++) {
      write,f,format=
        "Star#%2d   % 6.2f % 6.1f % 5.1f     %6.1f   %.3f     %6.1f\n",
        jt,(*target.lambda)(jl),(*target.xposition)(jt),(*target.yposition)(jt),
        fwhm(jt,jl),strehl(jt,jl),e50(jt,jl);
    }
    write,f,format="Field Avg % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(avg,jl),strehl(avg,jl),e50(avg,jl);
    write,f,format="Field rms % 5.2f                   %6.1f   %.3f     %6.1f\n",
      (*target.lambda)(jl),fwhm(rms,jl),strehl(rms,jl),e50(rms,jl);
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

  // create a header
  header = fitsBuildCard("WAVELENGTH", *target.lambda , "microns");
  grow, header, fitsBuildCard("PIXSIZE", psize , "milliarcsec");
  grow, header, fitsBuildCard("XPOSITION", *target.xposition , "arcsec");
  grow, header, fitsBuildCard("YPOSITION", *target.yposition , "arcsec");

  yao_fitswrite,YAO_SAVEPATH+parprefix+"-imav.fits",imav, header;

  // saved graphics
  window,7,display="",hcp=YAO_SAVEPATH+parprefix+".ps",wait=1,style="work.gs";
  fma;
  hcpon;

  disp2d,im,*target.xposition,*target.yposition,1,zoom=*target.dispzoom,init=1;
  for (jl=1;jl<=target._nlambda;jl++) {
    disp2d,imav(,,,jl),*target.xposition,*target.yposition,1,power=0.5;
  }

  for (j=1;j<=nwfs;j++) {
    plg,wfs(j).gspos(2),wfs(j).gspos(1),marker='\2',
      type="none",marks=1,color="red";
  }
  axisLegend,"arcsec","arcsec";
  mypltitle=parprefix+"/ Average PSF";
  if (strehlsp != []) {
    fma;
    limits;
    limits,square=0;
    plg,strehlsp,itv;
    if (strehllp!=[]) plg,strehllp,itv,color="red";
    myxytitles,"Iterations","Strehl";
    // hcp;
  }
  plt,escapechar(swrite(format="YAO %s | %s | %s",aoSimulVersion,parprefix,sim.name)),0.01,0.227,tosys=0;
  hcpoff;
  hcp_finish;

  gui_message,swrite(format="Dumping results in %s.res (ps,imav.fits)...DONE",YAO_SAVEPATH+parprefix);

  if (user_end_after_loop) status = user_end_after_loop();

  if (curw != -1) {window,curw;}
  //  if (is_set(disp)) {window,style="boxed.gs";}
  return imav;
}

dtor = pi/180.;
radeg = 1./dtor;


//===================
// API compatibility:
//===================
ZernikeWfs            = zernike_wfs;
PyramidWfs            = pyramid_wfs;
ShWfsInit             = shwfs_init;
multWfsIntMat         = mult_wfs_int_mat;
multWfs               = mult_wfs;
wfsCheckPixelSize     = wfs_check_pixel_size;
checkParameters       = check_parameters;
graphicConfig         = graphic_config;
modalGainOptimization = modal_gain_optimization;
MakePupil             = make_pupil;
FindDmModes           = build_dm_modes;
disp2D                = disp2d;
ftcbAoSimul           = ft_cb_ao_simul;
make_curv_wfsSubs     = make_curv_wfs_subs;
ssNoise               = ss_noise;
ShWfs                 = sh_wfs;
MakeCurvWfs           = make_curv_wfs;
CurvWfs               = curv_wfs;
doInter               = do_imat;
prepSVD               = prep_svd;
buildComMat           = build_cmat;
swapScreens           = swap_screens;
getTurbPhaseInit      = get_turb_phase_init;
getTurbPhase          = get_turb_phase;
getPhase2dFromDms     = get_phase2d_from_dms;
getPhase2dFromOptics  = get_phase2d_from_optics;
correctUpLinkTT       = correct_uplink_tt;
splitWfsVector        = split_wfs_vector;
splitDMCommandVector  = split_dm_vector;
MakePztIF             = make_pzt_dm;
MakeEltPztIF          = make_pzt_dm_elt;
MakeKLIF              = make_kl_dm;
MakeZernikeIF         = make_zernike_dm;
MakeDhIF              = make_dh_dm;
MakeBimorphIF         = make_curvature_dm;
MakeTipTiltIF         = make_tiptilt_dm;
projectAnisoIF        = project_aniso_dm;
MakeAnisoIF           = make_aniso_dm;
compDmShape           = comp_dm_shape;
controlScreen         = control_screen;
mcaoRayleigh          = mcao_rayleigh;
progressBar           = progress_bar;
PhaseStructFunc       = phase_struct_func;
CreatePhaseScreens    = create_phase_screens;
generateVKspectrum    = generate_von_karman_spectrum;
generatePhaseWithL0   = generate_phase_with_L0;
plotMTF               = plot_mtf;
plotDphi              = plot_dphi;
userPlot              = user_plot;
