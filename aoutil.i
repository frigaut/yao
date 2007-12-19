/* AOUTIL.I
 *
 * A collection of utility routines to go with yao.i
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: aoutil.i,v 1.4 2007-12-19 19:44:19 frigaut Exp $
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
 * $Log: aoutil.i,v $
 * Revision 1.4  2007-12-19 19:44:19  frigaut
 * - solved a number of bugs and inconsistencies between regular yao call and
 *   GUI calls.
 * - fixed misregistration for curvature systems
 * - change: misregistration entry from the GUI is now in pupil diameter unit,
 *   not in subaperture unit!
 * - changed default efd in c188-bench.par
 *
 * Revision 1.3  2007/12/19 15:45:32  frigaut
 * - implemented yao.conf which defines the YAO_SAVEPATH directory where
 * all temporary files and result files will be saved
 * - modified yao.i and aoutil.i to save in YAO_SAVEPATH
 * - bumped version to 4.2.0
 * - slight changes to GUI (edit conf file)
 *
 * Revision 1.2  2007/12/19 13:18:59  frigaut
 * - explicit message when screens are not present/found
 * - more messages in statusbar
 * - added statusbar1 (that can hide/show) for strehl status header
 *
 * Revision 1.1.1.1  2007/12/12 23:29:12  frigaut
 * Initial Import - yorick-yao
 *
 * Revision 1.9  2004/10/18 21:30:43  frigaut
 * added things relative to inter actuator coupling (MakePztIF),
 * for regular and elt configs.
 * added tests relative to that in checkparameters
 *
 * Revision 1.7  2004/09/14 04:32:56  frigaut
 * several modifs to do with the creation of turbulent phase screens
 * - implemented cosf and sinf which take and return floats to save space
 * - started coding generatePhaseWithL0 in turbulence.i, not finished
 * - modif YORICKEXE in makefiles, just "yorick" did not do it on osx
 * - modifs ok for both veclib and fftw implementations
 *
 * Revision 1.6  2004/08/02 07:10:53  frigaut
 * Added routine getTurbPhaseInitCheckOverflow, which checks for Y
 * interpolation indices larger than the max Y dimension, which would
 * cause problem in getTurbPhase.
 *
 * Revision 1.5  2004/07/29 04:06:50  frigaut
 * added cvs dollar Log in header
 *
 *
 * func graphicConfig(subsystemnum,dmnum)
 * func checkParameters(void)
 * func disp2D(ar,xpos,ypos,area,zoom=,power=,init=,nolimits=)
 * func hysteresis(v,n,first=) 
 * func modalGainOptimization(disp=,update=)
 * func ftcbAoSimul(FrameDelay,gain,dim)
 * func FindDmModes(disp=)
 * func wfsCheckPixelSize(ns,&binindices,&centroidw,printheader=,silent=)
 * func MakePztIF(nm,&def,disp=)
 * func MakeEltPztIF(nm,&def,disp=)
 * func MakeZernikeIF(nm,&def,disp=)
 * func projectAnisoIF(nmaniso,nmlow,nmhigh,disp=)
 * func MakeAnisoIF(nm,&def,disp=)
 * func MakeTipTiltIF(nm,&def,disp=)
 * func MakeBimorphIF(nm,&def,disp=,cobs=)
 * func MakeCurvWfsSubs(ns,dim,pupd,disp=,cobs=)
 * func _map2d(t,dim,cent)
 * func noll(ord)
 * func nollmat(ord)
 * func rotby90(image,rot)
 * func MakePupil(dim,pupd,xc=,yc=,real=,cobs=)
 * func fwhmStrehl(image,ps,lambda,teldiam,cobs,&strehl,&fwhm,&strehlab,&airy,&psf0,
 * func telftoh1(f,u,v)
 * func telftoh2(f,u)
 * func telftog(f,u)
 * func telftot0(f,cobs)
 * func telfto(lambda,dlambda,teldiam,cobs,pixsize,dim,freqc=,npt=,silent=,returnpsf=)
 * func ftcb(te,tcal,tmir,gain,dim,x=)
 * func encircled_energy(image,&ee50,xc=,yc=)
 * func findfwhm(image,psize)
 *
 */


//----------------------------------------------------

func print_struct(structure,strstr) {
  tmp = strParse(strParse(strParse(strJoin(print(structure)),"(")(2),")")(1),",");
  if (!am_subroutine()) {return tmp;} 
  for (i=1;i<=numberof(tmp);i++) {
    write,format=strstr+".%s\n",tmp(i);
  }
}

//----------------------------------------------------

//func xyxy2xxyy(void)
func xxyy2xyxy(void)
{
  offset=0;
  for (n=1;n<=nwfs;n++) {
    grow,reordered,transpose(reform(indgen(wfs(n)._nmes),
                                    [2,wfs(n)._nsub,2]))(*)+offset;
    offset+=wfs(n)._nmes;
  }
  return reordered;
}

//----------------------------------------------------
func getTurbPhaseInitCheckOverflow
/* DOCUMENT func getTurbPhaseInitCheckOverflow
   This routine has the sole purpose of checking the possible "overflow"
   of the Y index in the future calls to getTurbPhase. In any of the
   Y index at which the interpolation is to be done is larger than the
   phase screen(s) Y dimension, then an error is flagged.
   SEE ALSO: getTurbPhaseInit, getTurbPhase
 */
{
  extern screendim;

  dimy = screendim(2);

  if (max(wfsyposcub(max,,)+yposvec(max,)(,-)) > dimy) {
    write,"\nSome of the phase screens are too small (Y dimension) "+
      "for the specified system.";
    write,format="The following WFS/screens will cause a "+
      "index overflow (current Ydim is %d):\n",dimy;
    maxind = wfsyposcub(max,,)+yposvec(max,)(,-);
    w2     = where2(maxind > dimy);
    for (i=1;i<=dimsof(w2)(3);i++) {
      write,format="WFS#%d, screen#%d, Ymax=%f\n",w2(2,i),w2(1,i),
        maxind(w2(1,i),w2(2,i));
    }
    write,"To remedy this situation, you can either:";
    write,"   - use larger phase screens (Y axis) using 'CreatePhaseScreens'";
    write,"   - Modify the extremum Y offsets of your WFSs";
    write,"   - Lower the altitude of the offending atmospheric layer";
    exit;
  }

  if (max(gsyposcub(max,,)+yposvec(max,)(,-)) > dimy) {
    write,"\nSome of the phase screens are too small (Y dimension) "+
      "for the specified system.";
    write,format="The following Perf Star/screens will cause a "+
      "index overflow (current Ydim is %d):\n",dimy;
    maxind = gsyposcub(max,,)+yposvec(max,)(,-);
    w2     = where2(maxind > dimy);
    for (i=1;i<=dimsof(w2)(3);i++) {
      write,format="Perf.Star#%d, screen#%d, Ymax=%f\n",w2(2,i),w2(1,i),
        maxind(w2(1,i),w2(2,i));
    }
    write,"To remedy this situation, you can either:";
    write,"   - use larger phase screens (Y axis) using 'CreatePhaseScreens'";
    write,"   - Modify the Y position of your Perf. Star";
    write,"   - Lower the altitude of the offending atmospheric layer";
    exit;
  }
}
//----------------------------------------------------
func graphicConfig(subsystemnum,dmnum)
/* DOCUMENT func configGraphic(void)
   Plots a graphical representation of the system config,
   per subsystem and per level (altitude)
   subsystemnum and dmnum optional. If not set all subsystems
   and dms are displayed
   SEE ALSO:
 */
{
  t = span(0.,2*pi,100);
  c = ["red","blue","green","cyan","magenta","yellow","black"];
  markers = ['1','2','3','4','5','6','7'];
  //c: 7 elements: limit to 7 dm levels

  plsys,1;
  limits,square=1;
  psize = tel.diam/sim.pupildiam;

  for (nss=1;nss<=max(dm.subsystem);nss++) {
    if ( (subsystemnum != [] ) && (subsystemnum != nss) ) continue;

    fma;
    // one level by one level
    for (i=1;i<=ndm;i++) {
      if ( (dmnum != [] ) && (dmnum != i) ) continue;
      if ( dm(i).subsystem != nss ) continue;
      if ( dm(i).type == "aniso" ) continue;

      // radius of outmost actuator in meters:
      //      rad = (dm(i).nxact-1)/2.*dm(i).pitch*psize;
      // plots circle containing all actuators
      //      plg,rad*sin(t),rad*cos(t),type=2;

      // radius of outmost object beam in meter:
      tmp = sqrt((*target.xposition)^2.+(*target.yposition)^2.);
      maxobjrad = max(tmp);
      rad = max(tmp)*4.848e-6*dm(i).alt+tel.diam/2.;
      // plots circle containing all rays from objects:
      plg,rad*sin(t),rad*cos(t),type=3;

      if ( dm(i).type == "stackarray" ) {

        // build array containing (x,y) location of valid/controlled actuators:
        loc = ([*(dm(i)._x),*(dm(i)._y)]-sim._cent)*psize;
        // and plot:
        plg,loc(,2),loc(,1),type=0,marker=markers(i),color=c(i);

        if (numberof(*dm(i)._ex) != 0) {
          // build array containing (x,y) location of extrapolated actuators:
          loc = ([*(dm(i)._ex),*(dm(i)._ey)]-sim._cent)*psize;
          // and plot:
          plg,loc(,2),loc(,1),type=0,marker=markers(i);
        }
      }

      // loop on sensor:
      for (j=1;j<=nwfs;j++) {
        if ( wfs(j).subsystem != nss ) continue;
        // computes scaling factor due to limited altitude:
        if (wfs(j)._gsalt != 0) {
          fact = (wfs(j)._gsalt-dm(i).alt)/wfs(j)._gsalt;
        } else {
          fact=1;
        }
        offsets = wfs(j).gspos*4.848e-6*dm(i).alt;
        rad=tel.diam/2.*fact;
        // plotting beam shape for this wfs at this dm altitude:
        plg,offsets(2)+rad*sin(t),offsets(1)+rad*cos(t),color=c(i),width=3;
      }

      myxytitles,"Beam outline [m]","Beam outline [m]",[0.02,0.02];
      mypltitle,swrite(format="Patches of guide star beams on DM#%d, Subsystem %d",i,nss);
      comment = swrite(format="Configuration actuator/beams in a plane at altitude %.0f\n"+
                       "Solid lines: Outline of this subsystem GS beams\n"+
                       "Dotted line: Circle including outermost ray to specified science \n"+
                       "   objects (at %.1f arcsec off-center)\n"+
                       "Numbers for stackarrays actuators: colored:controlled, BW: extrapolated\n"+
                       "# of active actuators= %d, # of extrapolated actuators=%d",
                       dm(i).alt,maxobjrad,dm(i)._nact,dm(i)._enact);
      plt,comment,0.06,0.22,tosys=0,justify="LT";
      plt,sim.name,0.01,0.01,tosys=0;
      limits; lim=limits(); limits,lim(1)*1.1,lim(2)*1.1,lim(3)*1.1,lim(4)*1.1;
      hcp;
      if (sim.debug >=1) typeReturn;
      fma;
    }
    if ( subsystemnum != [] ) continue;


    // all levels
    for (i=1;i<=ndm;i++) {
      if ( dm(i).subsystem != nss ) continue;
      if ( dm(i).type == "aniso" ) continue;
      // radius of outmost actuator in meters:
      rad = (dm(i).nxact-1)/2.*dm(i).pitch*psize;
      // plots circle containing all actuators
      plg,rad*sin(t),rad*cos(t),type=i;
      if ( dm(i).type == "stackarray" ) {
        // build array containing (x,y) location of valid actuators:
        loc = ([*(dm(i)._x),*(dm(i)._y)]-sim._cent)*psize;
        // and plot:
        plg,loc(,2),loc(,1),type=0,marker=markers(i),color=c(i);
      }
      // loop on sensor:
      for (j=1;j<=nwfs;j++) {
        if ( wfs(j).subsystem != nss ) continue;
        // computes scaling factor due to limited altitude:
        if (wfs(j)._gsalt != 0) {
          fact = (wfs(j)._gsalt-dm(i).alt)/wfs(j)._gsalt;
        } else {
          fact=1;
        }
        offsets = wfs(j).gspos*4.848e-6*dm(i).alt;
        rad=tel.diam/2.*fact;
        // plotting beam shape for this wfs at this dm altitude:
        plg,offsets(2)+rad*sin(t),offsets(1)+rad*cos(t),color=c(i),width=3;
      }
    }

    myxytitles,"Beam outline [m]","Beam outline [m]",[0.02,0.02];
    mypltitle,swrite(format="Patches of guide star beams on DMs, Subsystem %d",nss);
    plt,sim.name,0.01,0.01,tosys=0;
    limits; lim=limits(); limits,lim(1)*1.1,lim(2)*1.1,lim(3)*1.1,lim(4)*1.1;
    hcp;
    if (sim.debug >=1) typeReturn;
  }
}

//----------------------------------------------------
func checkParameters(void)
/* DOCUMENT func checkParameters(void)
   Check the parameters in yao parfile
   - set defaults
   - value valid?
   - compatibility with other parameters
   SEE ALSO:
 */
{
  write,format="Checking parameters ... %s\n","";
  
  //==================================================================  
  // BASIC EXISTENCE AND CONSISTENCY CHECKINGS AND DEFAULT ASSIGNMENTS
  //==================================================================  

  if (nwfs == []) {exit,"nwfs has not been set";}
  if (ndm == []) {exit,"ndm has not been set";}
  
  // sim structure
  if (sim.pupildiam == 0) {exit,"sim.pupildiam has not been set";}

  // atm structure
  if ((*atm.screen) == []) {exit,"atm.screen has not been set";}
  if (typeof(*atm.screen) != "string") {exit,"*atm.screen is not a string";}

  ftmp = *atm.screen;
  for (i=1;i<=numberof(ftmp);i++) {
    if (!open(ftmp(i),"r",1)) { // file does not exist
      msg = swrite(format="Phase screen %s not found!\\n"+
                   "You need to generate phase screens for yao.\\n"+
                   "Go to \"Phase Screen -> Create phase Screen\"\\n"+
                   "If you already have phase screens, you can modify\\n"+
                   "the path in the parfile atm.screen definition",ftmp(i));                
      if (_pyk_proc) pyk,"set_cursor_busy(0)";
      pyk_warning,msg;
      msg = swrite(format="\nWARNING: %s not found!\nEdit the par file and change the "+
                   "path,\nand/or run \"Phase Screen -> Create phase Screen\"\n",ftmp(i));
      write,msg;
      break;
    }
  }

  if ((*atm.layerfrac) == []) {exit,"atm.layerfrac has not been set";}
  if ((*atm.layerspeed) == []) {exit,"atm.layerspeed has not been set";}
  if ((*atm.layeralt) == []) {exit,"atm.layeralt has not been set";}
  if ((*atm.winddir) == []) {exit,"atm.winddir has not been set";}

  if (nallof(_(numberof(*atm.screen),numberof(*atm.layerfrac),numberof(*atm.layerspeed), \
               numberof(*atm.layeralt)) == numberof(*atm.winddir))) {
    exit,"Some elements within atm.screen, layerfrac, layerspeed, layeralt \n"+\
      "or winddir do not have the same number of elements.";
  }

  // wfs structure
  for (ns=1;ns<=nwfs;ns++) {
    if (wfs(ns).type == string()) {
      exit,swrite(format="wfs(%d).type has not been set",ns);
    }
    if (wfs(ns).subsystem == 0) {wfs(ns).subsystem = 1;}
    if (wfs(ns).lambda == 0) {
      exit,swrite(format="wfs(%d).lambda has not been set",ns);
    }
    if (wfs(ns).dispzoom == 0) {wfs(ns).dispzoom = 1;}
    if ((wfs(ns).type == "curvature") && ((*wfs(ns).nsubperring) == [])) {
      write,format="wfs(%d).nsubperring has not been set.\n",ns;
      write,"Valid values are, for example:";
      exit,"wfs(n).nsubperring = &([6,12,18]);";
    }
    if ((wfs(ns).type == "curvature") && (wfs(ns).l == 0)) {
      exit,swrite(format="wfs(%d).l has not been set",ns);
    }
    if ((wfs(ns).type == "curvature") && (wfs(ns).gsdepth > 0)) {
      exit,swrite(format="wfs(%d): gsdepth not implemented for CWFS, sorry",ns);
    }
    if ((wfs(ns).type == "curvature") && (wfs(ns).rayleighflag == 1)) {
      exit,swrite(format="wfs(%d): Rayleigh not implemented for CWFS, sorry",ns);
    }
    if ((wfs(ns).gsalt == 0) && (wfs(ns).rayleighflag == 1)) {
      write,swrite(format="wfs(%d): gsalt = 0 and Rayleigh flag is set !",ns);
      exit,"Rayleigh is not computed for NGS sensors, sorry.";
    }
    if ((wfs(ns).gsalt > 0) && (wfs(ns).laserpower == 0)) {
      exit,swrite(format="wfs(%d): this is a LGS and you haven't set wfs.laserpower",ns);
    }
    if ((wfs(ns).type == "curvature") && (wfs(ns).fieldstopdiam == 0)) {
      wfs(ns).fieldstopdiam = 1.f;
    }
    if ((wfs(ns).type == "hartmann") && (wfs(ns).shmethod == 0)) {
      exit,swrite(format="wfs(%d).shmethod has not been set",ns);
    }
    if ((wfs(ns).type == "hartmann") && (wfs(ns).shnxsub == 0)) {
      exit,swrite(format="wfs(%d).shnxsub has not been set",ns);
    }
    if ((wfs(ns).type == "hartmann") && (wfs(ns).shmethod == 2)) {
      if ((wfs(ns).type == "hartmann") && (wfs(ns).pixsize == 0)) {
        exit,swrite(format="wfs(%d).pixsize has not been set",ns);
      }
      if ((wfs(ns).type == "hartmann") && (wfs(ns).npixels == 0)) {
        exit,swrite(format="wfs(%d).npixels has not been set",ns);
      }
    }
    if (wfs(ns).nintegcycles == 0) {wfs(ns).nintegcycles = 1;}
    if (wfs(ns).fracIllum == 0) {wfs(ns).fracIllum = 0.5;}
    if (wfs(ns).optthroughput == 0) {wfs(ns).optthroughput = 1.0;}
    wfs.ron = float(wfs.ron);
  }

  // dm structure
  for (nm=1;nm<=ndm;nm++) {
    if (dm(nm).type == string()) {
      exit,swrite(format="dm(%d).type has not been set",nm);
    }
    if (dm(nm).subsystem == 0) {dm(nm).subsystem = 1;}
    if (dm(nm).iffile == string()) {dm(nm).iffile = "";}
    if ((dm(nm).type=="stackarray") && (dm(nm).coupling==0)) {
      dm(nm).coupling = 0.2;
      write,format="dm(%d).coupling set to %f\n",nm,dm(nm).coupling;
    }
    if (dm(nm).ecmatfile == string()) {dm(nm).ecmatfile = "";}
    if (dm(nm).push4imat == 0) {dm(nm).push4imat = 20;}
    if (dm(nm).thresholdresp == 0) {dm(nm).thresholdresp = 0.3;}
    if (dm(nm).gain == 0) {dm(nm).gain = 1.;}
    if (dm(nm).unitpervolt == 0) {
      if ( (dm(nm).type == "tiptilt") || (dm(nm).type == "zernike") ){
        dm(nm).unitpervolt = 0.0005;
      } else if (dm(nm).type == "bimorph") {
        dm(nm).unitpervolt = 1.0;
      } else {
        dm(nm).unitpervolt = 0.01;
      }
    }
    if ( (dm(nm).type == "bimorph") && ((*dm(nm).nelperring) == []) ) {
      write,format="dm(%d).nelperring has not been set.\n",nm;
      write,"Valid values are, for example:";
      exit,"dm(n).nelperring = &([6,12,18]);";
    }
    if ( (dm(nm).type == "stackarray") && (dm(nm).nxact == 0) ) {
      exit,swrite(format="dm(%d).nxact has not been set",nm);
    }
    if ( (dm(nm).type == "stackarray") && (dm(nm).pitch == 0) ) {
      exit,swrite(format="dm(%d).picth has not been set",nm);
    }
    if (dm(nm).type=="stackarray") {
      if ((dm(nm).coupling<0.04) || (dm(nm).coupling>0.32)) {
        write,format="Invalid value for dm(%d).coupling -> %f\n",nm,dm(nm).coupling;
        exit,"Valid values from 0.04 to 0.30";
      }
    }
    if ( (dm(nm).type == "zernike") && (dm(nm).nzer == 0) ) {
      exit,swrite(format="dm(%d).nzer has not been set",nm);
    }    
  }

  // mat structure
  if ((*mat.condition) == []) {exit,"mat.condition has not been set";}
  if (numberof(*mat.condition) != max(_(wfs.subsystem,dm.subsystem)) ) {
    exit,"dimension of *mat.condition is not equal to the number of subsystems";
  }
  if (mat.file == string()) {mat.file = "";}

  // tel structure
  if (tel.diam == 0) {exit,"tel.diam has not been set";}

  // target structure
  if ((*target.lambda) == []) {exit,"target.lambda has not been set";}
  if ((*target.xposition) == []) {exit,"target.xposition has not been set";}
  if ((*target.yposition) == []) {exit,"target.yposition has not been set";}
  if ((*target.dispzoom) == []) {
    target.dispzoom = &(array(1.,numberof(*target.lambda)));
  }
  if (nallof(_(numberof(*target.xposition), \
               numberof(*target.yposition)) == numberof(*target.dispzoom) )) {
    exit,"Some elements within target.xposition, yposition, dispzoom "+\
      "do not have the same number of elements.";
  }

  // gs structure
  if (anyof(wfs.gsalt != 0) && (gs.lgsreturnperwatt == 0)) {
    gs.lgsreturnperwatt = 22.;
    write,format="gs.lgsreturnperwatt set to %f\n",gs.lgsreturnperwatt;
  }
  if (anyof(wfs.gsalt == 0) && (gs.zeropoint == 0)) {
    exit,"You have some NGS and gs.zeropoint has not been set";
  }
  
  // loop structure
  if (loop.gain == 0) {exit,"loop.gain has not been set";}
  if (loop.niter == 0) {exit,"loop.niter has not been set";}
  if (loop.ittime == 0) {exit,"loop.ittime has not been set";}
  if (loop.startskip == 0) { loop.startskip = 10; }
  if (loop.skipevery == 0) { loop.skipevery = loop.niter; }
  if (loop.skipby == 0) { loop.skipby = 10000; }
  if (loop.modalgainfile == string()) {loop.modalgainfile = "";}
  if (loop.stats_every==0) loop.stats_every=4;
  

  //============================================================================  
  // DONE WITH BASIC EXISTENCE AND CONSISTENCY CHECKINGS AND DEFAULT ASSIGNMENTS
  //============================================================================  

  // NOW GOING INTO MORE ELABORATE CHECKINGS
  // WFSs:
  for (ns=1;ns<=nwfs;ns++) {
    if (wfs(ns).lambda < 0.1) {
      write,format="WFS#%d: wfs.lambda < 0.1. That seems weird.\n",ns;
      write,"Remember: lambda should be in microns";
      exit;
    }
    if (wfs(ns).shthreshold < 0.) {
      write,format="WFS#%d: wfs.shthreshold < 0 does not make sense\n",ns;
      exit;
    }
    if ((wfs(ns).filtertilt == 1) && (wfs(ns).correctUpTT == 0)) {
      write,format="WARNING! WFS#%d: wfs.correctUpTT = 0 with wfs.filtertilt = 1\n",ns;
    }
    if ((wfs(ns).correctUpTT == 1) && (wfs(ns).uplinkgain == 0)) {
      write,format="WARNING! WFS#%d: wfs.correctUpTT = 1 but wfs.uplinkgain = 0\n",ns;
    }
    // for shmethod = 1, we still need to set npixels and pixsize to avoid crashing the
    // WFS initialization routine that is still used to set other parameters used with
    // method 1. So we put DUMMY VALUES.
    if (wfs(ns).shmethod == 1) {
      wfs(ns).pixsize = 0.1;
      wfs(ns).npixels = 2;
    }
    if ( (wfs(ns).type=="curvature") && (wfs(ns).nintegcycles != 1) ) {
      write,"\n I have not implemented yet nintegcycle > 1 for curvature WFS.";
      write,"I can see no reason to do it, as curvature sensors are usually";
      write,"read-out noise free, therefore reducing the gain should be";
      write,"prefered. If you have a compelling reason and want it, drop";
      write,"me an email.";
      exit;
    }
  }

  // Sets the Influence function file name if not set:
  for (nm=1;nm<=ndm;nm++) {
    if (dm(nm).iffile == "") {
      dm(nm).iffile = parprefix+swrite(format="-if%d",nm)+".fits";
      write,format="dm(%d).iffile set to %s\n",nm,dm(nm).iffile;
    }
    if (dm(nm).ecmatfile == "") {
      dm(nm).ecmatfile = parprefix+swrite(format="-ecmat%d",nm)+".fits";
    }
    if ((dm(nm).type != "stackarray") && (dm(nm).elt == 1)) {
      exit,swrite(format="DM %d: parameter dm.elt only used with stackarray mirrors\n",nm);
    }
    dm(nm)._eiffile = parprefix+swrite(format="-if%d",nm)+"-ext.fits";
  }

  if (anyof(dm.elt == 1) && anyof(dm.type == "aniso")) {
    exit,"You can not use currently dm.elt=1 with anisoplanatism modes";
  }
  
  for (i=1;i<=max(_(wfs.subsystem,dm.subsystem));i++) {
    if (noneof(wfs.subsystem == i)) {
      exit,swrite(format="There is no WFS in subsystem %d\n",i);
    }  
    if (noneof(dm.subsystem == i)) {
      exit,swrite(format="There is no DM in subsystem %d\n",i);
    }  
  }

  // Sets the interaction/command matrix file name if not set:
  if (mat.file == "") {mat.file = parprefix+"-mat.fits";}

  // # of targets (stars at which the performance are evaluated):
  target.xposition = &(float(*target.xposition));
  target.yposition = &(float(*target.yposition));
  target._ntarget = numberof(*target.xposition);
  target._nlambda = numberof(*target.lambda);
  if (anyof((*target.lambda) < 0.1)) {
      write,"Some or all target.lambda < 0.1. That seems weird.";
      write,"Remember: lambda should be in microns";
      exit;
    }

  wfs.type = strtolower(wfs.type);
  dm.type  = strtolower(dm.type);

  // Are we using a WFS we know?
  if (nallof((wfs.type == "curvature") | (wfs.type == "hartmann") |
             (wfs.type =="pyramid"))) {
    exit,"wfs.type : Unknown value";
  }

  // Are we using a DM we know?
  if (nallof((dm.type == "bimorph") | (dm.type == "stackarray") |
             (dm.type == "zernike") | (dm.type == "tiptilt") |
             (dm.type == "aniso"))) {
    exit,"dm.type : Unknown value";
  }

  if (opt!=[]) {
    if (opt.misreg==[]) {
      noptics = numberof(opt.phasemaps);
      for (i=1;i<=noptics;i++) opt.misreg = [0.,0.];
    }
    opt.misreg= float(opt.misreg);
  }
  
  write,"OK";
}
//----------------------------------------------------
func disp2D(ar,xpos,ypos,area,zoom=,power=,init=,nolimits=)
/* DOCUMENT func disp2D(arrayptr,xpos,ypos,area,zoom=,power=,init=)
   display several images in the same plsys, at position given
   by xpos and ypos.
   ar: ar can be either an array of pointers or an image cube
   xpos,ypos: the (X,Y) positions, in arbitrary coordinates
   area: plsys number
   zoom: keyword (vector, dim nim), additional zooming factor (on top of default)
   power: *arrayptr^power is displayed
   init: initialize, precompute stuff
   SEE ALSO:
 */
{
  extern basezoomptr;
  plsys,area;
  
  if (typeof(ar) == "pointer") {
    cas = "ptr";
    nim = numberof(ar);
    earlyExit = (*ar(1) == []); // in case we call init with ar undefined (legitimate)
  } else {
    cas = "cube";
    tmp = dimsof(ar);
    if (tmp(1) == 2) {nim=1;} else {nim=tmp(4);}
  }
  
  if (is_set(init)) {
    if (basezoomptr == []) {basezoomptr=array(pointer,10);}
    if ((zoom != []) && (numberof(zoom) != nim)) {zoom = array(zoom,nim);}
    if (zoom == []) {zoom = array(1.,nim);}
    xd = abs(xpos-xpos(-,));
    yd = abs(ypos-ypos(-,));
    di = sqrt(xd^2.+yd^2.);
    di = di+unit(nim)*max(di);
    basezoom = (1.+0.9*min(di)/2.)*zoom;
    basezoomptr(area) = &basezoom;
    if (!is_set(nolimits)) {
      limits,min(xpos-basezoom),max(xpos+basezoom),
        min(ypos-basezoom),max(ypos+basezoom),square=1;
    }
    if (earlyExit) return;
  }

  basezoom = *basezoomptr(area);

  if ( cas=="ptr") {
    if (!is_set(power)) {
      for (i=1;i<=nim;i++) {
        pli,*ar(i),xpos(i)-basezoom(i),ypos(i)-basezoom(i),
          xpos(i)+basezoom(i),ypos(i)+basezoom(i);
      }
    } else {
      for (i=1;i<=nim;i++) {
        pli,*ar(i)^power,xpos(i)-basezoom(i),ypos(i)-basezoom(i),
          xpos(i)+basezoom(i),ypos(i)+basezoom(i);
      }
    }
  } else {
    if (!is_set(power)) {
      for (i=1;i<=nim;i++) {
        pli,ar(,,i),xpos(i)-basezoom(i),ypos(i)-basezoom(i),
          xpos(i)+basezoom(i),ypos(i)+basezoom(i);
      }
    } else {
      for (i=1;i<=nim;i++) {
        pli,ar(,,i)^power,xpos(i)-basezoom(i),ypos(i)-basezoom(i),
          xpos(i)+basezoom(i),ypos(i)+basezoom(i);
      }
    }
  }
  
  if ((is_set(init)) & (!is_set(nolimits))) {
    limits,square=1;
    limits;
  }
}
//--------------------------------------------------------------------------
func hysteresis(v,n,first=) 
{

  // utilisation :
  // appel : pos = hysteresis(v) pour avoir la position des actuateurs
  // pour un vecteur de controle v, pris en compte l"hysteresis
  //
  // v = vecteur des voltages
  // hyst = fraction hysteresis (e.g. 0.05 = 5%)
  // first = vecteur v de depart
  //
  // a bit cumbersome. borrowed from IDL routine of same name in simullgs.pro

  if (!is_void(first)) {
    dm(n)._vold = &(first*1.f);
    dm(n)._posold = &(first*1.f);
    dm(n)._chpos = &(first*1.f);
    dm(n)._chv = &(first*1.f);
    dm(n)._dir = &(first*0.f);
    dm(n)._delta = &(first*0.f);
    return first;
  } else {
    vold = *dm(n)._vold;
    posold = *dm(n)._posold;
    chpos = *dm(n)._chpos;
    chv = *dm(n)._chv;
    dir = *dm(n)._dir;
    delta = *dm(n)._delta;
    hyst = dm(n).hyst;
  }

  chdir  = ((v-vold)*dir) <= 0;  // change direction ?
  wchd  = where(chdir);  // indices of which that changed dir
  if (anyof(chdir)) {
    delta(wchd) = -chv(wchd)+vold(wchd);  // delta v since last direction change
    chpos(wchd) = posold(wchd);  // new position where direction change occured
    chv(wchd) = vold(wchd);
  }

  pos = v;
  w = where(delta > 0.);  // was going up
  if (anyof(delta > 0.)) {pos(w) = min(v(w)+hyst/2.*delta(w),chpos(w));}
  w = where(delta < 0.);
  if (anyof(delta < 0.)) {pos(w) = max(v(w)+hyst/2.*delta(w),chpos(w));}


  dir = v-vold;
  vold = v;
  posold = pos;

  dm(n)._vold = &(vold);
  dm(n)._posold = &(posold);
  dm(n)._chpos = &(chpos);
  dm(n)._chv = &(chv);
  dm(n)._dir = &(dir);
  dm(n)._delta = &(delta);
  
  return pos;
}
//--------------------------------------------------------------------------

func modalGainOptimization(disp=,update=)
/* DOCUMENT func modalGainOptimization(disp=,update=)

   NOT UPGRADED TO VERSION 2.
   DO NOT USE.

   This routine optimizes the modal gains.
   Keywords:
   disp: set to display stuff
   update: set if this is a gain update (opposite to the first time run)

   This routine uses:
     - the saved error circular buffer (cberr.fits)

   This routine calls:
   ...
   This routine sets:
   ...
   
   SEE ALSO:
 */
{
  cberr        = fits_read("cberr.fits"); // CB of actuator error
  nGoodSamples = long(2^floor(log(ao.LoopNIter-ao.LoopStartSkip)/log(2)));
  cbmoderr     = atm(,+)*cberr(+,1-nGoodSamples:);// CB of mode coef errors
  modgains     = modalgain*ao.LoopGain; // overall mode gains
  bestGains    = modgains*0.;
  length       = 2048;
  gainNpt      = 30;
  toGains      = spanl(1e-2,0.9,gainNpt);
  errvect      = array(float,gainNpt);
  errTF        = array(float,length/2-1,gainNpt);

  fudgeFactor  = 0.85;

  // initialize Error Transfer functions:
  for (n=1;n<=gainNpt;n++)
    {errTF(,n) = ftcbAoSimul(ao.LoopFrameDelay,toGains(n),length/2)(,1);}
  
  
  // Loop over controlled modes:
  for (i=1;i<=NModesControlled;i++)
  {
    // Build PSD of the error
    psderr = psd(cbmoderr(i,),length,samp=ao.LoopItTime,filter=6,silent=1)(2:,2);
    pause,10;
    // ^^ psderr is in fact length/2 long

    // Build current gain transfer function
    curErrTF = ftcbAoSimul(ao.LoopFrameDelay,modgains(i),length/2)(,1);

    // Reconstruct open loop PSD
    OLmodPsd = psderr/curErrTF;
    
    // Loop over possible gains
    for (n=1;n<=gainNpt;n++) {errvect(n) = sum(OLmodPsd*errTF(,n));}
    if (is_set(disp)) {logxy,0,0; plot,errvect,toGains; pause,500;}

    // Determine best gain
    bestGains(i) = fudgeFactor*toGains(errvect(mnx));
    if (ao.verbose>=1) {write,format="Mode %i, gain to be updated from %f to %f\n",
                       i,modgains(i),bestGains(i);}
  }
  // Update modalgain
  if (is_set(update))
  {
    modalgain = bestGains/ao.LoopGain;
    fitsWrite,YAO_SAVEPATH+ao.LoopModalGainFile,modalgain;
    if (ao.verbose>=1) {write,format="Gains updated and saved in !",ao.LoopModalGainFile;}
  }
}

//----------------------------------------------------

func ftcbAoSimul(FrameDelay,gain,dim)
/* DOCUMENT func ftcbAoSimul(FrameDelay,gain,dim)

   NOT UPGRADED TO VERSION 2.
   DO NOT USE UNTIL UPGRADED.

   This routine simulates the time aspect of the numerical loop
   and compute the associated transfer functions (error, closeloop)
   for further use in modalGainOptimization().
   Inputs:
   FrameDelay: frame delay in close loop.
   gain: AO loop gain
   dim: desired linear size of the output transfer functions
   SEE ALSO: modalGainOptimization
 */
{
  input         = array(float,2*dim);
  mirvect       = array(float,2*dim);
  output        = array(float,2*dim);
  input(2)      = 1.;
  mir           = 0.;
  command       = 0.;
  wfsMesHistory = array(float,FrameDelay+1);

  for (i=1;i<=2*dim;i++)
  {
    in = input(i);
    // Subtract the mirror shape:
    out = in - mir;

    // Do the WFS:
    mes = out;

    // Handling frame delay (0 -> no frame delay at all, 1 -> regular
    // one frame delay of integrator, 2 -> integrator + one frame
    // computation, etc...)
    wfsMesHistory  = roll(wfsMesHistory,-1);
    wfsMesHistory(FrameDelay+1) = mes;
    usedMes        = wfsMesHistory(1);
    
    // computes the command vector (integrator with gain):
    err      = usedMes;
    command += gain*err;
    
    // Computes the mirror shape using influence functions:
    mir      = command;
    // Subtract the mirror shape:
    out = in - mir;

    mirvect(i) = command;
    output(i)  = err;
    //    write,input(i),command,err,out; pause,500;
  }
  errTF = abs((fft(output+1e-12,1))^2.)/abs((fft(input,1))^2.);
  errTF = errTF(2:dim);
  clTF  = abs((fft(mirvect,1))^2.)/abs((fft(input,1))^2.);
  clTF  = clTF(2:dim);
  
  return [errTF,clTF];
}

//----------------------------------------------------
func ssNoise(nss)
{

  write,format="Computing propagated noise for subsystem %d\n",nss;
  write,"Extracting cmat";
  sswfs = ssdm = [];
  for (ns=1;ns<=nwfs;ns++) {grow,sswfs,array(wfs(ns).subsystem,wfs(ns)._nmes);}
  for (nm=1;nm<=ndm;nm++)  {grow,ssdm,array(dm(nm).subsystem,dm(nm)._nact);}
  wsswfs = where(sswfs == nss);
  wssdm  = where(ssdm  == nss);
  cmat =  cMat(wssdm,wsswfs);

  write,"Computing subsystem modes";
  nmv = where(dm.subsystem == nss);
  if (numberof(nmv)==0) error,"No DM found in subsystems";
  modes = FindDmModes(nmv,actmodes,modvar);
  tmp = actmodes(+,)*cmat(+,);
  modcov = trace(tmp(,+)*tmp(,+));
  // now the measurements are in arcsec, so we got to convert to rd of
  // phase difference at the edge of the subaperture, which is what
  // the SH measure:
  nsv = where(wfs.subsystem == nss);
  if (numberof(nsv) != 1) error,"zero or more than one wfs in subsystem";
  ns = nsv(1);
  subsize = tel.diam/wfs(ns).shnxsub;
  arcsectord = 4.848e-6*subsize*2*pi/wfs(ns).lambda/1e-6;
  noise = modcov*modvar/arcsectord^2.;
  write,format="Total noise on phase for 1rd2 per subaperture = %f rd2\n",sum(noise);
  fma;plh,noise; limits,square=0; limits;
  return noise;
}
//----------------------------------------------------

func FindDmModes(nm,&u,&modvar,&eigenvalues,disp=)

  /* DOCUMENT unfinished, I think.
     NOT FINISHED, NOT UPGRADED TO VERSION 2.
     DO NOT USE.
   */
{
  if (anyof(dm(nm).elt)) exit,"Not implemented for dm.elt=1";
  def = [];
  for (i=1;i<=numberof(nm);i++) {
    grow,def,*dm(nm(i))._def;
  }
  defpup = ipupil(dm(nm(1))._n1:dm(nm(1))._n2,dm(nm(1))._n1:dm(nm(1))._n2);
  wpup = where(defpup);
  tmp = def(*,)(wpup,); // not really saving RAM, but...
  tpup	= sum(defpup);
  /*
    p	= (def*ipupil)(sum,sum,)/tpup;
    def	= def-p(-,-,);
    def	= reform(def,long(ao._size)*ao._size,ao._DmNAct);
    def	= def(where(ipupil),);
  */
  dd	= tmp(+,)*tmp(+,);
  eigenvalues = SVdec(dd,u,vt);
  modes = def(,,+)*u(+,);
  // mask with ipupil
  modes = modes*defpup(,,-);
  // compute mode variance:
  modvar = (modes(*,)(wpup,)(rms,))^2.;
  return modes;
  if (disp == 1) 
   {for (i=1;i<=ao._DmNAct;i++) {fma; pli,modes(,,i)*ipupil;pause,100;}}
  //normalize the modes:
  norm	= sqrt((modes^2.*ipupil)(avg,avg,));
  modes	= modes/norm(-,-,);

  ActIF	= modes(,,1:-1) //*(1./indgen(ao._DmNAct-1))(-,-,);
  ao._DmNAct = ao._DmNAct-1;
  doInter,disp=1;
  buildComMat,all=1,nomodalgain=1;
  mn	= cMat(,+)*cMat(,+);
  mn	= diag(mn)*indgen(ao._DmNAct-1)^2.
  if (disp == 1) {fma;plg,mn;}
  return modes;
}

//----------------------------------------------------
func wfsCheckPixelSize(ns,&binindices,&centroidw,printheader=,silent=)
{
/* DOCUMENT wfsCheckPixelSize(ns,&binindices,&centroidw)

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
   cause problems.
   This routine fills and returns binindices and centroidw
   SEE ALSO: ShWfs
 */
  
  if (odd(wfs(ns).npixels)) {
    write,format="WFS#%d: wfs.npixels odd values not supported.\n",ns;
    exit;
  }
  
  pupd	     = sim.pupildiam;
  nxsub	     = wfs(ns).shnxsub(0);
  subsize    = pupd/nxsub;
  sdim       = long(2^ceil(log(subsize)/log(2)+1));
  err        = 0;

  desiredPixelSize = wfs(ns).pixsize;
  desiredNpixels = wfs(ns).npixels;
  quantumPixelSize = wfs(ns).lambda/(tel.diam/sim.pupildiam)/4.848/sdim;
  rebinFactor = round(desiredPixelSize/quantumPixelSize);
  if ( rebinFactor == 0 ) {rebinFactor = 1l;}
  actualPixelSize  = rebinFactor*quantumPixelSize;

  wfs(ns).pixsize = actualPixelSize;

  //  write,format="WFS#%2d Pixel size: Desired = %f, Actual = %dx%f = %f\n",
  //    ns,desiredPixelSize,rebinFactor,quantumPixelSize,actualPixelSize;
  //  write,format="  (Max subaperture size = %f)\n",quantumPixelSize*sdim;

  while ( (rebinFactor*wfs(ns).npixels) > sdim ) {
    wfs(ns).npixels -= 2;
    err = 1;
  }

  if (!is_set(silent)) {
    f	= open(YAO_SAVEPATH+parprefix+".res","a+");
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


  
  //nbin = long(rebinFactor);
  rdim = rebinFactor*wfs(ns).npixels;
  xy = long(indices(rdim)-1.);
  //binxy = rdim/rebinFactor;
  tmp = xy(,,1)/rebinFactor+xy(,,2)/rebinFactor*wfs(ns).npixels;
  binindices = array(-1l,[2,sdim,sdim]);
  binindices(sdim/2-rdim/2+1:sdim/2+rdim/2,sdim/2-rdim/2+1:sdim/2+rdim/2) = tmp;
  binindices = int(eclat(binindices));

  centroidw = indgen(wfs(ns).npixels)-1.-(wfs(ns).npixels/2.-0.5);
  // we might as well express it in arcsec:
  centroidw = float(centroidw*actualPixelSize);

  return err;
}

//----------------------------------------------------

func MakePztIF(nm,&def,disp=)
  /* DOCUMENT function MakePztIF2(dm_structure,disp=)
     the influence functions are in microns per volt.
  */
{
  gui_progressbar_frac,0.;
  gui_progressbar_text,"Computing Influence Functions";
  coupling=dm(nm).coupling;

  // best parameters, as determined by a multi-dimensional fit
  // (see coupling3.i)
  a=[4.49469,7.25509,-32.1948,17.9493];
  p1 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;
  
  a = [2.49456,-0.65952,8.78886,-6.23701];
  p2 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [1.16136,2.97422,-13.2381,20.4395];
  irc = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;
  
  if (sim.debug==2) write,format="p1=%f  p2=%f  ir=%f\n",p1,p2,irc;

  dim   = dm(nm)._n2-dm(nm)._n1+1;
  size	= sim._size;
  nxact	= dm(nm).nxact;
  cobs	= tel.cobs;
  cent  = sim._cent;
  pitch = dm(nm).pitch;
  /*
    ir	= pitch*1.2;
    ir	= pitch*1.46;
    ir	= pitch*1.65;
    c  = 3.8; p1 = 3.9; p2 = 2.4; ir = pitch*1.20;  // good. no. coupling 8% too low
    c  = 3.8; p1 = 4; p2 = 2.4; ir = pitch*1.65;
    c  = 3.75; p1 = 4.2; p2 = 2.5; ir = pitch*1.25; //ok, coupling=13%
    c  = 3.74; p1 = 3.805; p2 = 2.451; ir = pitch*1.4; //good, coupling=17%
    c  = 4; p1 = 3.84; p2 = 2.5; ir = pitch*1.5; //good, coupling=20%
    c  = 3.74; p1 = 3.805; p2 = 2.451; ir = pitch*1.4; //good, coupling=17%
  */
  ir = irc*pitch;
    
  bord  = 0;
  cub   = array(float,nxact+bord*2,nxact+bord*2,4);

  // make X and Y indices array:
  xy    = indices(nxact+bord*2);

  // express "centered" coordinate of actuator in pixels:
  xy    = (xy-1.-bord-(nxact-1.)/2.)*pitch;

  // fill cub (X coord  and Y coord):
  cub(,,1) = xy(,,1); cub(,,2) = xy(,,2);
  dis      = sqrt(cub(,,1)^2.+cub(,,2)^2.);
  if (dm(nm).pitchMargin == 0) {
    pitchMargin = 1.44;
  } else {
    pitchMargin = dm(nm).pitchMargin;
  }
  rad      = ((nxact-1.)/2.+pitchMargin)*pitch; //+1.44 is the margin
  inbigcirc= where(dis < rad);
  // 1 if valid actuator, 0 if not:
  // selection is done after interaction matrix is done
  cub(,,3) = 1; 

  // 1 if valid guard ring actuator, 0 if not:
  //cub(,,4) = (dis >= (pupr+extent*pitch)) & (dis < (pupr+(1.+extent)*pitch));
 // I don't use extrapolation actuator anymore.
  cub(,,4) = 0.;

  // converting to array coordinates:
  cub(,,1) = cub(,,1)+cent;
  cub(,,2) = cub(,,2)+cent;

  cub      = cub(*,);
  // cub now has two indices: first one is actuator number (valid or extrap)
  // second one is: 1:Xcoord, 2:Ycoord, 3:valid?, 4:extrapolation actuator?

  // filtering actuators outside of a disk radius = rad (see above)
  cub      = cub(inbigcirc,);

  cubval   = cub(where(cub(,3)),);
  
  nvalid   = int(sum(cubval(,3)));
  
  xy    = indices(size);
  x     = xy(,,2); y = xy(,,1);
  def	= array(float,dim,dim,nvalid);

  dm(nm)._x  = &(cubval(,1));
  dm(nm)._y  = &(cubval(,2));

  x = x(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);
  y = y(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);

  if (sim.verbose != 0) {write,format="\nCreating Influence function for actuator #%s","";}

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  for (i=1;i<=nvalid;i++) {
    if (sim.verbose != 0) {write,format="%d ",i;}
    tmpx       = clip(abs((x-cubval(i,1))/ir),1e-8,2.);
    tmpy       = clip(abs((y-cubval(i,2))/ir),1e-8,2.);
    tmp        = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*
                 (1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
    def(,,i)   = tmp*(tmpx <= 1.)*(tmpy <= 1.);
    gui_progressbar_text,swrite(format="Computing Influence Functions %d/%d",i,nvalid);
    gui_progressbar_frac,float(i)/nvalid;
    if ((disp == 1) && (sim.debug == 2)) {fma; pli,def(,,i);}
  }
  if (sim.verbose)  write,"";

  tmp=pitch/abs(ir);
  coupling = 1.- tmp^p1 + c*log(tmp)*tmp^p2;

  if (sim.debug>=1) write,format="coupling=%.2f%%  ",coupling*100;
  // look for extrapolation actuator stuff in v1.0.8 if needed

  fact = dm(nm).unitpervolt/max(def);
  
  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  if (sim.debug>=1) {
    piston=def(,,sum)*ipupil(dm(1)._n1:dm(1)._n2,dm(1)._n1:dm(1)._n2);
    tv,piston;
  }
  
  clean_progressbar;
  return def;
}

//----------------------------------------------------
func MakeEltPztIF(nm,&def,disp=)

  /* DOCUMENT function MakeEltPztIF(dm_structure,disp=)
     the influence functions are in microns per volt.
     same as MakePztIF but returns only local IF and
     start indices
   */
{
  coupling=dm(nm).coupling;

  // best parameters, as determined by a multi-dimensional fit
  // (see coupling3.i)
  a=[4.49469,7.25509,-32.1948,17.9493];
  p1 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;
  
  a = [2.49456,-0.65952,8.78886,-6.23701];
  p2 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [1.16136,2.97422,-13.2381,20.4395];
  irc = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;
  
  if (sim.debug==2) write,format="p1=%f  p2=%f  ir=%f\n",p1,p2,irc;

  dim   = dm(nm)._n2-dm(nm)._n1+1;
  size	= sim._size;
  nxact	= dm(nm).nxact;
  cent  = sim._cent;
  pitch = dm(nm).pitch;
  ir = irc*pitch;

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  // compute IF on partial (local) support:
  smallsize = long(ceil(2*ir+10));
  dm(nm)._eltdefsize = smallsize;
  xy    = indices(smallsize)-smallsize/2-0.5;
  x     = xy(,,1); y = xy(,,2);
  tmpx  = clip(abs(x/ir),1e-8,2.);
  tmpy  = clip(abs(y/ir),1e-8,2.);
  tmp   = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*(1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
  def   = tmp*(tmpx <= 1.)*(tmpy <= 1.);

  /*
  //smooth out the edges
  mask = def > 0;
  w = where(mask);
  k=makegaussian(5,1.);
  k=k/sum(k);
  cdef=def;
  //  for (i=1;i<=3;i++) {cdef=convVE(cdef,k);cdef(w)=def(w);}
  // 2004jun07: changed convVE -> convol for compat with yao_fftw
  for (i=1;i<=3;i++) {cdef=convol2d(cdef,k);cdef(w)=def(w);}
  def = cdef;
  */
  
  // compute location (x,y and i,j) of each actuator:
  cub   = array(float,nxact,nxact,2);
  // make X and Y indices array:
  xy    = indices(nxact);

  // express "centered" coordinate of actuator in pixels:
  xy    = (xy-1.-(nxact-1.)/2.)*pitch;

  // fill cub (X coord  and Y coord):
  cub(,,1) = xy(,,1); cub(,,2) = xy(,,2);
  // the following determine if an actuator is to be considered or not
  // relative to the pitchmargin parameter.
  dis      = sqrt(cub(,,1)^2.+cub(,,2)^2.);
  if (dm(nm).pitchMargin == 0) {
    pitchMargin = 1.44;
  } else {
    pitchMargin = dm(nm).pitchMargin;
  }
  rad      = ((nxact-1.)/2.+pitchMargin)*pitch;
  inbigcirc= where(dis < rad);
  // 1 if valid actuator, 0 if not:

  // converting to array coordinates:
  cub += cent;

  cub      = cub(*,);
  // cub now has two indices: first one is actuator number
  // second one is: 1:Xcoord, 2:Ycoord

  // filtering actuators outside of a disk radius = rad (see above)
  cubval   = cub(inbigcirc,);

  dm(nm)._nact = dimsof(cubval)(2);
  // following 4 lines changed on 2007apr19 to be consistent with order
  // with elt=0 (and thus consistent with how subapertures are numbered)
  dm(nm)._x  = &(cubval(,2));
  dm(nm)._y  = &(cubval(,1));
  dm(nm)._i1  = &(int(long(cubval(,2)-smallsize/2+0.5)-dm(nm)._n1));
  dm(nm)._j1  = &(int(long(cubval(,1)-smallsize/2+0.5)-dm(nm)._n1));

  def	= def(,,-)*array(1.f,dm(nm)._nact)(-,-,);

  // look for extrapolation actuator stuff in v1.0.8 if needed

  fact = dm(nm).unitpervolt/max(def);
  def = float(def*fact);
  dm(nm)._def = &def;

  return def;
}
//----------------------------------------------------
func MakeZernikeIF(nm,&def,disp=)

  /* DOCUMENT function MakeZernikeIF,dm_structure,ActIF,disp=
     modified 2004jan22 to have scaled as tip-tilt (e.g.
     1 arcsec/volt).
   */
{
  gui_progressbar_frac,0.;
  gui_progressbar_text,"Computing Influence Functions";
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  nzer	= dm(nm).nzer;
  cobs	= tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;
  
  def	= array(float,dim,dim,nzer);

  for (i=1;i<=nzer;i++) {
    def(,,i) = zernike_ext(i);
    if (disp == 1) {fma; pli,def(,,i);}
    gui_progressbar_frac,float(i)/nzer;
  }
  if (sim.verbose>=1) {write,format="Number of zernike :%d\n",nzer;}

  // normalization factor: one unit of tilt gives 1 arcsec:
  current = def(dim/2,dim/2,2)-def(dim/2-1,dim/2,2);
  fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;
  
  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  clean_progressbar;
  return def;
}

//----------------------------------------------------
func projectAnisoIF(nmaniso,nmlow,nmhigh,disp=)
/* DOCUMENT func projectAnisoIF(nmaniso,nmlow,nmhigh,disp=)
   This function finds the actuator commands to apply on dmlow and dmhigh
   to produce the anisoplanatism modes (which upper part is in dm(nmaniso)).
   nmaniso: # indice of anisoplanatism DM
   nmlow: # indice of low DM (at 0 altitude)
   nmhigh: # indice of high DM (at non zero altitude)

   computes alow and ahigh, which are #actuator x #anisomode.
   also compute "comaniso", which is the alow and ahigh put into a
   global total_#_actuator x 3 matrix, which can be directly multiplied/added
   to the global system command vector (see aoloop).
   Store them in extern variables for future use.
   SEE ALSO: MakeAnisoIF
 */
{
  extern alow,ahigh,comaniso;
  
  // we address here the zero altitude layer. The pupil is well defined.
  // cut a ipupil of the appropriate size:
  puplow = ipupil(dm(nmlow)._n1:dm(nmlow)._n2,dm(nmlow)._n1:dm(nmlow)._n2);
  w = where(puplow);
  // this transform def into a #spatial_point x nact array and retains only
  // the spatial point inside the pupil
  def = (*dm(nmlow)._def)(*,)(w,);
  // compute the IF covariance matrix
  defcov = def(+,)*def(+,);
  // find the inverse:
  defcovi = LUsolve(defcov);

  // now look at the anisoplanatism modes:
  // same, extract ipupil of appropriate dimension
  pupaniso = ipupil(dm(nmaniso)._n1:dm(nmaniso)._n2,dm(nmaniso)._n1:dm(nmaniso)._n2);
  w = where(pupaniso);
  // this transform def into a #spatial_point x nact array and retains only
  // the spatial points inside the pupil
  defa = -(*dm(nmaniso)._def)(*,)(w,);
  // compute the product act * mode:
  anisoproj = def(+,)*defa(+,);

  // command vector (matrices, 3 modes) to apply to DM to get a given mode 
  alow = defcovi(+,)*anisoproj(+,);

  // display:
  if (disp) {
    for (i=1;i<=3;i++) {tv,(*dm(nmlow)._def)(,,+)*alow(+,i)*puplow; hitReturn;}
  }

  // Now, the altitude DM:
  // it's basically the same thing, except now there is no well-defined pupil.
  // So we define here the pupil as the area which is controllable by the actuators.
  // that should be perfectly acceptable as the is the only area which will be seen
  // by any beam.
  puphigh = (*dm(nmhigh)._def)(,,sum);
  puphigh = (puphigh > 0.8*max(puphigh));
  w = where(puphigh);
  def = (*dm(nmhigh)._def)(*,)(w,);
  defcov = def(+,)*def(+,);
  defcovi = LUsolve(defcov);

  pupaniso = array(float,[2,sim._size,sim._size]);
  pupaniso(dm(nmhigh)._n1:dm(nmhigh)._n2,dm(nmhigh)._n1:dm(nmhigh)._n2) = puphigh;
  pupaniso = pupaniso(dm(nmaniso)._n1:dm(nmaniso)._n2,dm(nmaniso)._n1:dm(nmaniso)._n2);
  w = where(pupaniso);
  defa = (*dm(nmaniso)._def)(*,)(w,);
  anisoproj = def(+,)*defa(+,);

  ahigh = defcovi(+,)*anisoproj(+,);

  if (disp || (sim.debug == 2)) {
    for (i=1;i<=3;i++) {tv,(*dm(nmhigh)._def)(,,+)*ahigh(+,i)*puphigh; hitReturn;}
  }

  indexDm       = array(long,2,ndm);
  indexDm(,1)   = [1,dm(1)._nact];
  for (nm=2;nm<=ndm;nm++) {
    indexDm(,nm) = [indexDm(2,nm-1)+1,sum(dm(1:nm)._nact)];
  }
  comaniso = array(float,[2,sum(dm._nact),3]);
  comaniso(indexDm(1,nmlow):indexDm(2,nmlow),) = alow;
  comaniso(indexDm(1,nmhigh):indexDm(2,nmhigh),) = ahigh;

}
//----------------------------------------------------
func MakeAnisoIF(nm,&def,disp=)

  /* DOCUMENT function MakeAnisoIF,dm_structure,ActIF,disp=
     2004jan22: implemented normalization as for zernikeIF,
     i.e. based on the same amplitude tip that gives 1"
  */
{
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  cobs	= tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;
  
  def	= array(float,dim,dim,3);

  for (i=1;i<=3;i++) {
      def(,,i) = zernike_ext(i+3);
      if (disp == 1) {fma; pli,def(,,i);}
    }
  if (sim.verbose>=1) {write,format="Number of Anisoplanatism modes :%d\n",3;}

  // normalization factor: see MakeZernikeIF and MakeTipTiltIF
  tip = zernike_ext(2);
  current = tip(dim/2,dim/2,1)-tip(dim/2-1,dim/2,1);
  fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;

  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  return def;
}

//----------------------------------------------------
func MakeTipTiltIF(nm,&def,disp=)

  /* DOCUMENT function MakeTipTiltIF,dm_structure,ActIF,disp=
     adapted from makeZernikeIF
     modified 2004jan22 to make it normalized at 1"
   */
{
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  nzer	= 2;
  cobs	= tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;
  
  def	= array(float,dim,dim,nzer);

  for (i=1;i<=nzer;i++) {
      def(,,i) = zernike_ext(i+1);
      if (disp == 1) {fma; pli,def(,,i);}
    }

  // normalization factor: one unit of tilt gives 1 arcsec:
  current = def(dim/2,dim/2,1)-def(dim/2-1,dim/2,1);
  fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;

  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  return def;
}

//----------------------------------------------------

func MakeBimorphIF(nm,&def,disp=,cobs=)

  /* DOCUMENT:
     func MakeBimorphIF,dm_structure,&def,disp=
     This function build up the curvature mirror influence functions
     dim = output dimension of arrays
     pupd = pupil diameter in pixels
     SupportRadius = Radius at which the support points are located
     CompDim = Array dimension used for Computations of IFs (usually 4-8 x pupd)
     NOT NORMALIZED IN ANY WAY (arbitrary divided by 50 so that one gets
     acceptable phase for imat with a few tens volts.
  */
{
  extern actNumIm;
  local WhichRing,ActThetaIn,ActThetaOut,ActRadiusIn,ActRadiusOut;

  gui_progressbar_frac,0.;
  gui_progressbar_text,"Computing Influence Functions";
  
  dimdef   = dm(nm)._n2-dm(nm)._n1+1;

  dim	= sim._size;
  pupd	= sim.pupildiam;
  psize = tel.diam/sim.pupildiam;  // pixel in meter

  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*(dm(nm).alt)/psize;

  SupportRadius = 2.2;
  CompDim	= dim*2;

  NRing	= sum(*(dm(nm).nelperring) != 0); // Number of Rings
  NActPerRing = (*(dm(nm).nelperring))(1:NRing);
  NAct = sum(NActPerRing);
  // Compute the internal and external radius of each rings
  // given the number of actuators per rings:
  SurfOneAct = pi/sum(NActPerRing); // Surface of one actuator
  RInRing = array(float,NRing);     // Internal Radius
  ROutRing = array(float,NRing);    // External Radius
  if (is_set(cobs)) {RInRing(1) = cobs;} 
  // loop on ring number
  for (i=1;i<=NRing;i++) {
    ROutRing(i) = sqrt(NActPerRing(i)*SurfOneAct/pi+RInRing(i)^2.);
    if (i != NRing) RInRing(i+1) = ROutRing(i);
  }
  if (is_set(cobs)) {RInRing(1) = 0.;} 
  ROutRing(NRing) = 1.6;
  //  RInRing(NRing)  = 1.05;
  // now we got to determine the inner and outer radius and angle for
  // each actuators:
  WhichRing = array(1,NActPerRing(1));  // Ring index per actuator
  for (i=2;i<=NRing;i++) {grow,WhichRing,array(i,NActPerRing(i));}

  // offset angle of first electrode in rings:
  if (*dm(nm).angleoffset==[]) angleoffset=array(0.,NRing);
  else angleoffset=(*dm(nm).angleoffset)*pi/180.;

  // if rint and rout are specified, use it instead:
  if ((*dm(nm).rint)!=[]) RInRing=*dm(nm).rint;
  if ((*dm(nm).rout)!=[]) ROutRing=*dm(nm).rout;
  if (dm(nm).supportRadius) SupportRadius=dm(nm).supportRadius;
  
  // loop to determine radiuses and angle:
  for (i=1;i<=NRing;i++) {
    dtheta = 2*pi/NActPerRing(i)
    for (j=1;j<=NActPerRing(i);j++) {
      t1 = (j-1.)*dtheta + angleoffset(i);
      t2 = t1+dtheta;
      grow,ActThetaIn,t1 ;
      grow,ActThetaOut,t2 ;
      grow,ActRadiusIn,RInRing(i) ; 
      grow,ActRadiusOut,ROutRing(i) ;
    }
  }
  // Now build the actuator images:
  x = span(1,CompDim,CompDim)(,-:1:CompDim)-CompDim/2.-1;
  y = transpose(x);
  ang = atan(y,x);
  //  ang = atan(x,y);
  ang1 = ang + (ang < 0)*2*pi;
  ang2 = (ang1+pi)%(2*pi)+pi;
  //  rad = dist(CompDim)/(pupd/2.);
  rad = dist(CompDim)/(patchDiam/2.);
  d2 = clip(eclat(dist(CompDim)^2.),1e-5,);
  cpupil = rad < 1.;

  supportOffset=90.;
  if (dm(nm).supportOffset!=[]) supportOffset=dm(nm).supportOffset;
  supportOffset *= (pi/180.);

  tmp = abs(rad-SupportRadius)*5 + abs(ang1-0*pi/3-supportOffset) ;
  Support1 = where(tmp == min(tmp))(1);
  tmp = abs(rad-SupportRadius)*5 + abs(ang1-2*pi/3-supportOffset) ;
  Support2 = where(tmp == min(tmp))(1);
  tmp = abs(rad-SupportRadius)*5 + abs(ang1-4*pi/3-supportOffset) ;
  Support3 = where(tmp == min(tmp))(1);

  def = array(float,dimdef,dimdef,NAct);
  i1 = CompDim/2-dim/2+1;
  i2 = CompDim/2+dim/2;
  tmp = array(1.,CompDim,CompDim);
  tmp = tmp-0.5*cpupil
  tmp(Support1) = 0;
  tmp(Support2) = 0;
  tmp(Support3) = 0;
  if (disp == 1) {fma; pli, tmp; limits;}

  for (i=1;i<=NAct;i++) {
    // the following to avoid issues due to discontinuity of ang array at 0=2pi
    if (ActThetaOut(i)>(2*pi)) ang=ang2; else ang=ang1;
    Act = (rad >= ActRadiusIn(i)) * (rad < ActRadiusOut(i)) * \
      (ang >= ActThetaIn(i)) * (ang < ActThetaOut(i));
    if (i==1) {
      actNumIm = Act;
    } else {
      actNumIm += Act*i;
    }
    aif = fft(fft(eclat(long(Act)),1)/d2,-1);
    aif.re = eclat(aif.re);
    aif.im = eclat(aif.im);
    //    aif = aif - (aif(Support2)-aif(Support3))*y/(y(Support2)-y(Support3));
    //    aif = aif - (aif(Support1)-aif(Support2))*x/(x(Support1)-x(Support2));
    xdif = x(Support2)-x(Support3);
    if (xdif!=0.) aif = aif - (aif(Support2)-aif(Support3))*x/xdif;
    ydif = y(Support1)-y(Support2);
    if (ydif!=0.) aif = aif - (aif(Support1)-aif(Support2))*y/ydif;
    aif = aif - aif(Support3);
    aif = float(aif);
    tdef = aif(i1:i2,i1:i2);
    def(,,i)   = tdef(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);
    if (disp == 1) {
      fma;
      mypltitle,swrite(format="Influence Function %d/%d",i,NAct),[0.,-0.005],height=12;
      pli,def(,,i)*ipupil(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);
    }
    gui_progressbar_text,swrite(format="Computing Influence Functions %d/%d",i,NAct);
    gui_progressbar_frac,float(i)/NAct;
  }
  def = def/max(def)*pi; //just to keep things within reasonable values.
  def *= dm(nm).unitpervolt;  // adjustable normalization factor

  def = float(def/50.);  // factor 50 arbitrary to scale roughtly as PZT IF
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  clean_progressbar;
  return def;
}

//----------------------------------------------------

func MakeCurvWfsSubs(ns,dim,pupd,disp=,cobs=)

  /* DOCUMENT func MakeCurvWfsSubs(ns,dim,pupd,disp=)
   */
{
  local WhichRing,SubThetaIn,SubThetaOut,SubRadiusIn,SubRadiusOut;

  NRing	 = sum(*(wfs(ns).nsubperring) != 0);  // Number of Rings
  NSubPerRing = (*(wfs(ns).nsubperring))(1:NRing);
  NSub = sum(NSubPerRing);
  // Compute the internal and external radius of each rings
  // given the number of subapertures per rings:
  SurfOneSub = pi/sum(NSubPerRing); // Surface of one actuator
  RInRing = array(float,NRing);     // Internal Radius
  ROutRing = array(float,NRing);    // External Radius
  if (is_set(cobs)) {RInRing(1) = cobs;} 
  // loop on ring number
  for (i=1;i<=NRing;i++) {
    ROutRing(i) = sqrt(NSubPerRing(i)*SurfOneSub/pi+RInRing(i)^2.);
    if (i != NRing) RInRing(i+1) = ROutRing(i);
  }
  if (is_set(cobs)) {RInRing(1) = 0.;} 
  ROutRing(NRing) = 1.6;
  // now we got to determine the inner and outer radius and angle for
  // each actuators:
  WhichRing = array(1,NSubPerRing(1));  // Ring index per actuator
  for (i=2;i<=NRing;i++) {grow,WhichRing,array(i,NSubPerRing(i));}

  // offset angle of first subaperture in rings:
  if (*wfs(ns).angleoffset==[]) angleoffset=array(0.,NRing);
  else angleoffset=(*wfs(ns).angleoffset)*pi/180.;

  // if rint and rout are specified, use it instead:
  if ((*wfs(ns).rint)!=[]) RInRing=*wfs(ns).rint;
  if ((*wfs(ns).rout)!=[]) ROutRing=*wfs(ns).rout;

  // loop to determine radiuses and angle:
  for (i=1;i<=NRing;i++) {
    dtheta = 2*pi/NSubPerRing(i);
    for (j=1;j<=NSubPerRing(i);j++) {
      t1 = (j-1.)*dtheta + angleoffset(i);
      t2 = t1+dtheta;
      grow,SubThetaIn,t1 ;
      grow,SubThetaOut,t2 ;
      grow,SubRadiusIn,RInRing(i) ; 
      grow,SubRadiusOut,ROutRing(i) ;
    }
  }
  // Now build the subapertures images:
  x = span(1,dim,dim)(,-:1:dim)-dim/2.-1;
  y = transpose(x);
  ang = atan(y,x);
  ang1 = ang + (ang < 0)*2*pi;
  ang2 = (ang1+pi)%(2*pi)+pi;
  rad = dist(dim)/(pupd/2.);
  d2 = eclat(dist(dim)^2.);
  Subs = array(short,dim,dim,NSub);
  for (i=1;i<=NSub;i++) {
    // the following to avoid issues due to discontinuity of ang array at 0=2pi
    if (SubThetaOut(i)>(2*pi)) ang=ang2; else ang=ang1;
    tmp = (rad >= SubRadiusIn(i)) * (rad < SubRadiusOut(i)) * \
      (ang >= SubThetaIn(i)) * (ang < SubThetaOut(i));
    if (i>1) tmp = tmp-Subs(,,i-1); // avoid double points
    Subs(,,i) = (tmp==1);
    if (disp == 1) {fma; pli,tmp;}
  }
  // ok, in the following, I have changed the way one sums
  // the signal over the subaperture area in CurvWFS:
  // now I am passing a vector of indices for each subaperture
  // and using them for the sum.
  // sind is the ensemble of indice vectors
  // nsind is the number of actual indices for a given subaperture
  tmp 	= Subs(sum,sum,);
  sind 	= array(long,max(tmp),NSub)*0+1;
  nsind	= array(long,NSub);
  for (i=1;i<=NSub;i++) 
    {
      sind(1:tmp(i),i) = where(Subs(,,i) == 1); 
      nsind(i)	= sum(Subs(,,i));
    }

  wfs(ns)._sind = &(int(sind));
  wfs(ns)._nsind = &(int(nsind));

  return Subs;
}
//----------------------------------------------------
func _map2d(t,dim,cent)
{
  // assumes XY square array as input
  curdim = dimsof(t);
  if (curdim(2) != curdim(3)) {exit,"_map2d only takes square arrays";}
  cdim = curdim(2);


  if (cdim == dim) {return t;}

  if (cdim > dim) {
   if (is_void(cent)) {cent = cdim/2+1.;}
    _p1 = long(ceil(cent-dim/2.));
    _p2 = _p1+dim-1;
    return t(_p1:_p2,_p1:_p2,..);
  }

  if (cdim < dim) {
   if (is_void(cent)) {cent = dim/2+1.;}
    curdim(2:3) = [dim,dim];
    tmp = array(structof(t),curdim);
    _p1 = long(ceil(cent-cdim/2.));
    _p2 = _p1+cdim-1;
    if (_p2 != long(floor(cent+cdim/2.-1))) {exit,"Problem in _map2d";}
    tmp(_p1:_p2,_p1:_p2,..) = t;
    return tmp;
  }
}
//----------------------------------------------------

func noll(ord)
{
/* DOCUMENT func noll(ord)
   Noll variance of zernike for D/r0 = 1
     
   SEE ALSO:
 */
  innp = 0.;
  cov  = array(double,ord);
  for (j=2;j<=ord+1;j++) {
      tmp  = zernumero(j); n = tmp(1);  m = tmp(2);
      innp = gamma(14./3.)*gamma((2*n-14./3.+3.)/2.)/
        (2.^(14./3.)*gamma((14./3.+1.)/2.)*gamma((14./3.+1.)/2.)*
         gamma((2*n+14./3.+3.)/2.));
      cov(j-1) = (0.046/pi*(n+1.)*innp*(-1.)^((2*n-2*m)/2.));
  }
  return cov*(2.*pi)^(11./3.)/pi;
}
//----------------------------------------------------
func nollmat(ord)
{
/* DOCUMENT func nollmat(ord)
   Noll covariance matrix for D/r0 = 1
     
   SEE ALSO:
 */
  innp = 0.;
  cov  = array(double,[2,ord,ord]);
  for (j=2;j<=ord+1;j++) {
    for (jp=2;jp<=ord+1;jp++) {
      tmp  = zernumero(j);
      tmpp = zernumero(jp);
      n    = tmp(1);  m  = tmp(2);
      np   = tmpp(1); mp = tmpp(2);

      print,(n+np-14./3.+3.)/2.,(np-n+14./3.+1.)/2.,
        (n-np+14./3.+1.)/2.,(n+np+14./3.+3.)/2.;

      innp = gamma(14./3.)*gamma((n+np-14./3.+3.)/2.)/
        (2.^(14./3.)*gamma((np-n+14./3.+1.)/2.)*gamma((n-np+14./3.+1.)/2.)*
         gamma((n+np+14./3.+3.)/2.));
      cov(j-1,jp-1) = (0.046/pi*sqrt((n+1.)*(np+1.))*innp*float(m == mp)*
               (-1.)^(float(m == mp)*(n+np-2*m)/2.))*long((abs(j-jp)%2) == 0);
    }
  }
  cov   = (cov+transpose(cov))/2.;
  return cov*(2.*pi)^(11./3.)/pi;
}

//----------------------------------------------------
func rotby90(image,rot)

  /* DOCUMENT rotby90(image,rot)
     with:
     rot = 0 -> no rotation
     rot = 1 -> 90  deg anticlockwise
     rot = 2 -> 180 deg anticlockwise
     rot = 3 -> 270 deg anticlockwise
     AUTHOR: F.Rigaut, June 11, 2002.
     SEE ALSO: rotate
  */
{
  if (rot == 0) {return image;}
  if (rot == 1) {return transpose(image)(::-1,);}
  if (rot == 2) {return image(::-1,::-1);}
  if (rot == 3) {return transpose(image)(,::-1);}
}

//----------------------------------------------------

func MakePupil(dim,pupd,xc=,yc=,real=,cobs=)

  /* DOCUMENT func MakePupil(dim,pupd,xc=,yc=,real=)
   */
{
  if (real == 1) {
    pup = exp(-(dist(dim,xc=xc,yc=yc)/(pupd/2.))^60.)^0.69314;
  } else {
    //    xc;yc;info,xc;
    //    tv,dist(dim,xc=xc,yc=yc);pause,2000;
    pup = dist(dim,xc=xc,yc=yc) < (pupd+1.)/2.;
  }
  if (is_set(cobs)) {
    if (real == 1) {
      pup -= exp(-(dist(dim,xc=xc,yc=yc)/(pupd*cobs/2.))^60.)^0.69314;
    } else {
      pup -= dist(dim,xc=xc,yc=yc) < (pupd*cobs+1.)/2.;
    }
  }
    
  return pup;
}

//----------------------------------------------------

func fwhmStrehl(image,ps,lambda,teldiam,cobs,&strehl,&fwhm,&strehlab,&airy,&psf0,
                fibre=,source=,rmask=,dlambda=,psfcomp=,autoback=,silent=)
{
  /* DOCUMENT function fwhmStrehl
     Syntax: fwhmStrehl(image,ps,lambda,teldiam,cobs,strehl,fwhm,strehlab,&airy,
                fibre=,source=,rmask=,dlambda=,psf0=,psfcomp=,autoback=)
     Cette routine calcule le strehl et la fwhm de l'image. Le strehl
     est calcule de la facon suivante :
     1) on calcule la fto et on filtre les frequences > coupure du telescope
     2) On corrige de l'effet de filtrage par la taille finie du pixel
     3) Le cas echeant ("source" is set), on divise par la fto de la source
     4) On synthetise une tache d'airy pour la longueur d'onde en question
     5) On rebinne par 4, normalise, et on compare les maxima -> Strehl
     La FWHM est simplement calcule a partir du nombre de pixel de valeur > au
     max de l'image/2.

     ps		= pixel size en arcsec
     lambda	= en um
     teldiam	= telescope diameter en m
     cobs	= fraction obstruction central/telescope diameter
     fibre	= cf ci-dessus
     source	= source dimension in arcsec.
     autoback   = automatic normalization of background by interpolation of
                  zero point in MTF
     By convention, real plane coordinate (0,0) -> [dim/2,dim/2]
                    fourier plane coordinates (0,0) -> [0,0]
  */
  
  if (!is_set(dlambda)) {dlambda=0.;}

  if (!is_set(fibre)) {fibre = "disk";}
  // valeur autorisee pour fibre = 'disk' ou 'gaussian'


  dim	= dimsof(image)(2);

  // CALCUL DU STREHL :
  // Image theorique :
  tfto	= telfto(lambda,dlambda,teldiam,cobs,ps,dim,silent=1);
  tfto  = roll(tfto,[1,1]); // bug fixed 2007jul28 (no impact)
  mask	= (tfto > 1e-4);
  airy	= roll(abs(fft(tfto,1)));
  // was fftrebin, caused ripples in ima below, 2007jul28
  airy	= spline2(airy,4);

  // Image experimentale :
  fcoup	= dim/2/(lambda*1e-6/2./teldiam/4.848e-6/ps);
  ifto	= fft(roll(image),1);
  mask	= roll((dist(dim) < fcoup));
  // on soustrait le bruit moyen : what's that? 2007jul28
  /*
  if (sum(1.-mask) != 0) {
    mfl	  = sum(ifto.re*(1.-mask))/sum(1.-mask);
    mim	  = sum(ifto.im*(1.-mask))/sum(1.-mask);
    // print,max(ifto.re)/mfl,max(ifto.im)/mim;
    ifto.re  = (ifto.re-mfl)*mask;
    ifto.im  = (ifto.im-mim)*mask;
  }
  */
  // calcul de l'attenuation de la FTM par le moyennage du pixel :
  tmp		= roll(dist(dim))/(dim/2.)*pi/2.;
  ifto		= ifto/sinc(tmp);
  if (is_set(source)) {
    // Calcul de la FTM de la source (etoile artificielle ) :
    if (fibre == "gaussian") {
      // fibre = gaussian
      imfib   = exp(-(clip((dist(dim)/(source/ps/1.66))^2.,,20)));
      fibfto  = abs(fft(imfib,-1));
    } else {
      // fibre = camember
      tmp     = roll(dist(dim))/(dim/2.)*pi/2.*source/ps;
      fibfto  = sinc(tmp);
    }
  } else {
    fibfto  = ifto.re*0.+1.;
  }

  ima	= roll(float(fft((ifto/fibfto)*(tfto > 1e-8),-1)));
  psf0	= ima/sum(ima);
  // was fftrebin, caused ripples, 2007jul28
  ima	= spline2(ima,4);

  if (is_set(rmask)) // diaphragme les deux images
  {
    diap  = roll(dist(4*dim)) <= 4*rmask;
    wm    = (where2(ima == max(ima)))(,1);
    ima	  = ima*roll(diap,wm);
    airy  = airy*roll(diap);
  }

  airy	 = airy/sum(airy);
  ima	 = ima/sum(ima);
  strehl = max(ima)/max(airy);

  // Calcul de la largeur a mi-hauteur :
  fwhm	= sqrt(4./pi*sum(ima >= max(ima)/2.))*ps/4.;

  if (!is_set(silent)) write,format="fwhm = %f, Strehl = %f\n",fwhm,strehl;

  if (is_set(autoback))
  {
    tmp	  = ifto.re;
    t1	  = (tmp(2,1)+tmp(1,2)+tmp(1,0)+tmp(0,1))/4.;
    t12	  = (tmp(2,1)-tmp(3,1)+tmp(1,2)-tmp(1,3)+
             tmp(1,0)-tmp(1,-1)+tmp(0,1)-tmp(-1,1))/4.;
    ifto.re(1,1) = t1+t12;
    ifto  = ifto/(t1+t12);

    ima	  = roll(float(fft((ifto/fibfto)*(tfto > 1e-8),-1)));
    psf0  = ima/sum(ima);
    ima	  = fftrebin(ima,4);

    if (is_set(rmask))   // diaphragme les deux images
    {
      diap  = dist(4*dim) <= 4*rmask;
      ima   = ima*roll(diap,[ima(mxx,),ima(,mxx)]);
    }

    airy   = airy/sum(airy);
    ima	   = ima/sum(ima);
    strehlab = max(ima)/max(airy);

    // Calcul de la largeur a mi-hauteur :
    fwhm   = sqrt(4./pi*sum(ima >= max(ima)/2.))*ps/4.;

    write,format="%s\n","With automatic determination of background from points 0,1,2 of MTF :";
    write,format="fwhm = %f, Strehl = %f\n",fwhm,strehlab; 
  }
}
//------------------------------------------------
func telftoh1(f,u,v)
{
  e  = -1;
  if (abs(1-v) < 1e-12) {e = 1.;}

  tmp= (v^2/pi)*acos((f/v)*(1.+e*(1-u^2.)/(4.*f^2.)));

  return tmp;
}
//------------------------------------------------
func telftoh2(f,u)
{return -1.*(f/pi)*(1.+u)*sqrt((1.-(2*f/(1+u))^2.)* (1-((1-u)/(2*f))^2));}
//------------------------------------------------
func telftog(f,u)
{
  tmp = f*0.;
  z1 = where(f <= (1-u)/2.);
  z2 = where(f >= (1+u)/2.);
  z3 = where( f > (1-u)/2. & f < (1+u)/2. );
  if (exist(z3))
    {tmp(z3) = telftoh1(f(z3),u,1.) + telftoh1(f(z3),u,u) + telftoh2(f(z3),u);}
  if (exist(z1)) {tmp(z1) = u^2.;}
  if (exist(z2)) {tmp(z2) = 0.;}

  return tmp;
}
//------------------------------------------------
func telftot0(f,cobs)
{return (telftog(f,1.)+cobs^2*telftog(f/cobs,1.)-2.*telftog(f,cobs))/(1.-cobs^2.);}
//------------------------------------------------
func telfto(lambda,dlambda,teldiam,cobs,pixsize,dim,freqc=,npt=,silent=,returnpsf=)
{
  /* DOCUMENT func telfto
     Syntax: telfto(lambda,dlambda,teldiam,cobs,pixsize,dim,freqc=,
       npt=,silent=,returnpsf=)
     Computes and returns the Modulation transfer function (FTO is
     the french equivalent...) of a telescope with central obstruction
     and perfect optics.
     Parameters and keywords:
     lambda in microns
     dlambda in microns
     teldiam in meters
     cobs in fraction of teldiam
     pixsize in arcsec
     dim is output array dimension
     freqc is the cut-off frequency ?
     npt the number of point in dlambda (for calculations in this routine)
     silent has to be set to one for suppressing printout comments
     setting returnpsf has the effect of returning the PSF, not the MTF
  */
     
  if (is_void(dim))
    {exit,"function telfto,lambda,dlambda,teldiam,cobs,pixsize,dim,freqc=,npt=,silent=";}

  if (is_set(freqc)) {
    dlamb = 0.;
    lamb  = 1.;
    teld  = 1.;
    pixs  = freqc*lamb/teld/4.848/dim;}
  else {
    dlamb = dlambda;
    lamb  = lambda;
    teld  = teldiam;
    pixs  = pixsize;}

  if (!is_set(npt)) {npt = 5;}
  if (dlambda == 0.) {npt = 1;}

  if (!is_set(silent))
    {write,format="Cut-off frequency in pixels : %7.4f\n",teld*4.848/lamb*dim*pixs;}

  mtf = array(float,dim,dim);
  dd = clip(roll(dist(dim)),1e-10,);

  for (i=0;i<=npt-1;i++) {
    if (npt > 1) {l = lamb - dlamb*(i-(npt-1)/2.)/(npt-1);} else {l=lamb;}
    fc = teld*4.848/l;
    fc = fc*dim*pixs;
    f  = dd/fc;
    mask = (f <= 1.);
    f  = f * mask + f(dim/2+1,dim/2+1) * (1-mask);
    mtf = mtf + telftot0(f,cobs)*mask/npt;}


  mtf = mtf*sin(pi*dd/2./dim)/(pi*dd/2./dim);

  // big bug detected on 2007feb13: was returning fft(mtf)^2. for PSF
  // fortunately, I was not using the flag in any other routine, so no harm done
  if (is_set(returnpsf)) {return abs(fft(mtf,1));}
  return mtf;
}

//-------------------------------------------------------

func ftcb(te,tcal,tmir,gain,dim,x=)
{
/* DOCUMENT ftcb(te,tcal,tmir,gain,dim,x=)
   returns [f,hbo,hcor,hbf]
   AUTHOR: F.Rigaut, way back in 1996?
   SEE ALSO: 
 */
  f = indgen(dim)/te/2./dim;
  if (!is_void(x)) { f = x;}
  p = 2i*pi*f;

  //  hzoh	= (1.-exp(-te*p))/(te*p);
  hzoh	= 1.;
  hmir	= 1./(1.+tmir*p);
  hwfs	= (1.-exp(-te*p))/(te*p);
  hcal	= gain*exp(-tcal*p);

  hbo	= hzoh*hmir*hwfs*hcal/p/te;

  hcor	= float((1./(1.+hbo))*conj(1./(1.+hbo)));
  hbf	= float((hbo/(1.+hbo))*conj(hbo/(1.+hbo)));
  hbo	= float(hbo*conj(hbo));

  return ([f,hbo,hcor,hbf]);
}

//---------------------------------------------------------

func encircled_energy(image,&ee50,xc=,yc=)
  /* DOCUMENT encircled_energy(image,ee50,xc=,yc=)
 * Computes encircled energy function for a 2D array. 
 * Keywords xc and yc optionaly specify center about which the encircled
 * energy profile is to be computed. Returns optionally the value of the 
 * 50% encircled energy *diameter* (output).
 * F.Rigaut, Nov 2001.
 * SEE ALSO: findfwhm
 */
{
  im	= float(image);
  dim	= (dimsof(image))(2);

  if (xc == []) {xc = dim/2+1;}
  if (yc == []) {yc = dim/2+1;}

  e	= 1.9;
  npt	= 20;
  rv	= span(1.,(dim/2.)^(1./e),npt)^e;
  ee	= rv*0.;
  for (i=1;i<=npt;i++) 
    {
      fil   = MakePupil(dim,rv(i),xc=xc,yc=yc,real=1);
      rv(i) = sqrt(sum(fil)*4/pi); // that's a diameter
      ee(i) = sum(fil*im);
    }
  rv	= grow(0.,rv);
  ee	= grow(0.,ee);
  ee	= ee/sum(im);
  //fma;plg,ee,rv;
  xp	= span(0.,dim/2.,2*dim);
  yp	= interp(ee,rv,xp);
  ayp	= abs(yp-0.5);
  ee50   = xp(ayp(mnx));
  return yp;
}

//---------------------------------------------------------

func findfwhm(image,psize,nreb=)
/* DOCUMENT findfwhm(image,psize)
 * Determine the FWHM of an image. It simply rebins the image by
 * a factor of 4 and counts the number of pixels above the maximum
 * of the image divided by 2. Not the best accurate method but it's
 * robust
 * F.Rigaut, 2001/11/10.
 * SEE ALSO:
 */
{
  if (psize == []) psize=1;
  if (nreb==[]) nreb=4;

  if (nreb==1) eq_nocopy,imreb,image;
  else imreb = fftrebin(image,nreb);
  
  fwhm	= sum(imreb >= (max(imreb)/2.))/nreb^2.;
  fwhm	= sqrt(4.*fwhm/pi)*psize;

  return fwhm;
}
