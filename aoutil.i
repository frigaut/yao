/* AOUTIL.I
 *
 * A collection of utility routines to go with yao.i
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: aoutil.i,v 1.8 2010-04-15 02:36:53 frigaut Exp $
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
 * Revision 1.8  2010-04-15 02:36:53  frigaut
 *
 *
 * final commit to upgrade this repo to yao 4.5.1
 *
 * Revision 1.7  2008/05/12 18:00:53  frigaut
 * fixed problem with zernike diameter for altitude DMs.
 *
 * Revision 1.6  2008/05/11 14:03:56  frigaut
 * - implemented zernike wfs
 * - gotten rid of faulty round function in yao_util
 *
 * Revision 1.5  2007/12/20 13:34:52  frigaut
 * - various bug fixes
 * - better handlng of default parfile path
 * - better handling of options menu (WFS and DM)
 *
 * Revision 1.4  2007/12/19 19:44:19  frigaut
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
 * added things relative to inter actuator coupling (make_pzt_dm),
 * for regular and elt configs.
 * added tests relative to that in checkparameters
 *
 * Revision 1.7  2004/09/14 04:32:56  frigaut
 * several modifs to do with the creation of turbulent phase screens
 * - implemented cosf and sinf which take and return floats to save space
 * - started coding generate_phase_with_L0 in turbulence.i, not finished
 * - modif YORICKEXE in makefiles, just "yorick" did not do it on osx
 * - modifs ok for both veclib and fftw implementations
 *
 * Revision 1.6  2004/08/02 07:10:53  frigaut
 * Added routine get_turb_phase_initCheckOverflow, which checks for Y
 * interpolation indices larger than the max Y dimension, which would
 * cause problem in get_turb_phase.
 *
 * Revision 1.5  2004/07/29 04:06:50  frigaut
 * added cvs dollar Log in header
 *
 *
 * func graphic_config(subsystemnum,dmnum)
 * func check_parameters(void)
 * func disp2d(ar,xpos,ypos,area,zoom=,power=,init=,nolimits=)
 * func hysteresis(v,n,first=) 
 * func modal_gain_optimization(disp=,update=)
 * func ft_cb_ao_simul(FrameDelay,gain,dim)
 * func build_dm_modes(disp=)
 * func wfs_check_pixel_size(ns,&binindices,&centroidw,printheader=,silent=)
 * func make_pzt_dm(nm,&def,disp=)
 * func make_pzt_dm_elt(nm,&def,disp=)
 * func make_zernike_dm(nm,&def,disp=)
 * func project_aniso_dm(nmaniso,nmlow,nmhigh,disp=)
 * func make_aniso_dm(nm,&def,disp=)
 * func make_tiptilt_dm(nm,&def,disp=)
 * func make_curvature_dm(nm,&def,disp=,cobs=)
 * func make_curv_wfs_subs(ns,dim,pupd,disp=,cobs=)
 * func _map2d(t,dim,cent)
 * func noll(ord)
 * func nollmat(ord)
 * func rotby90(image,rot)
 * func make_pupil(dim,pupd,xc=,yc=,real=,cobs=)
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
func get_turb_phase_initCheckOverflow(void)
/* DOCUMENT func get_turb_phase_initCheckOverflow
   This routine has the sole purpose of checking the possible "overflow"
   of the Y index in the future calls to get_turb_phase. In any of the
   Y index at which the interpolation is to be done is larger than the
   phase screen(s) Y dimension, then an error is flagged.
   SEE ALSO: get_turb_phase_init, get_turb_phase
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
    write,"   - use larger phase screens (Y axis) using 'create_phase_screens'";
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
    write,"   - use larger phase screens (Y axis) using 'create_phase_screens'";
    write,"   - Modify the Y position of your Perf. Star";
    write,"   - Lower the altitude of the offending atmospheric layer";
    exit;
  }
}
//----------------------------------------------------
func graphic_config(subsystemnum,dmnum)
/* DOCUMENT func configGraphic(void)
   Plots a graphical representation of the system config,
   per subsystem and per level (altitude)
   subsystemnum and dmnum optional. If not set all subsystems
   and dms are displayed
   SEE ALSO:
 */
{
  t       = span(0.,2*pi,100);
  col     = ["red","blue","green","cyan","magenta","yellow"];
  markers = ['1','2','3','4','5','6','7','8','9'];
  markers = char(indgen(36)+48); // 1,2,...A,B...
  maxcol  = numberof(col);

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

      for (j=1;j<=nwfs;j++) {
        if ( wfs(j).subsystem != nss ) continue;
        // overplot subaperture position and size for first wfs in this subsystem
        first_wfs_plotted = 0;
        if ( (first_wfs_plotted==0) && (dm(i).alt==0.) && (dm(i).type!="tiptilt") ) {
          for (jj=1;jj<=wfs(j)._nsub4disp;jj++) {
            if ((*wfs(j)._validsubs)(jj)==0) continue;
            x1 = (*wfs(j)._istart)(jj);
            y1 = (*wfs(j)._jstart)(jj);
            subsize    = sim.pupildiam/wfs(j).shnxsub;
            if (wfs(j).npixpersub) subsize = wfs(j).npixpersub;
            l1 = subsize;
            // convert pixels in meters
            x1 = (x1-sim._cent)*tel.diam/sim.pupildiam;
            y1 = (y1-sim._cent)*tel.diam/sim.pupildiam;
            l1 = l1*tel.diam/sim.pupildiam;
            plg,_(y1,y1,y1+l1,y1+l1,y1),_(x1,x1+l1,x1+l1,x1,x1),color=[50,50,50];
          }
          first_wfs_plotted = 1;
        }
      }
      
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
        plg,loc(,2),loc(,1),type=0,marker=markers(i),color=col(i%maxcol);

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
        plg,offsets(2)+rad*sin(t),offsets(1)+rad*cos(t),color=col(i%maxcol),width=3;

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
        plg,loc(,2),loc(,1),type=0,marker=markers(i),color=col(i%maxcol);
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
        plg,offsets(2)+rad*sin(t),offsets(1)+rad*cos(t),color=col(i%maxcol),width=3;
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
func check_parameters(void)
/* DOCUMENT func check_parameters(void)
   Check the parameters in yao parfile
   - set defaults
   - value valid?
   - compatibility with other parameters
   SEE ALSO:
 */
{
  extern sim,atm,wfs,dm,mat,tel,target,gs,loop;
  write,format="Checking parameters ... %s\n","";
  
  //==================================================================  
  // BASIC EXISTENCE AND CONSISTENCY CHECKINGS AND DEFAULT ASSIGNMENTS
  //==================================================================  

  if (nwfs == []) {exit,"nwfs has not been set";}
  if (ndm == []) {exit,"ndm has not been set";}
  
  // SIM STRUCTURE
  if (sim.pupildiam == 0) {exit,"sim.pupildiam has not been set";}

  // ATM STRUCTURE
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

  // WFS STRUCTURE
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
    if ((wfs(ns).type == "hartmann") && (wfs(ns).shmethod == 2)) {
      if ((strlen(wfs(ns).fsname)>0) && (strlen(wfs(ns).fstop)>0)) {
        write,format="You have set both wfs(%d).fsname and wfs(%d).fstop\n",ns,ns;
        write,format="I will ignore wfs(%d).fstop and use wfs(%d).fsname !\n",ns,ns;
      }
      if ((strlen(wfs(ns).fsname)==0) && (strlen(wfs(ns).fstop)==0)) {
        write,format="No field stop defined for wfs %d. Setting to 'square'\n",ns;
        wfs(ns).fstop = "square";
      }
      if (wfs(ns).fssize==0) {
        write,format="wfs(%d).fssize has not been set, will be forced to subap FoV\n",ns;
      }
    }    
    if (wfs(ns).nintegcycles == 0) {wfs(ns).nintegcycles = 1;}
    if (wfs(ns).fracIllum == 0) {wfs(ns).fracIllum = 0.5;}
    if (wfs(ns).optthroughput == 0) {wfs(ns).optthroughput = 1.0;}
    wfs.ron = float(wfs.ron);
  }

  // DM STRUCTURE
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
      if ((dm(nm).coupling<0.04) || (dm(nm).coupling>0.8)) {
        write,format="Invalid value for dm(%d).coupling -> %f\n",nm,dm(nm).coupling;
        exit,"Valid values from 0.04 to 0.80";
      }
    }
    if ( (dm(nm).type == "zernike") && (dm(nm).nzer == 0) ) {
      exit,swrite(format="dm(%d).nzer has not been set",nm);
    }    
    if ( (dm(nm).type == "kl") && (dm(nm).nkl == 0) ) {
      exit,swrite(format="dm(%d).nkl has not been set",nm);
    }    
    if (dm(nm).irexp==1) {
      if (dm(nm).irfact == 0.) {
        dm(nm).irfact = 1.;
        write,format="dm(%d).irfact set to %f\n",nm,dm(nm).irfact;
      }
    }
  }

  // MAT STRUCTURE
  if ((*mat.condition) == []) {exit,"mat.condition has not been set";}
  if (numberof(*mat.condition) != max(_(wfs.subsystem,dm.subsystem)) ) {
    exit,"dimension of *mat.condition is not equal to the number of subsystems";
  }
  if (mat.file == string()) {mat.file = "";}

  // TEL STRUCTURE
  if (tel.diam == 0) exit,"tel.diam has not been set";

  // TARGET STRUCTURE
  if ((*target.lambda) == []) exit,"target.lambda has not been set";
  if ((*target.xposition) == []) exit,"target.xposition has not been set";
  if ((*target.yposition) == []) exit,"target.yposition has not been set";
  if ((*target.dispzoom) == []) {
    target.dispzoom = &(array(1.,numberof(*target.lambda)));
  }
  if (nallof(_(numberof(*target.xposition), \
               numberof(*target.yposition)) == numberof(*target.dispzoom) )) {
    exit,"Some elements within target.xposition, yposition, dispzoom "+\
      "do not have the same number of elements.";
  }

  // GS STRUCTURE
  if (anyof(wfs.gsalt != 0) && (gs.lgsreturnperwatt == 0)) {
    gs.lgsreturnperwatt = 22.;
    write,format="gs.lgsreturnperwatt set to %f\n",gs.lgsreturnperwatt;
  }
  if (anyof(wfs.gsalt == 0) && (gs.zeropoint == 0)) {
    exit,"You have some NGS and gs.zeropoint has not been set";
  }
  
  // LOOP STRUCTURE
  if (loop.gain == 0) write,format="%s\n","Warning: loop.gain = 0";
  if ( (numberof(*loop.gainho)) != (numberof(*loop.leakho)) ) \
    exit,"*loop.gainho should have same number of element as *loop.leakho";
  if (loop.niter == 0) exit,"loop.niter has not been set";
  if (loop.ittime == 0) exit,"loop.ittime has not been set";
  if (loop.startskip == 0) loop.startskip = 10;
  if (loop.skipevery == 0) loop.skipevery = loop.niter;
  if (loop.skipby == 0) loop.skipby = 10000;
  if (loop.modalgainfile == string()) loop.modalgainfile = "";
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
    // Are we using a WFS we know?
    wfs_type = strtolower(wfs(ns).type);
    if ( (wfs_type != "curvature") && (wfs_type != "hartmann") &&
         (wfs_type != "zernike")   && (wfs_type !="kl") ) {
      // check if this is a user supplied function
      cmd = swrite(format="totype = typeof(%s)",wfs(ns).type);
      include,[cmd],1;
      if ( totype != "function") {
        error,swrite(format="wfs(%d).type : Unknown value or non-existing function \"%s\"",
                ns,wfs(ns).type)
      }
    } else wfs(ns).type = wfs_type;
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
    
    // Are we using a DM we know?
    dm_type = strtolower(dm(nm).type);
    if ( (dm_type != "bimorph") && (dm_type != "stackarray") &&
         (dm_type != "zernike") && (dm_type != "kl") &&
         (dm_type != "segmented") &&
         (dm_type != "tiptilt") && (dm_type != "aniso") ) {
      // check if this is a user supplied function
      cmd = swrite(format="totype = typeof(%s)",dm(nm).type);
      include,[cmd],1;
      if ( totype != "function") {
        error,swrite(format="dm(%d).type : Unknown value or non-existing function \"%s\"",
        nm,dm(nm).type)
      }
    } else dm(nm).type = dm_type;
    
    dm(nm)._eiffile = parprefix+swrite(format="-if%d",nm)+"-ext.fits";
  }


// now possible (2009oct6, v4.3.0)
//   if (anyof(dm.elt == 1) && anyof(dm.type == "aniso")) {
//     exit,"You can not use currently dm.elt=1 with anisoplanatism modes";
//   }
  
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


  if (opt!=[]) {
    if (opt.misreg==[]) {
      noptics = numberof(opt.phasemaps);
      for (i=1;i<=noptics;i++) opt.misreg = [0.,0.];
    }
    opt.misreg= float(opt.misreg);
  }
  
  write,format="%s\n","OK";
}
//----------------------------------------------------
func disp2d(ar,xpos,ypos,area,zoom=,power=,init=,nolimits=)
/* DOCUMENT func disp2d(arrayptr,xpos,ypos,area,zoom=,power=,init=)
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
    if (basezoomptr == []) {
      basezoomptr=array(pointer,10);
    }
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
        off = basezoom(i)/(dimsof(*ar(1))(2));
        pli,*ar(i),xpos(i)-basezoom(i)-off,ypos(i)-basezoom(i)-off,
          xpos(i)+basezoom(i)-off,ypos(i)+basezoom(i)-off;
      }
    } else {
      for (i=1;i<=nim;i++) {
        off = basezoom(i)/(dimsof(*ar(1))(2));
        pli,*ar(i)^power,xpos(i)-basezoom(i)-off,ypos(i)-basezoom(i)-off,
          xpos(i)+basezoom(i)-off,ypos(i)+basezoom(i)-off;
      }
    }
  } else {
    if (!is_set(power)) {
      for (i=1;i<=nim;i++) {
        off = basezoom(i)/(dimsof(ar(,,1))(2));
        pli,ar(,,i),xpos(i)-basezoom(i)-off,ypos(i)-basezoom(i)-off,
          xpos(i)+basezoom(i)-off,ypos(i)+basezoom(i)-off;
      }
    } else {
      for (i=1;i<=nim;i++) {
        off = basezoom(i)/(dimsof(ar(,,1))(2));
        pli,ar(,,i)^power,xpos(i)-basezoom(i)-off,ypos(i)-basezoom(i)-off,
          xpos(i)+basezoom(i)-off,ypos(i)+basezoom(i)-off;
      }
    }
  }
  
  if ((is_set(init)) & (!is_set(nolimits))) {
    limits,square=1;
    limits;
  }
}
//--------------------------------------------------------------------------

func modal_gain_optimization(disp=,update=)
/* DOCUMENT func modal_gain_optimization(disp=,update=)

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
  cberr        = fitsRead("cberr.fits"); // CB of actuator error
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
    {errTF(,n) = ft_cb_ao_simul(ao.LoopFrameDelay,toGains(n),length/2)(,1);}
  
  
  // Loop over controlled modes:
  for (i=1;i<=NModesControlled;i++)
  {
    // Build PSD of the error
    psderr = psd(cbmoderr(i,),length,samp=ao.LoopItTime,filter=6,silent=1)(2:,2);
    pause,10;
    // ^^ psderr is in fact length/2 long

    // Build current gain transfer function
    curErrTF = ft_cb_ao_simul(ao.LoopFrameDelay,modgains(i),length/2)(,1);

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

func ft_cb_ao_simul(FrameDelay,gain,dim)
/* DOCUMENT func ft_cb_ao_simul(FrameDelay,gain,dim)

   NOT UPGRADED TO VERSION 2.
   DO NOT USE UNTIL UPGRADED.

   This routine simulates the time aspect of the numerical loop
   and compute the associated transfer functions (error, closeloop)
   for further use in modal_gain_optimization().
   Inputs:
   FrameDelay: frame delay in close loop.
   gain: AO loop gain
   dim: desired linear size of the output transfer functions
   SEE ALSO: modal_gain_optimization
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
func ss_noise(nss)
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
  modes = build_dm_modes(nmv,actmodes,modvar);
  tmp = actmodes(+,)*cmat(+,);
  modcov = trace(tmp(,+)*tmp(,+));
  // now the measurements are in arcsec, so we got to convert to rd of
  // phase difference at the edge of the subaperture, which is what
  // the SH measure:
  nsv = where(wfs.subsystem == nss);
  if (numberof(nsv) != 1) error,"zero or more than one wfs in subsystem";
  ns = nsv(1);

  // subsize in pixel:
  subsize  = sim.pupildiam/wfs(ns).shnxsub(0);
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;

  // subsize in meters:
  subsize_m = subsize * tel.diam/sim.pupildiam;
  //subsize = tel.diam/wfs(ns).shnxsub;
  arcsectord = 4.848e-6*subsize_m*2*pi/wfs(ns).lambda/1e-6;
  noise = modcov*modvar/arcsectord^2.;
  write,format="Total noise on phase for 1rd2 per subaperture = %f rd2\n",sum(noise);
  fma;plh,noise; limits,square=0; limits;
  return noise;
}
//----------------------------------------------------

func build_dm_modes(nm,&u,&modvar,&eigenvalues,disp=)

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
  tpup  = sum(defpup);
  /*
    p = (def*ipupil)(sum,sum,)/tpup;
    def = def-p(-,-,);
    def = reform(def,long(ao._size)*ao._size,ao._DmNAct);
    def = def(where(ipupil),);
  */
  dd  = tmp(+,)*tmp(+,);
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
  norm  = sqrt((modes^2.*ipupil)(avg,avg,));
  modes = modes/norm(-,-,);

  ActIF = modes(,,1:-1) //*(1./indgen(ao._DmNAct-1))(-,-,);
  ao._DmNAct = ao._DmNAct-1;
  do_imat,disp=1;
  build_cmat,all=1,nomodalgain=1;
  mn  = cMat(,+)*cMat(,+);
  mn  = diag(mn)*indgen(ao._DmNAct-1)^2.
  if (disp == 1) {fma;plg,mn;}
  return modes;
}

//----------------------------------------------------

func make_curv_wfs_subs(ns,dim,pupd,disp=,cobs=)

  /* DOCUMENT func make_curv_wfs_subs(ns,dim,pupd,disp=)
   */
{
  local WhichRing,SubThetaIn,SubThetaOut,SubRadiusIn,SubRadiusOut;

  NRing  = sum(*(wfs(ns).nsubperring) != 0);  // Number of Rings
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
  tmp   = Subs(sum,sum,);
  sind  = array(long,max(tmp),NSub)*0+1;
  nsind = array(long,NSub);
  for (i=1;i<=NSub;i++) 
    {
      sind(1:tmp(i),i) = where(Subs(,,i) == 1); 
      nsind(i)  = sum(Subs(,,i));
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

func make_pupil(dim,pupd,xc=,yc=,real=,cobs=)

  /* DOCUMENT func make_pupil(dim,pupd,xc=,yc=,real=)
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

     ps   = pixel size en arcsec
     lambda = en um
     teldiam  = telescope diameter en m
     cobs = fraction obstruction central/telescope diameter
     fibre  = cf ci-dessus
     source = source dimension in arcsec.
     autoback   = automatic normalization of background by interpolation of
                  zero point in MTF
     By convention, real plane coordinate (0,0) -> [dim/2,dim/2]
                    fourier plane coordinates (0,0) -> [0,0]
  */
  
  if (!is_set(dlambda)) {dlambda=0.;}

  if (!is_set(fibre)) {fibre = "disk";}
  // valeur autorisee pour fibre = 'disk' ou 'gaussian'


  dim = dimsof(image)(2);

  // CALCUL DU STREHL :
  // Image theorique :
  tfto  = telfto(lambda,dlambda,teldiam,cobs,ps,dim,silent=1);
  tfto  = roll(tfto,[1,1]); // bug fixed 2007jul28 (no impact)
  mask  = (tfto > 1e-4);
  airy  = roll(abs(fft(tfto,1)));
  // was fftrebin, caused ripples in ima below, 2007jul28
  airy  = spline2(airy,4);

  // Image experimentale :
  fcoup = dim/2/(lambda*1e-6/2./teldiam/4.848e-6/ps);
  ifto  = fft(roll(image),1);
  mask  = roll((dist(dim) < fcoup));
  // on soustrait le bruit moyen : what's that? 2007jul28
  /*
  if (sum(1.-mask) != 0) {
    mfl   = sum(ifto.re*(1.-mask))/sum(1.-mask);
    mim   = sum(ifto.im*(1.-mask))/sum(1.-mask);
    // print,max(ifto.re)/mfl,max(ifto.im)/mim;
    ifto.re  = (ifto.re-mfl)*mask;
    ifto.im  = (ifto.im-mim)*mask;
  }
  */
  // calcul de l'attenuation de la FTM par le moyennage du pixel :
  tmp   = roll(dist(dim))/(dim/2.)*pi/2.;
  ifto    = ifto/sinc(tmp);
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

  ima = roll(float(fft((ifto/fibfto)*(tfto > 1e-8),-1)));
  psf0  = ima/sum(ima);
  // was fftrebin, caused ripples, 2007jul28
  ima = spline2(ima,4);

  if (is_set(rmask)) // diaphragme les deux images
  {
    diap  = roll(dist(4*dim)) <= 4*rmask;
    wm    = (where2(ima == max(ima)))(,1);
    ima   = ima*roll(diap,wm);
    airy  = airy*roll(diap);
  }

  airy   = airy/sum(airy);
  ima  = ima/sum(ima);
  strehl = max(ima)/max(airy);

  // Calcul de la largeur a mi-hauteur :
  fwhm  = sqrt(4./pi*sum(ima >= max(ima)/2.))*ps/4.;

  if (!is_set(silent)) write,format="fwhm = %f, Strehl = %f\n",fwhm,strehl;

  if (is_set(autoback))
  {
    tmp   = ifto.re;
    t1    = (tmp(2,1)+tmp(1,2)+tmp(1,0)+tmp(0,1))/4.;
    t12   = (tmp(2,1)-tmp(3,1)+tmp(1,2)-tmp(1,3)+
             tmp(1,0)-tmp(1,-1)+tmp(0,1)-tmp(-1,1))/4.;
    ifto.re(1,1) = t1+t12;
    ifto  = ifto/(t1+t12);

    ima   = roll(float(fft((ifto/fibfto)*(tfto > 1e-8),-1)));
    psf0  = ima/sum(ima);
    ima   = fftrebin(ima,4);

    if (is_set(rmask))   // diaphragme les deux images
    {
      diap  = dist(4*dim) <= 4*rmask;
      ima   = ima*roll(diap,[ima(mxx,),ima(,mxx)]);
    }

    airy   = airy/sum(airy);
    ima    = ima/sum(ima);
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

  hzoh  = (1.-exp(-te*p))/(te*p);
  //hzoh  = 1.;
  hmir  = 1./(1.+tmir*p);
  hwfs  = (1.-exp(-te*p))/(te*p);
  hcal  = gain*exp(-tcal*p);

  hbo = hzoh*hmir*hwfs*hcal/(1-exp(-p*te));

  hcor  = float((1./(1.+hbo))*conj(1./(1.+hbo)));
  hbf = float((hbo/(1.+hbo))*conj(hbo/(1.+hbo)));
  hbo = float(hbo*conj(hbo));

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
  im  = float(image);
  dim = (dimsof(image))(2);

  if (xc == []) {xc = dim/2+1;}
  if (yc == []) {yc = dim/2+1;}

  e = 1.9;
  npt = 20;
  rv  = span(1.,(dim/2.)^(1./e),npt)^e;
  ee  = rv*0.;
  for (i=1;i<=npt;i++) 
    {
      fil   = make_pupil(dim,rv(i),xc=xc,yc=yc,real=1);
      rv(i) = sqrt(sum(fil)*4/pi); // that's a diameter
      ee(i) = sum(fil*im);
    }
  rv  = grow(0.,rv);
  ee  = grow(0.,ee);
  ee  = ee/sum(im);
  //fma;plg,ee,rv;
  xp  = span(0.,dim/2.,2*dim);
  yp  = interp(ee,rv,xp);
  ayp = abs(yp-0.5);
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
  
  fwhm  = sum(imreb >= (max(imreb)/2.))/nreb^2.;
  fwhm  = sqrt(4.*fwhm/pi)*psize;

  return fwhm;
}


func phi2zer(i_num, phase,pup, nzer=, kl=)
/* DOCUMENT phi2zer(phase,nzer=,&zcoeff)
 * Do the decomposition of a phase screen on N zernike (or KL)
 * coefficients. To be used in yao_mcao.i, to save a circular buffer
 * with the Z (or KL) coefficients of the residual phase.
 * B.Neichel, 2009/05/13.
 * SEE ALSO:
 */
{
  //require,mcao_i_dir+"lib/libkl.i";
  //no need anymore, they are include in myst_init.i

  extern conv, w_def;
  
  //first, find the dimension of the phase
  dim_phi = dimsof(phase)(2); //in x and y

  if (dimsof(phase)(1) == 3) {ndir = dimsof(phase)(4);};
  //we are dealing with multiple directions
 

  //This should be computed only once (for the first iteration)
  //then, we only need "conv"
  if (i_num == 1) {
    if (kl) {

      kl_tab = array(float,[3,dim_phi,dim_phi,nzer]);
      kl_tab = make_kl(nzer,dim_phi,var,outbas,pup);
      w_def = where(kl_tab(,,1));
      def_kl =kl_tab(*,)(w_def,);
      conv = generalized_inverse(def_kl);
    
    } else {
      prepzernike,dim_phi,sim.pupildiam;
      tmp = zernike(1);
      w_def = where(tmp);
      zer_tab = array(float,[3,dimsof(tmp)(2),dimsof(tmp)(2),nzer]);
    
      for (z=1;z<=nzer;z++) {
        zer_tab(,,z) = zernike(z); };
      
      def_zer =zer_tab(*,)(w_def,);///norm;
      conv = generalized_inverse(def_zer);

    };
  };

  //do the loop on each directions here.
  if (ndir) { ztmp = array(float,[2,ndir,nzer]);}
  else { ztmp = array(float,[2,1,nzer]); ndir = 1;};

  for (nnn=1;nnn<=ndir;nnn++){
    ztmp(nnn,) = conv(,+)*(phase(,,nnn)(w_def))(+);}
  
  return ztmp;

}


func make_fieldstop(ns)
{
  extern wfs;

  // build field stop for WFS ns
  
  subsize    = sim.pupildiam/wfs(ns).shnxsub(0);
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;
  sdim       = long(2^ceil(log(subsize)/log(2)+1));  
  fftpixsize = wfs(ns).pixsize/wfs(ns)._rebinfactor; // fft pixel size in arcsec
  
  // if field stop size has not been set, set it to the subap fov (which seems
  // the best thing to do)
  if (wfs(ns).fssize==0) wfs(ns).fssize = wfs(ns).npixels * wfs(ns).pixsize;
  if (sim.verbose>=1) {
    write,format="WFS#%d Field stop size = %f\n",ns,wfs(ns).fssize;
  }
  
  fs_size_fftpix = wfs(ns).fssize/fftpixsize;
  
  if (wfs(ns).fstop=="none") {
    fs = array(1n,[2,sdim,sdim]);
  } else if (wfs(ns).fstop=="round") {
    fs = dist(sdim,xc=sdim/2+0.5,yc=sdim/2+0.5)<=(fs_size_fftpix/2.);    
  } else { // anything else -> square FS
    fs = array(0n,[2,sdim,sdim]);
    hsize = long(round(fs_size_fftpix/2.));
    hsize = clip(hsize,,sdim/2);
    fs(sdim/2-hsize+1:sdim/2+hsize,sdim/2-hsize+1:sdim/2+hsize) = 1n;
  }
  if (wfs(ns).fsoffset!=[]) {
    // only roll by an integer # of fft pixel:
    fsoffset = long(round(wfs(ns).fsoffset/fftpixsize));
    fs = roll(fs,fsoffset);
    // refresh value of fsoffset to actual value:
    wfs(ns).fsoffset = fsoffset * fftpixsize;
    if (sim.verbose>=2) {
      write,format="WFS#%d Field Stop offsets actual values: (%f,%f) arcsec\n",
        ns,(wfs(ns).fsoffset)(1),(wfs(ns).fsoffset)(2);
    }
  }
  wfs(ns)._submask = &(float(roll(fs)));
  wfs(ns)._domask = 1l;
}



func generate_vib_time_serie(sampling_time,length,white_rms,one_over_f_rms,peak,peak_rms,peak_width=)
/* DOCUMENT generate_vib_time_serie
   sampling_time  = loop sampling time [seconds]
   length         = number of point desired in generated time serie
   white_rms      = rms of white noise [arcsec]
   one_over_f_rms = rms of 1/f noise (from 1 Hz included up to cutoff 
                    frequency)
   peak           = vector containing the frequency [Hz] at which vibration 
                    peaks are to be generated
   peak_rms       = vector (same number of elements as peak) containing the 
                    rms of each peaks (arcsec)
   peak_width     = keyword optionaly containing a vector (same number of
                    elements as peak) of width of each peaks [FWHM in Hz,
                    default one Hz]. 
   SEE ALSO:
 */
{
  local psd;
  
  // construct power spectrum:
  npt = long(ceil(length/2.))+1; //+1 to include both 0 and freqmax
  freq = span(0.,1./sampling_time/2.,npt);
  psdall = array(0.,npt);
  
  // white noise:
  psd = array(1.,npt); 
  // normalize (power theorem)
  psd = psd/sum(psd)*white_rms^2.;
  psdall = psd;
  
  // one over f noise:
  psd = 1./clip(freq,freq(2),);
  onehz = long(ceil(1./freq(2)));
  psd = psd/sum(psd(onehz:))*one_over_f_rms^2.;
  psdall += psd;
  
  // peaks:
  if (peak!=[]) {
    if (peak_width==[]) peak_width=array(freq(2)/10.,numberof(peak));
    if ( (numberof(peak_rms) != numberof(peak)) ) \
      error,"numberof(peak_rms) != numberof(peak)";
    if ( (numberof(peak_width) != numberof(peak)) ) \
      error,"numberof(peak_width) != numberof(peak)";

    for (i=1;i<=numberof(peak);i++) {
      if ( (peak(i)<0) || (peak(i)>max(freq)) ) continue;
      peak_width(i) = clip(peak_width(i),freq(2)/10.,); 
      sigma = peak_width(i)/2.35;
      psd = exp( - (freq-peak(i))^2./(2*sigma^2.));
      psd = psd/sum(psd)*peak_rms(i)^2.;
      psdall += psd;
    } 
  }
  
  psdall(1) = 0.; // null zero freq component

  // add negative part:
  psdtot = _(psdall,psdall(2:-1)(::-1));
  norm = sum(psdtot);
  // amplitude
  amp = sqrt(psdtot);
  // generate random phase:
  pha = random(numberof(amp))*2*pi;
  // do the fourier transform:
  tmp = array(complex,numberof(amp));
  tmp.re = amp*cos(pha);
  tmp.im = amp*sin(pha);
  ts = float(fft(tmp,1));
  // post normalize:
  if ((normts=sum(ts^2))!=0.) ts = ts * sqrt(norm/normts);
  ts = ts*sqrt(length)/sqrt(2.); // don't ask, but it works and kinda make sense

  write,format="rms of time serie = %.3f milliarcsec\n",ts(rms)*1000.;
  
  if (sim.debug>=2) {
    fma; 
    plsys,2;
    logxy,0,0;
    plh,sqrt(psdall(2:)),freq(2:);
    xytitles,"Frequency [Hz]","sqrt(PSD) (one axis)";
    pltitle,swrite(format="Vibrations PSD (rms=%.1fmas)",ts(rms)*1000);
    limits,square=0;
    limits;
    limits,-10,freq(0)+10;
    hitReturn;
    logxy,0,0;
    limits;
  }

  return ts;    
}


func find_examples_path(void)
{
  parpath="./:"+pathform(_(Y_USER,Y_SITES,Y_SITE));
  tmp = find_in_path("sh6x6.par",takefirst=1,path=parpath);
  if (tmp==[]) tmp=find_in_path("data/sh6x6.par",takefirst=1,path=parpath);
  if (tmp==[]) tmp=find_in_path("share/yao/examples/sh6x6.par",takefirst=1,path=parpath);
  if (tmp==[]) {
    parpath="/usr/share/doc/yorick-yao/examples/";
    tmp = find_in_path("sh6x6.par",takefirst=1,path=parpath);
  }
  if (tmp==[]) {
    write,"Couldn't find the examples directory";
    return [];
  } else write,format="Found examples directory: %s\n",dirname(tmp);
  return dirname(tmp);
}

func find_doc_path(silent=)
{
  parpath="./:"+pathform(_(Y_USER,Y_SITES,Y_SITE));
  tmp = find_in_path("aosimul.html",takefirst=1,path=parpath);
  if (tmp==[]) tmp=find_in_path("doc/aosimul.html",takefirst=1,path=parpath);
  if (tmp==[]) tmp=find_in_path("share/yao/doc/aosimul.html",takefirst=1,path=parpath);
  if (tmp==[]) {
    parpath="/usr/share/doc/yorick-yao/doc/";
    tmp = find_in_path("aosimul.html",takefirst=1,path=parpath);
  }
  if (tmp==[]) {
    if (!silent) write,"Couldn't find the doc directory";
    return [];
  } else {
    if (!silent) write,format="Found doc directory: %s\n",dirname(tmp);
    return dirname(tmp);
  }
}


func yaodoc(void)
/* DOCUMENT
   will find the location of the documentation directory and open the
   manual in a browser (default browser for osx, firefox for linux).
   SEE ALSO:
 */

{
  docpath = find_doc_path(silent=1);
  if (os_env=="darwin") {
    system,"open "+docpath+"/manual.html";
  } else if (os_env=="linux") {
    system,"firefox "+docpath+"/manual.html &";
  }
}
