//----------------------------------------------------

func make_pzt_dm(nm,&def,disp=)
  /* DOCUMENT function make_pzt_dm2(dm_structure,disp=)
     the influence functions are in microns per volt.
  */
{
  gui_progressbar_frac,0.;
  gui_progressbar_text,swrite(format="Computing Influence Functions for DM#%d",nm);
  coupling=dm(nm).coupling;

  // best parameters, as determined by a multi-dimensional fit
  // (see coupling3.i)
  a=[4.49469,7.25509,-32.1948,17.9493];
  p1 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [2.49456,-0.65952,8.78886,-6.23701];
  p2 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [1.16136,2.97422,-13.2381,20.4395];
  irc = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  if (sim.debug>=2) write,format="p1=%f  p2=%f  ir=%f\n",p1,p2,irc;

  dim   = dm(nm)._n2-dm(nm)._n1+1;
  size  = sim._size;
  nxact = dm(nm).nxact;
  cobs  = tel.cobs;
  cent  = sim._cent;
  pitch = dm(nm).pitch;
  /*
    ir  = pitch*1.2;
    ir  = pitch*1.46;
    ir  = pitch*1.65;
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

  if (dm(nm).xflip) cub(,,1) = cub(::-1,,1);
  if (dm(nm).yflip) cub(,,2) = cub(,::-1,1);

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
//  x     = xy(,,2); y = xy(,,1);
  x     = xy(,,1); y = xy(,,2);
  def = array(float,dim,dim,nvalid);

  dm(nm)._x  = &(cubval(,1));
  dm(nm)._y  = &(cubval(,2));

  x = x(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);
  y = y(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);

  if (sim.verbose != 0) {write,format="\nCreating Influence function for actuator #%s","";}

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  for (i=1;i<=nvalid;i++) {
    if (sim.verbose != 0) write,format="%d ",i;
    if (dm(nm).irexp==1) {
      irfact = dm(nm).irfact;
      tmp = sqrt( ((x-cubval(i,1))/ir*irfact)^2.+((y-cubval(i,2))/ir*irfact)^2. );
      def(,,i)   = exp(-(tmp)^1.5);
    } else if (dm(nm).irexp==2) {
      //IF fitted from Hadamard experimental iMat:
      //      a_had = [0.2506,8.37,2.24497,26.2,0,0];//BETTER SET OF PARAM !!!
      a_had = [26.2,8.37]/8.*dm(nm).pitch;

      // make sure which sinc we're using:
      if (abs(sinc(1.))<1e-10) fact=1.; else fact=pi;

      def(,,i)= (sinc(fact * sqrt((x-cubval(i,1))^2.)/a_had(1))* \
                 sinc(fact * sqrt((y-cubval(i,2))^2.)/a_had(1))* \
                 exp(-((x-cubval(i,1))/a_had(2))^2.              \
                     -((y-cubval(i,2))/a_had(2))^2. ));
                     
      } else {
        tmpx       = clip(abs((x-cubval(i,1))/ir),1e-8,2.);
        tmpy       = clip(abs((y-cubval(i,2))/ir),1e-8,2.);
        tmp        = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*        \
                     (1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
        def(,,i)   = tmp*(tmpx <= 1.)*(tmpy <= 1.);
    }

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
    piston=def(,,sum)*ipupil(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2);
    tv,piston;
  }

  if (dm(nm)._puppixoffset!=[]) {
    // ok, so for now we'll do the following:
    // influence functions are actually shifted in comp_dm_shape,
    // so that the offset can be used for all types of DM.
    // here we also propagate on the IF coordinates, for completeness.
    // this may seem weird to someone that look at the defs later,
    // and compare to the coordinates, but, ok...
    *dm(nm)._x += dm(nm)._puppixoffset(1)
    *dm(nm)._y += dm(nm)._puppixoffset(2)
  }

  clean_progressbar;
  return def;
}

//----------------------------------------------------
func make_pzt_dm_elt(nm,&def,disp=)

  /* DOCUMENT function make_pzt_dm_elt(dm_structure,disp=)
     the influence functions are in microns per volt.
     same as make_pzt_dm but returns only local IF and
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

  if (sim.debug>=2) write,format="p1=%f  p2=%f  ir=%f\n",p1,p2,irc;

  dim   = dm(nm)._n2-dm(nm)._n1+1;
  size  = sim._size;
  nxact = dm(nm).nxact;
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
  if (dm(nm).xflip) cub(,,1) = cub(::-1,,1);
  if (dm(nm).yflip) cub(,,1) = cub(,::-1,1);

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
  dm(nm)._x  = &(cubval(,1));
  dm(nm)._y  = &(cubval(,2));
  dm(nm)._i1  = &(int(long(cubval(,1)-smallsize/2+0.5)-dm(nm)._n1));
  dm(nm)._j1  = &(int(long(cubval(,2)-smallsize/2+0.5)-dm(nm)._n1));

  def = def(,,-)*array(1.f,dm(nm)._nact)(-,-,);

  if (dm(nm)._puppixoffset!=[]) {
    // see comment above in make_pzt_dm
    *dm(nm)._x += dm(nm)._puppixoffset(1)
    *dm(nm)._y += dm(nm)._puppixoffset(2)
  }


  // look for extrapolation actuator stuff in v1.0.8 if needed

  fact = dm(nm).unitpervolt/max(def);
  def = float(def*fact);
  dm(nm)._def = &def;

  return def;
}
//----------------------------------------------------
func make_kl_dm(nm,&def,disp=)

  /* DOCUMENT function make_kl_dm,dm_number,ActIF,disp=
   */
{
  require,"yaokl.i";
  gui_progressbar_frac,0.;
  gui_progressbar_text,swrite(format="Computing Influence Functions for DM#%d: KL",nm);
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  //  nkllow = dm(nm).nklfiltered;
  nkllow = 1;
  nkl = dm(nm).nkl;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  gsdist = sqrt((abs(wfs.gspos)^2.)(sum,));
  patchDiam = long(ceil((sim.pupildiam+2*max(gsdist)*
                         4.848e-6*abs(dm(nm).alt)/psize)/2)*2);

  //  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;

  if (dm(nm).alt==0) {
    // enforce pupil to be system pupil.
    i1 = sim._size/2 - sim.pupildiam/2+1;
    i2 = sim._size/2 + sim.pupildiam/2;
    outpup = ipupil(i1:i2,i1:i2);
    //patchDiam should be good.
    write,format="KL: PatchDiam = %d, sim.pupildiam=%d\n",patchDiam,sim.pupildiam;
    cobs = tel.cobs;
  } else {
    outpup = [];
    cobs = 0.;
    patchDiam += 2; // margin
  }

  kl = float(make_kl(nkl,patchDiam,varkl,outbas,outpup,oc=cobs,nr=128));

  // order them in a similar order as zernike:
  kl = order_kls(kl,patchDiam,upto=20);

  def = array(float,dim,dim,nkl-nkllow+1);

  n1 = dim/2-patchDiam/2+1;
  n2 = n1+patchDiam-1;
  for (i=nkllow;i<=nkl;i++) {
    def(n1:n2,n1:n2,i-nkllow+1) = kl(,,i);
    if (disp == 1) {fma; pli,def(,,i-nkllow+1);}
  }
  if (sim.verbose>=1) {write,format="Number of KL :%d\n",nkl-nkllow+1;}

  // the KL are normalized so that rms over surface = 1 unit
  // meaning the TT go from -2 to 2.

  // we'll use the same normalization as zernike, just for consistency

  // current: over tel.diam, we have 4 units, and a rms of 1 we want a
  // rms of 957nm = 0.957microns

  def = def * 0.957f * float(dm(nm).unitpervolt);

  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  return def;
}


//----------------------------------------------------
func make_zernike_dm(nm,&def,disp=)

  /* DOCUMENT function make_zernike_dm,dm_structure,ActIF,disp=
     modified 2004jan22 to have scaled as tip-tilt (e.g.
     1 arcsec/volt).
   */
{
  gui_progressbar_frac,0.;
  gui_progressbar_text,swrite(format="Computing Influence Functions for DM#%d: Zernikes",nm);
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  nzer  = dm(nm).nzer;
  minzer = dm(nm).minzer;
  cobs  = tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  // below: bug discovered 2009mar24: << REDO mcao matrices
  // was using linear distance (abs(wfs.gspos), not working), not XY !!!
  gsdist = sqrt((abs(wfs.gspos)^2.)(sum,));
  patchDiam = sim.pupildiam+2*max(gsdist)*4.848e-6*abs(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;

  def = array(float,dim,dim,nzer-minzer+1);

  for (i=1;i<=(nzer-minzer+1);i++) {
    def(,,i) = zernike_ext(i+minzer-1);
    if (disp == 1) {fma; pli,def(,,i);}
    gui_progressbar_frac,float(i)/(nzer-minzer+1.);
  }
  if (sim.verbose>=1) {write,format="Number of zernike :%d\n",nzer;}

  // normalization factor: one unit of tilt gives 1 arcsec:
  z2 = zernike_ext(2);
  current = z2(dim/2,dim/2)-z2(dim/2-1,dim/2);
  fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;

  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  clean_progressbar;
  return def;
}

//----------------------------------------------------
func make_dh_dm(nm,&def,disp=)

  /* DOCUMENT function make_diskharm_dm,dm_structure
     adapted on 2010jun from the modal zernike dm function above.
   */
{

  gui_progressbar_frac,0.;
  gui_progressbar_text,swrite(format="Computing Influence Functions for DM#%d: disk harmonics",nm);
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  cobs  = tel.cobs;
  ndh   = dm(nm).ndh; // create this variable in dm structure
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  gsdist = sqrt((abs(wfs.gspos)^2.)(sum,));
  patchDiam = sim.pupildiam+2*max(gsdist)*4.848e-6*abs(dm(nm).alt)/psize;

  def = float(make_diskharmonic(dim,patchDiam,ndh,xc=cent-dm(nm)._n1+1,yc=cent-dm(nm)._n1+1));

  if (sim.verbose>=1) {write,format="Number of DH modes :%d\n",ndh;}

  // I am not sure if the normalization factor is correct for DH, but I am leaving it! (aurea)
  // normalization factor: one unit of tilt gives 1 arcsec:
  // current = def(dim/2,dim/2,3)-def(dim/2-1,dim/2,3);
  // fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;
  fact = dm(nm).unitpervolt;

  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);  // This is equal to ndh
  dm(nm)._def = &def;

  clean_progressbar;
  return def;
}

//----------------------------------------------------
func make_tiptilt_dm(nm,&def,disp=)

  /* DOCUMENT function make_tiptilt_dm,dm_structure,ActIF,disp=
     adapted from makeZernikeIF
     modified 2004jan22 to make it normalized at 1"
   */
{
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  nzer  = 2;
  cobs  = tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*abs(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;

  def = array(float,dim,dim,nzer);

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

func make_curvature_dm(nm,&def,disp=,cobs=)

  /* DOCUMENT:
     func make_curvature_dm,dm_structure,&def,disp=
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

  dim = sim._size;
  pupd  = sim.pupildiam;
  psize = tel.diam/sim.pupildiam;  // pixel in meter

  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*abs(dm(nm).alt)/psize;

  SupportRadius = 2.2;
  CompDim = dim*2;

  NRing = sum(*(dm(nm).nelperring) != 0); // Number of Rings
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
func make_aniso_dm(nm,&def,disp=)

  /* DOCUMENT function make_aniso_dm,dm_structure,ActIF,disp=
     2004jan22: implemented normalization as for zernikeIF,
     i.e. based on the same amplitude tip that gives 1"
  */
{
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  cobs  = tel.cobs;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  patchDiam = sim.pupildiam+2*max(abs(wfs.gspos))*
    4.848e-6*abs(dm(nm).alt)/psize;

  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;

  def = array(float,dim,dim,3);

  for (i=1;i<=3;i++) {
      def(,,i) = zernike_ext(i+3);
      if (disp == 1) {fma; pli,def(,,i);}
    }
  if (sim.verbose>=1) {write,format="Number of Anisoplanatism modes :%d\n",3;}

  // normalization factor: see make_zernike_dm and make_tiptilt_dm
  tip = zernike_ext(2);
  current = tip(dim/2,dim/2,1)-tip(dim/2-1,dim/2,1);
  fact = (dm(nm).unitpervolt*tel.diam/sim.pupildiam)*4.848/current;

  def = float(def*fact);
  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  return def;
}



//----------------------------------------------------
func project_aniso_dm(nmaniso,nmlow,nmhigh,disp=)
/* DOCUMENT func project_aniso_dm(nmaniso,nmlow,nmhigh,disp=)
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
   SEE ALSO: make_aniso_dm
 */
{
  extern alow,ahigh,comaniso;

  // we address here the zero altitude layer. The pupil is well defined.
  // cut a ipupil of the appropriate size:
  puplow = ipupil(dm(nmlow)._n1:dm(nmlow)._n2,dm(nmlow)._n1:dm(nmlow)._n2);
  w = where(puplow);
  // this transform def into a #spatial_point x nact array and retains only
  // the spatial point inside the pupil
  if (dm(nmlow).elt == 1) {
    n1 = dm(nmlow)._n1;
    n2 = dm(nmlow)._n2;
    sizedef=n2-n1+1;
    tabdef=array(float,sizedef,sizedef,sum(dm(nmlow)._nact));
    command = array(float,dm(nmlow)._nact);
    for (i=1;i<=dm(nmlow)._nact;i++) {
      command *= 0.0f;
      command(i) = 1.0f;
      tabdef(,,i) = comp_dm_shape(nmlow,&command);
    }
    def=tabdef(*,)(w,);
  }
  else {
    def = (*dm(nmlow)._def)(*,)(w,);
  }
  // compute the IF covariance matrix
  defcov = def(+,)*def(+,);

  // now look at the anisoplanatism modes:
  // same, extract ipupil of appropriate dimension
  pupaniso = ipupil(dm(nmaniso)._n1:dm(nmaniso)._n2,\
                    dm(nmaniso)._n1:dm(nmaniso)._n2);
  w = where(pupaniso);
  // this transform def into a #spatial_point x nact array and retains only
  // the spatial points inside the pupil
  defa = -(*dm(nmaniso)._def)(*,)(w,);
  // compute the product act * mode:
  anisoproj = def(+,)*defa(+,);

  // command vector (matrices, 3 modes) to apply to DM to get a given mode
  alow = LUsolve(defcov,anisoproj);

  // display:
  if (disp) {
    for (i=1;i<=3;i++) {
      if (dm(nmlow).elt == 1) {
        tv,tabdef(,,+)*alow(+,i)*puplow;
        }
      else tv,(*dm(nmlow)._def)(,,+)*alow(+,i)*puplow;
      hitReturn;
    }
  }

  // Now, the altitude DM:
  // it's basically the same thing, except now there is no well-defined pupil.
  // So we define here the pupil as the area which is controllable by the actuators.
  // that should be perfectly acceptable as the is the only area which will be seen
  // by any beam.

  if (dm(nmhigh).elt == 1) {
    n1 = dm(nmhigh)._n1;
    n2 = dm(nmhigh)._n2;
    sizedef=n2-n1+1;
    tabdef=array(float,sizedef,sizedef,sum(dm(nmhigh)._nact));
    command = array(float,dm(nmhigh)._nact);
    for (i=1;i<=dm(nmhigh)._nact;i++) {
      command *= 0.0f;
      command(i) = 1.0f;
      tabdef(,,i) = comp_dm_shape(nmhigh,&command);
    }
    puphigh = (tabdef)(,,sum);
    if (max(-puphigh)>max(puphigh)) puphigh = -puphigh;
    puphigh = (puphigh > 0.8*max(puphigh));
    w = where(puphigh);
    def=tabdef(*,)(w,);
  } else {
    //puphigh = (*dm(nmhigh)._def)(,,sum);
    puphigh = comp_dm_shape(nmhigh,&(array(1.0f,dm(nmhigh)._nact)));
    if (max(-puphigh)>max(puphigh)) puphigh = -puphigh;
    puphigh = (puphigh > 0.8*max(puphigh));
    w = where(puphigh);
    def = (*dm(nmhigh)._def)(*,)(w,);
  }

  defcov = def(+,)*def(+,);

  pupaniso = array(float,[2,sim._size,sim._size]);
  pupaniso(dm(nmhigh)._n1:dm(nmhigh)._n2,dm(nmhigh)._n1:dm(nmhigh)._n2)=puphigh;
  pupaniso = pupaniso(dm(nmaniso)._n1:dm(nmaniso)._n2,\
                      dm(nmaniso)._n1:dm(nmaniso)._n2);
  w = where(pupaniso);
  defa = (*dm(nmaniso)._def)(*,)(w,);
  anisoproj = def(+,)*defa(+,);

  ahigh = LUsolve(defcov,anisoproj);

  if (disp || (sim.debug == 2)) {
    for (i=1;i<=3;i++) {
      if (dm(nmhigh).elt == 1) {
        tv,tabdef(,,+)*ahigh(+,i)*puphigh;
        }
      else tv,(*dm(nmhigh)._def)(,,+)*ahigh(+,i)*puphigh;
      hitReturn;
    }
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


func make_segmented_dm(nm,&def,disp=)
{
  dim   = dm(nm)._n2-dm(nm)._n1+1;
//  dim = sim._size;
  cent = sim._cent-dm(nm)._n1+1;

  // make large, unfiltered segment map:
  map = make_seg_hexa_grid(dm(nm).pitch,dm(nm).nxseg,dim,x,y,cent=cent);
  // keep only the one within the area of interest:
  if (dm(nm).fradius) f_rad = dm(nm).fradius;
  else f_rad = dm(nm).pitch*(dm(nm).nxseg)/2.;

  map  = filt_seg_hexa_grid(map,x,y,f_rad,cent=cent);
  // repack/renumber
  map = renum_int_array(map);

  // number of segments
  nseg = max(map);
  // 3 degrees of freedom per segment (piston, tip and tilt)
  nact = 3*nseg;
  // allocate influence function data cube
  def = array(0.0f,[3,dim,dim,nact]);
  // build influence functions
  k = 1;
  for (i=1; i<=nseg; i++) {
    // piston
    def(,,k++) = float(map==i);
    // tip
    tmp        = float(def(,,k-1)*\
                    tip1arcsec(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2));
    w          = where(def(,,k-1));
    tmp(w)    -= avg(tmp(w));
    def(,,k++) = tmp;
    // tilt
    tmp        = float(def(,,k-2)*\
                    tilt1arcsec(dm(nm)._n1:dm(nm)._n2,dm(nm)._n1:dm(nm)._n2));
    w          = where(def(,,k-2));
    tmp(w)    -= avg(tmp(w));
    def(,,k++) = tmp;
    if (disp) tv,def(,,k-1);
  }

  dm(nm)._nact = nact;
  dm(nm)._def = &def;

  return map;
}


func make_seg_hexa_grid(pitch,nxseg,dim,&x,&y,cent=,rotby=)
/* DOCUMENT func make_seg_hexa_grid(pitch,nxseg,dim,&x,&y,cent=,rotby=)
   pitch = Segment center to segment center pitch [pixels]
   nxseg = Number of segment in long axis (X) diameter
           nxseg must be integer, odd or even
   dim   = Size of final maps (optional, default to (nxseg+2)*pitch)
   rotby = Angle in degrees to rotate maps by (optional)
   SEE ALSO: filt_seg_hexa_grid, renum_int_array
   Typical sequence:
   m   = make_seg_hexa_grid(20,9,,x,y) // make segment grid
   fm  = filt_seg_hexa_grid(m,x,y,20*9/2)
   you can call above as many time as you like until segment pattern
   outer figure is acceptable. then:
   map = renum_int_array(fm)
 */

{
  tic;
  if (!dim) dim = (long((nxseg+2)*pitch+1)/2)*2;
  if (odd(dim)) write,"Warning: make_segments_hexa(): dim is odd";
  if (!cent) cent=dim/2;

  // compute coordinates of segment centers
  // we need to oversize as Y dimension is compressed
  // (ypitch = xpitch / (sqrt(3)/2.)
  nyseg = long(nxseg / (sqrt(3.)/2.));
  if (odd(nyseg-nxseg)) nyseg++;
  // let's still oversize it a bit:
  nover = 4;
  xy = indices(nyseg+nover);
  // center it at (nxseg+nover)/2;
  xy = xy - 1 - (nyseg+nover)/2;
  // now let's shift X by 0.5 for all odd Y:
  xy(,,1) += odd(abs(xy(,,2)))*0.5;
  // adjust for the case nxseg = even (then no segment at exact center)
  xy += even(nxseg)*0.5;

  // scale X and Y by proper distances:
  x = xy(,,1)*pitch;
  y = xy(,,2)*pitch*sqrt(3.)/2.;

  // make vectors:
  x = x(*);
  y = y(*);

  if (debug) {
    write,format="nxseg = %d, nyseg = %d\n",nxseg,nyseg;
    fma;
    plp,y,x,symbol=4,size=0.3;
    limits,square=1;
    plmargin;
    hitReturn;
  }

  // rotate, if needed:
  if (rotby) {
    xy2 = mrot(rotby)(+,)*[x,y](,+);
    x = xy2(1,);
    y = xy2(2,);
  }

  // get rid of segments that will be completely, for sure,
  // outside of a circular aperture?
  //d = sqrt(x^2.+y^2.);
  //w = where(d<=(nxseg*pitch/2.+pitch/2.));
  //if (numberof(w)==0) error,"Problem when filtering segments";
  //x = x(w);
  //y = y(w);
  nseg = numberof(x);

  if (debug) {
    fma;
    plp,y,x,symbol=4,size=0.3;
    ang = span(0.,2*pi,100);
    plg,(nxseg*pitch)/2.*sin(ang),(nxseg*pitch)/2.*cos(ang);
    limits,square=1;
    plmargin;
    hitReturn;
  }
  if (debug) write,(t1=tac());

  x += cent;
  y += cent;

  // build the segments at requested locations:
  // To not go overboard with memory allocation, I chose
  // to do it by chunk. There might be other ways, but this
  // is cheap and quick:
  maxram = 256; // in MB
  maxdim = long(sqrt(maxram*1e6/(nseg*4)));
  if (debug) write,format="Max dim = %d\n",maxdim;

  if (maxdim>=dim) {
    cube = array(0.0f,[3,dim,dim,nseg]);
    for (i=1;i<=nseg;i++) {
      cube(,,i) = float(dist(dim,xc=x(i),yc=y(i)));
    }
    map = cube(,,mnx);
  } else {
    sdim = dim;
    while (sdim>maxdim) sdim = long(ceil(sdim/2.));
    ndim = long(ceil(dim*1./sdim));
    map = array(long,[2,ndim*sdim,ndim*sdim]);
    cube = array(0.0f,[3,sdim,sdim,nseg]);
    for (i=1;i<=ndim;i++) {
      if (debug) write,format="%d out of %d\n",i,ndim;
      for (j=1;j<=ndim;j++) {
        for (k=1;k<=nseg;k++) {
          cube(,,k) = float(dist(sdim,xc=-(i-1)*sdim+x(k),
                                      yc=-(j-1)*sdim+y(k)));
        }
//        error;
        map(1+(i-1)*sdim:i*sdim,1+(j-1)*sdim:j*sdim) = cube(,,mnx);
        if (debug) tv,cube(,,mnx);
      }
    }
    map = map(1:dim,1:dim);
  }

  if (debug) write,(t2=tac())-t1;

  return map;
}



func filt_seg_hexa_grid(map,&x,&y,filter_radius,cent=)
/* DOCUMENT
   func filt_seg_hexa_grid(map,x,y,filter_radius)
   Filter segments that are outside of area of interest

   map = map of segments (pixel value = segment number it
     belongs to). Generally built by make_seg_hexa_grid()
   x and y = coordinates of segment centers (in pixel, w.r.t
     center of map). Generally returned by make_seg_hexa_grid()
   filter_radius = segment will be filtered if distance to center
                   is greater than this value (pixels)
   cent = center of coodinates (pixels, optional, default dim/2)

   SEE ALSO: make_seg_hexa_grid, renum_int_array
 */
{
  if (!cent) cent=dimsof(map)(2)/2;
  fmap = map;
  if (filter_radius) {
    d = sqrt((x-cent)^2.+(y-cent)^2.);
    // keep only valid x and y
    w = where(d<=filter_radius);
    if (numberof(w)==0) error,"Problem: 0 segment selected";
    x = x(w);
    y = y(w);
    // zero out segments outside of defined radius:
    w = where(d>filter_radius);
    if (numberof(w)>0) {
      for (i=1;i<=numberof(w);i++) {
        ww = where(fmap==w(i));
        if (numberof(ww)) fmap(where(fmap==w(i))) = 0;
      }
    }
  }

  return fmap;
}



func renum_int_array(map)
/* DOCUMENT func renum_int_array(map)
  // re-assign segment numbers to get rid of missing subap/segments
  // in returned map, first segment/subap is 1, then 2, ... until Nmax
   SEE ALSO:
 */
{
  rmap = map;
  m = max(rmap);
  l = 1;
  for (i=1;i<=m;i++) {
    w = where(rmap==i);
    if (numberof(w)) rmap(w) = l++;
  }

  return rmap;
}



