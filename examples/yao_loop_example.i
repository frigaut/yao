require,"yao.i";
write,"CREATING PHASE SCREENS";
if (!open(Y_USER+"data/screen1.fits","r",1)) {
  CreatePhaseScreens,1024,256,prefix=YUSER+"data/screen";
 }

window,33,wait=1;

// read out parfile
aoread,"test3.par";
atm.dr0at05mic = 35; // be more gentle

// define vector on which we want to loop and final strehl array
gsmagv = [6,9,12];
gainv = [0.01,0.1,0.5,1.0];
strehlarray = array(0.,[2,numberof(gsmagv),numberof(gainv)]);

// loop on gsmag and gain
for (ii=1;ii<=numberof(gsmagv);ii++) {
  for (jj=1;jj<=numberof(gainv);jj++) {
    wfs(1).gsmag=gsmagv(ii);
    loop.gain=gainv(jj);
    // it's safer, but not always necessary, to call again
    // aoinit (here for gsmag). some parameters do not need it.
    aoinit,disp=1;
    aoloop,disp=1,controlscreen=10;
    plsys,1; animate,1; // double buffering: gives smoother displays
    for (kk=1;kk<loop.niter;kk++) go;
    plsys,1; animate,0; // turn off double buffering
    //after_loop;  // to wrap up the analysis and print out results
    strehlarray(ii,jj) = strehllp(0); // fill in result array
    // and display results as we go:
    window,33;
    fma;
    for (ll=1;ll<=ii;ll++) {
      plg,strehlarray(ll,),gainv,color=-ll-4;
      ylims=limits()(3:4); ymax=ylims(2); yspace=(ylims(2)-ylims(1))/15.;
      plt,swrite(format="gsmag=%d",gsmagv(ll)),0.011,ymax-yspace*(ll-1), \
        justify="LT",tosys=1,color=-ll-4;
    }
    logxy,1,0;
    xytitles,"Loop Gain",swrite(format="Strehl @ %.2fmicrons",(*target.lambda)(0));
    window,0;
  }
 }
