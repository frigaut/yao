require,"yao.i";
write,"CREATING PHASE SCREENS";
if (!open(Y_USER+"data/screen1.fits","r",1)) {
  mkdirp,Y_USER+"data/screen";
  CreatePhaseScreens,1024,256,prefix=Y_USER+"data/screen";
 }
f = findfiles("test*.par");
write,"LOOPING ON TEST*.PAR";
for (i=1;i<=numberof(f);i++) {
  write,format="\n\nTesting %s\n\n",f(i);
  pause,1000;
  aoread,f(i);
  sim.verbose = 1;
  loop.niter = 500;
  aoinit,disp=1,dpi=50,clean=1;
  aoloop,disp=10,controlscreen=10;
  for (k=1;k<=loop.niter;k++) go;
  after_loop;  // to wrap up the analysis and print out results
}
aoread,f(0); aoinit; aoloop;
write,"\nAll tests OK\n";
