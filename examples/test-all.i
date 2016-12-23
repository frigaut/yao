require,"yao.i";
write,format="YAO_SAVEPATH = \"%s\"\n",YAO_SAVEPATH;
write,format="cwd = \"%s\"\n",get_cwd();

write,"CREATING PHASE SCREENS";
if (!open(Y_USER+"data/screen1.fits","r",1)) {
  mkdirp,Y_USER+"data/screen";
  create_phase_screens,1024,256,prefix=Y_USER+"data/screen";
 }

f = findfiles("test*.par");
f = f(sort(f));
if (get_argv()(0)=="bench") f = findfiles("*-bench.par");

if (get_argv()(0)=="toy") f = findfiles("toy*.par");

struct perf_s{string parfile; string name; float itps; string display; float strehl; float lambda;};
//perf = array(perf_s,numberof(f)+1);
perf = [];

smdebug=0;

write,"LOOPING ON TEST*.PAR";
for (i=1;i<=numberof(f);i++) {
  write,format="\n\nTesting %s\n\n",f(i);
  aoread,f(i);
  sim.verbose = 1;
  sim.debug=0;
  // SVIPC settings:
  sim.svipc = 0;
  wfs.svipc = 0;
  // wfs.svipc = 2;
  // sim.svipc = 1+2;
  // sim.svipc=1+2+4;
  disp = 10;
  if (!strmatch(f(i),"fast")) loop.niter = 500;
  else disp=0; // let's not display for the fast demo
  aoinit,disp=disp,clean=1;
  random_seed,0.5; ran1init; 
  aoloop,disp=disp,no_reinit_wfs=1,controlscreen=10*(i==2)*(disp!=0);
  go,all=1;
  grow,perf,perf_s();
  perf(0).parfile = f(i);
  perf(0).name=sim.name;
  perf(0).itps = iter_per_sec;
  perf(0).display = (disp?"ON":"OFF");
  perf(0).strehl = strehl(1,0); // first position, last lambda
  perf(0).lambda = (*target.lambda)(0);
  // same without displays
  disp=0;
  random_seed,0.5; ran1init; 
  aoloop,disp=disp,no_reinit_wfs=1;
  go,all=1;
  grow,perf,perf_s();
  perf(0).parfile = f(i);
  perf(0).name=sim.name;
  perf(0).itps = iter_per_sec;
  perf(0).display = (disp?"ON":"OFF");
  perf(0).strehl = strehl(1,0); // first position, last lambda
  perf(0).lambda = (*target.lambda)(0);
  if (sim.svipc||anyof(wfs.svipc)) status = quit_forks();
  // hitReturn;
}

//~ write,format="\n\nTesting %s (no display)\n\n",f(0);
//~ aoread,f(0); 
//~ sim.verbose = 1;
//~ sim.debug=0;
//~ loop.niter = 500;
//~ aoinit,disp=0;
//~ aoloop,disp=0;
//~ go,all=1;
//~ perf(0).parfile = f(0);
//~ perf(0).name=sim.name; 
//~ perf(0).itps = iter_per_sec;
//~ perf(0).display = "OFF";
//~ perf(0).strehl = strehl(1,0); // first position, last lambda
//~ perf(0).lambda = (*target.lambda)(0);
write,"\nSUCCESS: ALL TESTS OK.\n\nPerformance Summary:";
f = open(swrite(format="test-all_v%s.dat",yaoVersion),"a");
write,f,format="%s\n","----------------------------------------------------------";
write,f,format="%s\n",timestamp();
write,format="%-12s %-34s %-8s%-10s%-10s\n","Parfile","Name","iter/s","Display?","Strehl";
write,f,format="%-12s %-34s %-8s%-10s%-10s\n","Parfile","Name","iter/s","Display?","Strehl";
for (i=1;i<=numberof(perf);i++) {
  if (strlen(perf(i).name)>34) name = strpart(perf(i).name,1:31)+"...";
  else name = perf(i).name;
  str = swrite(format="%-12s %-34s %-8.1f%-10s%.2f@%.2fmic",\
  strpart(strtok(perf(i).parfile,".")(1),1:12), \
  name,perf(i).itps,perf(i).display,perf(i).strehl,perf(i).lambda);
  write,format="%s\n",str;
  write,f,format="%s\n",str;
}
close,f;

