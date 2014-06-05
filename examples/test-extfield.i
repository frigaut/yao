require,"yao.i";
parfile = "sh6x6.par";
nois = 0;
skymag = 8;
gsmag = 10.5;
niter = 200;
disp = 10;
wfs_disp_crop_edges = 1;
alls = [];

aoread,parfile; loop.niter=niter;
wfs(1).fssize=1.6;
wfs(1).noise = nois; 
wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
aoinit,disp=1,clean=1;
fits_write,"imat1.fits",iMat,overwrite=1;
aoloop,disp=disp;
ran1init;
go,all=1;
grow,alls,avg(strehl);

// aoread,parfile; loop.niter=niter;
// wfs(1).fssize=1.6;
// wfs(1).spotpitch=8;
// wfs(1).noise = nois;
// wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
// aoinit,disp=1,clean=1;
// fits_write,"imat2.fits",iMat,overwrite=1;
// aoloop,disp=disp;
// ran1init;
// go,all=1;
// grow,alls,avg(strehl);

aoread,parfile; loop.niter=niter;
wfs(1).extfield=4.; 
wfs(1).spotpitch=20;
wfs(1).fssize=1.6;
wfs(1).noise = nois;
wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
aoinit,disp=1,clean=1;
fits_write,"imat3.fits",iMat,overwrite=1;
aoloop,disp=disp;
ran1init;
go,all=1;
grow,alls,avg(strehl);

aoread,parfile; loop.niter=niter;
wfs(1).extfield=8.; 
wfs(1).spotpitch=10;
wfs(1).fssize=1.6;
wfs(1).noise = nois;
wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
aoinit,disp=1,clean=1;
fits_write,"imat4.fits",iMat,overwrite=1;
aoloop,disp=disp;
ran1init;
go,all=1;
grow,alls,avg(strehl);

aoread,parfile; loop.niter=niter;
wfs(1).extfield=10.; 
wfs(1).spotpitch=24;
wfs(1).fssize=1.6;
wfs(1).noise = nois;
wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
aoinit,disp=1,clean=1;
fits_write,"imat5.fits",iMat,overwrite=1;
aoloop,disp=disp;
ran1init;
go,all=1;
grow,alls,avg(strehl);

// with sky
aoread,parfile; loop.niter=niter;
wfs(1).extfield=10.; 
wfs(1).spotpitch=24;
wfs(1).fssize=3.2;
wfs(1).noise = 1;
wfs(1).gsmag = gsmag; wfs(1).skymag = skymag;
aoinit,disp=1,clean=1;
aoloop,disp=disp;
ran1init;
go,all=1;


alls;

dev = (alls-avg(alls))/avg(alls); dev;
if (anyof(dev>0.001)) write,"There seems to be an issue with extfield."; \
  else write,"All great with extfield!";
