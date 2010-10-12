//======================================
// SVIPC functions for parallelization
//======================================

require,"svipc.i";
require,"yao_setnsync.i";

nshm    = 50;
nsem    = 50;


func init_keys(void)
{
  extern shmkey, semkey;
  
  if ( (sim.shmkey) || (sim.shmkey == 0) ) shmkey=0x0badcafe; \
  else shmkey = sim.shmkey;

  if ( (sim.semkey) || (sim.semkey == 0) ) semkey=0x0badbeef; \
  else semkey = sim.semkey;
}

status = init_keys();

func svipc_init(void)
{
  extern shm_init_done;

  if (!shm_init_done) {
    if (shm_init(shmkey,slots=nshm)==0) {
      sem_init,semkey,nums=nsem;
      if (sim.verbose>0) write,format="%s\n","SVIPC initialized";
    } // else already initialized by another process
  }
  shm_init_done = 1;
  shm_write,shmkey,"quit?",&([0]);
  
  // what's in semaphores? 
  // sem#          content
  // 0             trigger WFS child calculation
  // 1             ready from WFS child
  // 3             psf: trigger child PSFs calculation
  // 4             psf: notify parent PSFs ready
  // 20+ reserved for WFSs ||
}


func svipc_wfs_init(phase,ns)
{
  extern atm,wfs,dm,loop,mat,opt,tel,target,gs;
  extern pscreens,im,imav,mircube,cubphase;
  extern iMat,cMat,xposvec,yposvec;
  extern wfsxposcub,   wfsyposcub,   gsxposcub,   gsyposcub;
  extern dmwfsxposcub, dmwfsyposcub, dmgsxposcub, dmgsyposcub;
  extern optwfsxposcub,optwfsyposcub,optgsxposcub,optgsyposcub;
  extern statsokvec, sphase, bphase, imtmp, imphase;
  extern strehllp, strehlsp, itv, commb, errmb;
    
  //  extern nforks_per_wfs;

  if (sim.debug>=1) write,format="Entering svipc_wfs_init for WFS#%d\n",ns;
  
  // Just to make sure. no need for svipc in that case:
  if (wfs(ns).svipc<=1) return;
  
  // Bail out if init has already been done for this WFS
  if (wfs(ns)._svipc_init_done) return;

  // Make sure shm has been initialized
  if (!shm_init_done) status = svipc_init();

  //  if (nforks_per_wfs==[]) nforks_per_wfs = array(0,nwfs);
  //  if (numberof(nforks_per_wfs)!=nwfs)
  //    error,"You can't change the #wfs while using svipc, please restart yao";
  //  nforks_per_wfs(ns) = wfs(ns).svipc;
    
  // Initialize a few generic variable we're gonna need:
  // This on serves to indicate to the forks that a sync is needed
  shm_write,shmkey,swrite(format="sync_wfs%d_forks",ns),&([0]);
  // This on is used to quit the forks
  shm_write,shmkey,"quit_wfs_forks?",&([0]);

  // To avoid a shm_write/shm_read at each iter, I have opted to
  // share the variables (shm_var) fimage, phase and mesvec.
  // We need to create the slots (shm_write) first:
  // *wfs._fimage has necessarily been created in shwfs_init()
  shm_write,shmkey,swrite(format="wfs%d_fimage",ns),wfs(ns)._fimage;
  // phase exist as it has been passed as arg to sh_wfs()
  // write,format="%s: ","phase in svipc_wfs_init";  info,phase;
  shm_write,shmkey,swrite(format="wfs%d_phase",ns),&phase;
  // ... and we create mesvec:
  mesvec = array(float,2*wfs(ns)._nsub);
  shm_write,shmkey,swrite(format="wfs%d_mesvec",ns),&mesvec;

  
  // Compute which subapertures have to be dealt with by which forks:
  // Find subok for _shwfs_phase2spots()
  // The parallelization is for *spot* calculation, hence nsub4disp
  nsubsperforks = float(wfs(ns)._nsub4disp)/wfs(ns).svipc;
  tofork = long(floor((indgen(wfs(ns)._nsub4disp)-1)/nsubsperforks))+1;
  wfs(ns)._fork_subs = &array(int,[2,wfs(ns)._nsub4disp,wfs(ns).svipc]);
  
  // find yoffset, subok, ysize parameters for _shwfs_spots2slopes()
  wfs(ns)._fork_subs2 = &split_subok(ns,yoffset,ysize);
  wfs(ns)._fimny2 = &ysize;
  wfs(ns)._yoffset = &yoffset;

  if (sim.verbose>0) {
    write,format="WFS#%d, %d forks, %.1f subap/forks\n",
      ns,wfs(ns).svipc,nsubsperforks;
  }
  
  // usual thing: we can't fork() with windows open, so let's kill them
  wl = window_list();
  if (wl!=[]) for (i=1;i<=numberof(wl);i++) winkill,wl(i);

  extern svipc_wfs_ns;
  svipc_wfs_ns = ns; // we'll need that too

  // leave a trace that we've gone through here.
  wfs(ns)._svipc_init_done = 1;  
  
  for (nf=1;nf<=wfs(ns).svipc;nf++) {
    
    (*wfs(ns)._fork_subs)(,nf) = int(tofork==nf);

    if (nf<2) continue; // "nf=1" is main process.

    // fork:
    if (fork()==0) { // I'm the child

      // child: free resources
      atm = dm = mat = opt = tel = target = gs = [];
      pscreens = im = imav = mircube = cubphase = [];
      iMat = cMat = xposvec = yposvec = [];
      wfsxposcub = wfsyposcub = gsxposcub = gsyposcub = [];
      dmwfsxposcub = dmwfsyposcub = dmgsxposcub = dmgsyposcub = [];
      optwfsxposcub = optwfsyposcub = optgsxposcub = optgsyposcub = [];
      statsokvec = sphase = bphase = imtmp = imphase = [];
      strehllp = strehlsp = itv = commb = errmb = [];

      // start listening
      status = wfs_fork_listen(ns,nf);
    }
  }

  // restore windows if needed:
  if (wl!=[]) status = create_yao_window();
  return 0;
}

func wfs_fork_listen(ns,nf)
/* DOCUMENT wfs_fork_listen(ns,nf)
   Main loop for yao wfs forks.
   - Init all variables as in sh_wfs()
   - loop forever waiting for trigger from main
     and executing _shwfs_phase2spots and _shwfs_spots2slopes()
   - at each iter, check if sync is needed.
   SEE ALSO:
 */
{
  if (sim.debug>0) write,format="WFS#%d fork#%d listening\n",ns,nf;

  // init the internal variables just as in sh_wfs()
  pupd        = sim.pupildiam;
  size        = sim._size;
  nxsub       = wfs(ns).shnxsub(0);
  subsize     = int(pupd/nxsub);
  if (wfs(ns).npixpersub) subsize = wfs(ns).npixpersub;
  phasescale  = float(2*pi/wfs(ns).lambda);
  sdim        = long(2^ceil(log(subsize)/log(2)+1));
  sdimpow2    = int(log(sdim)/log(2));

  // to protect from a WFS sync, let's not put it in the wfs structure
  svipc_subok  = (*wfs(ns)._fork_subs)(,nf);
  svipc_subok2 = (*wfs(ns)._fork_subs2)(,nf);
  yoffset      = (*wfs(ns)._yoffset)(nf);
  fimny2       = (*wfs(ns)._fimny2)(nf);

  // if (phase!=[]) shm_unvar,phase;
  
  shm_var,shmkey,swrite(format="wfs%d_fimage",ns),ffimage;
  shm_var,shmkey,swrite(format="wfs%d_phase",ns),phase;
  shm_var,shmkey,swrite(format="wfs%d_mesvec",ns),mesvec;
  
  
  // then listen and execute ad libitum
  do {

    // wait for trigger:
    if (sim.debug>20) \
      write,format="fork: Waiting for trigger from main on sem ns=%d\n",20+4*(ns-1);
    sem_take,semkey,20+4*(ns-1);
    if (sim.debug>20) write,format="fork: gotten sem %d\n",20+4*(ns-1);
    
    // check if we have to quit:
    if (shm_read(shmkey,"quit_wfs_forks?")(1)) {
      if (sim.verbose>0) {
        write,format="WFS#%d, fork#%d quitting\n> ",ns,nf;
        //        maybe_prompt;
      }
      yorick_quit;
    }

    // sync if needed:
    status = sync_wfs_from_master(ns,nf);

    // write,format="%s ","in fork."; info,phase;
    
    // do our stuff: 
    err = _shwfs_phase2spots( pupsh, phase, phasescale,
             *wfs(ns)._tiltsh, int(size), *wfs(ns)._istart,
             *wfs(ns)._jstart, int(subsize), int(subsize), 
             wfs(ns)._nsub4disp, sdimpow2, wfs(ns)._domask, *wfs(ns)._submask, 
             *wfs(ns)._kernel, *wfs(ns)._kernels, *wfs(ns)._kerfftr,
             *wfs(ns)._kerffti, wfs(ns)._initkernels, wfs(ns)._kernelconv,
             *wfs(ns)._binindices, wfs(ns)._binxy, 
             wfs(ns)._rebinfactor, ffimage, svipc_subok,
             *wfs(ns)._imistart, *wfs(ns)._imjstart, 
             wfs(ns)._fimnx , wfs(ns)._fimny, 
             *wfs(ns)._fluxpersub, *wfs(ns)._raylfluxpersub,
             *wfs(ns)._skyfluxpersub, float(wfs(ns).darkcurrent*loop.ittime),
             int(wfs(ns).rayleighflag), 
             *wfs(ns)._rayleigh, wfs(ns)._bckgrdinit,
                              wfs(ns)._cyclecounter, wfs(ns).nintegcycles);

    // give trigger back:
    if (sim.debug>20) write,format="fork: giving trigger on sem %d\n",20+4*(ns-1)+1;
    sem_give,semkey,20+4*(ns-1)+1;
    
    sem_take,semkey,20+4*(ns-1)+2;
    if (sim.debug>20) write,format="fork: gotten sem %d\n",20+4*(ns-1)+2;
    
    threshold   = array(float,wfs(ns)._nsub4disp+1)+wfs(ns).shthreshold;
    
    err = _shwfs_spots2slopes( ffimage,
                *wfs(ns)._imistart2, *wfs(ns)._imjstart2,
                wfs(ns)._nsub4disp, wfs(ns).npixels,
                wfs(ns)._fimnx, fimny2, yoffset, 
                *wfs(ns)._centroidw, wfs(ns).shthmethod, threshold, *wfs(ns)._bias,
                *wfs(ns)._flat, wfs(ns).ron, wfs(ns).noise, 
                *wfs(ns)._bckgrdcalib, wfs(ns)._bckgrdinit, wfs(ns)._bckgrdsub,
                *wfs(ns)._validsubs, svipc_subok2, wfs(ns).nintegcycles,
                mesvec);
    

    sem_give,semkey,20+4*(ns-1)+3;

  } while (1);
}


func split_subok(ns,&yoffset,&ysize)
/* DOCUMENT split_subok(ns)
   Returns a matrix indicating which subap should be processed by
   which fork.
   wfs(ns).svipc rows.
   each rows has wfs._nsub4disp elements. In row N, 1 means
   it is to be processed by this process N, 0 means not process.
   The difficulty here is that in _shwfs_spots2slopes(), we add
   the noise to the image, so the subap for each forks have to
   span full rows.
   returns also the yoffsets and ysize (fimny) to be used in
   the call to _shwfs_spots2slopes()
   SEE ALSO:
 */
{
  nj = *wfs(ns)._jstart;
  //  d=_(ni(where(ni(dif))),ni(0));
  subsize  = sim.pupildiam/wfs(ns).shnxsub;
  nj = nj/subsize;
  nj = nj-min(nj);
  nt = int(floor(nj/(wfs(ns).shnxsub*1./wfs(ns).svipc)));

  subok = array(int,[2,wfs(ns)._nsub4disp,wfs(ns).svipc]);
  yoffset = ysize = array(int,wfs(ns).svipc);

  for (i=1;i<=wfs(ns).svipc;i++) {
    subok(,i)  = ((nt+1)==i);
    yoffset(i) = (*wfs(ns)._imjstart2)(where(subok(,i))(1));
    ysize(i)   = (*wfs(ns)._imjstart2)(where(subok(,i))(0)) - \
      yoffset(i)+wfs(ns).npixels;
  }
  return subok;
}

func svipc_wfs_profile(void)
{
  tfork = rdcols("/tmp/fork.res");
  tmain   = rdcols("/tmp/main.res");

  
}

func svipc_start_forks(void)
{
  extern iMat,cMat,dm,atm,wfs;
  extern fork_done;
  extern svipc_procname;
  extern all_svipc_procname;
  // sim.svipc bits (also see structure definition in yao_structure.i)
  //   bit 0 (1): parallelize WFS/DM
  //   bit 1 (2): parallelize PSFs

  if (fork_done) {
    write,"forks() already done, ignoring";
    return;
  }

  // usual thing: we can't fork() with windows open, so let's kill them
  wl = window_list();
  // if (anyof(wl==0)) { // that should be the yao window
  //   w0_dpi = window_geometry(0);
  // }
  if (wl!=[]) for (i=1;i<=numberof(wl);i++) winkill,wl(i);

  all_svipc_procname = [];
  
  // WFS CHILD
  // spawn just one process, the one taking care of wfsmes.
  // it could possibly spawn itself other process later (e.g.
  // for each of the WFS).
  if ((sim.svipc>>0)&1) {
    
    svipc_procname = "WFS";
    grow,all_svipc_procname,svipc_procname;
    
    if (fork()==0) { // I'm the child

      if (sim.verbose>0) \
        write,format="WFS child fork()ed with PID %d\n",getpid();
      
      // get rid of what we don't need
      iMat = cMat = [];
      for (i=1;i<=ndm;i++) dm(i)._def = &[];

      // start listening
      //set_idler,topwfs_listen;
      status = topwfs_listen();
      exit;
    }
  }

  // PSFs child
  if ((sim.svipc>>1)&1) {
    svipc_procname = "PSFs";
    grow,all_svipc_procname,svipc_procname;
    if (fork()==0) { // child
      write,format="PSFs child fork()ed with PID %d\n",getpid();
      // get rid of what we don't need
      //
      // start listening
      // set_idler,psf_listen;
      status = psf_listen();
      exit;
    }
  }
    
  fork_done = 1;
  
  // restore windows if needed:
  if (wl!=[]) status = create_yao_window();
  return 0;
}


yorick_quit = quit;

func quit
{
  shm_write,shmkey,"quit?",&([1]);
  shm_write,shmkey,"quit_wfs_forks?",&([1]);
// nforks = sum(clip(wfs.svipc-1,0,));
  for (i=0;i<=nsem;i++) sem_give,semkey,i,count=100;
  usleep,50;
  status = svipc_clean();
  yorick_quit;
}

func quit_forks(void)
{
  extern wfs,shm_init_done,fork_done;
  
  if ( (sim!=[]) && (sim.svipc) ) {
    // a previous yao run initialized this and
    // started forks. we can close them.
    shm_write,shmkey,"quit?",&([1]);
    for (i=1;i<=19;i++) sem_give,semkey,i,count=100;
    clean = 1;
  }
  
  if ( (wfs!=[]) && anyof(wfs.svipc>1) ) {
    shm_write,shmkey,"quit_wfs_forks?",&([1]);
    for (i=20;i<=nsem;i++) sem_give,semkey,i,count=100;
    wfs._svipc_init_done = 0;
    clean = 1;
  }
  if (clean) {
    usleep,500;
    status = svipc_clean();
    shm_init_done = 0;
    fork_done = 0;
  }
}

func quit_wfs_forks(void)
{
  extern wfs;
  extern prev_sync_counter;
  extern sync_init_done;
  //  extern nforks_per_wfs;

  if (!shm_init_done) return;
  
  shm_write,shmkey,"quit_wfs_forks?",&([1]);

  // FIXME: only set sem for wfs forks
  for (ns=1;ns<=nwfs;ns++) {
    if ((wfs(ns).type=="hartmann")&&(wfs(ns).svipc>1)) {
      sem_give,semkey,20+4*(ns-1),count=wfs(ns).svipc-1;
    }
  }
  
  wfs._svipc_init_done = 0;
  prev_sync_counter = [];
  //  nforks_per_wfs = [];
  usleep,200;
  //  sync_init_done=[];
}


func svipc_clean(void)
/* DOCUMENT svipc_clean(void)
   Clean and close the shared memory segment + semaphore allocation
   SEE ALSO:
 */
{
  shm_cleanup,shmkey;
  sem_cleanup,semkey;
}



func topwfs_svipc(void)
/* DOCUMENT topwfs_svipc(void)
   Routine called by *main* loop to get results (measurements)
   from parallelized WFS fork.
   SEE ALSO:
 */
{
  /* here:
     - sem_take sem_wfs_done
     - check that available data are indeed the one for iter we want
     - load back processed data in current session (mes and wfsimage)
     - load current phase, mircube, loopCounter, key into sharedmem
     - 1 sem_give sem_wfs_start to start wfs process
     - continue

     on the wfs side:
     - block on sem_wfs_start
     - once release, proceed to do calculations (mult_wfs)
     - store result data in sharemem, including mes, wfsimage
     - sem_give sem_wfs_done
   */

  // Block until slave done (skip for first iter)
  if (loopCounter>1) {
    if (smdebug) write,"topwfs_svipc: waiting for child to be done";
    sem_take,semkey,1;
    if (smdebug) write,"Child done, proceeding";
    // get data from shm
    svipc_mes = shm_read(shmkey,"svipc_mes");
    for (ns=1;ns<=nwfs;ns++) {
      wfs(ns)._fimage = &(shm_read(shmkey,swrite(format="wfs%d_image",ns)));
    }
  } else svipc_mes = array(0.0f,sum(wfs._nmes));
  // now prepare for next shot:
  if (smdebug) write,"topwfs_svipc: loading data in shm";
  shm_write,shmkey,"loop_counter",&[loopCounter];
  shm_write,shmkey,"mircube",&mircube;
  // etc?
  // give start signal to slave:
  if (smdebug) write,"topwfs_svipc: Giving trigger to Child";
  sem_give,semkey,0;

  return svipc_mes;
}



func psf_listen(void)
/* DOCUMENT psf_listen(void)
   The PSF forked process loop on this routine.
   It is triggered by a sem from the main process
   this routine execute the PSF calculation and
   then set a sem to notify the main process the PSF
   are ready.
   SEE ALSO:
 */
{
  extern mircube, loopCounter;
  extern im,imav;

  do {
    // wait for trigger from master
    if (smdebug) write,"psf_listen: waiting for trigger from parent";
    sem_take,semkey,3;
    if (smdebug) write,"psf_listen: got trigger";

    // check if we have to quit:
    if (shm_read(shmkey,"quit?")(1)) {
      write,format="%s child quitting\n",svipc_procname;
      yorick_quit;
    }

    // anything to sync?
    status = sync_child();
 
    // do the PSF calculations
    loopCounter = shm_read(shmkey,"loop_counter")(1);
    if (smdebug) write,"psf_listen: doing PSF calculations";

    mircube = shm_read(shmkey,"mircube");
 
    for (jl=1;jl<=target._nlambda;jl++) {
      for (jt=1;jt<=target._ntarget;jt++) {
        cubphase(,,jt)  = getPhase2dFromDms(jt,"target") +    \
          getPhase2dFromOptics(jt,"target") +                 \
          getTurbPhase(loopCounter,jt,"target");
      }
      // compute image cube from phase cube
      err = _calc_psf_fast(&pupil,&cubphase,&im,dimpow2,
                           target._ntarget,float(2*pi/(*target.lambda)(jl)));
       
      // Accumulate statistics:
      imav(,,,jl) = imav(,,,jl) + im;
    }
 
    // we're done.
    if (smdebug) write,"psf_listen: done, writing result in shm";
    shm_write,shmkey,"imsp",&im;
    shm_write,shmkey,"imlp",&imav;

    if (smdebug) write,"psf_listen: Giving trigger back to master";
    sem_give,semkey,4;

  } while (1);
  // if (smdebug) write,"psf_listen: calling oneself";
  // set_idler,psf_listen;
}


func topwfs_listen(void)
// for top WFS child
{
  extern mircube,loopCounter;

  do {
    // wait for trigger from master
    if (smdebug) write,"topwfs_listen: waiting for trigger from master";
    sem_take,semkey,0;
    if (smdebug) write,"topwfs_listen: got trigger";

    // check if we have to quit:
    if (shm_read(shmkey,"quit?")(1)) {
      write,format="%s child quitting\n",svipc_procname;
      yorick_quit;
    }
  
    // anything to sync?
    status = sync_child();
  
    // do the wfsing
    loopCounter = shm_read(shmkey,"loop_counter")(1);
    if (smdebug) write,"topwfs_listen: doing wfsing";

    mircube = shm_read(shmkey,"mircube");
    if ((sim.svipc>>2)&1) { // parallel WFSs
      mes = svipc_mult_wfs(loopCounter);
    } else { // single fork
      mes = mult_wfs(loopCounter);
    }
  
    // we're done.
    if (smdebug) write,"topwfs_listen: done, writing result in shm";
    // wfs._tt
    // wfs._lastvalidtt
    shm_write,shmkey,"svipc_mes",&mes;

    if (smdebug) write,"topwfs_listen: Giving trigger back to master";
    sem_give,semkey,1;

  } while (1);
  // if (smdebug) write,"topwfs_listen: calling oneself";
  // set_idler,topwfs_listen;
}



func wfs_listen(void)
// for WFSs children
{
  extern mircube;
  // wait for trigger from master
  if (smdebug) write,"wfs_listen: waiting for trigger from master";
  sem_take,semkey,sem4wfs(svipc_ns);
  if (smdebug) write,"wfs_listen: got trigger";

  // check if we have to quit:
  if (shm_read(shmkey,"quit?")(1)) {
    write,format="%s child quitting\n",svipc_procname;
    //    shm_free,shmkey,swrite(format="dm%d_shape",svipc_nm);
    yorick_quit;
    //    return;
  }
  
  // anything to sync?
  status = sync_child();
  
  // do something
  loopCounter = shm_read(shmkey,"loop_counter")(1);
  mircube = shm_read(shmkey,"mircube");
  if (smdebug) write,"doing wfsing of a single WFS";

  mes = svipc_single_wfs(loopCounter,svipc_ns);
  shm_write,shmkey,swrite(format="wfs%d_mes",svipc_ns),&mes;
  shm_write,shmkey,swrite(format="wfs%d_image",svipc_ns),wfs(svipc_ns)._fimage;
  
  if (smdebug) write,"wfs_listen: Giving trigger back to master";
  sem_give,semkey,sem4wfs(svipc_ns)+1;

  if (smdebug) write,"calling oneself";
  set_idler,wfs_listen;
}


func svipc_single_wfs(iter,ns)
{
  extern wfs;
  
  offsets = wfs(ns).gspos;
  phase   = getPhase2dFromDms(ns,"wfs");
  phase  += getPhase2dFromOptics(ns,"wfs");
  phase  += getTurbPhase(iter,ns,"wfs");

  if (wfs(ns).correctUpTT) {
    phase = correctUpLinkTT(phase,ns);
  }

  // get the measurements:
  if (wfs(ns).type == "hartmann" ) {smes = ShWfs(ipupil,phase,ns);}
  if (wfs(ns).type == "zernike")   {smes = ZernikeWfs(ipupil,phase,ns);}
  if (wfs(ns).type == "kl")        {smes = KLWfs(ipupil,phase,ns);}

  // subtract the reference vector for this sensor:
  if (wfs(ns)._cyclecounter == 1) {
    smes = smes - *wfs(ns)._refmes;
  }
    
  // compute the TT and subtract if required:
  // wfs._tt is computed inside this fork, so all good.
  wfs(ns)._tt(1) = sum( smes * (*wfs(ns)._tiprefvn) );
  wfs(ns)._tt(2) = sum( smes * (*wfs(ns)._tiltrefvn) );
  if (wfs(ns).filtertilt) {
    smes = smes - wfs(ns)._tt(1) * (*wfs(ns)._tiprefv) \
      - wfs(ns)._tt(2) * (*wfs(ns)._tiltrefv);
  }
  if (wfs(ns)._cyclecounter == 1) {
    wfs(ns)._lastvalidtt = wfs(ns)._tt;
  }
  return smes;
}



// func svipc_mult_wfs(iter,disp=)
// /* DOCUMENT func svipc_mult_wfs(iter,disp=)
//    svipc version of mult_wfs
//  */
// {
//   mes = [];

//   // parallelize main WFSs
//   for (ns=1;ns<=nwfs;ns++) {
//     // no data to send to slave (it will lookup mircube,
//     // which is the only data needed, with iter)
//     // send trigger
//     sem_give,semkey,sem4wfs(ns);
//   }

//   // wait for the children to be done:
//   for (ns=1;ns<=nwfs;ns++) sem_take,semkey,sem4wfs(ns)+1;

//   // all semaphores have been released, collect the pieces:
//   for (ns=1;ns<=nwfs;ns++) {
//     smes = shm_read(shmkey,swrite(format="wfs%d_mes",ns));
//     grow,mes,smes;
//   }
  
//   return mes;
// }


func cv2sv(charv)
/* DOCUMENT cv2sv(charv): char vector to string vector
   SEE ALSO: sv2cv
 */
{
  if ((allof(charv==[0x00]))||(charv==[])) return "";
  w = where(charv==0x00);
  i1 = _(1,w+1)(:-1);
  i2 = w;
  stringv = [];
  for (i=1;i<=numberof(w);i++) grow,stringv,string(&charv(i1(i):i2(i)));
  return stringv;
}

func sv2cv(stringv)
/* DOCUMENT sv2cv(charv): string vector to char vector
   SEE ALSO: cv2sv
 */
{
  if (stringv==[]) return char([0x00]);
  charv = [];
  for (i=1;i<=numberof(stringv);i++) grow,charv,*pointer(stringv(i));
  return charv;
}


// original_shm_read = shm_read;
// func shm_read(key,id,subscribe=)
// {
//   write,format="shm_read( %s )\n",id;
//   return original_shm_read(key,id,subscribe=subscribe);
// }

// original_shm_write = shm_write;
// func shm_write(key,id,a,publish=)
// {
//   write,format="shm_write, %s\n",id;
//   original_shm_write,key,id,a,publish=publish;
// }

// original_shm_free = shm_free;
// func shm_free(key,id)
// {
//   write,format="shm_free, %s\n",id;
//   original_shm_free,key,id;
// }
