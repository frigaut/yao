local yao_setnsync;
/* DOCUMENT yao_setnsync.i
   Implement setting of a selected number of important variables
   in an yao session, e.g., r0, wfs noise, wfs flux, etc...
   - Set the proper variables, and write in shared memory if need be.
   - Execute the needed action to effect the change (e.g. changing
     the Cn2 profile will need to re-run getTurbPhaseInit() )

   Implement syncing of the child, when in svipc mode.
   - read the shared memory structure
   - Execute the needed action to effect the change (e.g. changing
     the Cn2 profile will need to re-run getTurbPhaseInit() )
   SEE ALSO: yao, yao_svipc
*/


// UTILITY/GENERIC FUNCTIONS

func init_sync(void)
{
  for (i=1;i<=numberof(all_svipc_procname);i++) {
    shm_write,shmkey,"sync_"+all_svipc_procname(i),&sv2cv("");
  }
}


func sync_child(void)
{
  slotname = "sync_"+svipc_procname;
  mes = shm_read(shmkey,slotname);
  mes = cv2sv(mes)(1);  // one at a time. FIXME
  if (mes!="") {
    include,["status = "+mes+"()"],1;
    shm_free,shmkey,slotname;
    shm_write,shmkey,slotname,&sv2cv("");
  }
}

func broadcast_sync(targets,msg)
{
  for (i=1;i<=numberof(targets);i++) {
    slotname = "sync_"+targets(i);
    shm_free,shmkey,slotname;
    shm_write,shmkey,slotname,&sv2cv(msg);
  }
}

/*
struct setnsync { string prop; 

func set_generic(prop,value)
{
  extern atm,opt,sim,wfs,dm,mat,tel,target,gs,loop,phi;
  local var;
  
  wfs.noise = noise_flag;

  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3"];
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_noise";
  }
  // no action to take.
  
}
*/


// SETS & SYNCS FUNCTIONS

// for the WFS forks:
func sync_wfs_forks(pupsh)
/* DOCUMENT sync_wfs_forks(pupsh)
   Must be invoqued from main process to sync wfs forks after
   any modification to the wfs structure has been done.
   This will sync the fork on the next call to sh_wfs
   SEE ALSO:
 */
{
  extern sync_init_done;

  if (pupsh==[]) pupsh=ipupil;
  
  // save the wfs structures
  var = vsave(wfs);
  // write in shm
  // if doing a shm_free, we need to protect with sem
  // as the client could try to read in between the shm_free
  // and the shm_write. not done. bug still open.
  if (sync_init_done) shm_free,shmkey,"wfs_structs";
  shm_write,shmkey,"wfs_structs",&var;
  shm_write,shmkey,"pupsh",&ipupil;
  // broadcast message to children
  // increment shm variable each time a sync is done.
  for (ns=1;ns<=nwfs;ns++) {
    vname = swrite(format="sync_wfs%d_forks",ns);
    if ( (sync_init_done) && (wfs(ns)._svipc_init_done) ) {
      var = shm_read(shmkey,vname);
    } else {
      var = [0];
      sync_init_done = 1;
    }
    shm_write,shmkey,vname,&(var+1);
  }
}

func sync_wfs_from_master(ns,nf)
/* DOCUMENT sync_wfs_from_master(ns,nf)
   To be called from a wfs fork to sync itself with master process
   ns = wfs #
   nf = fork #
   SEE ALSO:
 */
{
  extern wfs;
  extern prev_sync_counter;
  extern pupsh;

  if (prev_sync_counter==[]) prev_sync_counter=0;
  
  vname = swrite(format="sync_wfs%d_forks",ns);
  sync_counter = shm_read(shmkey,vname)(1); 

  if (sync_counter == prev_sync_counter) return;
  
  pupsh = shm_read(shmkey,"pupsh");

  //  while (catch(-1)) {
  //    // protects the following shm_read in case
  //    // the read happens between the shm_free and shm_write
  //    // of the parent (see sync_wfs_forks() ) above.
  //    if (sim.debug>=1) write,"Catched error in shm_read(wfs_structs)";
  //    usleep,100;
  //  }

  // also to protect. 
  //  do {
  var = shm_read(shmkey,"wfs_structs");
  //    test = ((structof(var)==int)&&(var(1)==-1));
  //    if (test) write,"Catched error in shm_read(wfs_structs), var=-1";
  //    pause,10;
  //  } while (test);
  
  restore,openb(var);

  prev_sync_counter = sync_counter;
  
  if (sim.debug>1) \
    write,format="WFS sync'ed on WFS#%d fork#%d\n", ns, nf;
}




func set_noise(noise_flag)
{
  extern wfs;
  local var;
  
  wfs.noise = noise_flag;

  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3"];
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_noise";
  }
  // no action to take.
}

func sync_noise(void)
{
  extern wfs;

  if (sim.svipc) {
    var = shm_read(shmkey,"wfs_structs");
    restore,openb(var);
  }
  write,format="Noise sync'ed on child %s\n",svipc_procname;
  // no action to take.
}

func set_dr0(dr0)
{
  extern atm;
  local var;

  if (!is_scalar(dr0)) {
    write,"Usage: set_dr0, scalar value";
    return;
  }
  atm.dr0at05mic = dr0

  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3","PSFs"];
    // save the wfs structures
    var = vsave(atm);
    // write in shm
    shm_write,shmkey,"atm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_dr0";
  }
  // action to take.
  getTurbPhaseInit,skipReadPhaseScreens=1;
}

func sync_dr0(void)
{
  extern atm;

  if (sim.svipc) {
    var = shm_read(shmkey,"atm_structs");
    restore,openb(var);
  }
  write,format="D/r0 sync'ed on child %s\n",svipc_procname;
  // action to take.
  getTurbPhaseInit,skipReadPhaseScreens=1;  
}


func set_cn2(layerfrac)
{
  extern atm;
  local var;

  if (nallof(dimsof(layerfrac)==dimsof(*atm.layerfrac))) {
    write,format="Usage: set_cn2, layerfrac = %d element vector\n",\
      numberof(*atm.layerfrac);
    return;
  }

  layerfrac = layerfrac/sum(layerfrac);
  
  *atm.layerfrac = layerfrac;

  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3","PSFs"];
    // save the wfs structures
    var = vsave(atm);
    // write in shm
    shm_write,shmkey,"atm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_cn2";
  }
  // action to take.
  getTurbPhaseInit,skipReadPhaseScreens=1;
}

func sync_cn2(void)
{
  extern atm;

  if (sim.svipc) {
    var = shm_read(shmkey,"atm_structs");
    restore,openb(var);
  }
  write,format="Cn2 sync'ed on child %s\n",svipc_procname;
  // action to take.
  getTurbPhaseInit,skipReadPhaseScreens=1;  
}


func set_misreg(misreg,nm=)
{
  extern dm;
  local var;

  if (nm==[]) nm=1:3;
  dm(nm).misreg = misreg;

  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3"];
    if ((sim.svipc>>1)&1) grow,targets,["DM1","DM2"];
    // save the wfs structures
    var = vsave(dm);
    // write in shm
    shm_write,shmkey,"dm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_misreg";
  }
  // no action to take.
}

func sync_misreg(void)
{
  extern dm;

  if (sim.svipc) {
    var = shm_read(shmkey,"dm_structs");
    restore,openb(var);
  }
  write,format="Misreg sync'ed on child %s\n",svipc_procname;
  // no action to take.
}


func reset_strehl(void)
{
  extern imav, niterok,strehls_inter;

  if (sim.svipc) {
    targets = ["PSFs"];
    // nothing to pass
    // broadcast message to children
    broadcast_sync,targets,"sync_reset_strehl";
  }
  // no action to take.
  imav *= 0;
  niterok = 0;
  strehls_inter *= 0;
}

func sync_reset_strehl(void)
{
  extern imav, niterok;
  imav *= 0;
  niterok = 0;
  
}

func set_ngs_geometry(wfsxpos,wfsypos)
{
  extern wfs;

  ngs = numberof(wfsxpos);
  wfs(6:6+ngs-1).gspos(1,) = wfsxpos;
  wfs(6:6+ngs-1).gspos(2,) = wfsypos;
  
  if (sim.svipc) {
    targets = ["WFS","WFS1","WFS2","WFS3","PSFs"];
    // save the wfs structures
    var = vsave(atm);
    // write in shm
    shm_write,shmkey,"atm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_cn2";
  }
  getTurbPhaseInit,skipReadPhaseScreens=1;
}

func sync_ngs_geometry(void)
{
  extern wfs;

  if (sim.svipc) {
    var = shm_read(shmkey,"wfs_structs");
    restore,openb(var);
  }
  write,format="NGS geometry sync'ed on child %s\n",svipc_procname;
  // no action to take.
  getTurbPhaseInit,skipReadPhaseScreens=1;
}

