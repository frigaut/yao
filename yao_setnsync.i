/*
 * yao_setnsync.i
 *
 * Compilation of functions related to yao svipc process communication
 *
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2012, Francois Rigaut
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

local yao_setnsync;
/* DOCUMENT yao_setnsync.i
   Implement setting of a selected number of important variables
   in an yao session, e.g., r0, wfs noise, wfs flux, etc...
   - Set the proper variables, and write in shared memory if need be.
   - Execute the needed action to effect the change (e.g. changing
     the Cn2 profile will need to re-run get_turb_phase_init() )

   Implement syncing of the child, when in svipc mode.
   - read the shared memory structure
   - Execute the needed action to effect the change (e.g. changing
     the Cn2 profile will need to re-run get_turb_phase_init() )
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
    if (smdebug) write,format="%s executing %s\n",svipc_procname,mes;
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

  var = shm_read(shmkey,"wfs_structs");

  restore,openb(var);

  prev_sync_counter = sync_counter;
  
  if (sim.debug>1) \
    write,format="WFS sync'ed on WFS#%d fork#%d\n", ns, nf;
  if (sim.debug>20) \
    write,format="%d thmethod=%d  thres=%f\n",getpid(),
      wfs(ns).shthmethod, wfs(ns).shthreshold;
}




func set_noise(noise_flag)
{
  extern wfs;
  local var;
  
  wfs.noise = noise_flag;

  if (sim.svipc) {
    //~ targets = ["WFS","WFS1","WFS2","WFS3"];
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_noise";
  }
  if (anyof(wfs.svipc)) {
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    status = sync_wfs_forks();
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
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    targets = _(targets,"PSFs");
    // save the atm structures
    var = vsave(atm);
    // write in shm
    shm_write,shmkey,"atm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_dr0";
  }
  // action to take.
  get_turb_phase_init,skipReadPhaseScreens=0;
}

func sync_dr0(void)
{
  extern atm;

  if (sim.svipc) {
    var = shm_read(shmkey,"atm_structs");
    restore,openb(var);
  }
  write,format="D/r0 sync'ed on child %s (D/r0=%.2f)\n",\
    svipc_procname,atm.dr0at05mic;
  // action to take.
  get_turb_phase_init,skipReadPhaseScreens=0;
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
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    targets = _(targets,"PSFs");
    // save the wfs structures
    var = vsave(atm);
    // write in shm
    shm_write,shmkey,"atm_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_cn2";
  }
  // action to take.
  get_turb_phase_init,skipReadPhaseScreens=1;
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
  get_turb_phase_init,skipReadPhaseScreens=1;  
}


func set_misreg(misreg,nm=)
{
  extern dm;
  local var;

  if (nm==[]) nm=1:3;
  dm(nm).misreg = misreg;

  if (sim.svipc) {
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
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
  extern imav, niterok,strehls_inter,psf_child_started;

  if (sim.svipc) {
    targets = ["PSFs"];
    // nothing to pass
    // broadcast message to children
    broadcast_sync,targets,"sync_reset_strehl";
  }
  // no action to take.
  imav *= 0;
  niterok = 0;
  // shm_write,shmkey,"imlp",&imav;
  psf_child_started=0;
  // strehls_inter *= 0;
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
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    targets = _(targets,"PSFs");
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_ngs_geometry";
  }
  get_turb_phase_init,skipReadPhaseScreens=1;
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
  get_turb_phase_init,skipReadPhaseScreens=1;
}

func set_lgs_profile(amp,alt,lgs_focus_alt,ns)
{
  extern wfs;

  if (numberof(amp)!=numberof(alt)) {
    write,"numberof(amp)!=numberof(alt)";
    return;
  }
  
  for (i=1;i<=numberof(ns);i++) {
    if (amp!=[]) wfs(ns(i)).lgs_prof_amp = &float(amp);
    if (alt!=[]) wfs(ns(i)).lgs_prof_alt = &float(alt);
    if (lgs_focus_alt!=[]) wfs(ns(i)).lgs_focus_alt = lgs_focus_alt; \
    shwfs_comp_lgs_defocuses,ns(i);
  }
  
  if (sim.svipc) {
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"sync_lgs_profile";
  }
}

func sync_lgs_profile(void)
{
  extern wfs;
  
  if (sim.svipc) {
    var = shm_read(shmkey,"wfs_structs");
    restore,openb(var);
  }
  write,format="LGS profile sync'ed on child %s\n",svipc_procname;
  //~ for (i=1;i<=numberof(ns);i++) shwfs_comp_lgs_defocuses,ns(i);
}

// more generic functions:
func sync_wfs(void)
{
  if (sim.svipc) {
    targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
    // save the wfs structures
    var = vsave(wfs);
    // write in shm
    shm_free,shmkey,"wfs_structs";
    shm_write,shmkey,"wfs_structs",&var;
    // broadcast message to children
    broadcast_sync,targets,"child_sync_wfs";
  }
}
  
func child_sync_wfs(void)
{
  extern wfs;
  
  if (sim.svipc) {
    var = shm_read(shmkey,"wfs_structs");
    restore,openb(var);
  }
  write,format="WFS structure sync'ed on child %s\n",svipc_procname;
}


func sync_opt(void)
{
  if (opt==[]) return;
  if (!sim.svipc) return;

  targets = all_svipc_procname(where(strmatch(all_svipc_procname,"WFS")));
  grow,targets,all_svipc_procname(where(strmatch(all_svipc_procname,"PSFs")));
  // save the opt structures
  var = vsave(opt);
  // write in shm
  shm_free,shmkey,"opt_structs";
  shm_write,shmkey,"opt_structs",&var;      
  // broadcast message to children
  broadcast_sync,targets,"child_sync_opt";
}
  
func child_sync_opt(void)
{
  extern opt;

  if (opt==[]) return;
  if (!sim.svipc) return;
  
  var = shm_read(shmkey,"opt_structs");
  restore,openb(var);
  write,format="OPT structure sync'ed on child %s\n",svipc_procname;
}

