/* yao_disp
 * Proof of concept to display yao data from an external
 * yao process dedicated to display
 * This should have 2 main advantages:
 * - it doesn't take resources from the main yao process
 * - it allows for a clean separation of both function, code-wise and
 *   performance wise. Possibly allows for easiest coding of 
 *   display-on-demand, i.e. the user selects what he/she wants to 
 *   display (e.g. psf, phase)
 * Right now, this si just a test that display wfs image and residual
 * phase for wfs #1 (only for shwfs). But it works. Some sync issues,
 * plus I need to code how to launch this (probably with spawn within
 * a yao session, not sure)
*/

require,"svipc.i";
require,"yao.i";

shmkey=0x0badcafe;
semkey=0x0badbeef;
shm_init,shmkey;
max_disp_freq = 1.;
status = create_yao_window();

shm_var,shmkey,swrite(format="wfs%d_fimage",1),ffimage;
shm_var,shmkey,swrite(format="wfs%d_phase",1),phase;
loopCounter = shm_read(shmkey,"loop_counter")(1);

func yao_disp(void)
{
  extern loopCounter;
  
  lc = shm_read(shmkey,"loop_counter")(1);
  if (loopCounter==lc) {
    // no display, just wait some small amount of time
    after,0.1/max_disp_freq,yao_disp;
    return;
  }
  loopCounter = lc;
  write,format="Loop Counter = %d\n",loopCounter;
  fma;
  plsys,1;
  pli,phase;
  plsys,2;
  pli,ffimage;
  after,1./max_disp_freq,yao_disp;
}
