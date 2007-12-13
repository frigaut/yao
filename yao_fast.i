/*
 * yao_fast.i
 *
 * wrappers for the compiled functions in yao_fast.c
 *
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: yao_fast.i,v 1.2 2007-12-13 16:04:21 frigaut Exp $
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
 * $Log: yao_fast.i,v $
 * Revision 1.2  2007-12-13 16:04:21  frigaut
 * - modification to broken Makefile
 * - reshuffling of plug_in statement
 *
 * Revision 1.1.1.1  2007/12/12 23:29:11  frigaut
 * Initial Import - yorick-yao
 *
 *
 */


//write,"Yao FFTW version";

func calcPSFVE(pupil,phase,scale=)
/* DOCUMENT func calcPSFVE(pupil,phase,scale=)
   Similar to calcpsf, but way faster.
   This function calls the C routine _calcPSFVE that uses the vectorial
   library vDSP fft functions.
   Pupil and phase have to be float, this is insured in this wrapper
   routine. If you have any care for speed, I recommend that your
   input array are already floats. It takes time to cast from one
   type to another.
   Phase can be a data cube, in which case the return image is also
   a data cube of equal dimensions.
   Scale is a scaling factor on phase (the used phase is = to
   input phase * scaling factor)
   Warning: Works only for powers of 2 !
   SEE ALSO: _calcPSFVE
 */
{
  if (typeof(pupil) != "float") {pupil=float(pupil);}
  if (typeof(phase) != "float") {phase=float(phase);}
  if (is_set(scale)) {scale = float(scale);} else {scale=1.0f;}
  
  dims = dimsof(phase);
  if (dims(2) != dims(3)) { error,"X and Y dimension have to be the same"; }
  
  outimage = array(float,dims);
  n2       = int(log(dims(2))/log(2));
  if ((2^n2) != dims(2)) { error," Dimension has to be a power of 2"; }

  if (dims(1) == 3) {nplans = int(dims(4));} else {nplans = 1n;}
  
  err = _calcPSFVE(&pupil,&phase,&outimage,n2,nplans,scale);

  return outimage;
}
extern _calcPSFVE
/* PROTOTYPE
   int _calcPSFVE(pointer pupil, pointer phase, pointer image, int n2, int nplans, float scale)
*/

func fftw_wisdom(void)
/* DOCUMENT func fftw_wisdom(void)
   this function should be run at the start of each yorick session.
   It reads out the wisdom file, if any, or calls _init_fftw_plans
   to optimize the wisdow if it does not find the file.
   SEE ALSO: _init_fftw_plans
 */
{
  wisdom_file = Y_USER+"/fftw_wisdom.dat";
  if (open(wisdom_file,"r",1)) { //file exists
    _import_wisdom,wisdom_file;
    //    write,format="%s\n","fftw wisdow file loaded";
  } else { //file does not exist
    //    write,"I did not find a fftw_wisdom.dat file in ~/Yorick/";
    //    write,"When you have a minute (actually it takes several hours), run:\n";
    //    write,"init_fftw_wisdom; (default)\n";
    //    write,"or init_fftw_wisdom,nlimit; (run optimization up to n^nlimit)\n";
  }
}

func init_fftw_wisdom(nlimit)
/* DOCUMENT func init_fftw_wisdom(nlimit)
   this function should be run once on your hardware to optimize
   fftw and save the wisdom file
   nlimit: fft will be optimized for n = 2^[1,...,nlimit]
   SEE ALSO: _init_fftw_plans
 */
{
  if (nlimit == []) nlimit=11;
  
  wisdom_file = Y_USER+"/fftw_wisdom.dat";
  if (open(wisdom_file,"r",1)) { //file exists
    write,format="%s\n","fftw wisdow file already exist!";
    write,format="%s\n","If you wish to re-run init_fftw_wisdom,";
    write,format="%s\n","delete the existing fftw_wisdom.dat";
  } else { //file does not exist
    _init_fftw_plans,nlimit;
    _export_wisdom,wisdom_file;
  }
}

extern _init_fftw_plans
/* PROTOTYPE
   int _init_fftw_plans(long nlimit)
*/
extern _import_wisdom
/* PROTOTYPE
   int _import_wisdom(string wisdom_file)
*/
extern _export_wisdom
/* PROTOTYPE
   int _export_wisdom(string wisdom_file)
*/


func fftVE(realp,imagp,dir)
{
  if (typeof(realp) != "float") {realp=float(realp);}
  if (typeof(imagp) != "float") {imagp=float(imagp);}

  sub = am_subroutine();
  if (sub) {
    eq_nocopy,x,realp;
    eq_nocopy,y,imagp;
  } else {
    x = realp;
    y = imagp;
  }

  dims = dimsof(x);
  if (dims(2) != dims(3)) { error,"arrays should be square"; }

  n2       = int(log(dims(2))/log(2));
  if ((2^n2) != dims(2)) { error," Dimension has to be a power of 2"; }

  _fftVE,&x,&y,n2,dir;

  if (sub) return;
  
  return [x,y];
}

//==================================================================
extern _fftVE
/* PROTOTYPE
   int _fftVE(pointer realpart, pointer imagpart, int n2, int dir)
*/

// do not use _fftVE2, it was an experiment, but _fftVE is actually faster:
extern _fftVE2
/* PROTOTYPE
   int _fftVE2(pointer in, pointer out, int n, int dir)
*/


//==================================================================

extern _shwfs
/* PROTOTYPE
   int _shwfs(pointer pupil, pointer phase, float phasescale,
   pointer phaseoffset, int dimx, pointer istart, pointer jstart,
   int nx, int ny, int nsubs, int sdimpow2, long domask, pointer submask, pointer kernel,
   pointer kernels, pointer kerfftr, pointer kerffti, int initkernels,
   int kernelconv, pointer binindices, int binxy, pointer centroidw,
   pointer fimage, pointer imistart, pointer jmistart, int fimnx,
   int fimny, pointer flux, pointer rayleighflux, pointer skyflux,
   pointer threshold, pointer bias,
   pointer flat, float ron, float darkcurrent, int noise,
   int rayleighflag, pointer rayleigh, pointer bckgrdcalib, int bckgrdinit,
   int bckgrdsub, pointer mesvec, int counter, int niter)
*/

extern _shwfsSimple
/* PROTOTYPE
   int _shwfsSimple(pointer pupil, pointer phase, float phasescale,
   pointer phaseoffset, int dimx, int dimy, pointer istart, pointer jstart,
   int nx, int ny, int nsubs, float toarcsec, pointer mesvec)
*/

extern _cwfs
/* PROTOTYPE
   int _cwfs(pointer pupil, pointer phase, float phasescale, pointer phaseoffset,
   pointer cxdef,pointer sxdef, int dimpow2, pointer sind, pointer nsind,
   int nsubs, pointer fimage, pointer fimage2, float nphotons, float skynphotons,
   float ron, float darkcurrent, int noise, pointer mesvec)
*/

fftw_wisdom;
