/*
 * yao_utils.i
 *
 * wrapper for the compiled functions in yao_utils.c
 *
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2017, Francois Rigaut
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

require,"imutil.i";

extern usleep;

extern _mynoop2
/* PROTOTYPE
   int _mynoop2(pointer in, int nx, int ny, pointer out, int fx, int fy, int binfact)
*/

/* The following are wrapper to C routines that have no Yorick call functions.
   They are used in yao for specific operations */

extern _dmsum
/* PROTOTYPE
   void _dmsum(pointer def, int nxdef, int nydef, int nzdef, pointer coefs, pointer outphase)
*/

extern _dmsum2
/* PROTOTYPE
   void _dmsum2(pointer def, pointer inddef, long ninddef, long ndef, pointer coefs, pointer outphase, long ndmshape)
*/

extern _dmsumelt
/* PROTOTYPE
   void _dmsumelt(pointer def, int nxdef, int nydef, int nzdef, pointer i1, pointer j1,
   pointer coefs, pointer outphase, int outnx, int outny)
*/

extern _get2dPhase
/* PROTOTYPE
   int _get2dPhase(pointer pscreens, int psnx, int psny, int nscreens, pointer skip, pointer outphase, int phnx, int phny, pointer ishifts, pointer xshifts, pointer jshifts, pointer yshifts)
*/

func cosf(array)
/* DOCUMENT func cosf(array)
   Returns the cos of the argument.
   Input and output are float type.
   SEE ALSO:
 */
{
  if (typeof(array) != "float") {
    error,"This function is for float types!";
  }
  res = array;
  _cosf,&res,numberof(array);
  return res;
}
extern _cosf
/* PROTOTYPE
   int _cosf(pointer data, int size)
*/


func sinf(array)
/* DOCUMENT func sinf(array)
   Returns the sin of the argument.
   Input and output are float type.
   SEE ALSO:
 */
{
  if (typeof(array) != "float") {
    error,"This function is for float types!";
  }
  res = array;
  _sinf,&res,numberof(array);
  return res;
}
extern _sinf
/* PROTOTYPE
   int _sinf(pointer data, int size)
*/


// comment following line to have a deterministic random (!) start...
ran1init;  // init random function for poidev.
