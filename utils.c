/*
 * utils.c
 *
 * C utility functions for yao
 *
 * This file contains a number of utility functions, coded in C to gain
 * execution time. It addresses functionalities that are missing in
 * yorick, mostly concerning 2D image processing.
 * 
 * This file is part of the yao package, an adaptive optics
 * simulation tool.
 *
 * $Id: utils.c,v 1.5 2010/07/02 21:26:51 frigaut Exp $
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
 * $Log: utils.c,v $
 * Revision 1.5  2010/07/02 21:26:51  frigaut
 * - merged Aurea Garcia-Rissmann disk harmonic code
 * - implemented parallel extension (sim.svipc and wfs.svipc)
 * - a few bug fixes (and many more bug introduction with these major
 *   parallel changes (!). Fortunately, the svipc=0 behavior should be unchanged.
 *
 * Revision 1.4  2008/01/02 13:54:53  frigaut
 * - correct size for the graphic inserts (no black border)
 * - updated spec files
 *
 * Revision 1.3  2007/12/13 16:04:21  frigaut
 * - modification to broken Makefile
 * - reshuffling of plug_in statement
 *
 * Revision 1.2  2007/12/13 00:58:31  frigaut
 * added license and header
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "ydata.h"
#include "yapi.h"
#include "pstdlib.h"


/************************************************************************
 * noop. For testing and timing. with parameter passing                 *
 ************************************************************************/

int _mynoop2(float *in, int nx, int ny, float *out, int fx, int fy, int binfact)
{
  return(0);
}

void Y_usleep(int nArgs)
{
  long milliseconds = YGetInteger(sp-nArgs+1);
//  useconds_t us;
//  us = (useconds_t)(milliseconds*1000l);
  usleep(milliseconds*1000l);
}

int _cosf(float *x, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    x[i]=cos(x[i]);
  }
  return (0);
}

int _sinf(float *x, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    x[i]=sin(x[i]);
  }
  return (0);
}
