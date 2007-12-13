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
 * $Id: utils.c,v 1.2 2007-12-13 00:58:31 frigaut Exp $
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
 * Revision 1.2  2007-12-13 00:58:31  frigaut
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

float ran1()
{
  float norm;
  norm = 2147483647.f;

  return random()/norm;
}


void _eclat_float(float *ar, int nx, int ny)
{
  int i,j,k1,k2;
  float a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}



void _poidev(float *xmv, long n)
{ 
  float gammln(float xx); 
  /*  float ran1(long *idum);*/
  static float sq,alxm,g,oldm=(-1.0); 
  float xm,em,t,y,y1; 
  long i;

  for (i=0;i<n;i++) {
    xm = xmv[i];
    if (xm == 0.0f) continue;
    if (xm < 20.0) { /* Use direct method. */ 
      if (xm != oldm) { 
	oldm=xm; 
	g=exp(-xm);  /* If xm is new, compute the exponential. */
      } 
      em = -1; 
      t=1.0; 
      do { 
	++em; 
	t *= ran1(); 
      } while (t > g);
    } else {  /* Use rejection method. */ 
      if (xm != oldm) { 
	oldm=xm; 
	sq=sqrt(2.0*xm); 
	alxm=log(xm); 
	g=xm*alxm-gammln(xm+1.0); 
      } 
      do { 
	do { 
	  y=tan(3.141592654*ran1());
	  em=sq*y+xm; 
	} while (em < 0.0); 
	em=floor(em); 
	t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g); 
      } while (ran1() > t); 
    } 
    xmv[i] = em;
  }
} 


float gammln(float xx) 
{ 
  /* Returns the value ln[?(xx)] for xx>0. */
  double x,y,tmp,ser; 
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155, 
			0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j; 

  y=x=xx; 
  tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  ser=1.000000000190015; 
  for (j=0;j<=5;j++) ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x); 
} 


void _gaussdev(float *xmv, long n)
{ 
  /* Returns a normally distributed deviate with zero mean and unit variance, 
     using ran1() as the source of uniform deviates. */ 

  /*  float ran1(long *idum); */
  static int iset=0; 
  static float gset; 
  float fac,rsq,v1,v2; 
  long i;

  for (i=0;i<n;i++) {
    if (iset == 0) { 
      do { 
	v1=2.0*ran1()-1.0; 
	v2=2.0*ran1()-1.0; 
	rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0); 
      fac=sqrt(-2.0*log(rsq)/rsq); 
      gset=v1*fac; 
      iset=1; 
      xmv[i] = v2*fac; 
    } else { 
      iset=0; 
      xmv[i] = gset; 
    } 
  }
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
