/*
 * yaodh.i
 *
 * A collection of routines related to Disk Harmonic modes.
 *
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2013, Francois Rigaut
 * Original in MatLab: Norman Mark Milton (August 25, 2005)
 * Adapted to Yorick by Aurea Garcia Rissmann (May 23, 2010)
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

local yaodh;
/* dh: disk harmonic func Library

 dh_alt_elem       		- disk harmonic elements (alternate order)
 dh_alt_index      		- disk harmonic func index (alternate order)
 dh_alt_num        		- disk harmonic func numbers (alternate order)
 dh_bjprime_zero   		- disk harmonic BesselJ' zero
 dh_cos_coeff      		- disk harmonic func cos wave coefficient
 dh_dh             		- disk harmonic func evaluation
 dh_dhfast         		- disk harmonic func evaluation (fast)
 dh_dhindex        		- disk harmonic radial and azimuthal index
 dh_elem           		- disk harmonic elements
 dh_flip_x_coeff   		- disk harmonic func reflect (about y-axis) coefficient
 dh_flip_y_coeff   		- disk harmonic func reflect (about x-axis) coefficient
 dh_index          		- disk harmonic func index
 dh_norm           		- disk harmonic normalization constant
 dh_num            		- disk harmonic func numbers
 dh_rotate_coeff   		- disk harmonic func rotated coefficient
 dh_scale_coeff    		- disk harmonic func rescaled coefficient
 dh_shift_coeff    		- disk harmonic func shifted coefficient
 dh_sin_coeff      		- disk harmonic func sin wave coefficient
 dh_var           	 	- disk harmonic variance
 load_dh_bjprime_zero_tab	- load init file
 prepdiskharmonic		- prepares the frame with pupil


 Original in MatLab: Norman Mark Milton (August 25, 2005)
 Adapted to Yorick by Aurea Garcia Rissmann (May 23, 2010)

*/

  require,"bessel.i";
  // require,"hdf5.i";

/* how to use it:
   Example:

   prepdiskharmonic,128,100;
   load_dh_bjprime_zero_tab;
   p = dh_dhindex(2,1); // DH n=2,m=1, n always >= m.
   image = dh_dh(p(1),p(2),zr,ztheta);
   window,1;
   pli,image;

*/

func make_diskharmonic(size,diameter,ndhmodes,xc=,yc=,disp=)
/* DOCUMENT:
 * make_diskharmonic(size,diameter,ndhmodes,xc=,yc=,disp=)
 * shortcut to prepdiskharmonic. return data cube.
 */
{
  prepdiskharmonic,size,diameter,xc,yc;
  load_dh_bjprime_zero_tab;
  max_order = zernumero(ndhmodes)(1)+1;
  ndh=0;
  //ntmodes = sum(indgen(max_order+1));
  if (sim.debug) write,format="%s","#DH:  ";

  for (i=0;i<=max_order;i++) {
    for (k=0;k<=i;k++) {
      ndh = ndh+1;
      if (ndh>ndhmodes) break;
      if (sim.debug) write,format="\rDoing DH# %d",ndh;
      if (ndh == 1) {
        dh_tab = array(float,[3,size,size,ndhmodes]);
      }
      p = dh_dhindex(i,k);
      dh_tab(,,ndh) = dh_dh(p(1),p(2),zr,ztheta);
      if (disp == 1) {fma; pli,dh_tab(,,ndh);}
    }
  }
  write,"";

  return dh_tab(,,1:ndhmodes);

}


//===========================================================================

func prepdiskharmonic(size,diameter,xc,yc)

/* DOCUMENT prepdiskharmonic(size,diameter,xc,yc)
 * Call this function to set up the geometry for subsequent calls
 * to the diskharmonic functions. This is excatly the same function
 * as prepzernike; I just included it here with another name for
 * completeness of this disk harmonics file.
 * size : size of the 2d array on which future "diskharmonic" will be returned
 * diameter : diameter of the pupil in pixel in the array
 * xc, yc (optional) : Coordinates (in pixels of the center of the pupil)
 * Example:
 * > prepdiskharmonic,128,100
 * SEE ALSO:
 */
{
  extern zdim,zr,ztheta,zmask,zrmod,zmaskmod;

  if (xc == []) {xc = size/2+1;}
  if (yc == []) {yc = size/2+1;}

  radius= (diameter+1.)/2.;
  zdim	= size;
  zr	= dist(zdim,xc=xc,yc=yc)/radius;
  zmask	= (zr <= 1.);
  zmaskmod = (zr <= 1.2);
  zrmod	= zr*zmaskmod;
  zr	= zr*zmask;
  x	= float(span(1,zdim,zdim)(,-:1:zdim));
  y	= transpose(x);
  ztheta= atan(y-yc,x-xc);
}

//============================================================================

func dh_alt_elem(order)
{
/* DOCUMENT:
   dh_alt_elem: disk harmonic elements (alternate order)
   order - disk harmonic order number
   returns: number of disk harmonic basis funcs from
            order 0 through order o
*/
  if (order < 0) {
      error,"dh_alt_elem: invalid order number";
  } else {
      n = (2 * order^2) + order + 1;
  }
  return n;
  }

//======================================================================================

func dh_norm(dhn,dhm)
{
/* DOCUMENT:
   dh_norm: disk harmonic normalization constant
   dhn - Zero number
   dhm - Bessel order number
   returns: normalization constant for disk harmonic func
*/

  if (dhn < 0) error,"dh_norm: invalid zero number";
  if ((dhn == 0) && (dhm != 0)) error,"dh_norm: invalid Bessel number for order zero";
  if (dhn == 0){
      a = 1;
  } else {
      mabs = abs(dhm);
      l    = dh_bjprime_zero(dhn,dhm);
      k    = 2*pi *l;
      a    = sqrt(1/((1-(mabs/k)^2) * bessj(mabs,k)^2));
  }
  return a;
}

//========================================================================================

func dh_dh(dhn,dhm,r,theta)
{
/* DOCUMENT:
   dh_dh: disk harmonic func evaluation
   dhn   - Zero number
   dhm   - Bessel order number
   r     - radial coordinate
   theta - azimuthal angle coordinate

   returns: value of disk harmonic func at specified point
           (real part)
*/
  if (dhn < 0) error,"dh_dh: invalid zero number";

  d = dh_dhfast(dhm,dh_bjprime_zero(dhn,dhm),dh_norm(dhn,dhm),r,theta);
  return d;
}

//========================================================================================

func dh_bjprime_zero(dhn,dhm)
{
/* DOCUMENT:
   dh_bjprime_zero: disk harmonic BesselJ' zero
   dhn - Zero number
   dhm - Bessel order number
   returns: dhn'th zero to first derivative of dhm'th order BesselJ
          func (first kind)
*/
  extern dh_bjprime_zero_tab;

  if (dhn < 0) error,"dh_bjprime_zero: invalid zero number";
  mabs = abs(dhm);
  if (dhn == 0) {
    l = 0;
  } else {
    r = dimsof(dh_bjprime_zero_tab)(2);
    c = dimsof(dh_bjprime_zero_tab)(3);
    if (dhn>c) {
        error,"dh_bjprime_zero: invalid zero number";
    } else {
           if ((mabs+1) > r) {
               error,"dh_bjprime_zero: invalid order number";
           } else {
               l = dh_bjprime_zero_tab(mabs+1,dhn)
	   }
    }
  }
  return l;
}

//=====================================================================================

func dh_dhfast(dhm,l,a,r,theta)
{
/* DOCUMENT:
 dh_dhfast: disk harmonic func evaluation (fast)
 dhm   - Bessel order number
 l     - disk harmonic spatial frequency
 a     - disk harmonic normalization constant
 r     - radial coordinate
 theta - azimuthal angle coordinate
 returns: value of disk harmonic func at specified point
          (real part)
*/
   write,"teste2";
   mabs = abs(dhm);
   k    = 2*pi*l;
   d    = a*bessj(mabs,k*r);
   if (dhm > 0) {
    d = sqrt(2)*d*sin(mabs*theta);
   } else if (dhm < 0) d = sqrt(2)*d*cos(mabs*theta); // matrix elem-by-elem
                                                             // multiplication
   return d;
}

//========================================================================================

func dh_alt_index(o,n,m)
{
/*
   dh_alt_index: disk harmonic func index (alternate order)

   o - disk harmonic maximum order number
   n - Zero number
   m - Bessel order number

   returns: index of disk harmonic basis func with
           zero n and order m
*/
   if (o < 0) error,"dh_alt_index: invalid maximum order number";
   if (o < n) error,"dh_alt_index: invalid zero number";
   if (o < abs(m)) error,"dh_alt_index: invalid order number";
   if ((n == 0) && (m == 0)) {
       i = 1;
   } else {
       i    = ((2 * o) + 1) * (n - 1) + 2;
       mabs = abs(m);
       if (m != 0) {
          i = i + (2 * (mabs - 1));
          if (m > 0) {
              i = i + 1;
          } else {
              i = i + 2;
          }
       }
    }
   return i;
}

//========================================================================================

func dh_alt_num(o,i)
{
/* DOCUMENT:
  dh_alt_num: disk harmonic func numbers (alternate order)

  i - disk harmonic func index number

  returns: [n,m] func numbers disk harmonic basis func with
           index i with maximum order o
*/
   if (o < 0) error,"dh_alt_num: invalid order number";
   if (i < 1) error,"dh_alt_num: invalid func index";
   if (i == 1) {
       n = 0;
       m = 0;
   } else {
       i = i - 2;
       per_n = (2*o)+1;
       n     = floor(i/per_n) + 1;
       m     = i - ((n - 1) * per_n);
       if (m > 0) if (m%2==0) {  // arrumar esta funcao
                      m = -floor(m/2);
                  } else {
                      m = floor((m - 1)/2) + 1;
                  }
   }
   return [n,m];
}

//========================================================================================

func dh_cos_coeff(dhm, dhnb, dhmb, k)
{
/* DOCUMENT:
   dh_cos_coeff: disk harmonic func cos wave coefficient

   dhm  - Bessel order number
   dhnb - Zero number (basis)
   dhmb - Bessel order number (basis)
   k    - sin wave spatial frequency

   returns: value of disk harmonic coefficient for cos wave
*/
   if (dhm == dhmb) {
       if ((dhm <= 0) & (dhm%2 == 0)) {
           c = 2*cos(abs(dhm)*pi/2);
       } else {
           c = 0;
       }
       if (c != 0) {
           abnm = dh_norm(dhnb,dhm);
           kbnm = 2*pi*dh_bjprime_zero(dhnb,dhm);
           c    = c*abnm*bess_integral(dhm,k,kbnm);
           if (dhm != 0) {
               c = sqrt(2)*c;
           }
       }
   } else {
       c = 0;
   }
   return c;
}

//========================================================================================

func load_dh_bjprime_zero_tab(void)
{
/* DOCUMENT:
   dh_init: Initialize disk harmonic global variables
*/
extern dh_bjprime_zero_tab;
// require,"hdf5.i";

// dh_bjprime_zero_tab  = h5read("besseljprimezeros200.h5","/data");
dh_bjprime_zero_tab  = yao_fitsread(Y_SITE+"data/besseljprimezeros200.fits");

//clear besseljprimezeros200;
return;
}

//=====================================================================================

func dh_dhindex(n,k)
{
/* DOCUMENT:
   dh_dhindex: disk harmonic radial and azimuthal index

   n - disk harmonic order number
   k - disk harmonic mode number

   returns: return DH indices (dhn, dhm)
*/

/* assume n and k are valid
% if n < 0
%     error('dh_dhindex: invalid order number');
% end
% if (k < 0) | (k > n)
%     error('dh_dhindex: invalid mode number');
% end
*/

   dhm = n - (2 * k);
   if (dhm == 0) {
       dhn = floor(n/2);
   } else if (dhm > 0) {
             dhn = k + 1;
           } else {
             dhn = n - k + 1;
	   }
   return int([dhn,dhm]);
}

//======================================================================================

func dh_flip_x_coeff(n,m,nb,mb)
{
/*
  dh_flip_x_coeff: disk harmonic func reflect (about y-axis) coefficient

  n   - disk harmonic order number
  m   - disk harmonic azimuthal order number
  nb  - disk harmonic order  number (basis)
  mb  - disk harmonic azimuthal order number (basis)

  returns: value of disk harmonic coefficient with reflection (about y-axis)
*/
  c = dh_flip_y_coeff(n,m,nb,mb)*dh_rotate_coeff(n,m,nb,mb,pi);
  return c;
}

//======================================================================================

func dh_flip_y_coeff(n,m,nb,mb)
{
/*
  dh_flip_y_coeff: disk harmonic func reflect (about x-axis) coefficient

  n   - disk harmonic order number
  m   - disk harmonic azimuthal order number
  nb  - disk harmonic order  number (basis)
  mb  - disk harmonic azimuthal order number (basis)

  returns: value of disk harmonic coefficient with reflection (about x-axis)
*/
  if (n == nb) {
      if (m == mb) {
          if (mb > 0) {
              c = -1;
          } else { c = 1;}
      } else { c = 0;}
  } else { c = 0;}
  return c;
}

//=======================================================================================

func dh_rotate_coeff(dhn,dhm,dhnb,dhmb,theta)
{
/* DOCUMENT:
  dh_rotate_coeff: disk harmonic func rotated coefficient

  dhn   - Zero number
  dhm   - Bessel order number
  dhnb  - Zero number (basis)
  dhmb  - Bessel order number (basis)
  theta - rotation angle (radians)

  returns: value of disk harmonic coefficient with rotation
*/
  mabs = abs(dhmb);
  if (dhn == dhnb) {
    if (dhm == 0) {
        if (dhmb == 0) {
            c = 1;
        } else {
            c = 0;
        }
    } else if (dhm < 0) {
        // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
        if (dhm == dhmb) {
            c = cos(mabs*theta);
        } else if (dhm == -dhmb) {
            c = sin(mabs*theta);
        } else {
            c = 0;
        }
    } else {
        // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
        if (dhm == dhmb) {
            c = cos(mabs*theta);
        } else if (dhm == -dhmb) {
            c = -sin(mabs*theta);
        } else {
            c = 0;
        }
    }
  } else c = 0;
  return c;
}
//=============================================================================

func dh_scale_coeff(dhn,dhm,dhnb,dhmb,scale)
{
/* DOCUMENT:
  dh_scale_coeff: disk harmonic func rescaled coefficient

  dhn   - Zero number
  dhm   - Bessel order number
  dhnb  - Zero number (basis)
  dhmb  - Bessel order number (basis)
  scale - radial coordinate scale factor

  returns: value of disk harmonic coefficient with radial rescaling
*/
  if (dhm == dhmb) {
    anm  = dh_norm(dhn,dhm);
    abnm = dh_norm(dhnb,dhm);
    a    = 2*pi*dh_bjprime_zero(dhn,dhm)*scale;
    b    = 2*pi*dh_bjprime_zero(dhnb,dhm);
    c    = 2*abnm*anm*bess_integral(dhm,a,b);
  } else c = 0;
  return c;
}

//=================================================================================

func dh_shift_coeff(dhn,dhm,dhnb,dhmb,r0,theta0)
{
/*
  dh_shift_coeff: disk harmonic func shifted coefficient

  dhn    - Zero number
  dhm    - Bessel order number
  dhnb   - Zero number (basis)
  dhmb   - Bessel order number (basis)
  r0     - radial coordinate of shift (unit circle)
  theta0 - azimuthal coordinate of shift (rad)

  returns: value of disk harmonic coefficient with origin shift
*/
// dhnb = n, dhmb = m, dhn = n', dhm = m' where s + m' = m  or  s = m - m'

  mabs   = abs(dhm);
  mbabs  = abs(dhmb);
  sminus = mbabs - mabs;
  splus  = mbabs + mabs;
  if (dhm == 0) {
    if (dhmb == 0) d = dh_shift_term(mbabs,dhn,dhm,dhnb,dhmb,r0,theta0,0);
    if (dhmb < 0) d = sqrt(2) * dh_shift_term(mbabs,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
    if (dhmb > 0) d = sqrt(2) * dh_shift_term(mbabs,dhn,dhm,dhnb,dhmb,r0,theta0,1);
  }
  if (dhm < 0) {
    if (dhmb == 0) d = sqrt(2)*((-1)^mabs)* dh_shift_term(mabs,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
    if (dhmb < 0) {
        d1 = dh_shift_term(sminus,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
        d2 = dh_shift_term(splus,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
        d  = d1 + (((-1)^mabs)*d2);
    }
    if (dhmb > 0) {
        d1 = dh_shift_term(sminus,dhn,dhm,dhnb,dhmb,r0,theta0,1);
        d2 = dh_shift_term(splus,dhn,dhm,dhnb,dhmb,r0,theta0,1);
        d  = d1 + (((-1)^mabs)*d2);
    }
  }
  if (dhm > 0) {
    if (dhmb == 0) d = sqrt(2)*((-1)^mabs)* dh_shift_term(mabs,dhn,dhm,dhnb,dhmb,r0,theta0,1);
    if (dhmb < 0) {
        d1 = dh_shift_term(sminus,dhn,dhm,dhnb,dhmb,r0,theta0,1);
        d2 = dh_shift_term(splus,dhn,dhm,dhnb,dhmb,r0,theta0,1);
        d  = - d1 + (((-1)^mabs) * d2);
    }
    if (dhmb > 0) {
        d1 = dh_shift_term(sminus,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
        d2 = dh_shift_term(splus,dhn,dhm,dhnb,dhmb,r0,theta0,-1);
        d  = d1 - (((-1)^mabs)*d2);
    }
  }
  return d;
}

//=========================================================================================

func dh_shift_term(s,dhn,dhm,dhnb,dhmb,r0,theta0,trig)
{
   anm   = dh_norm(dhn,dhm);
   abnm  = dh_norm(dhnb,dhmb);
   knm   = 2*pi*dh_bjprime_zero(dhn, dhm);
   kbnm  = 2*pi*dh_bjprime_zero(dhnb, dhmb);
   d     = 2*abnm*anm*bess_integral(abs(dhmb),knm,kbnm) * bessj(s,knm*r0);
   if (trig > 0) {
    d = d*sin(s*theta0);
   } else if (trig < 0) d = d*cos(s*theta0);
   return d;
}

//=========================================================================================

func dh_sin_coeff(dhm,dhnb,dhmb,k)
{
/* DOCUMENT:
 dh_sin_coeff: disk harmonic func sin wave coefficient

 dhm  - Bessel order number
 dhnb - Zero number (basis)
 dhmb - Bessel order number (basis)
 k    - sin wave spatial frequency

 returns: value of disk harmonic coefficient for sin wave
*/
  if (dhm == dhmb) {
    if ((dhm < 0) & (dhm%2 != 0)) {
        c = 2*sin(abs(dhm)*pi/2);
    } else {
        c = 0;
    }
    if (c != 0) {
        abnm = dh_norm(dhnb,dhm);
        kbnm = 2*pi*dh_bjprime_zero(dhnb,dhm);
        c    = c*abnm*bess_integral(dhm,k,kbnm);
        if (dhm != 0) c = sqrt(2.)*c;
    }
  } else { c = 0;}
  return c;
}

//===============================================================================

func bess_integral(m,f1,f2)
{
/* DOCUMENT:
   calculates the Bessel integral over unit circle
*/
  if (f1==f2) {
     c = 0.5*(bessj(m,f1)^2-bessj(m-1,f1)*bessj(m+1,f1))
  } else {
     c = (f2*bessj(m-1,f2)*bessj(m,f1)-f1*bessj(m-1,f1)*bessj(m,f2))/(f1^2-f2^2);
  }
  return c;
}
//==============================================================================================

func dh_var(dhn,dhm,D,r0,oscale,vonKarman)
{
/* DOCUMENT:
  dh_var: disk harmonic variance

  dhm       - bessel order number
  dhn       - bessel zero number
  D         - aperture diameter (m)
  r0        - coherence length (m)
  oscale    - outer scale (m)
  vonKarman - use von Karman spectrum instead of Kolmogorov

 returns: variance of mode (dhn, dhm)
          for aperture D and coherence length r0
*/

if (dhn < 0) error,"dh_var: invalid zero number";

/*
% dh_alpha = 0.023 / (2.0^(5.0 / 3.0));
% dh_alpha = 0.023 / (2.0^(3.0 / 2.0));
% dh_alpha = 0.023 / pi;

% dh_alpha = 0.023 / (pi^4);
% dh_alpha = 2.8 * 10^(-4);
*/
  dh_alpha = 2.6*10^(-4);

  lnm      = dh_bjprime_zero(dhn,dhm);
  anm      = dh_norm(dhn,dhm);

  if (vonKarman == 1) {
    f0     = 1.0/(2.0*oscale);      // scaled to radius not diameter
    lscale = ((lnm^2 + f0^2)^(-11.0/6.0))/lnm;
  } else {
    lscale = lnm^(-14.0/3.0);
  }
  // v = dh_alpha / (pi^3) * ((D / r0)^(5.0 / 3.0)) * anm^2 * lscale;
  v = dh_alpha*((D/r0)^(5.0/3.0))*anm^2*lscale;

  return v;
}

func dh_dhfast(dhm,l,a,r,theta)
{
/* DOCUMEMNT:
  dh_dhfast: disk harmonic function evaluation (fast)

  dhm   - Bessel order number
  l     - disk harmonic spatial frequency
  a     - disk harmonic normalization constant
  r     - radial coordinate
  theta - azimuthal angle coordinate

  returns: value of disk harmonic function at specified point
           (real part)

*/
  mabs = abs(dhm);
  k    = 2*pi*l;
  d    = a * bessj(mabs,k*r);
  if (dhm > 0) {
    d = sqrt(2)*d*sin(mabs*theta);
  } else if (dhm < 0) d = sqrt(2)*d*cos(mabs*theta);
  return d;
}

func dh_bjprime_zero(dhn,dhm)
{
/* DOCUMENT:
  dh_bjprime_zero: disk harmonic BesselJ' zero

  dhn - Zero number
  dhm - Bessel order number

  returns: dhn'th zero to first derivative of dhm'th order BesselJ
           function (first kind)

*/
  extern dh_bjprime_zero_tab;
  if (dhn < 0) error,"dh_bjprime_zero: invalid zero number";
  mabs = abs(dhm);
  if (dhn == 0) {
    l = 0;
  } else {
    r = dimsof(dh_bjprime_zero_tab)(2);
    c = dimsof(dh_bjprime_zero_tab)(3);
    if (dhn > c) {
        error,"dh_bjprime_zero: invalid zero number";
    } else if (mabs+1 > r) {
              error,"dh_bjprime_zero: invalid order number";
           } else {
              l = dh_bjprime_zero_tab(mabs+1,dhn);
           }
  }
  return l;
}
