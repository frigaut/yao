/*
 * yao_util.i
 *
 * A collection of routines for general purpose for yao
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

require,"style.i";
require,"pkg_mngr.i"; // kinput
require,"util_fr.i"; // nprint, wheremin, wheremax, typeReturn, exist
require,"linalg.i";

func escapechar(s)
{
  s=streplace(s,strfind("_",s,n=20),"!_");
  s=streplace(s,strfind("^",s,n=20),"!^");
  return s;
}


//---------------------------------------------------------

func zernumero(zn)
/* DOCUMENT zernumero(zn)
 * Returns the radial degree and the azimuthal number of zernike
 * number zn, according to Noll numbering (Noll, JOSA, 1976)
 * SEE ALSO: prepzernike, zernike
 */
{
  j	= 0;
  for (n=0;n<=100;n++)
   {
    for (m=0;m<=n;m++)
     {
      if (even(n-m))
       {
        j	= j+1;
        if (j == zn) {return [n,m];}
        if (m != 0)
         {
          j	= j+1;
          if (j == zn) {return [n,m];}
         }
       }
     }
   }
}

//---------------------------------------------------------

func gamma(arg)
/* DOCUMENT gamma(arg)
 * Gamma function.
 * SEE ALSO: gamma.i in yorick/i/
 */
{
  if (arg == 0.) {
    return 1.;
  } else {
    return exp(ln_gamma(arg));
  }
}

//---------------------------------------------------------

func factoriel(arg)
/* DOCUMENT factoriel(arg)
 * Return factoriel of the argument
 * SEE ALSO:
 */
{
  if (arg == 0) {
    return 1.;
  } else {
    res = 1.;
    for (i=1;i<=arg;i++) res = res*i;
    return res;
   }
}

//---------------------------------------------------------

func zernike(zn)
/* DOCUMENT zernike(zn)
 * Returns the zernike number zn, defined on a 2D array as per
 * the prepzernike function.
 * These zernikes follow the Noll (JOSA, 1976) numbering and
 * definition (rms of 1 over the pupil)
 * Example:
 * > prepzernike,128,100
 * > pli,zernike(6)
 * SEE ALSO: prepzernike, zernumero
 */
{
  extern zdim,zr,zteta,zmask,zrmod,zmaskmod;

  z	= array(float,zdim,zdim);
  znm	= zernumero(zn) ; n=znm(1) ; m=znm(2);

  for (i=0;i<=(n-m)/2;i++) {
    z = z + (-1.)^i*zr^(n-2.*i)*factoriel(n-i)/
      (factoriel(i)*factoriel((n+m)/2-i)*factoriel((n-m)/2-i));
  }
  if (odd(zn)) {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*sin(m*zteta);
    }
  } else {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*cos(m*zteta);
    }
  }

  return z*zmask;
}

//---------------------------------------------------------

func zernike_ext(zn)
/* DOCUMENT zernike_ext(zn)
 * Same as the zernike function, except that the function is not masked
 * at R=1. This might be useful for some WFS codes where the derivative
 * of the wavefront is needed (and therefore a pixel outside of the
 * pupil is used to compute the derivatives).
 * These zernikes follow the Noll (JOSA, 1976) numbering and
 * definition (rms of 1 over the pupil)
 * Example:
 * > prepzernike,128,100
 * > pli,zernike(6)
 * SEE ALSO: zernike, prepzernike, zernumero
 */
{
  extern zdim,zr,zteta,zmask,zrmod,zmaskmod;

  z	= array(float,zdim,zdim);
  znm	= zernumero(zn) ; n=znm(1) ; m=znm(2);

  for (i=0;i<=(n-m)/2;i++) {
    z = z + (-1.)^i*zrmod^(n-2.*i)*factoriel(n-i)/
      (factoriel(i)*factoriel((n+m)/2-i)*factoriel((n-m)/2-i));
  }
  if (odd(zn)) {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*sin(m*zteta);
    }
  } else {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*cos(m*zteta);
    }
  }

  return z*zmaskmod;
}

//---------------------------------------------------------

func prepzernike(size,diameter,xc,yc)
/* DOCUMENT prepzernike(size,diameter,xc,yc)
 * Call this function to set up the geometry for subsequent calls
 * to the zernike function.
 * size : size of the 2d array on which future "zernike" will be returned
 * diameter : diameter of the pupil in pixel in the array
 * xc, yc (optional) : Coordinates (in pixels of the center of the pupil)
 * Example:
 * > prepzernike,128,100
 * > pli,zernike(6)
 * SEE ALSO: zernike,zernike_ext,zernumero
 */
{
  extern zdim,zr,zteta,zmask,zrmod,zmaskmod;

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
  zteta	= atan(y-yc,x-xc);
}

//---------------------------------------------------------

func jpg_write_color(im,filename,cmin=,cmax=,quality=,noflip=)

/* DOCUMENT jpg_write_color(im,filename,cmin=,cmax=,quality=,noflip=)
 * Wrapper for the jpg_write procedure.
 * Reads out the current palette and uses it as color table to write
 * the image "im" in a jpg file.
 * Flags and keywords as in jpeg_write.
 * SEE ALSO: jpg_write, jpg_read, jpg_info
 */
{
  dimx = dimsof(im)(2);
  dimy = dimsof(im)(3);
  palette,r,g,b,query=1;
  bytim = bytscl(im);
  cubim= array(char,[3,3,dimx,dimy]);
  cubim(1,,) = r(bytim+1);
  cubim(2,,) = g(bytim+1);
  cubim(3,,) = b(bytim+1);
  jpg_write,cubim,filename,cmin=cmin,cmax=cmax,quality=quality,noflip=noflip;
}

//---------------------------------------------------------

func log00(void){logxy,0,0;}
func log11(void){logxy,1,1;}
/* DOCUMENT log00() and log11()
 * Shortcuts for logxy,0,0 and logxy,1,1
 * SEE ALSO: logxy
 */

//---------------------------------------------------------

func tbget(fh,dataptr,keyword)
{
  nfields = numberof(dataptr);
  fields  = array(string,nfields);

  for (i=1;i<=nfields;i++) {
    fields(i) = fits_get(fh,swrite(format="TTYPE%d",i));
  }

  t = transpose(*dataptr(where(fields == keyword)(1)));
  if (numberof(t) == max(dimsof(t)(2:))) { t = t(*);}
  return t;
}

//---------------------------------------------------------

func mrot(ang)
/* DOCUMENT mrot(angle)
 * returns the matrix of rotation for a given angle.
 * It has to be used as follow:
 * If you want to rotate a vector of two coefficients xy=[x,y],
 * You should do rotated vector = mrot(+,)*xy(+,);
 * Angle is in degrees.
 * SEE ALSO:
 */
{
  dtor=pi/180.;
  return [[cos(ang*dtor),-sin(ang*dtor)],[sin(ang*dtor),cos(ang*dtor)]];
}

//---------------------------------------------------------

func clmfit(y,x,&a,function,&yfit)
/* DOCUMENT clmfit(y,x,&a,function,&yfit)
 * Useful wrapper for the lmfit procedure.
 * y = the data to fit vs x
 * a = the output coefficients (may have initial value on input)
 * function = a string containing the function definition where
 * x and a must be used as variable and coefficients name
 * e.g. "a(1)+a(2)*cos(x)"
 * yfit = optional output. Best fit.
 * SEE ALSO: lmfit
 */
{
  system,"rm /tmp/foo.i";

  f	= open("/tmp/foo.i","w");
  write,f,"func foo(x,a) {return "+function+";}";
  close,f;
  include,"/tmp/foo.i",10;
  require,"lmfit.i";
  r= lmfit(foo,x,a,y);

  yfit = foo(x,a);
  return a;
}

//---------------------------------------------------------

func bin2(image)
/* DOCUMENT bin2(image)
 * Returns the input image, binned by a factor of 2.
 * That is, a 512x512 image is transformed in a 256x256 image.
 * one output pixel is the average of the 4 corresponding input ones,
 * so that it conserves the total intensity.
 * SEE ALSO:
 */
{
  d = dimsof(image);
  if (d(1) != 2) {
    error,"Bin only accepts images";
  }
  if (((d(2) % 2) != 0) || ((d(3) % 2) != 0)) {
    error,"Bin only accepts dimensions with even # of pixels";
  }
  sim= image+roll(image,[-1,-1])+roll(image,[-1,0])+roll(image,[0,-1]);
  return sim(::2,::2);
}

//---------------------------------------------------------

func extractImage(image,dimx,dimy,method=)
/* DOCUMENT extractImage(image,dimx,dimy,method=)
 * Interactively extract a subimage from a larger one.
 * dimx[,dimy] are the dimensions of the subimage (not used for method 2)
 * method = 1 center is selected with mouse
 * method = 2 use mouse to click and drag to define subimage
 * NOTE: The image has to be displayed before calling the rountine.
 * This allows the user to arrange cuts and zoom level.
 * SEE ALSO:
 */
{
  if (!is_set(method)) {method = 1;}
  if (is_void(dimy)) {
    if (is_void(dimx)) {
      method = 2; // all is done with mouse
    } else {
      dimy = dimx;
      method = 1; // only center is selected with mouse
    }
  }
  s  = dimsof(image);
  if (exist(dimx) && dimx > s(2)) {write,"X dimension too large";}
  if (exist(dimy) && dimy > s(3)) {write,"Y dimension too large";}

  if (method == 1) {
    beg: co = mouse(1,0,"Click on center of image to extract");
    co1 = long(round(co(1:2)-[dimx/2.,dimy/2.]+1));
    co2 = co1 + [dimx,dimy]-1;
    xc  = [co1(1),co2(1)];
    yc  = [co1(2),co2(2)];
    if (anyof(xc <=0) || anyof(xc > s(2))) {
      write,"Out of bound in X";
      goto beg;
    }
    if (anyof(yc <=0) || anyof(yc > s(2))) {
      write,"Out of bound in Y";
      goto beg;
    }
    sim = image(xc(1):xc(2),yc(1):yc(2));
  }

  if (method == 2) {
    co = long(round(mouse(1,1,"Click and drag to select image to extract")+0.5));
    xc  = [co(1),co(3)];
    yc  = [co(2),co(4)];
    xc  = xc(sort(xc));
    yc  = yc(sort(yc));
    sim = image(xc(1):xc(2),yc(1):yc(2));
  }
  return sim;
}

//---------------------------------------------------------

func surface(image,shade=)
/* DOCUMENT surface(image,shade=)
 * Simple wrapper to get a simple mesh (or shaded) of the input array.
 * SEE ALSO:
 */
{
  if (!exist(image))
  {
    text = \
    "For the equivalent of the IDL 'surface' routine, you should use : \n\
     #include \"plwf.i\" \n\
     window,style=\"nobox.gs\";  // Possibly \n\
     orient3,-35*pi/180.,25*(pi/180.); // for instance \n\
     light3, diffuse=.5, specular=1., sdir=[1,.5,1];  // defines lighting\n\
     xy = indices(dimsof(image)); \n\
     plwf,image,xy(,,1),xy(,,2); \n\
     or \n\
     plwf,image,xy(,,1),xy(,,2),edges=0,shade=1; for 'shade_surf'\n\n\
     Or Just call this function \n\
     \"> surface,image,shade=0/1\" and this function will do the job";
    write,text;
  } else
  {
    require,"plwf.i";
    window,style="nobox.gs";
    orient3,-35*pi/180.,25*(pi/180.);
    light3, diffuse=.5, specular=1., sdir=[1,.5,1];
    xy = indices(dimsof(image));
    if (is_set(shade))
      {plwf,image,xy(,,1),xy(,,2),edges=0,shade=1;}
    else
      {plwf,image,xy(,,1),xy(,,2);}
  }
}

//---------------------------------------------------------

func __sinc(ar)
/* DOCUMENT sinc(ar)
 * Return the sinus cardinal of the input array
 * F.Rigaut, 2002/04/03
 * SEE ALSO: Eric Thiebault wrote a sinc which is probably better.
 */
{
  local ar;
  ar = double(ar);
  w  = where(abs(ar) < 1e-8);
  if (exist(w)) {ar(w) = 1e-8;}
  return sin(ar)/ar;
}
if (!is_func(sinc)) sinc=__sinc;

//---------------------------------------------------------

func window2(dummy)
/* DOCUMENT window2()
 * Create a window with style "myboxed.gs"
 * F.Rigaut 2002/04/03
 * SEE ALSO: window, window3
 */
{window,wait=1,style="/Users/frigaut/Yorick/Francois/myboxed.gs";}

//---------------------------------------------------------

func uint(arg)
/* DOCUMENT uint(arg)
 * Return the unsigned version of an integer of long argument.
 * Sorry, I neede this for some fits files. There is a uint type
 * in yorick-mb, beware.
 * F.Rigaut, 2001/10
 * SEE ALSO:
 */
{
  arg = long(arg);
  tmp = where(arg < 0);
  if (numberof(tmp) > 0) {arg(tmp) = arg(tmp) + 65536;}
  return arg;
}

//---------------------------------------------------------

func decimal_time(str,delim)
/* DOCUMENT decimal_time(str,delim)
 *  Returns the decimal time (in hours) from string like "20:33:12"
 *  or "21&32&01" or "06 55 32". You can specify a delimiter.
 *  F.Rigaut 2001/10
 *  OBSOLETE. THIS ROUTINE IS SUPERSEEDED BY "ParseTime".
 *  SEE ALSO:
 */
{
  local res;
  for (i=1;i<=numberof(str);i++) {
    v	= grow(strtok(str(i),delim)(1),strtok(strtok(str(i),delim)(2),delim));
    hh = mm = ss = 0.;
    sread,v,hh,mm,ss;
    grow,res,hh+mm/60.+ss/3600.;
  }
  return res;
}

//---------------------------------------------------------

func fftfit(yin,fraccut,nsig)
/* DOCUMENT fftfit(yin,fraccut,nsig)
 * routine of iterative fit by FT, discarding aberrant points
 * yin = input vector to fit
 * fraccut = cut in the fourier plane in fraction of cut-off frequency
 * nsig = number of sigma for rejection of aberrant points
 * F.Rigaut 2001/10
 * SEE ALSO:
 */
{
  n	= numberof(yin);
  v	= grow(yin,yin(::-1));
  np	= 1;
  iter	= 0;
  while ((np != 0) && (iter <= 20))
   {
    iter= iter+1;
    f	= fft(v,1);
    mask	= float(f)*0.;
    mask(1:long(n*fraccut)) = 1;
    mask	= mask(::-1);
    mask(1:long(n*fraccut+1)) = 1;
    vfit	= float(fft(f*mask,-1))/(2*n);
    // fma;plg,v;plg,vfit,type=2; pause,500;
    sig	= (v-vfit)(rms);
    w	= where(abs(v-vfit) > nsig*sig);
    np	= numberof(w); // print,np;
    if (is_array(w) == 0) {np=0;}
    if (np != 0) {v(w) = vfit(w);}
   }
  if (iter >= 20) {print,"Max number of iteration reach. Exiting.";}
  yout	= vfit(1:n);
  return yout;
}

//---------------------------------------------------------

func rfftconvol(image,kernel)
/* DOCUMENT rfftconvol(image,kernel)
 * Not specialy optimized fft convolution.
 * Inputs are two real positive arrays (images) of identical size.
 * Output is a real array. The output is Normalized in flux.
 * F.Rigaut, 2002/04/04
 * SEE ALSO: convol
 */
{
  sz = dimsof(image);
  if (anyof(dimsof(kernel)-sz)) {
    error,"Image and kernel should have the same size !";
  }
  return float(fft(fft(image,1)*fft(kernel,1),-1))/sz(2)^2.;
}

//---------------------------------------------------------

func convol2d(image,kernel)
/* DOCUMENT convol2d(image,kernel)
 * For small kernel. Not FT based. Really slow ! use only for small images.
 * you can use convVE if you are running yao on a mac (see yao_veclib/yao_fast.i)
 * F.Rigaut, 2001
 * NOTE: name changed as of yorick-1.6.01. used to be named convol
 * SEE ALSO: rfftconvol
 */
{
  local im;
  s	= dimsof(image);
  sk	= dimsof(kernel);
  im	= array(float,s(2)+sk(2),s(3)+sk(3));
  im(1:s(2),1:s(3)) = image;
  imac	= im*0.;
  mask	= im*0.;
  mask(1:s(2),1:s(3)) = 1.;
  maskac= mask*0.;
  for (i=1;i<=sk(2);i++) {
      for (j=1;j<=sk(3);j++) {
          imac	 = imac+roll(im,[i,j])*kernel(i,j);
          maskac = maskac+roll(mask,[i,j])*kernel(i,j);
        }
    }
  imac	= imac/clip(maskac,1e-4,);
  //  return maskac;
  //  return imac;
  return imac(sk(2)/2+2:sk(2)/2+s(2)+1,sk(3)/2+2:sk(3)/2+s(3)+1);
}

//---------------------------------------------------------
extern _nowtime;
 _nowtime = array(double,10);

//---------------------------------------------------------
func myxytitles(xtitle,ytitle,xyoff,font=,height=)
{
  l = viewport();
  xdim = l(2)-l(1);
  ydim = l(4)-l(3);
  if (is_void(xyoff)) {xyoff=[0.,0.];}
  // coordinate of X axis label
  xx = l(1)+xdim/2.;
  xy = l(3)-ydim/9.+xyoff(2);
  plt,xtitle,xx,xy,height=height,font=font,tosys=0,justify="CC";
  // coordinate of Y axis label
  yx = l(1)-xdim/6.+xyoff(1);
  yy = l(3)+ydim/2.;
  plt,ytitle,yx,yy,height=height,font=font,tosys=0,orient=1,justify="CC";
}
//---------------------------------------------------------
func mypltitle(title,xyoff,font=,height=)
{
  l = viewport();
  xdim = l(2)-l(1);
  ydim = l(4)-l(3);
  if (is_void(xyoff)) {xyoff=[0.,0.];}
  // coordinate of title
  xx = l(1)+xdim/2.+xyoff(1);
  xy = l(4)+ydim/15.+xyoff(2);
  plt,title,xx,xy,height=height,font=font,tosys=0,justify="CC";
}
//---------------------------------------------------------

func axisLegend(xtext,ytext,xyoff=,yxoff=)
/* DOCUMENT axisLegend(xtext,ytext,xyoff=,yxoff=)
 * plot the axis captions. works for myboxed.gs graphic style.
 * F.Rigaut, 2001/11/10.
 */
{
  if (!is_set(xyoff)) {xyoff=0.;}
  if (!is_set(yxoff)) {yxoff=0.;}
  plt,xtext,0.42,0.39+xyoff,justify="CC",font="helvetica";
  plt,ytext,0.12+yxoff,0.65,justify="CC",font="helvetica",orient=1;
}

//---------------------------------------------------------

func fftshift(image,xs,ys)
/* DOCUMENT fftshift(image,xs,ys)
 * Shift the input array by an arbitrary amount (xs,ys) in pixel units.
 * Of course xs and ys can be fractional. This rountine shift the
 * image by passing in the Fourier plane. All the usual restrictions
 * apply:
 *    - The image should be well sampled (Nyquist)
 *    - There should not be discontinuities at the edges, etc...
 * The input array can be 1D or 2D
 * SEE ALSO: fftrebin
 */
{
  s = dimsof(image);

  if (s(1) == 2)                                   // case 2D
  {
    xy     = indices(s);
    x      = xy(,,1)-s(2)/2.-1.;
    y      = xy(,,2)-s(3)/2.-1.;
    // normalisation factor:
    tilt   = 64.*0.098174773*(xs*x/s(2)+ys*y/s(3));
    fsh    = array(complex,s);
    fsh.re = roll(cos(tilt));
    fsh.im = roll(sin(tilt));

    sim    = float(fft(fft(image,-1)*fsh,1));
    sim    = sim/sum(sim)*sum(image);
  }
  if (s(1) == 1)                                   // case 1D
  {
    x      = indgen(s(2))-s(2)/2.-1.;
    tilt   = 64*0.09817*xs*x/s(2);
    fsh    = array(complex,s);
    fsh.re = roll(cos(tilt));
    fsh.im = roll(sin(tilt));

    sim    = float(fft(fft(image,-1)*fsh,1));
    sim    = sim/sum(sim)*sum(image);
  }

  return sim;
}
//---------------------------------------------------------

func fftrebin(image,nreb)
/* DOCUMENT fftrebin(image,nreb)
 * Returns "image" rebinned nreb times (nreb should be an integer,
 * power of 2, i.e. 2, 4, 8, ...) using a Fourier technique (basically,
 * extention of the support in the Fourier plane by zero values.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: fft, fft_setup, fft_inplace
 */
{
  dim	= (dimsof(image))(2);
  ndim	= nreb*dim;
  imout	= array(complex,ndim,ndim);

  imfft = fft(eclat(image),1);
  imfft.re = eclat(imfft.re);
  imfft.im = eclat(imfft.im);
  imout(1:dim,1:dim) = imfft;
  imout = roll(imout,[-dim/2,-dim/2]);
  return eclat(float(fft(imout,-1))/dim/dim/nreb^2);
}

//---------------------------------------------------------
func yao_wfs_rotate_shift(ar,rot,sh,integer=)
/* DOCUMENT yao_rotate_shift
   Rotate and shift an input array
   This is meant to be used to rotate and shift the phase at the input
   of the wfs routines.
   Uses FFT shift and FFT rotate
   Restiction: rotation has to be [-90,+90]
   shift has to be reasonable in pixels (a few)
   rot is rotation angle in degrees (CW)
   shift is 2 element vector [x,y] of shift in pixels
   when integer is set, input is expected to be integer array and is trated as such
   (output is also integer array)
   rot shift it/s
   0   0     103
   with bilinear
   1   0     58
   0   1     63
   1   1     61
   with fft*
   1   0     32
   0   1     44
   1   1     25
   */
{
  extern ywrs_dim,ywrs_disk,ywrs_mask,ywrs_sdisk;

  if (sum(ar)==0) return ar;
  if (allof((rot==0.)&(sh==[0,0]))) return ar; 
  // First let's deal with rotation:
  if ((abs(rot)>90)&&(use_fftrotate)) error,"WFS rotation has to be within [-90,90]";
  // this is especially to rotate eiuther a phase or a pupil, but likely something 
  // that is defined over a disk of diameter < array size (allegedly < array size/2)
  // let's build a apodising mask to avoid ringing:
  marg = 3;
  if (ywrs_dim != sim._size) {
    ywrs_dim = sim._size;
    ywrs_disk = dist(sim._size) < (sim._size/2-marg);
    ywrs_sdisk = sum(ywrs_disk)*1.0f;
    ywrs_mask = smooth(ywrs_disk,marg);
  }
  if (integer) {
    tmp = indgen(sim._size);
    ar = bilinear(ar,tmp+sh(1),tmp+sh(2),grid=1);
    ar = rotate2(ar,-rot); // rotate2() rot parameter is opposite fft_rotate()
    return ar;
  }
  // else we are dealing with phase, float arrays
  cst = sum(ar*ywrs_disk)/ywrs_sdisk;
  // prepare array to rotate/shift. The following should be approximately 
  // zero-centred and rolled off at the edge to avoid/minimise ringing
  ar2rs = (ar-cst)*ywrs_mask;
  // Now we can shift:
  if (nallof(sh==[0,0])) {
    if (use_fftshift) ar2rs = fftshift(ar2rs,-sh(1),-sh(2));
    else {
      tmp = indgen(sim._size);
      ar2rs = bilinear(ar2rs,tmp+sh(1),tmp+sh(2),grid=1);
    } 
  }
  // and rotate:
  if (rot!=0) {
    if (use_fftrotate) ar2rs = fftrotate(ar2rs,rot);
    else ar2rs = rotate2(ar2rs,-rot);
  }
  return ar2rs;
}


func fftrotate(im,angle,xc=,yc=)
/* DOCUMENT fftrotate(im,angle)
    im    : square image
    angle : rotation angle in degrees (clockwise)
    xc,yc : x,y positions of the rotation center

    high precision image rotation using fft
    no information loss if : image is shannon sampled
                             image has sufficiently large guard band
    using the algorithm from :
    "Fast Fourier method for the accurate rotation of sampled images"
    Kieran G. Larkin, Michael A. Oldfield, Hanno Klemm
    Optics Communications 139 (1997) 99-106

    routine by d. gratadour 31/05/2011
    SEE ALSO:
  */

{
   size = dimsof(im);
   if (size(2) != size(3)) error,"works only on square images";
   nx = size(2);


   if (xc == []) xc = (nx%2 == 0 ? nx/2+0.5 : nx/2+1);
   if (yc == []) yc = (nx%2 == 0 ? nx/2+0.5 : nx/2+1);

   theta = angle * pi / 180.;

   stepx = tan(theta/2);
   stepy = -1.*sin(theta);

   if (nx%2) {
     tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-0.5));
     tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-0.5));
   } else {
     tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-1));
     tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-1));
   }

   compltiltx=array(complex,nx);
   compltiltx.im=roll(tiltx);
   compltiltx = compltiltx(,-::nx-1);

   compltilty=array(complex,nx);
   compltilty.im=roll(tilty);
   compltilty = compltilty(-::nx-1,);

   col  = span(1,nx,nx)(,-:1:nx);
   lig  = transpose(col);

   tmpc=array(complex,nx,nx);

   tmpc = fft(exp(compltiltx*(lig-xc))*fft(im,[1,0]),[-1,0]);
   tmpc = fft(exp(compltilty*(col-yc))*fft(tmpc,[0,1]),[0,-1]);
   tmpc = fft(exp(compltiltx*(lig-xc))*fft(tmpc,[1,0]),[-1,0]);

   return tmpc.re/nx/nx/nx;
}


// December 2003:This function is superseeded by the clip function that call the
// clip C rountine. Much faster than this one. see yorickfr.i
// 2004mar03: I realized that the custom yorick version may not be available
// on all machines I have installed yorick on. So here I provide a substitute
// if clip is not defined.
func oldclip(arg,lt,ht)
/* DOCUMENT pli, clip(arg, mini, maxi);
 * Returns the argument, which has been "clipped" to mini
 * and maxi, i.e. in which all elements lower than "mini"
 * have been replaced by "mini" and all elements greater
 * than "maxi" by "maxi". Array is converted to float.
 * Either "mini" and "maxi" can be ommited, in which case
 * the corresponding mini or maxi is not clipped.
 * Equivalent to the IDL ">" and "<" operators.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO:
 */
{
  require,"utils.i";  // from yutils
  local imo;
  imo = double(arg);
  if (lt == []) lt=min(imo);
  if (ht == []) ht=max(imo);
  if (is_scalar(imo)) {
    if (max(imo) > ht) imo=ht;
    if (min(imo) < lt) imo=lt;
    return imo;
  }
  if (max(imo) > ht) imo(where(imo > ht)) = ht;
  if (min(imo) < lt) imo(where(imo < lt)) = lt;
  return imo;
}

if (clip == []) {clip = oldclip;}

//---------------------------------------------------------

func strInt(ivec,nchar)
/* DOCUMENT strInt(ivec,nchar);
 * Create a string array which elements are the string
 * equivalent of each elements of "ivec", with as many
 * heading "0" added to fill a string of length nchar.
 * example:
 * print,strInt(indgen(10:12),4)
 * ["0010","0011","0012"]
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: str routines in string.i
 */
{return strpart("00000000000000"+swrite(ivec,format="%i"),-(nchar-1):0);}

//---------------------------------------------------------

func medianCube(cube)
/* DOCUMENT medianCube(cube)
 * Returns a 2D array which elements are the median along the
 * 3rd dimension of the input variable "cube".
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: median.
 */
{return median(cube,3);}

//---------------------------------------------------------

// this function is superseeded by the call to the C function distraw
// see yorickfr.i
// 2004mar03: same remark as for clip (see above). Superseeded by the faster
// dist in yorickfr.i, but as this custom version is not always available,
// I provide this one as backup.
func olddist(dim,xc=,yc=)
/* DOCUMENT dist(size,xc=,yc=)
 * Returns an array which elements are the distance to (xc,yc). xc and
 * yc can be omitted, in which case they are defaulted to size/2+1.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: indices
 */
{
  dim	= long(dim);
  if (xc == []) xc = int(dim/2)+1;
  if (yc == []) yc = int(dim/2)+1;
  x	= float(span(1,dim,dim)(,-:1:dim));
  y	= transpose(x);
  d	= float(sqrt((x-xc)^2.+(y-yc)^2.));
  d	= clip(d,1e-5,);
  return d;
}

if (dist == []) {dist = olddist;}

//---------------------------------------------------------
// 2004mar03: same remark as for clip (see above). Superseeded by the faster
// eclat in yorickfr.i, but as this custom version is not always available,
// I provide this one as backup.
func oldeclat(image)
/* DOCUMENT eclat(image)
 * Equivalent, but slightly faster (?) than roll. Transpose the four main
 * quadrants of a 2D array. Mostly used for FFT applications.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: roll.
 */
{
  d	= dimsof(image);
  dx	= d(2);
  dy	= d(3);
  x1=1; x2=dx/2 ; x3=x2+1 ; x4=dx;
  y1=1; y2=dy/2 ; y3=y2+1 ; y4=dy;
  out	= image*0.;
  out(x1:x2,y1:y2) = image(x3:x4,y3:y4);
  out(x3:x4,y1:y2) = image(x1:x2,y3:y4);
  out(x1:x2,y3:y4) = image(x3:x4,y1:y2);
  out(x3:x4,y3:y4) = image(x1:x2,y1:y2);
  return out;
}

if (eclat == []) {eclat = oldeclat;}

//---------------------------------------------------------

func calcpsf(pupil,phase,init=)
/* DOCUMENT calcpsf(pupil,phase,init=)
 * Compute psfs from pupil and phase using FFT.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: fft, rfftconvol
 */
{
  extern calcPsfWorkSpace;

  if ((calcPsfWorkSpace = []) || (is_set(init)))
  {
    // write,"\nSetting up FFT workspace for calcpsf";
    workspace= fft_setup(dimsof(pupil),1);
  }

  dim	= (dimsof(pupil))(2);
  p	= array(complex,dim,dim);
  p.re	= pupil*cos(phase);
  p.im	= pupil*sin(phase);
  psf = eclat(abs(fft(p,1,setup=calcPsfWorkSpace)))^2.;
  return psf;
}

//---------------------------------------------------------

func apod(length,degree)
/* DOCUMENT apod(length,degree)
 * Returns apodization functions for 1D Fourier transforms.
 * degree = 1 to 3 : apodization functions for the FTS
 * degree = 4      : sinc
 * degree = 5      : Bartlett filter (cf NR p 547)
 * degree = 6      : Hann filter
 * degree = 7      : Welch filter
 * SEE ALSO:
 */
{
  if ((degree < 0) || (degree > 7))
    {
      print,"Apodization degree should be 0, 1, 2, 3, 4, 5, 6 or 7";
      return -1;
    }

  if ((degree >= 0) && (degree <= 3))
    {
      c = transpose([[1.,0.,0.,0.],  // matrix of the coefficients of
       [0.548,-0.0833,0.5353,0.],    // the apodization fonction.
       [0.26,-0.154838,0.894838,0.], // c(2,*) = vector of coef for weak apod
       [0.09,0.,0.5875,0.3223]]);     // c(4,*) = "      "    "   "  strong "

      u = 2.*(indgen(length)-1.)/length-1.;

      ffil=u*0.;

      for (i=1;i<=4;i++) {ffil = ffil + c(degree,i)*(1-u^2.)^i;}
    }

  if (degree == 4)
    {
      x    = 2.*pi*(((indgen(length)-1.)/(length-1.))-0.5);
      ffil = 1.4914*sinc(x);
    }

  n  = indgen(length)-1.;
  if (degree == 5) {ffil = 1. - abs((n-0.5*length)/(0.5*length));}
  if (degree == 6) {ffil = 0.5*(1 - cos(2*pi*n/length));}
  if (degree == 7) {ffil = 1. - ((n-0.5*length)/(0.5*length))^2.;}

  return ffil;
}

//---------------------------------------------------------

func psd(s, length, step=, filter=, samp=, db=,noplot=,overplot=,
         sqroot=,roddier=,xtitre=,ytitre=,silent=,color=,type=,hist=)
/* DOCUMENT psd(s, length, step=, filter=, samp=, db=,noplot=,overplot=,
 *    sqrt=,roddier=,xtitre=,ytitre=)
 * Procedure PSD : Compute the Power Spectral Density of a vector
 * s       = variable on which the PSD has to be computed
 * length  = length of the subsample for FFTs
 * step    = shift in pixels between subsamples
 * filter  = apodization function as in apod.pro (usually 6)
 * samp    = sampling time
 * db      = plots in dB :10*alog10(dsp)
 * noplot  = do not plot
 * overplot= over plot
 * sqroot  = returns the sqrt of the dsp
 * roddier = plots the psd roddier style
 * SEE ALSO:
 */
{
  extern psdnumberofoverplots;
  if (is_void(s)) {
    write,"psd,data,fftlength,step=,filter=,samp=,db=,noplot=,overplot=,sqrt=,roddier=,xtitre=,ytitre=,silent=";
    return;
  }
  if (is_void(step)) { step = length/2; }
  if (is_void(filter)) { filter = 0; }
  if (is_void(samp)) { samp = 1; }
  if (is_void(noplot)) { noplot = 0; }
  if (is_void(xtitre)) { xtitre = "Frequency"; }
  if (is_void(ytitre)) { ytitre = ""; }
  if (is_set(roddier)) { ytitre = ytitre+"[freq * PSD]";} else {ytitle = ytitre+"[PSD]"; }

  if (length > numberof(s)) {
    error,"length > number of element in vector!!!";
  }

  if (!is_set(silent)) {
    print,"Stdev of input PSD vector (in psd.i) : ",s(rms);
  }

  nb = long((numberof(s) - length)/step)+1;

  if (!is_set(silent)) {
    write,format="Averaging %2d sample of length %5d shifted by %5d\n", \
              nb,length,step;
  }

  dsp  = array(double,length);
  if (filter == -1) {
    fil = array(1.,length);
  } else {
    fil = apod(length,filter);
  }
  fn   = sum(abs(fft(fil,1))^2.)/length^2.;

  for (i=0;i<=nb-1;i++) {
    tmp  = s(i*step+1:i*step+length);
    tmp  = tmp-avg(tmp);
    dsp  = dsp + abs(fft(tmp*fil,1))^2./length/fn;
  }

  dsp = dsp/nb;
  dsp = dsp(1:length/2);

  f = (indgen(length/2)-1.)/samp/2./(length/2-1.);
  f = (indgen(length/2)-1.)/samp/2./(length/2);

  if (!is_set(silent)) {
    write,format="Freq. Max = %8.6e\n",max(f);
  }

  dsp = dsp*2.; // _x2 because negative part omitted
  dsp = dsp/length*(length/2./max(f)); // to get in unit^2/Hz

  if (!is_void(sqroot)) { dsp = sqrt(dsp); }

  if (db) dsp = 10.*log10(dsp);
  if (roddier) dsp = f*dsp;

  if (!noplot) {
    if (!overplot) {
      if (hist) plh,dsp(2:),f(2:),color=color,type=type;
      else plg,dsp(2:),f(2:),color=color,type=type;
      psdnumberofoverplots=0;
    }

    if (overplot) {
      psdnumberofoverplots++;
      cols = ["red","blue","green","magenta"];
      if (hist) {
        plh,dsp(2:),f(2:),color=cols(psdnumberofoverplots),color=color,type=type;
      } else {
        plg,dsp(2:),f(2:),color=cols(psdnumberofoverplots),color=color,type=type;
      }
    }
    myxytitles,xtitre,ytitre;

    limits;
    if (roddier) { logxy,1,0; }
    else if (db) { logxy,1,0; }
    else { logxy,1,1; }
  }

  return [f,dsp];
}

//---------------------------------------------------------

man=help;
hitReturn=typeReturn;
encircledEnergy=encircled_energy;
findFiles = findfiles;
