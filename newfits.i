/*
 * NEWFITS.I
 *
 * $Id: newfits.i,v 1.2 2008-11-19 00:53:19 frigaut Exp $
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
 * Yorick procedures to read a FITS file (including extensions) and
 * mechanism to parse header cards into keywords, value and comments
 * fields.  The goal of this package is only to provide these basic
 * functionalities.  No elaborate mechanisms is provided to work with
 * the header keywords, such as handling of world coordinate system.
 *
 * A good set of documents on fits standards can be found at:
 * http://fits.gsfc.nasa.gov
 *
 * Fits Standard:
 * ftp://nssdcftp.gsfc.nasa.gov/standard/fits/fits_standard.pdf
 *
 * newfits.i,v 0.1 2002/09/27 rigaut
 *
 * Routines included:
 * fitsRead(name,&phdr,&ehdr,extension=,onlyheader=) Reads a fits file
 * fitsHdrValue(hdr,keyword,&position,default=) Get a keyword value
 * fitsHdrComment(hdr) Return string vector with all comment lines
 * fitsHdrHistory(hdr)Return string vector with all history lines
 *
 * fitsWrite(name,data,header,exttype=,append=,rescale=) still in development
 * fitsBuildCard(keyword,value,comment,type=,form=) built a Card line
 * fitsAddCard(hdr,keyword,value,comment,type=,form=,after=,before=,first=,last=)
 * fitsDeleteCard(hdr,keyword) delete an entry in a header
 * fitsAddComment(hdr,comment,after=,before=,first=,last=)
 * fitsAddHistory(hdr,comment,after=,before=,first=,last=)
 * 
 * _fGetKeyword(line) Return keyword associated with the input line
 * _fGetValue(line,&comment) Return value associated with input line
 * _fCheckPrimaryHeader(hdr) Checks header is conform to fits standards
 * _fCheckExtensionHeader(hdr,&XTtype) Checks header is conform to Xtension std
 * _fGetNaxis(hdr) determine # of axis and dimensions
 * _fBuildArray(hdr) return array confirming info of TYPE and NAXISn
 * _fReadHeader(file,&address) Directly access and reads header
 * exitInError(message) Close currently open file and exits in error
 * _fBuildCoreHeader(exttype,data) Build the bear bone header for a given Xtension
 * _fStripCoreFromHeader(hdr) Strip a header from the core keywords
 * is_set(arg) test if the keyword "arg" is set (=1)
 * _fRescale(data, bitpix, &bscale, &bzero, data_min=, data_max=) Rescale array
 * 
 * Known bugs:
 * * fitsHdrValue does not parse properly the card image string values which
 *   contain a "/" before the comment.
 * Restrictions:
 * * Does not handle binary table extension yet. Only Image extensions.
 *
 * 2004apr30
 *   - fixed format bug in func fitsBuildCard (%20#g -> %#20g)
 *
*/

require,"string.i";

func fitsHead(name, &phdr, &ehdr, extension=, hdu=)
{
  fitsRead,name,phdr,ehdr,extension=,hdu=,onlyheader=1;
  return phdr;
}

func fitsRead(name, &phdr, &ehdr, extension= , hdu=, onlyheader= )
/* DOCUMENT data = fitsRead(filename, PrimaryHDR, ExtensionHDR,
                                          extension=, onlyheader=)
     Returns the Primary Data Array or an Extension (one based)
     Data Array of the FITS file FILENAME.
     Mandatory argument:
       - filename     : Input, Name of the fits file
     Optional arguments:
       - PrimaryHDR   : Output, Primary Header
       - ExtensionHDR : Output, Requested extension Header
     Optional keywords:
       - extension or hdu: Input, extension number (1 based, i.e. 1rst ext=1)
       - onlyheader   : Input, if set, data is not read, only header.
     The header is returned as a string array. Values can be parsed using the
     routine fitsHdrValue.
  EXAMPLES:
  Read an image from the Primary Data Array:
  im = fitsRead("hrwfs.1.fits")
  Read the Primary Header and Primary Data Array:
  im = fitsRead("hrwfs.1.fits",hdr);
  Read only the Primary Header:
  im = fitsRead("hrwfs.1.fits",hdr,onlyheader=1);
  Read the Data in extension#2, primary and extension(#1) header:
  im = fitsRead("~/N20020926S0050.fits",hdr,ehdr,extension=1);
  SEE ALSO: fitsHdrValue.
*/
{
  if (!is_void(hdu)) extension=hdu;
  if (is_void(extension)) extension=0;
  /* (extension number starts at 1) */
  
  file= open(name, "rb");
  sun_primitives, file;
  address= 0L;

  /* First read the primary header: */
  hdr= _fReadHeader(file,address);
  phdr= hdr;

  /* Some error checking */
  if (!_fCheckPrimaryHeader(hdr))
    exitInError,"Some standard fits keywords are missing in the primary header";
  if ((extension >= 1) && (!fitsHdrValue(hdr,"EXTEND",default=1))) {
    write,"WARNING: This FITS file is not supposed to contain fits Extensions";
    write,"(Keyword EXTEND not present or set to 'F'). I will try anyway...";
  }

  extn= 1;
  while (extn <= extension) {
    /* We have to deal with a fits extension. There is no point in
       reading the PDA (or previous extension data array), so we
       skip it and find the start address of the next Extension header. */
    naxis= _fGetNaxis(hdr);
    dim= (naxis(1) > 0);  /* In case NAXIS=0 */
    dim= dim*abs(fitsHdrValue(hdr,"BITPIX"))/8; /* # of byte in type */
    for (i=1;i<=naxis(1);i++) dim= dim*naxis(i+1);
    skip= long(ceil(dim/2880.)*2880);
    address+= skip;
    /* Now read the next extension header: */
    hdr= _fReadHeader(file,address);
    if (!_fCheckExtensionHeader(hdr,XTtype))
      exitInError,"Some standard fits keywords are missing in the Extension header"+extn;
    ehdr= hdr;
    extn++;
  }

  /* return without data if this keyword has been set */
  if (onlyheader) return;
  /* That's also a legal value and there might be no data */
  if (fitsHdrValue(hdr,"NAXIS") == 0) {return 0.;}

  /* Build data array (including determination of array type) */
  //  data= _fBuildArray(hdr);
  // The above call was replaced by putting the function directly inside
  // the main program as it is slitghtly faster.

  naxis = _fGetNaxis(hdr);
  type  = fitsHdrValue(hdr,"BITPIX");
  /* Valid types: */
  x= [];
  if (type == 8)   {x = char();}   /* 8 pixel values are unsigned bytes */
  if (type == 16)  {x = short();}  /* 16 pixel values are signed 2-byte integers */
  if (type == 32)  {x = long();}   /* 32 pixel values are signed 4-byte integers */
  if (type == -32) {x = float();}  /* -32 pixel values are 4-byte floating points */
  if (type == -64) {x = double();} /* -64 pixel values are 8-byte floating points */
  if (is_void(x)) exitInError,"Unrecognized BITPIX value";
  
  data = array(x,naxis);
  
  if (_read(file,address,data) != numberof(data))
    exitInError,"EOF encountered before all data were read";

  /* Scaling of data */
  bzero= fitsHdrValue(hdr,"BZERO",default=0.);
  bscale= fitsHdrValue(hdr,"BSCALE",default=1.);
  if (bscale != 1.) {data *= bscale;}
  if (bzero  != 0.) {data += bzero;}
  /* We're done */
  close,file;
  return data;
}

func fitsHdrValue2(hdr,keyword,&position,default=)
/* DOCUMENT val = fitsHdrValue(hdr,keyword,&position,default=)
   Returns the value associated with a keyword in a fits header.
   Returns "Not a Keyword" if keyword is not found and no default
   was given.
   Mandatory parameters:
   hdr:         (input) header, as returned by fitsRead (string array)
   keyword:     (input) keyword name, string
   Optional parameters:
   position:    (output) position of keyword in header
   Optional keyword:
   default:     (input) default return value if keyword is not found
                in the header.
   SEE ALSO: fitsRead, _fGetKeyword, _fGetValue
   This was an attempt to find a faster procedure that the following
   function, but it's slightly slower (not by much though).
 */
{
  key= strtrim(strtoupper(keyword));
  w=where(strtrim(_fGetKeyword(hdr)) == key);
  if (is_void(hdr(w)) & !is_void(default)) return default;
  if (is_void(hdr(w)) & is_void(default)) return [];
  position = w(1);
  return _fGetValue(hdr(position));
}

func fitsHdrValue(hdr,keyword,&position,default=)
/* DOCUMENT val = fitsHdrValue(hdr,keyword,&position,default=)
   Returns the value associated with a keyword in a fits header.
   Returns "Not a Keyword" if keyword is not found and no default
   was given.
   Mandatory parameters:
   hdr:         (input) header, as returned by fitsRead (string array)
   keyword:     (input) keyword name, string
   Optional parameters:
   position:    (output) position of keyword in header
   Optional keyword:
   default:     (input) default return value if keyword is not found
                in the header.
   SEE ALSO: fitsRead, _fGetKeyword, _fGetValue
 */
{
  key= strtrim(strtoupper(keyword));
  position= -1;
  
  for (i=1;i<=numberof(hdr);i++) {
    if (strtrim(_fGetKeyword(hdr(i))) == key) {
      position = i;
      return _fGetValue(hdr(i));
    }
  }
  if (!is_void(default)) return default;
  return "Not a Keyword";
}

sxpar = fitsHdrValue;
    
func _fGetKeyword(line) {return strtoupper(strpart(line,1:8));}

func _fGetValue(line,&comment)
{

  if (strpart(line,9:9) != "=") {return [];}
  /* This is not a line that contains a valid keyword+value */
  
  vline   = strpart(line,10:80);
  vline   = strtok(vline,"/");
  val     = strtrim(vline(1));
  comment = strtrim(vline(2));

  value   = val;
  /* Case logical */
  if (value == "T" || value == "t") {return int(1);}
  if (value == "F" || value == "f") {return int(0);}

  /* Case string */
  if (strpart(val,1:1) == "'") {return strpart(val,2:-1);}
  
  x = long();
  if (sread(val,format="%i",x)) {value = x;}

  if (strmatch(val,".")) {
    x = double();
    if (sread(val,format="%g",x)) {value = x;}
  }
  
  return value;
}

func fitsHdrComment(hdr)
{
  /* Return a string array containing the comment cards in hdr */
  comment = [];
  for (i=1;i<=numberof(hdr);i++) {
    if (strtrim(_fGetKeyword(hdr(i))) == "COMMENT")
      grow,comment,strpart(hdr(i),10:80);
  }
  return comment;
}

func fitsHdrHistory(hdr)
{
  /* Return a string array containing the history cards in hdr */
  history = [];
  for (i=1;i<=numberof(hdr);i++) {
    if (strtrim(_fGetKeyword(hdr(i))) == "HISTORY")
      grow,history,strpart(hdr(i),10:80);
  }
  return history;
}

func _fCheckPrimaryHeader(hdr)
{
  /* Is the header conform to FITS standard ?
     Some mandatory keywords have to be there */
  return (strtrim(_fGetKeyword(hdr(1))) == "SIMPLE" && 
          strtrim(_fGetKeyword(hdr(2))) == "BITPIX" && 
          strtrim(_fGetKeyword(hdr(3))) == "NAXIS");
}

func _fCheckExtensionHeader(hdr,&XTtype)
{
  /* Is the header conform to FITS standard ?
     Some mandatory keywords have to be there */

  if (strtrim(_fGetKeyword(hdr(1))) != "XTENSION") return 0;

  XTtype = strtrim(fitsHdrValue(hdr,"XTENSION"));

  if (XTtype == "IMAGE") {
    /* fits image extension */
    return (strtrim(_fGetKeyword(hdr(1))) == "XTENSION" && 
            strtrim(_fGetKeyword(hdr(2))) == "BITPIX" && 
            strtrim(_fGetKeyword(hdr(3))) == "NAXIS" &&
            (fitsHdrValue(hdr,"PCOUNT") == 0) &&
            (fitsHdrValue(hdr,"GCOUNT") == 1));
  }

  if (XTtype == "TABLE") {
    /* fits table extension */
    return (strtrim(_fGetKeyword(hdr(1))) == "XTENSION" && 
            (fitsHdrValue(hdr,"BITPIX") == 8) &&
            (fitsHdrValue(hdr,"NAXIS") == 2) &&
            strtrim(_fGetKeyword(hdr(4))) == "NAXIS1" &&
            strtrim(_fGetKeyword(hdr(5))) == "NAXIS2" &&
            (fitsHdrValue(hdr,"PCOUNT") == 0) &&
            (fitsHdrValue(hdr,"GCOUNT") == 1) &&
            strtrim(_fGetKeyword(hdr(8))) == "TFIELDS");
  }

  if (XTtype == "BINTABLE") {
    /* fits binary table extension */
    exitInError,"Binary table not yet handle by fitsRead";
    return (strtrim(_fGetKeyword(hdr(1))) == "XTENSION" && 
            strtrim(_fGetKeyword(hdr(2))) == "BITPIX" && 
            strtrim(_fGetKeyword(hdr(3))) == "NAXIS" &&
            strtrim(_fGetKeyword(hdr(4))) == "NAXIS1" &&
            strtrim(_fGetKeyword(hdr(5))) == "NAXIS2" &&
            strtrim(_fGetKeyword(hdr(6))) == "PCOUNT" &&
            strtrim(_fGetKeyword(hdr(7))) == "GCOUNT" &&
            strtrim(_fGetKeyword(hdr(8))) == "TFIELDS");
  }

  exitInError,"XTENSION type not recognized";
}

func _fGetNaxis(hdr)
{
  /* determine # of axis and dimensions */
  nax     = fitsHdrValue(hdr,"NAXIS",position);
  naxis   = nax;
  
  for (i=1;i<=nax;i++) {
    grow,naxis,_fGetValue(hdr(position+i));
  }
  return naxis;
}
  
func _fBuildArray(hdr)
{
  naxis = _fGetNaxis(hdr);

  type  = fitsHdrValue(hdr,"BITPIX");
  /* Valid types: */
  x= [];
  if (type == 8)   {x = char();}   /* 8 pixel values are unsigned bytes */
  if (type == 16)  {x = short();}  /* 16 pixel values are signed 2-byte integers */
  if (type == 32)  {x = long();}   /* 32 pixel values are signed 4-byte integers */
  if (type == -32) {x = float();}  /* -32 pixel values are 4-byte floating points */
  if (type == -64) {x = double();} /* -64 pixel values are 8-byte floating points */
  if (is_void(x)) exitInError,"Unrecognized BITPIX value";
  
  data = array(x,naxis);
  return data;
}

func _fReadHeader(file,&address)
{
  /* Internal function to read a header (HDU or extension)
     the atomic fits block is 2880 bytes long. */
  hdr           = [];
  buffer        = array(char, 80, 36);
  buffer_size   = 80*36;
  end_not_found = 1;

  do {
    if (_read(file, address, buffer) != buffer_size)
      exitInError, "cannot read header";

    address  += buffer_size;

    for (i=1; i<=dimsof(buffer)(3) && end_not_found; i++) {
      line    = string(&buffer(,i));
      grow,hdr,line;
      keyword = strtoupper(strpart(line,1:3));

      if (keyword == "END") {
        end_not_found=0;
        break;
      }
    }
  } while (end_not_found);
  return hdr;
}

func exitInError(message)
{
  close,file;
  error,message;
}
/*----------------------------------------------------------*/

func fitsWrite(filename,data,header,exttype=,append=,rescale=)
/* DOCUMENT fitsWrite(filename,data,header,exttype=,append=,rescale=)
   Write a fits file. Primary HDU and possible append extensions.
   Current valid extension = only image.
   Usage:
   fitsWrite,filename,data,header,exttype=,append=,rescale=
     filename: input; filename as a string
     data: data. Can be void, scalar or multidimensionnal
     header: header. mandatory keywords are overwritten.
     exttype: Extension type. "PRIMARY" (default), "IMAGE"
     rescale: rescale data array
   SEE ALSO:
 */
{
  if (is_void(exttype)) exttype="PRIMARY";

  exttype= strtoupper(strtrim(exttype));

  if (exttype != "PRIMARY" && exttype != "IMAGE" && exttype != "BINTABLE")
    error,"Not a valid Extension type. Type help,fitsWrite for valid type";

  if (exttype != "PRIMARY") append=1;

  // Rescale the data if requested:
  bzero= 0.;
  bscale= 1.;
  if (is_set(rescale) & !is_void(data)) {
    if (structof(data) == float || structof(data) == double) {
      if (structof(data) == float) {
        bitpix= 16;     // float  -> short
      } else {
        bitpix= 32;     // double -> long
      }
      data= _fRescale(data, bitpix, bscale, bzero);
    } else {
      error,"Can only use rescale on float and double";
    }
  }

  // Build a core header from data type/size
  hdr = _fBuildCoreHeader(exttype,data);
  // append the newly determined bzero and bscale
  grow,hdr,fitsBuildCard("BZERO",bzero,);
  grow,hdr,fitsBuildCard("BSCALE",bscale,);
  // test provided header to append to core header:
  if ((!is_void(header)) && (typeof(header) != "string"))
    error,"Wrong header type.";
  // strip input header from core keywords:
  if (!is_void(header)) grow,hdr,_fStripCoreFromHeader(header);
  // Add the END keyword
  grow,hdr,_fEndCard();

  // Open the file for writes:
  local address;
  if (append) {
    stream = open(filename, "ab");
  } else {
    stream = open(filename, "wb");
  }
  sun_primitives, stream;
  address= 0L;

  // write header:
  for (i=1;i<=numberof(hdr);i++) {
    line = *pointer(hdr(i));
    _write, stream, address, line(1:80);
    address +=80;
  }
  // pad for address = % 2880
  _fPadFits,stream,address;

  // writes data:
  if (data != []) {
    _write,stream,address,data;
    address += long(abs(numberof(data)*fitsHdrValue(hdr,"BITPIX"))/8.);
    _fPadFits,stream,address;
  }
  
  close, stream;

  // Yorick uses a Contents Log file, which should be removed...
  if (open(filename+"L", "", 1)) {remove, filename+"L";}
}  

func _fPadFits(stream,&address)
{
/* DOCUMENT _fPadFits(stream,address)
   This function pads with blanks up to the next multiple of 2880 and
   update address.
   SEE ALSO:
 */

  toad = long(ceil(address/2880.)*2880);
  lstr = toad-address;
  if (lstr != 0) {
    pad  = string(&array('\x20',lstr));
    _write, stream, address, (*pointer(pad))(1:lstr);
  }
  address = toad;
}

func _fBuildCoreHeader(exttype,data)

/* DOCUMENT header = _fBuildCoreHeader(exttype,data)
   Build the bearbone header for a given type (primary,
   image xtension, table xtension.
   Valid type are "primary","image","bintable"
   I don't handle binary Extension yet
   SEE ALSO:
*/
{
  upxtype = strtoupper(exttype);
  
  bitpix= [];
  if (typeof(data) == "void") bitpix=8;
  if (typeof(data) == "char") bitpix=8;
  if (typeof(data) == "short") bitpix=16;
  if (typeof(data) == "long") bitpix=32;
  if (typeof(data) == "int") bitpix=32;
  if (typeof(data) == "float") bitpix=-32;
  if (typeof(data) == "double") bitpix=-64;
  if (is_void(bitpix)) error,"Unsupported data type";
  if (data != []) {dim= dimsof(data);} else {dim=[0];}
  hdr= [];
  
  if (upxtype == "PRIMARY") {
    grow,hdr,fitsBuildCard("SIMPLE",int(1),"file conforms to FITS standard");
  }

  if (upxtype == "IMAGE") {
    grow,hdr,fitsBuildCard("XTENSION","IMAGE","IMAGE extension");
  }
  
  if (upxtype == "BINTABLE") {
    grow,hdr,fitsBuildCard("XTENSION","BINTABLE","BINARY TABLE extension");
  }

  grow,hdr,fitsBuildCard("BITPIX",bitpix,"number of bits per data pixel");
  grow,hdr,fitsBuildCard("NAXIS",dim(1),"number of data axes");
  for (i=1;i<=dim(1);i++)
    grow,hdr,fitsBuildCard("NAXIS"+swrite(format="%i",i),dim(1+i),
                         "length of axis"+swrite(format="%i",i));

  if (upxtype == "PRIMARY") {
    grow,hdr,fitsBuildCard("EXTEND",int(1),"FITS dataset may contain extensions");
  }
  
  if (upxtype == "IMAGE") {
    grow,hdr,fitsBuildCard("PCOUNT",0l,"required keyword; must = 0");
    grow,hdr,fitsBuildCard("GCOUNT",1l,"required keyword; must = 1");
  }
  
  return hdr;
}

func _fStripCoreFromHeader(hdr)

/* DOCUMENT header = _fStripCoreFromHeader(hdr)
   Strip (existing) header from the core card (simply,naxis,...)
   to be able to combine it with the core generated with
   _fBuildCoreHeader.
   SEE ALSO:
*/
{
  mask = (_fGetKeyword(hdr) != "SIMPLE  ") &
    (strpart(_fGetKeyword(hdr),1:5) != "NAXIS") & 
    (_fGetKeyword(hdr) != "BITPIX  ") & 
    (_fGetKeyword(hdr) != "EXTEND  ") & 
    (_fGetKeyword(hdr) != "PCOUNT  ") & 
    (_fGetKeyword(hdr) != "GCOUNT  ") & 
    (_fGetKeyword(hdr) != "XTENSION") & 
    (_fGetKeyword(hdr) != "BZERO   ") & 
    (_fGetKeyword(hdr) != "BSCALE  ") & 
    (_fGetKeyword(hdr) != "END     ");

  return hdr(where(mask));
}

func fitsDeleteCard(&hdr,keyword)
/* DOCUMENT fitsDeleteCard(hdr,keyword)
     Delete the card label "keyword" in hdr
     Returns hdr
   SEE ALSO: fitsAddCard
 */
{
  mask = (strtrim(_fGetKeyword(hdr)) != strtoupper(keyword));
  hdr = hdr(where(mask));
  return hdr;
}

func fitsAddCard(&hdr,keyword,value,comment,type=,form=,
                 after=,before=,first=,last=)
/* DOCUMENT fitsAddCard(hdr,keyword,value,comment,type=,form=,
                 after=,before=,first=,last=)
     Add a new entry (card) in existing header "hrd".
     All keywords are optional.
     type (string): to force type of value (e.g. for logical)
     form (string): format, e.g "%13.8f"
     after (string): keyword before which to insert new card
       in header (if it doesn't exist, new card will be
       inserted at the end)
     before (string): see after
     first (long, boolean): if set, new card will be inserted
       as first card in header
     last (long, boolean): if set, new card will be inserted
       as last card in header
   SEE ALSO: fitsBuildCard
 */
{
  line = fitsBuildCard(keyword,value,comment,type=type,form=form);

  if (noneof([is_set(first),is_set(last),is_set(before),
              is_set(after)])) last=1;

  if (is_set(first)) { hdr = _(line,hdr); return hdr; }
  if (is_set(last)) { hdr = _(hdr,line); return hdr; }

  // after or before are set.
  // set key to which ever is set
  key = (is_set(before) ? before : after);
  // find position in current header
  fitsHdrValue,hdr,key,pos;
  // if requested after/before keyword does not exist, append
  // at the end and return:
  // note: it is not a problem to append at the end, possibly
  // after the "END" keyword, as this will be cleaned up by
  // when it's written up (fitsWrite).
  if (pos == -1) return _(hdr,line);
  
  ncard = numberof(hdr);
  if ((is_set(before)) & (pos == 1)) { hdr = _(line,hdr); return hdr; } 
  if ((is_set(after)) & (pos == ncard)) { hdr = _(hdr,line); return hdr; }
  if (is_set(before)) { hdr = _(hdr(1:pos-1),line,hdr(pos:)); return hdr; }
  if (is_set(after)) { hdr = _(hdr(1:pos),line,hdr(pos+1:)); return hdr; }
}

func fitsAddComment(&hdr,comment,after=,before=,first=,last=)
/* DOCUMENT fitsAddComment(hdr,comment,after=,before=,
                           first=,last=)
     Add a new comment entry in existing header "hrd".
     All keywords are optional.
     after (string): keyword before which to insert new card
       in header (if it doesn't exist, new card will be
       inserted at the end)
     before (string): see after
     first (long, boolean): if set, new card will be inserted
       as first card in header
     last (long, boolean): if set, new card will be inserted
       as last card in header
   SEE ALSO: fitsAddHistory, fitsAddCard
 */
{
  line = "COMMENT "+comment+string(&(array(' ',80)));
  line = strpart(line,1:80);
  
  if (noneof([is_set(first),is_set(last),is_set(before),
              is_set(after)])) last=1;

  if (is_set(first)) { hdr = _(line,hdr); return hdr; }
  if (is_set(last)) { hdr = _(hdr,line); return hdr; }

  // after or before are set.
  // set key to which ever is set
  key = (is_set(before) ? before : after);
  // find position in current header
  fitsHdrValue,hdr,key,pos;
  // if requested after/before keyword does not exist, append
  // at the end and return:
  // note: it is not a problem to append at the end, possibly
  // after the "END" keyword, as this will be cleaned up by
  // when it's written up (fitsWrite).
  if (pos == -1) { hdr = _(hdr,line);return hdr; }
  
  ncard = numberof(hdr);
  if ((is_set(before)) & (pos == 1)) { hdr = _(line,hdr) ; return hdr; }
  if ((is_set(after)) & (pos == ncard)) { hdr = _(hdr,line) ; return hdr; }
  if (is_set(before)) { hdr = _(hdr(1:pos-1),line,hdr(pos:)) ; return hdr; }
  if (is_set(after)) { hdr = _(hdr(1:pos),line,hdr(pos+1:)) ; return hdr; }
}

func fitsAddHistory(&hdr,comment,after=,before=,first=,last=)
/* DOCUMENT fitsAddHistory(hdr,comment,after=,before=,
                           first=,last=)
     Add a new history entry in existing header "hrd".
     All keywords are optional.
     after (string): keyword before which to insert new card
       in header (if it doesn't exist, new card will be
       inserted at the end)
     before (string): see after
     first (long, boolean): if set, new card will be inserted
       as first card in header
     last (long, boolean): if set, new card will be inserted
       as last card in header
   SEE ALSO: fitsAddHistory, fitsAddCard
 */
{
  line = "HISTORY "+comment+string(&(array(' ',80)));
  line = strpart(line,1:80);
  
  if (noneof([is_set(first),is_set(last),is_set(before),
              is_set(after)])) last=1;

  if (is_set(first)) { hdr = _(line,hdr); return hdr; }
  if (is_set(last)) { hdr = _(hdr,line); return hdr; }

  // after or before are set.
  // set key to which ever is set
  key = (is_set(before) ? before : after);
  // find position in current header
  fitsHdrValue,hdr,key,pos;
  // if requested after/before keyword does not exist, append
  // at the end and return:
  // note: it is not a problem to append at the end, possibly
  // after the "END" keyword, as this will be cleaned up by
  // when it's written up (fitsWrite).
  if (pos == -1) { hdr = _(hdr,line); return hdr; }
  
  ncard = numberof(hdr);
  if ((is_set(before)) & (pos == 1)) { hdr = _(line,hdr);  return hdr; }
  if ((is_set(after)) & (pos == ncard)) { hdr = _(hdr,line); return hdr; }
  if (is_set(before)) { hdr = _(hdr(1:pos-1),line,hdr(pos:)); return hdr; }
  if (is_set(after)) { hdr = _(hdr(1:pos),line,hdr(pos+1:)); return hdr; }
}

func fitsBuildCard(keyword,value,comment,type=,form=)
  /* type optional, to force type of value (e.g. for logical) */
  /* form= format, optional */
{
  if (is_void(comment)) comment="";
  pad= string(&array('\x20',80));
  line= strpart(swrite(format="%-8s",strtoupper(keyword)),1:8);
  if (is_void(value)) {
    line = strpart(line+pad,1:80);
    return line;
  }
  
  line+= "= ";
  if (!is_void(form)) {
    line+= swrite(format=form,value);
  } else {
    if (typeof(value) == "string")
      line+= swrite(format="%-20s",swrite(format="'%-8s'",value));
    if (typeof(value) == "int") {
      if (value == 0) line+= swrite(format="%20s","F");
      if (value == 1) line+= swrite(format="%20s","T");
    }
    if (typeof(value) == "long") line+= swrite(format="%20i",value);
    if (typeof(value) == "float" || typeof(value) == "double")
      line+= swrite(format="%#20g",value);
  }
  
  line+= " / "+comment;
  line= strpart(line+pad,1:80);
  
  return line;
}

func is_set(arg)
  /* DOCUMENT is_set(arg)
     Returns 0 if element is void or equal to zero, 1 otherwise
     F.Rigaut 2002/06/03
     SEE ALSO: is_void
  */
{
  if (is_void(arg) | (arg == 0)) {return 0;}
  else {return 1;}
}

func _fRescale(data, bitpix, &bscale, &bzero, data_min=, data_max=)
/* DOCUMENT rescaled_data= _fRescale(data, bitpix, bscale, bzero, data_min=, data_max=)

        Linearly rescale the values of input array DATA to fit into
        integers with BITPIX bits per value (i.e., `char', `short' or
        `long' for BITPIX being 8, 16 and 32 respectively).

        Arguments BSCALE and BZERO are optional and purely outputs passed
        by address.  Their value will be set so that:
                DATA(i) = BZERO + BSCALE * RESCALED_DATA(i)
        where the equality may not be exact due to rounding errors.  The
        difference is however the smallest possible, i.e., less than
        BSCALE / 2.

        Keywords DATA_MIN and DATA_MAX may be used to supply the maximum
        and minimum data values or to set brightness cutoffs.

        Originally from fitsRescale in the Yorick fits.i distrib.
   SEE ALSO:
*/
{
  if (bitpix == 8) {
    // assume ``char'' is 8 bits unsigned
    file_type= char;
    file_unsigned= 1N;
  } else if (bitpix == 16) {
    // assume ``short'' is 16 bits signed
    file_type= short;
    file_unsigned= 0N;
  } else if (bitpix == 32) {
    // assume ``long'' is 32 bits signed
    file_type= long;
    file_unsigned= 0N;
  } else {
    error, "bad BITPIX (should be 8, 16 or 32)";
  }

  data_min= double((is_void(data_min)? min(data) : data_min));
  data_max= double((is_void(data_max)? max(data) : data_max));
  if (data_max < data_min) error, "bad DATA_MAX and DATA_MIN";

  if (file_unsigned) {
    file_min= 0.;
    file_max= 2.^bitpix - 1.;
  } else {
    file_min= -2.^(bitpix - 1);
    file_max=  2.^(bitpix - 1) - 1.;
  }

  if (data_max == data_min) {
    bzero= double(data_min);
    return array(file_type(0), dimsof(data));
  }

  bscale= (data_max - data_min) / (file_max - file_min);
  bzero= data_min - bscale * file_min;
  return file_type(min(file_max, max(file_min, (data - bzero) / bscale)));
}

func _fEndCard(void)
{return "END                                                                              ";}

// lowercase equivalent:
fitshead       = fitsHead;      
fitsread       = fitsRead;     
fitshdrvalue   = fitsHdrValue;  
fitshdrcomment = fitsHdrComment;
fitshdrhistory = fitsHdrHistory;
fitswrite      = fitsWrite;     
fitsdeletecard = fitsDeleteCard;
fitsaddcard    = fitsAddCard;   
fitsaddcomment = fitsAddComment;
fitsaddhistory = fitsAddHistory;
fitsbuildcard  = fitsBuildCard; 
