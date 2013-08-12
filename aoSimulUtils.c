/*
 * aoSimulUtils.c
 *
 * C utility functions for yao
 *
 * This file contains a number or C routines used in conjunction with
 * yao.i, the AO simulation tool written in Yorick. This codes
 * efficiently, in C, the most critical elements of the AO loop (func-
 * tion aoloop in yao.i, for yao versions > 2.0.
 * 
 * This file is part of the yao package, an adaptive optics simulation tool.
 *
 * Copyright (c) 2002-2013, Francois Rigaut
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


/************************************************************************
 * Function int _get2dPhase                                             *
 * Computes the integrated phase along one given direction (angle)      *
 * All init data are pre-computed to accelerate execution time within   *
 * the time critical aoloop() in yao.i.                                 *
 * Last modified: Dec 12, 2003.                                         *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

int _get2dPhase(float *pscreens, /* dimension [psnx,psny,nscreens] */
                int psnx, 
                int psny, 
                int nscreens, 
                float *outphase, /* dimension [phnx,phny] */
                int phnx, 
                int phny, 
                int *ishifts,    /* array of X integer shifts dimension [phnx,nscreens] */
                float *xshifts,  /* array of X fractional shifts dimension [phnx,nscreens] */
                int *jshifts,    /* array of Y integer shifts dimension [phnx,nscreens] */
                float *yshifts)  /* array of Y fractional shifts dimension [phnx,nscreens] */
     /* ishifts[k,nscreens] and jshifts[k,nscreens] are the integer shifts for screen[k]
        xshifts[k,nscreens] and yshifts[k,nscreens] are the fractional shifts for screen[k],
     */
{
  int i,j,k,ips,jps;
  int firstel;
  float wx1,wx2,wy1,wy2;

  /* Loop on phase screens */
  for (k=0;k<nscreens;++k) {

    /* first indice for this screen */
    firstel = k*(psnx*psny);
    
    /* Loop on indices of output (integrated) phase */
    for (j=0;j<phny;++j) {
      for (i=0;i<phnx;++i) {

        ips = ishifts[i+k*phnx];
        jps = jshifts[j+k*phny];
        
        /* Computes the weights for the 4 surrounding pixels */
        wx1 = (1.0f - xshifts[i+k*phnx]);
        wx2 = (xshifts[i+k*phnx]);
        wy1 = (1.0f - yshifts[j+k*phny]);
        wy2 = (yshifts[j+k*phny]);
        
        
        /* Safety net: don't access elements outside of pscreens memory space */
        if ( (firstel+ips+1+(jps+1)*psnx) >= (psnx*psny*nscreens) ) {return (1);}
        
        /* Finaly, compute and integrate outphase */
        outphase[i+j*phnx] += ( wx1*wy1*pscreens[firstel+ips+jps*psnx]
                              + wx2*wy1*pscreens[firstel+ips+1+jps*psnx]
                              + wx1*wy2*pscreens[firstel+ips+(jps+1)*psnx]
                              + wx2*wy2*pscreens[firstel+ips+1+(jps+1)*psnx]);

      }
    }
  }
  return(0);
}


/************************************************************************
 * Function void _dmsum                                                 *
 * This routine simply loop on the number of actuator and computes the  *
 * DM global shape, given a serie of influence functions and a command  *
 * vector. It can also accomodate more than one DM.                     *
 * It saves a bunch of time in aoloop() for such a simple function.     *
 * Last modified: Dec 14, 2003.                                         *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

void _dmsum(float *def,     // pointer to dm influence functions 
            int  nxdef,     // X dim
            int  nydef,     // Y dim
            int  nzdef,     // Z (3rd) dim = # IFs
            float *coefs,   // command coefficients
            float *dmshape) // pointer to output phase.
{
  /* Declarations */
  int i, k;
  int n = (nxdef*nydef);
  float co;

  /* Zero out dmshape */
  for ( i=0 ; i<n ; i++) { dmshape[i] = 0.0f; }


  /* Loop over influence functions and commands */
  for ( k=0 ; k<nzdef ; k++ ) {
    co = coefs[k];
    for ( i=0 ; i<n ; i++) {
      dmshape[i] += co * def[i+k*n];
    }
  }
}

void _dmsum2(float *def,     // pointer to dm influence functions
             long  *inddef,   // vector of indices where def != 0
             long  ninddef,   // # of elements in inddef
             long  ndef,      // Z (3rd) dim = # IFs
             float *coefs,   // command coefficients
             float *dmshape, // pointer to output phase.
             long  ndmshape) // # of elements in dmshape
{
  /* Declarations */
  int i, k;
  float co;

  /* Zero out dmshape */
  for ( i=0 ; i<ndmshape ; i++) { dmshape[i] = 0.0f; }

  /* Loop over influence functions and commands */
  for ( k=0 ; k<ndef ; k++ ) {
    co = coefs[k];
    for ( i=0 ; i<ninddef ; i++) {
      dmshape[inddef[i]] += co * def[inddef[i]+k*ndmshape];
    }
  }
}

/************************************************************************
 * Function void _dmsumelt                                                 *
 * This routine simply loop on the number of actuator and computes the  *
 * DM global shape, given a serie of influence functions and a command  *
 * vector. It can also accomodate more than one DM.                     *
 * It saves a bunch of time in aoloop() for such a simple function.     *
 * Written 2004apr06                                                    *
 * Last modified 2004apr06
 * Author: F.Rigaut                                                     *
 ************************************************************************/

void _dmsumelt(float *def,     // pointer to dm influence functions 
               int  nxdef,     // X dim
               int  nydef,     // Y dim
               int  nzdef,     // Z (3rd) dim = # IFs
               int  *i1,       // i indice at which to put def # k in output phase
               int  *j1,       // j indice at which to put def # k in output phase
               float *coefs,   // command coefficients
               float *dmshape, // pointer to output phase.
               int  outnx,     // X dim of output array
               int  outny)     // Y dim of output array
{
  /* Declarations */
  int i, j, k, ioff;
  int n = (nxdef*nydef);
  float co;

  /* Zero out dmshape */
  for ( i=0 ; i<(outnx*outny) ; i++) { dmshape[i] = 0.0f; }

  /* Loop over influence functions and commands */
  for ( k=0 ; k<nzdef ; k++ ) {
    co = coefs[k];
    for ( i=0 ; i<nxdef ; i++) {
      // if index out of final image, continue.
      if ( ((i1[k]+i) < 0) | ((i1[k]+i) > (outnx-1)) ) continue;
      for ( j=0 ; j<nydef ; j++) {
        // if index out of final image, continue.
        if ( ((j1[k]+j) < 0) | ((j1[k]+j) > (outny-1)) ) continue;
        // now we map the def indices into the output array indices
        ioff = (i1[k]+i) + outnx*(j1[k]+j);
        dmshape[ioff] += co * def[i+nxdef*j+n*k];
      }
    }
  }
}

