
//----------------------------------------------------

/* DOCUMENT wfs, dm, atmospheric, etc.. structures:
   Main structures for the AO simul package parameters
   If additional parameters are needed, they should be entered in these
   structures definition and changes reflected in the parameter file
   (e.g. sh12.par)

   There are several type of entries:
   - long, float, string scalars: e.g.
   > atm.dr0at05mic = 33.;
   - pointers: These are pointers to variables that can have arbitrary
   number of elements. You have to define them in the following way:
   > wfs.nsubperring = &( [6,12,18] );

   Structure members can be accessed in the following way:
   var = dm.type;
   var = *atm.layerfrac;  // "*" dereference a pointer

   If there are several instance of a given structure: for instance, it is
   common to have a system with several dm:
   dm.type is thus a vector of the types of all instance of structure "dm"
   dm(1).type is a scalar
   *wfs(1).nsubperring is a vector.
   (*wfs(1).nsubperring)(1) is the first element of this vector.

   The variables with a "_" in frnt of them are internal variables. they are set
   and used while the program is running. You can still access them, and possibly
   modify them to reach a particular purpose.

   The following generic structures are instanced into structures of the same
   name, without the "_struct" appended, when the parameter file is read. For
   instance, "atm" will be the structure containing the atmospheric parameters
   as defined in the generic type atm_struct below.


   SYNTAX OF THE COMMENTS BELOW:
   For each entries, we give the type (scalar, vectorptr -vector pointer-, etc),
   what the parameter is, possible comments, whether the parameter is
   required or optional, and the default between [].
   As a general comment: when the structures are instanciated, all their elements
   get a default value. This is zero (0) for a float or long (scalar or vector),
   empty string for a string, and 0x0 for a pointer.

   modified 2004oct18 for aosimulv3.3 to 3.5
   modified 2004july for aosimulv3.0 to 3.3
   modified 2004jan-mar for aosimulv2.4 to 3.0
   modified 2003dec20-24 for aosimul-v2.3
   modified 2003dec19 for aosimul-v2.2
   modified 2003feb19 for aosimul-v1.2
   AUTHOR: F.Rigaut, beginning 2002
   SEE ALSO: All ao simul functions (aoread, aoinit, aoloop).
*/

struct sim_struct
{
  string  name;           // A name for this simulation run. Optional [none]
  long    pupildiam;      // Pupil diameter in pixels. Required [none]
  long    pupilapod;      // whether to use an apodized pupil (rolled off @ edges) for
                          // some image calculations. Optional [0]
  long    debug;          // set the debug level (0:no debug, 1:some, 2: most). Optional [0]
  long    verbose;        // set the verbose level (0:none, 1: some, 2: most). Optional [0]
  long    svipc;          // set to the number of process for parallelization
                          // 0 = no parallelization
                          // bits    effect
                          // 0(1)    WFS/DM global split (2 process)
                          // 1(2)    PSFs calculation parallelization
                          // e.g. sim.svipc = 1 -> split WFS/DM
                          // e.g. sim.svipc = 2 -> parallelize PSFs
                          // e.g. sim.svipc = 3 -> WFS/DM & PSFs
  long    svipc_wfs_nfork;// nb of forks when splitting WFSs (one or more WFS per
                          // fork). if not set, will be min([nwfs,nprocessors]);
  pointer svipc_wfs_forknb;// sim.svipc_wfs_forknb is a vector, with nb of
                          // elements = # of WFS, and which contains
                          // the fork# for each WFS; e.g.:
                          // sim.svipc_wfs_forknb = &([1,1,2,2,3]);
                          // means run WFS 1 and 2 in WFS fork 1
                          //       run WFS 3 and 4 in WFS fork 2
                          //   and run WFS 5 in WFS fork 3.
                          // as a special case, 0 means run in WFS
                          // parent (the one from which all WFS
                          // children were forked). To be implemented
                          // at a later stage.  If this vector is not
                          // specified, we will spread the WFS evenly
                          // (as best as possible) within the N
                          // (=sim.svipc_wfs_nfork) WFS forks.
  long    shmkey;         // shared memory key (there's a default).
                          // Change to run multiple simul in parallel.
  long    semkey;         // shared memory key (there's a default)
  // Internal keywords:
  long    _size;          // Internal. Size of the arrays [pixels]
  float   _cent;          // Internal. Pupil is centered on (_cent,_cent)
};

struct atm_struct
{
  float   dr0at05mic;     // Dr0 at sensing wavelength, at zenith. Required [none].
  pointer screen;         // string vectorptr. Screen file names. Required [none].
  pointer layerfrac;      // float vectorptr. Layer fraction. Sum to one is insured
  // in aoinit. Required [none]
  pointer layerspeed;     // float vectorptr. Layer speed. Required [none]
  pointer layeralt;       // float vectorptr. Layer altitude (m). Specified at Zenith.
  // Required [none]
  pointer winddir;        // Wind dir (use 0 for now)
  // Internal variables
  pointer _layeralt;      // float vectorptr. Actual layer altitude (m), from atm.alt & zen.angle
};

struct opt_struct
{
  string   phasemaps;      // filename of phasemap. Z scale should be in microns
  float    alt;            // float. equivalent altitude in m.
  float    misreg(2);      // float vector. misreg. (similar to DM, see below)
  float    _cent;          // center of the phase maps arrays (similar to sim._cent)
};

struct wfs_struct
{
  string  type;           // valid type are "curvature", "hartmann", "dh",
                          // "pyramid", "zernike" or "user_function" where user_function
                          // is the name of a function defined by the user (see doc).
                          // Required [none]
  long    subsystem;      // Subsystem this WFS belongs to. Optional [1]
  float   lambda;         // WFS wavelength in microns. Required [none]
  long    noise;          // Enable noise (photon noise/read out noise). Optional [0=no].
  float   ron;            // Read out noise in electrons. Optional [0]
  float   darkcurrent;    // dark current in e-/sec/pixel or APD. Optional [0]
  float   gspos(2);       // This WFS guide star position (x<y) in arcsec. Optional [0,0]
  float   gsalt;          // This WFS guide star altitude in meter. 0 for infinity.
                          // ...Specified at zenith. Optional [0]
  float   gsdepth;        // This WFS GS depth in meter (e.g. Na layer thickness).
                          // Specified at zenith. Optional [0]
  float   laserpower;     // this wfs laser power (Na laser only), in Watts projected on sky.
                          // Required when using lasers. Exclusive with gsmag; i.e.
                          // define one OR the other.
  float   gsmag;          // This WFS guide star magnitude. Optional [0]. For LGSs, see above.
  float   skymag;         // This WFS sky mag. Optional [no sky]
  long    filtertilt;     // Filter TT on this sensor? Optional [0=no]
  long    correctUpTT;    // Correct up link tp-tilt ? Optional [0=no]
  float   uplinkgain;     // Up link TT loop gain. Optional [0]
  float   dispzoom;       // Zoom factor for the display (typically around 1). Optional [1]
  float   optthroughput;  // optical throughput to WFS. Optional [1.0]
  long    disjointpup;    // boolean. if set, wfs(n) will be filtered by an array
                          // disjointpup(,,n) that has to be defined by the user
                          // see user_pupil(). Allow for GMT-type topology.
  long    svipc;          // number of parallel process to use for this WFS.
                          // (0 or 1: don't parallelize)
  float   zeropoint;      // zeropoint for the wavefront sensor. Optional [0.]

  // Curvature WFS only keywords:
  pointer nsubperring;    // Long vectorptr. # subapertures per ring. Required [none]
  pointer angleoffset;    // float vectorptr. offset angle for first subaperture of ring.
  float   l;              // Extra focal distance in a F/60 beam. Required [none]
  pointer rint;           // float vectorptr. if set, specify the inner radius for each ring
  pointer rout;           // float vectorptr. if set, specify the outer radius for each ring
  float   fieldstopdiam;  // diameter of field stop in arcsec. Optional [1]. used only
                          // to compute sky contribution (with skymag).

  // Pyramid WFS only keywords:
  float   pyr_mod_ampl;   // pyramid wfs modulation amplitude radius [arcsec]
  long    pyr_mod_npts;   // total number of point along modulation circle [unitless]
  pointer pyr_mod_pos;    // positions for modulation, overwrites ampl and npts [arcsec]
  long    pyr_padding;    // Pad the pupil image to reduce spatial aliasing [unitless]
                          // A pad of 1 means adding wfs.npixpersub pixels
                          // on each side of the pupil image. Typical 0 to 4.
  string  pyr_mod_loc;    // Location of modulation (before/after the field stop.
                          // valid value are "before" or "after"

  // Shack-Hartmann WFS only keywords:
  long    shmethod;       // 1 = simple gradient average, 2=full propagation. Required [none]
  long    shnxsub;        // # of subaperture in telescope diameter. Required [none]
                          // shnxsub also for pyramid.
  long    npixpersub;     // to force npixpersub and bypass constraint that
                          // pupildiam should be a multiple of this number
                          // e.g. to investigate lenslet larger than pupildiam (or mask inpupil)
  float   pixsize;        // Subaperture pixel size in arcsec. Required [none]
  int     npixels;        // Final # of pixels in subaperture. Required [none]
  float   pupoffset(2);   // offset of the whole wfs subs w.r.t telescope aperture [meter]
                          // allow misregistration w.r.t tel pupil and funky configurations
  long    shthmethod;     // 1: yao default, 2: podium, 3: brightest pixels. Required [1]
  float   shthreshold;    // Threshold in computation of subaperture signal, >=0. Optional [0]
  float   biasrmserror;   // rms error on WFS bias in electron. Optional [0]
  float   flatrmserror;   // rms error on WFS flat, referenced to 1. Optional [0]
                          // Typical value can be 0.01
  string  fsname;         // fits file with subaperture amplitude mask. It should have
                          // dimension 2^sdimpow2 square. can be float or long.
  string  fstop;          // "none", "square" or "round" are allowable values
  float   fssize;         // side (square) or diameter (round) of field stop [arcsec]
  float   fsoffset(2);    // offset (arcsec) in x and y (2 elements vector)
  float   kernel;         // FWHM in arcsec of WFS gaussian kernel. Optional.
                          // Default is computed as a function of D/r0
  int     nintegcycles;   // # of cycles/iterations over which to integrate. Optional [1]
  float   fracIllum;      // fraction illuminated to be valid. Optional [0.5]
  float   LLTxy(2);       // 2 element vector with (x,y) of laser projector [m]
  long    centGainOpt;    // Centroid Gain optimization flag. only for LGS (correctupTT and
                          // filtertilt must be set). Optional [0]
  int     rayleighflag;   // set to one to take rayleigh into account
                          // zernike wfs only
  int     nzer;           // # of zernike sensed

  // DH wfs only
  int     ndh;            // # of dh sensed
  int     ndhfiltered;    // # of dh filtered. 2 would filter tip and tilt.

  // Internal keywords:
  int     _initkernels;   // put in wfs struct for svipc sync 2010jun25
  int     _svipc_init_done; // svipc init done for this wfs
  pointer _svipc_subok;   // vector (length=nsub4disp) 0-> skip comput. of spot. 1->do it.
  pointer _fork_subs;   // matrix signing which subap have to be processed by the various
                          // wfs fork (dim=nforkxnsub). 0-> skip, 1-> process.
  pointer _fork_subs2;  // same, but for call to _shwfs_spots2slopesa
  pointer _validsubs;     // 0/1 mark invalid/valid, out of the ones selected for display
  float   _origpixsize;   // Internal.
  int     _rebinfactor;   // fft pixels to big pixels
  float   _gsalt;         // This WFS guide star altitude in meter. 0 for infinity.
  float   _gsdepth;       // This WFS GS depth in meter (e.g. Na layer thickness).
  int     _nsub;          // Internal. Tot # of valid subs.
  int     _nsub4disp;     // Internal. Tot # of subs to display.
  long    _nmes;          // internal. Tot # of measurements.
  pointer _sind;          // Internal: see CurvWFS
  pointer _nsind;         // Internal: see CurvWFS
  pointer _cxdef;         // Internal: see CurvWFS
  pointer _sxdef;         // Internal: see CurvWFS
  pointer _tiltsh;        // Internal: see sh_wfs
  pointer _masks;         // Internal: see sh_wfs
  pointer _fluxpersub;    // Internal: see sh_wfs
  pointer _raylfluxpersub;// Internal: see sh_wfs
  pointer _skyfluxpersub; // Internal: see sh_wfs
  float   _nphotons;      // Internal: see WFS routines
  float   _skynphotons;   // Internal: see WFS routines
  float   _tt(2);         // Internal: WFS measured Tip and tilt
  float   _lastvalidtt(2);// Internal: WFS measured Tip and tilt
  float   _upttcommand(2);// Internal:
  pointer _refmes;        // internal: reference measurement vector
  pointer _tiprefv;       // internal: tip reference meas. vector
  pointer _tiltrefv;      // internal: tilt reference meas. vector
  pointer _tiprefvn;      // internal: tip reference meas. vector. normalized (norm=1)
  pointer _tiltrefvn;     // internal: tilt reference meas. vector. normalized.
  pointer _istart;        //
  pointer _jstart;        //
  pointer _binindices;    //
  int     _binxy;         //
  pointer _centroidw;     //
  pointer _fimage;        //
  pointer _fimage2;       //
  pointer _imistart;      //
  pointer _imjstart;      //
  pointer _imistart2;     //
  pointer _imjstart2;     //
  int     _fimnx;         // x dim of wfs._fimage
  int     _fimny;         // y dim of wfs._fimage
  pointer _fimny2;        // y dim of (possibly split) wfs._fimage for _shwfs_spots2slopes (svipc)
  pointer _yoffset;       // y offset of wfs._fimage for _shwfs_spots2slopes (svipc)
  pointer _bias;          // internal: array of bias error
  pointer _flat;          // internal: array of flat error
  long    _domask;        // internal. flag to do submask calculations in _shwfs
  pointer _submask;       // internal: array. subaperture amplitude mask.
  pointer _kernel;        // internal: kernel for _shwfs. use: dointer or LGS uplink im.
  pointer _kernels;       // internal: subaperture dependant image kernel
  pointer _kerfftr;       // internal: storage of FFTs of kernels
  pointer _kerffti;       // internal: storage of FFTs of kernels
  int     _kernelconv;    // interal: convolve with kernel in _shwfs?
  int     _cyclecounter;  // counter in integration sequence (see nintegcycles above)
  pointer _dispimage;     // image to display (same as fimage except if nintegcycles!=1)
  pointer _x;             // shwfs: X positions of subaperture centers
  pointer _y;             // shwfs: Y positions of subaperture centers
  float   _centroidgain;  // internal: centroid gain if dithering on
  pointer _rayleigh;      // pointer to rayleigh images array for this sensor
  pointer _bckgrdcalib;   // pointer to background array calibration
  int     _bckgrdinit;    // set to one to fill calibration array
  int     _bckgrdsub;     // set to one to subtract background (default)
  pointer _meashist;      // measurement history, useful for nintegcycles > 1
  float   _zeropoint;     // zeropoint for the wavefront sensor.
  pointer _pha2dhc;       // projection matrix phase to DH coefs for this wfs
  pointer _wpha2dhc;      // valid phase points indices
  int     _n12(2);        //
};

struct dm_struct
{
  string  type;           // valid types are "bimorph", "stackarray" "tiptilt",
                          // "dh", "zernike", "kl", "segmented", "aniso" or
                          // "user_function", where user_function is the name of
                          // a function provided by the user. Required [none]
  long    subsystem;      // Subsystem this DM belongs to. Optional [1]
  long    virtual;        // virtual DMs for tomography, don't correct wavefront
  pointer fitvirtualdm;   // which tomographic virtual DMs are used to drive this DM
  string  iffile;         // Influence function file name. Leave it alone.
  long    pitch;          // Actuator pitch (pixel). stackarray/segmented only. Required [none]
  float   alt;            // Conjugation altitude in meter. Specified @ zenith! Optional [0]
  float   hyst;           // DM actuator hysteresis (0. to 1.). Optional [0]
  float   push4imat;      // Voltage to apply for imat calibration. Optional [20].
                          // Note: the default is not OK for many configs. Change at will.
  float   thresholdresp;  // Normalized response threshold for an act. to be valid. Optional [0.3]
  float   unitpervolt;    // Influence function sensitivity in unit/volt. Optional [0.01]
                          // Stackarray: mic/volt, Tip-tilt: arcsec/volt.
  float   maxvolt;        // Saturation voltage (- and +) in volt. Optional [none if not set]
  float   gain;           // loop gain for this DM (total = this x loop.gain). Optional [1]
  float   misreg(2);      // dm misregistration (pixels). optional [0,0]
  long    xflip;          // flip influence functions left/right
  long    yflip;          // flip influence functions up/down
  float   pupoffset(2);   // global offset in pupil of whole actuator pattern [m]
  long    disjointpup;    // boolean. if set, dm(n) will be filtered by an array
                          // disjointpup(,,n) that has to be defined by the user
                          // see user_pupil(). Allow for GMT-type topology.

  // Bimorph-only keywords:
  pointer nelperring;     // long vectorptr. # of elec. per ring, e.g &([6,12,18]). Required [none]
  pointer angleoffset;    // float vectorptr. offset angle for first electrode in ring.
  pointer rint;           // float vectorptr. if set, specify the inner radius for each ring
  pointer rout;           // float vectorptr. if set, specify the outer radius for each ring
  float   supportRadius;  // radius of support points, normalized in pupil radius
  float   supportOffset;  // angle offset of first support points, in degree (default=90)

  // Stackarray-only keywords:
  long    nxact;          // # of act. in pupil diameter. Required [none]
  float   pitchMargin;    // margin to include more corner actuators when creating inf.functions
                          // optional [1.44]
  float   coupling;       // coupling coef in influence function. optional [0.2].
                          // valid values from 0.04 to 0.30.
  string  ecmatfile;      // valid to extrap. actuator matrix (extrap_com). Optional.
  long    noextrap;       // set to 1 if no extrapolated actuators whatsoever are to be used [0]
  long    elt;            // set to 1 if fast dmsum to be used
  long    irexp;          // use old/regular form (irexp=0) or
                          // exp(-(d/irfact)^1.5) model (irexp=1) or
                          // sinc*gaussian (irexp=2)
  float   irfact;         // use when irexp=1 (see above)
  long    filtertilt;     // Filter TT on this DM? Optional [0=no]

  // Zernike-only keywords:
  long    nzer;           // Number of modes, including piston. Required [none]
  long    minzer;         // lowest order zernike, default=1 (piston)

  // Disk-Harmonic only keywords
  long    ndh;            // number of DH modes

  // KL-only keywords:
  long    nkl;            // Number of modes, including piston. Required [none]

  // Segmented only keywords:
  long    nxseg;          // number of segments in long axis (X)
  float   fradius;        // segments are created over a wider area than the
                          // nxseg defined above. Only segments which distance
                          // to the (0,0) pupil coordinates is <= fradius
                          // will be kept (pixels). default dm.pitch*dm.nxseg/2.

  // MMSE and sparse MMSE matrix reconstructor parameters:
  float   regparam;       // regularization parameter
  string  regtype;        // regulatization matrix generation method.
  pointer regmatrix;      // matrix used in the regularization
  // Internal keywords:
  long    _puppixoffset(2);
  float   _alt;           // Actual conjugation altitude in meter, from dm.alt and zen.
  long    _nact;          // Internal. Tot # of actuators.
  pointer _def;           // Internal: Pointer to IF data
  pointer _x;             // Internal: x position of actuator in pixels
  pointer _y;             // Internal: x position of actuator in pixels
  pointer _i1;            //
  pointer _j1;            //
  pointer _ei1;           //
  pointer _ej1;           //
  string  _eiffile;       // Influence function file name for extrap. actuators
  pointer _edef;          // Internal: Pointer to IF data for extrap. actuators
  pointer _ex;            // Internal: x position of extrap. actuator in pixels
  pointer _ey;            // Internal: x position of extrap. actuator in pixels
  long    _enact;         // Internal. Tot # of extrap. actuators.
  long    _n1;            // Internal: position of leftmost pixel in ao._size^2 array
  long    _n2;            // Internal: position of leftmost pixel in ao._size^2 array
  pointer _pupil;         // Internal. Mask for display.
  pointer _vold;          // internal: hysteresis
  pointer _posold;        // internal: hysteresis
  pointer _chpos;         // internal: hysteresis
  pointer _chv;           // internal: hysteresis
  pointer _dir;           // internal: hysteresis
  pointer _delta;         // internal: hysteresis
  pointer _command;       // pointer to command vector
  pointer _extrapcmat;    // extrapolation matrix: extrap_com = extrapmat(,+)*valid_com(+)
  int     _eltdefsize;    // size of def in case elt=1
  pointer _regmatrix;     // regularization matrix used, if any
  pointer _fMat;          // fitting matrix for tomography
};

struct mat_struct
{
  string  method;         // reconstruction method: "svd" (default), "mmse", "mmse-sparse"
  pointer condition;      // float vecorptr. Condition numbers for SVD, per subsystem. Required [none]
  long    sparse_MR;      // maximum number of rows for sparse method
  long    sparse_MN;      // maximum number of elements for sparse method
  float   sparse_thresh;  // threshold for non-zero sparse elements
  float   sparse_pcgtol;  // tolerance for reconstruction, default = 1e-3
  string  file;           // iMat and cMat filename. Leave it alone.
  // fitting parameters for tomographic reconstruction
  long  fit_simple;            // 0 or 1, default = 0. Simple optimizes on the optical axis and only works if the tomographic DM is the same as the corresponding virtual DMs, but is faster.
  // the following parameters only apply to "mmse" fitting
  long    fit_subsamp;    // subsampling the phase for fitting matrix (set to larger than 1 for speed), default = 1
  long    fit_target;     // which target to optimize fitting for, default = 1
  float   fit_minval;     // minimum value for sparse fitting matrix, default = 1e-2
};

struct tel_struct
{
  float   diam;               // Telescope diameter in meters. Required [none]
  float   cobs;               // Central obstruction / telescope diameter ratio. Optional [0]

  // TIP vibrations parameters
  float   tipvib_white_rms;   // rms [arcsec] of vibration white noise
  float   tipvib_1overf_rms;  // rms [arcsec] of vibration 1/f noise (from 1 Hz to cutoff)
  pointer tipvib_peaks;       // positions [Hz] of vibration peak in PSD
  pointer tipvib_peaks_rms;   // rms [arcsec] of each vibration peaks (defined above)
  pointer tipvib_peaks_width; // width [Hz] of each vibration peaks (default 1 freq bin)

  // TILT vibrations parameters
  float   tiltvib_white_rms;  // rms [arcsec] of vibration white noise
  float   tiltvib_1overf_rms; // rms [arcsec] of vibration 1/f noise (from 1 Hz to cutoff)
  pointer tiltvib_peaks;      // positions [Hz] of vibration peak in PSD
  pointer tiltvib_peaks_rms;  // rms [arcsec] of each vibration peaks (defined above)
  pointer tiltvib_peaks_width;// width [Hz] of each vibration peaks (default 1 freq bin)
};

struct target_struct
{
  pointer lambda;         // float vectorptr. Image wavelengths in micron. Required [none]
  pointer xposition;      // float vectorptr. X positions in arcsec. Required [none]
  pointer yposition;      // float vectorptr. Y positions in arcsec. Required [none]
  pointer dispzoom;       // float vectorptr. Display zoom (typically around 1.). Optional [1.]

  // Internal keywords
  long    _ntarget;       // Internal: # of target
  long    _nlambda;       // Internal: # of lambda
};

struct gs_struct
{
  float   zeropoint;      // Photometric zero point (#photons@pupil/s/full_aper, mag0 star).
                          // Required [none]
  float   zenithangle;    // zenith angle. Optional [0.]. The zenith angle is used to compute:
                          // - r0 off-zenith
                          // - atmopheric turbulence layer altitude
                          // - LGS altitude and thickness of Na Layer
                          // - LGS brighness
                          // note that dm altitude is unchanged.
  float  lgsreturnperwatt;// Sodium LGS return in photons/cm2/s at entrance pupil.
                          // Specified at zenith. Modified by gs.zenithangle. Optional [22.]
                          // basically, you have to fold in this the sodium density
                          // and your model of return.
};

struct loop_struct
{
  float   gain;            // Loop gain. Optional, but important! [0]
  float   leak;            // leak term (0 means no leak) [0]
  pointer gainho;          // vector of higher order gains (starting at 2nd order)
  pointer leakho;          // vector of higher order leaks (starting at 2nd order)
  long    framedelay;      // loop delay (# of frames). Optional [0]
                           // Regular CCD 1 frame integration -> framedelay=1
                           // + readout & Calculation -> framedelay=2
  long    niter;           // # of total iteration. Required [none]
  float   ittime;          // Iteration time in seconds. Required [none]
  long    startskip;       // # iter to skip before collecting statistics. Optional [10]
  long    skipevery;       // skip by "skipby" every "skipevery" iterations. Optional [0=none]
  long    skipby;          // see above. this is to get better statistical
                           // coverage. Optional [10000]
  long    stats_every;     // compute stats every so many iteration (default 4)
  long    jumps2swapscreen;//number of jumps (i.e. niter/skipevery) after which screens
                           // will be swapped (rotation, 2->1, 3->2... 1->last
  string  modalgainfile;   // Name of file with mode gains. Optional.
  //float   dithering;     // TT dithering for centroid gain (volts).
  string  method;          // "closed-loop", "open-loop", "pseudo open-loop"
};
