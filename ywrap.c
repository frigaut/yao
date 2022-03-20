/* codger-generated yorick package wrapper file */
#include "play.h"
#include "ydata.h"

/*----------------begin yao_fast.i */
extern BuiltIn Y__calc_psf_fast;

extern int _calc_psf_fast(void *, void *, void *, int , int , 
  float , int );
void
Y__calc_psf_fast(int n)
{
  if (n!=7) YError("_calc_psf_fast takes exactly 7 arguments");
  PushIntValue(_calc_psf_fast(yarg_sp(6), yarg_sp(5), yarg_sp(4), 
    yarg_si(3), yarg_si(2), yarg_sf(1), yarg_si(0)));
}

extern BuiltIn Y__init_fftw_plans;

extern int _init_fftw_plans(int );
void
Y__init_fftw_plans(int n)
{
  if (n!=1) YError("_init_fftw_plans takes exactly 1 arguments");
  PushIntValue(_init_fftw_plans(yarg_si(0)));
}

extern BuiltIn Y__init_fftw_plan;

extern int _init_fftw_plan(int );
void
Y__init_fftw_plan(int n)
{
  if (n!=1) YError("_init_fftw_plan takes exactly 1 arguments");
  PushIntValue(_init_fftw_plan(yarg_si(0)));
}

extern BuiltIn Y__import_wisdom;

extern int _import_wisdom(char *);
void
Y__import_wisdom(int n)
{
  if (n!=1) YError("_import_wisdom takes exactly 1 arguments");
  PushIntValue(_import_wisdom(yarg_sq(0)));
}

extern BuiltIn Y__export_wisdom;

extern int _export_wisdom(char *);
void
Y__export_wisdom(int n)
{
  if (n!=1) YError("_export_wisdom takes exactly 1 arguments");
  PushIntValue(_export_wisdom(yarg_sq(0)));
}

extern BuiltIn Y__set_sincos_approx;

extern void _set_sincos_approx(int );
void
Y__set_sincos_approx(int n)
{
  if (n!=1) YError("_set_sincos_approx takes exactly 1 arguments");
  _set_sincos_approx(yarg_si(0));
}

extern BuiltIn Y__get_sincos_approx;
extern BuiltIn Y__sinecosinef;

extern void _sinecosinef(float , void *, void *);
void
Y__sinecosinef(int n)
{
  if (n!=3) YError("_sinecosinef takes exactly 3 arguments");
  _sinecosinef(yarg_sf(2), yarg_sp(1), yarg_sp(0));
}

extern BuiltIn Y__fftVE;

extern int _fftVE(void *, void *, int , int );
void
Y__fftVE(int n)
{
  if (n!=4) YError("_fftVE takes exactly 4 arguments");
  PushIntValue(_fftVE(yarg_sp(3), yarg_sp(2), yarg_si(1), 
    yarg_si(0)));
}

extern BuiltIn Y__fftVE2;

extern int _fftVE2(void *, void *, int , int );
void
Y__fftVE2(int n)
{
  if (n!=4) YError("_fftVE2 takes exactly 4 arguments");
  PushIntValue(_fftVE2(yarg_sp(3), yarg_sp(2), yarg_si(1), 
    yarg_si(0)));
}

extern BuiltIn Y_embed_image;

extern int embed_image(float *, int , int , float *, int , int , 
  int , int );
void
Y_embed_image(int n)
{
  if (n!=8) YError("embed_image takes exactly 8 arguments");
  PushIntValue(embed_image(yarg_f(7,0), yarg_si(6), yarg_si(5), 
    yarg_f(4,0), yarg_si(3), yarg_si(2), yarg_si(1), yarg_si(0)));
}

extern BuiltIn Y__shwfs_phase2spots;

extern int _shwfs_phase2spots(float *, float *, float , float *, 
  int , int *, int *, int , int , int , int , long , float *, 
  float *, int , float *, float *, float *, int , int , int *, 
  int , int , int , float *, float *, float *, float *, int , 
  float *, float *, int *, int *, int *, int , int , float *, 
  float *, float *, float , int , float *, int , int , int );
void
Y__shwfs_phase2spots(int n)
{
  if (n!=45) YError("_shwfs_phase2spots takes exactly 45 arguments");
  PushIntValue(_shwfs_phase2spots(yarg_f(44,0), yarg_f(43,0), yarg_sf(42), 
    yarg_f(41,0), yarg_si(40), yarg_i(39,0), yarg_i(38,0), yarg_si(37), 
    yarg_si(36), yarg_si(35), yarg_si(34), yarg_sl(33), 
    yarg_f(32,0), yarg_f(31,0), yarg_si(30), yarg_f(29,0), yarg_f(28,0), 
    yarg_f(27,0), yarg_si(26), yarg_si(25), yarg_i(24,0), yarg_si(23), 
    yarg_si(22), yarg_si(21), yarg_f(20,0), yarg_f(19,0), yarg_f(18,0), 
    yarg_f(17,0), yarg_si(16), yarg_f(15,0), yarg_f(14,0), yarg_i(13,0), 
    yarg_i(12,0), yarg_i(11,0), yarg_si(10), yarg_si(9), yarg_f(8,0), 
    yarg_f(7,0), yarg_f(6,0), yarg_sf(5), yarg_si(4), yarg_f(3,0), 
    yarg_si(2), yarg_si(1), yarg_si(0)));
}

extern BuiltIn Y__shwfs_spots2slopes;

extern int _shwfs_spots2slopes(float *, int *, int *, int , int , 
  int , int , int , float *, long , float *, float *, float *, 
  float , float , long , float *, int , int , int *, int *, int , 
  float *);
void
Y__shwfs_spots2slopes(int n)
{
  if (n!=23) YError("_shwfs_spots2slopes takes exactly 23 arguments");
  PushIntValue(_shwfs_spots2slopes(yarg_f(22,0), yarg_i(21,0), yarg_i(20,0), 
    yarg_si(19), yarg_si(18), yarg_si(17), yarg_si(16), yarg_si(15), 
    yarg_f(14,0), yarg_sl(13), yarg_f(12,0), yarg_f(11,0), yarg_f(10,0), 
    yarg_sf(9), yarg_sf(8), yarg_sl(7), yarg_f(6,0), 
    yarg_si(5), yarg_si(4), yarg_i(3,0), yarg_i(2,0), yarg_si(1), 
    yarg_f(0,0)));
}

extern BuiltIn Y__shwfs_simple;

extern int _shwfs_simple(float *, float *, float , float *, int , 
  int , int *, int *, int , int , int , float , float *);
void
Y__shwfs_simple(int n)
{
  if (n!=13) YError("_shwfs_simple takes exactly 13 arguments");
  PushIntValue(_shwfs_simple(yarg_f(12,0), yarg_f(11,0), yarg_sf(10), 
    yarg_f(9,0), yarg_si(8), yarg_si(7), yarg_i(6,0), yarg_i(5,0), 
    yarg_si(4), yarg_si(3), yarg_si(2), yarg_sf(1), 
    yarg_f(0,0)));
}

extern BuiltIn Y__cwfs;

extern int _cwfs(float *, float *, float , float *, float *, 
  float *, int , int *, int *, int , float *, float *, float , 
  float , float , float , float , int , float *);
void
Y__cwfs(int n)
{
  if (n!=19) YError("_cwfs takes exactly 19 arguments");
  PushIntValue(_cwfs(yarg_f(18,0), yarg_f(17,0), yarg_sf(16), yarg_f(15,0), 
    yarg_f(14,0), yarg_f(13,0), yarg_si(12), yarg_i(11,0), yarg_i(10,0), 
    yarg_si(9), yarg_f(8,0), yarg_f(7,0), yarg_sf(6), yarg_sf(5), 
    yarg_sf(4), yarg_sf(3), yarg_sf(2), yarg_si(1), 
    yarg_f(0,0)));
}


/*----------------begin yao_utils.i */
extern BuiltIn Y_usleep;
extern BuiltIn Y__mynoop2;

extern int _mynoop2(void *, int , int , void *, int , int , int );
void
Y__mynoop2(int n)
{
  if (n!=7) YError("_mynoop2 takes exactly 7 arguments");
  PushIntValue(_mynoop2(yarg_sp(6), yarg_si(5), yarg_si(4), 
    yarg_sp(3), yarg_si(2), yarg_si(1), yarg_si(0)));
}

extern BuiltIn Y__dmsum;

extern void _dmsum(void *, int , int , int , void *, void *);
void
Y__dmsum(int n)
{
  if (n!=6) YError("_dmsum takes exactly 6 arguments");
  _dmsum(yarg_sp(5), yarg_si(4), yarg_si(3), yarg_si(2), 
    yarg_sp(1), yarg_sp(0));
}

extern BuiltIn Y__dmsum2;

extern void _dmsum2(void *, void *, long , long , void *, void *, 
  long );
void
Y__dmsum2(int n)
{
  if (n!=7) YError("_dmsum2 takes exactly 7 arguments");
  _dmsum2(yarg_sp(6), yarg_sp(5), yarg_sl(4), yarg_sl(3), 
    yarg_sp(2), yarg_sp(1), yarg_sl(0));
}

extern BuiltIn Y__dmsumelt;

extern void _dmsumelt(void *, int , int , int , void *, void *, 
  void *, void *, int , int );
void
Y__dmsumelt(int n)
{
  if (n!=10) YError("_dmsumelt takes exactly 10 arguments");
  _dmsumelt(yarg_sp(9), yarg_si(8), yarg_si(7), yarg_si(6), 
    yarg_sp(5), yarg_sp(4), yarg_sp(3), yarg_sp(2), yarg_si(1), 
    yarg_si(0));
}

extern BuiltIn Y__get2dPhase;

extern int _get2dPhase(void *, int , int , int , void *, void *, 
  int , int , void *, void *, void *, void *);
void
Y__get2dPhase(int n)
{
  if (n!=12) YError("_get2dPhase takes exactly 12 arguments");
  PushIntValue(_get2dPhase(yarg_sp(11), yarg_si(10), yarg_si(9), 
    yarg_si(8), yarg_sp(7), yarg_sp(6), yarg_si(5), yarg_si(4), 
    yarg_sp(3), yarg_sp(2), yarg_sp(1), yarg_sp(0)));
}

extern BuiltIn Y__cosf;

extern int _cosf(void *, int );
void
Y__cosf(int n)
{
  if (n!=2) YError("_cosf takes exactly 2 arguments");
  PushIntValue(_cosf(yarg_sp(1), yarg_si(0)));
}

extern BuiltIn Y__sinf;

extern int _sinf(void *, int );
void
Y__sinf(int n)
{
  if (n!=2) YError("_sinf takes exactly 2 arguments");
  PushIntValue(_sinf(yarg_sp(1), yarg_si(0)));
}


/*----------------list include files */

static char *y0_includes[] = {
  "yao_fast.i",
  "yao_utils.i",
  0,
  0
};

/*----------------collect pointers and names */

static BuiltIn *y0_routines[] = {
  &Y__calc_psf_fast,
  &Y__init_fftw_plans,
  &Y__init_fftw_plan,
  &Y__import_wisdom,
  &Y__export_wisdom,
  &Y__set_sincos_approx,
  &Y__get_sincos_approx,
  &Y__sinecosinef,
  &Y__fftVE,
  &Y__fftVE2,
  &Y_embed_image,
  &Y__shwfs_phase2spots,
  &Y__shwfs_spots2slopes,
  &Y__shwfs_simple,
  &Y__cwfs,
  &Y_usleep,
  &Y__mynoop2,
  &Y__dmsum,
  &Y__dmsum2,
  &Y__dmsumelt,
  &Y__get2dPhase,
  &Y__cosf,
  &Y__sinf,
  0
};

static void *y0_values[] = {
  0
};

static char *y0_names[] = {
  "_calc_psf_fast",
  "_init_fftw_plans",
  "_init_fftw_plan",
  "_import_wisdom",
  "_export_wisdom",
  "_set_sincos_approx",
  "_get_sincos_approx",
  "_sinecosinef",
  "_fftVE",
  "_fftVE2",
  "embed_image",
  "_shwfs_phase2spots",
  "_shwfs_spots2slopes",
  "_shwfs_simple",
  "_cwfs",
  "usleep",
  "_mynoop2",
  "_dmsum",
  "_dmsum2",
  "_dmsumelt",
  "_get2dPhase",
  "_cosf",
  "_sinf",
  0
};

/*----------------define package initialization function */

PLUG_EXPORT char *yk_yao(char ***,
                         BuiltIn ***, void ***, char ***);
static char *y0_pkgname = "yao";

char *
yk_yao(char ***ifiles,
       BuiltIn ***code, void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}
