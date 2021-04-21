/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__BK
#define _nrn_initial _nrn_initial__BK
#define nrn_cur _nrn_cur__BK
#define _nrn_current _nrn_current__BK
#define nrn_jacob _nrn_jacob__BK
#define nrn_state _nrn_state__BK
#define _net_receive _net_receive__BK 
#define rates rates__BK 
#define states states__BK 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define conductance _p[1]
#define c_instances _p[2]
#define c_forwardRate_TIME_SCALE _p[3]
#define c_forwardRate_VOLT_SCALE _p[4]
#define c_forwardRate_CONC_SCALE _p[5]
#define c_reverseRate_TIME_SCALE _p[6]
#define c_reverseRate_VOLT_SCALE _p[7]
#define c_reverseRate_CONC_SCALE _p[8]
#define c_q10Settings_q10Factor _p[9]
#define c_q10Settings_experimentalTemp _p[10]
#define c_q10Settings_TENDEGREES _p[11]
#define gion _p[12]
#define c_forwardRate_V _p[13]
#define c_forwardRate_ca_conc _p[14]
#define c_forwardRate_r _p[15]
#define c_reverseRate_V _p[16]
#define c_reverseRate_ca_conc _p[17]
#define c_reverseRate_r _p[18]
#define c_q10Settings_q10 _p[19]
#define c_rateScale _p[20]
#define c_alpha _p[21]
#define c_beta _p[22]
#define c_fcond _p[23]
#define c_inf _p[24]
#define c_tau _p[25]
#define conductanceScale _p[26]
#define fopen0 _p[27]
#define fopen _p[28]
#define g _p[29]
#define c_q _p[30]
#define temperature _p[31]
#define ek _p[32]
#define ik _p[33]
#define cai _p[34]
#define cao _p[35]
#define rate_c_q _p[36]
#define Dc_q _p[37]
#define v _p[38]
#define _g _p[39]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ik	*_ppvar[2]._pval
#define _ion_dikdv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_BK", _hoc_setdata,
 "rates_BK", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_BK", "S/cm2",
 "conductance_BK", "uS",
 "c_forwardRate_TIME_SCALE_BK", "ms",
 "c_forwardRate_VOLT_SCALE_BK", "mV",
 "c_forwardRate_CONC_SCALE_BK", "mM",
 "c_reverseRate_TIME_SCALE_BK", "ms",
 "c_reverseRate_VOLT_SCALE_BK", "mV",
 "c_reverseRate_CONC_SCALE_BK", "mM",
 "c_q10Settings_experimentalTemp_BK", "K",
 "c_q10Settings_TENDEGREES_BK", "K",
 "gion_BK", "S/cm2",
 "c_forwardRate_r_BK", "kHz",
 "c_reverseRate_r_BK", "kHz",
 "c_alpha_BK", "kHz",
 "c_beta_BK", "kHz",
 "c_tau_BK", "ms",
 "g_BK", "uS",
 0,0
};
 static double c_q0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"BK",
 "gmax_BK",
 "conductance_BK",
 "c_instances_BK",
 "c_forwardRate_TIME_SCALE_BK",
 "c_forwardRate_VOLT_SCALE_BK",
 "c_forwardRate_CONC_SCALE_BK",
 "c_reverseRate_TIME_SCALE_BK",
 "c_reverseRate_VOLT_SCALE_BK",
 "c_reverseRate_CONC_SCALE_BK",
 "c_q10Settings_q10Factor_BK",
 "c_q10Settings_experimentalTemp_BK",
 "c_q10Settings_TENDEGREES_BK",
 0,
 "gion_BK",
 "c_forwardRate_V_BK",
 "c_forwardRate_ca_conc_BK",
 "c_forwardRate_r_BK",
 "c_reverseRate_V_BK",
 "c_reverseRate_ca_conc_BK",
 "c_reverseRate_r_BK",
 "c_q10Settings_q10_BK",
 "c_rateScale_BK",
 "c_alpha_BK",
 "c_beta_BK",
 "c_fcond_BK",
 "c_inf_BK",
 "c_tau_BK",
 "conductanceScale_BK",
 "fopen0_BK",
 "fopen_BK",
 "g_BK",
 0,
 "c_q_BK",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 40, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	conductance = 1e-05;
 	c_instances = 1;
 	c_forwardRate_TIME_SCALE = 1;
 	c_forwardRate_VOLT_SCALE = 1;
 	c_forwardRate_CONC_SCALE = 1e+06;
 	c_reverseRate_TIME_SCALE = 1;
 	c_reverseRate_VOLT_SCALE = 1;
 	c_reverseRate_CONC_SCALE = 1e+06;
 	c_q10Settings_q10Factor = 3;
 	c_q10Settings_experimentalTemp = 303.15;
 	c_q10Settings_TENDEGREES = 10;
 	_prop->param = _p;
 	_prop->param_size = 40;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 prop_ion = need_memb(_k_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _BK_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", 2.0);
 	ion_reg("k", 1.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 40, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 BK /home/jovyan/work/NeuroML_Examples/Rocher/LEMSFILES3/x86_64/BK.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Mod file for component: Component(id=BK type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dc_q = rate_c_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dc_q = Dc_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    c_q = c_q - dt*(- ( rate_c_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
   double _lcaConc ;
 _lcaConc = cai ;
   c_forwardRate_V = v / c_forwardRate_VOLT_SCALE ;
   c_forwardRate_ca_conc = _lcaConc / c_forwardRate_CONC_SCALE ;
   c_forwardRate_r = ( 7.0 / ( 1.0 + ( 0.0015 * ( exp ( c_forwardRate_V / - 11.765 ) ) / ( c_forwardRate_ca_conc * 1e6 ) ) ) ) / c_forwardRate_TIME_SCALE ;
   c_reverseRate_V = v / c_reverseRate_VOLT_SCALE ;
   c_reverseRate_ca_conc = _lcaConc / c_reverseRate_CONC_SCALE ;
   c_reverseRate_r = ( 1.0 / ( 1.0 + ( c_reverseRate_ca_conc * 1e6 ) / ( 0.00015 * ( exp ( c_reverseRate_V / - 11.765 ) ) ) ) ) / c_reverseRate_TIME_SCALE ;
   c_q10Settings_q10 = pow( c_q10Settings_q10Factor , ( ( temperature - c_q10Settings_experimentalTemp ) / c_q10Settings_TENDEGREES ) ) ;
   c_rateScale = c_q10Settings_q10 ;
   c_alpha = c_forwardRate_r ;
   c_beta = c_reverseRate_r ;
   c_fcond = pow( c_q , c_instances ) ;
   c_inf = c_alpha / ( c_alpha + c_beta ) ;
   c_tau = 1.0 / ( ( c_alpha + c_beta ) * c_rateScale ) ;
   rate_c_q = ( c_inf - c_q ) / c_tau ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  c_q = c_q0;
 {
   ek = - 88.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   c_q = c_inf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   conductanceScale = 1.0 ;
   fopen0 = c_fcond ;
   fopen = conductanceScale * fopen0 ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ik = gion * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(c_q) - _p;  _dlist1[0] = &(Dc_q) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/jovyan/work/NeuroML_Examples/Rocher/LEMSFILES3/BK.mod";
static const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=BK type=ionChannelHH)\n"
  "\n"
  "COMMENT\n"
  "\n"
  "    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)\n"
  "         org.neuroml.export  v1.7.0\n"
  "         org.neuroml.model   v1.7.0\n"
  "         jLEMS               v0.10.2\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX BK\n"
  "    USEION ca READ cai,cao VALENCE 2\n"
  "    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE c_instances                       : parameter\n"
  "    \n"
  "    RANGE c_alpha                           : exposure\n"
  "    \n"
  "    RANGE c_beta                            : exposure\n"
  "    \n"
  "    RANGE c_tau                             : exposure\n"
  "    \n"
  "    RANGE c_inf                             : exposure\n"
  "    \n"
  "    RANGE c_rateScale                       : exposure\n"
  "    \n"
  "    RANGE c_fcond                           : exposure\n"
  "    RANGE c_forwardRate_TIME_SCALE          : parameter\n"
  "    RANGE c_forwardRate_VOLT_SCALE          : parameter\n"
  "    RANGE c_forwardRate_CONC_SCALE          : parameter\n"
  "    \n"
  "    RANGE c_forwardRate_r                   : exposure\n"
  "    RANGE c_reverseRate_TIME_SCALE          : parameter\n"
  "    RANGE c_reverseRate_VOLT_SCALE          : parameter\n"
  "    RANGE c_reverseRate_CONC_SCALE          : parameter\n"
  "    \n"
  "    RANGE c_reverseRate_r                   : exposure\n"
  "    RANGE c_q10Settings_q10Factor           : parameter\n"
  "    RANGE c_q10Settings_experimentalTemp    : parameter\n"
  "    RANGE c_q10Settings_TENDEGREES          : parameter\n"
  "    \n"
  "    RANGE c_q10Settings_q10                 : exposure\n"
  "    RANGE c_forwardRate_V                   : derived variable\n"
  "    RANGE c_forwardRate_ca_conc             : derived variable\n"
  "    RANGE c_reverseRate_V                   : derived variable\n"
  "    RANGE c_reverseRate_ca_conc             : derived variable\n"
  "    RANGE conductanceScale                  : derived variable\n"
  "    RANGE fopen0                            : derived variable\n"
  "    \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    \n"
  "    (nA) = (nanoamp)\n"
  "    (uA) = (microamp)\n"
  "    (mA) = (milliamp)\n"
  "    (A) = (amp)\n"
  "    (mV) = (millivolt)\n"
  "    (mS) = (millisiemens)\n"
  "    (uS) = (microsiemens)\n"
  "    (molar) = (1/liter)\n"
  "    (kHz) = (kilohertz)\n"
  "    (mM) = (millimolar)\n"
  "    (um) = (micrometer)\n"
  "    (umol) = (micromole)\n"
  "    (S) = (siemens)\n"
  "    \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    \n"
  "    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!\n"
  "    \n"
  "    conductance = 1.0E-5 (uS)\n"
  "    c_instances = 1 \n"
  "    c_forwardRate_TIME_SCALE = 1 (ms)\n"
  "    c_forwardRate_VOLT_SCALE = 1 (mV)\n"
  "    c_forwardRate_CONC_SCALE = 1000000 (mM)\n"
  "    c_reverseRate_TIME_SCALE = 1 (ms)\n"
  "    c_reverseRate_VOLT_SCALE = 1 (mV)\n"
  "    c_reverseRate_CONC_SCALE = 1000000 (mM)\n"
  "    c_q10Settings_q10Factor = 3 \n"
  "    c_q10Settings_experimentalTemp = 303.15 (K)\n"
  "    c_q10Settings_TENDEGREES = 10 (K)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    ek (mV)\n"
  "    ik (mA/cm2)\n"
  "    \n"
  "    cai (mM)\n"
  "    \n"
  "    cao (mM)\n"
  "    \n"
  "    \n"
  "    c_forwardRate_V                        : derived variable\n"
  "    \n"
  "    c_forwardRate_ca_conc                  : derived variable\n"
  "    \n"
  "    c_forwardRate_r (kHz)                  : derived variable\n"
  "    \n"
  "    c_reverseRate_V                        : derived variable\n"
  "    \n"
  "    c_reverseRate_ca_conc                  : derived variable\n"
  "    \n"
  "    c_reverseRate_r (kHz)                  : derived variable\n"
  "    \n"
  "    c_q10Settings_q10                      : derived variable\n"
  "    \n"
  "    c_rateScale                            : derived variable\n"
  "    \n"
  "    c_alpha (kHz)                          : derived variable\n"
  "    \n"
  "    c_beta (kHz)                           : derived variable\n"
  "    \n"
  "    c_fcond                                : derived variable\n"
  "    \n"
  "    c_inf                                  : derived variable\n"
  "    \n"
  "    c_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopen0                                 : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_c_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    c_q  \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    ek = -88.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    c_q = c_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=BK type=ionChannelHH), from conductanceScaling; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    conductanceScale = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=BK type=ionChannelHH), from gates; Component(id=c type=gateHHrates)\n"
  "    ? multiply applied to all instances of fcond in: <gates> ([Component(id=c type=gateHHrates)]))\n"
  "    fopen0 = c_fcond ? path based, prefix = \n"
  "    \n"
  "    fopen = conductanceScale  *  fopen0 ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ik = gion * (v - ek)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    c_q' = rate_c_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    LOCAL caConc\n"
  "    \n"
  "    caConc = cai\n"
  "    \n"
  "    c_forwardRate_V = v /  c_forwardRate_VOLT_SCALE ? evaluable\n"
  "    c_forwardRate_ca_conc = caConc /  c_forwardRate_CONC_SCALE ? evaluable\n"
  "    c_forwardRate_r = (7/(1 + (0.0015 * (exp ( c_forwardRate_V /-11.765))/( c_forwardRate_ca_conc  * 1e6)))) /  c_forwardRate_TIME_SCALE ? evaluable\n"
  "    c_reverseRate_V = v /  c_reverseRate_VOLT_SCALE ? evaluable\n"
  "    c_reverseRate_ca_conc = caConc /  c_reverseRate_CONC_SCALE ? evaluable\n"
  "    c_reverseRate_r = (1/(1 + ( c_reverseRate_ca_conc  * 1e6)/(0.00015* (exp ( c_reverseRate_V /-11.765)) ))) /  c_reverseRate_TIME_SCALE ? evaluable\n"
  "    c_q10Settings_q10 = c_q10Settings_q10Factor ^((temperature -  c_q10Settings_experimentalTemp )/ c_q10Settings_TENDEGREES ) ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=c type=gateHHrates), from q10Settings; Component(id=null type=q10ExpTemp)\n"
  "    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10ExpTemp)]))\n"
  "    c_rateScale = c_q10Settings_q10 ? path based, prefix = c_\n"
  "    \n"
  "    ? DerivedVariable is based on path: forwardRate/r, on: Component(id=c type=gateHHrates), from forwardRate; Component(id=null type=Golgi_KC_c_alpha_rate)\n"
  "    c_alpha = c_forwardRate_r ? path based, prefix = c_\n"
  "    \n"
  "    ? DerivedVariable is based on path: reverseRate/r, on: Component(id=c type=gateHHrates), from reverseRate; Component(id=null type=Golgi_KC_c_beta_rate)\n"
  "    c_beta = c_reverseRate_r ? path based, prefix = c_\n"
  "    \n"
  "    c_fcond = c_q ^ c_instances ? evaluable\n"
  "    c_inf = c_alpha /( c_alpha + c_beta ) ? evaluable\n"
  "    c_tau = 1/(( c_alpha + c_beta ) *  c_rateScale ) ? evaluable\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_c_q = ( c_inf  -  c_q ) /  c_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "}\n"
  "\n"
  ;
#endif
