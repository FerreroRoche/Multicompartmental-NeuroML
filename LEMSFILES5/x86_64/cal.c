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
 
#define nrn_init _nrn_init__cal
#define _nrn_initial _nrn_initial__cal
#define nrn_cur _nrn_cur__cal
#define _nrn_current _nrn_current__cal
#define nrn_jacob _nrn_jacob__cal
#define nrn_state _nrn_state__cal
#define _net_receive _net_receive__cal 
#define rates rates__cal 
#define states states__cal 
 
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
#define k_instances _p[2]
#define k_steadyState_rate _p[3]
#define k_steadyState_midpoint _p[4]
#define k_steadyState_scale _p[5]
#define k_timeCourse_tau _p[6]
#define l_instances _p[7]
#define l_steadyState_rate _p[8]
#define l_steadyState_midpoint _p[9]
#define l_steadyState_scale _p[10]
#define l_timeCourse_TIME_SCALE _p[11]
#define l_timeCourse_VOLT_SCALE _p[12]
#define gion _p[13]
#define k_steadyState_x _p[14]
#define k_timeCourse_t _p[15]
#define k_rateScale _p[16]
#define k_fcond _p[17]
#define k_inf _p[18]
#define k_tauUnscaled _p[19]
#define k_tau _p[20]
#define l_steadyState_x _p[21]
#define l_timeCourse_V _p[22]
#define l_timeCourse_t _p[23]
#define l_rateScale _p[24]
#define l_fcond _p[25]
#define l_inf _p[26]
#define l_tauUnscaled _p[27]
#define l_tau _p[28]
#define conductanceScale _p[29]
#define fopen0 _p[30]
#define fopen _p[31]
#define g _p[32]
#define k_q _p[33]
#define l_q _p[34]
#define temperature _p[35]
#define eca _p[36]
#define ica _p[37]
#define rate_k_q _p[38]
#define rate_l_q _p[39]
#define Dk_q _p[40]
#define Dl_q _p[41]
#define v _p[42]
#define _g _p[43]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_dicadv	*_ppvar[1]._pval
 
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
 "setdata_cal", _hoc_setdata,
 "rates_cal", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_cal", "S/cm2",
 "conductance_cal", "uS",
 "k_steadyState_midpoint_cal", "mV",
 "k_steadyState_scale_cal", "mV",
 "k_timeCourse_tau_cal", "ms",
 "l_steadyState_midpoint_cal", "mV",
 "l_steadyState_scale_cal", "mV",
 "l_timeCourse_TIME_SCALE_cal", "ms",
 "l_timeCourse_VOLT_SCALE_cal", "mV",
 "gion_cal", "S/cm2",
 "k_timeCourse_t_cal", "ms",
 "k_tauUnscaled_cal", "ms",
 "k_tau_cal", "ms",
 "l_timeCourse_t_cal", "ms",
 "l_tauUnscaled_cal", "ms",
 "l_tau_cal", "ms",
 "g_cal", "uS",
 0,0
};
 static double delta_t = 0.01;
 static double k_q0 = 0;
 static double l_q0 = 0;
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
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cal",
 "gmax_cal",
 "conductance_cal",
 "k_instances_cal",
 "k_steadyState_rate_cal",
 "k_steadyState_midpoint_cal",
 "k_steadyState_scale_cal",
 "k_timeCourse_tau_cal",
 "l_instances_cal",
 "l_steadyState_rate_cal",
 "l_steadyState_midpoint_cal",
 "l_steadyState_scale_cal",
 "l_timeCourse_TIME_SCALE_cal",
 "l_timeCourse_VOLT_SCALE_cal",
 0,
 "gion_cal",
 "k_steadyState_x_cal",
 "k_timeCourse_t_cal",
 "k_rateScale_cal",
 "k_fcond_cal",
 "k_inf_cal",
 "k_tauUnscaled_cal",
 "k_tau_cal",
 "l_steadyState_x_cal",
 "l_timeCourse_V_cal",
 "l_timeCourse_t_cal",
 "l_rateScale_cal",
 "l_fcond_cal",
 "l_inf_cal",
 "l_tauUnscaled_cal",
 "l_tau_cal",
 "conductanceScale_cal",
 "fopen0_cal",
 "fopen_cal",
 "g_cal",
 0,
 "k_q_cal",
 "l_q_cal",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 44, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	conductance = 1e-05;
 	k_instances = 3;
 	k_steadyState_rate = 1;
 	k_steadyState_midpoint = -61;
 	k_steadyState_scale = 4.2;
 	k_timeCourse_tau = 1;
 	l_instances = 1;
 	l_steadyState_rate = 1;
 	l_steadyState_midpoint = -85.5;
 	l_steadyState_scale = -8.5;
 	l_timeCourse_TIME_SCALE = 1;
 	l_timeCourse_VOLT_SCALE = 1;
 	_prop->param = _p;
 	_prop->param_size = 44;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
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

 void _cal_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", 2.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 44, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cal /home/jovyan/work/NeuroML_Examples/Rocher/LEMSFILES5/x86_64/cal.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Mod file for component: Component(id=cal type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dk_q = rate_k_q ;
   Dl_q = rate_l_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dk_q = Dk_q  / (1. - dt*( 0.0 )) ;
 Dl_q = Dl_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    k_q = k_q - dt*(- ( rate_k_q ) ) ;
    l_q = l_q - dt*(- ( rate_l_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
   k_steadyState_x = k_steadyState_rate / ( 1.0 + exp ( 0.0 - ( v - k_steadyState_midpoint ) / k_steadyState_scale ) ) ;
   k_timeCourse_t = k_timeCourse_tau ;
   k_rateScale = 1.0 ;
   k_fcond = pow( k_q , k_instances ) ;
   k_inf = k_steadyState_x ;
   k_tauUnscaled = k_timeCourse_t ;
   k_tau = k_tauUnscaled / k_rateScale ;
   l_steadyState_x = l_steadyState_rate / ( 1.0 + exp ( 0.0 - ( v - l_steadyState_midpoint ) / l_steadyState_scale ) ) ;
   l_timeCourse_V = v / l_timeCourse_VOLT_SCALE ;
   l_timeCourse_t = l_timeCourse_TIME_SCALE * ( ( 20.0 * exp ( ( l_timeCourse_V + 160.0 ) / 30.0 ) / ( 1.0 + exp ( ( l_timeCourse_V + 84.0 ) / 7.3 ) ) ) + 35.0 ) ;
   l_rateScale = 1.0 ;
   l_fcond = pow( l_q , l_instances ) ;
   l_inf = l_steadyState_x ;
   l_tauUnscaled = l_timeCourse_t ;
   l_tau = l_tauUnscaled / l_rateScale ;
   rate_k_q = ( k_inf - k_q ) / k_tau ;
   rate_l_q = ( l_inf - l_q ) / l_tau ;
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
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  k_q = k_q0;
  l_q = l_q0;
 {
   eca = 120.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   k_q = k_inf ;
   l_q = l_inf ;
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
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   conductanceScale = 1.0 ;
   fopen0 = k_fcond * l_fcond ;
   fopen = conductanceScale * fopen0 ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ica = gion * ( v - eca ) ;
   }
 _current += ica;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(k_q) - _p;  _dlist1[0] = &(Dk_q) - _p;
 _slist1[1] = &(l_q) - _p;  _dlist1[1] = &(Dl_q) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/jovyan/work/NeuroML_Examples/Rocher/LEMSFILES5/cal.mod";
static const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=cal type=ionChannelHH)\n"
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
  "    SUFFIX cal\n"
  "    USEION ca WRITE ica VALENCE 2 ? Assuming valence = 2 (Ca ion); TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE k_instances                       : parameter\n"
  "    \n"
  "    RANGE k_tau                             : exposure\n"
  "    \n"
  "    RANGE k_inf                             : exposure\n"
  "    \n"
  "    RANGE k_rateScale                       : exposure\n"
  "    \n"
  "    RANGE k_fcond                           : exposure\n"
  "    RANGE k_steadyState_rate                : parameter\n"
  "    RANGE k_steadyState_midpoint            : parameter\n"
  "    RANGE k_steadyState_scale               : parameter\n"
  "    \n"
  "    RANGE k_steadyState_x                   : exposure\n"
  "    RANGE k_timeCourse_tau                  : parameter\n"
  "    \n"
  "    RANGE k_timeCourse_t                    : exposure\n"
  "    RANGE l_instances                       : parameter\n"
  "    \n"
  "    RANGE l_tau                             : exposure\n"
  "    \n"
  "    RANGE l_inf                             : exposure\n"
  "    \n"
  "    RANGE l_rateScale                       : exposure\n"
  "    \n"
  "    RANGE l_fcond                           : exposure\n"
  "    RANGE l_steadyState_rate                : parameter\n"
  "    RANGE l_steadyState_midpoint            : parameter\n"
  "    RANGE l_steadyState_scale               : parameter\n"
  "    \n"
  "    RANGE l_steadyState_x                   : exposure\n"
  "    RANGE l_timeCourse_TIME_SCALE           : parameter\n"
  "    RANGE l_timeCourse_VOLT_SCALE           : parameter\n"
  "    \n"
  "    RANGE l_timeCourse_t                    : exposure\n"
  "    RANGE k_tauUnscaled                     : derived variable\n"
  "    RANGE l_timeCourse_V                    : derived variable\n"
  "    RANGE l_tauUnscaled                     : derived variable\n"
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
  "    k_instances = 3 \n"
  "    k_steadyState_rate = 1 \n"
  "    k_steadyState_midpoint = -61 (mV)\n"
  "    k_steadyState_scale = 4.2 (mV)\n"
  "    k_timeCourse_tau = 1 (ms)\n"
  "    l_instances = 1 \n"
  "    l_steadyState_rate = 1 \n"
  "    l_steadyState_midpoint = -85.5 (mV)\n"
  "    l_steadyState_scale = -8.5 (mV)\n"
  "    l_timeCourse_TIME_SCALE = 1 (ms)\n"
  "    l_timeCourse_VOLT_SCALE = 1 (mV)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    eca (mV)\n"
  "    ica (mA/cm2)\n"
  "    \n"
  "    \n"
  "    k_steadyState_x                        : derived variable\n"
  "    \n"
  "    k_timeCourse_t (ms)                    : derived variable\n"
  "    \n"
  "    k_rateScale                            : derived variable\n"
  "    \n"
  "    k_fcond                                : derived variable\n"
  "    \n"
  "    k_inf                                  : derived variable\n"
  "    \n"
  "    k_tauUnscaled (ms)                     : derived variable\n"
  "    \n"
  "    k_tau (ms)                             : derived variable\n"
  "    \n"
  "    l_steadyState_x                        : derived variable\n"
  "    \n"
  "    l_timeCourse_V                         : derived variable\n"
  "    \n"
  "    l_timeCourse_t (ms)                    : derived variable\n"
  "    \n"
  "    l_rateScale                            : derived variable\n"
  "    \n"
  "    l_fcond                                : derived variable\n"
  "    \n"
  "    l_inf                                  : derived variable\n"
  "    \n"
  "    l_tauUnscaled (ms)                     : derived variable\n"
  "    \n"
  "    l_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopen0                                 : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_k_q (/ms)\n"
  "    rate_l_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    k_q  \n"
  "    l_q  \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    eca = 120.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    k_q = k_inf\n"
  "    \n"
  "    l_q = l_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=cal type=ionChannelHH), from conductanceScaling; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    conductanceScale = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=cal type=ionChannelHH), from gates; Component(id=k type=gateHHtauInf)\n"
  "    ? multiply applied to all instances of fcond in: <gates> ([Component(id=k type=gateHHtauInf), Component(id=l type=gateHHtauInf)]))\n"
  "    fopen0 = k_fcond * l_fcond ? path based, prefix = \n"
  "    \n"
  "    fopen = conductanceScale  *  fopen0 ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ica = gion * (v - eca)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    k_q' = rate_k_q \n"
  "    l_q' = rate_l_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    \n"
  "    k_steadyState_x = k_steadyState_rate  / (1 + exp(0 - (v -  k_steadyState_midpoint )/ k_steadyState_scale )) ? evaluable\n"
  "    k_timeCourse_t = k_timeCourse_tau ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=k type=gateHHtauInf), from q10Settings; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    k_rateScale = 1 \n"
  "    \n"
  "    k_fcond = k_q ^ k_instances ? evaluable\n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=k type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)\n"
  "    k_inf = k_steadyState_x ? path based, prefix = k_\n"
  "    \n"
  "    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=k type=gateHHtauInf), from timeCourse; Component(id=null type=fixedTimeCourse)\n"
  "    k_tauUnscaled = k_timeCourse_t ? path based, prefix = k_\n"
  "    \n"
  "    k_tau = k_tauUnscaled  /  k_rateScale ? evaluable\n"
  "    l_steadyState_x = l_steadyState_rate  / (1 + exp(0 - (v -  l_steadyState_midpoint )/ l_steadyState_scale )) ? evaluable\n"
  "    l_timeCourse_V = v /  l_timeCourse_VOLT_SCALE ? evaluable\n"
  "    l_timeCourse_t = l_timeCourse_TIME_SCALE *((20 * exp(( l_timeCourse_V  + 160) / 30) / (1 + exp(( l_timeCourse_V  + 84) / 7.3))) +35) ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=l type=gateHHtauInf), from q10Settings; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    l_rateScale = 1 \n"
  "    \n"
  "    l_fcond = l_q ^ l_instances ? evaluable\n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=l type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)\n"
  "    l_inf = l_steadyState_x ? path based, prefix = l_\n"
  "    \n"
  "    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=l type=gateHHtauInf), from timeCourse; Component(id=null type=CaL_tau)\n"
  "    l_tauUnscaled = l_timeCourse_t ? path based, prefix = l_\n"
  "    \n"
  "    l_tau = l_tauUnscaled  /  l_rateScale ? evaluable\n"
  "    \n"
  "     \n"
  "    rate_k_q = ( k_inf  -  k_q ) /  k_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_l_q = ( l_inf  -  l_q ) /  l_tau ? Note units of all quantities used here need to be consistent!\n"
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
