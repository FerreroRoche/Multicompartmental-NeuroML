/* Created by Language version: 6.2.0 */
/* VECTORIZED */
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
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
 
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
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
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
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
  hoc_register_dparam_size(_mechtype, 3);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cal /cygdrive/c/Users/roche/Thesis code nml/Cells/cal.mod\n");
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
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
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
 }}

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
 
}}

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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
 double _break, _save;
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
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
 _break = t + .5*dt; _save = t;
 v=_v;
{
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
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
