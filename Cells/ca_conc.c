/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
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
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 
#define _threadargsprotocomma_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define restingConc _p[0]
#define decayConstant _p[1]
#define rho _p[2]
#define concentration _p[3]
#define extConcentration _p[4]
#define cai _p[5]
#define cao _p[6]
#define ica _p[7]
#define rate_concentration _p[8]
#define Dconcentration _p[9]
#define DextConcentration _p[10]
#define _g _p[11]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _style_ca	*((int*)_ppvar[3]._pvoid)
#define diam	*_ppvar[4]._pval
#define area	*_ppvar[5]._pval
 
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
 /* external NEURON variables */
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
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ca_conc", _hoc_setdata,
 "rates_ca_conc", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define iCa iCa_ca_conc
 double iCa = 0;
#define initialExtConcentration initialExtConcentration_ca_conc
 double initialExtConcentration = 0;
#define initialConcentration initialConcentration_ca_conc
 double initialConcentration = 0;
#define surfaceArea surfaceArea_ca_conc
 double surfaceArea = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "surfaceArea_ca_conc", "um2",
 "iCa_ca_conc", "nA",
 "initialConcentration_ca_conc", "mM",
 "initialExtConcentration_ca_conc", "mM",
 "restingConc_ca_conc", "mM",
 "decayConstant_ca_conc", "ms",
 "rho_ca_conc", "mM",
 "concentration_ca_conc", "mM",
 "extConcentration_ca_conc", "mM",
 0,0
};
 static double concentration0 = 0;
 static double delta_t = 0.01;
 static double extConcentration0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "surfaceArea_ca_conc", &surfaceArea_ca_conc,
 "iCa_ca_conc", &iCa_ca_conc,
 "initialConcentration_ca_conc", &initialConcentration_ca_conc,
 "initialExtConcentration_ca_conc", &initialExtConcentration_ca_conc,
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
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_synonym(int, double**, Datum**);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"ca_conc",
 "restingConc_ca_conc",
 "decayConstant_ca_conc",
 "rho_ca_conc",
 0,
 0,
 "concentration_ca_conc",
 "extConcentration_ca_conc",
 0,
 0};
 static Symbol* _morphology_sym;
 extern Node* nrn_alloc_node_;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	restingConc = 0;
 	decayConstant = 13.3333;
 	rho = 300000;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[4]._pval = &prop_ion->param[0]; /* diam */
 	_ppvar[5]._pval = &nrn_alloc_node_->_area; /* diam */
 prop_ion = need_memb(_ca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for ca */
 
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

 void _ca_conc_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", 2.0);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 7);
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ca_conc /cygdrive/c/Users/roche/Thesis code nml/Cells/ca_conc.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Mod file for component: Component(id=ca_conc type=fixedFactorConcentrationModel)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates();
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargs_ ) ;
   Dconcentration = rate_concentration ;
   cai = concentration ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargs_ ) ;
 Dconcentration = Dconcentration  / (1. - dt*( 0.0 )) ;
 return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargs_ ) ;
    concentration = concentration - dt*(- ( rate_concentration ) ) ;
   cai = concentration ;
   }
  return 0;
}
 
static int  rates (  ) {
   surfaceArea = area ;
   iCa = - 1.0 * ( 0.01 ) * ica * surfaceArea ;
   rate_concentration = ( iCa / surfaceArea ) * rho - ( ( concentration - restingConc ) / decayConstant ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
  ica = _ion_ica;
  cai = _ion_cai;
     _ode_spec1 ();
  _ion_cai = cai;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 static void _ode_synonym(int _cnt, double** _pp, Datum** _ppd) { 
 	int _i; 
	for (_i=0; _i < _cnt; ++_i) {_p = _pp[_i]; _ppvar = _ppd[_i];
 _ion_cai =  concentration ;
 }}
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
  ica = _ion_ica;
  cai = _ion_cai;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  concentration = concentration0;
  extConcentration = extConcentration0;
 {
   initialConcentration = cai ;
   initialExtConcentration = cao ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   concentration = initialConcentration ;
   extConcentration = initialExtConcentration ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  ica = _ion_ica;
  cai = _ion_cai;
 initmodel();
  _ion_cai = cai;
  nrn_wrote_conc(_ca_sym, (&(_ion_cai)) - 1, _style_ca);
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  cai = _ion_cai;
  cao = _ion_cao;
  ica = _ion_ica;
  cai = _ion_cai;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 85 in file ca_conc.mod:\n    \n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } {
   if ( concentration < 0.0 ) {
     concentration = 0.0 ;
     }
   }
  _ion_cai = cai;
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(concentration) - _p;  _dlist1[0] = &(Dconcentration) - _p;
_first = 0;
}
