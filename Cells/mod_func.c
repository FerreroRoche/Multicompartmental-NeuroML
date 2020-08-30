#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ca_conc_reg();
extern void _cacc_reg();
extern void _cah_reg();
extern void _cal_reg();
extern void _h_reg();
extern void _iClamp0_reg();
extern void _k_reg();
extern void _kca_reg();
extern void _kdr_reg();
extern void _leak_reg();
extern void _na_a_reg();
extern void _na_s_reg();
extern void _vClamp0_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ca_conc.mod");
fprintf(stderr," cacc.mod");
fprintf(stderr," cah.mod");
fprintf(stderr," cal.mod");
fprintf(stderr," h.mod");
fprintf(stderr," iClamp0.mod");
fprintf(stderr," k.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," kdr.mod");
fprintf(stderr," leak.mod");
fprintf(stderr," na_a.mod");
fprintf(stderr," na_s.mod");
fprintf(stderr," vClamp0.mod");
fprintf(stderr, "\n");
    }
_ca_conc_reg();
_cacc_reg();
_cah_reg();
_cal_reg();
_h_reg();
_iClamp0_reg();
_k_reg();
_kca_reg();
_kdr_reg();
_leak_reg();
_na_a_reg();
_na_s_reg();
_vClamp0_reg();
}
