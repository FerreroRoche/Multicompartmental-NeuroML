#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _BK_reg(void);
extern void _cacc_reg(void);
extern void _ca_conc_reg(void);
extern void _cah_reg(void);
extern void _cal_reg(void);
extern void _h_reg(void);
extern void _iclamp1_reg(void);
extern void _kca_reg(void);
extern void _kdr_reg(void);
extern void _k_reg(void);
extern void _leak_reg(void);
extern void _na_a_reg(void);
extern void _na_s_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," BK.mod");
    fprintf(stderr," cacc.mod");
    fprintf(stderr," ca_conc.mod");
    fprintf(stderr," cah.mod");
    fprintf(stderr," cal.mod");
    fprintf(stderr," h.mod");
    fprintf(stderr," iclamp1.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kdr.mod");
    fprintf(stderr," k.mod");
    fprintf(stderr," leak.mod");
    fprintf(stderr," na_a.mod");
    fprintf(stderr," na_s.mod");
    fprintf(stderr, "\n");
  }
  _BK_reg();
  _cacc_reg();
  _ca_conc_reg();
  _cah_reg();
  _cal_reg();
  _h_reg();
  _iclamp1_reg();
  _kca_reg();
  _kdr_reg();
  _k_reg();
  _leak_reg();
  _na_a_reg();
  _na_s_reg();
}
