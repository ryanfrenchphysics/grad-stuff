#ifndef _GET_OP_FE_H
#define _GET_OP_FE_H

#include "complex.h"

void parameters_on_FS(int AN, double complex **g_dos);
void FSsetup(int *ANFS, double complex **g_dos);
void domain_OP(int SelfCnst, int AN, double temp, double Dthick,
    double *xgrid, double *zgrid, double **OP11, double **OP22,
    double **OP33, double **OP13, double **OP31, int argc, char *argv[]);
void next(double *z, double znext, double *err, char *error, char *where);
void angle_integrator_OP_fake(double complex en);
void angle_integrator(double complex en, int AN);
void trajectory_integrator(double complex en, int an,
    double **DEL11, double **DEL22, double **DEL33,
    double **DEL13, double **DEL31);
void nextc(double complex *z, double complex znext,
    double *err, char *error, char *where);
void get_bc(double complex en, int nxb, int an,
    double complex **a_bc, double complex **aa_bc);
void OP_carpet_DW(void);

#endif  /* _GET_OP_FE_H */