#ifndef _OP_GUESS_H
#define _OP_GUESS_H

#include <complex.h>

void op_guess(double dfilm, double t, double *xgr, double *zgr, 
	double **A11, double **A22, double **A33, double **A13, double **A31,
    int rescale_op, int read_op_film);

void gapBfilm_specular(double dfilm, double z_tt, double complex **Delta,
    double *zgr); 

#endif /* _OP_GUESS_H */