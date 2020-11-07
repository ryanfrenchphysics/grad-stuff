#ifndef _OP_UPDATE_H
#define _OP_UPDATE_H
#include <complex.h>

void OP_update(int Nvector, double complex *vec_next, 
	double complex *Fvec_next, double complex *vector[], 
	double complex *Fvec[], double complex *Fmin_prev, double Norm[]);
void ludcmp_z(double complex **a, int z_n, int *indx, double complex *det);
void lubksb_z(double complex **a, int z_n, int *indx, double complex z_b[]);

#endif /* _OP_UPDATE_H */