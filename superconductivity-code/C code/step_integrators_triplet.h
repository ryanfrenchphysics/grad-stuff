#ifndef _STEP_INTEGRATORS_TRIPLET_H
#define _STEP_INTEGRATORS_TRIPLET_H

#include <complex.h>

/*
 * Step integrators for Riccati equations
 * 
 * Updated July 2013 - ABV
 * Updated Jan 2019 - Ryan French
 */

typedef double complex matrix2x2[2][2];

void step_integratorc_a(double ds, double complex tm, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);
void step_integratorc_b(double ds, double complex tm, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);


/* 
* These step integrators are valid for NON-UNITARY phases with general Delta[0,1,2,3]
* and scalar self-energy that is included in en
* 		en = I * t(m + 0.5) - \nu  or en = energy + I * 0 - \nu 
* 		( but no magnetic field or vector self-energy )
*/

void step_integratornu_a(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);
void step_integratornu_b(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);



void step_integratorm_a(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);
void step_integratorm_b(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf);


void cmatrixsum(int sgn, matrix2x2 A, matrix2x2 B, matrix2x2 C);
void cmatrixmult(matrix2x2 A, matrix2x2 B, matrix2x2 C);
void c2x2matrixinv(matrix2x2 C, matrix2x2 C_inv);


#endif  /* _STEP_INTEGRATORS_TRIPLET_H */