#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrutil.h"
#include "ricparam.h"
#include "step_integrators_triplet.h"

/* Coherence lengths - basically uniform along this direction */
#define MAXds       1e+5


/*
 * step_integratorc_a(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for a-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * This integrator is valid only for UNITARY phases with
 * \vDelta \times \vDelta^* = 0.
 * 
 * Arguments:
 *      double ds
 *      double complex tm       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 * 
 */

void step_integratorc_a(double ds, double complex tm, double complex *Deltai, 
		    double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j, sgn_tm;
	double complex c0[4], Delta[4], delta_c[4];
	double complex ee, eexp, Deltac, x, A, dd, Delta2; 
	double complex denom;

	sgn_tm = (creal(tm) >= 0.0) ? 1 : -1;
	
	if (ds > MAXds) {
        /* Homogeneous solution along this direction */
		denom = 0.0;
		for (j = 1; j <= 3; j++) {
			Delta[j] = Deltaf[j];
			denom += pow(cabs(Delta[j]), 2);
		}

		denom = tm + I * csqrt(denom - tm * tm);
		cf[0] = ci[0];
		for (j = 1; j <= 3; j++)
            cf[j] = c0[j] = -Delta[j] / denom;
		return;
	}

	denom = 0.0;
	for (j = 1; j <= 3; j++) {
		/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
		Delta[j] = Deltaf[j];
		denom += pow(cabs(Delta[j]), 2);
	}

	denom = tm + I * csqrt(denom - tm * tm);

	ee = tm;
	Deltac = Delta2 = dd = 0.0;

	for (j = 1; j <= 3; j++) {
		c0[j] = -Delta[j] / denom;
		ee += c0[j] * conj(Delta[j]);
		delta_c[j] = ci[j] - c0[j];
		dd += delta_c[j] * delta_c[j];
		Delta2 += conj(Delta[j]) * conj(Delta[j]);
		Deltac += delta_c[j] * conj(Delta[j]);
	}

	eexp = cexp(2.0 * I * ee * ds);
	A = (1.0 - eexp) / 2.0 / ee;
	x = 1.0 + A * Deltac;
	denom = x * x + A * A * (Delta2 * dd - Deltac * Deltac);

	cf[0] = ci[0];

	for (j = 1; j <= 3; j++) {
		cf[j] = c0[j] + (delta_c[j] + A * conj(Delta[j]) * dd) * eexp / denom;

		if (isnan(cabs(cf[j])) != 0.0) {
			printf("NAN in a-Riccati! Exiting.\n");
			printf("Deltai[%d]=%e+I(%e)\n", j, creal(Deltai[j]),
                cimag(Deltai[j]));
			printf("Deltaf[%d]=%e+I(%e)\n", j, creal(Deltaf[j]),
                cimag(Deltaf[j]));
			printf("delta_c[%d]=%e+I(%e)\n", j, creal(delta_c[j]),
                cimag(delta_c[j]));
			printf("ds = %f \n", ds);
			printf("ee = %e+I(%e) \n", creal(ee), cimag(ee));
			printf("dd = %e+I(%e) \n", creal(dd), cimag(dd));
			printf("eexp = %e+I(%e) \n", creal(eexp), cimag(eexp));
			printf("A = %e+I(%e) \n", creal(A), cimag(A));
			printf("x = %e+I(%e) \n", creal(x), cimag(x));
			printf("denom = %e+I(%e) \n", creal(denom), cimag(denom));
			exit(1);
		}
	}
}


/*
 * step_integratorc_b(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for b-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * This integrator is valid only for UNITARY phases with
 * \vDelta \times \vDelta^* = 0.
 * 
 * Arguments:
 *      double ds
 *      double complex tm       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 * 
 */

void step_integratorc_b(double ds, double complex tm, double complex *Deltai, 
		    double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j, sgn_tm;
	double complex c0[4], Delta[4], delta_c[4];
	double complex ee, eexp, dd, x, B, Deltac, Delta2; 
	double complex denom;

	sgn_tm = (creal(tm) >= 0.0) ? 1 : -1;
    denom = 0.0;
	
	if (ds > MAXds) {
        /* Homogeneous solution along this direction */

		for (j = 1; j <= 3; j++) {
			Delta[j] = Deltaf[j];
			denom += pow(cabs(Delta[j]), 2);
		}

		denom = tm + I * csqrt(denom - tm * tm);
		cf[0] = ci[0];

		for (j = 1; j <= 3; j++)
            cf[j] = c0[j] = conj(Delta[j]) / denom;
		return;
	}

	for (j = 1; j <= 3; j++) {
		/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
		Delta[j] = Deltaf[j];
		denom += pow(cabs(Delta[j]), 2);
	}

	denom = tm + I * csqrt(denom - tm * tm);

	ee = tm;
	Deltac = Delta2 = dd = 0.0;

	for (j = 1; j <= 3; j++) {
		c0[j] = conj(Delta[j]) / denom;
		ee -= c0[j] * Delta[j];
		delta_c[j] = ci[j] - c0[j];
		Delta2 += Delta[j] * Delta[j];
		Deltac += delta_c[j] * Delta[j];
		dd += delta_c[j] * delta_c[j];
	}

	eexp = cexp(2.0 * I * ee * ds);
	B = (1.0 - eexp) / 2 / ee;
	x = 1.0 - B * Deltac;
	denom = x * x + B * B * (Delta2 * dd - Deltac * Deltac);

	cf[0] = ci[0];

	for (j = 1; j <= 3; j++) {
		cf[j] = c0[j] + (delta_c[j] - B * Delta[j] * dd) * eexp / denom;

		if (isnan(cabs(cf[j])) != 0.0) {
			printf("NAN in b-Riccati! Exiting.\n");
			exit(1);
		}
	}
}


/*
 * step_integratornu_a(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for a-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * Arguments:
 *      double ds
 *      double complex en       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 * 
 */

void step_integratornu_a(double ds, double complex en, double complex *Deltai, 
		    double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j;
	double complex sqvec; 
	double complex a_h[4], dum[4], delc[4], q[4], sqrtq[4], ex, ch, sh;
	double complex Delta[4], DeltaC[4], aux_g[4], aux_h[4], aux_f[4]; 

	if (ds > MAXds) {
		for (j = 0; j <= 3; j++) {
			Delta[j] = Deltaf[j];
		}
		ab_bulk_non_unitary(a_h, dum, en, Delta);
		for (j = 0; j <= 3; j++)
            cf[j] = a_h[j];

		return;
	} else {
		for (j = 0; j <= 3; j++) {
			/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
			Delta[j] = Deltaf[j];
			DeltaC[j] = conj(Delta[j]);
		}
		ab_bulk_non_unitary(a_h, dum, en, Delta);
		for (j = 0; j <= 3; j++)
            delc[j] = ci[j] - a_h[j];
	}

	ric_multiply(Delta, DeltaC, q);
	q[0] -= en * en; 
	ric_sqrt(q, sqrtq); 
	sqvec = csqrt(sqrtq[1] * sqrtq[1] + sqrtq[2] * sqrtq[2] +
        sqrtq[3] * sqrtq[3]);

	if (cabs(sqvec) == 0)
        sqvec = 1e-10;      /* Avoid NAN */

	ex = cexp(-sqrtq[0] * ds);
	ch = ccosh(sqvec * ds); 
	sh = -csinh(sqvec * ds);

	aux_g[0] = ex * ch; 
	aux_g[1] = ex * sh * sqrtq[1] / sqvec; 
	aux_g[2] = ex * sh * sqrtq[2] / sqvec; 
	aux_g[3] = ex * sh * sqrtq[3] / sqvec; 

	ex = cexp(-2 * sqrtq[0] * ds);
	ch = ccosh(2 * sqvec * ds); 
	sh = -csinh(2 * sqvec * ds); 
	aux_f[0] = -I * 0.5 * (ex * ch - 1); 
	aux_f[1] = -I * 0.5 * ex * sh * sqrtq[1] / sqvec; 
	aux_f[2] = -I * 0.5 * ex * sh * sqrtq[2] / sqvec; 
	aux_f[3] = -I * 0.5 * ex * sh * sqrtq[3] / sqvec;

	ric_invert(sqrtq, dum); 
	ric_multiply(dum, aux_f, aux_f);
	ric_multiply(DeltaC, aux_f, aux_f);

	if (cabs(Delta[0]) == 0) {
		aux_h[0] = aux_g[0]; 
		aux_h[1] = -aux_g[1]; 
		aux_h[2] = -aux_g[2]; 
		aux_h[3] = -aux_g[3]; 
	} else {
		ric_multiply(DeltaC, Delta, q);
		q[0] -= en * en; 
		ric_sqrt(q, sqrtq); 
		sqvec = csqrt(sqrtq[1] * sqrtq[1] + sqrtq[2] * sqrtq[2] +
            sqrtq[3] * sqrtq[3]);

		if (cabs(sqvec) == 0)
            sqvec = 1e-10;      /* Avoid NAN */

		ex = cexp(-sqrtq[0] * ds);
		ch = ccosh(sqvec * ds); 
		sh = -csinh(sqvec * ds);

		aux_h[0] = ex * ch; 
		aux_h[1] = ex * sh * sqrtq[1] / sqvec; 
		aux_h[2] = ex * sh * sqrtq[2] / sqvec; 
		aux_h[3] = ex * sh * sqrtq[3] / sqvec; 
	}

	ric_multiply(aux_f, delc, dum);

	dum[0] = 1 - dum[0];
    dum[1] = -dum[1];
    dum[2] = -dum[2];
    dum[3] = -dum[3];

	ric_invert(dum, cf); 
	ric_multiply(cf, aux_h, cf);
	ric_multiply(delc, cf, cf);
	ric_multiply(aux_g, cf, cf);

	for (j = 0;j <= 3;j++)
        cf[j] += a_h[j];

	for (j = 0; j <= 3; j++) {
		if (isnan(cabs(cf[j])) != 0.0) {
			printf("\n !!! NAN in a-Riccati! Exiting.\n");
			printf("Delta[1]=%f Delta[2]=%f Delta[3]=%f \n", cabs(Delta[1]),
                cabs(Delta[2]), cabs(Delta[3])); 
			printf("en=%f+I(%f)\n", creal(en), cimag(en)); 
			printf("q[0]=%f+I(%f)\n", creal(q[0]), cimag(q[0])); 
			printf("q[1]=%f+I(%f)\n", creal(q[1]), cimag(q[1])); 
			printf("q[2]=%f+I(%f)\n", creal(q[2]), cimag(q[2])); 
			printf("q[3]=%f+I(%f)\n", creal(q[3]), cimag(q[3])); 
			printf("sqrtq[0]=%f+I(%f)\n", creal(sqrtq[0]), cimag(sqrtq[0])); 
			printf("sqrtq[1]=%f+I(%f)\n", creal(sqrtq[1]), cimag(sqrtq[1])); 
			printf("sqrtq[2]=%f+I(%f)\n", creal(sqrtq[2]), cimag(sqrtq[2])); 
			printf("sqrtq[3]=%f+I(%f)\n", creal(sqrtq[3]), cimag(sqrtq[3])); 
			printf("sqvec=%f+I(%f)\n", creal(sqvec), cimag(sqvec)); 
			printf("ex=%f ch=%f sh=%f \n", cabs(ex), cabs(ch), cabs(sh)); 
			exit(1);
		}
	}
}


/*
 * step_integratornu_b(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for b-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * Arguments:
 *      double ds               ->
 *      double complex en       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 * 
 */

void step_integratornu_b(double ds, double complex en, double complex *Deltai, 
		    double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j;
	double complex sqvec; 
	double complex a_h[4], dum[4], delc[4], q[4], sqrtq[4], ex, ch, sh;
	double complex Delta[4], DeltaC[4], aux_g[4], aux_h[4], aux_f[4]; 

	if (ds > MAXds) {
		for (j = 0; j <= 3; j++)
			Delta[j] = Deltaf[j];

		ab_bulk_non_unitary(dum, a_h, en, Delta);

		for (j = 0; j <= 3; j++)
            cf[j] = a_h[j];

		return;
	} else {
		for (j = 0; j <= 3; j++) {
			/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
			Delta[j] = Deltaf[j];
			DeltaC[j] = conj(Delta[j]);
		}
		ab_bulk_non_unitary(dum, a_h, en, Delta);
		for (j = 0; j <= 3; j++)
            delc[j] = ci[j] - a_h[j];
	}

	ric_multiply(Delta, DeltaC, q);
	q[0] -= en * en; 
	ric_sqrt(q, sqrtq); 
	sqvec = csqrt(sqrtq[1] * sqrtq[1] + sqrtq[2] * sqrtq[2] +
        sqrtq[3] * sqrtq[3]);

	if (cabs(sqvec) == 0)
        sqvec = 1e-10;          /* Avoid NAN */

	ex = cexp(-sqrtq[0] * ds);
	ch = ccosh(sqvec * ds); 
	sh = -csinh(sqvec * ds);

	aux_h[0] = ex * ch; 
	aux_h[1] = ex * sh * sqrtq[1] / sqvec; 
	aux_h[2] = ex * sh * sqrtq[2] / sqvec; 
	aux_h[3] = ex * sh * sqrtq[3] / sqvec; 

	ex = cexp(-2 * sqrtq[0] * ds);
	ch = ccosh(2 * sqvec * ds); 
	sh = - csinh(2 * sqvec * ds);

	aux_f[0] = I * 0.5 * (ex * ch - 1); 
	aux_f[1] = I * 0.5 * ex * sh * sqrtq[1] / sqvec; 
	aux_f[2] = I * 0.5 * ex * sh * sqrtq[2] / sqvec; 
	aux_f[3] = I * 0.5 * ex * sh * sqrtq[3] / sqvec;

	ric_invert(sqrtq, dum); 
	ric_multiply(dum, aux_f, aux_f);
	ric_multiply(aux_f, Delta, aux_f);

	if (cabs(Delta[0]) == 0) {
		aux_g[0] = aux_h[0]; 
		aux_g[1] = -aux_h[1]; 
		aux_g[2] = -aux_h[2]; 
		aux_g[3] = -aux_h[3]; 
	} else {
		ric_multiply(DeltaC, Delta, q);
		q[0] -= en * en; 
		ric_sqrt(q, sqrtq); 
		sqvec = csqrt(sqrtq[1] * sqrtq[1] + sqrtq[2] * sqrtq[2] +
            sqrtq[3] * sqrtq[3]);

		if (cabs(sqvec) == 0)
            sqvec = 1e-10;          /* Avoid NAN */

		ex = cexp(-sqrtq[0] * ds);
		ch = ccosh(sqvec * ds); 
		sh = -csinh(sqvec * ds);

		aux_g[0] = ex * ch; 
		aux_g[1] = ex * sh * sqrtq[1] / sqvec; 
		aux_g[2] = ex * sh * sqrtq[2] / sqvec; 
		aux_g[3] = ex * sh * sqrtq[3] / sqvec; 
	}

	ric_multiply(aux_f, delc, dum);

	dum[0] = 1.0 - dum[0];
    dum[1] = -dum[1];
    dum[2] = -dum[2];
    dum[3] = -dum[3];

	ric_invert(dum, cf); 
	ric_multiply(cf, aux_h, cf);
	ric_multiply(delc, cf, cf);
	ric_multiply(aux_g, cf, cf);

	for (j = 0; j <= 3; j++)
        cf[j] += a_h[j]; 

	for (j = 0; j <= 3; j++) {
		if (isnan(cabs(cf[j])) != 0.0) {
			printf("\n !!! NAN in b-Riccati! Exiting.\n");
			printf("Delta[1]=%f Delta[2]=%f Delta[3]=%f \n", cabs(Delta[1]),
                cabs(Delta[2]), cabs(Delta[3])); 
			printf("en=%f+I(%f)\n", creal(en), cimag(en)); 
			printf("q[0]=%f+I(%f)\n", creal(q[0]), cimag(q[0])); 
			printf("q[1]=%f+I(%f)\n", creal(q[1]), cimag(q[1])); 
			printf("q[2]=%f+I(%f)\n", creal(q[2]), cimag(q[2])); 
			printf("q[3]=%f+I(%f)\n", creal(q[3]), cimag(q[3])); 
			printf("sqrtq[0]=%f+I(%f)\n", creal(sqrtq[0]), cimag(sqrtq[0])); 
			printf("sqrtq[1]=%f+I(%f)\n", creal(sqrtq[1]), cimag(sqrtq[1])); 
			printf("sqrtq[2]=%f+I(%f)\n", creal(sqrtq[2]), cimag(sqrtq[2])); 
			printf("sqrtq[3]=%f+I(%f)\n", creal(sqrtq[3]), cimag(sqrtq[3])); 
			printf("sqvec=%f+I(%f)\n", creal(sqvec), cimag(sqvec)); 
			printf("ex=%f ch=%f sh=%f \n", cabs(ex), cabs(ch), cabs(sh)); 
			exit(1);
		}
	}
}


/*
 * step_integratorm_a(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for a-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * This integrator is valid for NON-UNITARY phases with 
 * 		\vq = \vDelta\times\vDelta^* != 0  (No scalar Delta[0] component)
 * and scalar self energy that is included in en
 * 		en=I*t(m+0.5)-\nu  or en=energy+I*0-\nu 
 * Still no magnetic field or vector self-energy, but maybe this is the way 
 * to implement them LATER.
 * 
 * This integrator is slower than the above  Non-Unitary integrator by
 * 1.5-2 times.
 * 
 * Arguments:
 *      double ds               ->
 *      double complex en       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 */

void step_integratorm_a(double ds, double complex en, double complex *Deltai, 
		    double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j;
	double ampl, phi, Delta2;
	double complex E, q[4], q2, alpha, Norm; 
	double complex g_h[4], f_h[4], ff_h[4], a_h[4], f2, ff2, f_ff;
	double complex Delta[4], m[4], Delta_a, denom, lambda;
	double complex ee, eexp, eexp2, f[4]; 
	matrix2x2 S, S_inv, SS, SS_inv;
	matrix2x2 LAMBDA, UNITY, G0, H0, F0;
	matrix2x2 DELTA, delta, af;

	UNITY[0][0] = UNITY[1][1] = 1.0; 
	UNITY[0][1] = UNITY[1][0] = 0.0;

	Delta2 = 0.0;
	if (ds > MAXds) {
		for (j = 1; j <= 3; j++) {
			Delta[j] = Deltaf[j];
			Delta2 += pow(cabs(Delta[j]), 2);
		}
	} else {
		for (j = 1; j <= 3; j++) {
			/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
			Delta[j] = Deltaf[j];
			Delta2 += pow(cabs(Delta[j]), 2);
		}
	}

	E = Delta2 - en * en;
	q[1] = -2 * I * (creal(Delta[2]) * cimag(Delta[3]) -
        cimag(Delta[2]) * creal(Delta[3]));
	q[2] = -2 * I * (-creal(Delta[1]) * cimag(Delta[3]) +
        cimag(Delta[1]) * creal(Delta[3]));
	q[3] = -2 * I * (creal(Delta[1]) * cimag(Delta[2]) -
        cimag(Delta[1]) * creal(Delta[2]));
	q2 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	alpha = 1.0 / (E + csqrt(E * E + q2));
	Norm = -1.0 / csqrt(2.0 * alpha) / csqrt(E * E + q2); 

	/* Bulk (homogeneous) propagator */
	g_h[0] = Norm * en;
	g_h[1] = -I * Norm * alpha * en * q[1];
	g_h[2] = -I * Norm * alpha * en * q[2];
	g_h[3] = -I * Norm * alpha * en * q[3];

	f_h[0] = ff_h[0] = 0.0;
	f_h[1] = -Norm * (Delta[1] - alpha * (Delta[2] * q[3] - Delta[3] * q[2]));
	f_h[2] = -Norm * (Delta[2] - alpha * (-Delta[1] * q[3] + Delta[3] * q[1]));
	f_h[3] = -Norm * (Delta[3] - alpha * (Delta[1] * q[2] - Delta[2] * q[1]));

	ff_h[1] = -Norm * (conj(Delta[1]) + alpha * (conj(Delta[2]) * q[3] -
        conj(Delta[3]) * q[2]));
	ff_h[2] = -Norm * (conj(Delta[2]) + alpha * (-conj(Delta[1]) * q[3] +
        conj(Delta[3]) * q[1]));
	ff_h[3] = -Norm * (conj(Delta[3]) + alpha * (conj(Delta[1]) * q[2] -
        conj(Delta[2]) * q[1]));

	f2 = f_h[1] * f_h[1] + f_h[2] * f_h[2] + f_h[3] * f_h[3];
	ff2 = ff_h[1] * ff_h[1] + ff_h[2] * ff_h[2] + ff_h[3] * ff_h[3];
	f_ff = f_h[1] * ff_h[1] + f_h[2] * ff_h[2] + f_h[3] * ff_h[3];
	denom = (I - g_h[0]) * (I - g_h[0]) +
        (f2 * ff2 - (f_ff * f_ff)) / 4 / g_h[0] / g_h[0];
	
	/* Bulk (homogeneous) Riccati a-amplitude */
	a_h[0] = 0.0;
	a_h[1] = -1.0 / denom * ((I - g_h[0]) * f_h[1] +
        1.0 / 2 / g_h[0] * (f_h[1] * f_ff - ff_h[1] * f2));
	a_h[2] = -1.0 / denom * ((I - g_h[0]) * f_h[2] +
        1.0 / 2 / g_h[0] * (f_h[2] * f_ff - ff_h[2] * f2));
	a_h[3] = -1.0 / denom * ((I - g_h[0]) * f_h[3] +
        1.0 / 2 / g_h[0] * (f_h[3] * f_ff - ff_h[3] * f2));

	if (ds > MAXds) {
		for (j = 0; j <= 3; j++)
            cf[j] = a_h[j];

		return;
	}

	Delta_a = a_h[1] * conj(Delta[1]) + a_h[2] * conj(Delta[2]) +
        a_h[3] * conj(Delta[3]);
	ee = en + Delta_a;
	eexp = cexp(I * ee * ds);
	eexp2 = eexp * eexp;

	/* Diagonalization and so on... */
	m[1] = a_h[2] * conj(Delta[3]) - a_h[3] * conj(Delta[2]);
	m[2] = a_h[3] * conj(Delta[1]) - a_h[1] * conj(Delta[3]);
	m[3] = a_h[1] * conj(Delta[2]) - a_h[2] * conj(Delta[1]);

	if (cabs(m[3]) != 0.0)
		lambda = I * m[3] * csqrt(1.0 +
            cpow(m[1] / m[3], 2) + cpow(m[2] / m[3], 2));
	else
		lambda = I * csqrt(m[1] * m[1] + m[2] * m[2] + m[3] * m[3]);

	if (cabs(lambda) == 0.0) {
		G0[0][0] = eexp * (1.0 - m[3] * ds);
		G0[1][1] = eexp * (1.0 + m[3] * ds);
		G0[0][1] = eexp * (-m[1] + I * m[2]) * ds;
		G0[1][0] = eexp * (-m[1] - I * m[2]) * ds;

		H0[0][0] = eexp * (1.0 - m[3] * ds);
		H0[1][1] = eexp * (1.0 + m[3] * ds);
		H0[0][1] = eexp * (-m[1] - I * m[2]) * ds;
		H0[1][0] = eexp * (-m[1] + I * m[2]) * ds;

		f[1] = (eexp2 - 1.0) / 2 / ee *
			(conj(Delta[1]) + (conj(Delta[2]) * m[3] - conj(Delta[3]) * m[2]) / ee) -
			I * eexp2 / ee * ds * (conj(Delta[2]) * m[3] - conj(Delta[3]) * m[2]);
		f[2] = (eexp2 - 1.0) / 2 / ee *
			(conj(Delta[2]) + (conj(Delta[3]) * m[1] - conj(Delta[1]) * m[3]) / ee) -
			I * eexp2 / ee * ds * (conj(Delta[3]) * m[1] - conj(Delta[1]) * m[3]);
		f[3] = (eexp2 - 1.0) / 2 / ee *
			(conj(Delta[3]) + (conj(Delta[1]) * m[2] - conj(Delta[2]) * m[1]) / ee) -
			I * eexp2 / ee * ds * (conj(Delta[1]) * m[2]-conj(Delta[2]) * m[1]);

		F0[0][0] = f[1] + I * f[2];
		F0[1][1] = -f[1] + I * f[2];
		F0[0][1] = F0[1][0] = -f[3];
	} else {
		LAMBDA[0][0] = lambda;
		LAMBDA[1][1] = -lambda;
		LAMBDA[0][1] = LAMBDA[1][0] = 0.0;

		ampl = sqrt(0.5 * cabs(1 + I * m[3] / lambda));
		phi = -0.5 * carg(2.0 / (1 + I * m[3] / lambda));

		S[0][0] = ampl * cexp(I * phi);
		S[1][1] = ampl * cexp(-I * phi);
		S[0][1] = -S[1][1] * (m[1] - I * m[2]) / (m[3] - I * lambda);
		S[1][0] = S[0][0] * (m[1] + I * m[2]) / (m[3] - I * lambda);

		S_inv[0][0] = S[1][1];
		S_inv[1][1] = S[0][0];
		S_inv[0][1] = -S[0][1];
		S_inv[1][0] = -S[1][0];
		
		SS[0][0] = ampl * cexp(I * phi);
		SS[1][1] = ampl * cexp(-I * phi);
		SS[0][1] = -SS[1][1] * (m[1] + I * m[2]) / (m[3] - I * lambda);
		SS[1][0] = SS[0][0] * (m[1] - I * m[2]) / (m[3] - I * lambda);

		SS_inv[0][0] = SS[1][1];
		SS_inv[1][1] = SS[0][0];
		SS_inv[0][1] = -SS[0][1];
		SS_inv[1][0] = -SS[1][0];

		DELTA[0][0] = conj(Delta[1]) + I * conj(Delta[2]);
		DELTA[1][1] = -conj(Delta[1]) + I * conj(Delta[2]);
		DELTA[0][1] = -conj(Delta[3]);
		DELTA[1][0] = -conj(Delta[3]);

		cmatrixmult(DELTA, S, DELTA);
		cmatrixmult(SS_inv, DELTA, DELTA);

		F0[0][0] = (eexp2 * cexp(I * (LAMBDA[0][0] + LAMBDA[0][0]) * ds) - 1) /
			(2 * ee + LAMBDA[0][0] + LAMBDA[0][0]) * DELTA[0][0];
		F0[1][1] = (eexp2 * cexp(I * (LAMBDA[1][1] + LAMBDA[1][1]) * ds) - 1) /
			(2 * ee + LAMBDA[1][1] + LAMBDA[1][1]) * DELTA[1][1];
		F0[0][1] = (eexp2 * cexp(I * (LAMBDA[0][0] + LAMBDA[1][1]) * ds) - 1) /
			(2 * ee + LAMBDA[0][0] + LAMBDA[1][1]) * DELTA[0][1];
		F0[1][0] = (eexp2 * cexp(I * (LAMBDA[1][1] + LAMBDA[0][0]) * ds) - 1) /
			(2 * ee + LAMBDA[1][1] + LAMBDA[0][0]) * DELTA[1][0];

		cmatrixmult(F0, S_inv, F0);
		cmatrixmult(SS, F0, F0);
		/* found F0, now find G0, H0 */

		G0[0][0] = H0[0][0] = eexp * cexp(I * LAMBDA[0][0] * ds);
		G0[1][1] = H0[1][1] = eexp * cexp(I * LAMBDA[1][1] * ds);
		G0[0][1] = H0[0][1] = 0.0;
		G0[1][0] = H0[1][0] = 0.0;

		cmatrixmult(S, G0, G0);
		cmatrixmult(G0, S_inv, G0);
		cmatrixmult(SS, H0, H0);
		cmatrixmult(H0, SS_inv, H0);
	}
	
	af[0][0] = delta[0][0] = -(ci[1] - a_h[1]) + I * (ci[2] - a_h[2]);
	af[1][1] = delta[1][1] = (ci[1] - a_h[1]) + I * (ci[2] - a_h[2]);
	af[0][1] = delta[0][1] = (ci[0] - a_h[0]) + (ci[3] - a_h[3]);
	af[1][0] = delta[1][0] = -(ci[0] - a_h[0]) + (ci[3] - a_h[3]);

	cmatrixmult(F0, af, af);
	cmatrixsum(+1, UNITY, af, af);
	c2x2matrixinv(af, af);
	cmatrixmult(af, H0, af);
	cmatrixmult(delta, af, af);
	cmatrixmult(G0, af, af);

	cf[0] = a_h[0] + 0.5 * (af[0][1] - af[1][0]);
	cf[3] = a_h[3] + 0.5 * (af[0][1] + af[1][0]);
	cf[1] = a_h[1] + 0.5 * (-af[0][0] + af[1][1]);
	cf[2] = a_h[2] - 0.5 * I * (af[0][0] + af[1][1]);

	for (j = 0; j <= 3; j++) {
		if (isnan(cabs(cf[j])) != 0.0) {
			printf("NAN in a-Riccati! Exiting.\n");
			exit(1);
		}
	}
}


/*
 * step_integratorm_a(double, double complex, double complex*, double complex*,
 *                  double complex*, double complex*)
 * 
 * This is step integrator for b-Riccati amplitude from point 1 with 
 * (Deltai, ci-known) to point 2 with (Deltaf, cf-unknown) separated by ds.
 * 
 * This integrator is valid for NON-UNITARY phases with 
 * 		\vq = \vDelta\times\vDelta^* != 0  (No scalar Delta[0] component)
 * and scalar self energy that is included in en
 * 		en=I*t(m+0.5)-\nu  or en=energy+I*0-\nu 
 * Still no magnetic field or vector self-energy, but maybe this is the way 
 * to implement them LATER.
 * 
 * This integrator is slower than the above  Non-Unitary integrator by
 * 1.5-2 times.
 * 
 * Arguments:
 *      double ds               ->
 *      double complex en       ->
 *      double complex *Deltai  ->
 *      double complex *Deltaf  ->
 *      double complex *ci      ->
 *      double complex *cf      ->
 * 
 * Returns:
 *      void
 */

void step_integratorm_b(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf)
{
	int j;
	double ampl, phi, Delta2;
	double complex E, q[4], q2, alpha, Norm; 
	double complex g_h[4], f_h[4], ff_h[4], aa_h[4], f2, ff2, f_ff;
	double complex Delta[4], m[4], Delta_aa, denom, lambda;
	double complex ee, eexp, eexp2, f[4]; 
	matrix2x2 S, S_inv, SS, SS_inv;
	matrix2x2 LAMBDA, UNITY, G0, H0, F0;
	matrix2x2 DELTA, delta, aaf;

	UNITY[0][0] = UNITY[1][1] = 1.0; 
	UNITY[0][1] = UNITY[1][0] = 0.0;
	ds = -ds;

	Delta2 = 0.0;
	if (ds > MAXds) {
		for (j = 1; j <= 3; j++) {
			Delta[j] = Deltaf[j];
			Delta2 += pow(cabs(Delta[j]), 2);
		}
	} else {
		for (j = 1; j <= 3; j++) {
			/* Delta[j]=0.5*(Deltai[j]+Deltaf[j]); */
			Delta[j] = Deltaf[j];
			Delta2 += pow(cabs(Delta[j]), 2);
		}
	}

	E = Delta2 - en * en;

    /*
     * q[1] = Delta[2]*conj(Delta[3]) - conj(Delta[2])*Delta[3];
     * q[2] = -Delta[1]*conj(Delta[3]) + conj(Delta[1])*Delta[3];
     * q[3] = Delta[1]*conj(Delta[2]) - conj(Delta[1])*Delta[2];
     */

	q[1] = -2 * I * (creal(Delta[2]) * cimag(Delta[3]) -
        cimag(Delta[2]) * creal(Delta[3]));
	q[2] = -2 * I * (-creal(Delta[1]) * cimag(Delta[3]) +
        cimag(Delta[1]) * creal(Delta[3]));
	q[3] = -2 * I * (creal(Delta[1]) * cimag(Delta[2]) -
        cimag(Delta[1]) * creal(Delta[2]));

	q2 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	alpha = 1.0 / (E + csqrt(E * E + q2));
	Norm = -1.0 / csqrt(2.0 * alpha) / csqrt(E * E + q2); 

	/* Bulk (homogeneous) propagator */
	g_h[0] = Norm * en;
	g_h[1] = -I * Norm * alpha * en * q[1];
	g_h[2] = -I * Norm * alpha * en * q[2];
	g_h[3] = -I * Norm * alpha * en * q[3];

	f_h[0] = ff_h[0] = 0.0;
	f_h[1] = -Norm * (Delta[1] - alpha * (Delta[2] * q[3] - Delta[3] * q[2]));
	f_h[2] = -Norm * (Delta[2] - alpha * (-Delta[1] * q[3] + Delta[3] * q[1]));
	f_h[3] = -Norm * (Delta[3] - alpha * (Delta[1] * q[2] - Delta[2] * q[1]));

	ff_h[1] = -Norm * (conj(Delta[1]) + alpha *
        (conj(Delta[2]) * q[3] - conj(Delta[3]) * q[2]));
	ff_h[2] = -Norm * (conj(Delta[2]) + alpha *
        (-conj(Delta[1]) * q[3] + conj(Delta[3]) * q[1]));
	ff_h[3] = -Norm * (conj(Delta[3]) + alpha *
        (conj(Delta[1]) * q[2] - conj(Delta[2]) * q[1]));

	f2 = f_h[1] * f_h[1] + f_h[2] * f_h[2] + f_h[3] * f_h[3];
	ff2 = ff_h[1] * ff_h[1] + ff_h[2] * ff_h[2] + ff_h[3] * ff_h[3];
	f_ff = f_h[1] * ff_h[1] + f_h[2] * ff_h[2] + f_h[3] * ff_h[3];
	denom = (I - g_h[0]) * (I - g_h[0]) +
        (f2 * ff2 - f_ff * f_ff) / 4 / g_h[0] / g_h[0];
	
	/* Bulk (homogeneous) Riccati a-amplitude */
	aa_h[0] = 0.0;
	aa_h[1] = 1.0 / denom * ((I - g_h[0]) * ff_h[1] +
        1.0 / 2 / g_h[0] * (ff_h[1] * f_ff - f_h[1] * ff2));
	aa_h[2] = 1.0 / denom * ((I - g_h[0]) * ff_h[2] +
        1.0 / 2 / g_h[0] * (ff_h[2] * f_ff - f_h[2] * ff2));
	aa_h[3] = 1.0 / denom * ((I - g_h[0]) * ff_h[3] +
        1.0 / 2 / g_h[0] * (ff_h[3] * f_ff - f_h[3] * ff2));

	if (ds > MAXds) {
		for (j = 0; j <= 3; j++)
            cf[j] = aa_h[j];

		return;
	}

	Delta_aa = aa_h[1] * (Delta[1]) + aa_h[2] * (Delta[2]) + aa_h[3] * (Delta[3]);
	ee = -en + Delta_aa;
	eexp = cexp(I * ee * ds);
	eexp2 = eexp * eexp;

	/* Diagonalization and so on... */
	m[1] = aa_h[2] * (Delta[3]) - aa_h[3] * (Delta[2]);
	m[2] = aa_h[3] * (Delta[1]) - aa_h[1] * (Delta[3]);
	m[3] = aa_h[1] * (Delta[2]) - aa_h[2] * (Delta[1]);

	if (cabs(m[3]) != 0.0)
		lambda = I * m[3] * csqrt(1.0 + cpow(m[1] / m[3], 2) +
            cpow(m[2] / m[3], 2));
	else
		lambda = I * csqrt(m[1] * m[1] + m[2] * m[2] + m[3] * m[3]);

	if (cabs(lambda) == 0.0) {
		G0[0][0] = eexp * (1.0 + m[3] * ds);
		G0[1][1] = eexp * (1.0 - m[3] * ds);
		G0[0][1] = eexp * (m[1] + I * m[2]) * ds;
		G0[1][0] = eexp * (m[1] - I * m[2]) * ds;

		H0[0][0] = eexp * (1.0 + m[3] * ds);
		H0[1][1] = eexp * (1.0 - m[3] * ds);
		H0[0][1] = eexp * (m[1] - I * m[2]) * ds;
		H0[1][0] = eexp * (m[1] + I * m[2]) * ds;

		f[1] = (eexp2 - 1.0) / 2 / ee *
			((Delta[1]) + ((Delta[2]) * m[3] - (Delta[3]) * m[2]) / ee) -
			I * eexp2 / ee * ds * ((Delta[2]) * m[3] - (Delta[3]) * m[2]);
		f[2] = (eexp2 - 1.0) / 2 / ee *
			((Delta[2]) + ((Delta[3]) * m[1] - (Delta[1]) * m[3]) / ee) -
			I * eexp2 / ee * ds * ((Delta[3]) * m[1] - (Delta[1]) * m[3]);
		f[3] = (eexp2 - 1.0) / 2 / ee *
			((Delta[3]) + ((Delta[1]) * m[2] - (Delta[2]) * m[1]) / ee) -
			I * eexp2 / ee * ds * ((Delta[1]) * m[2] - (Delta[2]) * m[1]);

		F0[0][0] = -f[1] + I * f[2];
		F0[1][1] = f[1] + I * f[2];
		F0[0][1] = F0[1][0] = f[3];
	} else {
		LAMBDA[0][0] = lambda;
		LAMBDA[1][1] = -lambda;
		LAMBDA[0][1] = LAMBDA[1][0] = 0.0;

		ampl = sqrt(0.5 * cabs(1 + I * m[3] / lambda));
		phi = -0.5 * carg(2.0 / (1 + I * m[3] / lambda));

		S[0][0] = ampl * cexp(I * phi);
		S[1][1] = ampl * cexp(-I * phi);
		S[0][1] = -S[1][1] * (m[1] + I * m[2]) / (m[3] - I * lambda);
		S[1][0] = S[0][0] * (m[1] - I * m[2]) / (m[3] - I * lambda);

		S_inv[0][0] = S[1][1];
		S_inv[1][1] = S[0][0];
		S_inv[0][1] = -S[0][1];
		S_inv[1][0] = -S[1][0];
		
		SS[0][0] = ampl * cexp(I * phi);
		SS[1][1] = ampl * cexp(-I * phi);
		SS[0][1] = -SS[1][1] * (m[1] - I * m[2]) / (m[3] - I * lambda);
		SS[1][0] = SS[0][0] * (m[1] + I * m[2]) / (m[3] - I * lambda);

		SS_inv[0][0] = SS[1][1];
		SS_inv[1][1] = SS[0][0];
		SS_inv[0][1] = -SS[0][1];
		SS_inv[1][0] = -SS[1][0];

		DELTA[0][0] = -(Delta[1]) + I * (Delta[2]);
		DELTA[1][1] = (Delta[1]) + I * (Delta[2]);
		DELTA[0][1] = (Delta[3]);
		DELTA[1][0] = (Delta[3]);

		cmatrixmult(DELTA, S, DELTA);
		cmatrixmult(SS_inv, DELTA, DELTA);

		F0[0][0] = (eexp2 * cexp(-I * (LAMBDA[0][0] + LAMBDA[0][0]) * ds) - 1) /
			(2 * ee - LAMBDA[0][0] - LAMBDA[0][0]) * DELTA[0][0];
		F0[1][1] = (eexp2 * cexp(-I * (LAMBDA[1][1] + LAMBDA[1][1]) * ds) - 1) /
			(2 * ee - LAMBDA[1][1] - LAMBDA[1][1]) * DELTA[1][1];
		F0[0][1] = (eexp2 * cexp(-I * (LAMBDA[0][0] + LAMBDA[1][1]) * ds) - 1) /
			(2 * ee - LAMBDA[0][0] - LAMBDA[1][1]) * DELTA[0][1];
		F0[1][0] = (eexp2 * cexp(-I * (LAMBDA[1][1] + LAMBDA[0][0]) * ds) - 1) /
			(2 * ee - LAMBDA[1][1] - LAMBDA[0][0]) * DELTA[1][0];

		cmatrixmult(F0, S_inv, F0);
		cmatrixmult(SS, F0, F0);
		/* found F0, now find G0, H0 */

		G0[0][0] = H0[0][0] = eexp * cexp(-I * LAMBDA[0][0] * ds);
		G0[1][1] = H0[1][1] = eexp * cexp(-I * LAMBDA[1][1] * ds);
		G0[0][1] = H0[0][1] = 0.0;
		G0[1][0] = H0[1][0] = 0.0;
		cmatrixmult(S, G0, G0);
		cmatrixmult(G0, S_inv, G0);
		cmatrixmult(SS, H0, H0);
		cmatrixmult(H0, SS_inv, H0);
	}
	
	aaf[0][0] = delta[0][0] = (ci[1] - aa_h[1]) + I * (ci[2] - aa_h[2]);
	aaf[1][1] = delta[1][1] = -(ci[1] - aa_h[1]) + I * (ci[2] - aa_h[2]);
	aaf[0][1] = delta[0][1] = (ci[0] - aa_h[0]) - (ci[3] - aa_h[3]);
	aaf[1][0] = delta[1][0] = -(ci[0] - aa_h[0]) - (ci[3] - aa_h[3]);

	cmatrixmult(F0, aaf, aaf);
	cmatrixsum(+1, UNITY, aaf, aaf);
	c2x2matrixinv(aaf, aaf);
	cmatrixmult(aaf, H0, aaf);
	cmatrixmult(delta, aaf, aaf);
	cmatrixmult(G0, aaf, aaf);

	cf[0] = aa_h[0] + 0.5 * (aaf[0][1] - aaf[1][0]);
	cf[3] = aa_h[3] - 0.5 * (aaf[0][1] + aaf[1][0]);
	cf[1] = aa_h[1] + 0.5 * (aaf[0][0] - aaf[1][1]);
	cf[2] = aa_h[2] - 0.5 * I*(aaf[0][0] + aaf[1][1]);

	for (j = 0; j <= 3; j++) {
		if (isnan(cabs(cf[j])) != 0.0) {
			printf("NAN in b-Riccati! Exiting.\n");
			exit(1);
		}
	}
}


/*
 * cmatrixsum(int, matrix2x2, matrix2x2, matrix2x2)
 * 
 * (?)
 * 
 * Arguments:
 *      int sgn         ->
 *      matrix2x2 A     ->
 *      matrix2x2 B     ->
 *      matrix2x2 C     ->
 * 
 * Returns:
 *      void
 */

void cmatrixsum(int sgn, matrix2x2 A, matrix2x2 B, matrix2x2 C)
{
	int m, n;

	for (m = 0; m < 2; m++)
		for (n = 0; n < 2; n++)
			C[m][n] = A[m][n] + sgn * B[m][n];
}


/*
 * cmatrixmult(matrix2x2, matrix2x2, matrix2x2)
 * 
 * (?)
 * 
 * Arguments:
 *      matrix2x2 A     ->
 *      matrix2x2 B     ->
 *      matrix2x2 C     ->
 * 
 * Returns:
 *      void
 */
	
void cmatrixmult(matrix2x2 A, matrix2x2 B, matrix2x2 C)
{
	int m, n, k;
	matrix2x2 a, b;

	for (m = 0; m < 2; m++) {
		for (n = 0; n < 2; n++) {
			a[m][n] = A[m][n];
			b[m][n] = B[m][n];
		}
    }
	for (m = 0; m < 2; m++) {
		for (n = 0; n < 2; n++) {
			C[m][n] = 0.0;
			for (k = 0; k < 2; k++) {
				C[m][n] += a[m][k] * b[k][n];
            }
		}
    }
}


/*
 * c2x2matrixinv(matrix2x2, matrix2x2)
 * 
 * 
 * Arguments:
 *      matrix2x2 C         ->
 *      matrix2x2 C_inv     ->
 * 
 * Returns:
 *      void
 */

void c2x2matrixinv(matrix2x2 C, matrix2x2 C_inv)
{
	double complex DetC, tmp;

	DetC = C[0][0] * C[1][1] - C[0][1] * C[1][0];

	if (cabs(DetC) < 1.0e-10) {
		printf("Cannot invert matrix, |DetC|=%e\n", cabs(DetC));
		exit(1);
	}

	tmp = C[0][0];                  /* in case **C_inv = **C */
	C_inv[0][0] = C[1][1] / DetC;
	C_inv[1][1] = tmp / DetC;
	C_inv[0][1] = -C[0][1] / DetC;
	C_inv[1][0] = -C[1][0] / DetC;
}
