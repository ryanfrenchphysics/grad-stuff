#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ricparam.h"

#define TINYsr      1.0e-20

/* TODO: Add info on args / function information */

/*
 * Riccati matrices are defined through components as
 * ra = (ra[0] + ra[1,2,3]*sigma_1,2,3) * (i sigma_2)
 * rb = (i sigma_2) * (rb[0] + rb[1,2,3]*sigma_1,2,3)
 * 
 * Green's functions are normalized to G^2 = -1
 * 
 *     [ g    f ]
 * G = [        ]
 *     [ ff  gg ]
 * 
 * g  = g[0] + g[1, 2, 3] * sigma_1,2,3
 * gg = gg[0] + gg[1, 2, 3] * sigma^tr_1,2,3
 * f = (f[0] + f[1, 2, 3] * sigma_1,2,3) * (i sigma_2)
 * ff = (i sigma_2) * (ff[0] + ff[1, 2, 3] * sigma_1,2,3)
 * 
 * and in terms of Riccati amplitudes
 * g = -I / (1 - ra * rb) * (1 + ra * rb)
 * gg = I / (1 - rb * ra) * (1 + rb * ra)
 * f = -I / (1 - ra * rb) * 2 * ra
 * ff = I / (1 - rb * ra) * 2 * rb
 */


/*
 * void ab_bulk(double complex*, double complex*,
 *              double complex, double complex*)
 * 
 * Unitary Triplet
 * Im(en) should be > 0
 * 
 * Args:
 *      double complex *ra      ->
 *      double complex *rb      ->
 *      double complex en       ->
 *      double complex *Delta   ->
 * 
 * Returns:
 *      void
 */

void ab_bulk(double complex *ra, double complex *rb, 
		double complex en, double complex *Delta)
{
	int j;
	double complex denom;

	/*
     * for (j = 1; j <= 3; j++)
     *      printf("Delta[%d] = %f + I(%f)\n", j, creal(Delta[j]), cimag(Delta[j]));
     */
	/* fflush(stdout); */

	denom = pow(cabs(Delta[1]), 2) + pow(cabs(Delta[2]), 2) +
            pow(cabs(Delta[3]), 2);
	denom = en + I * csqrt(denom - en * en);

	ra[0] = rb[0] = 0.0;
	for (j = 1; j <= 3; j++) {
		ra[j] = -Delta[j] / denom;
		rb[j] = conj(Delta[j]) / denom;
	}
}

/* 
 * ab_bulk_non_unitary(double complex*, double complex*,
 *                      double complex, double complex*)
 * 
 * non_Unitary uniform:  no (i\sigma_2) factors already
 * a = -1 / [en + I * sqrt(Delta * _Delta_ - en2)] * Delta 
 * b =  _Delta_ * 1 / [en + I * sqrt(Delta * _Delta_ - en2)]
 * 
 * Delta can have singlet [0] and triplet [1,2,3] components
 * !! - No Magnetic field
 * 
 * Im(en) should be > 0
 * 
 * Args:
 *      double complex *ra      ->
 *      double complex *rb      ->
 *      double complex en       ->
 *      double complex *Delta   ->
 * 
 * Returns:
 *      void
 */
void ab_bulk_non_unitary(double complex *ra, double complex *rb, 
		double complex en, double complex *Delta)
{
	int j;
	double Delta2, S[4];
	double complex DeltaC[4], q[4], sqrtq[4], den[4]; 

	Delta2 = 0.0;
	for (j = 0; j <= 3; j++) {
		Delta2 += pow(cabs(Delta[j]), 2);
		DeltaC[j] = conj(Delta[j]); 
	}
	
	/* Spin of the Cooper pairs */
	S[1] = 2 * creal(Delta[0] * DeltaC[1]) - 2 * cimag(Delta[2] * DeltaC[3]);
	S[2] = 2 * creal(Delta[0] * DeltaC[2]) - 2 * cimag(Delta[3] * DeltaC[1]);
	S[3] = 2 * creal(Delta[0] * DeltaC[3]) - 2 * cimag(Delta[1] * DeltaC[2]);

    /* q = gap2 - en2 + \vsigma * \vS */
	q[0] = Delta2 - en * en;  
	q[1] = S[1]; 
	q[2] = S[2];
	q[3] = S[3]; 

    /* den = en + I * sqrt(q) */
	ric_sqrt(q, sqrtq);
	den[0] = en + I * sqrtq[0];
	den[1] = I * sqrtq[1];
	den[2] = I * sqrtq[2];
	den[3] = I * sqrtq[3];

    /* q = 1 / den */
	ric_invert(den, q);	
	ric_multiply(q, Delta, ra);

	for (j = 0; j <= 3; j++)
        ra[j] *= -1;    /* ra = -q * Delta */

    /* rb = _Delta_ * q */
	ric_multiply(DeltaC, q, rb); 
}


/*
 * a_from_g(double complex*, double complex*, double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *ra      ->
 *      double complex *g       ->
 *      double complex *f       ->
 *      double complex *ff      ->
 * 
 * Returns:
 *      void
 */

void a_from_g(double complex *ra, double complex *g, 
		double complex *f, double complex *ff)
{
	double complex gf, gxf[4], g2, denom;

	gf = g[1] * f[1] + g[2] * f[2] + g[3] * f[3];
	g2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
	gxf[1] = g[2] * f[3] - g[3] * f[2];
	gxf[2] = -g[1] * f[3] + g[3] * f[1];
	gxf[3] = g[1] * f[2] - g[2] * f[1];

	denom = (I - g[0]) * (I - g[0]) - g2;

	ra[0] = -1.0 / denom * (f[0] * (I - g[0]) + gf);
	ra[1] = -1.0 / denom * (f[1] * (I - g[0]) + f[0] * g[1] + I * gxf[1]);
	ra[2] = -1.0 / denom * (f[2] * (I - g[0]) + f[0] * g[2] + I * gxf[2]);
	ra[3] = -1.0 / denom * (f[3] * (I - g[0]) + f[0] * g[3] + I * gxf[3]);
}


/*
 * b_from_g(double complex*, double complex*, double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *rb      ->
 *      double complex *gb      ->
 *      double complex *f       ->
 *      double complex *ff      ->
 * 
 * Returns:
 *      void
 */

void b_from_g(double complex *rb, double complex *gg, 
		double complex *f, double complex *ff)
{
	double complex ggff, ggxff[4], gg2, denom;

	ggff = gg[1] * ff[1] + gg[2] * ff[2] + gg[3] * ff[3];
	gg2 = gg[1] * gg[1] + gg[2] * gg[2] + gg[3] * gg[3];
	ggxff[1] = gg[2] * ff[3] - gg[3] * ff[2];
	ggxff[2] = -gg[1] * ff[3] + gg[3] * ff[1];
	ggxff[3] = gg[1] * ff[2] - gg[2] * ff[1];

	denom = cpow(I + gg[0], 2) - gg2;
    /* denom = (I+gg[0])*(I+gg[0]) - gg2; */

	rb[0] = 1.0 / denom * (ff[0] * (I + gg[0]) + ggff);
	rb[1] = 1.0 / denom * (ff[1] * (I + gg[0]) + ff[0] * gg[1] + I * ggxff[1]);
	rb[2] = 1.0 / denom * (ff[2] * (I + gg[0]) + ff[0] * gg[2] + I * ggxff[2]);
	rb[3] = 1.0 / denom * (ff[3] * (I + gg[0]) + ff[0] * gg[3] + I * ggxff[3]);
}


/*
 * g_from_ric(double complex*, double complex*, double complex*)
 * 
 * g = -I / (1 - ra * rb) * (1 + ra * rb)
 * 
 * Args:
 *      double complex *g       ->
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      void
 */

void g_from_ric(double complex *g, double complex *ra, double complex *rb)
{
	double complex ab, axb[4], axb2, denom;

	ab = ra[0] * rb[0];
	ab += ra[1] * rb[1];
	ab += ra[2] * rb[2];
	ab += ra[3] * rb[3];

	axb[1] = ra[0] * rb[1] + ra[1] * rb[0] + I * (ra[2] * rb[3] - ra[3] * rb[2]);
	axb[2] = ra[0] * rb[2] + ra[2] * rb[0] + I * (ra[3] * rb[1] - ra[1] * rb[3]);
	axb[3] = ra[0] * rb[3] + ra[3] * rb[0] + I * (ra[1] * rb[2] - ra[2] * rb[1]);

	axb2 = axb[1] * axb[1];
	axb2 += axb[2] * axb[2];
	axb2 += axb[3] * axb[3];

	denom = (1 + ab) * (1 + ab) - axb2;

	g[0] = -I / denom * (1 - ab * ab + axb2);
	g[1] = 2.0 * I / denom * axb[1];
	g[2] = 2.0 * I / denom * axb[2];
	g[3] = 2.0 * I / denom * axb[3];
}


/*
 * gg_from_ric(double complex*, double complex*, double complex*)
 * 
 * gg = I / (1 - rb * ra) * (1 + rb * ra)
 * 
 * Args:
 *      double complex *gg      ->
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      void
 */

void gg_from_ric(double complex *gg, double complex *ra, double complex *rb)
{
	double complex ab, axb[4], axb2, denom;

	ab = ra[0] * rb[0];
	ab += ra[1] * rb[1];
	ab += ra[2] * rb[2];
	ab += ra[3] * rb[3];

	axb[1] = ra[0] * rb[1] + ra[1] * rb[0] - I * (ra[2] * rb[3] - ra[3] * rb[2]);
	axb[2] = ra[0] * rb[2] + ra[2] * rb[0] - I * (ra[3] * rb[1] - ra[1] * rb[3]);
	axb[3] = ra[0] * rb[3] + ra[3] * rb[0] - I * (ra[1] * rb[2] - ra[2] * rb[1]);

	axb2 = axb[1] * axb[1];
	axb2 += axb[2] * axb[2];
	axb2 += axb[3] * axb[3];

	denom = (1 + ab) * (1 + ab) - axb2;

	gg[0] = I / denom * (1 - ab * ab + axb2);
	gg[1] = 2.0 * I / denom * axb[1];
	gg[2] = 2.0 * I / denom * axb[2];
	gg[3] = 2.0 * I / denom * axb[3];
}


/*
 * f_from_ric(double complex*, double complex*,
 *              double complex*, double complex*)
 * 
 * f = -I / (1 - ra * rb) * 2 * ra
 * ff = I / (1 - rb * ra) * 2 * rb
 * 
 * Args:
 *      double complex *f       ->
 *      double complex *ff      ->
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      void
 */

void f_from_ric(double complex *f, double complex *ff,
		double complex *ra, double complex *rb)
{
	double complex ab, av2, bv2, axb[4], axb2, AXB[4], AXB2, denom, DENOM;

	ab = ra[0] * rb[0];
	ab += ra[1] * rb[1];
	ab += ra[2] * rb[2];
	ab += ra[3] * rb[3];

	av2 = ra[1] * ra[1] + ra[2] * ra[2] + ra[3] * ra[3];
	bv2 = rb[1] * rb[1] + rb[2] * rb[2] + rb[3] * rb[3];

	axb[1] = ra[0] * rb[1] + ra[1] * rb[0] + I * (ra[2] * rb[3] - ra[3] * rb[2]);
	axb[2] = ra[0] * rb[2] + ra[2] * rb[0] + I * (ra[3] * rb[1] - ra[1] * rb[3]);
	axb[3] = ra[0] * rb[3] + ra[3] * rb[0] + I * (ra[1] * rb[2] - ra[2] * rb[1]);

	axb2 = axb[1] * axb[1];
	axb2 += axb[2] * axb[2];
	axb2 += axb[3] * axb[3];

	denom = (1 + ab) * (1 + ab) - axb2;

	AXB[1] = ra[0] * rb[1] + ra[1] * rb[0] - I * (ra[2] * rb[3] - ra[3] * rb[2]);
	AXB[2] = ra[0] * rb[2] + ra[2] * rb[0] - I * (ra[3] * rb[1] - ra[1] * rb[3]);
	AXB[3] = ra[0] * rb[3] + ra[3] * rb[0] - I * (ra[1] * rb[2] - ra[2] * rb[1]);

	AXB2 = AXB[1] * AXB[1];
	AXB2 += AXB[2] * AXB[2];
	AXB2 += AXB[3] * AXB[3];

	DENOM = (1 + ab) * (1 + ab) - AXB2;

	f[0] = -2.0 * I / denom * (ra[0] + ra[0] * ra[0] * rb[0] - rb[0] * av2);
	f[1] = -2.0 * I / denom * (ra[1] - ra[0] * ra[0] * rb[1] + rb[1] * av2);
	f[2] = -2.0 * I / denom * (ra[2] - ra[0] * ra[0] * rb[2] + rb[2] * av2);
	f[3] = -2.0 * I / denom * (ra[3] - ra[0] * ra[0] * rb[3] + rb[3] * av2);
	
	ff[0] = 2.0 * I / DENOM * (rb[0] + rb[0] * rb[0] * ra[0] - ra[0] * bv2);
	ff[1] = 2.0 * I / DENOM * (rb[1] - rb[0] * rb[0] * ra[1] + ra[1] * bv2);
	ff[2] = 2.0 * I / DENOM * (rb[2] - rb[0] * rb[0] * ra[2] + ra[2] * bv2);
	ff[3] = 2.0 * I / DENOM * (rb[3] - rb[0] * rb[0] * ra[3] + ra[3] * bv2);
}


/*
 * Manipulations with matrices of the kind: x = x[0] + x[1,2,3] * sigma_1,2,3
 * vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 */


/*
 * ric_multiply(double complex*, double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *x           ->
 *      double complex *y           ->
 *      double complex *result      ->
 * 
 * Returns:
 *      void
 */

void ric_multiply(double complex *x, double complex *y, double complex *result)
{
	double complex tmp[4]; 

	tmp[0] = x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]; 
	tmp[1] = x[0] * y[1] + x[1] * y[0] + I * (x[2] * y[3] - x[3] * y[2]);
	tmp[2] = x[0] * y[2] + x[2] * y[0] + I * (x[3] * y[1] - x[1] * y[3]);
	tmp[3] = x[0] * y[3] + x[3] * y[0] + I * (x[1] * y[2] - x[2] * y[1]);

	result[0] = tmp[0]; 
	result[1] = tmp[1]; 
	result[2] = tmp[2]; 
	result[3] = tmp[3]; 
}


/*
 * ric_invert(double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *x           ->
 *      double complex *x_inv       ->
 * 
 * Returns:
 *      void
 */

void ric_invert(double complex *x, double complex *x_inv)
{
	double complex den;

	den = x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3]; 
	x_inv[0] =  x[0] / den; 
	x_inv[1] = -x[1] / den; 
	x_inv[2] = -x[2] / den; 
	x_inv[3] = -x[3] / den;
}

/*
 * ric_sqrt(double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *x           ->
 *      double complex *x_sqrt      ->
 * 
 * Returns:
 *      void
 */

void ric_sqrt(double complex *x, double complex *x_sqrt)
{
	double complex den;

	/* den = x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3]; */
	/* x_sqrt[0] = 0.5 * (x[0] + csqrt(den)); */

	den = 1 - (x[1] * x[1] + x[2] * x[2] + x[3] * x[3]) / (x[0] * x[0]); 
	x_sqrt[0] = 0.5 * (1.0 + csqrt(den));
	x_sqrt[0] = csqrt(x[0]) * csqrt(x_sqrt[0]);
	x_sqrt[1] = x[1] / 2 / x_sqrt[0]; 
	x_sqrt[2] = x[2] / 2 / x_sqrt[0]; 
	x_sqrt[3] = x[3] / 2 / x_sqrt[0];
}


/*
 * Manipulations with Greens Function Matrices
 * vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * 
 * For this purpose I define a new parametrization
 * of GFs to avoid sigma^tr in gg:
 * 
 * g = G
 * gg = (-i sigma_2) GG (i sigma_2)
 * f = F (i sigma_2)
 * ff = (i sigma_2) FF
 * 
 * All matrices G, GG, F, FF are of the form x[0] + x[1,2,3] * sigma_1,2,3
 * 
 * gg-GG connection: gg[0] = GG[0], gg[1,2,3] = -GG[1,2,3]
 * 
 * C = AB  in terms of this parametrization
 * CG = AG * BG - AF * BFF
 * CGG= AGG * BGG - AFF * BF
 * CF = AG * BF + AF * BGG
 * CFF= AFF * BG + AGG * BFF
 */


/*
 * G_multiply(double complex*, double complex*, double complex*,
 *              double complex*, double complex*, double complex*,
 *              double complex*, double complex*, double complex*,
 *              double complex*, double complex*, double complex*)
 * 
 * Description?
 * 
 * Args:
 *      double complex *Ag      ->
 *      double complex *Agg     ->
 *      double complex *Af      ->
 *      double complex *Aff     ->
 *		double complex *Bg      ->
 *      double complex *Bgg     ->
 *      double complex *Bf      ->
 *      double complex *Bff     ->
 *		double complex *Cg      ->
 *      double complex *Cgg     ->
 *      double complex *Cf      ->
 *      double complex *Cff     ->
 * 
 * Returns:
 *      void
 */

void G_multiply(double complex *Ag, double complex *Agg, double complex *Af,
            double complex *Aff, double complex *Bg, double complex *Bgg,
            double complex *Bf, double complex *Bff, double complex *Cg,
            double complex *Cgg, double complex *Cf, double complex *Cff)
{
	int i;
	double complex tmpG1[4], tmpG2[4], tmpF1[4], tmpF2[4];
    double complex tmpFF1[4], tmpFF2[4], tmpGG1[4], tmpGG2[4];

	/* Get new parametrization of gg-functions */
	for (i = 1; i <= 3; i++) {
        Agg[i] *= -1;
        Bgg[i] *= -1;
    }

	/* Now use matrix multiplication and inversion to get the results */
	ric_multiply(Ag, Bg, tmpG1);
    ric_multiply(Af, Bff, tmpG2); 
	ric_multiply(Ag, Bf, tmpF1);
    ric_multiply(Af, Bgg, tmpF2); 
	ric_multiply(Aff, Bg, tmpFF1);
    ric_multiply(Agg, Bff, tmpFF2); 
	ric_multiply(Agg, Bgg, tmpGG1);
    ric_multiply(Aff, Bf, tmpGG2); 

	/* Return back to traditional parametrization of gg-functions */
	for (i = 1; i <= 3; i++) {
        Agg[i] *= -1;
        Bgg[i] *= -1;
    }

	/* 
     * Get the C-matrix. Should be done here in case *C
     * points to either A or B again
     */
	for (i = 0; i <= 3; i++)
        Cg[i] = tmpG1[i] - tmpG2[i]; 
	for (i = 0; i <= 3; i++)
        Cf[i] = tmpF1[i] + tmpF2[i]; 
	for (i = 0; i <= 3; i++)
        Cff[i] = tmpFF1[i] + tmpFF2[i]; 
	for (i = 0; i <= 3; i++)
        Cgg[i] = tmpGG1[i] - tmpGG2[i]; 
	for (i = 1; i <= 3; i++)
        Cgg[i] *= -1;     /* Return to traditional parametrization */


	/*Cg[0]=Ag[0]*Bg[0]+Ag[1]*Bg[1]+Ag[2]*Bg[2]+Ag[3]*Bg[3] -Af[0]*Bff[0]-Af[1]*Bff[1]-Af[2]*Bff[2]-Af[3]*Bff[3];
	Cg[1]=Ag[0]*Bg[1]+Ag[1]*Bg[0]+I*(Ag[2]*Bg[3]-Ag[3]*Bg[2]) -Af[0]*Bff[1]-Af[1]*Bff[0]-I*(Af[2]*Bff[3]-Af[3]*Bff[2]); 
	Cg[2]=Ag[0]*Bg[2]+Ag[2]*Bg[0]+I*(Ag[3]*Bg[1]-Ag[1]*Bg[3]) -Af[0]*Bff[2]-Af[2]*Bff[0]-I*(Af[3]*Bff[1]-Af[1]*Bff[3]); 
	Cg[3]=Ag[0]*Bg[3]+Ag[3]*Bg[0]+I*(Ag[1]*Bg[2]-Ag[2]*Bg[1]) -Af[0]*Bff[3]-Af[3]*Bff[0]-I*(Af[1]*Bff[2]-Af[2]*Bff[1]); 

	Cgg[0]=Agg[0]*Bgg[0]+Agg[1]*Bgg[1]+Agg[2]*Bgg[2]+Agg[3]*Bgg[3] -Af[0]*Bff[0]-Af[1]*Bff[1]-Af[2]*Bff[2]-Af[3]*Bff[3];
	Cgg[1]=Agg[0]*Bgg[1]+Agg[1]*Bgg[0]-I*(Agg[2]*Bgg[3]-Agg[3]*Bgg[2]) +Aff[0]*Bf[1]+Aff[1]*Bf[0]+I*(Aff[2]*Bf[3]-Aff[3]*Bf[2]); 
	Cgg[2]=Agg[0]*Bgg[2]+Agg[2]*Bgg[0]-I*(Agg[3]*Bgg[1]-Agg[1]*Bgg[3]) +Aff[0]*Bf[2]+Aff[2]*Bf[0]+I*(Aff[3]*Bf[1]-Aff[1]*Bf[3]); 
	Cgg[3]=Agg[0]*Bgg[3]+Agg[3]*Bgg[0]-I*(Agg[1]*Bgg[2]-Agg[2]*Bgg[1]) +Aff[0]*Bf[3]+Aff[3]*Bf[0]+I*(Aff[1]*Bf[2]-Aff[2]*Bf[1]); 

	Cf[0]= Ag[0]*Bf[0]+Ag[1]*Bf[1]+Ag[2]*Bf[2]+Ag[3]*Bf[3] +Af[0]*Bgg[0]-Af[1]*Bgg[1]-Af[2]*Bgg[2]-Af[3]*Bgg[3];
	Cf[1]= Ag[0]*Bf[1]+Ag[1]*Bf[0]+I*(Ag[2]*Bf[3]-Ag[3]*Bf[2]) +Af[1]*Bgg[0]-Af[0]*Bgg[1]-I*(Af[2]*Bgg[3]-Af[3]*Bgg[2]);
	Cf[2]= Ag[0]*Bf[2]+Ag[2]*Bf[0]+I*(Ag[3]*Bf[1]-Ag[1]*Bf[3]) +Af[2]*Bgg[0]-Af[0]*Bgg[2]-I*(Af[3]*Bgg[1]-Af[1]*Bgg[3]);
	Cf[3]= Ag[0]*Bf[3]+Ag[3]*Bf[0]+I*(Ag[1]*Bf[2]-Ag[2]*Bf[1]) +Af[3]*Bgg[0]-Af[0]*Bgg[3]-I*(Af[1]*Bgg[2]-Af[2]*Bgg[1]);

	Cf[0]= Agg[0]*Bff[0]-Agg[1]*Bff[1]-Agg[2]*Bff[2]-Agg[3]*Bff[3] +Aff[0]*Bg[0]+Aff[1]*Bg[1]+Aff[2]*Bg[2]+Aff[3]*Bg[3];
	Cf[1]= Agg[0]*Bff[1]-Agg[1]*Bff[0]-I*(Agg[2]*Bff[3]-Agg[3]*Bff[2]) +Aff[1]*Bg[0]+Aff[0]*Bg[1]+I*(Aff[2]*Bg[3]-Aff[3]*Bg[2]);
	Cf[2]= Agg[0]*Bff[2]-Agg[2]*Bff[0]-I*(Agg[3]*Bff[1]-Agg[1]*Bff[3]) +Aff[2]*Bg[0]+Aff[0]*Bg[2]+I*(Aff[3]*Bg[1]-Aff[1]*Bg[3]);
	Cf[3]= Agg[0]*Bff[3]-Agg[3]*Bff[0]-I*(Agg[1]*Bff[2]-Agg[2]*Bff[1]) +Aff[3]*Bg[0]+Aff[0]*Bg[3]+I*(Aff[1]*Bg[2]-Aff[2]*Bg[1]);*/
}


/*
 * G_multiply(double complex*, double complex*, double complex*,
 *              double complex*, double complex*, double complex*,
 *              double complex*, double complex*)
 * 
 * B = inverse(A)
 * 
 * Args:
 *      double complex *Ag      ->
 *      double complex *Agg     ->
 *      double complex *Af      ->
 *      double complex *Aff     ->
 *		double complex *Bg      ->
 *      double complex *Bgg     ->
 *      double complex *Bf      ->
 *      double complex *Bff     ->
 * 
 * Returns:
 *      void
 */

void G_invert(double complex *Ag, double complex *Agg, double complex *Af,
            double complex *Aff, double complex *Bg, double complex *Bgg,
            double complex *Bf, double complex *Bff)
{
	int i;
	double complex tmp[4], g[4], gg[4], f[4], ff[4];

	/* Get new parametrization of gg-functions */
	for (i = 1; i <= 3; i++)
        Agg[i] *= -1;
	
	ric_invert(Agg, g); 
	ric_multiply(g, Aff, tmp);
	ric_multiply(Af, tmp, g);

	for (i = 0; i <= 3; i++)
        g[i] += Ag[i];

	ric_invert(g, g);
	ric_multiply(tmp, g, ff);
	ric_invert(Ag, gg); 
	ric_multiply(gg, Af, tmp);
	ric_multiply(Aff, tmp, gg);

	for (i = 0; i <= 3; i++)
        gg[i] += Agg[i];

	ric_invert(gg,gg);
	ric_multiply(tmp,gg,f);

	/* Get traditional parametrization of gg-functions */
	for (i = 1; i <= 3; i++)
        Agg[i] *= -1;

	for (i = 0; i <= 3; i++) {
        Bg[i] = g[i];
        Bf[i] = -f[i];
        Bff[i] = -ff[i];
        Bgg[i] = gg[i];
    }

	for (i = 1; i <= 3; i++)
        Bgg[i] *= -1;
}





/*
 * Some Greens Functions
 * vvvvvvvvvvvvvvvvvvvvv
 */


/*
 * g_bulk(double complex*, double complex*, double complex*, double complex*,
 *          double complex, double complex*)
 * 
 * Non-unitary phases, NO magnetic field.
 * Delta can have singlet [0] and triplet [1,2,3] components
 * 
 * Args:
 *      double complex *g       ->
 *      double complex *gg      ->
 *      double complex *f       ->
 *      double complex *ff      ->
 *      double complex en       ->
 *      double complex *Delta   ->
 * 
 * Returns:
 *      void
 */

void g_bulk(double complex *g, double complex *gg, double complex *f,
            double complex *ff, double complex en, double complex *Delta)
{
	int j;
	double Delta2, S, Svec[4];
	double complex Epls, Emns, Ep, Em;
	
	Delta2 = 0.0; /* gap */

	for (j = 0; j <= 3; j++) {
		Delta2 += pow(cabs(Delta[j]), 2);
	}
	
	/* Spin of the Cooper pairs */
	Svec[1] = 2 * creal(Delta[0] * conj(Delta[1])) -
        2 * cimag(Delta[2] * conj(Delta[3]));
	Svec[2] = 2 * creal(Delta[0] * conj(Delta[2])) -
        2 * cimag(Delta[3] * conj(Delta[1]));
	Svec[3] = 2 * creal(Delta[0] * conj(Delta[3])) -
        2 * cimag(Delta[1] * conj(Delta[2]));
	
	S = sqrt(pow(Svec[1], 2) + pow(Svec[2], 2) + pow(Svec[3], 2));

	/* Unit spin vector */
	Svec[1] /= (S + 1e-10);
	Svec[2] /= (S + 1e-10);
	Svec[3] /= (S + 1e-10);

	Epls = csqrt(Delta2 + S - en * en); 
	Emns = csqrt(Delta2 - S - en * en); 

	Ep = 0.5 * (1 / Epls + 1 / Emns); 
	Em = 0.5 * (1 / Epls - 1 / Emns); 
	
	/* Uniform propagator */
	g[0] = -en * Ep; 
	g[1] = -en * Em * Svec[1]; 
	g[2] = -en * Em * Svec[2]; 
	g[3] = -en * Em * Svec[3];

	gg[0] =  en * Ep;
	gg[1] =  en * Em * Svec[1]; 
	gg[2] =  en * Em * Svec[2]; 
	gg[3] =  en * Em * Svec[3];

	f[0] = Delta[0] * Ep; 
	f[1] = Delta[1] * Ep - Delta[0] * Em * Svec[1] -
        I * Em * (Delta[2] * Svec[3] - Delta[3] * Svec[2]);
	f[2] = Delta[2] * Ep - Delta[0] * Em * Svec[2] -
        I * Em * (Delta[3] * Svec[1] - Delta[1] * Svec[3]);
	f[3] = Delta[3] * Ep - Delta[0] * Em * Svec[3] -
        I * Em * (Delta[1] * Svec[2] - Delta[2] * Svec[1]);

	ff[0] = conj(Delta[0]) * Ep; 
	ff[1] = conj(Delta[1]) * Ep + conj(Delta[0]) * Em * Svec[1] +
        I * Em * (conj(Delta[2]) * Svec[3] - conj(Delta[3]) * Svec[2]);
	ff[2] = conj(Delta[2]) * Ep + conj(Delta[0]) * Em * Svec[2] +
        I * Em * (conj(Delta[3]) * Svec[1] - conj(Delta[1]) * Svec[3]);
	ff[3] = conj(Delta[3]) * Ep + conj(Delta[0]) * Em * Svec[3] +
        I * Em * (conj(Delta[1]) * Svec[2] - conj(Delta[2]) * Svec[1]);
}


/*
 * gS_bulk(double complex*, double complex*, double complex*, double complex,
 *          double, double complex)
 * 
 * Singlet S,D-wave in magnetic field B.
 * 
 * Scalar parts in variable[0] and length of vector parts in variable[1],
 * remember that vector parts are directed along B
 * 
 * Args:
 *      double complex *g       ->
 *      double complex *f       ->
 *      double complex *ff      ->
 *      double complex en       ->
 *      double B                ->
 *      double complex Delta    ->
 * 
 * Returns:
 *      void
 */

void gS_bulk(double complex *g, double complex *f, double complex *ff, 
		    double complex en, double B, double complex Delta) /* (?) Is Delta pointer, like the other functions? */
{
	double Delta2;
	double complex eplus, eminus, Dplus, Dminus;
	
	Delta2 = pow(cabs(Delta), 2) + TINYsr;
	eplus = en - B;
	eminus = en + B;
	Dplus = csqrt(Delta2 - eplus * eplus);
	Dminus = csqrt(Delta2 - eminus * eminus);

	/* Bulk (homogeneous) propagator */
	g[0] = -1.0 / 2.0 * (eplus / Dplus + eminus / Dminus);
	g[1] = -1.0 / 2.0 * (eplus / Dplus - eminus / Dminus);
	g[2] = g[3] = 0.0;

	f[0] = 1.0 / 2.0 * Delta * (1 / Dplus + 1 / Dminus);
	f[1] = 1.0 / 2.0 * Delta * (1 / Dplus - 1 / Dminus);
	f[2] = f[3] = 0.0;

	ff[0] = 1.0 / 2.0 * conj(Delta) * (1 / Dplus + 1 / Dminus);
	ff[1] = 1.0 / 2.0 * conj(Delta) * (1 / Dplus - 1 / Dminus);
	ff[2] = ff[3] = 0.0;
}


/*
 * afromg(double complex*, double complex, double complex*, double complex*)
 * 
 * (?)
 * 
 * Args:
 *      double complex *ra      ->
 *      double complex gs       ->
 *      double complex *fv      ->
 *      double complex *ffv     ->
 * 
 * Returns:
 *      void
 */

void afromg(double complex *ra, double complex gs, double complex *fv,
            double complex *ffv)
{
	int j;
	double complex sf2, sf3, sf4;

	sf2 = sf3 = sf4 = 0.0;

	for (j = 1; j <= 3; j++) {
		sf2 += fv[j] * fv[j];
		sf3 += fv[j] * ffv[j];
		sf4 += ffv[j] * ffv[j];
	}

	for (j = 1; j <= 3; j++) {
		ra[j] = 
            ((fv[j] * (gs - I) - (fv[j] * sf3 - ffv[j] * sf2) / gs / 2.0) /
			(2.0 * gs * (gs - I) - sf3));
    }
}


/*
 * bfromg(double complex*, double complex, double complex*, double complex*)
 * 
 * (?)
 * 
 * Args:
 *      double complex *rb      ->
 *      double complex gs       ->
 *      double complex *fv      ->
 *      double complex *ffv     ->
 * 
 * Returns:
 *      void
 */

void bfromg(double complex *rb, double complex gs, 
		double complex *fv, double complex *ffv)
{
	int j;
	double complex sf2, sf3, sf4;

	sf2 = sf3 = sf4 = 0.0;

	for (j = 1; j <= 3; j++) {
		sf2 += fv[j] * fv[j];
		sf3 += fv[j] * ffv[j];
		sf4 += ffv[j] * ffv[j];
	}

	for (j=1; j<=3; j++) {
		rb[j] = 
            (-(ffv[j] * (gs - I) + (fv[j] * sf4 - ffv[j] * sf3) / gs / 2.0) /
			(2.0 * gs * (gs - I) - sf3));
    }
}


/*
 * gsfromric(double complex*, double complex*)
 * 
 * (?)
 * 
 * Args:
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      double complex          ->
 */

double complex gsfromric(double complex *ra, double complex *rb)
{
	int j;
	double complex saa, sbb, sab, gs;

	saa = sbb = sab = 0.0;

	for (j = 1; j <= 3; j++) {
		saa += ra[j] * ra[j];
		sbb += rb[j] * rb[j];
		sab += ra[j] * rb[j];
	}

	gs = -I * (1.0 - saa * sbb) / (1.0 + 2.0 * sab + saa * sbb);
	return gs;
}


/*
 * ffromric(double complex*, double complex*, double complex*, double complex*)
 * 
 * (?)
 * 
 * Args:
 *      double complex *fv      ->
 *      double complex *ffv     ->
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      void
 */

void ffromric(double complex *fv, double complex *ffv,
		    double complex *ra, double complex *rb)
{
	int j;
	double complex saa, sbb, sab, den;

	saa = sbb = sab = 0.0;
	for (j = 1; j <= 3; j++) {
		saa += ra[j] * ra[j];
		sbb += rb[j] * rb[j];
		sab += ra[j] * rb[j];
	}

	den = 1.0 + 2.0 * sab + saa * sbb;
	
	for (j = 1; j <= 3; j++) {
		fv[j] = -2.0 * I / den * (ra[j] + rb[j] * saa);
		ffv[j] = 2.0 * I / den * (rb[j] + ra[j] * sbb);
	}
}


/*
 * gfvfromric(double complex*, double complex*, double complex*,
 *          double complex*, double complex*)
 * 
 * (?)
 * 
 * Args:
 *      double complex *gv      ->
 *      double complex *fv      ->
 *      double complex *ffv     ->
 *      double complex *ra      ->
 *      double complex *rb      ->
 * 
 * Returns:
 *      void
 */

void gfvfromric(double complex *gv, double complex *fv, double complex *ffv,
		double complex *ra, double complex *rb)
{
	int j;
	double complex saa, sbb, sab, den;

	saa = sbb = sab = 0.0;

	for (j = 1; j <= 3; j++) {
		saa += ra[j] * ra[j];
		sbb += rb[j] * rb[j];
		sab += ra[j] * rb[j];
	}

	den = 1.0 + 2.0 * sab + saa * sbb;
	
	gv[1] = -2.0 / den * (ra[2] * rb[3] - ra[3] * rb[2]);
	gv[2] = -2.0 / den * (ra[3] * rb[1] - ra[1] * rb[3]);
	gv[3] = -2.0 / den * (ra[1] * rb[2] - ra[2] * rb[1]);

	for (j = 1; j <= 3; j++) {
		fv[j] = -2.0 * I / den * (ra[j] + rb[j] *saa);
		ffv[j] = 2.0 * I / den * (rb[j] + ra[j] *sbb);
	}
}