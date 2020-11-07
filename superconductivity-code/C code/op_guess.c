#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "op_guess.h"
#include "mpi_timer.h"
#include "nrutil.h"
#include "ricparam.h"
#include "subroutines.h"
#include "integrators.h"
#include "constvar.h"

#ifndef M_PI
    /* In case it isn't imported from math.h */
    #define M_PI 3.14159265358979323846
#endif

#define NZinit  0
#define Nza     0
#define Nzb     NZ

/* Functions used to calculate initial guess for the gap in the film */


/*
 * op_guess(double, double, double*, double*, double*, double**, double**,
 *          double**, double**, double**, int, int)
 * 
 * (?) Insert description of this function
 * 
 * Args:
 *      double dfilm    -> (?)
 *      double z_tt     -> (?)
 *      double *xgr     ->  (?)
 *      double *zgr     ->  (?)
 *      double **A11        ->
 *      double **A22        ->
 *      double **A33        ->
 *      double **A13        ->
 *      double **A31        ->
 *      int rescale_op      ->
 *      int read_op_film    ->
 * 
 * Returns:
 *      void
 */
void op_guess(double dfilm, double z_tt, double *xgr, double *zgr, 
	    double **A11, double **A22, double **A33, double **A13, double **A31,
        int rescale_op, int read_op_film)
{
	int nx, nz;
	double x0, gap, delta1, delta2, delta3, delta13, delta31, xx, zz;
    double scale11, scale22, scale33;
	double complex **Delta;
	char string[LL];
	FILE *outgap, *grid;
    /* (?) What do all of these variables mean? */


	Delta = Cmatrix(-NZ, NZ, 1, 3);
	gap = bulkgap(z_tt, 'A');

	/* initial guess for the 'far from domain-wall' OP */
	for (nz = -NZ; nz <= NZ; nz++) {
		Delta[nz][2] = Delta[nz][1] = gap;
		Delta[nz][3] = gap * cos(M_PI * zgr[nz] / dfilm); 
	}
	grid = fopen("xzgrid.dat", "w");

	for (nz = -NNZ; nz <= NNZ; nz++)
        for (nx = -NNX; nx <= NNX; nx++)
            fprintf(grid, "%f %f\n", xgr[nx], zgr[nz]);

	fclose(grid);

	switch (read_op_film) {
        case 0:
            gapBfilm_specular(dfilm, z_tt, Delta, zgr);

            sprintf(string, "gapBspec_D%.3f_t%.3f.dat", dfilm, z_tt);
            outgap = fopen(string, "w");
            for (nz = -NZ; nz <= NZ; nz++)
                fprintf(outgap, "%f \t%e \t%e \t%e\n", zgr[nz], 
                        creal(Delta[nz][1]), creal(Delta[nz][2]),
                        creal(Delta[nz][3]));
            fclose(outgap);
            break;

        case 1:
            sprintf(string, "gapBspec_D%.3f_t%.3f.dat", dfilm, z_tt);
            printf("reading OP from <%s>\n", string);
            outgap = fopen(string, "r");
            for (nz = -NZ; nz <= NZ; nz++) {
                fgets(string, LL, outgap); 
                sscanf(string, "%*s %lf %lf %lf", &delta1, &delta2, &delta3); 
                Delta[nz][1] = delta1;
                Delta[nz][2] = delta2;
                Delta[nz][3] = delta3;
            }
            fclose(outgap);
            break;

        case 2:
            sprintf(string, "gapFilm_t%.3f_D%.3f_SDW.dat", z_tt, dfilm);
            printf("reading OP from <%s>\n", string);
            outgap = fopen(string, "r");
            
            for (nx = -NNX; nx <= NNX; nx++) {
                for (nz = -NZ; nz <= NZ; nz++) {
                    fgets(string, LL, outgap); 
                    sscanf(string, "%lf %lf %lf %lf %lf %lf %lf", 
                        &xx, &zz, &delta1, &delta2, &delta3, &delta13,
                        &delta31); 

                    /*
                     * fprintf(outgap, "%f\t%f \t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
                     * xgr0[nx], zgr0[nz], A11[nx][nz], A22[nx][nz], A33[nx][nz],
                     *      A13[nx][nz], A31[nx][nz]);
                     */

                    if (fabs(xgr[nx] - xx) > 1e-10)
                        printf("Read xgr[%d] != set up xgr!\n", nx); 
                    if (fabs(zgr[nz] - zz) > 1e-10)
                        printf("Read zgr[%d] != set up zgr!\n", nz);

                    A11[nx][nz] = delta1;
                    A22[nx][nz] = delta2;
                    A33[nx][nz] = delta3;
                    A13[nx][nz] = delta13;
                    A31[nx][nz] = delta31;
                }
            }
		    fclose(outgap);
		    return;     /* (?) Is return intentional? It won't run the rest. */
	}

	x0 = 5.0;
	if (!rescale_op) {
		for (nx = NNX; nx >= -NNX; nx--) {
            for (nz = NZ; nz >= -NZ; nz--) {
                A11[nx][nz] = creal(Delta[nz][1]);
                A22[nx][nz] = creal(Delta[nz][1]); 
                A33[nx][nz] = creal(Delta[nz][3])*tanh(xgr[nx]/x0); 
                A13[nx][nz] = (0.2 * (Delta[0][3]) *
                    sin(2 * M_PI * zgr[nz] / dfilm) * 
                    tanh(xgr[nx] / x0) / cosh(xgr[nx] / x0)); 
                /* A13[nx][nz]=0.0; */
                A31[nx][nz] = -(Delta[0][3]) * 
                    sin(M_PI * zgr[nz] / dfilm) / cosh(xgr[nx] / x0); 
                /* A31[nx][nz]=0.0; */

                /*
                 * The '-' sign for 31. It prefers it for some reason.
                 * Even if I start first with a '+' it changes to '-' after
                 * several iterations. There is certain 'chirality' in a sense.
                 */ 
            }
        }
	} else {
		for (nz = NZ; nz >= -NZ; nz--) {
			scale11 = cabs(Delta[nz][1]) / (fabs(A11[NNX][nz]) + abs_prec);
			scale22 = cabs(Delta[nz][2]) / (fabs(A22[NNX][nz]) + abs_prec);
			scale33 = cabs(Delta[nz][3]) / (fabs(A33[NNX][nz]) + abs_prec);

            for (nx = NNX; nx >= -NNX; nx--) {
                A11[nx][nz] *= scale11; 
                A22[nx][nz] *= scale22;
                A33[nx][nz] *= scale33;
                /* leave the same */ 
                /* A13[nx][nz] */
                /* A31[nx][nz] */
            }
		}
	}
	free_Cmatrix(Delta, -NZ, NZ, 1, 3);
	return;
}

/* (?) Does this comment go with the next function? */
/** ******************************************************************
*
* 	GAP for He3-B in the film with 
* 		specular wall and specular free surface
*
* ********************************************************************
* 	
*  OP has the following structure \hat{Delta}=\vec{d} i\vec{\sigma}\sigma_y
*  with \vec{d} given by:
*  ( Delta[1]*p_x, Delta[2]*p_y, Delta[3]*p_z ) 
*
* ******************************************************************** */

/*
 * gapBfilm_specular(double, double, double complex**, double*)
 * 
 * This is a sub-routine for calculation of the Green Function and the 
 * OP far from the domain wall. Green functions and Ricatti amplitudes are 
 * in most general form (zero component of Ric_amp is for singlet case) but 
 * the integrator is for UNITARY case only (for speed). So one has to 
 * remember to put Ric_amp[0]=0 in some places.
 * 
 * Args:
 *      double dfilm                -> (?)
 *      double z_tt                 ->
 *      double complex **Delta      ->
 *      double *zgr                 ->
 * 
 * Returns:
 *      void
 */

void gapBfilm_specular(double dfilm, double z_tt, double complex **Delta,
        double *zgr)
{
	time_t time_of_start;
	
	int done, z_M, nz, nzerr, pan, panm, z_k, z_m, z_j, z_v, iter;
	double theta, dtheta, thetam, z_dz, z_ds, z_em;  
	double rerr, err1, err, err3; 
	double complex **sum, **Sum, **D, z_en; 

	double complex **Deltanext;
	double complex *Deltai, *Deltaf;
    double complex dum[4]={0,0,0,0}, Deltac[4];
	double complex **a[2*PAN], **b[2*PAN];
	double complex *aprev[2*PAN], *bprev[2*PAN];
	double complex **g[2*PAN], **f[2*PAN], **ff[2*PAN]; 
	/* double complex *f0, *ff0, *g0, *a0, *b0; //Ovchin BC variables */

    double complex *vec_next, *Fvec_next, *Fmin_prev;
    double *Normf; 
    double complex **vector, **Fvec;
	FILE *ric, *gap;
    int Nvector = 3 * (2 * NZ + 1);
    /* (?) What do all of these variables mean? */

	/* allocation of memory for the arrays */
	sum = Cmatrix(-NZ, NZ, 1, 3);
	Sum = Cmatrix(-NZ, NZ, 1, 3);
	D = Cmatrix(-NZ, NZ, 1, 3);
	Deltanext = Cmatrix(-NZ, NZ, 1, 3);
	Deltai = Cvector(1, 3);
	Deltaf = Cvector(1, 3);
	
	for (pan = 0; pan < 2 * PAN; pan++) {
		a[pan] = Cmatrix(-NZ, NZ, 0, 3);
		b[pan] = Cmatrix(-NZ, NZ, 0, 3);
		aprev[pan] = Cvector(0, 3);
		bprev[pan] = Cvector(0, 3);
		g[pan] = Cmatrix(-NZ, NZ, 0, 3);
		f[pan] = Cmatrix(-NZ, NZ, 0, 3);
		ff[pan] = Cmatrix(-NZ, NZ, 0, 3);
	}
	
	/*a0 = Cvector(0,3);
	b0 = Cvector(0,3);
	g0=Cvector(0,3);
	f0=Cvector(0,3);
	ff0=Cvector(0,3);*/

    vec_next = Cvector(1, Nvector); 
    Fvec_next = Cvector(1, Nvector); 
    Fmin_prev = Cvector(1, Nvector); 
	Normf = dvector(1, Mk);
    vector = Cmatrix(1, Mk, 1, Nvector); 
    Fvec  = Cmatrix(1, Mk, 1, Nvector);

	time_of_start = starttime();

	dtheta = M_PI / 2 / PAN;

	z_M = Mats(z_tt);      /* M - cutoff for Matsubara sum */
	
	for (z_k = 1; z_k <= Mk; z_k++)
        Normf[z_k] = 0.0;
	
	iter = done = 0;
	do {
        /* iterations for the OP begin */
		iter++;

		for (z_j = 1; z_j <= 3; z_j++)
		    for (nz = -NZ; nz <= NZ; nz++)
			    Sum[nz][z_j] = 0.0;

		gap = fopen("gap_current.dat", "w");
		for (nz = -NZ; nz <= NZ; nz++) 
		    fprintf(gap, "%f %e %e %e\n", zgr[nz], 
		        creal(Delta[nz][1]), creal(Delta[nz][2]), creal(Delta[nz][3]));
		fclose(gap);

	/* Start doing sum over Matsubara's, M - cutoff */
	for (z_m = 0; z_m <= z_M; z_m++) {
		z_em = z_tt * (z_m + 0.5);

		/* calculate initial b_out and a_in */
		for (pan = 0; pan < PAN; pan++) {
			panm = 2 * PAN - 1 - pan;
			theta = dtheta * (pan + 0.5);
			thetam = dtheta * (panm + 0.5);

			Deltac[1] = Delta[NZinit][1] * sin(theta);
			Deltac[2] = 0.0;
			Deltac[3] = Delta[NZinit][3] * cos(theta);
			z_en = I * z_em;

			/*
             * g_bulk(g[pan][NZinit], f[pan][NZinit], ff[pan][NZinit], en, Deltai);
             * a_from_g(a[pan][NZinit], g[pan][NZinit], f[pan][NZinit], ff[pan][NZinit]);
             */

			ab_bulk_non_unitary(a[pan][NZinit], dum, z_en, Deltac);

			for (nz = NZinit; nz < NZ; nz++) {
				z_en = I * z_em;
				Deltai[1] = Delta[nz][1] * sin(theta);
				Deltai[2] = 0.0;
				Deltai[3] = Delta[nz][3] * cos(theta);
				Deltaf[1] = Delta[nz + 1][1] * sin(theta);
				Deltaf[2] = 0.0;
				Deltaf[3] = Delta[nz + 1][3] * cos(theta);
				z_dz = zgr[nz + 1] - zgr[nz];
				z_ds = z_dz / cos(theta);
					
				step_integratorc_a(z_ds, z_en, dum, Deltaf, a[pan][nz],
                    a[pan][nz + 1]);

				/*
                 * step_integratorc_a(ds, en, Deltai, Deltaf, a[pan][nz], a[pan][nz+1]);
                 * step_integratorm_a(ds, en, Deltai, Deltaf, a[pan][nz], a[pan][nz+1]);
                 */
			}

			for (z_j = 0; z_j <= 3; z_j++)
				a[panm][NZ][z_j] = a[pan][NZ][z_j];

			b[pan][NZ][0] = conj(a[panm][NZ][0]);
			b[pan][NZ][1] = conj(a[panm][NZ][1]);
			b[pan][NZ][2] = conj(a[panm][NZ][2]);
			b[pan][NZ][3] = -conj(a[panm][NZ][3]);
		    } 

		do {
			for (pan = 0; pan < PAN; pan++) {
                    for (z_j = 0; z_j <= 3; z_j++) {
                        bprev[pan][z_j] = b[pan][Nzb][z_j];
                        aprev[pan][z_j] = a[pan][Nza][z_j];
                    }
            }
			
			for (pan = 0; pan < PAN; pan++) {
				panm = 2 * PAN - 1 - pan;
				theta = dtheta * (pan + 0.5);
				thetam = M_PI - theta;
	
				for (nz = NZ; nz > -NZ; nz--) {
                    z_dz = zgr[nz] - zgr[nz - 1];
                    z_ds = z_dz / cos(theta);
					
                    z_en = I * z_em;
                        
                    /* find b_out[nz = -NZ..NZ] */
                    Deltai[1] = Delta[nz][1] * sin(theta);
                    Deltai[2] = 0.0;
                    Deltai[3] = Delta[nz][3] * cos(theta);
                    Deltaf[1] = Delta[nz - 1][1] * sin(theta);
                    Deltaf[2] = 0.0;
                    Deltaf[3] = Delta[nz - 1][3] * cos(theta);
                    
                    step_integratorc_b(z_ds, z_en, dum, Deltaf, b[pan][nz],
                        b[pan][nz - 1]);
                    /* step_integratorc_b(ds, en, Deltai, Deltaf, b[pan][nz],
                     *      b[pan][nz-1]);
                     */
        
                    /* Find a_in[nz = -NZ...NZ] */
                    Deltai[1] = Delta[nz][1] * sin(thetam);
                    Deltai[2] = 0.0;
                    Deltai[3] = Delta[nz][3] * cos(thetam);
                    Deltaf[1] = Delta[nz - 1][1] * sin(thetam);
                    Deltaf[2] = 0.0;
                    Deltaf[3] = Delta[nz - 1][3] * cos(thetam);
                    
                    step_integratorc_a(z_ds, z_en, dum, Deltaf, a[panm][nz],
                        a[panm][nz-1]);
                    /*
                     * step_integratorc_a(ds, en, Deltai, Deltaf, a[panm][nz],
                     *   a[panm][nz-1]);
                     */
				}
				
				/* assume a_out[nz = -NZ] */
				for (z_j = 0; z_j <= 3; z_j++)
					a[pan][-NZ][z_j] = a[panm][-NZ][z_j];
			}

			/* calculate a_out[nz=-NZ...NZ] */
			for (pan = 0; pan < PAN; pan++) {
				theta = dtheta * (pan + 0.5);
	
				for (nz = -NZ; nz < NZ; nz++) {
        
                    z_dz = zgr[nz + 1] - zgr[nz];
                    z_ds = z_dz / cos(theta);
                    z_en = I * z_em;
                    Deltai[1] = Delta[nz][1] * sin(theta);
                    Deltai[2] = 0.0;
                    Deltai[3] = Delta[nz][3] * cos(theta);
                    Deltaf[1] = Delta[nz + 1][1] * sin(theta);
                    Deltaf[2] = 0.0;
                    Deltaf[3] = Delta[nz + 1][3] * cos(theta);
                        
                    step_integratorc_a(z_ds, z_en, dum, Deltaf, a[pan][nz],
                        a[pan][nz+1]);
                    /*
                     * step_integratorc_a(ds, en, Deltai, Deltaf, a[pan][nz], a[pan][nz+1]);
                     * step_integratorm_a(ds, en, Deltai, Deltaf, a[pan][nz], a[pan][nz+1]);
                     */
				}
	
				panm = 2 * PAN - 1 - pan;
				for (z_j = 0; z_j <= 3; z_j++) 
					a[panm][NZ][z_j] = a[pan][NZ][z_j];

				b[pan][NZ][0] = conj(a[panm][NZ][0]);
				b[pan][NZ][1] = conj(a[panm][NZ][1]);
				b[pan][NZ][2] = conj(a[panm][NZ][2]);
				b[pan][NZ][3] = -conj(a[panm][NZ][3]);
			} /* pan = 0; an < PAN; pan++ */
		
			err = 0.0;
				
			for (pan = 0; pan < PAN; pan++) {
                for (z_j = 0; z_j <= 3; z_j++) {
                    rerr = cabs(b[pan][Nzb][z_j] - bprev[pan][z_j]) / 
                        (cabs(b[pan][Nzb][z_j]) + abs_prec);
                    if (err < rerr)
                        err = rerr;
                    rerr = cabs(a[pan][Nza][z_j] - aprev[pan][z_j]) /
                        (cabs(a[pan][Nza][z_j]) + abs_prec);
                    if (err < rerr)
                        err = rerr;
                }
            }
		} while (err > 0.1 * prec); 
		
		/* 
         * Find phys propagator on the z-axis.
         * Average over Fermi surface
         */
		for (nz = -NZ; nz <= NZ; nz++)
			sum[nz][1] = sum[nz][2] = sum[nz][3] = 0.0;
		for (pan = 0; pan < PAN; pan++) {
			theta = dtheta * (pan + 0.5); 
			
			for (nz = -NZ; nz <= NZ; nz++) {
                /* Make sure there is only the triplet part */
				a[pan][nz][0] = b[pan][nz][0] = 0.0;
				g_from_ric(g[pan][nz], a[pan][nz], b[pan][nz]);
				f_from_ric(f[pan][nz], ff[pan][nz], a[pan][nz], b[pan][nz]);

				/* if needed - below are the Green's function compnents */
				/*g00[pan][nz][m] = g[pan][nz][0];
				f01[pan][nz][m]= f[pan][nz][1];
				f03[pan][nz][m]= f[pan][nz][3];
				ff01[pan][nz][m]= ff[pan][nz][1];
				ff03[pan][nz][m]= ff[pan][nz][3];*/
			}
			
			/* Average over Fermi surface for a certain 
			 * matsubara energy */
			for (nz = -NZ; nz <= NZ; nz++) {
				sum[nz][1] += 0.5 * dtheta * sin(theta) *
		            (f[pan][nz][1] + conj(ff[pan][nz][1]) -
				    2.0 * Delta[nz][1] * sin(theta) / z_em) * sin(theta);
				sum[nz][3] += dtheta * sin(theta) *
				    (f[pan][nz][3] + conj(ff[pan][nz][3]) -
				    2.0 * Delta[nz][3] * cos(theta) / z_em) * cos(theta);
			}
		}

		for (nz = -NZ; nz <= NZ; nz++) {
			Sum[nz][1] += 3.0 * z_tt / 2 * sum[nz][1];
			Sum[nz][3] += 3.0 * z_tt / 2 * sum[nz][3];
			Sum[nz][2] = Sum[nz][1];
		}
			
		if (z_m==0) {
			ric = fopen("ric.check.dat", "w");
			pan = 20;
			for (nz = -NZ; nz <= NZ;nz++)
                fprintf(ric, "%f %f %f %f %f\n", z_dz * nz,
                cimag(a[pan][nz][1]), cimag(a[pan][nz][3]),
                cimag(b[pan][nz][1]), cimag(b[pan][nz][3]));
			fclose(ric);
		}
	}

	z_v = 0;
	for (z_j = 1; z_j <= 3; z_j++) {
        for (nz = -NZ; nz <= NZ; nz++) {
                z_v++;
                vec_next[z_v] = Delta[nz][z_j];
                Fvec_next[z_v] = Delta[nz][z_j] - 1.0 / log(z_tt) *
                    Sum[nz][z_j];
            }
    }

	/* 
     * After 40 iterations do simple relaxation again, for this just make
	 * Normf[k]=0
     */

	if (iter < 12)
        for (z_j = 2; z_j < Mk; z_j++)
            Normf[z_j] = 0.0;
	if (iter % 40 == 0)
        for (z_j = 2; z_j < Mk; z_j++)
            Normf[z_j] = 0.0;
			
	OP_update(Nvector, vec_next, Fvec_next, vector, Fvec, Fmin_prev, Normf);

	z_v = 0;
	for (z_j = 1; z_j <= 3; z_j++) {
        for (nz = -NZ; nz <= NZ; nz++) {
            z_v++;
            Deltanext[nz][z_j] = vec_next[z_v];
        }
    }

	err = err1 = err3 = 0.0;
    nzerr = 0;

	for (nz = -NZ; nz <= NZ; nz++) {
		err = cabs(Deltanext[nz][1] - Delta[nz][1]) /
            (cabs(Delta[nz][1]) + abs_prec);

		if (err > err1) {
            err1 = err;
            nzerr = nz;
            err3 = 1.0;
        }
		err = cabs(Deltanext[nz][3] - Delta[nz][3]) / 
            (cabs(Delta[nz][3]) + abs_prec);
		if(err > err1) {
            err1 = err;
            nzerr = nz;
            err3 = 3.0;
        }
	}

	if (err1 < prec)
        done=1;

	for (nz = -NZ; nz <= NZ; nz++) {
		Delta[nz][1] = Deltanext[nz][1];
		Delta[nz][2] = Deltanext[nz][1];
		Delta[nz][3] = Deltanext[nz][3];
	}
	
	printf("%d iteration is over. Error = %e (nz=%d in component %.0f ) \n \n", 
			iter, err1, nzerr, err3); 
	
	} while (!done); /* check if we need to continue iterations */

	printf("Initial guess OP is determined for t=%f D=%f.\n", z_tt, dfilm);
	
	/* free_dvector(nu, -NZ,NZ);
	free_dvector(D_nu, -NZ,NZ); */
	free_Cmatrix(sum, -NZ, NZ, 1, 3);
	free_Cmatrix(Sum, -NZ, NZ, 1, 3);
	free_Cmatrix(D, -NZ, NZ, 1, 3);
	free_Cmatrix(Deltanext, -NZ, NZ, 1, 3);
	free_Cvector(Deltai, 1, 3);
	free_Cvector(Deltaf, 1, 3);
	
	for (pan = 0; pan < 2 * PAN; pan++) {
		free_Cmatrix(a[pan], -NZ, NZ, 0, 3);
		free_Cmatrix(b[pan], -NZ, NZ, 0, 3);
		free_Cvector(aprev[pan], 0, 3);
		free_Cvector(bprev[pan], 0, 3);
		free_Cmatrix(g[pan], -NZ, NZ, 0, 3);
		free_Cmatrix(f[pan], -NZ, NZ, 0, 3);
		free_Cmatrix(ff[pan], -NZ, NZ, 0, 3);
	}
	
	/*free_Cvector(a0, 0,3);
	free_Cvector(b0, 0,3);
	free_Cvector(g0, 0,3);
	free_Cvector(f0, 0,3);
	free_Cvector(ff0, 0,3);*/
		
    free_Cvector(vec_next, 1, Nvector);
    free_Cvector(Fvec_next, 1, Nvector);
    free_Cvector(Fmin_prev, 1, Nvector);
	free_dvector(Normf, 1, Mk);
    free_Cmatrix(vector, 1, Mk, 1, Nvector);
    free_Cmatrix(Fvec, 1, Mk, 1, Nvector);	
	endtime(time_of_start);
}