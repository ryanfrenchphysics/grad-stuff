#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include"nrutil.h"

#include"constvar.h"

/* No need to undef this */
#define TINYop  1.0e-20

/*
 * This contains methods for updating self-consistent order parameter 
 * and self-energies.
 */

/*
 * It requires the following variables to be included in 
 * the main program:
 * 
 * definitions are: 
 *	vector[v] - arranged column of all self-consistent fields,
 *		for example: 
 *		vector[v] = (Delta[0...NZ][1], Delta[0...NZ][2], ..., 
 *			     Sigma[0..NZ][1], Sigma[0..NZ][2], ...);
 *	Fvec[v] - vector of deviation of the self-consistency from zero,
 *	Fvec_next[v] - latest iteration outcome,
 *		for example, I always define for the OP part:
 * 		Fvec_next = Delta - 
 *			    1.0/log(t)*< T sum_{epsilon_m} (f-Delta/epsilon) >_{fs};
 *	Normf - norm of the deviation vector, 
 *		Normf = sum_v |Fvec[v]|^2; the self-consistency is reached when 
 *		Normf = 0;
 */


/* (?) Is this all necessary? */
/*	
 v = 1 ...... Nvector;
 #define Mk 10  - number of remembered previous steps

 int v, Nvector;
 double complex *vec_next, *Fvec_next;
 double complex **vector, **Fvec, *Fmin_prev;
 double *Normf; 


 Memory allocation and freeing:

        vec_next= Cvector(1, Nvector);
        Fvec_next=Cvector(1, Nvector);
        Fmin_prev= Cvector(1, Nvector); 
	Norm = dvector(1,Mk);
       	vector = Cmatrix(1,Mk, 1,Nvector); 
       	Fvec  = Cmatrix(1,Mk, 1,Nvector);
        
	Call to this subroutine in the main program has to be like this, 
       for example,	
	================================================================
	
	v=0;
	for (j=1; j<=3; j++)
		for (nz=0; nz<=NZ; nz++)
		{
			v++;
			vec_next[v] = Delta[nz][j];
			Fvec_next[v] = Delta[nz][j]-1.0/log(t)*Sum[nz][j];
		}
			
	OP_update(Nvector, vec_next, Fvec_next, vector, Fvec, Fmin_prev, Norm);

	v=0;
	for (j=1; j<=3; j++)
		for (nz=0; nz<=NZ; nz++)
		{
			v++;
			Deltanext[nz][j] = vec_next[v];
		}
	

	==================================================================
	
        free_Cvector(vec_next, 1, Nvector);
        free_Cvector(Fvec_next, 1, Nvector);
        free_Cvector(Fmin_prev, 1, Nvector);
	free_dvector(Norm, 1,Mk);
       	free_Cmatrix(vector, 1,Mk, 1,Nvector);
       	free_Cmatrix(Fvec, 1,Mk, 1,Nvector);
*/

/*
 * OP_update(int, double complex*, double complex*, double complex*[],
 *          double complex *[], double complex*, double[])
 * 
 * (?) Add description
 * 
 * Args:
 *      int Nvector                 -> (?)
 *      double complex *vec_next    ->
 *      double complex *Fvec_next   ->
 *      double complex *vector[]    ->
 *      double complex *Fvec[]      ->
 *      double complex *Fmin_prev   ->
 *      double Norm[]               ->
 * 
 * Returns:
 *      void
 */

void OP_update(int Nvector, double complex *vec_next, 
	double complex *Fvec_next, double complex *vector[], 
	double complex *Fvec[], double complex *Fmin_prev, double Norm[])
{
	int z_k, z_k0, z_l, z_v, mnz, mk_throw, *indx, zeros;
	double Normnext, Normaux, z_p, z_pp, z_q, divider, min;
	double complex **MM, *Mv, *g, det; 
	double complex *vec_min, *Fvec_min;
    /* (?) What do these mean? */
 
    vec_min = Cvector(1, Nvector);
    Fvec_min = Cvector(1, Nvector);

	/* Throwing away procedure */
	
	Normnext = 0.0;
	for (z_v = 1; z_v <= Nvector; z_v++)
		Normnext += pow(cabs(Fvec_next[z_v]), 2);

	Normnext = sqrt(Normnext);
	printf("OP update: Normnext=%e\n", Normnext);

	mnz = 1;
	for (z_k = 1; z_k <= Mk; z_k++){
		if (Norm[z_k] != 0.0)
            mnz++;
		/* printf("%d-%e ", k, Norm[k]); */
	}
	/* printf("\n"); */

    int if_is_done = 0;        

	if (mnz == Mk + 1) {
		mk_throw = 5;           /* number of vectors to throw away */

		if (Normnext < Norm[1] / 10.0) {
			for (z_k = Mk - mk_throw + 1; z_k <= Mk; z_k++) {
				Norm[z_k] = 0.0;
			}
			if_is_done = 1;
		}

		if (!if_is_done || Normnext > Norm[Mk]) {
			for (z_k = Mk - mk_throw + 1; z_k <= Mk; z_k++) {
				Norm[z_k] = 0.0;
			}
			if_is_done = 1;
		}

		if (!if_is_done || Normnext < Norm[Mk]) {	
			Norm[Mk] = 0.0;
			if_is_done = 1;
		}
	}	

    /* throwing away is over */
	mnz = 1;
	for (z_k = 1; z_k <= Mk; z_k++)
		if (Norm[z_k] != 0.0)
            mnz++;
	if (mnz == Mk + 1)
        printf("Throwing away didn't work!\n");

	/* add new element */
	Norm[mnz] = (Normnext != 0.0) ? Normnext : TINYop;

	for (z_v = 1; z_v <= Nvector; z_v++) {
		vector[mnz][z_v] = vec_min[z_v] = vec_next[z_v]; 
		Fvec[mnz][z_v] = Fvec_min[z_v] = Fvec_next[z_v];
	}

	if (mnz <= 4) {
        /*
         * change '==0' to '<= Nsi' (Nsi<5?) if you want to have Nsi
         * simple relaxations in the begining.
         */
        /* (?) What is this talking about? */

		z_p = 0.10;
		printf("Relaxation p=%.9f\n", z_p);
		goto next;          /* Skip all the following steps */
	}

	/* Solution for set of generators g[l], mnz-1 - number of eq-s */

	indx = ivector(1, mnz - 1);
    g = Cvector(1, mnz - 1);
    Mv = Cvector(1, mnz - 1);
    MM = Cmatrix(1, mnz - 1, 1, mnz - 1);
	
	divider = 0.0;
	zeros = 0;
	for (z_k = 1; z_k < mnz; z_k++) {
        for (z_l = 1; z_l < mnz; z_l++) {
            MM[z_k][z_l] = 0.0;
            for (z_v = 1; z_v <= Nvector; z_v++) {
                MM[z_k][z_l] += conj(Fvec[z_k][z_v] - Fvec[mnz][z_v]) *
                    (Fvec[z_l][z_v] - Fvec[mnz][z_v]);
            }
            if (cabs(MM[z_k][z_l]) != 0.0)
                divider += log10(cabs(MM[z_k][z_l]));
            if (cabs(MM[z_k][z_l]) == 0.0)
                zeros++;
        }
    }
	/* printf("divider=%e 1/power=%d \n", divider, (mnz-1)*(mnz-1)); */
	divider = divider/((mnz-1)*(mnz-1)-zeros);
	divider=pow(10,divider);
	/* printf("divider=%e \n", divider); */

	for (z_k = 1; z_k < mnz; z_k++) {
		Mv[z_k] = 0.0;
		for (z_v = 1; z_v <= Nvector; z_v++)
			Mv[z_k] += -Fvec[mnz][z_v] *
				conj(Fvec[z_k][z_v] - Fvec[mnz][z_v]);
	}

	min = 1e10;
	for (z_k = 1; z_k < mnz; z_k++) {
		Mv[z_k] /= divider; 
		if (isnan(cabs(Mv[z_k])))
            printf("Mv[%d]=%f\n", z_k, cabs(Mv[z_k]));

		for (z_l = 1; z_l < mnz; z_l++) {
			MM[z_k][z_l] /= divider; 
			if (isnan(cabs(MM[z_k][z_l]))) 
                printf("MM[%d][%d]=%f\n", z_k, z_l, cabs(MM[z_k][z_l]));
			if (cabs(MM[z_k][z_l]) < min && cabs(MM[z_k][z_l]) != 0.0)
                min = cabs(MM[z_k][z_l]);
		} 
	}

	/*
     * Now we have everything ready to find generators {g[l=1...mnz-1]}
	 * from equation MM[k][l]*g[l]=Mv[k] .
	 *
	 * use functions ludcmp_z and lubksb_z
	 */

	ludcmp_z(MM, mnz - 1, indx, &det);

	printf("mnz=%d, det = %e+I(%e)\n", mnz, creal(det), cimag(det));
	/* printf("zeros=%d, min = (%e)\n", zeros, min); */

	if (cabs(det) == 0.0) {
		z_p = 0.1 * rand() / (RAND_MAX) + 0.1;
		printf("Singular matrix: Random Relaxation p=%.9f\n", z_p);
		goto next;          /* Skip the following */
	}

	lubksb_z(MM, mnz - 1, indx, Mv);
	for (z_k = 1; z_k < mnz; z_k++) {
        g[z_k] = Mv[z_k];
		if (isnan(cabs(g[z_k])))
            printf("g[%d]=%f\n", z_k,cabs(g[z_k]));
    }
	
	for (z_v = 1; z_v <= Nvector; z_v++) {
		for (z_k = 1; z_k < mnz; z_k++) {
			Fvec_min[z_v] += g[z_k] *
				(Fvec[z_k][z_v] - Fvec[mnz][z_v]);
			vec_min[z_v] += g[z_k] *
				(vector[z_k][z_v] - vector[mnz][z_v]);
		}
	}

	Normaux = 0.0;
	for (z_v = 1; z_v <= Nvector; z_v++)
		Normaux += pow(cabs(vec_min[z_v] - vector[mnz][z_v]), 2);

	Normaux = sqrt(Normaux);

	z_pp = Normaux / Norm[mnz];
	if (isnan(z_pp))
        printf("pp=%f\n", z_pp);	
	Normaux = 0.0;

	for (z_v = 1; z_v <= Nvector; z_v++)
		Normaux += pow(cabs(Fmin_prev[z_v]), 2);
	Normaux = sqrt(Normaux);

	z_q = (1.0 < Normaux / Norm[mnz]) ? 1.0 : Normaux / Norm[mnz];
	if (isnan(z_q))
        printf("q=%f\n", z_q);

	z_p = z_pp * z_q;

	printf("Jump p=%.9e\n", z_p);
    free_ivector(indx, 1, mnz - 1);
    free_Cvector(g, 1, mnz - 1);
    free_Cvector(Mv, 1, mnz - 1);
    free_Cmatrix(MM, 1, mnz - 1, 1, mnz - 1);

next:
	/* Ordering of Norms and vectors */
	if (mnz != 1 && (Norm[mnz] < Norm[mnz - 1])) {
		z_k0 = mnz;
		for (z_k = 1; z_k < mnz; z_k++) {
			if (Norm[mnz] <= Norm[z_k]) {
                z_k0 = z_k;
                break;
            }
        }
		for (z_k = mnz - 1; z_k >= z_k0; z_k--) {
			Norm[z_k + 1] = Norm[z_k];
			for (z_v = 1; z_v <= Nvector; z_v++) {
				vector[z_k + 1][z_v] = vector[z_k][z_v];
				Fvec[z_k + 1][z_v] = Fvec[z_k][z_v];
			}
		}
		Norm[z_k0] = Normnext;
		for (z_v = 1; z_v <= Nvector; z_v++) {
			vector[z_k0][z_v] = vec_next[z_v];
			Fvec[z_k0][z_v] = Fvec_next[z_v];
		}
	}

	/* Next Order Parameter */
	for (z_v = 1; z_v <= Nvector; z_v++) {
		vec_next[z_v] = vec_min[z_v] + z_p * Fvec_min[z_v]; /* (?) Parentheses? */
		Fmin_prev[z_v] = Fvec_min[z_v];	
	}
	
    free_Cvector(vec_min, 1, Nvector);
    free_Cvector(Fvec_min, 1, Nvector);
}


/*
 * ludcmp_z(double complex**, int, int*, double complex*)
 * 
 * Given a matrix a[1..n][1..n], this routine replaces it by the LU 
 * decomposition of a rowwise permutation of itself. a and n are input. a 
 * is output, arranged as in equation NR(2.3.14) above; indx[1..n] is an 
 * output vector that records the row permutation effected by the 
 * partial pivoting; d is output as +/-1 depending on whether the number 
 * of row interchanges was even or odd, respectively. (I assign 
 * the value of determinant to it in the end. - A.V.) This routine is 
 * used in combination with lubksb_z to solve linear equations or 
 * invert a matrix
 * 
 * Args:
 *      double complex **a      ->
 *      int z_n                 ->
 *      int *indx               ->
 *      double complex *det     ->
 */

void ludcmp_z(double complex **a, int z_n, int *indx, double complex *det)
{ 
	int z_i, imax, z_j, z_k;
	double big, dum, temp; 
	double *z_vv;
	double complex sum, dum_z;
    /* (?) What do these represent? */

	z_vv = dvector(1, z_n);
	*det = 1.0;
	for (z_i = 1; z_i <= z_n; z_i++) {
		big = 0.0;
		for (z_j = 1; z_j <= z_n; z_j++) {
			if ((temp = cabs(a[z_i][z_j])) > big)
                big = temp;
        }
		/* print the matrix */
		// for (z_j = 1; z_j <= z_n; z_j++) {
        //    printf("%.2e\t", creal(a[z_i][z_j]));
        //     printf("\n");
        // }

		if (big == 0.0) {
			/* printf("Singular matrix!\n"); */
			*det = 0.0;
            return;
		}
		z_vv[z_i] = 1.0 / big;
	}

	imax = 1;
	for (z_j = 1; z_j <= z_n; z_j++) {
		for (z_i = 1; z_i < z_j; z_i++) {
			sum = a[z_i][z_j];
			for (z_k = 1; z_k < z_i; z_k++)
                sum -= a[z_i][z_k] * a[z_k][z_j];
			a[z_i][z_j] = sum;
		}
		big = 0.0;
		for (z_i = z_j; z_i <= z_n; z_i++) {
			sum = a[z_i][z_j];
			for (z_k = 1; z_k < z_j; z_k++)
				sum -= a[z_i][z_k] * a[z_k][z_j];
			a[z_i][z_j] = sum;
			if ((dum = z_vv[z_i] * cabs(sum)) >= big) {
				big = dum;
				imax = z_i;
			}
		}
		if (z_j != imax) {
			for (z_k = 1; z_k <= z_n; z_k++) {
				dum_z = a[imax][z_k];
				a[imax][z_k] = a[z_j][z_k];
				a[z_j][z_k] = dum_z;
			}
			*det = -(*det);
			z_vv[imax] = z_vv[z_j];
		}
		indx[z_j] = imax;
		if (cabs(a[z_j][z_j]) == 0.0)
            a[z_j][z_j] = TINYop;
		if (z_j != z_n) {
			dum_z = 1.0 / (a[z_j][z_j]);
			for (z_i = z_j + 1; z_i <= z_n; z_i++)
                a[z_i][z_j] *= dum_z;
		}
	}
	// determinant calculation
	for (z_j = 1; z_j <= z_n; z_j++)
        *det *= a[z_j][z_j];

	free_dvector(z_vv, 1, z_n);
}


/*
 * lubksb_z(double complex**, int, int*, double complex[])
 * 
 * Solves the set of n linear equations A�X = B. Here a[1..n][1..n] is input, 
 * not as the matrix A but rather as its LU decomposition, determined by the
 * routine ludcmp_z. indx[1..n] is input as the permutation vector returned by
 * ludcmp_z. b[1..n] is input as the right-hand side vector B, and returns
 * with the solution vector X. a, n, and indx are not modifed by this routine 
 * and can be left in place for successive calls with different right-hand
 * sides b. This routine takes into account the possibility that b will begin 
 * with many zero elements, so it is efficient to use in matrix inversion.
 * 
 * args:
 *      double complex **z_a    -> (?)
 *      int n                   ->
 *      int *indx               ->
 *      double complex z_b[]    ->
 * 
 * returns:
 *      void
 */

void lubksb_z(double complex **z_a, int n, int *indx, double complex z_b[])
{
	int z_i, z_ip, z_j;
    int z_ii = 0;
	double complex sum;

	for (z_i = 1; z_i <= n; z_i++) {
		z_ip = indx[z_i];
		sum = z_b[z_ip];
		z_b[z_ip] = z_b[z_i];
		if (z_ii) {
			for (z_j = z_ii; z_j <= z_i - 1; z_j++)
                sum -= z_a[z_i][z_j] * z_b[z_j];
        } else if (sum) {
            z_ii = z_i;
        }
		z_b[z_i] = sum;
	}
	for (z_i = n; z_i >= 1; z_i--) {
		sum = z_b[z_i];
		for (z_j = z_i + 1; z_j <= n; z_j++)
            sum -= z_a[z_i][z_j] * z_b[z_j];
		z_b[z_i] = sum / z_a[z_i][z_i];
	}
}