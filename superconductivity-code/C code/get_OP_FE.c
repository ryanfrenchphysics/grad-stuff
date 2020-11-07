#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "get_OP_FE.h"
#include "FSsphere_lebedev.h"
#include "nrutil.h"
#include "ricparam.h"
#include "mpi_timer.h"
#include "subroutines.h"
#include "constvar.h"

#ifndef M_PI
    /* In case it isn't imported from math.h */
    #define M_PI 3.14159265358979323846
#endif


/* 
 * (?) Find a better name for these variables,
 * or explanations
 */
static 	double **A11, **A13, **A22, **A31, **A33, D, t;
static 	double **sum11, **sum13, **sum22, **sum31, **sum33;
static double complex **Deltai, **Deltaf, *f, *ff, ***a, ***aa;
/* static double complex **a_tr, **aa_tr; */
static double complex **ampl_prev, **a1_bc, **aa1_bc, **a2_bc, **aa2_bc;
static double *zgr, *xgr;

/*
 * static Fermi surface variables that are set once after
 * the input.dat file is read (for each node)
 */
static double *dOmg, *Vx, *Vy, *Vz, **YY, **fsY2;

static void master(int AN, int self_cnst);
static void slave(int AN);

/*
 * parameters_on_FS(int, double complex **)
 * 
 * set up OP parameters, etc, for FS points
 * 
 * Args:
 *      int AN  ->  TODO: description
 *      double complex **g_dos  ->  info how many trajectories to calculate,
 *                                  resolved dos for Ang = creal(g_dos[0][0])
 * 
 * Returns:
 *      void
 * 
 * TODO: Refine function description to describe vvv
 *
 * creal(g_dos[0][m=1..Ang])=angles[m] - angles for trajectories 
 * (in degrees, measured from the surface normal in normal(x)-superflow(y)-plane)
 * on output, cimag(g_dos[0][m])='an_dos' the number of FS point, closest to angle[m] 
 * angle[m] is also updated
 * 
 * ATTENTION: g_dos is just a place holder. Consult ~/qccalcul/pgap/slab/srcNonUnit/get_OPetc.c 
 * on how to set it up
 * 
 */

void parameters_on_FS(int AN, double complex **g_dos)
{
    int an;
	double px, py, pz, angle;

	for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
            fsY2[a][b]=0.0;

	for (an = 1; an <= AN; an++) {
		px = Vx[an]; py=Vy[an]; pz=Vz[an]; 
		YY[1][an] = YY1(px,py,pz); 
		YY[2][an] = YY2(px,py,pz); 
		YY[3][an] = YY3(px,py,pz);

		for (int a = 1; a <= 3; a++)
            for(int b = 1; b <= 3; b++)
                fsY2[a][b] += YY[a][an] * YY[b][an] * dOmg[an]; 
	}
}


/*
 * FSsetup(int *, double complex **)
 * 
 * the 2D sphere of the FS is put into a 1-d array 1...AN
 * On input, ANFS is approximate number of angles over 4Pi sphere.
 * On output, ANFS = AN is actual number of points over 4Pi sphere.
 * 
 * Args:
 *      int *ANFS   ->  TODO: add description
 *      double complex **g_dos  ->  TODO: add description
 * 
 * Returns:
 *      void
 */

void FSsetup(int *ANFS, double complex **g_dos)
{
	int i,j, AN, order; 

	AN = *ANFS;
	FSgrid_numbers(&order, &AN); 

	/* set up the static FS parameters */
    Vx = dvector(1, AN); 
    Vy = dvector(1, AN); 
    Vz = dvector(1, AN); 
    dOmg = dvector(1, AN); 
	fsY2 = dmatrix(1, 3, 1, 3);
	YY = dmatrix(1, 3, 1, AN); 

	/* spherical coords: (x,y,z) -> (Vx,Vy,Vz) */
	FSlebedev_map(Vx, Vy, Vz, dOmg, order, AN); 

	parameters_on_FS(AN, g_dos);

	/*
     * printf(" The <Y[a]*Y[b]>  matrix: \n");
     * for (int i = 1; i <= 3; i++)
     *      for(int j = 1; j <= 3; j++)
     *          printf("\t%.10f /n", fsY2[i][j]);
     */

	*ANFS = AN;
}

/*
 * domain_OP(int, int, double, double, double *, double *, double **,
 *          double **, double **, double **, double **, int, char **)
 * 
 * TODO: Description
 * 
 * Args:
 *      int self_cnst    ->
 *      int AN          ->
 *	    double temp     ->
 *	    double Dthick   ->
 *	    double *xgrid   ->
 *	    double *zgrid   ->
 *	    double **OP11   ->
 *	    double **OP22   ->
 *	    double **OP33   ->
 *	    double **OP13   ->
 *	    double **OP31   ->
 *	    int argc        ->
 *	    char *argv[]    ->
 *
 * Returns:
 *      void
 */

void domain_OP(int self_cnst, int AN, double temp, double Dthick,
        double *xgrid, double *zgrid, double **OP11, double **OP22,
        double **OP33, double **OP13, double **OP31, int argc, char *argv[])
{
    /* (?) Do we need to include these for MPI's sake? */
	/*
     * int self_cnst, AN;
     * double temp, Dthick; 
	 * double *xgrid, *zgrid;
     * double **OP11, **OP22, **OP33, **OP13, **OP31;
	 * int argc;
	 * char *argv[];
     */

	int myrank; 

	f = Cvector(0,3);
	ff = Cvector(0,3);
	a = C3tensor(-NNX, NNX, -NNZ, NNZ, 0, 3);
	aa = C3tensor(-NNX, NNX, -NNZ, NNZ, 0, 3);
	ampl_prev = Cmatrix(-NMAX, NMAX, 0, 3);

	Deltai = Cmatrix(1, 2*NMAX, 0, 3);
	Deltaf = Cmatrix(1, 2*NMAX, 0, 3);
	a1_bc = Cmatrix(-NMAX, NMAX, 0, 3);
	aa1_bc = Cmatrix(-NMAX, NMAX, 0, 3);
	a2_bc = Cmatrix(-NMAX, NMAX, 0, 3);
	aa2_bc = Cmatrix(-NMAX, NMAX, 0, 3);
	
	sum11 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	sum13 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	sum22 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	sum31 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	sum33 = dmatrix(-NNX, NNX, -NNZ, NNZ);

	A11 = OP11;
	A13 = OP13;
	A22 = OP22;
	A31 = OP31;
	A33 = OP33;

	xgr = xgrid;
	zgr = zgrid;

	D = Dthick;
	t = temp;

	MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("rank %d received:\tt = %f D = %f \n", myrank, temp, Dthick);
	fflush(stdout);

	if (myrank == 0)
        master(AN, self_cnst);
	else
        slave(AN);

	free_Cvector(f, 0,3);
	free_Cvector(ff, 0,3);
	free_C3tensor(a, -NNX, NNX, -NNZ, NNZ, 0, 3);
	free_C3tensor(aa, -NNX, NNX, -NNZ, NNZ, 0, 3);
	free_Cmatrix(ampl_prev, -NMAX, NMAX, 0, 3);

	free_Cmatrix(Deltai, 1, 2*NMAX, 0, 3);
	free_Cmatrix(Deltaf, 1, 2*NMAX, 0, 3);
	free_Cmatrix(a1_bc, -NMAX, NMAX, 0, 3);
	free_Cmatrix(aa1_bc, -NMAX, NMAX, 0, 3);
	free_Cmatrix(a2_bc, -NMAX, NMAX, 0, 3);
	free_Cmatrix(aa2_bc, -NMAX, NMAX, 0, 3);

	free_dmatrix(sum11, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(sum13, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(sum22, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(sum31, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(sum33, -NNX, NNX, -NNZ, NNZ);
    /* (?) Do we need to return int: MPI error? */
}


/*
 * next(double *, double, double *, char *, char *)
 * 
 * TODO: Description
 * 
 * Args:
 *      double *z       ->
 *      double znext    ->
 *      double *err     ->
 *      char *error     ->
 *      char *where     ->
 * 
 * Returns:
 *      void
 */

void next(double *z, double znext, double *err, char *error, char *where)
{
	double errc = 0.0; 

	if (fabs(*z - znext) < abs_prec)
        errc = fabs(*z - znext);
	else
        errc = fabs(*z - znext) / (fabs(znext) + abs_prec);

	if (errc > *err) {
        *err = errc;
        sprintf(error, "%s", where);
    }
	*z = znext; 
}


/*
 * angle_integrator_OP_fake(double complex)
 * 
 * TODO: Description
 * 
 * Args:
 *      double complex en   ->
 * 
 * Returns:
 *      void
 */

void angle_integrator_OP_fake(double complex en)
{

	for (int nx = -NNX; nx <= NNX; nx++) {
        for (int nz = -NZ; nz <= NZ; nz++) {
            sum11[nx][nz] = sum13[nx][nz] = sum22[nx][nz] = sum31[nx][nz]
            = sum33[nx][nz] = 0.0;
        }
    }
}


/*
 * angle_integrator(double complex, int)
 * 
 * TODO: Description
 * 
 * Args:
 *      double complex en   ->
 *      int AN              ->
 * 
 * Returns:
 *      void
 */

void angle_integrator(double complex en, int AN)
{
	int an, j, nx, nz;
	double vx, vy, vz, em;
	
	em = cimag(en);
	for (nx = -NNX; nx <= NNX; nx++) {
        for (nz = -NZ; nz <= NZ; nz++) {
            sum11[nx][nz] = sum13[nx][nz] = sum22[nx][nz] = sum31[nx][nz]
            = sum33[nx][nz] = 0.0;
        }
    }

	for (nx = -NNX; nx <= NNX; nx++)
        for (nz = -NNZ; nz <= NNZ; nz++)
            a[nx][nz][0] = aa[nx][nz][0] = 0.0;

    /* Angular integration */
    for (an = 1; an <= AN; an++) {
        get_bc(en, -NNX, an, a1_bc, aa1_bc);
        get_bc(en, NNX, an, a2_bc, aa2_bc);
        
            /* (initial) boundary condition on the perimeter */
            for (nz = -NNZ; nz <= NNZ; nz++) {
                /* calculate on Left and Right */
                for (j = 1; j <= 3; j++) { 
                    a[-NNX][nz][j] = a1_bc[nz][j];
                    aa[-NNX][nz][j] = aa1_bc[nz][j];
                    a[NNX][nz][j] = a2_bc[nz][j];
                    aa[NNX][nz][j] = aa2_bc[nz][j];
                }
            }

            for (nx = -NNX; nx <= NNX; nx++) {
                /* guess on Top and Bottom */
                for (j = 1; j <= 3; j++) { 
                    a[nx][NNZ][j] = a1_bc[NNZ][j];
                    aa[nx][NNZ][j] = aa2_bc[NNZ][j];
                    a[nx][-NNZ][j] = a1_bc[-NNZ][j];
                    aa[nx][-NNZ][j] = aa2_bc[-NNZ][j];
                }
            }

        /* Integration along trajectories */
        trajectory_integrator(en, an, A11, A22, A33, A13, A31);

        /* err = 0.0; */
        for (nx= -NX; nx<=NX; nx++) {
            for (nz= -NZ; nz<=NZ; nz++){
                f_from_ric(f, ff, a[nx][nz], aa[nx][nz]);
                
                /*
                 * if(nx==nxpr && nz==nzpr) fprintf(outf, "%f %e %e %e %e\n", 
                 * theta, creal(f[1]+ff[1]), creal(f[3]+ff[3]), px, pz);
                 */

                /*
                 * note that sum f and ff reflect em -> -em symmetry,
                 * p -> -p summetry is taken care of in dOmega
                 * creal - just to make the OP real
                 */

                sum11[nx][nz] += dOmg[an] * YY[1][an] * 
                    (0.5 * creal(f[1] + ff[1]) - A11[nx][nz] * YY[1][an] / em);
                /* sum22[nx][nz] = sum11[nx][nz]; */
                sum22[nx][nz] += dOmg[an] * YY[2][an] * 
                    (0.5 * creal(f[2] + ff[2]) - A22[nx][nz] * YY[2][an] / em);
                sum33[nx][nz] += dOmg[an] * YY[3][an] * 
                    (0.5 * creal(f[3] + ff[3]) - A33[nx][nz] * YY[3][an] / em);
                
                sum13[nx][nz] += dOmg[an] * YY[3][an] * 
                    (0.5 * creal(f[1] + ff[1]) - A13[nx][nz] * YY[3][an] / em);
                sum31[nx][nz] += dOmg[an] * YY[1][an] * 
                    (0.5 * creal(f[3] + ff[3]) - A31[nx][nz] * YY[1][an] / em);
            }
        }
    }

	printf("em = %f completed\n", em);
	fflush(stdout);
}


/*
 * trajectory_integrator(double complex, int, double **, double **, double **,
 *                      double **, double **)
 * 
 * We find Riccati amplitudes on the square grid, same grid as the OP grid. 
 * For square grid to find Riccati amplitude in one corner (x) we start
 * with extrapolated amplitude (e) from known 3 other corners (k), which
 * have been found in previous steps:
 *
 *	 k2------x
 *	  |     /|
 *	  |    / |
 *	  |   /  |
 *	  |  /   |                               Grid points notation:
 *	 k1-e----k3  (e) from {k1,k3}          [nx,nzz]---------[nxx,nzz]
 *	                                          |                 |
 *	 OR                                       |                 | 
 *	                                          |                 |
 *	 k2------x                                |                 |
 *	  |     /|                                |                 |
 *	  |    / |                                |                 |
 *	  |   /  |                             [nx,nz]----------[nxx,nz]
 *	  |  /   |
 *	  | /    |
 *	  |/     |
 *	  e      |
 *	  |      |
 *	 k1------k3  (e) from {k1,k2}
 *
 *	This is FORWARD integration. 
 *	Advance: Given nz, go over nx: nxi -> nxf; then increment nz; 
 *	Similarly for REVERSE integration
 *
 *	The convergence for Riccati is checked on the boundary:
 *  must be periodic (-NNZ = NNZ.)
 * 
 * Args:
 *      double complex en   ->
 *      int an              ->
 *      double **DEL11      ->
 *      double **DEL22      ->
 *      double **DEL33      ->
 *      double **DEL13      ->
 *      double **DEL31      ->
 * 
 * Returns:
 *      void
 */

void trajectory_integrator(double complex en, int an, double **DEL11,
        double **DEL22, double **DEL33, double **DEL13, double **DEL31)
{
	int j, nx, nxx, nz, nzz, nxi, nxf, nzi, nzf, dirx, dirz, converged; 
	int trip_across;
	double dx, dz, ds, ratio, ds_x, ds_z, d1, d2, err, rerr; 
	double complex ai[4], aai[4], Delta[4], dum[4];
	char name[LL], where[LL];

	for (j = 0;j <= 3;j++)
        dum[j] = Delta[0] = 0.0;

	if (Vx[an] >= 0) {
		dirx = 1; 
		nxi = -NNX;
		nxf = NNX;
	} else {
		dirx = -1; 
		nxi = NNX;
		nxf = -NNX;
	}

	if (Vz[an] >= 0) {
		dirz = 1; 
		nzi = -NNZ;
		nzf = NNZ;
	} else {
		dirz = -1; 
		nzi = NNZ;
		nzf = -NNZ;
	}

	/* 
     * printf("traject. integration for an=%d, V=(%f,%f,%f)\n",
     *      an, Vx[an], Vy[an], Vz[an]); 
     */

	trip_across = converged = 0;

	do {
		/* Forward integration */
		nz = nzi;
		do {
			nzz = nz + dirz;
			dz = zgr[nzz] - zgr[nz];
			nx = nxi;
			do {
				nxx = nx + dirx;
				dx = xgr[nxx] - xgr[nx];
				/* 
                 * printf("nz=%d(%d) dz=%f  |  nx=%d(%d), dx=%f\n",
                 *      nz, dirz, dz, nx, dirx, dx); 
                 */
				if (Vx[an] == 0 && Vz[an] == 0) {
					for( j = 0; j <= 3; j++)
                        ai[j] = 0.0;
					ds = 1e+10;
				} else if ((ratio = dx / Vx[an] * Vz[an] / dz) < 1) {
					for(j=0; j<=3; j++)
                        ai[j] = a[nx][nzz][j] * (1 - ratio) +
                            a[nx][nz][j] * ratio;
					ds = dx / Vx[an]; 
				} else {
                    ratio = 1.0 / ratio;
					for(j=0; j<=3; j++)
                        ai[j] = a[nxx][nz][j] * (1 - ratio) +
                            a[nx][nz][j] * ratio;
					ds = dz / Vz[an]; 
				}

				Delta[1] = DEL11[nxx][nzz] * YY[1][an] +
                    DEL13[nxx][nzz] * YY[3][an];
				Delta[2] = DEL22[nxx][nzz] * YY[2][an];
				Delta[3] = DEL31[nxx][nzz] * YY[1][an] +
                    DEL33[nxx][nzz] * YY[3][an];

				step_integratorc_a(ds, en, dum, Delta, ai, a[nxx][nzz]);
				nx = nxx;
			} while ((nxf - nx) * dirx > 0);
			nz = nzz;
		} while ((nzf - nz) * dirz > 0);

		/* Backward integration */
		nz = nzf;
		do {
			nzz = nz - dirz;
			dz = zgr[nz] - zgr[nzz];
			nx = nxf;
			do {
				nxx = nx - dirx;
				dx = xgr[nx] - xgr[nxx];
				if (Vx[an] == 0 && Vz[an] == 0) {
					for (j = 0; j <= 3; j++)
                        aai[j] = 0.0;
					ds = 1e+10;
				} else if ((ratio = dx / Vx[an]*Vz[an] / dz) < 1) {
					for(j = 0; j <= 3; j++)
                        aai[j] = aa[nx][nzz][j] * (1 - ratio) +
                            aa[nx][nz][j] * ratio;
					ds = dx / Vx[an]; 
				} else {
                    ratio = 1.0 / ratio;
					for (j = 0; j <= 3; j++)
                        aai[j] = aa[nxx][nz][j] * (1 - ratio) +
                            aa[nx][nz][j] * ratio;
					ds = dz / Vz[an]; 
				}

				Delta[1] = DEL11[nxx][nzz] * YY[1][an] +
                    DEL13[nxx][nzz] * YY[3][an];
				Delta[2] = DEL22[nxx][nzz] * YY[2][an];
				Delta[3] = DEL31[nxx][nzz] * YY[1][an] +
                    DEL33[nxx][nzz] * YY[3][an];

				step_integratorc_b(ds, en, dum, Delta, aai, aa[nxx][nzz]);
				nx = nxx;
			} while ((nx - nxi) * dirx > 0);
			nz = nzz;
		} while ((nz - nzi) * dirz > 0);

		trip_across++;
		err = 0;
        sprintf(where, "%s", "??????");

		if (!converged) {
            /* for future use with BULK */
			for (nx = -NNX; nx <= NNX; nx++) {
                for(j = 0; j <= 3; j++) {
                    sprintf(name, "a_mid[nx = %d]|", nx); 
                    nextc(&(a[nx][nzi][j]), a[nx][nzf][j], &err, where, name); 
                    sprintf(name, "aa_mid[nx = %d]|", nx); 
                    nextc(&(aa[nx][nzf][j]), aa[nx][nzi][j], &err, where, name); 
			    }
            }
		}

		if (err<precR || trip_across>100)
            converged = 1;

		if (!converged) {
			for (nx = -NNX; nx <= NNX; nx++)
                for(j = 0;j <= 3;j++)
                    a[nx][nzi][j] = a[nx][nzf][j];

			for (nx = -NNX; nx <= NNX; nx++)
                for (j = 0; j <= 3; j++)
                    aa[nx][nzf][j] = aa[nx][nzi][j];
		}

	} while (!converged);
}


/*
 * nextc(double complex *, double complex, double *, char *, char *)
 * 
 * TODO: Description
 * 
 * Args:
 *      double complex *z       ->
 *      double complex znext    ->
 *      double *err             ->
 *      char *error             ->
 *      char *where             ->
 * 
 * Returns:
 *      void
 */

void nextc(double complex *z, double complex znext, double *err,
        char *error, char *where)
{
	double errc = 0.0; 

	if (cabs(*z - znext) < abs_prec)
        errc = cabs(*z - znext);
	else
        errc = cabs(*z - znext) / (cabs(znext) + abs_prec);
	if (errc > *err) {
        *err = errc;
        sprintf(error, "%s", where);
    }
	*z = znext;
}


/*
 * get_bc(double complex, int, int, double complex **, double complex **)
 * 
 * Boundary conditions at xgr=\pm infinity
 * 
 * Args:
*       double complex en        ->
*       int nxb                  ->
*       int an                   ->
*       double complex **a_bc    ->
*       double complex **aa_bc   ->
* 
* Returns:
*       a_bc(xgr =-infty, theta, phi = 0)
*       aa_bc(xgr = +infty, theta, phi = 0)
*/

void get_bc(double complex en, int nxb, int an,
        double complex **a_bc, double complex **aa_bc) 
{
	int dirz, nz, nzz, nzi, nzf, done, converged, nzconverged, j;
	double err, dz, ds;
	double complex amplnext[4], Delta[4];
    double dum[4] = {0, 0, 0, 0};
	char name[LL], where[LL];

	/* 
     * for(j=0;j<=3; j++)
     *      dum[j]=0;
     */

	if (Vz[an] >= 0) {
		dirz = 1; 
		nzi = -NNZ;
		nzf = NNZ;
	} else {
		dirz = -1; 
		nzi = NNZ;
		nzf = -NNZ;
	}

	/* Forward */
	done = converged = nzconverged = 0;
	nz = nzi; 
	Delta[0] = 0.0;
	Delta[1] = A11[nxb][nz] * YY[1][an] + A13[nxb][nz] * YY[3][an];
	Delta[2] = A22[nxb][nz] * YY[2][an];
	Delta[3] = A31[nxb][nz] * YY[1][an] + A33[nxb][nz] * YY[3][an];
	ab_bulk_non_unitary(a_bc[nz], dum, en, Delta);
	/*
     * printf("get_bc at nxb=%d YY=(%f,%f,%f) nz=\n",
     *      nxb, YY[1][an], YY[2][an], YY[3][an]);
     */

	do {
		/* printf("%d  ", nz); */
		nzz = nz + dirz;
		dz = zgr[nzz] - zgr[nz];
		ds = dz / Vz[an];

		Delta[1] = A11[nxb][nzz] * YY[1][an] + A13[nxb][nzz] * YY[3][an];
		Delta[2] = A22[nxb][nzz] * YY[2][an];
		Delta[3] = A31[nxb][nzz] * YY[1][an] + A33[nxb][nzz] * YY[3][an];

		step_integratorc_a(ds, en, dum, Delta, a_bc[nz], a_bc[nzz]);

		/* if( abs(nzz)<=NZ ){} */
		if (nzz == nzf) {
			if (converged)
                done = 1;
			err = 0.0;
            sprintf(where, "%s", "?????????");
			for (j = 1; j <= 3; j++) {
				sprintf(name, "a_bc[%d]|", j); 
				nextc(&(a_bc[nzi][j]), a_bc[nzf][j], &err, where, name); 
			}
			if (err < precR)
                converged = 1;
			nz = nzi;
		/* printf("\n"); */
		} else {
            nz = nzz;
        } 
	} while (!done);

	/* printf("\n"); */

	/* Backward */
	done = converged = 0;
	nz = nzf; 
	Delta[0] = 0.0;
	Delta[1] = A11[nxb][nz] * YY[1][an] + A13[nxb][nz] * YY[3][an];
	Delta[2] = A22[nxb][nz] * YY[2][an];
	Delta[3] = A31[nxb][nz] * YY[1][an] + A33[nxb][nz] * YY[3][an];
	ab_bulk_non_unitary(dum, aa_bc[nz], en, Delta);

	do {
		nzz = nz-dirz;
		dz = zgr[nz] - zgr[nzz];
		ds = dz / Vz[an];

		Delta[1] = A11[nxb][nzz] * YY[1][an] + A13[nxb][nzz] * YY[3][an];
		Delta[2] = A22[nxb][nzz] * YY[2][an];
		Delta[3] = A31[nxb][nzz] * YY[1][an] + A33[nxb][nzz] * YY[3][an];

		step_integratorc_b(ds, en, dum, Delta, aa_bc[nz], aa_bc[nzz]);

		if (nzz == nzi) {
			if (converged)
                done = 1;
			err = 0.0;
            sprintf(where, "%s", "?????????");
			for (j = 1; j <= 3; j++) {
				sprintf(name, "aa_bc[%d]|", j); 
				nextc(&(aa_bc[nzf][j]), aa_bc[nzi][j], &err, where, name); 
			}
			if (err < precR)
                converged = 1;
			nz = nzf;
		} else {
            nz = nzz;
        } 
	} while (!done);
		
	/* printf("%s\n", "get_bc DONE"); */
}


/*
 * OP_carpet_DW(void)
 * 
 * TODO: Description
 * 
 * Args:
 *      void
 * 
 * Returns:
 *      void
 */

void OP_carpet_DW(void)
{
	int sgn, nz, nx, nz1;

	for (nz = -NNZ; nz <= NNZ; nz++) {
		if (nz < -NZ) {
            nz1 = -(nz + 2 * NZ);
            sgn = -1;
        } else {
			if (nz > NZ) {
                nz1 = 2 * NZ - nz;
                sgn = -1;
            } else {
                nz1=nz;
                sgn=1;
            }
		}

		for (nx = -NNX; nx <= NNX; nx++) {
			A11[nx][nz] = A11[nx][nz1];
			A22[nx][nz] = A22[nx][nz1];
			A33[nx][nz] = sgn*A33[nx][nz1];
			A13[nx][nz] = sgn*A13[nx][nz1];
			A31[nx][nz] = A31[nx][nz1];
		}
	}	
}


/*
 * master(int, int)
 * 
 * TODO: Description
 * 
 * Args:
 *      int AN          ->
 *      int self_cnst   ->
 * 
 * Returns:
 *      void
 */

static void master(int AN, int self_cnst)
{
	time_t time_of_start;

	int done, nx, nz, m, mats_cutoff, iter, k, j, nxerr, nzerr;
    int component, corner;
	double em, *emgrid, ln, err, err1, convpoints;
	double complex en;
	double **Sum11, **Sum13, **Sum22, **Sum31, **Sum33;
	char string[LL], name[LL], where[LL];

	FILE *outgap, *initguess, *gapcurrent;

    int v;
    int Svector, Nvector;
    Svector = Nvector = 5 * (2 * NX + 1) * (2 * NZ + 1);
    double complex *vec_next, *Fvec_next, *Fmin_prev;
    double *Normf; 
    double complex **vector, **Fvec;

	/* richardson extrapolation variables for Matsubara sums */
	int enf, i_rich, period_rich, count_rich, shift, extrapol;
    int rich_converg, en_rich[N_rich + 1];
	double *grid_rich, x_rich[N_rich + 1];
	double complex sum_element_rich;
	double complex **S_rich, *S_rich_prev, *S_rich_new, *S_rich_direct_sum;

	int tag, m_slave, n_tasks, n_proc, who; 
	MPI_Status status; 
	MPI_Datatype MPIvect, MPImatrix;

	MPI_Type_vector(2 * NNX + 1, 1, 1, MPI_DOUBLE, &MPIvect);
	MPI_Type_commit(&MPIvect);
	MPI_Type_vector(2 * NNX + 1, 2 * NNZ + 1, 2 * NNZ + 1,
                    MPI_DOUBLE, &MPImatrix);
	MPI_Type_commit(&MPImatrix);

	/* outtest=fopen("Sum.dat","w"); */
	time_of_start = starttime();
   
	Sum11 = dmatrix(-NNX, NNX, -NZ, NZ);
	Sum13 = dmatrix(-NNX, NNX, -NZ, NZ);
	Sum22 = dmatrix(-NNX, NNX, -NZ, NZ);
	Sum31 = dmatrix(-NNX, NNX, -NZ, NZ);
	Sum33 = dmatrix(-NNX, NNX, -NZ, NZ);

    vec_next = Cvector(1, Nvector); 
    Fvec_next = Cvector(1, Nvector); 
    Fmin_prev = Cvector(1, Nvector); 
    Normf = dvector(1, Mk);
    vector = Cmatrix(1, Mk, 1, Nvector); 
    Fvec = Cmatrix(1, Mk, 1, Nvector);

	ln = log(t);
	mats_cutoff = Mats(t);                  /* cutoff for Matsubara */
	emgrid = dvector(0, mats_cutoff);
	for (m = 0; m <= mats_cutoff; m++)
        emgrid[m] = t * (m + 0.5);

	S_rich = Cmatrix(1, Svector, 0, N_rich);
	S_rich_prev = Cvector(1, Svector);
	S_rich_new = Cvector(1, Svector);
	S_rich_direct_sum = Cvector(1, Svector);
	grid_rich = emgrid;
	period_rich = 6;                        /* Sampling period in Matsubaras */

	if (period_rich * (N_rich + 1) > mats_cutoff)
        period_rich = 1.0 * mats_cutoff / (N_rich + 1);


	printf("%s",
        "\nOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
	printf("   OP calculation: t=%f  D=%f   NNZ=%d\n", t, D, NNZ);
	printf("%s",
        "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");

	/* how many nodes? */
	 MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
	 printf("there are %d slaves\n", n_tasks-1);

	for ( k = 1; k <= Mk; k++)
        Normf[k] = 0.0;

	done = iter = 0;

    /* iterations for the OP begin */
	do {
		iter++;
		for (nx = - NNX; nx <= NNX; nx++) {
		    for (nz= -NZ; nz<=NZ; nz++) {
			    Sum11[nx][nz] = Sum13[nx][nz] = Sum22[nx][nz] =
				Sum31[nx][nz] = Sum33[nx][nz] = 0.0;
            }
        }

	printf("\niter=%d\n", iter);
	OP_carpet_DW();                   /* make new carpet from the OP */

	if (iter == 1){
	  sprintf(string, "initguess_carpet.dat");
	  initguess = fopen(string, "w");
	  /* fprintf(initguess, "# dx = %f dz = %f t=%f\n", dx, dz, t); */
	  for (nx = - NNX; nx <= NNX; nx++) {
	    for (nz = - NNZ; nz <= NNZ; nz++) {
	    	fprintf(initguess, "%f\t%f \t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
	  	        xgr[nx], zgr[nz], A11[nx][nz], A22[nx][nz], A33[nx][nz],
                A13[nx][nz], A31[nx][nz]);
	    }
      }
	  fflush(initguess);
	  fclose(initguess);
	}

	sprintf(string, "gap_curcarp.dat");
	gapcurrent = fopen(string, "w");
	/* fprintf(gapcurrent, "# dx = %f dz = %f t=%f\n", dx, dz, t); */
	for (nx = - NNX; nx <= NNX; nx++) {
        for (nz = - NNZ; nz <= NNZ; nz++){
            fprintf(gapcurrent, "%f\t%f \t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
                xgr[nx], zgr[nz], A11[nx][nz], A22[nx][nz], A33[nx][nz],
                A13[nx][nz], A31[nx][nz]);
        }
    }
	fflush(gapcurrent);
	fclose(gapcurrent);

    MPI_Bcast(&(A11[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(A13[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(A22[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(A31[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(A33[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
	
	/* printf("\nBroadcasted new OP for \titer=%d\n", iter); */

    /* go over Matsubara energies */
	m = 0; 
	for (n_proc = 1; n_proc < n_tasks; n_proc++) { 
		em = emgrid[m]; 
		en = I * em;
		MPI_Send(&m, 1, MPI_INT, n_proc, tag_work, MPI_COMM_WORLD); 
		MPI_Send(&en, 1, MPI_DOUBLE_COMPLEX, n_proc, tag_work, MPI_COMM_WORLD); 
		m++; 
	}

	while (m <= mats_cutoff) {
		em = emgrid[m];
		en = I * em;
		
		/* calculates sums[nx] for given em */
		MPI_Recv(&m_slave, 1, MPI_INT, MPI_ANY_SOURCE, 30, 
				MPI_COMM_WORLD, &status); 
		who = status.MPI_SOURCE;
		MPI_Recv(&(sum11[-NNX][-NNZ]), 1, MPImatrix, who, 11,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum13[-NNX][-NNZ]), 1, MPImatrix, who, 13,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum22[-NNX][-NNZ]), 1, MPImatrix, who, 22,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum31[-NNX][-NNZ]), 1, MPImatrix, who, 31,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum33[-NNX][-NNZ]), 1, MPImatrix, who, 33,
				MPI_COMM_WORLD, &status); 

		MPI_Send(&m, 1, MPI_INT, who, tag_work, MPI_COMM_WORLD); 
		MPI_Send(&en, 1, MPI_DOUBLE_COMPLEX, who, tag_work, MPI_COMM_WORLD);

		for (nx = -NNX; nx <= NNX; nx++) {
            for (nz = -NZ; nz <= NZ; nz++){
                Sum11[nx][nz] += t * sum11[nx][nz];
                Sum13[nx][nz] += t * sum13[nx][nz];
                Sum22[nx][nz] += t * sum22[nx][nz];
                Sum31[nx][nz] += t * sum31[nx][nz];
                Sum33[nx][nz] += t * sum33[nx][nz];
            }
            m++;
        }
	}
	/* receive the remainings */
	for (n_proc = 1; n_proc < n_tasks; n_proc++) { 
		MPI_Recv(&m_slave, 1, MPI_INT, MPI_ANY_SOURCE, 30, 
				MPI_COMM_WORLD, &status); 
		who = status.MPI_SOURCE;
		MPI_Recv(&(sum11[-NNX][-NNZ]), 1, MPImatrix, who, 11,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum13[-NNX][-NNZ]), 1, MPImatrix, who, 13,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum22[-NNX][-NNZ]), 1, MPImatrix, who, 22,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum31[-NNX][-NNZ]), 1, MPImatrix, who, 31,
				MPI_COMM_WORLD, &status); 
		MPI_Recv(&(sum33[-NNX][-NNZ]), 1, MPImatrix, who, 33,
				MPI_COMM_WORLD, &status); 

		for (nx = -NNX; nx <= NNX; nx++) {
            for (nz = -NZ; nz <= NZ; nz++){
                Sum11[nx][nz] += t * sum11[nx][nz];
                Sum13[nx][nz] += t * sum13[nx][nz];
                Sum22[nx][nz] += t * sum22[nx][nz];
                Sum31[nx][nz] += t * sum31[nx][nz];
                Sum33[nx][nz] += t * sum33[nx][nz];
            }
        }
	}

	v = 0;
	for (nx = -NX; nx <= NX; nx++) {
        for (nz = -NZ; nz <= NZ; nz++) {
            vec_next[++v] = A11[nx][nz];
            Fvec_next[v] = fsY2[1][1] * A11[nx][nz]
                        - (1.0 / ln * Sum11[nx][nz]);
            
            vec_next[++v] = A13[nx][nz];
            Fvec_next[v] = fsY2[3][3] * A13[nx][nz]
                        - (1.0 / ln * Sum13[nx][nz]);
            
            vec_next[++v] = A22[nx][nz];
            Fvec_next[v] = fsY2[2][2] * A22[nx][nz]
                        - (1.0 / ln * Sum22[nx][nz]);
            
            vec_next[++v] = A31[nx][nz];
            Fvec_next[v] = fsY2[1][1] * A31[nx][nz]
                        - (1.0 / ln * Sum31[nx][nz]);
            
            vec_next[++v] = A33[nx][nz];
            Fvec_next[v] = fsY2[3][3] * A33[nx][nz]
                        - (1.0 / ln * Sum33[nx][nz]);
        }
    }

	/* 
     * After 40 iterations do simple relaxation again,
     * for this just make Normf[k] = 0
     */

	if (iter < 10)
        for (j = 4; j <= Mk; j++)
            Normf[j] = 0.0; 


	OP_update(Nvector, vec_next, Fvec_next, vector, Fvec, Fmin_prev, Normf);


	sprintf(where, "????");
	err = 0.0;
	v = 0;

	for (nx = -NX; nx <= NX; nx++) {
        for (nz = -NZ; nz <= NZ; nz++)
        { 
            sprintf(name, "A11[%d][%d]", nx, nz);
            next(&(A11[nx][nz]), creal(vec_next[++v]), &err, where, name); 

            sprintf(name, "A13[%d][%d]", nx, nz);
            next(&(A13[nx][nz]), creal(vec_next[++v]), &err, where, name); 

            sprintf(name, "A22[%d][%d]", nx, nz);
            next(&(A22[nx][nz]), creal(vec_next[++v]), &err, where, name); 

            sprintf(name, "A31[%d][%d]", nx, nz);
            next(&(A31[nx][nz]), creal(vec_next[++v]), &err, where, name); 

            sprintf(name, "A33[%d][%d]", nx, nz);
            next(&(A33[nx][nz]), creal(vec_next[++v]), &err, where, name); 
        }
    }

	printf("After iter %d,  MAX error is %e in %s \n", iter, err, where); 
	
	sprintf(string, "gap_current.dat");
	gapcurrent = fopen(string, "w");
	fprintf(gapcurrent, "\n"); //# dx = %f dz = %f t=%f\n", dx, dz, t);
	for (nx = -NNX; nx <= NNX; nx++) {
        for (nz = -NNZ; nz <= NNZ; nz++){
            fprintf(gapcurrent, "%f\t%f \t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
            xgr[nx], zgr[nz], A11[nx][nz], A22[nx][nz], A33[nx][nz],
            A13[nx][nz], A31[nx][nz]);
        }
    }
	fflush(gapcurrent);
	fclose(gapcurrent);

	if (err < prec || iter > 50)
        done = 1;
	else
        done = 0;

	if (!done)
        tag = tag_relax;
	else
        tag = tag_stop;

	for (n_proc = 1; n_proc < n_tasks; n_proc++) { 
		MPI_Send(0, 0, MPI_INT, n_proc, tag, MPI_COMM_WORLD); 
		MPI_Send(0, 0, MPI_INT, n_proc, tag, MPI_COMM_WORLD);
	}

	} while(!done); 

	sprintf(string, "gapFilm_t%.3f_D%.3f_DW.dat", t, D);
	outgap = fopen(string, "w");
	fprintf(outgap, "%s",
        "#  OP components are in units of Tc \n");
	fprintf(outgap, "%s",
        "#  x\t\t\tz \t\tAxx \t\t\tAyy \t\t\tAzz \t\t\tAxz \t\t\tAzx \n");
	
	for (nx = -NNX; nx <= NNX; nx++) {
        for (nz = -NZ; nz <= NZ; nz++){
            fprintf(outgap, "%s",
                "%f\t%f \t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
                xgr[nx], zgr[nz], A11[nx][nz] * PI2, A22[nx][nz] * PI2,
                A33[nx][nz] * PI2, A13[nx][nz] * PI2, A31[nx][nz] * PI2);
	    }
    }
	fprintf(outgap, "%s", "\n");
	fflush(outgap);

	free_dmatrix(Sum11, -NNX, NNX, -NZ, NZ);
	free_dmatrix(Sum13, -NNX, NNX, -NZ, NZ);
	free_dmatrix(Sum22, -NNX, NNX, -NZ, NZ);
	free_dmatrix(Sum31, -NNX, NNX, -NZ, NZ);
	free_dmatrix(Sum33, -NNX, NNX, -NZ, NZ);

    free_Cvector(vec_next, 1, Nvector);
    free_Cvector(Fvec_next, 1, Nvector);
    free_Cvector(Fmin_prev, 1, Nvector);
    free_dvector(Normf, 1, Mk);
    free_Cmatrix(vector, 1, Mk, 1, Nvector);
    free_Cmatrix(Fvec, 1, Mk, 1, Nvector);

	fclose(outgap);

	MPI_Type_free(&MPIvect);
	MPI_Type_free(&MPImatrix);

	endtime(time_of_start);
	/* Return int: MPI error code ? */
}


/*
 * slave(int)
 * 
 * TODO: Description
 * 
 * Args:
 *      int AN  ->
 * 
 * Returns:
 *      void
 */

static void slave(int AN)
{
	int tag, myrank, m_slave;
	double complex en;
	/* 
     * FILE *out;
     * out = fopen("Slave_OP.dat","w");
     */

	MPI_Status status; 
	MPI_Datatype MPIvect, MPImatrix;

	MPI_Type_vector(2 * NNX + 1, 1, 1,MPI_DOUBLE, &MPIvect);
	MPI_Type_commit(&MPIvect);
	MPI_Type_vector(2 * NNX + 1, 2 * NNZ + 1, 2 * NNZ + 1,
                    MPI_DOUBLE, &MPImatrix);
	MPI_Type_commit(&MPImatrix);

	while (1) {
		MPI_Bcast(&(A11[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(A13[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(A22[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(A31[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(A33[-NNX][-NNZ]), 1, MPImatrix, 0, MPI_COMM_WORLD);
		/*for (nz= -NNZ;nz<=NNZ;nz++)
			fprintf(out, "%d\t%e\t%e\t%e\n", nz, A33[-NNX][nz], A33[6][nz], A33[NNX][nz]);
		fflush(out);*/
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		/*printf("A33[-NNX][-NNZ]=%f\n", A33[-NNX][-NNZ]);
		printf("rank %d received Broadcast\n", myrank);
		fflush(stdout);*/
		while (1) {
			MPI_Recv(&m_slave, 1, MPI_INT, 0, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);
			MPI_Recv(&en, 1, MPI_DOUBLE_COMPLEX, 0, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
			if (tag == tag_relax)
                break;
			if (tag == tag_stop)        /* (?) Should be else? */
                return;
		//printf("rank %d received Send\n", myrank);

			angle_integrator(en, AN);
			
		    printf("rank \t%d  finished \tm=%d \n", myrank, m_slave);
            fflush(stdout);

			/* 
             * Sends tag - garbage, one could send Matsubara. 
			 * This is to tell Master to listen to this Slave only
             */
			MPI_Send(&m_slave, 1, MPI_INT, 0, 30, MPI_COMM_WORLD); 
			
			MPI_Send(&(sum11[-NNX][-NNZ]), 1, MPImatrix, 0, 11, MPI_COMM_WORLD); 
			MPI_Send(&(sum13[-NNX][-NNZ]), 1, MPImatrix, 0, 13, MPI_COMM_WORLD); 
			MPI_Send(&(sum22[-NNX][-NNZ]), 1, MPImatrix, 0, 22, MPI_COMM_WORLD); 
			MPI_Send(&(sum31[-NNX][-NNZ]), 1, MPImatrix, 0, 31, MPI_COMM_WORLD); 
			MPI_Send(&(sum33[-NNX][-NNZ]), 1, MPImatrix, 0, 33, MPI_COMM_WORLD); 
		}
	}

	MPI_Type_free(&MPIvect);
	MPI_Type_free(&MPImatrix);
	/* (?) Should return int: MPI error? */
}