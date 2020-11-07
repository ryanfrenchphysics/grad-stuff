#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "nrutil.h"
#include "ricparam.h"
#include "mpi_timer.h"
#include "subroutines.h"
#include "constvar.h"
#include "get_OP_FE.h"
#include "return_values.h"

/* (?) Where are these functions defined? */
void FSsetup_simple(int *AN, int ANzmax, double complex **g_dos); 
void FSsetup_circle(int *ANmax, double complex **g_dos); 

/*
	bcd,bcu = 1(specular), 0(ovchinnikov)
  z
  |
  bc                 
 > - - - - - - - - - D/2 - - - - - - - - - < bcu
 >                    & theta            <
 >                    &   /              <
 >                    &--/               <
 >                    & /                <
 -LX ---------------- 0------------------LX---x
 >                    &                  <
 >                    &                  <
 >                    &                  <
 >                    &                  <
 >                    &                  <
 >  - - - - - - - - -D/2 - - - - - - - - < bcd

FS is sphere 4Pi - Lebedev parametrization
*/

int main(int argc, char *argv[])
{
	int nproc, myrank, Angles, nx, nz, done, namelength, selfcnst;
    int freeenergy, rescaleOP, readOPfilm; 
	double LX, D, Di, Df, dD, t, Dprev, Dcrit, FEd, FEb, FEdp, FEbp, dF, dFp;
	double *xgr, *zgr; /*, dx, dz; */
	double **A11, **A13, **A22, **A31, **A33, **FrEn;
	double complex **g_dos; 
	char string[LL], data[LL];
	FILE *in, *out;
	MPI_Datatype MPIxvect, MPIzvect;

	MPI_Init(&argc, &argv);
	MPI_Type_vector(2 * NNX + 1, 1, 1, MPI_DOUBLE, &MPIxvect);
	MPI_Type_commit(&MPIxvect);
	MPI_Type_vector(2 * NNZ + 1, 1, 1, MPI_DOUBLE, &MPIzvect);
	MPI_Type_commit(&MPIzvect);

	xgr = dvector(-NNX, NNX);
	zgr = dvector(-NNZ, NNZ);
	A11 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	A13 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	A22 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	A31 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	A33 = dmatrix(-NNX, NNX, -NNZ, NNZ);
	FrEn = dmatrix(-NNX, NNX, -NNZ, NNZ);
	
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0) {
		printf("Given %d nodes.\n", nproc);
		if (nproc < 2) {
			printf("Need more than 1 node! \n");
			MPI_Finalize();
			exit(1);
		}
	}

	MPI_Get_processor_name(string, &namelength);
	printf("Connection set up with: \t(p%d) \t%s\n", myrank, string);
    fflush(stdout);

	if (myrank == 0) {
		if ((in = fopen("input.dat", "r")) == NULL) {
			printf("File <<input.dat>> does not exist!\n");
			MPI_Finalize();
			exit(1);
		} else {
			fgets(data, LL, in);
            sscanf(data, "%d  %d", &selfcnst, &freeenergy);
            printf("self-consistency=%d\n", selfcnst);

			fgets(data, LL, in);
            sscanf(data, "%lf", &t);
            printf("t =%f\n", t);

			fgets(data, LL, in);
            sscanf(data, "%lf", &LX);
            printf("LX =%f\n", LX);

			fgets(data, LL, in);
            sscanf(data, "%lf  %lf  %lf", &Di, &Df, &dD); 
            printf("Di=%f\t", Di);
            printf("Df=%f\t", Df);
            printf("dD=%f\n", dD);

			fgets(data, LL, in);
            sscanf(data, "%d", &Angles);
            printf("Approximete number of FS angles: Angles=%d\n", Angles);

			fgets(data, LL, in);
            sscanf(data, "%d", &readOPfilm);
            printf("readOP=%d\n", readOPfilm);
		}
	}

	MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&LX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Di, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Df, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Angles, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		sprintf(string, "fe_dens_t%.3f.dat", t);
		out = fopen(string, "a");
	}

	/*
     * FS grid is set up in get_OP_FE.c as static variables
     * for FS integration in that file
     */
	FSsetup(&Angles, g_dos);        /* Lebedev sphere */

	D = Di;
	Dprev = D - dD;
	Dcrit = D;
	FEdp = FEbp = 10000.0;
	rescaleOP = 0;

	do {
		done = 0;
		xgr[0] = 0.0;
		for (nx = -NNX; nx <= NNX; nx++) {
			xgr[nx] = LX * nx / NNX * pow(fabs(1.0 * nx / NNX), 0.5); 
		}

		zgr[0] = 0.0;
		for (nz = 1; nz <= NZ; nz++) {
			zgr[nz] = 0.5 * D * (1 - pow(fabs(nz - NZ) / 1.0 / NZ, 1.5)); 
			zgr[-nz] = -zgr[nz]; 
		}

		for(nz=NZ+1;nz<=NNZ;nz++) {
            zgr[nz] = zgr[NZ] + (zgr[NZ] - zgr[2 * NZ - nz]);
            zgr[-nz] = -zgr[nz];
        }

		if (myrank == 0)
			op_guess(D, t, xgr, zgr, A11, A22, A33, A13, A31,
                rescaleOP, readOPfilm);

		rescaleOP = 1;
        readOPfilm = 0;

		domain_OP(selfcnst, Angles, t, D, xgr, zgr,
            A11, A22, A33, A13, A31, argc, argv);
		printf("Exited OP - p%d\n", myrank);

		/*domainFE(&FEd, &FEb, Angles, t, D, xgr, zgr,
         *  A11, A22, A33, A13, A31, argc, argv); */

		/* printf("Exited FE - p%d\n", myrank); */

		if (myrank == 0) {
			fprintf(out, "t=%f D=%f FEd=%e FEdunif=%e\n", t, D, FEd, FEb);
			dF = FEd - FEb;
			dFp = FEdp - FEbp;
			if (dF < 0) {
				Dcrit = (dF * Dprev - dFp * D) / (dF - dFp);
				fprintf(out, "#SDW transition for t=%f at Dcrit=%f\n",
                    t, Dcrit);
				done = 1;
			}
			fflush(out);
		}

		FEdp = FEd;
        FEbp = FEb;
		Dprev = D;
		D += dD;
		if ((Df - D) * dD < -1e-5)
            done = 1;

		MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} while (!done);

	free_dvector(xgr, -NNX, NNX);
	free_dvector(zgr, -NNZ, NNZ);
	
	free_dmatrix(A11, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(A13, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(A22, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(A31, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(A33, -NNX, NNX, -NNZ, NNZ);
	free_dmatrix(FrEn, -NNX, NNX, -NNZ, NNZ);

	MPI_Type_free(&MPIxvect);
	MPI_Type_free(&MPIzvect);
	MPI_Finalize();
	return NO_ERROR;
}