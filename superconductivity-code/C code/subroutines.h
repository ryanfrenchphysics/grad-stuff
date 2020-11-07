#include <complex.h>

void OP_update(int Nvector, double complex *vec_next, 
	double complex *Fvec_next, double complex *vector[], 
	double complex *Fvec[], double complex *Fmin_prev, double Norm[]);

double bulkgap(double t, char phase); 

void step_integratorc_a(double ds, double complex tm, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 
void step_integratorc_b(double ds, double complex tm, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 
void step_integratorm_a(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 
void step_integratorm_b(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 
void step_integratornu_a(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 
void step_integratornu_b(double ds, double complex en, double complex *Deltai, 
		double complex *Deltaf, double complex *ci, double complex *cf); 

//---------------- problem specific subroutines --------------------------------------------------------
void op_guess(double dfilm, double t, double *xgr, double *zgr, 
	double **A11, double **A22, double **A33, double **A13, double **A31, int rescaleOP, int readOPfilm);

void domainOP(int SelfCons, int AN, double t, double D, double *xgr, double *zgr, 
	double **A11, double **A22, double **A33, double **A13, double **A31, int argc, char *argv[]);

void domainFE(double *FEd, double *FEb, double t, double D, double *xgr, double *zgr, 
	double **A11, double **A22, double **A33, double **A13, double **A31, int argc, char *argv[]);

void domainDOS( double Dthick, double *xgrid, double *zgrid, int *nx_dos, int *nz_dos, int *pan_dos, int *an_dos,
	double **OP11, double **OP22, double **OP33, double **OP13, double **OP31, int argc, char *argv[]);

//Richardson extrapolation
void en_placement_rich(int i_new, int *i_rich, int next_i, int *count_rich, int *shift,  int *extrapol);
void shiftX(double *x_rich, double next_x);
void sumupdate(double complex *Sum, double complex sum_new_element, int i_new, int *i_rich, int shift);
double complex sum_extrapol(int N, double complex *sum, double *x);
void se_extrapol(int ni, int period, int nf, /*double complex *se,*/ double complex **SE, int Nch, double *x);

// ------subroutins from the slab code ----------- // 
void get_op(int *SelfCons, int AN, double bc, double bcf, double t, double As1, double *Ps, double *xgr, double complex **g_dos, 
		double complex ***Amatrix, 
		double complex **a_restr, double *Nu_prl, double *Nu_prp, double *jparal, double *jperp, double **jspin, char *string); 
void	get_dos(int ENi, int ENf, int AN, double imp_gamma, double imp_shift, double *Ps, double *xgr, double complex **g_dos, 
			double complex ***Amatrix, double *Nu_prl, double *Nu_prp, char *string); 
void	get_free_energy(double t, double A1s, int AN, double imp_gamma, double imp_shift, double *Ps, double *xgr, 
			double complex **g_dos, double *free_en, 
			double complex ***Amatrix, double *Nu_prl, double *Nu_prp, char *string);
void dFErecursive(double *dfeAB, double emA, double emB, double *feA, double *feB, int AN, double bc, double bcf, double *xgr, 
			double complex ***Amatrix, double *Ps, double *Nu_prl, double *Nu_prp, double *gap2, 
			double complex ***dumfsF, double complex *fsG, double complex **dumfsG, double complex **g_dos, 
			double prec_lambda_int, int *depth); 

