#ifndef _RICPARAM_H
#define _RICPARAM_H

#include <complex.h>

/* UNITARY. 'ra, rb, fv, ffv, Delta' all have 3 components: [1,2,3] */
void ab_bulk(double complex *ra, double complex *rb, 
		double complex en, double complex *Delta);
void ab_bulk_non_unitary(double complex *ra, double complex *rb, 
		double complex en, double complex *Delta);
void a_from_g(double complex *ra, double complex *g, 
		double complex *f, double complex *ff);
void b_from_g(double complex *rb, double complex *gg, 
		double complex *f, double complex *ff);
void g_from_ric(double complex *g, double complex *ra, double complex *rb);
void gg_from_ric(double complex *gg, double complex *ra, double complex *rb);
void f_from_ric(double complex *f, double complex *ff,
		double complex *ra, double complex *rb);

// extras: matrix operations + some Ovchinnikov BC subroutines
void ric_multiply(double complex *x, double complex *y, double complex *result); 
void ric_invert(double complex *x, double complex *x_inv);
void ric_sqrt(double complex *x, double complex *x_sqrt);
void G_multiply(double complex *Ag, double complex *Agg, double complex *Af, double complex *Aff, 
		    double complex *Bg, double complex *Bgg, double complex *Bf, double complex *Bff, 
		    double complex *Cg, double complex *Cgg, double complex *Cf, double complex *Cff); 
void G_invert(double complex *Ag, double complex *Agg, double complex *Af, double complex *Aff, 
		    double complex *Bg, double complex *Bgg, double complex *Bf, double complex *Bff);

// NON-UNITARY g_bulk from spin-vector Delta[0,1,2,3], 'g,gg,f,ff' [0,1,2,3] 
void g_bulk(double complex *g, double complex *gg, double complex *f, double complex *ff, 
		double complex en, double complex *Delta);
void gS_bulk(double complex *g, double complex *f, double complex *ff, 
		double complex en, double B, double complex Delta);
void afromg(double complex *ra, double complex gs, 
		double complex *fv, double complex *ffv);
void bfromg(double complex *rb, double complex gs, 
		double complex *fv, double complex *ffv);
double complex gsfromric(double complex *ra, double complex *rb);
void ffromric(double complex *fv, double complex *ffv,
		double complex *ra, double complex *rb);
void gfvfromric(double complex *gv, double complex *fv, double complex *ffv,
		double complex *ra, double complex *rb);


/* (?) What are these functions for? Where are they declared? */
double complex ovchinBg(double complex g, double complex g0, double complex *fv,
                double complex *ffv, double complex *fv0, double complex *ffv0);
void ovchinBgv(double complex *Bgv, double complex *gv, double complex *gv0,
                double complex *fv, double complex *ffv,
                double complex *fv0, double complex *ffv0);
void ovchinBfv(double complex *Bfv, double complex g, double complex g0,
                double complex *gv, double complex *gv0,
                double complex *fv, double complex *fv0);
void ovchinBffv(double complex *Bfv, double complex g, double complex g0,
                double complex *gv, double complex *gv0,
                double complex *ffv, double complex *ffv0);

#endif /* _RICPARAM_H */
