#ifndef _INTEGRATORS_H
#define _INTEGRATORS_H

void integratora(double ds, int nsinit, int nsfinal , double tm, 
                double complex **Deltas, double complex **c);
                                                                                  
void integratorb(double ds, int nsinit, int nsfinal, double tm,
		double complex **Deltas, double complex **c);

void integratorac(double ds, int nsinit, int nsfinal , double complex tm, 
                double complex **Deltas, double complex **c);
                                                                                  
void integratorbc(double ds, int nsinit, int nsfinal, double complex tm,
		double complex **Deltas, double complex **c);

void gapalongtraj(int nsinit, int nsfinal, int nz,
	double theta, double complex **Delta, double complex **Deltas, char phase);

void step_integrator_a(double ds, double tm, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

void step_integrator_b(double ds, double tm, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

void step_integratorc_a(double ds, double complex tm, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

void step_integratorc_b(double ds, double complex tm, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

void step_integratorm_a(double ds, double complex en, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

void step_integratorm_b(double ds, double complex en, double complex *Deltai,
                double complex *Deltaf, double complex *ci, double complex *cf);

#endif /* _INTEGRATORS_H */

