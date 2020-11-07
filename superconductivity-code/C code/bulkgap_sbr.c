/*
 *
 * 	This function calculates bulk value of the gap 
 * 	for temperature z_tt and phase 'A' or 'P' or 'B' or 'D'
 *
*/

#include <stdio.h>
#include <math.h>

#ifndef M_PI
    /* In case it isn't imported from math.h */
    #define M_PI 3.14159265358979323846
#endif

#define DPREC 1.0e-6



/*
 * bulkgap (double, double)
 * 
 * 		Calculates the bulk value of the gap for temperature
 * 
 * Args:
 * 		double z_tt 	->	temperature (?)
 * 		char phase	->	(?)
 * 
 * Returns:
 * 		double	->	Delta (?)
 * 
 */
double bulkgap(double z_tt, char phase)
{
	double z_ID[2][10], z_EI3, z_EI4, ei4;
	int z_M, z_k, pan_aux, pan, z_PAN;
	double z_em, z_srt, z_D, z_Dder, z_Delta;
	double z_Deltanext, err, theta, dtheta, z_Y;

	double ec3[] = {
		0.0,
		0.0,
		1.0/640.0,
		-90.0/640.0,
		729.0/640.0
	};
	double ec4[] = {
		0.0,
		-1.0/465920.0,
		9.0/5120.0,
		-729.0/5120.0,
		531441.0/465920.0
	};
	
	/* (?) Purpose? */
	if (z_tt == 0.0)
		z_tt = 1.0e-3;


	z_M = 25.0 / z_tt;
	z_Delta = 0.2;

	if (phase =='B') {
		/* (?) What does this mean? */
	  do {
		z_D = z_Delta * log(z_tt);
		z_Dder = log(z_tt);
		for (int m = 0; m < z_M; m++) {
			z_em = z_tt * (m + 0.5);
			z_srt = sqrt((z_em * z_em) + (z_Delta * z_Delta));
			
			z_D -= z_tt * z_Delta * (1.0 / z_srt - 1.0 / z_em);
			z_Dder -= z_tt * ((z_em * z_em) / (pow(z_srt, 3) - 1.0 / z_em));
		}
			
		z_Deltanext = z_Delta - z_D / z_Dder;
		err = fabs(z_Deltanext - z_Delta) / (fabs(z_Delta) + DPREC);
		z_Delta = z_Deltanext;
	  } while(err > DPREC);
	}
	
	if (phase == 'A' || phase == 'P') {
		/* (?) What does this mean? */
	  do {
		z_k = 0;
		z_EI4 = ei4 = 0.0;
		do {
			/* Do an angle integral */
			z_ID[0][z_k] = z_ID[1][z_k] = 0.0;
			dtheta = M_PI / 2 / pow(3, z_k); /* (?) Same as M_PI / (2 * pow)? */
			z_PAN = 2 * pow(3, z_k - 1); /* number of new angles */

			if (z_k == 0)
				z_PAN = 1;

			/* (?) Create one single loop here? */
			for (pan_aux = 1; pan_aux <= z_PAN; pan_aux++) {
				pan = pan_aux + pan_aux / 2;
				theta = (pan - 0.5) * dtheta;
				for (int m = 0; m < z_M; m++){
					z_em = z_tt * (m + 0.5);
					z_srt = sqrt((z_em * z_em) + pow(z_Delta * sin(theta), 2));

					/* (?) Shorten the math here (make intermediate vals) */
					z_ID[0][z_k] += (
						1.5 * pow(sin(theta), 3) * dtheta * z_tt * z_Delta * (1.0 / z_srt - 1.0 / z_em)
					);
					z_ID[1][z_k] += (
						1.5 * pow(sin(theta), 3) * dtheta * z_tt * (z_em * z_em / pow(z_srt, 3) - 1.0 / z_em)
					);
				}
			}
			if (z_k != 0) {
				/* (?) What are we doing? */
				z_ID[0][z_k] += z_ID[0][z_k - 1] / 3.0;
				z_ID[1][z_k] += z_ID[1][z_k - 1] / 3.0;
			}
			err = 1.0;
			if ( z_k >= 3) {
				/* (?) What are we doing? */
				/* (?) Shorten the math, perhaps use a loop? */
				z_EI3 = (
					(ec3[4] * z_ID[0][z_k]) + (ec3[3] * z_ID[0][z_k - 1]) + (ec3[2] * z_ID[0][z_k - 2])
				);
				z_EI4 = (
					(ec4[4] * z_ID[0][z_k]) + (ec4[3] * z_ID[0][z_k - 1]) + (ec4[2] * z_ID[0][z_k - 2]) + (ec4[1] * z_ID[0][z_k - 3])
				);
				ei4 = (
					(ec4[4] * z_ID[1][z_k]) + (ec4[3] * z_ID[1][z_k - 1]) + (ec4[2] * z_ID[1][z_k - 2]) + (ec4[1] * z_ID[1][z_k - 3])
				);
				err = fabs(z_EI3 - z_EI4) / (fabs(z_EI4) + DPREC);
			}
			z_k++;
		} while(err > DPREC);

		z_D = z_Delta * log(z_tt) - z_EI4;
		z_Dder = log(z_tt) - ei4;
		z_Deltanext = z_Delta - z_D / z_Dder; /* (?) Do we want subtraction before division? */
		err = fabs(z_Deltanext - z_Delta) / (fabs(z_Delta) + DPREC);
		z_Delta = z_Deltanext;
	  } while(err > DPREC);
	}

	if (phase == 'D') {
		/* (?) What does this mean? */
	  do {
		z_k = 0;
		z_EI4 = ei4 = 0.0;
		do {
			/* angle integral */
			z_ID[0][z_k] = z_ID[1][z_k] = 0.0;
			dtheta = M_PI / 2 / pow(3, z_k); /* (?) Is this M_PI / (2 * pow) ? */
			z_PAN = 2 * pow(3, z_k - 1); /* number of new angles */
			if (z_k == 0)
				z_PAN=1;

			for (pan_aux = 1; pan_aux <= z_PAN; pan_aux++) {
				pan = pan_aux + pan_aux / 2; /* (?) Order of operations? */
				theta = (pan - 0.5) * dtheta;
				z_Y = sin(2 * theta);
				for (int m = 0; m < z_M; m++) {
					z_em = z_tt * (m + 0.5);
					z_srt = sqrt((z_em * z_em) + pow(z_Delta * z_Y, 2));

					z_ID[0][z_k] += (
						4 / M_PI * pow(z_Y, 2) * dtheta * z_tt * z_Delta * (1.0 / z_srt - 1.0 / z_em)
					);
					z_ID[1][z_k] += (
						4 / M_PI * pow(z_Y, 2) * dtheta* z_tt * ((z_em * z_em) / pow(z_srt, 3) - 1.0 / z_em) /* (?) Order of operations?? */
					);
				}
			}

			if (z_k != 0) {
				z_ID[0][z_k] += z_ID[0][z_k - 1] / 3.0;
				z_ID[1][z_k] += z_ID[1][z_k - 1] / 3.0;
			}
			err = 1.0;
			if (z_k >= 3) {
				z_EI3 = (
					(ec3[4] * z_ID[0][z_k]) + (ec3[3] * z_ID[0][z_k - 1]) + (ec3[2] * z_ID[0][z_k - 2])
				);
				z_EI4 = (
					(ec4[4] * z_ID[0][z_k]) + (ec4[3] * z_ID[0][z_k - 1]) + (ec4[2] * z_ID[0][z_k - 2]) + (ec4[1] * z_ID[0][z_k - 3])
				);
				ei4 = (
					(ec4[4] * z_ID[1][z_k]) + (ec4[3] * z_ID[1][z_k - 1]) + (ec4[2] * z_ID[1][z_k - 2]) + (ec4[1] * z_ID[1][z_k - 3])
				);

				err = fabs(z_EI3 - z_EI4) / (fabs(z_EI4) + DPREC);
			}
			z_k++;
		} while(err > DPREC);

		z_D = z_Delta * log(z_tt) - z_EI4;
		z_Dder = log(z_tt) - ei4;
		z_Deltanext = z_Delta - z_D / z_Dder; /* (?) OOO? */
		err = fabs(z_Deltanext - z_Delta) / (fabs(z_Delta) + DPREC);
		z_Delta = z_Deltanext;
	  } while(err > DPREC);
	}
	return z_Delta;
}
