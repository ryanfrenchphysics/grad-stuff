#ifndef _CONSTVAR_H
#define _CONSTVAR_H

/*
 * TODO: Insert an introduction to this file
 */

/* Preprocessor Constants */

#define Mk                  10  /* Number of stored vectors for OP update */  
#define N_rich              4   /* Richardson extrapolation number */
#define FErecursive_depth   17  /* recursion depth of F.E. lambda-integral */
#define PAN                 40
/*#define	LX                  30.0 */
#define NX                  50
#define NZ                  30
#define prec                1.0e-3
#define precM               1.0e-4
#define abs_prec            1.0e-4
#define LL                  140
#define tag_work            0
#define tag_relax           1
#define tag_stop            2
#define ENmax               1000
#define precR       (0.01 * prec)   /* Riccati amplitude convergence prec. */
#define precRfe     (0.01 * prec)   /* R.A. convergence prec. F.E. */
#define RfrG        (M_PI / 180)
#define GfrR        (180 / M_PI)
#define PI2         (2 * M_PI)
#define NNX         (NX) 
#define NNZ         (2 * NZ)
#define NMAX        ( (NNX > NNZ) ? NNX : NNZ )


/* Preprocessor datatypes */

#define YY1(x,y,z)      x       /*      Triplet      */
#define YY2(x,y,z)      y       /*      Order        */
#define YY3(x,y,z)      z       /*      Parameter    */


/* Preprocessor functions */

#define Mats(t) ( (int) 35.0 / t ) 

#endif /* _CONSTVAR_H */
