#include <stdio.h>
#ifndef _poly_h_
#define _poly_h_

/* $Id: poly.h,v 1.2 2009/06/03 22:54:20 glk Exp glk $ */

#define loop(A,B) for(A=0;A<B;A++)

typedef struct poly
        {
        double  x0;
        double* coef;         /* y = coef[0] + coef[1]*x + coef[2]*x^2 ... coef[order]*x^order  */
        int32_t   order;          /* number of elements of v */
        }  poly;

typedef struct freq
        {
        double x0;
        double* freq;         /* frequency array */
        int32_t   nfreq;          /* number of frequencies */
        double* fcoef;        /* y = fcoef[0]*cos(freq[0]*x    + fcoef[1]*sin(freq[0]*x  + ... 
                                     fcoef[2n-2]*cos(freq[n/2]*x fcoef[2n-1]*sin(freq[n/2]*x */
       
        }  freq;

poly poly_create(int32_t order);
freq freq_create(int32_t nfreq);
int32_t poly_lsq_fit(double *x, double *y, double *w,int32_t nobs,poly *poly_fit,
                  double *rms, double *rmsfit, freq *freq_fit);
double eval_poly(poly *poly_fit,double x, freq *freq_fit);
double print_poly(FILE *dst,poly *poly_fit,freq *freq_fit);


#endif /* _poly_h_ */
