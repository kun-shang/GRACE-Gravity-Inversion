/* $Id: stats.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

#ifndef _stats_h_
#define _stats_h_

#include <stdio.h>              /* def NULL */
#include <math.h>               /* fabs */
#include <inttypes.h>
#include <stdlib.h>


typedef struct stats_t{
  int32_t n;
  double v[4];		        /* acccumulated statistic values */
                               /*  stats_add */
    
} stats_t;

stats_t* stats_init();

void stats_reset( stats_t *x );

void stats_free( stats_t *x );

void stats_add( stats_t *x, double val );

void stats_calc( 
/* Input */
		stats_t *x ,	/* unchanged, more stats maybe accumulated by calling
                                   stats add  */
		int32_t *n,		/* number of points in stats */
		double s[5]	/* mean sigma min max rms  */
		);

int32_t stats_array( 
/* Input */

		double *a, 
		int32_t a_size,	/* 1 greater than the last index, number elements */ 
/* Output */
		int32_t *n,		/* a_size+1 */
		double s[5]	/* stats mean min max std_dev rms */
		);


int32_t SigEdit(
            /* Input */
            double *a,         /* array of doubles to edit */
            int32_t la,            /* number of elements in a last element
                                  a[la-1] */
            double SigFac,     /* outlier is more than SigFac*sigma
                                  away from mean */
            int32_t maxiter,       /* max iterations for sigma edit */
            stats_t *st,       /* if NULL an internal stats structure
                                 will be allocated and freed, if
                                 non-null it will be used and
                                 overwritten */
            /* output  */
            int32_t *edit_flgs);
     /* return number good points , 0 if all data is deleted, -1 if
        maxiter reached before converging */


#endif
