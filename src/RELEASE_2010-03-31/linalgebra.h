#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#
#ifndef _linalgebra_h_
#define _linalgebra_h_

/* $Id: linalgebra.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

#define loop(A,B) for(A=0;A<B;A++)

typedef struct vector
        {
        double* value;
        int32_t   size;          /* number of elements of v */
        }  vector;

typedef struct vint{
  int*   value;
  int32_t    size;
} vint;

typedef struct matrix
        {
        int32_t row_dim;            /* number of elements in a column */
        int32_t col_dim;            /* number of elements in a row */
        double** value;
        }  matrix;       


/* prototypes */

matrix matrix_create( 
/* Input:    */
          int32_t row_dim,
          int32_t col_dim
/* Output: function value */
 ); 

vector vector_create( 
/* Input:    */
          int32_t size
/* Output: function value */
 );

vint vint_create( 
/* Input:    */
          int32_t size
/* Output: function value */
 );

void matrix_destroy( matrix v );
void vector_destroy( vector v );
void vint_destroy( vint v );
void transpose(matrix *a,matrix *at);
void ab(matrix *a,matrix *b,matrix *c);
void atb(matrix *a,matrix *b,matrix *c);
void lcab(matrix *a,matrix *b,matrix *c,double fa,double fb);
void accum_normal_eq(matrix *H,double y,double weight,matrix *HTWH, matrix *HTWy);
void solve_normal_eq(matrix *HTWH, matrix *HTWy, matrix *covariance, matrix *solution);
void print_matrix(FILE *dst,matrix *a,int8_t *label);
void compute_correlation(matrix *cov,matrix *corr);

#endif /* _linalgebra_h_ */
