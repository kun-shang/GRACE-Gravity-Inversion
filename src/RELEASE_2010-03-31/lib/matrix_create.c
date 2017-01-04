/* $Id: matrix_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $
	Purpose:
           Creates storage for matrix structure. Storage is 
           _not_ initialized to 0.

Initial coding:
   09/21/95 Willy Bertiger
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

#include <stdlib.h>		/* malloc */
#include <stdio.h>		/* printf */

matrix matrix_create( 
/* Input:    */
          int32_t row_dim,
          int32_t col_dim
/* Output: function value */
 ) 
{

/* Local Variables: */
 static int8_t SccsId[] = "$Id: matrix_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $";

  matrix r;			/* return value */
  int32_t    i;			/* loop index  */

  r.row_dim = row_dim;
  r.col_dim = col_dim;

/* first allocate the set of row pointers */

  if ( row_dim == 0 || col_dim == 0 ) {
    r.value = NULL;
    return r;
  }

  r.value = (real* * ) malloc( row_dim * sizeof(real *));

  r.value[0] = (real*) malloc( row_dim * col_dim * sizeof(real));
  if ( !r.value[0] ) fprintf( stderr ,"Error in matrix_create\n" ) ;

  for (i = 1; i < row_dim; i++) r.value[i] = r.value[0] + i * col_dim;

  return r;
}
