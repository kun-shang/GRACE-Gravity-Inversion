/* $Id: matrix_destroy.c,v 1.3 2009/06/06 22:28:26 glk Exp $
	Purpose:
           Frees storage for matrix structure. 
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

#include <stdlib.h>

void matrix_destroy( matrix v ) 
{
 static int8_t SccsId[] = "$Id: matrix_destroy.c,v 1.3 2009/06/06 22:28:26 glk Exp $";
  free ( v.value[0] ) ;
  free ( v.value ) ;
}
