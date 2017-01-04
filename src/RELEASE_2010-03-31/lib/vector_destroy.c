/* $Id: vector_destroy.c,v 1.3 2009/06/06 22:28:26 glk Exp $

	Purpose: Frees storage for vector structure.  */

#include "GRACEdefs.h"
#include "GRACEprototype.h"

#include <stdlib.h>

void vector_destroy( vector v ) 
{
 static int8_t SccsId[] = "$Id: vector_destroy.c,v 1.3 2009/06/06 22:28:26 glk Exp $";
  free ( v.value ) ;
}
