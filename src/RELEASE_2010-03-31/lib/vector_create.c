/* $Id: vector_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $
	Purpose:
           Creates storage for vector structure. Storage is 
           _not_ initialized to 0.

Initial coding:
   09/21/95 Willy Bertiger
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

#include <stdlib.h>		/* malloc */
#include <stdio.h>		/* printf */

vector vector_create( 
/* Input:    */
          int32_t size
/* Output: function value */
 ){

 static int8_t SccsId[] = "$Id: vector_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $";
  vector r; 

  r.size = size;
  r.value = (real * ) malloc( size * sizeof(real) );
  if ( !r.value )  fprintf( stderr ,"Error in vector_create\n" ) ;

  return r;
}
