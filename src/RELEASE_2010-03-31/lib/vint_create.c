/* $Id: vint_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $
	Purpose:
           Creates storage for vint structure. Storage is 
           _not_ initialized to 0.

Initial coding:
   09/29/95 RJM 
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

#include <stdlib.h>		/* malloc */
#include <stdio.h>		/* printf */

vint vint_create( 
/* Input:    */
          int32_t size
/* Output: function value */
 ){

 static int8_t SccsId[] = "$Id: vint_create.c,v 1.3 2009/06/06 22:28:26 glk Exp $";
  vint vi; 

  vi.size = size;
  vi.value = (int32_t * ) malloc( size * sizeof(int32_t) );
  if ( ! vi.value )  fprintf( stderr ,"Error in vint_create\n" ) ;

  return vi;
}
