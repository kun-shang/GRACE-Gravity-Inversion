/* Purpose: normalizing a vector */

#include "GRACEdefs.h"
#include "GRACEprototype.h"

void normalize ( vector a )
{
  static int8_t SccsId[] = "$Id: normalize.c,v 1.3 2009/06/06 22:28:26 glk Exp $";
  real norm ;
  int32_t  i ;

  norm = sqrt( dot( a, a ) );

  if (norm) for (i=0; i<a.size; i++) a.value[i] /= norm ;

  return ;
}
