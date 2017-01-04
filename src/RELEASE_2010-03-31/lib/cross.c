/* $Id: cross.c,v 1.3 2009/06/06 22:28:26 glk Exp $
   Purpose: cross product of two vectors

Initial coding:
   02/06/96	Sien Wu
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

void cross ( vector v1, vector v2, vector v3 ) 
{
  static int8_t SccsId[] = "$Id: cross.c,v 1.3 2009/06/06 22:28:26 glk Exp $";

  v3.value[0]=v1.value[1]*v2.value[2] - v1.value[2]*v2.value[1] ;
  v3.value[1]=v1.value[2]*v2.value[0] - v1.value[0]*v2.value[2] ;
  v3.value[2]=v1.value[0]*v2.value[1] - v1.value[1]*v2.value[0] ;
}
