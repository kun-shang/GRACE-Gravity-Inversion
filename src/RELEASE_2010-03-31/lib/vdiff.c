/* $Id: vdiff.c,v 1.3 2009/06/06 22:28:26 glk Exp $

   Purpose: Difference of two vectors.
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

int32_t vdiff( vector v1, vector v2, vector v3 )
{
 static int8_t SccsId[] = "$Id: vdiff.c,v 1.3 2009/06/06 22:28:26 glk Exp $";

 int32_t  i;                       /* index */
 for (i = 0 ; i < v3.size ; i++ ) v3.value[i] = v1.value[i] - v2.value[i];
 return;
}
