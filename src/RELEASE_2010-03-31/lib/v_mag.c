/* $Id: v_mag.c,v 1.4 2009/06/06 22:28:26 glk Exp $

   Purpose: compute magnitude of a vector
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

real v_mag ( vector a )
{
 static int8_t SccsId[] = "$Id: v_mag.c,v 1.4 2009/06/06 22:28:26 glk Exp $";

 real d = 0;		/* accumulates dot product */
 int32_t  i;		/* index */
 for (i = 0 ; i < a.size ; i++ ) d += a.value[i] * a.value[i];
 return sqrt( d );
}
