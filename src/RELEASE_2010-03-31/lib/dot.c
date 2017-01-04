/* $Id: dot.c,v 1.5 2009/06/06 22:28:26 glk Exp $

   Purpose: Return dot product of two vectors.
*/

#include "GRACEdefs.h"
#include "GRACEprototype.h"

real dot( vector x, vector y )
{
 static int8_t SccsId[] = "$Id: dot.c,v 1.5 2009/06/06 22:28:26 glk Exp $";

 real d = 0;		/* accumulates dot product */
 int32_t  i;		/* index */
 for (i = 0 ; i < x.size ; i++ ) d += x.value[i] * y.value[i];
 return d;
}
