#include "GRACEiolib.h"

#define NBITSMAX 16

static int8_t SccsId[] = "$Id: GetShortBits.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

void GetShortBits(int16_t value, int8_t *bits)
/*----------------------------------------------------------------------------->
/ purpose:  return bits in int16_t value
/
/ coded by: Gerhard L.H. Kruizinga                11/13/98
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  for ( i = 15 ; i >= 0; i-- )
  {
    bits[i] = (value & (1<<i)) ? 1 : 0;
  }
}
