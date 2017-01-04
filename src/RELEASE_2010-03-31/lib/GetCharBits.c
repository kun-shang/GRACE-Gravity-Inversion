#include "GRACEiolib.h"

static int8_t SccsId[] = "$Id: GetCharBits.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

void GetCharBits(uint8_t value, int8_t *bits)
/*----------------------------------------------------------------------------->
/ purpose:  return bits in int16_t value
/
/ coded by: Gerhard L.H. Kruizinga                11/13/98
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  for ( i = 7 ; i >= 0; i-- )
  {
    bits[i] = (value & (1<<i)) ? 1 : 0;
  }
}
