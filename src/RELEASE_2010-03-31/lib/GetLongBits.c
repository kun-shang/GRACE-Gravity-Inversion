#include "GRACEiolib.h"

#define NBITSMAX 32 

static int8_t SccsId[] = "$Id: GetLongBits.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

void GetLongBits(int32_t value, int8_t *bits)
/*----------------------------------------------------------------------------->
/ purpose:  return bits in int32_t value
/
/ coded by: Gerhard L.H. Kruizinga                11/13/98
/ borrowed from Gerhard and modified by Jean      06/22/00
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  for ( i = 31 ; i >= 0; i-- )
  {
    bits[i] = (value & (1<<i)) ? 1 : 0;
  }
}
