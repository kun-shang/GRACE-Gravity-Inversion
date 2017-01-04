#include "GRACEiolib.h"
#include "GRACEio_utils.h"

static int8_t SccsId[] = "$Id: UnSetCharBit.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

void UnSetCharBit(uint8_t *value, int32_t nth_bit)
/*----------------------------------------------------------------------------->
/ purpose:  set n-th bit to zero in value
/
/ coded by: Gerhard L.H. Kruizinga                07/13/00
/ modified: Gerhard L.H. Kruizinga                10/08/01
/
<-----------------------------------------------------------------------------*/
{
  uint16_t x;

  int32_t           px,py,n;

  x  = 0;
  px = (int32_t) nth_bit;
  py = px;
  n  = 1;

  if (nth_bit < 0 || nth_bit > 7) 
  {
    fprintf(stderr,"\n Invalid bit location in UnSetShortBit specified = %d \n",
                    nth_bit);
    fprintf(stderr," Bit Range should be 0 <= nth_bit <= 7\n\n");
  }

  if (PutCharBits(x,px,value,py,n) == -1)
  {
    fprintf(stderr,"\n Invalid call to PutShortBits in UnSetShortBit\n",
                    nth_bit);
    fprintf(stderr," Check call statement !!\n\n");
  }
}
