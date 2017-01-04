#include "GRACEiolib.h"
#include "GRACEio_utils.h"

static int8_t SccsId[] = "$Id: UnSetLongBit.c,v 1.5 2009/06/06 22:15:38 glk Exp $";

void UnSetLongBit(uint32_t *value, int32_t nth_bit)
/*----------------------------------------------------------------------------->
/ purpose:  set n-th bit to zero in value
/
/ coded by: Gerhard L.H. Kruizinga                07/13/00
/ modified: Gerhard L.H. Kruizinga                08/23/00
/
<-----------------------------------------------------------------------------*/
{
  uint32_t x;

  int32_t           px,py,n;

  x  = 0;
  px = (int32_t) nth_bit;
  py = px;
  n  = 1;

  if (nth_bit < 0 || nth_bit > 31) 
  {
    fprintf(stderr,"\n Invalid bit location in UnSetLongBit specified = %d \n",
                    nth_bit);
    fprintf(stderr," Bit Range should be 0 <= nth_bit <= 31\n\n");
  }

  if (PutLongBits(x,px,value,py,n) == -1)
  {
    fprintf(stderr,"\n Invalid call to PutLongBits in UnSetLongBit\n",
                    nth_bit);
    fprintf(stderr," Check call statement !!\n\n");
  }
}
