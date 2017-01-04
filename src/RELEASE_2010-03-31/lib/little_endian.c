#include "GRACEiolib.h"

static int8_t SccsId[] = "$Id: little_endian.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

int32_t little_endian()
/*----------------------------------------------------------------------------->
/ purpose:  return if current architecture uses little endian convetion
/
/ coded by: Gerhard L.H. Kruizinga                08/24/01
/           Larry Romans
/
/ return: 1 Little endian architecture (e.g. Linux86 Vax) 
/         0 Big endian architecture (e.g. sun hpux)
/
/ In big-endian architectures, the leftmost bytes (those with a lower address) 
/ are most significant. In little-endian architectures, the rightmost bytes 
/ are most significant. 
/
/ The terms big-endian and little-endian are derived from the Lilliputians of 
/ Gulliver's Travels, whose major political issue was whether soft-boiled eggs 
/ should be opened on the big side or the little side. Likewise, the 
/ big-/little-endian computer debate has much more to do with political issues
/ than technological merits. 
<-----------------------------------------------------------------------------*/
{
  static int32_t retval = -1;
  static int8_t endi[] = "ENDI";
  int32_t        i, ldum;
  int8_t        *cp;

  if (retval < 0L) {

    cp = (int8_t *) &ldum;
    loop(i,4) cp[i] = endi[i];

    if (ldum != 1162757193) retval = 1L;
    else retval = 0L;
  }

  return(retval);

}
