/* 

Purpose: 
 
   This subroutine (Julian date 2 calendar DATE) converts the input inte
   Julian date to the corresponding Gregorian calendar date. Since the J
   date is an integer, this correspondence is exact for noon of the cale
   date.
 
   The algorithm for this conversion is taken from the following article
   Tantzen,R.T., "Communications of the ACM", Volume 6, Number 8, August
   Algorithm 199, page 444.  

   Adapted from timetrans J2DATE.F
   Bruce Haines 3/5/96

*/

#include "TimeLib.h"

void j2date(

/* Input: */
		int32_t jd,      /* is the integer julian data */
/* Output */

                int32_t *year,
		int32_t *month,
		int32_t *day
	   )
{

/* Local: */
  static int8_t SccsId[] = "$Id: j2date.c,v 1.3 2009/06/06 22:25:01 glk Exp $";

  int32_t j;
  int32_t y;
  int32_t m;
  int32_t d;

  j = jd;
  j = j - 1721119;
  y = (4*j-1)/146097;
  j = 4*j - 1 - 146097*y;
  d = j/4;
  j = (4*d+3)/1461;
  d = 4*d + 3 -1461*j;
  d = (d+4)/4;
  m = (5*d-3)/153;
  d = 5*d - 3 - 153*m;
  d = (d+5)/5;
  y = 100*y + j;
  if (m < 10) {
     m = m + 3;
  } else {
     m = m - 9;
     y = y + 1;
  }
  *year   = y;
  *month  = m;
  *day    = d;

}
