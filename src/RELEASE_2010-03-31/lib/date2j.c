#include "TimeLib.h"

int32_t date2j( int32_t year, int32_t month, int32_t day )
{
  static int8_t SccsId[] = "$Id: date2j.c,v 1.3 2009/06/06 22:25:01 glk Exp $";

  /* Purpose:
  This integer function (calendar DATE 2 Julian date) takes the input
  Gregorian calendar date and returns as its functional value the
  corresponding integer Julian date. Since the Julian date is an integer
  this correspondence is exact for noon of the output calendar date.

  The algorithm for this conversion is taken from the following article
  Tantzen,R.T., "Communications of the ACM", Volume 6, Number 8, August
  Algorithm 199, page 444. */

  /* Input_Arguments:

  year     is the year number.
  month    is the month number.
  day      is the day number.  */

  /* Output_Arguments */

  /* Declarations_of_External_Functions */
  /* none */

  /* Declarations_of_Local_Variables */

  int32_t Y;
  int32_t M;
  int32_t D;

  int32_t C;
  int32_t YA;
      
  Y = year;
  M = month;
  D = day;

  if( M > 2 )
    M = M - 3;
  else
  {
    M = M + 9;
    Y = Y - 1;
  }

  C  = Y/100;
  YA = Y - 100*C;
  return ((146097*C)/4 + (1461*YA)/4 + (153*M+2)/5 + D + 1721119);

}
