/* 

   Purpose:
          Convert seconds past J2000.0 to Julian Date

   11/07/95 Willy Bertiger

*/

#include "TimeLib.h" 

real sec2jd( 

/* Input: */
            real sec 
/* Output: function value */
            )

{

/* Local: */
  static int8_t SccsId[] = "$Id: sec2jd.c,v 1.3 2009/06/06 22:25:01 glk Exp $";

  return sec/sec_per_day + jdref ;

}






