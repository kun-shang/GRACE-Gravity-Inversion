#include "TimeLib.h"

  /* Purpose:  This function will return the tdt time corresponding to
     the UTC time */

real utc2gps(
		real utc,             /* Input: utc is in seconds past J2000 */
                table  leap)          /* Input: leap second table */
  /* Output: real (return val)   GPS in seconds */
{
  /* leap seconds table contains TAI-UTC */
   static int8_t SccsId[] = "$Id: utc2gps.c,v 1.3 2009/06/06 22:25:01 glk Exp $";

   real r;
 
   r = (taiutc(utc, leap) + utc - 19.0);   /* This returns gps time. */
   return r;

} 


