/****************************************************************************
*      RTG Source Code,                                                     *
*      Copyright (C) 1996, California Institute of Technology               *
*      U.S. Government Sponsorship under NASA Contract NAS7-1260            *
*                    (as may be time to time amended)                       *
*                                                                           *
*      RTG is a trademark of the California Institute of Technology.        *
*                                                                           *
*                                                                           *
*      written by Yoaz Bar-Sever, Willy Bertiger, Bruce Haines,             *
*                 Angelyn Moore, Ron Muellerschoen, Tim Munson,             *
*                 Larry Romans, and Sien Wu                                 *
****************************************************************************/
/* 
 
Purpose:
 
   This subroutine (SEConds to CALendar date) takes an input seconds past
   the Julian reference date (JDREF) for this library and returns the 
   components of the corresponding calendar date. The components of the
   calendar date are all returned as numbers to allow for use in computations
   For instance, the month is returned as the integer month number rather
   than as a character string.

   Adapted from timetrans SECCAL.F
   Bruce Haines  3/5/95

*/
#include "TimeLib.h"
#include <math.h>

void seccal(

/* Input: */

		real sec,   /* seconds past ref date */

/* Output */ 

                int32_t *year,  
      		int32_t *month,
		int32_t *day,
		int32_t *hour,
		int32_t *minute,
		int32_t *second,
		real *frac 
            )
{

/* Local Variables */

   static int8_t SccsId[] = "$Id: seccal.c,v 1.3 2009/06/06 22:25:01 glk Exp $";
   real jd;
   real temp; 

/* Method:
 
   Extract the fraction seconds (frac) from sec. Note that care must
   taken if sec < 0 since frac must be a non-negative number < 1. The o
   fraction seconds of the calendar representation is simply the frac
   part of the input seconds and can be computed immediately. The rem
   integral seconds is converted to the rest of the calendar date. Th
   done to avoid round off error that could be introduced by the
   intermediate conversion to Julian date. */

   *frac = sec - ((int32_t) sec);
   if ( *frac < 0.0 ) {
      *frac += 1;
   }

/* call sec2jd to convert the integral seconds to the Julian date. */

   jd = sec2jd( sec - *frac + 0.5);

/* call jd2cal to convert the Julian date to calandar date. */

   jd2cal( jd, year, month, day, hour, minute, second, &temp);
   
}
