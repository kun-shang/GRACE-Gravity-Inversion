/****************************************************************************
*      RTG Source Code,                                                     *
*      Copyright (C) 1996, California Institute of Technology               *
*      U.S. Government Sponsorship under NASA Contract NAS7-1260            *
*                    (as may be time to time amended)                       *
*                                                                           *
*      RTG is a trademark of the California Institute of Technology.        *
*                                                                           *
*      $Id: calsec.c,v 1.4 2009/06/06 22:25:01 glk Exp $                                                                 *
*                                                                           *
*                                                                           *
*      written by Yoaz Bar-Sever, Willy Bertiger, Bruce Haines,             *
*                 Angelyn Moore, Ron Muellerschoen, Tim Munson,             *
*                 Larry Romans, and Sien Wu                                 *
****************************************************************************/
#include "TimeLib.h"  

real calsec(int32_t year, int32_t month, int32_t day, int32_t hour, int32_t minute, int32_t second, 
  real frac )
{
  /* Purpose:  This double precision function (CALendar date to SEConds) 
     takes the components of a calendar date and time
            Year / Month / Day , Hour : Minute : Second.Frac
     and returns the corresponding seconds past the reference date (jdref)
     this library. */

  /* Input:
     year                    is the year.
     month                   is the month number.
     day                     is the day.
     hour                    is the hour.
     minute                  is the minute.
     second                  is the second.
     frac                    is the fractional seconds. 

     Output: 
     real (return value)     is the seconds past J2000 (including fractions) */
  

  /* Declarations_of_Local_Variables */
  static int8_t SccsId[] = "$Id: calsec.c,v 1.4 2009/06/06 22:25:01 glk Exp $";

  real jd;                   /* Julian date */
  
  /* External functions */
  int32_t date2j(int32_t year, int32_t month, int32_t day);

  /* Call DATE2J to compute the double precision Julian date at the start
     the current day. */
  jd = date2j( year, month, day ) - 0.5;

  /* Call jd2sec to compute the seconds past the reference date at the st
  of the current day and then add in the seconds in the remaining frac
  day. */

  return ((jd - jdref)*sec_per_day + hour*3600.0 + minute*60.0 + second + frac);

}

