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

/*  $Id: TimeLib.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $  */ 

#include "GRACEdefs.h"

#ifndef _TimeLib_h_
#define _TimeLib_h_

#define jdref (2451545.0)      /* J2000.0 (January 1, 2000, 12 hours) */
#define sec_per_day (86400.0)

#define gps_strt (-630763200.0) /* seconds past J2000 
                                   on 06-Jan-1980 00:00:00.000 
                                   gps time is counted from here*/


void cal2ch(
/* Inputs*/ int32_t year,
            int32_t month,
            int32_t day,
            int32_t hour,
            int32_t minute,
            int32_t second,
            real frac_second,
            int16_t option_code, 
/* Outputs */
            int8_t *cal_str);

real calsec(int32_t year, int32_t month, int32_t day, int32_t hour, int32_t minute, int32_t second, 
  real frac );

void ch2cal(int8_t *string, int32_t *year, int32_t *month, int32_t *day, int32_t *hour, 
  int32_t *minute, int32_t *second, real *frac);

real ch2sec(int8_t *string);

int32_t date2j( int32_t year, int32_t month, int32_t day );

int32_t dayoyr(

/* Input: */
		int32_t year, 	/* The Year Number      */
		int32_t month, 	/* The Month Number     */
		int32_t day		/* The Day of the Month */
/* Output:  function return                             */
          );

real gps2tdt( real gps );

real gps2utc(real gps_time, table table_leap);

void gpslpsec(real utc_time, int32_t *lpsec, real *before, table leapsecs);

double gpsmin2sec( unsigned int   gps_minutes, uint16_t gps_millisecs );

double gpsws2sec(
  /* Input: */
	    int32_t gpsweek,	/* the gps week number */
	    double tow)	        /* seconds of GPS week */;

void j2date(

/* Input: */
		int32_t jd,      /* is the integer julian data */
/* Output */

                int32_t *year,
		int32_t *month,
		int32_t *day
	   );

void jd2cal (

/* Input */
		real jd,	/* Julian date       */
/* Output */
		int32_t *year,  	/* Year               */
		int32_t *month,	/* Month Number       */
		int32_t *day,	/* Day                */
		int32_t *hour,	/* Hour               */
		int32_t *minute,	/* Minute             */
		int32_t *second,	/* Second             */
		real *frac 	/* Fractional Seconds */
            );

real jd2sec(

/* Input */
		real jd	/* Julian date       */
/* Output */
/* Return value:
   real   - seconds past J2000 */
            );

void sec2ch(
/* Inputs */
           real tsec,        /* time, in seconds */
           int16_t opt_code,   /*0x0000: string is "DD-MMM HH:MM:SS.SSS" */
                             /*0x0001: string is "DD-MMM-YYYY HH:MM:SS.SSS"*/
/* Outputs */
           int8_t *cal_str)  /* calendar character string */;

void sec2gpsmin( double gpstime, 
                 unsigned int   *gps_minutes, uint16_t *gps_millisecs );

void sec2gpsws ( double sec, int32_t *week, double *tow );

real sec2jd( 
/* Input: */
            real sec 
/* Output: function value */
            );

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
            );

int8_t *stamper( 

/* Input: */
	      int16_t option_code,   /* option code */
/* Output: */
              int8_t *timestamp_str  /* note that stamper() returns a ptr*/
                                   /* to this string as well */ 
	      );

real taiutc(
  /* Input: */
	    real utc,		/* The value at which TAI-UTC is to be */
				/* determined.  utc is in seconds past J2000*/
	    table  leap)	/* leap second table */;

real tdt2et( real tdt );

real tdt2gps( real tdt );

real tdt2utc(
		real tdt,             /* Input: tdt is in seconds past J2000 */
                table  leap)          /* Input: leap second table */;

real utc2gps(
		real utc,             /* Input: utc is in seconds past J2000 */
                table  leap)          /* Input: leap second table */;

real utc2tdt(
		real utc,             /* Input: utc is in seconds past J2000 */
                table  leap)          /* Input: leap second table */;

table table_create(
                int32_t size );           /*                                     */

table LoadLeapSeconds();              /* load leap second table              */

#endif /* _TimeLib_h_ */
