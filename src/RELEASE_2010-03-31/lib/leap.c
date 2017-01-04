#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TimeLib.h" 


#define CFILNMAX   80
#define MAXCHAR  1000

static int8_t SccsId[] = "$Id: leap.c,v 1.4 2009/11/18 21:29:59 glk Exp $";

table LoadLeapSeconds()

{
/*----------------------------------------------------------------------------->
/
/   purpose: load all leap seconds into table from file
/
/   output:   leap.x[]   time of leap second occurance [seconds past J2000 UTC]
/             leap.y[]   leap seconds offset
/
<-----------------------------------------------------------------------------*/

  FILE    *leap_file;

  double   sec_frac,seconds;

  int8_t     filename[CFILNMAX],line[MAXCHAR];
  int8_t     *ref_month[] = {"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG",
                           "SEP","OCT","NOV","DEC"};
  int8_t     month[MAXCHAR],dummy[MAXCHAR];

  int32_t      nleap;
  int32_t     mon,day,year,hour,min;
  table    leap;

/*----------------------------------------------------------------------------->
/ Open leap second file
<-----------------------------------------------------------------------------*/

/*
 strcpy(filename,"/time-pole/LEAPSECS");
*/

 strcpy(filename,"/GIPSY_source/goa-var/time-pole/LEAPSECS");

 leap_file = fopen(filename, "rb");

 if (leap_file == NULL)
 {
  fprintf(stderr,"\n Leap Second file %s cannot be opened !! \n\n",filename);
  exit(1);
 }

/*----------------------------------------------------------------------------->
/ Load leap seconds into table
<-----------------------------------------------------------------------------*/

 nleap = 0;

 while (fgets(line,MAXCHAR,leap_file) != NULL) nleap++;

 rewind(leap_file);

 leap = table_create(nleap);

 nleap = 0;

 while (fgets(line,MAXCHAR,leap_file) != NULL) 
 {
  sscanf((line),"%2d%1s%3s%1s%4d%3d%1s%2d%1s%7lf",&day,&dummy,month,&dummy,
                    &year,&hour,&dummy,&min,&dummy,&seconds);

  mon      = 1;
  while ( strncmp(month,ref_month[mon-1],3) != 0) mon++;

  sec_frac = seconds - (double)( (int32_t)seconds);

  leap.x[nleap] = calsec (year,mon,day,hour,min,(int32_t) seconds,sec_frac);
  leap.y[nleap] = -atof((line+38));

  nleap++;
 }

 fclose(leap_file);
 return leap;

}
