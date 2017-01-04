#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "GRACEiolib.h"
#include "GRACEgpslib.h"
#include "TimeLib.h"

static int8_t SccsId[] = "$Id: prn2gpsh.c,v 1.6 2009/07/13 18:49:24 wib Exp $";

#define MAXCHAR   1000
#define MAXPRN2GPS  10 
#define MAXPRNSAT  101 
#define NOGPSID     -1

int32_t prn2gpsh(int32_t prn, double T2000)
/*----------------------------------------------------------------------------->
/ purpose: return gps number based on prn and time tag 2000
/
/ coded by: G.L.H. Kruizinga        04/10/01
/
/ input:  int32_t prn     prn number 
/         double T2000 time in sec past 2000 for which gps number is requested
/
/ output: Return gps number if valid requestion otherwise 0 is returned.
/
/ Note: Enviroment variable GOA_VAR must be set for PRN_GPS
/-----------------------------------------------------------------------------*/
{
  FILE             *prngps;

  static int32_t      first = 1;

  int8_t             filename[MAXCHAR],line[MAXCHAR];
  int8_t             start_date[MAXCHAR],deact_date[MAXCHAR];
  int8_t             *p;

  int32_t             prn_id,gps_id;
  int32_t             i,j,gps_out,ndx;

  int32_t              year,month,day;

  static double    sec_frac = 0.0,No_time = -999999999.0;

  static double    start_time[MAXPRNSAT][MAXPRN2GPS];
  static double    final_time[MAXPRNSAT][MAXPRN2GPS];

  static int32_t      prn2gps_count[MAXPRNSAT]; 
  static int32_t      prn2gps[MAXPRNSAT][MAXPRN2GPS];

  static int32_t       hour= 0,minute = 0,second = 0;

  if (first == 1)
  {
    loop(i,MAXPRNSAT)
    {
      loop(j,MAXPRN2GPS) 
      { 
        start_time[i][j] = No_time;
        final_time[i][j] = No_time;
        prn2gps[i][j] = NOGPSID;
      }
      prn2gps_count[i] = 0;
    }
    strcpy(filename,"");
    strcat(filename,getenv("GOA_VAR"));
    strcat(filename,"/etc/PRN_GPS");

    prngps = fopen(filename,"r");
    if (prngps == NULL)
    {
     fprintf(stderr,"\n PRN_GPS file  %s cannot be opened !! \n\n", filename);
     return 0L;
    }

    if (fgets(line,MAXCHAR,prngps) == NULL)
    {
     fprintf(stderr,"\n PRN_GPS file  %s is empty !! \n\n", filename);
     return 0L;
    }

    while (fgets(line,MAXCHAR,prngps) != NULL)
    {
      sscanf(line,"%s %s %d %d",start_date,deact_date,&gps_id,&prn_id);

      if (prn_id < 0 || prn_id > MAXPRNSAT-1)
      {
        fprintf(stderr,"\n Prn = %d is out of range. should be 0<->100\n\n",
                       prn_id);
        exit(1);
      }
      if (prn2gps_count[prn_id] + 1 > MAXPRN2GPS) 
      {
        fprintf(stderr,"\n Number of gps2prn's exceeds %d in file %s !!\n",
                        MAXPRN2GPS,filename);
        fprintf(stderr,
                " Change MAXPRN2GPS in prn2gpsh.c to accomodate entries\n\n");
        exit(1);
      }

      sscanf(start_date,"%4d-%2d-%2d",&year,&month,&day);
 
      prn2gps[prn_id][prn2gps_count[prn_id]] = gps_id;

      start_time[prn_id][prn2gps_count[prn_id]] = 
                         calsec(year,month,day,hour,minute,second,sec_frac);
      
      if (strcmp(deact_date,"0000") != 0) 
      {
        sscanf(deact_date,"%4d-%2d-%2d",&year,&month,&day);
        final_time[prn_id][prn2gps_count[prn_id]] = 
                         calsec(year,month,day,hour,minute,second,sec_frac);
      }

      prn2gps_count[prn_id]++;
    }

    first = 0;
    fclose(prngps);

/*
    loop(i,MAXPRNSAT)
    {
      if (prn2gps_count[prn_id] != 0)
      {
        loop(j,prn2gps_count[i])
        {
          fprintf(stderr,"prn ts tf : %d %f %f %d\n",i,start_time[i][j],
                  final_time[i][j],prn2gps[i][j]);
        }
        fprintf(stderr,"\n");
      }
    }
*/
  }

  if (prn < 0 || prn > MAXPRNSAT-1) 
  {
    fprintf(stderr,"\n Prn = %d is out of range. should be 0<->100\n\n",
                    prn_id);
    return 0L;
  }

  gps_out = 0;

  ndx = 0;

  while (T2000 > start_time[prn][ndx] && ndx < prn2gps_count[prn])
  {
    gps_out = prn2gps[prn][ndx]; 
    if (final_time[prn][ndx] != No_time)
      if (T2000 > final_time[prn][ndx]) gps_out = 0;
    ndx++;
  }

  return gps_out;
}
