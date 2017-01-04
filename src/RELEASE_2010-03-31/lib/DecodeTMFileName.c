#include <stdio.h>
#include <string.h>
#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"
#include "TimeLib.h"

static int8_t SccsId[] = "$Id: DecodeTMFileName.c,v 1.5 2009/06/06 22:15:38 glk Exp $";

#define MAXCHAR 1000

int32_t DecodeTMFileName(int8_t *TMFilename_in , int8_t *Satellite, int8_t *Version,
                      int8_t *DateString, int32_t *DateT2000, int8_t *Station,
                      int8_t *ProcCenter,int8_t *DataStream,int8_t *DataType)
/*----------------------------------------------------------------------------->
/ purpose:  Decode Telemetry file and return contents
/
/ coded by: Gerhard L.H. Kruizinga             06/06/01
/
/ input: TMFilename           Standard Telemetry Filename
/ output:
/        int8_t satellite       Satellite indicator (1 = GRACE A and 2 = GRACE B"
/        int32_t Version         Version string of file for Date
/        int8_t DateString      Timestamp in yyyy-mm-dd hh:mm:ss
/        int32_t DateT2000       Timestamp in seconds past 01/01/2000 12:00:00
/        int8_t Station         Downlink Station Name (eg NZ = Neustrelitz)
/        int8_t ProcCenter      Processing center Name (eg RDC)
/        int8_t DataStream      Satellite Data Stream (eg RT)
/        int8_t DataType        Date type (eg HK SC)
/
/ return 0L   if succesful decoding
/        1L   if invalid satnam
/        2L   if invalid product level
/        3L   if invalid year
/        4L   if invalid day of year
/        5L   if invalid hour
/        6L   if invalid minute
/        7L   if not all items have been found in tm filename
/
/ example of Telemetry filename : GR1-0-RDC-RT-HK+NZ_2003_158_10_24_1_2
/ GRACE-A application packets provided by RDC in Level-0 format, extracted from
/ Telemetry Packets of RT-HK received at Neustrelitz ground station on 158th
/ day, 2003 (dump start at 10:24), version 1.2
<-----------------------------------------------------------------------------*/
{
   int32_t  year,month,day,hour,minute,second,doy;

   int32_t i,level,item_count,word_count;

   double frac,Time,Tjan01;

   int8_t SatName[MAXCHAR],TMFilename[MAXCHAR];
   int8_t temp[MAXCHAR];

   second = 0;
   frac   = 0.0;

   strcpy(TMFilename,TMFilename_in);
   strcat(TMFilename,"_na");

   item_count = 1;
   word_count = 0;

   loop(i,strlen(TMFilename)-1)
   {
     if (strncmp((TMFilename+i),"-",1) == 0 || 
         strncmp((TMFilename+i),"_",1) == 0 ||
         strncmp((TMFilename+i),".",1) == 0 ||
         strncmp((TMFilename+i),"+",1) == 0)  

     {
       switch (item_count)
       {
         case 1:
                temp[word_count] = '\0';
                strcpy(SatName,temp);
                if (strncmp(SatName,"GR1",3)!=0 && strncmp(SatName,"GR2",3)!=0)
                return 1L;
                break;
         case 2:
                temp[word_count] = '\0';
                if (sscanf(temp,"%d",&level) != 1) return 2L;
                break;
         case 3:
                temp[word_count] = '\0';
                strcpy(ProcCenter,temp);
                break;
         case 4:
                temp[word_count] = '\0';
                strcpy(DataStream,temp);
                break;
         case 5:
                temp[word_count] = '\0';
                strcpy(DataType,temp);
                break;
         case 6:
                temp[word_count] = '\0';
                strcpy(Station,temp);
                break;
         case 7:
                temp[word_count] = '\0';
                if (sscanf(temp,"%d",&year) != 1) return 3L;
                break;
         case 8:
                temp[word_count] = '\0';
                if (sscanf(temp,"%d",&doy) != 1) return 4L;
                break;
         case 9:
                temp[word_count] = '\0';
                if (sscanf(temp,"%d",&hour) != 1) return 5L;
                break;
         case 10:
                temp[word_count] = '\0';
                if (sscanf(temp,"%d",&minute) != 1) return 6L;
                break;
       }
       if (item_count < 11)
       {
         item_count++;
         word_count = 0;
       }
       else
       {
         temp[word_count] = TMFilename[i];
         word_count++;
       }
     }
     else
     {
       temp[word_count] = TMFilename[i];
       word_count++;
     }
   }

   if (item_count < 11) return 7L;

   temp[word_count] = '\0';
   strcpy(Version,temp);

   Tjan01 = calsec(year,1,1,0,0,0,0.0);
   Time   = Tjan01 + (double) (doy - 1)*86400.0 + (double) hour * 3600.0 +
                     (double) minute * 60.0;

   seccal(Time,&year,&month,&day,&hour,&minute,&second,&frac);

   sprintf(DateString,"%4d-%02d-%02d %02d:%02d:%02d",year,month,day,hour,
           minute,second);

   *DateT2000 = (int32_t) Time;

   strncpy(Satellite,(SatName+2),1);
   if (Satellite[0] == '1') Satellite[0] = 'A';
   if (Satellite[0] == '2') Satellite[0] = 'B';
   Satellite[1] = '\0';

   loop(i,strlen(Version)-1) if (strncmp((Version+i),".",1)==0) Version[i]='\0';

   return 0L;
}
