#include <stdio.h>
#include <string.h>
#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

static int8_t SccsId[] = "$Id: ConstructFileName.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

#define MAXCHAR 1000

void ConstructFileName(int8_t *Satellite , double Time, int32_t VersionNumber,
                       int8_t *FileTypeName, int8_t *Filename, int8_t *ext)
/*----------------------------------------------------------------------------->
/ purpose:  Construct filename based on input information
/
/ coded by: Gerhard L.H. Kruizinga             08/22/00
/
/ input: int8_t satellite       Satellite indicator (A = GRACE A and B = GRACE B"
/        double Time          time for filename time tag (sec past 2000)
/        int32_t VersionNumber   Version number of file for Date
/        int32_t FileTypeName    Standard file type name(aka filekey,eg ACC1B);
/        int8_t ext             file extension
/ output:
/        int8_t Filename        Filename constructed based on input according to
/                             TTTTT_yyyy_mm_dd_S_vn.ext
/                             where
/                             TTTTT file label based on File Type Pointer
/                             yyyy  year
/                             mm    month
/                             dd    day of month
/                             S     satellite indicator A=(GRACE A) B=(GRACE B)
/                             vn    version number
<-----------------------------------------------------------------------------*/
{
   int8_t Date[MAXCHAR];

   int32_t  year,month,day,hour,minute,second;

   double frac;

   seccal(Time,&year,&month,&day,&hour,&minute,&second,&frac);

   sprintf(Date,"%4d-%02d-%02d",year,month,day);

   sprintf(Filename,"%c%c%c%c%c_%s_%s_%02d.%s",
           FileTypeName[0],FileTypeName[1],FileTypeName[2],
           FileTypeName[3],FileTypeName[4],Date,Satellite,VersionNumber,ext);


}
