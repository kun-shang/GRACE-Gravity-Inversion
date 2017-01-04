#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


static int8_t SccsId[] = "$Id: PrintXXXVOFRecord.c,v 1.6 2009/06/10 21:27:51 glk Exp $";


void PrintXXXVOFRecord(FILE *dst, XXXVO_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print an detailed ascii description of Vector Orientation 
/          Data Format record to file pointer dst
/
/ coded by: J. E. Patterson                  06/21/00
/
/ modified name of routine and structure     02/15/01                           
/
/ input:  *dst    Pointer to XXXVO Data Format File
/         *record Pointer to XXXVO Data struct (XXXVO_t)
<-----------------------------------------------------------------------------*/
{
 int8_t string[3];

 strcpy(string,"-");

/*----------------------------------------------------------------------------->
/ Write Header to dst
<-----------------------------------------------------------------------------*/

 fprintf(dst," %-20s = %d\n","record->gps_time",record->gps_time);
 string[0] = record->GRACE_id;
 fprintf(dst," %-20s = %s\n","record->GRACE_id",string);
 fprintf(dst," %-20s = %le\n","record->mag",record->mag);
 fprintf(dst," %-20s = %le\n","record->cosx",record->cosx);
 fprintf(dst," %-20s = %le\n","record->cosy",record->cosy);
 fprintf(dst," %-20s = %le\n","record->cosz",record->cosz);
 fprintf(dst," %-20s = %d\n","record->qualflg",record->qualflg);


}
