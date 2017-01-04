#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

static int8_t SccsId[] = "$Id: PrintIOA1BFRecord.c,v 1.6 2009/06/10 21:27:49 glk Exp $";


void PrintIOA1BFRecord(FILE *dst, IOA1B_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of the Inertial Orientation 
/          of ACC Data Format record to file pointer dst
/
/ coded by: J. E. Patterson                  06/21/00
/
/ modified name of routine and structure     02/15/01
/
/ input:  *dst    Pointer to IOA1B Data Format File
/         *record Pointer to IOA1B Data struct (IOA1B_t)
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
 fprintf(dst," %-20s = %lf\n","record->quatangle",record->quatangle);
 fprintf(dst," %-20s = %lf\n","record->quaticoeff",record->quaticoeff);
 fprintf(dst," %-20s = %lf\n","record->quatjcoeff",record->quatjcoeff);
 fprintf(dst," %-20s = %lf\n","record->quatkcoeff",record->quatkcoeff);
 fprintf(dst," %-20s = %d\n","record->qualflg",record->qualflg);

}
