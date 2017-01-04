#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


static int8_t SccsId[] = "$Id: PrintTHR1XFRecord.c,v 1.8 2009/06/10 21:27:51 glk Exp $";

void PrintTHR1XFRecord(FILE *dst, THR1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of Thruster activation Level 1B 
/          Data Format record to file pointer dst
/
/ coded by: J. E. Patterson                  06/21/00
/
/ modified name of routine and structure     02/15/01
/
/ input:  *dst    Pointer to THR1X Data Format File
/         *record Pointer to THR1X Data struct (THR1X_t)
<-----------------------------------------------------------------------------*/
{

  int32_t i; 

  int8_t label[100];
  int8_t string[3];


/*----------------------------------------------------------------------------->
/ Write Header to dst
<-----------------------------------------------------------------------------*/

  fprintf(dst," %-21s = %d\n","record->time_intg",record->time_intg);
  fprintf(dst," %-21s = %d\n","record->time_frac",record->time_frac);
  
  strcpy(string,"-");
  string[0] = record->time_ref;
  fprintf(dst," %-21s = %s\n","record->time_ref",string);

  strcpy(string,"-");
  string[0] = record->GRACE_id;
  fprintf(dst," %-21s = %s\n","record->GRACE_id",string);

  loop(i,MAXTHRSTRS)
  {
    fprintf(dst," %-18s%02d%1s = %d\n","record->thrust_count[",i+1,"]",
                record->thrust_count[i]);
  }

  loop(i,MAXTHRSTRS)
  {
    fprintf(dst," %-18s%02d%1s = %d\n","record->on_time[",i+1,"]",
                record->on_time[i]);
  }

  loop(i,MAXTHRSTRS)
  {
    fprintf(dst," %-18s%02d%1s = %d\n","record->accum_dur[",i+1,"]",
                record->accum_dur[i]);
  }

  fprintf(dst," %-21s = %d\n","record->qualflg",record->qualflg);

}
