#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiOSCFQFRecord.c,v 1.8 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiOSCFQFRecord(FILE *dst, OSCFQ_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records from the Ultra Stable Oscillator 
/          data from the file pointed to by dst 
/
/ coded by: J. E. Patterson                  09/07/00
/
/ modified name of routine and structure     02/15/01
/
/ input:  *dst    Pointer to OSCFQ Data Format File
/         *record Pointer to OSCFQ Data struct (OSCFQ_t)
<-----------------------------------------------------------------------------*/
{
  int8_t string[3];
  int8_t bits8[8];

  int32_t i;

  GetCharBits(record->qualflg,bits8);

  strcpy(string,"-");
  string[0] = record->GRACE_id;                                                          

  fprintf(dst,"%d %s %d %.16g %.16g %.16g", record->gps_time,
      string, record->uso_id, record->uso_freq, record->K_freq,record->Ka_freq);

  fprintf(dst,"  ");
  loop(i,8)fprintf(dst,"%d",bits8[7-i]);
  fprintf(dst,"\n");
 
  return Success;
}
