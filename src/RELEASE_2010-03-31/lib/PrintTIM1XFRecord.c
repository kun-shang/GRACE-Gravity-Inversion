#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"
#include "GRACEio_utils.h"

static int8_t SccsId[] = "$Id: PrintTIM1XFRecord.c,v 1.4 2009/06/10 21:27:51 glk Exp $";

void PrintTIM1XFRecord(FILE *dst, TIM1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of the OBDH mapping to GPS Time
/          record to file pointer dst
/
/ coded by: Gerhard Kruizinga                   10/19/01
/
/ input:  *dst    Pointer to TIM1X Data Format File
/         *record Pointer to TIM1X Data struct (TIM1X_t)
<-----------------------------------------------------------------------------*/
{
  int8_t bits[8];
  
  int32_t i;
 
/*----------------------------------------------------------------------------->
/ Write record elements to dst
<-----------------------------------------------------------------------------*/

  fprintf(dst," %-25s = %d\n","record->obdh_time",record->obdh_time);
  fprintf(dst," %-25s = %lc\n","record->GRACE_id",record->GRACE_id);
  fprintf(dst," %-25s = %d [","record->TS_suppid",record->TS_suppid);
/*----------------------------------------------------------------------------->
/ Decode TS_suppid 
<-----------------------------------------------------------------------------*/
  GetCharBits(record->TS_suppid,bits);
  loop(i,8)fprintf(dst,"%d",bits[7-i]);  
  fprintf(dst,"]\n");

  fprintf(dst," %-25s = %d\n","record->gpstime_intg",record->gpstime_intg);
  fprintf(dst," %-25s = %d\n","record->gpstime_frac",record->gpstime_frac);
  fprintf(dst," %-25s = %d\n","record->first_icu_blknr",record->first_icu_blknr);
  fprintf(dst," %-25s = %d\n","record->final_icu_blknr",record->final_icu_blknr);
  fprintf(dst," %-25s = %d [","record->qualflg",record->qualflg);
/*----------------------------------------------------------------------------->
/ Decode qualflag 
<-----------------------------------------------------------------------------*/
  GetCharBits(record->qualflg,bits);
  loop(i,8)fprintf(dst,"%d",bits[7-i]);  
  fprintf(dst,"]\n");

}
