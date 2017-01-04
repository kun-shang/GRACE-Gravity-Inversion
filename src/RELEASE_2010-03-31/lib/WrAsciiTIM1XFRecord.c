#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"
#include "GRACEio_utils.h"

#define Success 1
#define Failure 0

static int8_t SccsId[] = "$Id: WrAsciiTIM1XFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

boolean WrAsciiTIM1XFRecord(FILE *dst, TIM1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records from the OBDH time mapping to GPS time records to
/          file pointed to by dst 
/
/ coded by: Gerhard Kruizinga                       10/19/01
/
/ input:  *dst    Pointer to TIM1X Data Format File
/         *record Pointer to TIM1X Data struct (TIM1X_t)
<-----------------------------------------------------------------------------*/
{
 int8_t bits8[8];                                                                         

 int32_t i;

 GetCharBits(record->qualflg,bits8);                                                    

 fprintf(dst,"%d %c %d %d %d %d %d ", 
         record->obdh_time, record->GRACE_id, record->TS_suppid, 
         record->gpstime_intg, record->gpstime_frac, 
         record->first_icu_blknr, record->final_icu_blknr);

 loop(i,8)fprintf(dst,"%d",bits8[7-i]);                                                 
 fprintf(dst,"\n");
 
 return Success ;
}
