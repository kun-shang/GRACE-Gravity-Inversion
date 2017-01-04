#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Failure 0

static int8_t SccsId[] = "$Id: WrAsciiIHK1XFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

boolean WrAsciiIHK1XFRecord(FILE *dst, IHK1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Write Ascii line the IPU HK 1X Data Format record to file pointer dst
/
/ coded by: Gerhard L.H. Kruizinga           09/25/01
/
/ input:  *dst    Pointer to IHK1XF Data Format File
/         *record Pointer to IHK1XF Data struct (IHK1X_t)
<-----------------------------------------------------------------------------*/
{
  int32_t n_len,i;
  int8_t bits[8];
 
/*----------------------------------------------------------------------------->
/ Decode Product flag 
<-----------------------------------------------------------------------------*/
  GetCharBits(record->qualflg,bits);
 
  fprintf(dst,"%d %d %c %c ",record->time_intg,record->time_frac, 
                            record->time_ref, record->GRACE_id);

  loop(i,8)fprintf(dst,"%d",bits[7-i]);  
    
  fprintf(dst," %c %.16g %s\n",record->sensortype, record->sensorvalue,
                               record->sensorname);

  return True;
}
