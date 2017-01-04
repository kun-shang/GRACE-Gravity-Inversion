#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1 
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiCLK1BFRecord.c,v 1.10 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiCLK1BFRecord(FILE *dst, CLK1B_t *record)
/*-----------------------------------------------------------------------------'
/ purpose: Dump ascii records of Receiver Clock and GPS Clock data 
/          from the file pointed to by dst 
/
/ coded by: J. E. Patterson                         07/18/00
/
/ input:  *dst    Pointer to CLK1B Data Format File
/         *record Pointer to CLK1B Data struct (CLK1B_t)
'-----------------------------------------------------------------------------*/
{
 int8_t bits8[8];

 int32_t i;

 GetCharBits(record->qualflg,bits8);

 fprintf(dst,"%d %c %d %.16g %.16g %.16g %.16g",record->rcv_time, 
             record->GRACE_id,record->clock_id, record->eps_time, 
             record->eps_err, record->eps_drift,record->drift_err);

 fprintf(dst,"  ");
 loop(i,8)fprintf(dst,"%d",bits8[7-i]);
 fprintf(dst,"\n");
 
 return(Success);
}
