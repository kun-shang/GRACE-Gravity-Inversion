#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


#define Success 1 
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiPCI1AFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiPCI1AFRecord(FILE *dst, PCI1A_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records of the PCI Level 1A data 
/          from file pointer dst 
/
/ coded by: Gerhard Kruizinga                       03/18/03
/
/ input:  *dst    Pointer to PCI1A Data Format File
/         *record Pointer to PCI1A Data struct (PCI1A_t)
<-----------------------------------------------------------------------------*/
{
 int8_t bits8[8];
 
 int32_t i;

 GetCharBits(record->qualflg,bits8);


 fprintf(dst,"%d %c %.16g %.16g %.16g", record->gps_time, record->GRACE_id,
         record->ant_centr_corr, record->ant_centr_rate, record->ant_centr_accl);

 fprintf(dst,"  ");
 loop(i,8)fprintf(dst,"%d",bits8[7-i]);
 fprintf(dst,"\n");
 
 return Success;
}
