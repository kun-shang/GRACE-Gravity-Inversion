#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiXXXVOFRecord.c,v 1.6 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiXXXVOFRecord(FILE *dst, XXXVO_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records of Vector Orientation data 
/          from the file pointed to by dst 
/
/ coded by: J. E. Patterson                  07/18/00
/
/ modified name of routine and structure     02/15/01                           
/
/ input:  *dst    Pointer to XXXVO Data Format File
/         *record Pointer to XXXVO Data struct (XXXVO_t)
<-----------------------------------------------------------------------------*/
{
 int8_t string[3];
 int8_t bits8[8];                                                                         

 int32_t i;

 GetCharBits(record->qualflg,bits8);                                                    

 strcpy(string,"-");
 string[0] = record->GRACE_id;                                                          

 fprintf(dst,"%d %s %.16g %.16g %.16g %.16g", 
         record->gps_time, string, record->mag, record->cosx, 
         record->cosy, record->cosz);

 fprintf(dst," ");                                                                     
 loop(i,8)fprintf(dst,"%d",bits8[7-i]);                                                 
 fprintf(dst,"\n");
 
  return Success;
}
