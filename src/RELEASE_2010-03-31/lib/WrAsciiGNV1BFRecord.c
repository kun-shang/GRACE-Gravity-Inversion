#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1 
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiGNV1BFRecord.c,v 1.8 2009/06/06 22:15:38 glk Exp $";     


boolean WrAsciiGNV1BFRecord(FILE *dst, GNV1B_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records of GPS Navigation Level 1B data 
/          from the file pointed to by dst 
/
/ coded by: J. E. Patterson                         07/18/00
/
/ input:  *dst    Pointer to GNV1BF Data Format File
/         *record Pointer to GNV1BF Data struct (GNV1B_t)
<-----------------------------------------------------------------------------*/
{
 int8_t string1[3];
 int8_t string2[3];
 int8_t bits8[8];

 int32_t i;

 GetCharBits(record->qualflg,bits8);

 strcpy(string1,"-");
 string1[0] = record->GRACE_id;                                                          

 strcpy(string2,"-");
 string2[0] = record->coord_ref;                                                          

 fprintf(dst,"%d %s %s %.16g %.16g %.16g %.16g %.16g"
        " %.16g %.16g %.16g %.16g %.16g %.16g %.16g",
         record->gps_time, string1, string2, record->xpos, 
         record->ypos, record->zpos, record->xpos_err, 
         record->ypos_err, record->zpos_err, record->xvel, 
         record->yvel, record->zvel, record->xvel_err, record->yvel_err, 
         record->zvel_err);

 fprintf(dst,"  ");
 loop(i,8)fprintf(dst,"%d",bits8[7-i]);
 fprintf(dst,"\n");
 
  return Success;
}
