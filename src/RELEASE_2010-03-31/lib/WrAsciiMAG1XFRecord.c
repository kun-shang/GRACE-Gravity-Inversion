#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1 
#define Failure 0


static int8_t SccsId[] = "$Id: WrAsciiMAG1XFRecord.c,v 1.7 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiMAG1XFRecord(FILE *dst, MAG1X_t *record)
/*-----------------------------------------------------------------------------'
/ purpose: Dump ascii records of Magnetometer/Magnettorquer data 
/          from the file pointed to by dst 
/
/ coded by: J. E. Patterson                  07/18/00
/
/ modified name of routine and structure     02/15/01
/
/ input:  *dst    Pointer to MAG1XF Data Format File
/         *record Pointer to MAG1XF Data struct (MAG1X_t)
-----------------------------------------------------------------------------*/
{
 int8_t string1[3];
 int8_t string2[3];
 int8_t bits8[8];

 int32_t i;

 GetCharBits(record->qualflg,bits8);

 strcpy(string1,"-");
 strcpy(string2,"-");
 string1[0] = record->GRACE_id;                                                          
 string2[0] = record->time_ref;

 fprintf(dst,"%d %d %s %s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g", 
         record->time_intg, record->time_frac, 
         string1, string2, record->MfvX_RAW, record->MfvY_RAW,
         record->MfvZ_RAW, record->torque1A, record->torque2A,
         record->torque3A, record->torque1B, record->torque2B,
         record->torque3B, record->MF_BCalX, record->MF_BCalY,
         record->MF_BCalZ, record->torque_cal);

 fprintf(dst,"  ");
 loop(i,8)fprintf(dst,"%d",bits8[7-i]);
 fprintf(dst,"\n");
 
  return(Success);
}
