#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define NBITSMAX 16

static int8_t SccsId[] = "$Id: PrintGFD1XFRecord.c,v 1.7 2009/06/10 21:27:48 glk Exp $";



void PrintGFD1XFRecord(FILE *dst, GFD1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of GPS Flight Data Format record 
/          to file pointer dst
/
/ coded by: Gerhard L.H. Kruizinga                  06/08/00
/ modified: Gerhard L.H. Kruizinga                  06/12/00
/ modified: J.E. Patterson                          08/29/00
/        changed definition of receiver time from
/        one double to 2 uint32_ts
/ modified: J.E. Patterson                          03/02/01
/        change record type name 
/
/ input:  *dst    pointer to GPS Flight Data Format File
/         *record Pointer to GPS Flight Data struct (GFD1X_t)
<-----------------------------------------------------------------------------*/
{

 int32_t i;

 int8_t bits[NBITSMAX];
 
/*----------------------------------------------------------------------------->
/ Write Header to dst
<-----------------------------------------------------------------------------*/
  fprintf(dst," %-20s = %d\n","record->rcvtime_intg",record->rcvtime_intg);
  fprintf(dst," %-20s = %d\n","record->rcvtime_frac",record->rcvtime_frac);
  fprintf(dst," %-20s = %c\n","record->GRACE_id",record->GRACE_id);
  fprintf(dst," %-20s = %d\n","record->prn_id",record->prn_id);
  fprintf(dst," %-20s = %d\n","record->ant_id",record->ant_id);
  fprintf(dst," %-20s = %d [","record->prod_flag",record->prod_flag);
/*----------------------------------------------------------------------------->
/ Decode Product flag 
<-----------------------------------------------------------------------------*/
  GetShortBits(record->prod_flag,bits);
  loop(i,NBITSMAX)fprintf(dst,"%d",bits[15-i]);  
  fprintf(stderr,"]\n");

  fprintf(dst," %-20s = %d\n","record->qualflg",record->qualflg);

/*----------------------------------------------------------------------------->
/ Write all product specified by prod_flag to dst
<-----------------------------------------------------------------------------*/

  loop(i,NBITSMAX)
  {
    if (bits[i] == 1)
    {
      switch(i)
      {
       case  0:
               fprintf(dst," %-20s = %f\n","record->CA_range",record->CA_range);
               break;
       case  1:
               fprintf(dst," %-20s = %f\n","record->L1_range",record->L1_range);
               break;
       case  2:
               fprintf(dst," %-20s = %f\n","record->L2_range",record->L2_range);
               break;
       case  3:
               fprintf(dst," %-20s = %f\n","record->CA_phase",record->CA_phase);
               break;
       case  4:
               fprintf(dst," %-20s = %f\n","record->L1_phase",record->L1_phase);
               break;
       case  5:
               fprintf(dst," %-20s = %f\n","record->L2_phase",record->L2_phase);
               break;
       case  6:
               fprintf(dst," %-20s = %d\n","record->CA_SNR",record->CA_SNR);
               break;
       case  7:
               fprintf(dst," %-20s = %d\n","record->L1_SNR",record->L1_SNR);
               break;
       case  8:
               fprintf(dst," %-20s = %d\n","record->L2_SNR",record->L2_SNR);
               break;
       case  9:
               fprintf(dst," %-20s = %d\n","record->CA_chan",record->CA_chan);
               break;
       case 10:
               fprintf(dst," %-20s = %d\n","record->L1_chan",record->L1_chan);
               break;
       case 11:
               fprintf(dst," %-20s = %d\n","record->L2_chan",record->L2_chan);
               break;
       case 12:
               fprintf(dst," %-20s = %f\n","record->K_phase",record->K_phase);
               break;
       case 13:
               fprintf(dst," %-20s = %f\n","record->Ka_phase",record->Ka_phase);
               break;
       case 14:
               fprintf(dst," %-20s = %d\n","record->K_SNR",record->K_SNR);
               break;
       case 15:
               fprintf(dst," %-20s = %d\n","record->Ka_SNR",record->Ka_SNR);
               break;
       default:
               fprintf(stderr,"\n Product Flag index %d is invalid!!!\n\n",i);
               exit(0);
      }
    }
  }
}
