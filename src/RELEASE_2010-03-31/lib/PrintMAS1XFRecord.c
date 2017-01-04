#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define NBITSMAX 8

static int8_t SccsId[] = "$Id: PrintMAS1XFRecord.c,v 1.5 2009/06/10 21:27:50 glk Exp $";

void PrintMAS1XFRecord(FILE *dst, MAS1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of SC Mass
/          Data Format record to file pointer dst
/
/ coded by: Gerhard L.H. Kruizinga                  09/25/01
/
/ input:  *dst    Pointer to MAS1X Data Format File
/         *record Pointer to MAS1X Data struct (MAS1X_t)
<-----------------------------------------------------------------------------*/
{

  int32_t i;

  int8_t bits[NBITSMAX];


/*----------------------------------------------------------------------------->
/ Write Header to dst
<-----------------------------------------------------------------------------*/

  fprintf(dst," %-21s = %d\n","record->time_intg",record->time_intg);
  fprintf(dst," %-21s = %d\n","record->time_frac",record->time_frac);
  fprintf(dst," %-21s = %c\n","record->time_ref",record->time_ref);
  fprintf(dst," %-21s = %c\n","record->GRACE_id",record->GRACE_id);
  fprintf(dst," %-21s = %d [","record->qualflg",record->qualflg);
  GetCharBits(record->qualflg,bits);
  loop(i,NBITSMAX)fprintf(dst,"%d",bits[NBITSMAX-1-i]);
  fprintf(stderr,"]\n");

  fprintf(dst," %-21s = %d [","record->prod_flag",record->prod_flag);
/*----------------------------------------------------------------------------->
/ Decode Product flag
<-----------------------------------------------------------------------------*/
  GetCharBits(record->prod_flag,bits);
  loop(i,NBITSMAX)fprintf(dst,"%d",bits[NBITSMAX-1-i]);
  fprintf(stderr,"]\n");

/*----------------------------------------------------------------------------->
/ priting all products specified by prod_flag 
/-----------------------------------------------------------------------------*/

  loop(i,NBITSMAX)
  {
    if (bits[i] == 1)
    {
      switch(i)
      {
         case 0:
                fprintf(dst," %-21s = %.16g\n","record->mass_thr",record->mass_thr);
                break;
         case 1:
                fprintf(dst," %-21s = %.16g\n","record->mass_thr_err",record->mass_thr_err);
                break;
         case 2:      
                fprintf(dst," %-21s = %.16g\n","record->mass_tnk",record->mass_tnk);
                break;
         case 3:      
                fprintf(dst," %-21s = %.16g\n","record->mass_tnk_err",record->mass_tnk_err);
                break;
         case 4:
                fprintf(dst," %-21s = %.16g\n","record->gas_mass_thr1",record->gas_mass_thr1);
                break;
         case 5:
                fprintf(dst," %-21s = %.16g\n","record->gas_mass_thr2",record->gas_mass_thr2);
                break;
         case 6:      
                fprintf(dst," %-21s = %.16g\n","record->gas_mass_tnk1",record->gas_mass_tnk1);
                break;
         case 7:      
                fprintf(dst," %-21s = %.16g\n","record->gas_mass_tnk2",record->gas_mass_tnk2);
                break;
          default:
          fprintf(stderr,"Product Flag index %d in PrintMAS1X is invalid!!! \n\n",i);
          exit(0);
      }
    }
  }
 
}
