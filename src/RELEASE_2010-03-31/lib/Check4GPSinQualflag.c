#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"
#include "GRACEio_utils.h"

static int8_t SccsId[] = "$Id: Check4GPSinQualflag.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

int32_t Check4GPSinQualflag(uint8_t qualflg)
/*----------------------------------------------------------------------------->
/ purpose:  Return if time in qualflag is GPS time
/
/ coded by: Gerhard L.H. Kruizinga                10/08/01
/
/ return 1L if bit 0 = 0 (GPS time)
/        0L if bit 1 = 1 (Space Craft Elapsed Time SCET) 
<-----------------------------------------------------------------------------*/
{
  int8_t           bits8[8];

  GetCharBits(qualflg,bits8);

  if (bits8[0] == 0) return 1L;
  if (bits8[1] == 0) return 0L;
}
