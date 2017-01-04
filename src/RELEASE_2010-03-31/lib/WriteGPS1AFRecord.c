#include "GRACEiolib.h"
#include "GRACEio_utils.h"
#include "GRACEio_prototypes.h"

#define Failure 0


#define NBITSMAX 16

static int8_t SccsId[] = "$Id: WriteGPS1AFRecord.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

boolean WriteGPS1AFRecord(FILE *dst, GPS1A_t *record)
/*----------------------------------------------------------------------------->
/ purpose: write GPS Flight Data Format level 1A record to file pointer dst
/          using generic routine WriteGFD1XFRecord.c
/
/ input:  *dst    pointer to GPS Flight Data Format File
/ output: *record Pointer to GPS Flight Data struct (GPS1A_t)
/
/ return:      1       normal return
/              0       End Of File reached
<-----------------------------------------------------------------------------*/
{

  if (WriteGFD1XFRecord(dst, record))
  {
    return True;
  }
  else
  {
    return False;
  }
}
