#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


static int8_t SccsId[] = "$Id: WrAsciiGPS1BFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";


boolean WrAsciiGPS1BFRecord(FILE *dst, GPS1B_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Dump ascii records of flight data from the file pointed to by dst 
/
/ coded by: J. E. Patterson                         07/18/00
/
/ input:  *dst    Pointer to GPS1B Data Format File
/         *record Pointer to GPS1B Data struct (GPS1B_t)
<-----------------------------------------------------------------------------*/
{

  if (WrAsciiGFD1XFRecord(dst, record))
  {
    return True;
  }
  else
  {
    return False;
  }
 
}
