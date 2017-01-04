#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


static int8_t SccsId[] = "$Id: PrintGPS1AFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";


boolean PrintGPS1AFRecord(FILE *dst, GPS1A_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of GPS 1A Flight Data Format 
/          record to file pointer dst
/
/ input:  *dst    pointer to GPS 1A Flight Data Format File
/         *record Pointer to GPS 1A Flight Data struct (GPS1A_t)
<-----------------------------------------------------------------------------*/
{

  if (PrintGFD1X(dst, record))
  {
    return True;
  }
  else
  {
    return False;
  }
}
