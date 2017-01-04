#include "GRACEiolib.h"
#include "GRACEio_utils.h"
#include "GRACEio_prototypes.h"


#define NBITSMAX 16

static int8_t SccsId[] = "$Id: ReadAHK1XFRecord.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

boolean ReadAHK1XFRecord(FILE *src, AHK1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: read Accelerometer housekeeping Format record from file pointer src
/          using generic routine ReadACC1AFRecord.c
/
/ input:  *src    pointer to ACC housekeeping file        
/ output: *record Pointer to AHK struct (AHK1X_t)
/
/ return:      1       normal return
/              0       End Of File reached
<-----------------------------------------------------------------------------*/
{

  if (ReadACC1AFRecord(src, record))
  {
    return True;
  }
  else
  {
    return False;
  }
}
