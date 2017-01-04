#include "GRACEiolib.h"
#include "GRACEio_utils.h"
#include "GRACEio_prototypes.h"

#define Failure 0


#define NBITSMAX 16

static int8_t SccsId[] = "$Id: WriteKBR1AFRecord.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

boolean WriteKBR1AFRecord(FILE *dst, KBR1A_t *record)
/*----------------------------------------------------------------------------->
/ purpose: write KBR data format level 1A record to file pointer dst
/          using generic routine WriteGFD1XFRecord.c
/
/ input:  *dst    pointer to GPS Flight Data Format File
/ output: *record Pointer to GPS Flight Data struct (KBR1A_t)
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
