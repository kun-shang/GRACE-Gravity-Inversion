#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"


static int8_t SccsId[] = "$Id: PrintKBR1AFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";


boolean PrintKBR1AFRecord(FILE *dst, KBR1A_t *record)
/*----------------------------------------------------------------------------->
/ purpose: Print a detailed ascii description of KBR 1A Flight Data Format 
/          record to file pointer dst
/
/ input:  *dst    pointer to KBR 1A Flight Data Format File
/         *record Pointer to KBR 1A Flight Data struct (KBR1A_t)
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
