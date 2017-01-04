#include <stdio.h>
#include <string.h>
#include "GRACEiolib.h"
#include "GRACEgpslib.h"

static int8_t SccsId[] = "$Id: ReadRecByte.c,v 1.5 2009/06/06 22:14:36 glk Exp $";

void ReadRecByte(FILE* src, int32_t *RecBytes)
/*----------------------------------------------------------------------------->
/ purpose: read leading or trailing 4 bytes in fortran unformatted binary file
/          and return number of bytes
/
/ coded by: G.L.H. Kruizinga        04/06/01
/
/ input:  *src      Pointer to fortran unformatted binary file
/ output: *RecBytes Pointer to Integer to return number of bytes in record
/
/-----------------------------------------------------------------------------*/
{
  fread(RecBytes,sizeof(int32_t),1,src);
}
