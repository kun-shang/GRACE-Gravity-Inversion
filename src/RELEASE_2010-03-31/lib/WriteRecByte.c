#include <stdio.h>
#include <string.h>
#include "GRACEiolib.h"
#include "GRACEgpslib.h"

static int8_t SccsId[] = "$Id: WriteRecByte.c,v 1.5 2009/06/06 22:14:36 glk Exp $";

void WriteRecByte(FILE* dst, int32_t RecBytes)
/*----------------------------------------------------------------------------->
/ purpose: write leading or trailing 4 bytes in fortran unformatted binary file
/
/ coded by: G.L.H. Kruizinga        04/07/01
/
/ input:  *dst      Pointer to fortran unformatted binary file
/ output: *RecBytes Pointer to Integer to return number of bytes in record
/
/-----------------------------------------------------------------------------*/
{
  if (fwrite(&RecBytes,sizeof(int32_t),1,dst) != 1)
  {
    fprintf(stderr,"Error writing RecBytes =%d\n",RecBytes);
    exit(1);
  }
}
