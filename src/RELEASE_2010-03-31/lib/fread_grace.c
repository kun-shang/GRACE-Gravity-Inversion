#include <stdio.h>
#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

static int8_t SccsId[] = "$Id: fread_grace.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

size_t fread_grace(void *ptr, size_t  size,  size_t  nitems,  FILE *stream)
/*----------------------------------------------------------------------------->
/ purpose:  identical function as standard c-function fread with the inclusion
/           of endian architecture check. If little-endian then bytes will
/           be swapped
/           endian architectures
/
/ coded by: Gerhard L.H. Kruizinga                08/27/01
/
<-----------------------------------------------------------------------------*/
{

   size_t  out;

   out = fread(ptr,size,nitems,stream); 

   if (little_endian()) swapbyte((int8_t *)ptr,size);

   return out;
}
