#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define MAXLINECHAR 1000;

static int8_t SccsId[] = "$Id: WriteFileHeader.c,v 1.7 2009/11/16 21:33:30 glk Exp $";


boolean WriteFileHeader(FILE *dst,FileHeader_t *header)
/*----------------------------------------------------------------------------->
/ purpose: Write Header information from to file dst from information
/          in FileHeader_t struct
/
/ coded by: Gerhard L.H. Kruizinga                           08/15/2000
/
/ input:  *dst    Pointer to data file
/ output: *header Pointer to header struct
/
/ return:      1       normal return
/              0       error writing header
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  int8_t write_format[HEADERMAXCHAR];
  int8_t print_char[HEADERMAXCHAR+10];

  sprintf(write_format,"%%-%ds\n",HEADERMAXCHAR);

  rewind(dst);

  loop(i,header->nrecord)
  {
    strncpy(print_char,&header->HeaderCards[i][0],HEADERMAXCHAR-1);
    print_char[HEADERMAXCHAR-1] = '\0';
    fprintf(dst,write_format,print_char);
  }

  fprintf(dst,write_format,"END OF HEADER");

  return True;
}
