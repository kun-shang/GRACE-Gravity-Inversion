#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "GRACEfiletype.h"

#define loop(A,B) for(A=0;A<B;A++)

static int8_t SccsId[] = "$Id: PointerLib.c,v 1.4 2009/11/16 21:24:18 glk Exp $";

int32_t GetFileType (int8_t *txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return file type pointer based on text version of pointername
/
/ coded by:  Gerhard L.H. Kruizinga                      03/06/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  loop(i,NRFILETYPES)
  {
    if (strcmp(txt_pointer,FileTypeName[i]) == 0) return FileTypePointer[i];
  }

  fprintf(stderr,"\n FileTypePointer = %s has not been defined!!\n",txt_pointer);
  fprintf(stderr," Check input to GetFileType.\n\n");
  exit(1);
}
int32_t GetAcc1aProd (int8_t *txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return acc1a product pointer based on text version of pointername
/
/ coded by:  Gerhard L.H. Kruizinga                      07/26/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  loop(i,NRACC1APODS)
  {
    if (strcmp(txt_pointer,Acc1aProdName[i]) == 0) return Acc1aProdPointer[i];
  }

  fprintf(stderr,"\n Acc1aProdName = %s has not been defined!!\n",txt_pointer);
  fprintf(stderr," Check input to GetAcc1aProd.\n\n");
  exit(1);
}
int32_t GetHeaderLabel (int8_t *txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return file header label pointer based on text version of header
/          label name
/
/ coded by:  Gerhard L.H. Kruizinga                      03/06/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  loop(i,NRHEADERLABELS)
  {
    if (strcmp(txt_pointer,HeaderLabelName[i]) == 0) return HeaderLabelPointer[i];
  }

  fprintf(stderr,"\n HeaderLabelPointer = %s has not been defined!!\n",txt_pointer);
  fprintf(stderr," Check input to GetHeaderLabel.\n\n");
  exit(1);
}
int32_t GetFileTypeName (int32_t SectorPointer, char* txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return file type name based on file type pointer based 
/
/ coded by:  Gerhard L.H. Kruizinga                      03/06/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  strcpy(txt_pointer,"NOT_DEFINED");

  loop(i,NRFILETYPES)
  {
    if (SectorPointer == FileTypePointer[i]) 
    {
      strcpy(txt_pointer,FileTypeName[i]);
      return 0;
    }
  }

  fprintf(stderr,"\n SectorPointer = %d has not been defined!!\n",SectorPointer);
  fprintf(stderr," Check input to GetFileTypeName.\n\n");
  exit(1);
}
int32_t GetAcc1aProdName (int32_t SectorPointer, char* txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return acc1a product name based on pointer based  acc1a pointer
/
/ coded by:  Gerhard L.H. Kruizinga                      07/26/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  strcpy(txt_pointer,"NOT_DEFINED");

  loop(i,NRACC1APODS)
  {
    if (SectorPointer == Acc1aProdPointer[i]) 
    {
      strcpy(txt_pointer,Acc1aProdName[i]);
      return 0;
    }
  }

  fprintf(stderr,"\n Acc1aProdPointer = %d has not been defined!!\n",SectorPointer);
  fprintf(stderr," Check input to GetAcc1aProdName.\n\n");
  exit(1);
}
int32_t GetHeaderLabelName (int32_t HeaderRecPointer, int8_t *txt_pointer)
/*----------------------------------------------------------------------------->
/ purpose: return header label name based on file header label pointer 
/
/ coded by:  Gerhard L.H. Kruizinga                      03/06/01
/
<-----------------------------------------------------------------------------*/
{
  int32_t i;

  strcpy(txt_pointer,"NOT_DEFINED");

  loop(i,NRHEADERLABELS)
  {
    if (HeaderRecPointer == HeaderLabelPointer[i]) 
    {
      strcpy(txt_pointer,HeaderLabelName[i]);
      return 0;
    }
  }

  fprintf(stderr,"\n HeaderRecPointer = %d has not been defined!!\n",HeaderRecPointer);
  fprintf(stderr," Check input to GetHeaderLabelName.\n\n");
  exit(1);
}
