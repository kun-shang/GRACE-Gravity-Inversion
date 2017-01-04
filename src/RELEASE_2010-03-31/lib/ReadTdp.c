#include <stdio.h>
#include <string.h>
#include "GRACEiolib.h"
#include "GRACEgpslib.h"

static int8_t SccsId[] = "$Id: ReadTdp.c,v 1.5 2009/06/06 22:14:36 glk Exp $";

#define MAXCHAR 1000

int32_t ReadTdp(FILE* src, tdp_t *TdpRec)
/*----------------------------------------------------------------------------->
/ purpose: Read a tdp record from file src and put data in TdpRec struct
/
/ coded by: G.L.H. Kruizinga        04/16/01
/
/ input:  *src      Pointer to tdp file
/ output: *TdpRec   Pointer to TdpRec Struct
/
/ when end of file return 1L;
/-----------------------------------------------------------------------------*/
{

  int32_t i;

  int8_t line[MAXCHAR];

  if (fgets(line,MAXCHAR,src) == NULL) return 1L;

  sscanf(line,"%lf %lf %lf %lf %s",&TdpRec->time,&TdpRec->apriori,
                                   &TdpRec->value,&TdpRec->sigma,TdpRec->name);

  /* find pointer in line that points to TdpRec-> name */

  i = 0;
  while (strncmp(line+i,TdpRec->name,strlen(TdpRec->name)) != 0 && 
         i < strlen(line) ) i++;

  strncpy(TdpRec->name,line+i,TDPNAMELENGTH);

  TdpRec->name[TDPNAMELENGTH] = '\0';
 
  return 0L;
}
