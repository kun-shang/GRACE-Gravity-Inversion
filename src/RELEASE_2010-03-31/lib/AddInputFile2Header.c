#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"
#include "GRACEfiletype.h"
#include <string.h>

static int8_t SccsId[] = "$Id: AddInputFile2Header.c,v 1.8 2009/11/16 21:25:19 glk Exp $";

#define MAXCHAR 1000

boolean AddInputFile2Header(FileHeader_t *header,int8_t *filekey, 
                            int8_t *input_filename,int8_t *input_file_ttag,
                            int8_t *version, int8_t *linktime)
/*----------------------------------------------------------------------------->
/ purpose: Change Contents in header struct for RecPointer with Modstring
/
/ coded by: Gerhard L.H. Kruizinga                           08/28/2001
/
/ input:   *header          Pointer to header struct
/          *filekey         Filekey name of inputfile     
/          *input_filename  Input filename
/          *input_file_ttag Input filename time tag
/ output:  *header Pointer to header struct
/
/ return:      1       normal return
/              0       error modifying header
/
<-----------------------------------------------------------------------------*/
{
  int32_t i,RecPointer,ndx;

  int8_t write_format[HEADERMAXCHAR],String[HEADERMAXCHAR+10];
  int8_t Base_Filename[MAXCHAR];

  GetBaseFilename (input_filename, Base_Filename);

  /*>>>> sanity checks <<<<*/

  if ((header->nrecord+2) > NMAXHEADERREC)
  {
    fprintf(stderr,"\nNumber of total Header labels exceeds %d",NMAXHEADERREC);
    fprintf(stderr," in AddInputFile2Header\n\n");
    return False;
  }

  if (header->NinputFileLabel+2 > MAXINPUTFILELABEL)
  {
    fprintf(stderr,"\nNumber of Header input file labels exceeds %d",
                   MAXINPUTFILELABEL);
    fprintf(stderr," in AddInputFile2Header\n\n");
    return False;
  }

  if (header->nrecord < NRHEADERLABELS)
  {
    header->NinputFileLabel = 0;
    loop(i,NRHEADERLABELS) 
         strcpy(&header->HeaderCards[i][0],"NOT DEFINED");
    header->nrecord = NRHEADERLABELS;
  }


  RecPointer = header->nrecord;
  ndx        = header->NinputFileLabel;

  strcpy(header->InputFileLabel[ndx].filekey,filekey);
  strcpy(header->InputFileLabel[ndx].name,Base_Filename);
  strcpy(header->InputFileLabel[ndx].time_tag,input_file_ttag);
  if (version[0]) strcpy(header->InputFileLabel[ndx].software_version,version);
  if (linktime[0]) strcpy(header->InputFileLabel[ndx].linktime,linktime);

  sprintf(write_format,"%%-%ds: %%-%ds",HEADERLABELMAXCHAR,
          HEADERMAXCHAR-HEADERLABELMAXCHAR-2);

  sprintf(String,"%s<-%s",filekey,Base_Filename);
  sprintf(&header->HeaderCards[RecPointer][0],write_format,
          "INPUT FILE NAME               ",String);

  RecPointer++;
  
  sprintf(String,"%s<-%s",filekey,input_file_ttag);
  sprintf(&header->HeaderCards[RecPointer][0],write_format,
          "INPUT FILE TIME TAG (UTC)     ",String);

  header->nrecord += 2;

  if (version[0])
  {
    RecPointer++;
  
    sprintf(String,"%s<-%s",filekey,&version[5]);
    sprintf(&header->HeaderCards[RecPointer][0],write_format,
            "INPUT FILE SOFTWARE VERSION   ",String);
    header->nrecord++;
  }
  if (linktime[0])
  {
    RecPointer++;
  
    sprintf(String,"%s<-%s",filekey,&linktime[5]);
    sprintf(&header->HeaderCards[RecPointer][0],write_format,
            "INPUT FILE LINKTIME TAG       ",String);
    header->nrecord++;
  }

  header->NinputFileLabel++;

  return True;
}
