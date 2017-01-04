#include <stdio.h>
#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1
#define Failure 0
#define MAXCHAR 1000

static int8_t SccsId[] = "$Id: ReadILG1XFRecord.c,v 1.4 2009/06/10 21:27:55 glk Exp $";

boolean ReadILG1XFRecord(FILE *src, ILG1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: read an IPU log message file from 
/          pointer src
/
/ coded by: Gerhard L.H. Kruizinga 10/15/00
/
/ input:  *src    pointer to ILG1X Data Format File
/ output: *record Pointer to ILG1X Data struct (ILG1X_t)
/
/ return:      1       normal return
/              0       End Of File reached
<-----------------------------------------------------------------------------*/
{
  int8_t          line[MAXCHAR], *c;
  int32_t          i, writing;

/*----------------------------------------------------------------------------->
/ Test for EOF
<-----------------------------------------------------------------------------*/

  if (fgets(line,MAXCHAR,src) == NULL) return False;

  if (sscanf(line,"%d %d %c >",&record->rcv_time,&record->pkt_count,
                                 &record->GRACE_id) != 3) return False;

  c = line;
  i = writing = 0;
  while ( (*c != '\0') && (*c != '\n') ) {
    if (writing) {
      if (i >= MAXLOGNAME-1) {
        fprintf(stderr,
          "Log packet line has > %d characters in ReadILG1XFRecord, truncated\n",
          MAXLOGNAME);
        record->logpacket[MAXLOGNAME-1] = '\0';
        return True;
      }
      record->logpacket[i++] = *c;
    }
    else if (*c == '>') writing = 1;
    c++;
  }
  record->logpacket[i] = '\0';
  
  return True;

}
