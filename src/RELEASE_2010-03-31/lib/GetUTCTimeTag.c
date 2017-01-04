#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "GRACEiolib.h"

static int8_t SccsId[] = "$Id: GetUTCTimeTag.c,v 1.5 2009/06/06 22:15:38 glk Exp $";


void GetUTCTimeTag(int8_t *time_tag)
/*----------------------------------------------------------------------------->
/ purpose: produce UTC time tag at time of evaluation of this routine
/          in FileHeader_t struct
/
/ coded by: Gerhard L.H. Kruizinga                           08/15/2000
/
/ output: *time_tag  Pointer to time tag string
/
/
<-----------------------------------------------------------------------------*/
{
  time_t   now;
  int32_t len;

  int8_t time_tag1[HEADERMAXCHAR];
 
  now = time(NULL);

  strftime(time_tag,HEADERMAXCHAR,"%Y-%m-%d %H:%M:%S",gmtime(&now));
}
