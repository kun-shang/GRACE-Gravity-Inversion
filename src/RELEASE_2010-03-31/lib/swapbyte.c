#include "GRACEiolib.h"

static int8_t SccsId[] = "$Id: swapbyte.c,v 1.4 2009/06/06 22:15:38 glk Exp $";

int32_t swapbyte(int8_t buf[],size_t num_bytes)
/*----------------------------------------------------------------------------->
/ purpose:  swap bytes depinding on specification to accomadate different
/           endian architectures
/
/ coded by: Gerhard L.H. Kruizinga                08/24/01
/
<-----------------------------------------------------------------------------*/
{
   int8_t            temp;

   int32_t            i,half_num_bytes;

   half_num_bytes = (int32_t) num_bytes/2;
   
   loop(i,half_num_bytes)
   {
     temp               = buf[i];
     buf[i]             = buf[num_bytes-i-1];
     buf[num_bytes-i-1] = temp;
   }

   return 0L;
} 
