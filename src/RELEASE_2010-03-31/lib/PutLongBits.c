#include "GRACEiolib.h"

static int8_t SccsId[] = "$Id: PutLongBits.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

/*******************************************************************
/Entry Point: PutLongBits.c
/
/Description:
/   This function takes n bits in the input variable starting
/   at position px, and puts them into the output variable starting
/   at position py.  The call is PutLongBits(x, px, &y, py, n).
/   For example, PutLongBits(x, 3, &y, 6, 2) puts
/   2 bits starting at position 3 of x and puts them into the starting
/   position 6 of y.  In other words, bits 3 and 2 of x are put into
/   bits 6 and 5 of y.  Note that the ADDRESS of y has to be passed.
/
/Input variables:
/   x  input int32_t from which the bits are to be taken
/   px starting position of the bits, with 0 being the lsb.
/   n  number of bits to put.
/Output variables:
/   y  output int32_t where the bits are put
/   py starting position of the bits, with 0 being the lsb.
/
/Written by: R. Berwin
/
/Modification history:
/
/Author     Date       Revision
/--------------------------------------------------------------------
/
/
/
/********************************************************************/
int32_t            PutLongBits(uint32_t x, int32_t px,
                   uint32_t *y, int32_t py, int32_t n)
{
   uint32_t   a, b, c;

   if ((px + 1 - n) < 0 || (py + 1 - n) < 0)
      return -1;
   else {
      a = b = c = 0;
      /* Move n bits of x starting at px to right justified of a */
      a = x >> (px + 1 - n) & ~(~0 << n);
      a = a << (py + 1 - n);    /* Move bits to position py */
      /* Move n bits of y starting at py to right justified of b */
      b = *y >> (py + 1 - n) & ~(~0 << n);
      b = b << (py + 1 - n);    /* Move bits to position py */
      c = b ^ *y;               /* exclusive OR */
      *y = c | a;               /* insert bits into y */
   }
   return 0;
}

