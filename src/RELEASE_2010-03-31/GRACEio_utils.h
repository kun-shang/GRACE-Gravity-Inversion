/* $Id: GRACEio_utils.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

#ifndef _GRACEio_utils_h_
#define _GRACEio_utils_h_

#include "GRACEdefs.h"

/*----------------------------------------------------------------------------->
/  Structure definitions
<-----------------------------------------------------------------------------*/
/* None at the moment                                                         */
      
/*----------------------------------------------------------------------------->
/ Function Prototypes 
<-----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
void SetCharBit(
                uint8_t  *value,   /* set nth-bit to 1 in value         */
                int32_t           nth_bit
               );

/*----------------------------------------------------------------------------*/
int32_t PutCharBits(
                uint8_t x,         /* This function takes n bits        */
                int32_t px,                  /* in the input variable x starting  */
                uint8_t *y,        /* at position px, and puts them into*/
                int32_t py,                  /* the output variable y starting at */
                int32_t n                    /* position py                       */
               );
/* return:	0	normal return
               -1	incorrect arguments used for subroutine
*/

/*----------------------------------------------------------------------------*/
void SetShortBit(
                uint16_t  *value,  /* set nth-bit to 1 in value         */
                int32_t           nth_bit
               );
/*----------------------------------------------------------------------------*/
void UnSetShortBit(
                uint16_t  *value,  /* set nth-bit to 0 in value         */
                int32_t           nth_bit
               );
/*----------------------------------------------------------------------------*/
int32_t PutShortBits(
                uint16_t x,        /* This function takes n bits        */
                int32_t px,                  /* in the input variable x starting  */
                uint16_t *y,       /* at position px, and puts them into*/
                int32_t py,                  /* the output variable y starting at */
                int32_t n                    /* position py                       */
               );
/* return:	0	normal return
               -1	incorrect arguments used for subroutine
*/
/*----------------------------------------------------------------------------*/
void SetLongBit(
                uint32_t  *value,   /* set nth-bit to 1 in value         */
                int32_t           nth_bit
               );
/*----------------------------------------------------------------------------*/
void UnSetLongBit(
                uint32_t  *value,   /* set nth-bit to 0 in value         */
                int32_t           nth_bit
               );
/*----------------------------------------------------------------------------*/
int32_t PutLongBits(
                uint32_t x,         /* This function takes n bits        */
                int32_t px,                  /* in the input variable x starting  */
                uint32_t *y,        /* at position px, and puts them into*/
                int32_t py,                  /* the output variable y starting at */
                int32_t n                    /* position py                       */
               );
/* return:	0	normal return
               -1	incorrect arguments used for subroutine
*/

/*----------------------------------------------------------------------------*/
void GetCharBits
               (
                uint8_t   value,

                int8_t            *bits
               );

/*----------------------------------------------------------------------------*/
void GetShortBits
               (
                int16_t           value,

                int8_t            *bits
               );

/*----------------------------------------------------------------------------*/
void GetLongBits
               (
                int32_t		value,

                int8_t            *bits
               );

/*----------------------------------------------------------------------------*/
int32_t GetLongArgv(                         /* return int32_t after string in args */
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str                  /* search string                    */
               );
/* return:     int32_t after "str" in argument list to program
*/

/*----------------------------------------------------------------------------*/
void Get2LongArgv(                        /* return 2 longs after string      */
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str,                 /* search string                    */
               int32_t *quan1,               /* first int32_t after string          */
               int32_t *quan2                /* second int32_t after string         */
               );
/*----------------------------------------------------------------------------*/
int32_t GetIntArgv(
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str                  /* search string                    */
               );
/* return:     int32_t after "str" in argument list to program
*/

/*----------------------------------------------------------------------------*/
int8_t *GetStrArgv(
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str                  /* search string                    */
               );
/* return:     chars after "str" in argument list to program
*/

/*----------------------------------------------------------------------------*/
double GetDoubleArgv(
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str                  /* search string                    */
               );
/* return:     double after "str" in argument list to program
*/
/*----------------------------------------------------------------------------*/
int32_t FindOpt (
/* input */
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[],              /* arguments array                  */
               int8_t *str                  /* search string                    */
               );
/* return:     index of string "str" in argv[] arguments array
*/

int32_t GetArgEnvGraceVersion(
               int32_t argc,                  /* number of arguments              */
               int8_t *argv[]               /* arguments array                  */
               );
/* return:     return version number 0<= <= 99 or -1 upon errror
*/


#undef _mk_extern_

#endif	/* _GRACEio_utils_h_ */
