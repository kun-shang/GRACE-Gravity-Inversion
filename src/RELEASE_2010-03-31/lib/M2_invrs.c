/* $Id: M2_invrs.c,v 1.3 2009/06/06 22:28:26 glk Exp $

   Purpose:
        Inversion of a 2x2 matrix

   06/18/2001   Sien Wu         Created
*/
#include "GRACEdefs.h"
#include "GRACEprototype.h"

int32_t M2_invrs   (
/* input */
		real	A[2][2],	/* input 2 x 2 matrix */
/* output */
		real    B[2][2]         /* inverse of A */

               ) /* return 1 if matrix is singular; 0 otherwise */
{
 static int8_t SccsId[] = "$Id: M2_invrs.c,v 1.3 2009/06/06 22:28:26 glk Exp $";

 real	det;	/* matrix determinant */

 det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
 if( det == 0 ) return 1;

 B[0][0] =  A[1][1] / det;
 B[0][1] = -A[0][1] / det;
 B[1][0] = -A[1][0] / det;
 B[1][1] =  A[0][0] / det;

return 0;
}
