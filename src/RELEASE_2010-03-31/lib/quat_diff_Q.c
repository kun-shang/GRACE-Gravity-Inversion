/*   $Id: quat_diff_Q.c,v 1.4 2009/06/06 22:28:26 glk Exp $"


   Purpose:
          output the difference of quaternion Q from quaternion P
	  in terms of difference quaternion

   11/30/2001	Sien Wu		created 
   03/12/2004	Sien Wu		swap sequence of quat_prod wo that output is equivalent to
				(P)inv x (Q) instead of (Q) x (P)inv

*/
#include "GRACEdefs.h"
#include "GRACEprototype.h"

void quat_diff_Q(

/* Input: */
        quaternion	Q,	/* target quaternion */
        quaternion	P,	/* reference quaternion */

/* Output: */
        quaternion     *D       /* difference quaterion of Q from P */
           )
{
/* Local: */
  static int8_t SccsId[] = "$Id: quat_diff_Q.c,v 1.4 2009/06/06 22:28:26 glk Exp $";


  quaternion PN;	/* inverse quaternion of P */

  PN.q0 = -P.q0;
  PN.q1 =  P.q1;
  PN.q2 =  P.q2;
  PN.q3 =  P.q3;

  quat_prod( PN, Q, D);

  return ;
}
