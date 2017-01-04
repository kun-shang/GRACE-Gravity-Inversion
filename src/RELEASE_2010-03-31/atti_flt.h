/* $Id: atti_flt.h,v 1.8 2004/11/01 22:33:35 wib Exp $ */

#ifndef atti_flt_h_
#define atti_flt_h_

#include "matrix.h"
#include "matrix2.h"

#include "atti.h"
#include "GRACEiolib.h"
#include "atti_utils.h" 
#include "atti_reg_prototypes.h"

#define MaxBatchMeas 6

typedef struct edit_criteria_t {
  double SigmaFactor;
  int32_t    SigWindow;
  double NominalSigTest_cam[3];
} edit_criteria_t;

typedef struct output_requests_t {
  int32_t out_every_n_sec;
  int32_t cov_cor;                  /* 0 = cov non-zero = correlation */
  int32_t iter_out;                 /* output only on iteration = iter ( start numbering at 1) */
} output_requests_t;

typedef struct proc_noise {

  int32_t Np;                       /* number of process noise states
                                   internal order*/
  double *M;                    /* array Np time correlation
                                   coefficients, from 0 to 1, function
                                   of tau ; M = exp( -dT / tau ) */
  double *Q;                     /* sqrt( 1 - M**2 ) * sigP; sigP =
                                   steady-state sigma P(j+1) = M P(j)
                                   + Q */

  double *tau;                  /* time correlation , tau = -1 =
                                   infinite = random walk; tau = 0 =>
                                   white noise , M = 0*/
  double *sigp;

} proc_noise;

typedef struct atti_flt {

  int32_t nDyn;                     /* number of dynamic states, was hard
                                   coded to 12 */
  int32_t n_NonDyn;                 /* number of non-dynamic states
                                   partials computed in regres */

  int32_t maxNonDyn;                /* size of estNonDyn and NonDynState */

  int32_t *estNonDyn;               /* array of len maxnondyn is 1 for est
                                   0 for not estimated, the number of
                                   1's = n_NonDyn */

  double *NonDynState;           /* maxNonDyn int32_t containing Val */

  int32_t maxDyn;                   /* size of estDyn */
  int32_t *estDyn;                 /* array of length ndyn either 1 or 0 
                                   1 = estimate ; exclude from filter */

  double *apsig;              /* apriori sigmas, dynamic then
                                   non-dynamic currently length is
                                   maxDyn + maxNonDyn*/

  double *apval;                /* currently not used in atti_lsq but
                                   it is in the SRIF */

  int32_t    *Param_order;          /* integer array of size maxNonDyn +
                                   maxDyn ( 18 )mapping the standard
                                   input order pq,omega, scale, bias,
                                   cam1 align, cam2 align to
                                   filter/smoother order which will
                                   always place the process noise
                                   paramters first (Np ) , a value >=
                                   Np + Nx indicates the parameter is
                                   not esitmated */

  int32_t   *Param_order_inv;       /* size Np+Nx maps these parmaters
                                   from silter smoother order back to
                                   input */

  int32_t Nx;                       /* number of non-stochastic states,
                                   this should always be >= 6 since
                                   the pointing (pq) and angular
                                   velocity are never stochastic */

  proc_noise *PN;                /* process noise characteristics */

} atti_flt;


#endif
