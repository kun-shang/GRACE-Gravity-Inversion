/*  $Id: atti_prototypes.h,v 1.8 2004/08/30 20:14:57 glk Exp $  */ 

#ifndef _atti_prototypes_h_
#define _atti_prototypes_h_

#include "atti_sim.h"

VEC *atti_rhs(                    /* return the RHS of 1st order system
                                    of ODE's */

/* Input: need to pass S, B as VEC? */
            double t,
            VEC *x,             /* currently a pointer to a 6 vector
                                   proj quat, angular velocity omega
                                   or proj quat , omega, partials wrt
                                   initial proj quat, 3 scale factors
                                   and 3 bias parameters = dim 6 +
                                   6*12 + 3*6 , Parameters for camera alignment could 
                                   go in the model par of the code,O-C

                                   Order of *x = pq,omega,
                                   dpq/dParm[0], domega/dParm[0],...
                                   for each parmater where the
                                   parameters are pq0, omega0
                                   [epoch state ] and S, a scale
                                   vector and B is a Bias vector for
                                   the angular accelerations
                                   domega/dParam = 0 unless Param = S
                                   or B so these are not necessary

                                   Also note that most of domega/dS
                                   and domega/dB are 0 except for the
                                   diagonal, but I've left those in
                                   for now. Might simplify indexing


                                */
/* Output pointer to a vector of values */
            VEC *x_dot,

            double *SB           /* accelerometer scale, bias , scale then bias */
);

int32_t delta_ij(                   /* dirac integer delta functiion */

/* Input: */
            int32_t i,    
            int32_t j
/* Output function value 0 or +-1 Levi-Cevita permutation , -999999 on error*/
            );

int32_t eps_ijk(

/* Input: */
            int32_t i,              /* 0, 1, or 2 c indexing */
            int32_t j,
            int32_t k
/* Output function value 0 or +-1 Levi-Cevita permutation , -999999 on error*/
            );

double rk4_m(
             VEC *(*f)(),
             double t,
             VEC *x,
             double h,
             double *aux
             );

int32_t atti(                       /* returns 0 unless error */
 
/* Input */
         double t,              /* desired value at output */
 
/* Input/Output */
         atti_state *state,     /* contains both initial condition,
                                   and previous time point on input;
                                   current values on output */
         atti_vari *partl       /* same comments as state , see struct */
 
         );


int32_t vari_nm( char* nout[][12]);

int32_t guess_epoch_state (
/* Input: */
           double epoch_time,   /* epoch to return p,w in GPS J2000 time */
           atti_sim_t *SimInfo, /* for passing to next_quat in sim mode */
/* Input: */
           double *p,           /* projective quaternion at epoch_time */ 
           double *w,           /* angular velocity                    */
           int8_t *GRACE_id_out,  /* GRACE S/C id (A || B)               */
           int8_t sca_desig       /* sca designator (P || S)             */
           );

int32_t ang_acc(                    /* return angular measured/interpolater acceleration at t */

/* Input: */
            double t,           /* input time GPS J2000 sec
                                   assumes t strictly increasing for now */
                                /* ENV ANG_ACC_FILE environement variable with name of file */
            double *SB,         /* accelerometer scale, bias , scale then bias */
/* Output */
            double *ang_acc,      /* radians/sec**2 */
            double *MB,           /* mean angular acc bias of all measurements
                                     after 3-sigma editing for x,y,z axis */
            int8_t *Grace_id        /* GRACE S/C id (A || B) */

            );

int32_t write_state2memory(
/* input*/
         atti_state *state      /* contains both initial condition,
                                   and previous time point on input;
                                   current values on output */

/* return 0  successful write
          1  time tag out of range
*/
         );

int32_t read_state_memory(
/* input*/
         atti_state *state      /* contains both initial condition,
                                   and previous time point on input;
                                   current values on output */

/* return 0  successful read
          1  time tag out of range
*/
         );


#endif
