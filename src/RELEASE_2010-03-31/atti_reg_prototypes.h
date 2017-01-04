/*  $Id: atti_reg_prototypes.h,v 1.12 2004/08/09 16:21:38 scw Exp $ */

#ifndef _atti_reg_prototypes_h_
#define _atti_reg_prototypes_h_

#include "atti_sim.h"

int32_t regres(

/* Input: */
	int32_t	restart,	/* if = 1, restart from beginning of data */
	int32_t	dT,		/* measurement time interval (sec) */
	double	tstart,		/* ignore any measurements before this time */
	int32_t     mode,		/* = 0 if epoch-state dynamic partials
				   = 1 if current-state dynamic partials */
/* Input & Output */
	atti_state  *x,		/* state vector; initial conditions at first call */
	atti_vari   *vp,	/* variational partials; initial conditions at first call */
	double	    *cam_off,	/* camera rotation angle offset (3 angles x 2 cameras) */
        atti_sim_t        *SimInfo,
/* Output */
	double	*tc,		/* current time */
	int8_t	*GRACE_id,	/* GRACE satellit id ('A' or 'B')*/
	int32_t	*sca_id,	/* star camera id number */
				/* = 0 if no measurements at current time */
	double	weight[4],	/* data weights to be applied */
	double	res[4],		/* residuals (O - C) quaternion */
	double	prtl[][4],	/* measurement partials */
	double	dQis_dP[4][3],	/* partial matrix of Qis wrt. Projection components */
	double	Qis[4],		/* model Qis quaternion */
	double	Qic[4]		/* model Qic quaternion */
/* Return 0 if no more measurement; 1 otherwise */
	  );

int32_t next_quat(
 
/* Input: */
	int32_t	RW_flg,		/* if = 1 then rewind SCA1A file to first meas. and return */
	int32_t	dT,		/* data interval (only data at N*dT will be output) */
        atti_sim_t        *SimInfo,
/* Output */
	double	*t,		/* time tag of quat */
	double	*q,		/* SCA1A quaternion components */
	int8_t	*GRACE_id,	/* GRACE id ('A' or 'B') */
	int32_t	*camera,	/* star camera id number */
	double	*noise_scale 	/* some function of nstars, residuals to be determined */
/* Return 0 if no more measurement; 1 otherwise */
         );

int32_t dQi2c_dQi2s
        (
/* Input: */
        double  t,              /* time past J2000 (sec) */
        int8_t    GRACE_id,       /* GRACE satellite id ('A' or 'B')*/
        int32_t     sca_id,         /* star camera id number */
	double	Qc_cp[2][4],	/* quaternion accounting for rotational offset of camera frame */
				/* from true to pseudo camera frames, for 2 cameras */
/* Output */
        double  Qs2c[4],        /* srf to star camera frame quaternion */
        double  prtl[4][4]      /* dQic_dQis partials */
/* Return 0 if normal; 1 otherwise */
        );

int32_t write_quat2sca1a(
/* Input: */
        double  time,             /* time past J2000 (sec) */
        double *Qic,              /* quaternion inertial to camera frame */
        int32_t     sca_id,           /* star camera id number */
        atti_sim_t *AttiSimSetup, /* Struct with simulation setup parameters includes GRACE_id */
        int32_t CloseFileFlag        /* =1 update header and close SCA1A file =0 write record */
        );
#endif
