/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _GRVTS_H_
    #define _GRVTS_H_


    #ifndef __STDIO_H__
        #include <stdio.h>
    #endif

    #ifndef __MATH_H__
        #include <math.h>
    #endif

    #ifndef __STRING_H__
        #include <string.h>
    #endif

    #ifndef __STDLIB_H__
        #include <stdlib.h>
    #endif

    #ifndef __CTYPE_H__
        #include <ctype.h>
    #endif

    #ifndef __TIME_H__
        #include <time.h>
    #endif

    #ifndef _NOVAS_H_
        #include "novas.h"
    #endif

    #ifndef _NOVASCON_H_
        #include "novascon.h"
    #endif

    #ifndef _NUTATION_H_
        #include "nutation.h"
    #endif

    #ifndef _EPHMAN_H_
        #include "eph_manager.h"
    #endif

    #ifndef _COORD_H_
        #include "coord.h"
    #endif


        typedef struct
        {
            double ds;
            int argn[6];
            char name[5];
            int n;
            int m;
            double cp;
            double sp;
            double cm;
            double sm;
            double fnm;
        }OTSperturb;

        typedef struct 
        {
            double ds;
            int argn[6];
            double ang;
            double cosf;
            double sinf;
        }OTSconstit;


        extern int NMAJ_PER, NMAJ_OTS, NALL_OTS;
        extern double *OTADM;
        extern OTSperturb *MAJ_PER;
        extern OTSconstit *ALL_OTS;




        int ots_open (char *ots_name, int nmax);
    
        double otidecs_admit(double jdt, double gmst, int nmax, int minor, 
                double *coef);

        int adm_open (char *adm_name);

        double arg2theta (double jdt, double gmst, int *n, 
                double *ang);



        double aod_open (char *file_aod, int nmax_aod, double *aod_eph);


        OTSperturb *ts_open (char *ts_name, int nmax, int *nper)
    
        int ts_read_minor(double jdt, double gmst, int nmax, int minor, 
                int nper, OTSperturb *per, double *coef)




#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


