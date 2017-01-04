
/*

------------------------------------------------------------------------

    Purpose: transformation between time and space coordinates
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

        double aod_open (char *file_aod, int nmax_aod, double *aod_eph);

     Global variables:

        extern int NMAJ_PER, NMAJ_OTS, NALL_OTS;
        extern double *OTADM;
        extern OTSperturb *MAJ_PER;
        extern OTSconstit *ALL_OTS;



------------------------------------------------------------------------



*/




#ifndef _GRVTS_H_

#include "grvts.h"
    #define _GRVTS_H_
#endif



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double aod_open (char *file_aod, int nmax_aod, double *aod_eph)
{
    FILE *fp_aod;
    double c00, s00, c06, s06, c12, s12, c18, s18;
    int n,m, l, ind, par_aod;
    char string[400];

    if ((fp_aod = fopen (file_aod,"r")) == NULL)
    {
        printf ("Cannot open AOD file?\n");
        exit (0);
    }

    par_aod = (nmax_aod + 1) * (nmax_aod + 1);

    aod_eph[0 * (par_aod + 1)] = gps2tt(0);
    aod_eph[1 * (par_aod + 1)] = gps2tt(21600);
    aod_eph[2 * (par_aod + 1)] = gps2tt(21600 * 2);
    aod_eph[3 * (par_aod + 1)] = gps2tt(21600 * 3);

    while (1)
    {
        if (fgets (string, 400, fp_aod) == NULL) break;
        sscanf (string, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf", 
            &n, &m, &c00, &s00, &c06, &s06, &c12, &s12, &c18, &s18);    

        if (n > nmax_aod )  continue;
        else if (m == 0)
        {
            aod_eph[0 * (par_aod + 1) + 1 + n] = c00;
            aod_eph[1 * (par_aod + 1) + 1 + n] = c06;
            aod_eph[2 * (par_aod + 1) + 1 + n] = c12;
            aod_eph[3 * (par_aod + 1) + 1 + n] = c18;
        }
        else 
        {
            l = nmax_aod - m + 1;
            ind = nmax_aod + 1 + (2 * nmax_aod - m + 2) * (m - 1);
            aod_eph[0 * (par_aod + 1) + 1 + ind + n - m] = c00;
            aod_eph[1 * (par_aod + 1) + 1 + ind + n - m] = c06;
            aod_eph[2 * (par_aod + 1) + 1 + ind + n - m] = c12;
            aod_eph[3 * (par_aod + 1) + 1 + ind + n - m] = c18;
            aod_eph[0 * (par_aod + 1) + 1 + ind + n - m + l] = s00;
            aod_eph[1 * (par_aod + 1) + 1 + ind + n - m + l] = s06;
            aod_eph[2 * (par_aod + 1) + 1 + ind + n - m + l] = s12;
            aod_eph[3 * (par_aod + 1) + 1 + ind + n - m + l] = s18;
        }
    }

    fclose(fp_aod);
    return 0;
}


