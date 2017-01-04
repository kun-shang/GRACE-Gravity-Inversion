
/*

------------------------------------------------------------------------

    Purpose: transformation between time and space coordinates
    Notes: 
    Programmer: Kun Shang @ 4.29.2014
    Functions:

        OTSperturb *ts_open (char *ts_name, int nmax, int *nper)
    
        int ts_read_minor(double jdt, double gmst, int nmax, int minor, 
                int nper, OTSperturb *per, double *coef)

        int adm_open (char *adm_name);

        double arg2theta (double jdt, double gmst, int *n, 
                double *ang);

     Global variables:

        extern int nper, NMAJ_OTS, NALL_OTS;
        extern double *OTADM;
        extern OTSperturb *per;
        extern OTSconstit *ALL_OTS;



------------------------------------------------------------------------



*/




#ifndef _GRVTS_H_

#include "grvts.h"
    #define _GRVTS_H_
#endif




    int NMAJ_OTS = 18, NALL_OTS  = 256;
    double *OTADM;
    OTSconstit *ALL_OTS;



OTSperturb *ts_open (char *ts_name, int nmax, int *nper)
{
    FILE *fp_ot;
    int i, n;
    char string[200];

    OTSperturb *per;

    if ((fp_ot = fopen (ts_name,"r")) == NULL)
    {
        printf ("Cannot open tide file %s\n", ts_name);
        exit (0);
    }

    i = 0;
    while (1)
    {
        if (fgets (string, 200, fp_ot) == NULL) break;
        sscanf (string, "%*f%*s%d", &n);
        if (n > nmax) continue;
        i ++;
    }
    rewind(fp_ot);

    *nper = i;
    printf ("NPER = %d\n", *nper);
    per = (OTSperturb *) calloc ( *nper, sizeof(OTSperturb));
    if ( per == NULL )    
    {
        printf("Malloc error!\n");
    }

    i = 0;
    while (1)
    {
        if (fgets (string, 200, fp_ot) == NULL) break;
        sscanf (string, "%*f%*s%d", &n);
        if (n > nmax) continue;
        sscanf (string, "%lf%s%d%d%lf%lf%lf%lf", 
            &per[i].ds, per[i].name, &per[i].n, &per[i].m, 
            &per[i].cp, &per[i].sp, &per[i].cm, &per[i].sm);    

        per[i].argn[0] = (int)(per[i].ds/100)%10;
        per[i].argn[1] = (int)(per[i].ds/10)%10 - 5;
        per[i].argn[2] = (int)(per[i].ds/1)%10 - 5;
        per[i].argn[3] = (int)(per[i].ds*10)%10 - 5;
        per[i].argn[4] = (int)(per[i].ds*100)%10 - 5;
        per[i].argn[5] = (int)(per[i].ds*1000)%10 - 5;
        per[i].fnm = 1e-11;

        i++;
    }

    fclose(fp_ot);

    return per;

}






int adm_open (char *adm_name)
{
    FILE *fp_ad;    
    int i, j;

    double doodson[256] = {
    55.565,55.575,56.554,56.556,57.355,57.553,57.555,57.565,57.575,
    58.554,59.553,62.656,63.645,63.655,63.665,64.456,64.555,65.445,65.455,
    65.465,65.655,65.665,65.675,66.454,67.455,67.465,71.755,72.556,73.545,
    73.555,73.565,74.556,75.345,75.355,75.365,75.555,75.565,75.575,76.554,
    77.355,77.365,81.655,82.656,83.445,83.455,83.655,83.665,83.675,84.456,
    85.255,85.455,85.465,85.475,86.454,91.555,92.556,93.355,93.555,93.565,
    93.575,95.355,95.365,107.755,109.555,115.845,115.855,117.645,117.655,118.654,
    119.455,125.745,125.755,126.556,126.754,127.545,127.555,128.554,129.355,133.855,
    134.656,135.435,135.635,135.645,135.655,135.855,136.555,136.654,137.445,137.455,
    137.655,137.665,138.454,139.455,143.535,143.745,143.755,144.546,144.556,145.535,
    145.545,145.555,145.755,145.765,146.554,147.355,147.555,147.565,148.554,153.645,
    153.655,154.656,155.435,155.445,155.455,155.645,155.655,155.665,155.675,156.555,
    156.654,157.445,157.455,157.465,158.454,161.557,162.556,163.545,163.555,163.755,
    164.554,164.556,165.545,165.555,165.565,165.575,166.554,167.355,167.555,167.565,
    168.554,172.656,173.445,173.645,173.655,173.665,174.456,174.555,175.445,175.455,
    175.465,175.655,175.665,175.675,176.454,182.556,183.545,183.555,183.565,185.355,
    185.365,185.555,185.565,185.575,191.655,193.455,193.465,193.655,193.665,195.255,
    195.455,195.465,195.475,207.855,209.655,215.955,217.755,219.555,225.855,227.645,
    227.655,228.654,229.455,234.756,235.745,235.755,236.556,236.655,236.754,237.545,
    237.555,238.554,239.355,243.635,243.855,244.656,245.435,245.645,245.655,246.456,
    246.555,246.654,247.445,247.455,247.655,248.454,253.535,253.755,254.556,255.535,
    255.545,255.555,255.557,255.755,255.765,256.554,257.355,257.555,257.565,257.575,
    262.656,263.645,263.655,264.555,265.445,265.455,265.655,265.665,265.675,267.455,
    267.465,271.557,272.556,273.545,273.555,273.557,274.554,274.556,275.545,275.555,
    275.565,275.575,276.554,277.555,283.655,283.665,285.455,285.465,285.475,293.555,
    293.565,295.355,295.365,295.555,295.565,295.575,455.555};


    if ((fp_ad = fopen (adm_name,"r")) == NULL)
    {
        printf ("Cannot open admittance file?\n");
        exit (0);
    }
 

    ALL_OTS = (OTSconstit *) calloc ( NALL_OTS, sizeof(OTSconstit));
    for (i = 0; i < NALL_OTS; i ++)
    {
        ALL_OTS[i].ds = doodson[i];

        ALL_OTS[i].argn[0] = (int)(ALL_OTS[i].ds/100)%10;
        ALL_OTS[i].argn[1] = (int)(ALL_OTS[i].ds/10)%10 - 5;
        ALL_OTS[i].argn[2] = (int)(ALL_OTS[i].ds/1)%10 - 5;
        ALL_OTS[i].argn[3] = (int)(ALL_OTS[i].ds*10)%10 - 5;
        ALL_OTS[i].argn[4] = (int)(ALL_OTS[i].ds*100)%10 - 5;
        ALL_OTS[i].argn[5] = (int)(ALL_OTS[i].ds*1000)%10 - 5;
    }
 
   
    OTADM = (double *) calloc ( NALL_OTS * NMAJ_OTS, sizeof(double));
    for (j = 0; j < NMAJ_OTS; j++)
    {
        for (i = 0; i < NALL_OTS; i ++)
        {
            fscanf (fp_ad, "%lf", &OTADM[NALL_OTS * j + i]);
        }
    }

    fclose (fp_ad);    
        
    return 0; 
}





int ts_read_minor(double jdt, double gmst, int nmax, int minor, 
        int nper, OTSperturb *per, double *coef)
{
    double ang, cp, sp, cm, sm, fnm, ds, 
        cosang, sinang, fcos[18], fsin[18];
    int i,j, n,m, l, ind;



    for (i = 0; i < (nmax + 1) * (nmax + 1); i++)
    {
        coef[i] = 0;
    }

    if (minor == 1)
    {
        for (i = 0; i < NALL_OTS; i++)
        {
            arg2theta (jdt, gmst, ALL_OTS[i].argn, &ALL_OTS[i].ang);
            ALL_OTS[i].cosf = cos(ALL_OTS[i].ang);
            ALL_OTS[i].sinf = sin(ALL_OTS[i].ang);

        }

        for (j = 0; j < NMAJ_OTS; j++)
        {
            fcos[j] = 0;
            fsin[j] = 0;
            for (i = 0; i < NALL_OTS; i ++)
            {
                fcos[j] = fcos[j] + ALL_OTS[i].cosf * OTADM[NALL_OTS * j + i];
                fsin[j] = fsin[j] + ALL_OTS[i].sinf * OTADM[NALL_OTS * j + i];
            }
        }
    }

    

    for (i = 0; i < nper; i++)
    {
        if (per[i].n > nmax) 
//        if (per[i].n > nmax || per[i].n < 2) 
        {
            continue;
        }

        n = per[i].n;
        m = per[i].m;
        cp = per[i].cp;
        sp = per[i].sp;
        cm = per[i].cm;
        sm = per[i].sm;
        fnm = per[i].fnm;
        ds = per[i].ds;


        if (minor == 1)
        {
            if (ds == 55.565) j = 0;
            else if (ds == 55.575) j = 1;
            else if (ds == 56.554) j = 2;
            else if (ds == 57.555) j = 3;
            else if (ds == 65.455) j = 4;
            else if (ds == 75.555) j = 5;
            else if (ds == 85.455) j = 6;
            else if (ds == 93.555) j = 7;

            else if (ds == 135.655) j = 8;
            else if (ds == 145.555) j = 9;
            else if (ds == 163.555) j = 10;
            else if (ds == 165.555) j = 11;

            else if (ds == 235.755) j = 12;
            else if (ds == 245.655) j = 13;
            else if (ds == 255.555) j = 14;
            else if (ds == 273.555) j = 15;
            else if (ds == 275.555) j = 16;

            else if (ds == 455.555) j = 17;
            else 
            {
                printf ("error in otidecs: ds = %f not 18 major tides\n", ds);   
                exit (1);
            }
            cosang = fcos[j];
            sinang = fsin[j];
        }
        else if (minor == 0)
        {       
            arg2theta (jdt, gmst, per[i].argn, &ang);
            cosang = cos(ang);
            sinang = sin(ang);
        }
        else 
        {
            printf ("error in otidecs: minor = %d != 0 or 1\n", minor);   
            exit (1);
        }
    
      
        if (m == 0)
        {
            coef[n] = coef[n] + fnm * ((cp+cm) * cosang + (sp+sm) * sinang);
        }
        else
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = coef[ind + n - m] + fnm * ((cp+cm) * cosang + (sp+sm) * sinang);
            coef[ind + n - m + l] = coef[ind + n - m + l] + fnm * ((sp-sm) * cosang - (cp-cm) * sinang);
        }
    }


    return 0;

}









double arg2theta (double jdt, double gmst, int *n, 
        double *ang)
{
    double t, a[5], b[6], theta, l, lp, F, D, Om;
    int i;
      
    t = (jdt - T0) / 36525.0;
    fund_args (t, a);
    l  = a[0];
    lp = a[1];
    F  = a[2];
    D  = a[3];
    Om = a[4];
//    Om = 0;

    b[1] = F + Om;
    b[0] = gmst + TWOPI / 2.0 - b[1];
    b[2] = b[1] - D;
    b[3] = b[1] - l;
//    b[4] = b[2];
    b[4] = - Om;
//    b[4] = 0;
    b[5] = b[1] - D - lp;

    theta = 0;
    for (i = 0; i < 6; i ++)
    {
        theta = theta + b[i] * n[i];
    }


    *ang = theta;

    return 0;
}








