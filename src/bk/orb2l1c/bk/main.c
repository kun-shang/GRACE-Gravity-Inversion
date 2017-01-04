/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE L1B to L1C & L2

  Version: 20 Dec 2011 
           fixed and stochastic constraint solution

  Copyright (c) 2011 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef __OMP_H__
    #include <omp.h>
#endif

#include "l1b2l1c.h"

//#define LFIX 4

short int ACCURACY = 1;
    
short int INTERP = 1;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
main (int argc, char *argv[])
{
    FILE *fp_stdin, *fpout_sim, *fp_otbin, *fp_atbin, *fpout_dbg;  
    int i, i_a, i_b, iter, n, nmax, NOMAX, NAMAX, NPMAX, DIM_OT, DIM_AT,  NDMAX, ncut,mmax, mcut, 
        nply, ncpr, nplycpr, ndt, order_poly, order_cpr, overlapcyc = 0,
        n_gnva, n_gnvb, n_kbr, n_acca, n_accb, n_scaa, simlabel, kbrr,
        n_scab, factor_b1, factor_b2, solvefor, solvegrv, m, l, k, ind, bias, 
        noinv, t, n_ot, LFIX = 0, ifix[4000], LSTK = 0, istk[4000],
        year, month, day, days, dbg;
    short int error, de_num;
    double *fkn, *fint, rpv[3], rpvt[3], rpt[3], rvt[3], wi[3], we[3], tjd[2], 
        mjd0, jd0, GMR2, r2R, vl2, vl1, p1e[3], p2e[3], *coefcsr, 
        *coefjpl, *coefgfz, *coefggm, vrb, vrn, vrs, outlier, fltpar, endcut,
        vf, ef2, ef1, *pt1, *pt2, *ai, *coefb, *coefbl, *coefrfl, *coefrfh, dv12b, dv12n, *COEFA, *COEFD, step_ot, step_at, *COEFO, *COEFP,
        dv12s, df12, acc1[3], acc2[3], agr1a[3] = {0}, agr2b[3] = {0},
        *ptdum, factor_n, *ptdum1, *ptdum2, vgn, ve2, ve1, va1, va2, dv1, dv2, 
        *lambda, *lambdam, *ksiupd, *ksiupdm, *k0, *aksi, *factor_a, *factor_b, factor_p, factor_r,
        *ksi, *ksim, *dksi, *dksiupd, r0[1], sigma0, pweight, *p0, *p0b0, *knk0, 
        *knk, *kksi, *kksim, *id, *kfix, *kfixt,*kn, *nk, *nkknk,
        *atpa, *atpy, *atpym, *cest, *cestm, efbias, efbdot, zero = 0.0, jd_beg, jd_end, 
        omega, ctksi[1], ytpy[1], btpb[1], etpe, e0p0e0, af1[3], af2[3],
        c, s, vfix[4000], vstk[4000];
    char line[500], card[20], f_gnv1c_a[200], f_gnv1c_b[200], f_acc1b_a[200], f_par_a[200], f_par_b[200],
        f_eph[200], f_acc1b_b[200], f_sca1b_a[200], f_sca1b_b[200], f_kbr1b_x[200], f_aot[200],
        f_grv[2][200]={"\0", "\0"}, 
        f_ref[2][200]={"\0", "\0"}, 
        f_csr[2][200]={"\0", "\0"}, 
        f_jpl[2][200]={"\0", "\0"}, 
        f_gfz[2][200]={"\0", "\0"}, 
        f_ot[200],f_pt[200], f_eop[200], stdname[200];    
    time_t s0, s1, s2, s3, s4, s5, s34;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// read input cards
    LFIX = 0;
    if (argc == 1)
    {
        printf ("input file name: ");
        scanf ("%s", stdname);
    }
    else if (argc == 2)
    {
        sscanf (argv[1], "%s", stdname);
    }
    else
    {
        printf ("input argument error!\n");
        exit (0);
    }
    if ( (fp_stdin = fopen (stdname,"r")) == NULL)
    {
        printf ("Cannot open stdin file!\n");
        exit (0);
    }
//    printf ("reading input file \"%s\"...\n", stdname);

    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        if (fgets (line, 500, fp_stdin) == NULL) break;
        sscanf (line, "%s", card);  
        
        if (strcmp (card,"GNV1C") ==0)  
        {
            sscanf (line, "%*s%s%s", f_gnv1c_a, f_gnv1c_b); 
        }
        if (strcmp (card,"PAR1B") ==0)  
        {
            sscanf (line, "%*s%s%s", f_par_a, f_par_b); 
        }
        if (strcmp (card,"ACC1B") ==0)  
        {
            sscanf (line, "%*s%s%s", f_acc1b_a, f_acc1b_b); 
        }
        if (strcmp (card,"SCA1B") ==0)  
        {
            sscanf (line, "%*s%s%s", f_sca1b_a, f_sca1b_b); 
        }
        if (strcmp (card,"KBR1B") ==0)  
        {
            sscanf (line, "%*s%s", f_kbr1b_x);  
        }
        if (strcmp (card,"AOT") ==0)
        {
            sscanf (line, "%*s%s%d", f_aot, &NDMAX);
        }

        if (strcmp (card,"GRV") ==0)    
        {
            sscanf (line, "%*s%s%s", f_grv[0], f_grv[1]);   
        }

        if (strcmp (card,"REF") ==0)    
        {
            sscanf (line, "%*s%s%s", f_ref[0], f_ref[1]);   
        }

        if (strcmp (card,"L2MCSR") ==0)
        {
            sscanf (line, "%*s%s%s", f_csr[0], f_csr[1]);
        }
        if (strcmp (card,"L2MJPL") ==0)
        {
            sscanf (line, "%*s%s%s", f_jpl[0], f_jpl[1]);
        }
        if (strcmp (card,"L2MGFZ") ==0)
        {
            sscanf (line, "%*s%s%s", f_gfz[0], f_gfz[1]);
        }

        if (strcmp (card,"OTIDE") ==0)
        {
            sscanf (line, "%*s%d%d", &NOMAX, &OTIDE);
        }

        if (strcmp (card,"ATIDE") ==0)
        {
            sscanf (line, "%*s%d%d", &NAMAX, &ATIDE);
        }

        if (strcmp (card,"PTIDE") ==0)
        {
            sscanf (line, "%*s%s%d", f_pt, &NPMAX);
        }



        if (strcmp (card,"EPH") ==0)
        {
            sscanf (line, "%*s%s", f_eph);
        }
        if (strcmp (card, "NMAX") ==0)  
        {
            sscanf (line, "%*s%d%d", &nmax, &mmax);
        }
        if (strcmp (card, "NCUT") ==0)  
        {
            sscanf (line, "%*s%d%d", &ncut, &mcut);
        }
        if (strcmp (card, "DT") ==0)    
        {
            sscanf (line, "%*s%d", &DT);
        }
        if (strcmp (card, "YEAR") ==0)
        {
            sscanf (line, "%*s%d", &year);
        }
        if (strcmp (card, "MONTH") ==0)
        {
            sscanf (line, "%*s%d", &month);
        }
        if (strcmp (card, "DAY") ==0)
        {
            sscanf (line, "%*s%d", &day);
        }
        if (strcmp (card, "DAYS") ==0)
        {
            sscanf (line, "%*s%d", &days);
        }
        if (strcmp (card, "KBRR") ==0)
        {
            sscanf (line, "%*s%d", &kbrr);
        }
        
        if (strcmp (card, "FITRES") ==0)
        {
//            sscanf (line, "%*s%d%d%d", &order_poly, &order_cpr, &overlapcyc);
            sscanf (line, "%*s%d%d%d%d", &nply, &ncpr, &nplycpr, &ndt);
            order_poly = nply - 1;
            order_cpr = (int)(ncpr / 2);
        }

        if (strcmp (card, "FILTER") ==0)
        {
            sscanf (line, "%*s%lf%lf", &fltpar, &endcut);
        }
        

        if (strcmp (card, "OUTLIER") ==0)    
        {
            sscanf (line, "%*s%lf", &outlier);
        }
        if (strcmp (card,"EOP") ==0)    
        {
            sscanf (line, "%*s%s", f_eop);  
        }
        if (strcmp (card, "STIDE") ==0)
        {
            sscanf (line, "%*s%d%lf%d", &PERMT, &C20PERM, &STIDE);
        }

        if (strcmp (card, "ACC_BIAS") ==0)
        {
            sscanf (line, "%*s%d%d", &MACC_BIAS, &MACC_DTBS);
        }

        if (strcmp (card, "ACC_SCAL") ==0)
        {
            sscanf (line, "%*s%d%d", &MACC_SCAL, &MACC_DTSL);
        }

        if (strcmp (card, "DEBUG") ==0)
        {
            sscanf (line, "%*s%d", &dbg);
        }


    }


    jd0 = julian_date (year,month,day,0);
    NDATA = days * 86400 / DT;
    GPS_S = (int)((jd0 - T0) * 86400 + 0.5);
    mjd0 = (int)(jd0 - 2400000.5);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open output files 



    if ( (fpout_sim = fopen ("sim.asc","w")) == NULL)
    {
        printf ("Cannot open fpout_sim file!\n");
        exit (0);
    }

    if ( (fpout_dbg = fopen ("dbg.asc","w")) == NULL)
    {
        printf ("Cannot open fpout_dbg file!\n");
        exit (0);
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open ephemeris file

    if ((error = ephem_open (f_eph, &jd_beg,&jd_end,&de_num)) != 0)
    {
      if (error == 1)
         printf ("JPL ephemeris file not found.\n");
       else
         printf ("Error reading JPL ephemeris file header.\n");
      return (error);
    }
    else
    {
//      printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
//         de_num, jd_beg, jd_end);
//      printf ("\n");
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// allocate memory

    solvegrv = (nmax + 1) * (nmax + 1);
    solvefor = (nmax + 1) * (nmax + 1) + bias;

    kbrx = (KBR1B *) calloc ( NDATA, sizeof(KBR1B));
    gnva = (GNV1C *) calloc ( NDATA, sizeof(GNV1C));
    gnvb = (GNV1C *) calloc ( NDATA, sizeof(GNV1C));
    ACA_EPH  = (double *) calloc (86400 * 4 * days, sizeof(double));
    SCA_EPH  = (double *) calloc (17280 * 5 * days, sizeof(double));
    ACB_EPH  = (double *) calloc (86400 * 4 * days, sizeof(double));
    SCB_EPH  = (double *) calloc (17280 * 5 * days, sizeof(double));  
  
    info = (InfStruct *) calloc ( NDATA, sizeof(InfStruct));

    sat1a = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));
    sat2b = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));
    sat12 = (DATAL1C *) calloc ( NDATA, sizeof(DATAL1C));

    fkn = (double *) calloc ( NDATA, sizeof(double));
    fint = (double *) calloc ( NDATA, sizeof(double));


    coefbl  = (double *) calloc (solvegrv, sizeof(double));
    coefb  = (double *) calloc ((ncut + 1) * (ncut + 1), sizeof(double));
    coefrfl  = (double *) calloc (solvegrv, sizeof(double));
    coefrfh  = (double *) calloc ((ncut + 1) * (ncut + 1), sizeof(double));
    coefcsr  = (double *) calloc (solvegrv, sizeof(double));
    coefjpl  = (double *) calloc (solvegrv, sizeof(double));
    coefgfz  = (double *) calloc (solvegrv, sizeof(double));


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// read L1B data

    s0 = time(NULL);  

    for (i = 0; i < NDATA; i ++)
    {
        sat1a[i].t = GPS_S + i * DT;
        sat2b[i].t = GPS_S + i * DT;
        sat12[i].t = GPS_S + i * DT;

        sat1a[i].error = 0;
        sat2b[i].error = 0;
        sat12[i].error = 0;
    }

    readkbr(f_kbr1b_x, &n_kbr);
    readgnv(f_gnv1c_a, f_gnv1c_b, &n_gnva, &n_gnvb);
    readsca(f_sca1b_a, f_sca1b_b);
    readacc(f_acc1b_a, f_acc1b_b);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open ocean tide file and EOP file: OTIDE & INF1B

//    NPMAX = 100;

    OPTM1 = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
    OPTM2 = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
    COEFP = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));
//    openopt (f_pt, NPMAX);




    step_ot = 60;
    DIM_OT = (int)(86400/step_ot);
    OT_EPH  = (double *) calloc ( ( (NOMAX + 1) * (NOMAX + 1) + 1 ) * DIM_OT , sizeof(double));
//    COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));

    if((fp_otbin=fopen(f_aot,"rb"))==NULL)
    {
        printf("Cannot write oteph.bin!\n");
        exit(0);
    }


    fread (OT_EPH, sizeof(double) * ( (NOMAX + 1) * (NOMAX + 1) + 1 ) * DIM_OT, 1, fp_otbin);


    step_at = 60;
    DIM_AT = (int)(86400/step_at);
    AT_EPH  = (double *) calloc ( ( (NAMAX + 1) * (NAMAX + 1) + 1 ) * DIM_AT , sizeof(double));
//    COEFA = (double *) calloc ( (NAMAX + 1) * (NAMAX + 1), sizeof(double));

//    if((fp_atbin=fopen("ateph.bin","rb"))==NULL)
//    {
//        printf("Cannot write ateph.bin!\n");
//        exit(0);
//    }


//    fread (AT_EPH, sizeof(double) * ( (NAMAX + 1) * (NAMAX + 1) + 1 ) * DIM_AT, 1, fp_atbin);












/*        
    if (OTIDE == 1 || OTIDE == 2)    
        openotcs_fes (f_ot);
    if (OTIDE == 3 || OTIDE == 4)
        openotcs_csr (f_ot);

*/



    AOD_EPH  = (double *) calloc ( ( (NDMAX + 1) * (NDMAX + 1) + 1 ) * 5 , sizeof(double));
//        COEFA  = (double *) calloc ( (NAMAX + 1) * (NAMAX + 1) , sizeof(double));
//    openaod (f_aod, NDMAX, NDMAX);


    getinfo(f_eop);
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// open gravity field file

    opengrav (f_grv, coefbl, GMA, nmax, mmax, 0);
    opengrav (f_grv, coefb, GMA, ncut, mcut, 0);
    opengrav (f_ref, coefrfl, GMA, nmax, mmax, 0);
    opengrav (f_ref, coefrfh, GMA, ncut, mcut, 0);
    opengrav (f_csr, coefcsr, GMA, nmax, mmax, 0);
    opengrav (f_jpl, coefjpl, GMA, nmax, mmax, 0);
    opengrav (f_gfz, coefgfz, GMA, nmax, mmax, 0);
    coefgfz[2] = coefgfz[2] + C20PERM;
//    printf ("PERM = %d\n", PERM);
//    if (PERMT == 0)
//    {
//        coefa[2] = coefa[2] + C20PERM;
//        coefe[2] = coefe[2] + C20PERM;
//    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
// calibrate accelerometer data


    cal_acc_01();
//    cal_acc_00();
    calbiaseph(f_par_a, f_par_b);
//    bias_acc(4);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    for (i = 0; i < NDATA; i ++)
    {    
        sat1a[i].rp[0] = gnva[i].xpos; 
        sat1a[i].rp[1] = gnva[i].ypos; 
        sat1a[i].rp[2] = gnva[i].zpos;
        sat1a[i].rv[0] = gnva[i].xvel; 
        sat1a[i].rv[1] = gnva[i].yvel; 
        sat1a[i].rv[2] = gnva[i].zvel;
        sat2b[i].rp[0] = gnvb[i].xpos; 
        sat2b[i].rp[1] = gnvb[i].ypos; 
        sat2b[i].rp[2] = gnvb[i].zpos;
        sat2b[i].rv[0] = gnvb[i].xvel; 
        sat2b[i].rv[1] = gnvb[i].yvel; 
        sat2b[i].rv[2] = gnvb[i].zvel;

        sat1a[i].pos = modvect(sat1a[i].rp);
        sat1a[i].vel = modvect(sat1a[i].rv);
        sat2b[i].pos = modvect(sat2b[i].rp);
        sat2b[i].vel = modvect(sat2b[i].rv);

        for (n = 0; n < 3; n ++)
        {
            sat12[i].rp[n] = sat2b[i].rp[n] - sat1a[i].rp[n];
            sat12[i].rv[n] = sat2b[i].rv[n] - sat1a[i].rv[n];
        }
        sat12[i].pos = modvect(sat12[i].rp);
        sat12[i].vel = modvect(sat12[i].rv);


        sat12[i].range = sat12[i].pos;
        sat12[i].rate = dotvect(sat12[i].rp, sat12[i].rv) / sat12[i].pos;

        if (kbrx[i].gps_time != 0)
        {
            sat12[i].rate = kbrx[i].rate;
        }
        else
        {
            sat12[i].error = 9;
        }

        if (kbrr == 2)
        {
            kbrvel2pos (i);
        }
        else if (kbrr == 1)
        {
            kbrpos2vel (i);
        }

        brmul(info[i].c_ie, sat1a[i].rp, 3, 3, 1, p1e);  
        brmul(info[i].c_ie, sat2b[i].rp, 3, 3, 1, p2e);

        xyz2llh(p1e, sat1a[i].llr);
        xyz2llh(p2e, sat2b[i].llr);

    }

    s1 = time(NULL);  
    printf("\n%5ld: seconds of reading L1B data\n", s1-s0);
    fflush(stdout);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    
    we[0] = 0; we[1] = 0; we[2] = ANGVEL;

    for (i = 0; i < NDATA; i ++)
    {    

        tjd[0] = info[i].jd0;
        tjd[1] = info[i].tt/86400.0;

        crsvect(info[i].wi, sat1a[i].rp, sat1a[i].wr);
        crsvect(info[i].wi, sat2b[i].rp, sat2b[i].wr);

/////////////////
//vk: 1/2*v^2 - GM/r
        sat1a[i].vk = dotvect (sat1a[i].rv, sat1a[i].rv) / 2.0 - GMA[0]/sat1a[i].pos;
        sat2b[i].vk = dotvect (sat2b[i].rv, sat2b[i].rv) / 2.0 - GMA[0]/sat2b[i].pos;
        sat12[i].vk = sat2b[i].vk - sat1a[i].vk;

////////////////
//vr:
        crsvect(sat1a[i].rp, sat1a[i].rv, rpv);
        sat1a[i].vr = dotvect(rpv, info[i].wi);
        sat1a[i].vr3 = dotvect(rpv, we);

        crsvect(sat2b[i].rp, sat2b[i].rv, rpv);
        sat2b[i].vr = dotvect(rpv, info[i].wi);
        sat2b[i].vr3 = dotvect(rpv, we);

        sat12[i].vr = sat2b[i].vr - sat1a[i].vr;
        sat12[i].vr3 = sat2b[i].vr3 - sat1a[i].vr3;

        sat1a[i].vr12 = sat1a[i].vr - sat1a[i].vr3;
        sat2b[i].vr12 = sat2b[i].vr - sat2b[i].vr3;
        sat12[i].vr12 = sat12[i].vr - sat12[i].vr3;

///////////////
//vf:
//        disse_a (sat1a[i].rv, i, &sat1a[i].dvaf, sat1a[i].af);
//        disse_b (sat2b[i].rv, i, &sat2b[i].dvaf, sat2b[i].af);



        accia (tjd, 1, sat1a[i].af);
        accia (tjd, 2, sat2b[i].af);
            
        sat1a[i].dvaf = dotvect(sat1a[i].af, sat1a[i].rv);
        sat2b[i].dvaf = dotvect(sat2b[i].af, sat2b[i].rv);
        sat12[i].dvaf = sat2b[i].dvaf - sat1a[i].dvaf;

        sat1a[i].dvrf =  dotvect(sat1a[i].wr, sat1a[i].af);
        sat2b[i].dvrf =  dotvect(sat2b[i].wr, sat2b[i].af);
        sat12[i].dvrf = sat2b[i].dvrf - sat1a[i].dvrf;

        for (n=0;n<3;n++)
            sat12[i].af[n] = sat2b[i].af[n] - sat1a[i].af[n];


////////////////
//vg:


        accgr (tjd, sat1a[i].rp, sat1a[i].rv, sat1a[i].ag);
        accgr (tjd, sat2b[i].rp, sat2b[i].rv, sat2b[i].ag);


        sat1a[i].dvag = dotvect(sat1a[i].ag, sat1a[i].rv);
        sat2b[i].dvag = dotvect(sat2b[i].ag, sat2b[i].rv);
        sat12[i].dvag = sat2b[i].dvag - sat1a[i].dvag;

        sat1a[i].dvrg =  dotvect(sat1a[i].wr, sat1a[i].ag);
        sat2b[i].dvrg =  dotvect(sat2b[i].wr, sat2b[i].ag);
        sat12[i].dvrg = sat2b[i].dvrg - sat1a[i].dvrg;

        for (n=0;n<3;n++)
            sat12[i].ag[n] = sat2b[i].ag[n] - sat1a[i].ag[n];





///////////////
//vn:
        nbodyv (tjd, sat1a[i].rp, &sat1a[i].vn, &sat1a[i].dvtn, sat1a[i].an);
        nbodyv (tjd, sat2b[i].rp, &sat2b[i].vn, &sat2b[i].dvtn, sat2b[i].an);

        sat12[i].vn = sat2b[i].vn - sat1a[i].vn;
        sat12[i].dvtn = sat2b[i].dvtn - sat1a[i].dvtn;

        sat1a[i].dvan = dotvect(sat1a[i].an, sat1a[i].rv);
        sat2b[i].dvan = dotvect(sat2b[i].an, sat2b[i].rv);
        sat12[i].dvan = sat2b[i].dvan - sat1a[i].dvan;

        sat1a[i].dvrn =  dotvect(sat1a[i].wr, sat1a[i].an);
        sat2b[i].dvrn =  dotvect(sat2b[i].wr, sat2b[i].an);
        sat12[i].dvrn = sat2b[i].dvrn - sat1a[i].dvrn;

        for (n=0;n<3;n++)
            sat12[i].an[n] = sat2b[i].an[n] - sat1a[i].an[n];



//////////////
//vs:
        stidev (i, sat1a[i].llr, &sat1a[i].vs, &sat1a[i].dvts, sat1a[i].as);
        stidev (i, sat2b[i].llr, &sat2b[i].vs, &sat2b[i].dvts, sat2b[i].as);

        sat12[i].vs = sat2b[i].vs - sat1a[i].vs;
        sat12[i].dvts = sat2b[i].dvts - sat1a[i].dvts;

        sat1a[i].dvas = dotvect(sat1a[i].as, sat1a[i].rv);
        sat2b[i].dvas = dotvect(sat2b[i].as, sat2b[i].rv);
        sat12[i].dvas = sat2b[i].dvas - sat1a[i].dvas;

        sat1a[i].dvrs =  dotvect(sat1a[i].wr, sat1a[i].as);
        sat2b[i].dvrs =  dotvect(sat2b[i].wr, sat2b[i].as);
        sat12[i].dvrs = sat2b[i].dvrs - sat1a[i].dvrs;

        for (n=0;n<3;n++)
            sat12[i].as[n] = sat2b[i].as[n] - sat1a[i].as[n];

    }

    s2 = time(NULL);  
    printf("\n%5ld: seconds of processing tides correction\n", s2-s1);
    fflush(stdout);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//*simulation potential difference*/


#pragma omp parallel private(pt1, pt2, COEFA, COEFD, COEFO, COEFP, i, n, m)
    {
        pt1  = (double *) calloc ( solvegrv, sizeof(double));
        pt2  = (double *) calloc ( solvegrv, sizeof(double));
        COEFA = (double *) calloc ( (NAMAX + 1) * (NAMAX + 1), sizeof(double));
        COEFD = (double *) calloc ( (NDMAX + 1) * (NDMAX + 1), sizeof(double));
        COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));
        COEFP = (double *) calloc ( (NPMAX + 1) * (NPMAX + 1), sizeof(double));

        #pragma omp for
    for (i = 0; i < NDATA; i ++)
    {

//////////
//vp:


        poletide (&info[i], NPMAX, COEFP);


        cs2acc (i, sat1a[i].llr, COEFP, GMA[0], GMA[1], NPMAX,
            &sat1a[i].vp, &sat1a[i].dvtp, sat1a[i].ap);
        cs2acc (i, sat2b[i].llr, COEFP, GMA[0], GMA[1], NPMAX,
            &sat2b[i].vp, &sat2b[i].dvtp, sat2b[i].ap);
        sat12[i].vp = sat2b[i].vp - sat1a[i].vp;

        sat1a[i].dvrp =  dotvect(sat1a[i].wr, sat1a[i].ap);
        sat2b[i].dvrp =  dotvect(sat2b[i].wr, sat2b[i].ap);
        sat12[i].dvrp = sat2b[i].dvrp - sat1a[i].dvrp;

        for (n=0;n<3;n++)
            sat12[i].ap[n] = sat2b[i].ap[n] - sat1a[i].ap[n];    // background model n=0~60

        sat1a[i].dvap = dotvect(sat1a[i].ap, sat1a[i].rv);
        sat2b[i].dvap = dotvect(sat2b[i].ap, sat2b[i].rv);
        sat12[i].dvap = sat2b[i].dvap - sat1a[i].dvap;







//////////////
//vo:
//        otidev12 (i, nomax, sat1a[i].llr, sat2b[i].llr);

//        otidecs_csr(&info[i], NOMAX, COEFO);

        lgr_order (OT_EPH, DIM_OT, (NOMAX + 1) * (NOMAX + 1) + 1 , info[i].tt, COEFO, 4);

        cs2acc (i, sat1a[i].llr, COEFO, GMA[0], GMA[1], NOMAX, 
            &sat1a[i].vo, &sat1a[i].dvto, sat1a[i].ao);
        cs2acc (i, sat2b[i].llr, COEFO, GMA[0], GMA[1], NOMAX, 
            &sat2b[i].vo, &sat2b[i].dvto, sat2b[i].ao);
        sat12[i].vo = sat2b[i].vo - sat1a[i].vo;

        sat1a[i].dvro =  dotvect(sat1a[i].wr, sat1a[i].ao);
        sat2b[i].dvro =  dotvect(sat2b[i].wr, sat2b[i].ao);
        sat12[i].dvro = sat2b[i].dvro - sat1a[i].dvro;

        for (n=0;n<3;n++)
            sat12[i].ao[n] = sat2b[i].ao[n] - sat1a[i].ao[n];    // background model n=0~60

        sat1a[i].dvao = dotvect(sat1a[i].ao, sat1a[i].rv);
        sat2b[i].dvao = dotvect(sat2b[i].ao, sat2b[i].rv);
        sat12[i].dvao = sat2b[i].dvao - sat1a[i].dvao;


//////////
//va: air tide
//
        lgr_order (AT_EPH, DIM_AT, (NAMAX + 1) * (NAMAX + 1) + 1 , info[i].tt, COEFA, 4);

        cs2acc (i, sat1a[i].llr, COEFA, GMA[0], GMA[1], NAMAX, 
            &sat1a[i].va, &sat1a[i].dvta, sat1a[i].aa);
        cs2acc (i, sat2b[i].llr, COEFA, GMA[0], GMA[1], NAMAX, 
            &sat2b[i].va, &sat2b[i].dvta, sat2b[i].aa);
        sat12[i].va = sat2b[i].va - sat1a[i].va;

        sat1a[i].dvra =  dotvect(sat1a[i].wr, sat1a[i].aa);
        sat2b[i].dvra =  dotvect(sat2b[i].wr, sat2b[i].aa);
        sat12[i].dvra = sat2b[i].dvra - sat1a[i].dvra;

        for (n=0;n<3;n++)
            sat12[i].aa[n] = sat2b[i].aa[n] - sat1a[i].aa[n];    // background model n=0~60

        sat1a[i].dvaa = dotvect(sat1a[i].aa, sat1a[i].rv);
        sat2b[i].dvaa = dotvect(sat2b[i].aa, sat2b[i].rv);
        sat12[i].dvaa = sat2b[i].dvaa - sat1a[i].dvaa;








//vh:
        lgr_order (AOD_EPH, 5, (NDMAX + 1) * (NDMAX + 1) + 1 , (info[i].gps - GPS_S) / 86400.0, COEFD, 1);

        cs2acc (i, sat1a[i].llr, COEFD, GMA[0], GMA[1], NDMAX, 
            &sat1a[i].vh, &sat1a[i].dvth, sat1a[i].ah);
        cs2acc (i, sat2b[i].llr, COEFD, GMA[0], GMA[1], NDMAX, 
            &sat2b[i].vh, &sat2b[i].dvth, sat2b[i].ah);
        sat12[i].vh = sat2b[i].vh - sat1a[i].vh;

//        sat12[i].dvth = sat2b[i].dvth - sat1a[i].dvth;

        sat1a[i].dvrh =  dotvect(sat1a[i].wr, sat1a[i].ah);
        sat2b[i].dvrh =  dotvect(sat2b[i].wr, sat2b[i].ah);
        sat12[i].dvrh = sat2b[i].dvrh - sat1a[i].dvrh;

        for (n=0;n<3;n++)
            sat12[i].ah[n] = sat2b[i].ah[n] - sat1a[i].ah[n];    // background model n=0~60

        sat1a[i].dvah = dotvect(sat1a[i].ah, sat1a[i].rv);
        sat2b[i].dvah = dotvect(sat2b[i].ah, sat2b[i].rv);
        sat12[i].dvah = sat2b[i].dvah - sat1a[i].dvah;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//reference l1c and los for RL05
//
        cs2acc (i, sat1a[i].llr, coefrfl, GMA[0], GMA[1], nmax, 
            &sat1a[i].vrfl, &sat1a[i].dvtb, sat1a[i].arfl);
        cs2acc (i, sat2b[i].llr, coefrfl, GMA[0], GMA[1], nmax, 
            &sat2b[i].vrfl, &sat2b[i].dvtb, sat2b[i].arfl);
        sat12[i].vrfl = sat2b[i].vrfl - sat1a[i].vrfl;    // background model n=0~60
        for (n=0;n<3;n++)
            sat12[i].arfl[n] = sat2b[i].arfl[n] - sat1a[i].arfl[n];    // background model n=0~60
        sat12[i].grfl = dotvect (sat12[i].arfl, sat12[i].rp) / sat12[i].pos;


        cs2acc (i, sat1a[i].llr, coefrfh, GMA[0], GMA[1], ncut,
            &sat1a[i].vrfh, &sat1a[i].dvtb, sat1a[i].arfh);
        cs2acc (i, sat2b[i].llr, coefrfh, GMA[0], GMA[1], ncut,
            &sat2b[i].vrfh, &sat2b[i].dvtb, sat2b[i].arfh);
        sat12[i].vrfh = sat2b[i].vrfh - sat1a[i].vrfh;    // background model n=0~60
        for (n=0;n<3;n++)
            sat12[i].arfh[n] = sat2b[i].arfh[n] - sat1a[i].arfh[n];    // background model n=0~60
        sat12[i].grfh = dotvect (sat12[i].arfh, sat12[i].rp) / sat12[i].pos;



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


        cs2acc (i, sat1a[i].llr, coefcsr, GMA[0], GMA[1], nmax, 
            &sat1a[i].vcsr, &sat1a[i].dvtb, sat1a[i].acsr);
        cs2acc (i, sat2b[i].llr, coefcsr, GMA[0], GMA[1], nmax, 
            &sat2b[i].vcsr, &sat2b[i].dvtb, sat2b[i].acsr);
        sat12[i].vcsr = sat2b[i].vcsr - sat1a[i].vcsr;  
        for (n=0;n<3;n++)
            sat12[i].acsr[n] = sat2b[i].acsr[n] - sat1a[i].acsr[n];  
        sat12[i].gcsr = dotvect (sat12[i].acsr, sat12[i].rp) / sat12[i].pos;

        cs2acc (i, sat1a[i].llr, coefgfz, GMA[0], GMA[1], nmax, 
            &sat1a[i].vgfz, &sat1a[i].dvtb, sat1a[i].agfz);
        cs2acc (i, sat2b[i].llr, coefgfz, GMA[0], GMA[1], nmax, 
            &sat2b[i].vgfz, &sat2b[i].dvtb, sat2b[i].agfz);
        sat12[i].vgfz = sat2b[i].vgfz - sat1a[i].vgfz;  
        for (n=0;n<3;n++)
            sat12[i].agfz[n] = sat2b[i].agfz[n] - sat1a[i].agfz[n];  
        sat12[i].ggfz = dotvect (sat12[i].agfz, sat12[i].rp) / sat12[i].pos;

        cs2acc (i, sat1a[i].llr, coefjpl, GMA[0], GMA[1], nmax, 
            &sat1a[i].vjpl, &sat1a[i].dvtb, sat1a[i].ajpl);
        cs2acc (i, sat2b[i].llr, coefjpl, GMA[0], GMA[1], nmax, 
            &sat2b[i].vjpl, &sat2b[i].dvtb, sat2b[i].ajpl);
        sat12[i].vjpl = sat2b[i].vjpl - sat1a[i].vjpl;  
        for (n=0;n<3;n++)
            sat12[i].ajpl[n] = sat2b[i].ajpl[n] - sat1a[i].ajpl[n];  
        sat12[i].gjpl = dotvect (sat12[i].ajpl, sat12[i].rp) / sat12[i].pos;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//background model of orbit 
//

        cs2acc (i, sat1a[i].llr, coefbl, GMA[0], GMA[1], nmax,
            &sat1a[i].vbl, &sat1a[i].dvtb, sat1a[i].abl);
        cs2acc (i, sat2b[i].llr, coefbl, GMA[0], GMA[1], nmax,
            &sat2b[i].vbl, &sat2b[i].dvtb, sat2b[i].abl);
        sat12[i].vbl = sat2b[i].vbl - sat1a[i].vbl;    // background model n=0~60
        for (n=0;n<3;n++)
            sat12[i].abl[n] = sat2b[i].abl[n] - sat1a[i].abl[n];    // background model n=0~60
        sat12[i].gbl = dotvect (sat12[i].abl, sat12[i].rp) / sat12[i].pos;


        cs2acc (i, sat1a[i].llr, coefb, GMA[0], GMA[1], ncut, 
            &sat1a[i].vb, &sat1a[i].dvtb, sat1a[i].ab);
        cs2acc (i, sat2b[i].llr, coefb, GMA[0], GMA[1], ncut, 
            &sat2b[i].vb, &sat2b[i].dvtb, sat2b[i].ab);

        sat12[i].dvtb = sat2b[i].dvtb - sat1a[i].dvtb; 

        sat1a[i].dvrb =  dotvect(sat1a[i].wr, sat1a[i].ab);
        sat2b[i].dvrb =  dotvect(sat2b[i].wr, sat2b[i].ab);
        sat12[i].dvrb = sat2b[i].dvrb - sat1a[i].dvrb; 

        for (n=0;n<3;n++)
            sat12[i].ab[n] = sat2b[i].ab[n] - sat1a[i].ab[n];    // background model n=0~60

        sat12[i].vb = sat2b[i].vb - sat1a[i].vb; 
//        sat1a[i].vbh = sat1a[i].vb - sat1a[i].vbl;   
//        sat2b[i].vbh = sat2b[i].vb - sat2b[i].vbl;
//        sat12[i].vbh = sat2b[i].vbh - sat1a[i].vbh;   

        sat1a[i].dvab = dotvect(sat1a[i].ab, sat1a[i].rv);
        sat2b[i].dvab = dotvect(sat2b[i].ab, sat2b[i].rv);
        sat12[i].dvab = sat2b[i].dvab - sat1a[i].dvab;

        
        for (n=0;n<3;n++)
        {
            sat12[i].ra[n] = - GMA[0]/pow(sat2b[i].pos,3) * sat2b[i].rp[n] 
                             + GMA[0]/pow(sat1a[i].pos,3) * sat1a[i].rp[n]
//                           + sat12[i].ab[n] 
//                           + sat12[i].arfh[n] 
                           + sat12[i].an[n] 
                           + sat12[i].as[n]
                           + sat12[i].af[n]
                           + sat12[i].ao[n]
//                           + sat12[i].ah[n]
//                           + sat12[i].aa[n]
//                           + sat12[i].ap[n]
                           + sat12[i].ag[n];
//            data1c[i].a12m[n] = data1c[i].ac[n];
        }
        sat12[i].dvcorr = 
                    - sat12[i].dvaf + sat12[i].dvrf
                    - sat12[i].dvan + sat12[i].dvrn 
                    - sat12[i].dvas + sat12[i].dvrs 
                    - sat12[i].dvao + sat12[i].dvro 
//                    - sat12[i].dvah + sat12[i].dvrh
//                    - sat12[i].dvaa + sat12[i].dvra
//                    - sat12[i].dvap + sat12[i].dvrp
                    - sat12[i].dvag + sat12[i].dvrg;


        sat12[i].gb = dotvect (sat12[i].ab, sat12[i].rp) / sat12[i].pos;

        sat12[i].glos = kbrx[i].accl 
                     + (kbrx[i].rate * kbrx[i].rate - sat12[i].vel * sat12[i].vel) / sat12[i].pos
//                     - dotvect (sat12[i].ra, sat12[i].rp) / sat12[i].pos -  sat12[i].grfh;
                     - dotvect (sat12[i].ra, sat12[i].rp) / sat12[i].pos -  sat12[i].gb;
        
    }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


        free(pt1);
        free(pt2);
        free(COEFA);
        free(COEFD);
        free(COEFO);

    }

    s3 = time(NULL);  
    printf("\n%5ld: seconds of processing earth gravity and partial\n", s3-s2);
    fflush(stdout);


///////////////////
//vint: all the potential needed to be intergrated
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvcorr;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vcorr = fint[i];



    s4 = time(NULL);  
    printf("\n%5ld: seconds of finshing\n", s4-s3);
    fflush(stdout);





    for (i = 0; i < NDATA; i ++)
    {


/*
        sat12[i].vc = sat12[i].vk  - sat12[i].vr
                    - sat12[i].vaf + sat12[i].vrf
                    - sat12[i].van + sat12[i].vrn 
                    - sat12[i].vas + sat12[i].vrs 
                    - sat12[i].vao + sat12[i].vro 
                    - sat12[i].vah + sat12[i].vrh
                    - sat12[i].vag + sat12[i].vrg
                    - sat12[i].vap + sat12[i].vrp
                    - sat12[i].vaa + sat12[i].vra;
*/
        sat12[i].vc = sat12[i].vk  - sat12[i].vr + sat12[i].vcorr;

        sat12[i].vl1c = sat12[i].vc - sat12[i].vb;
//        sat12[i].vl1c = sat12[i].vc - sat12[i].vrfh;


        sat12[i].vl1cm = sat12[i].vk - sat12[i].vrb
                     - sat12[i].vaf 
                     - sat12[i].van 
                     - sat12[i].vas
                     - sat12[i].vao 
                     - sat12[i].vah
                     - sat12[i].vag
                     - sat12[i].vap
                     - sat12[i].vaa;
        sat12[i].vl1cm = sat12[i].vcm - sat12[i].vb;


        sat12[i].vrm = sat12[i].vrb + sat12[i].vrf 
                     + sat12[i].vrn + sat12[i].vrs 
                     + sat12[i].vro + sat12[i].vrh 
                     + sat12[i].vrg + sat12[i].vrp + sat12[i].vra;
        
//        sat12[i].vref = (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl) / 3.0;
//        sat12[i].gref = (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl) / 3.0;
//        sat12[i].vref = (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0;
//        sat12[i].gref = (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl - 3 * sat12[i].gbl) / 3.0;

        sat12[i].vref = 0;
        sat12[i].gref = 0;

//        sat12[i].vref = (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0;
//        sat12[i].gref = (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl - 3 * sat12[i].gbl) / 3.0;
//        sat12[i].vref = sat12[i].vb;
//        sat12[i].gref = sat12[i].gb;

    }





    fflush(stdout);

//    estacc ();
    
    endcut = endcut / DT;

    for (i = 0; i < NDATA; i ++)
//        fkn[i] = sat12[i].vl1c - sat12[i].vref;
        fkn[i] = sat12[i].vl1c;
    daydft (fkn, fint, NDATA,  DT, 5400 / fltpar);
    for (i = endcut; i < NDATA - endcut; i ++)
//        sat12[i].vl1cflt = fint[i] + sat12[i].vref;
        sat12[i].vl1cflt = fint[i];

    for (i = 0; i < NDATA; i ++)
//        fkn[i] = sat12[i].glos - sat12[i].gref;
        fkn[i] = sat12[i].glos;
    daydft (fkn, fint, NDATA,  DT, 5400 / fltpar);
    for (i = endcut; i < NDATA - endcut; i ++)
//        sat12[i].glosflt = fint[i] + sat12[i].gref;
        sat12[i].glosflt = fint[i];



// HARD CODE for fitting!!!!
// HARD CODE for fitting!!!!
//    fitres (order_poly, order_cpr, overlapcyc);

//    if (order_poly != 0 || order_cpr !=0)
    fitl1c_pw (order_poly, order_cpr);
    fitlos_pw (order_poly, order_cpr);

    fitl1c_mv (nply, ncpr, nplycpr, ndt);
    fitlos_mv (nply, ncpr, nplycpr, ndt);

//    printf ("!!!!!!!!!\n");

    for (i = 0; i < NDATA; i ++)
    {

        if (sat12[i].error != 0) continue;
// OUTPUT SIM
// OUTPUT SIM
        sat12[i].vmvfitstd = fabs(sat12[i].vl1c - sat12[i].vmvfit - sat12[i].vref);
//            - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
        sat12[i].vpwfitstd = fabs(sat12[i].vl1c - sat12[i].vpwfit - sat12[i].vref);
//            - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
        sat12[i].gmvfitstd = fabs(sat12[i].glos - sat12[i].gmvfit - sat12[i].gref);
//            - (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl - 3 * sat12[i].gbl) / 3.0);
        sat12[i].gpwfitstd = fabs(sat12[i].glos - sat12[i].gpwfit - sat12[i].gref);
//            - (sat12[i].gcsr + sat12[i].ggfz + sat12[i].gjpl - 3 * sat12[i].gbl) / 3.0);
    
        sat12[i].vrfh -= sat12[i].vb;
        sat12[i].grfh -= sat12[i].gb;


        fprintf (fpout_sim, "%10d ", gnva[i].gps_time);                     //1: gps_time
        fprintf (fpout_sim, "%25.15f ", sat1a[i].llr[0]);                   //2
        fprintf (fpout_sim, "%25.15f ", sat1a[i].llr[1]);                   //3
        fprintf (fpout_sim, "%25.12f ", sat1a[i].llr[2]);                   //4
        fprintf (fpout_sim, "%25.15f ", sat2b[i].llr[0]);                   //5
        fprintf (fpout_sim, "%25.15f ", sat2b[i].llr[1]);                   //6
        fprintf (fpout_sim, "%25.12f ", sat2b[i].llr[2]);                   //7
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh);                     //8
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos - sat12[i].grfh);                     //9
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh - sat12[i].vmvfit);   //10
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh - sat12[i].vpwfit);   //11
        fprintf (fpout_sim, "%20.12e ", sat12[i].vl1cflt - sat12[i].vrfh);                  //12
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos - sat12[i].grfh - sat12[i].gmvfit);   //13
        fprintf (fpout_sim, "%20.12e ", sat12[i].glos - sat12[i].grfh - sat12[i].gpwfit);   //14
        fprintf (fpout_sim, "%20.12e ", sat12[i].glosflt - sat12[i].grfh);                  //15
        fprintf (fpout_sim, "%20.12e ", sat12[i].vcsr - sat12[i].vrfl);      //16
        fprintf (fpout_sim, "%20.12e ", sat12[i].vgfz - sat12[i].vrfl);      //17
        fprintf (fpout_sim, "%20.12e ", sat12[i].vjpl - sat12[i].vrfl);      //18
        fprintf (fpout_sim, "%20.12e ", sat12[i].gcsr - sat12[i].grfl);      //19
        fprintf (fpout_sim, "%20.12e ", sat12[i].ggfz - sat12[i].grfl);      //20
        fprintf (fpout_sim, "%20.12e ", sat12[i].gjpl - sat12[i].grfl);      //21

        fprintf (fpout_sim, "%d\t", sat12[i].error);                        //22
        fprintf (fpout_sim, "%20.12f ", sat12[i].vmvfitres);                //23
        fprintf (fpout_sim, "%20.12f ", sat12[i].vpwfitres);                //24
        fprintf (fpout_sim, "%20.12f ", sat12[i].gmvfitres);                //25
        fprintf (fpout_sim, "%20.12f ", sat12[i].gpwfitres);                //26
        fprintf (fpout_sim, "%20.12f ", sat12[i].vmvfitstd);                //27
        fprintf (fpout_sim, "%20.12f ", sat12[i].vpwfitstd);                //28
        fprintf (fpout_sim, "%20.12f ", sat12[i].gmvfitstd);                //29
        fprintf (fpout_sim, "%20.12f ", sat12[i].gpwfitstd);                //30

        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp2);                      //34
        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp1);                      //33
        fprintf (fpout_sim, "%20.12f ", gnva[i].Tp);                        //31
        fprintf (fpout_sim, "%20.12f ", gnvb[i].Tp);                        //32

/*

        fprintf (fpout_sim, "%20.12f\t", (sat1a[i].llr[1] + sat2b[i].llr[1])/2.0);  //2:
        fprintf (fpout_sim, "%20.12f\t", (sat1a[i].llr[0] + sat2b[i].llr[0])/2.0);  //3:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vcm);                             //4:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vc);                              //5:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vl1cm );                           //6:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitm );                           //7:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vl1c );                            //8:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfit );                            //9:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vb);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vbl);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vbh);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vk);   //14:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vr);   //15:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vaf);   //16:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrf);   //17:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].van);   //18:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrn);   //19:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vas);   //20:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrs);   //21:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vao);   //22:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vro);   //23:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vah);   //24:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrh);   //25:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vag);   //26:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrg);   //27:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vap);   //28:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrp);   //29:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vaa);   //28:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vra);   //29:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vab);   //30:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrb);   //31:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vtb);   //32:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vn);   //33:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vs);   //34:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vo);   //35:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vh);   //36:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vp);  //37:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].va); //38:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vr12 ); //39:

        fprintf (fpout_sim, "%25.15f\t", sat1a[i].llr[0]);  //40:
        fprintf (fpout_sim, "%25.15f\t", sat1a[i].llr[1]);  //41:
        fprintf (fpout_sim, "%25.12f\t", sat1a[i].llr[2]);  //42:
        fprintf (fpout_sim, "%25.15f\t", sat2b[i].llr[0]);  //43:
        fprintf (fpout_sim, "%25.15f\t", sat2b[i].llr[1]);  //44:
        fprintf (fpout_sim, "%25.12f\t", sat2b[i].llr[2]);  //45:
        fprintf (fpout_sim, "%20.12e\t", sat12[i].vl1c - sat12[i].vpwfit);   //46:
        fprintf (fpout_sim, "%20.12e\t", sat12[i].vl1c - sat12[i].vmvfit);   //47:
        fprintf (fpout_sim, "%d\t", sat12[i].error);   //48:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfit);                       //49
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitres);                       //50
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitstd);                       //51
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfit);                       //52
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfitres);  //53:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfitstd);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].glos - sat12[i].gpwfit);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].glos - sat12[i].gmvfit);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfit);                       //49
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfitres);                       //50
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfitstd);                       //51
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfit);                       //52
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfitres);  //53:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfitstd);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gcsr - sat12[i].gbl);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].ggfz - sat12[i].gbl);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gjpl - sat12[i].gbl);  //54:
*/
        fprintf (fpout_sim, "\n");
    }



    s5 = time(NULL);  
    printf("\n%5ld: seconds of process L1C data\n", s5-s4);
    fflush(stdout);


    fflush(fpout_sim);
//        exit(1);





if (dbg != 0)
{

    if ( (fpout_dbg = fopen ("dbg.asc","w")) == NULL)
    {
        printf ("Cannot open fpout_dbg file!\n");
        exit (0);
    }



/////////////////
//vrf:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvrf;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vrf = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvrf;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vrf = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vrf = sat2b[i].vrf - sat1a[i].vrf;


/////////////////
//vaf:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvaf;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vaf = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvaf;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vaf = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vaf = sat2b[i].vaf - sat1a[i].vaf;



///////////////
//vab:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvab;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vab = fint[i];




///////////////////
//vrb:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvrb;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vrb = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvrb;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vrb = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vrb = sat2b[i].vrb - sat1a[i].vrb;


///////////////
//vtb:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvtb;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vtb = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvtb;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vtb = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vtb = sat2b[i].vtb - sat1a[i].vtb;




//////////////////
//vrn:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvrn;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vrn = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvrn;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vrn = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vrn = sat2b[i].vrn - sat1a[i].vrn;



///////////////
//vtn:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvtn;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vtn = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvtn;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vtn = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vtn = sat2b[i].vtn - sat1a[i].vtn;



//////////////////
//van:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvan;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].van = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvan;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].van = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].van = sat2b[i].van - sat1a[i].van;




//////////////////
//vrs:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvrs;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vrs = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvrs;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vrs = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vrs = sat2b[i].vrs - sat1a[i].vrs;



///////////////
//vts:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvts;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vts = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvts;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vts = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vts = sat2b[i].vts - sat1a[i].vts;



///////////////
//vas:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvas;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vas = fint[i];
    
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvas;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vas = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vas = sat2b[i].vas - sat1a[i].vas;




////////////////////
//vro:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvro;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vro = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvro;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vro = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vro = sat2b[i].vro - sat1a[i].vro;


////////////////////
//vto:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvto;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vto = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvto;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vto = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vto = sat2b[i].vto - sat1a[i].vto;


////////////////////
//vao:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvao;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vao = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvao;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vao = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vao = sat2b[i].vao - sat1a[i].vao;


//////////
//vrh
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvrh;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vrh = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvrh;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vrh = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vrh = sat2b[i].vrh - sat1a[i].vrh;


/////////////////////////
//vah
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat1a[i].dvah;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat1a[i].vah = fint[i];

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat2b[i].dvah;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat2b[i].vah = fint[i];

    for (i = 0; i < NDATA; i ++)
        sat12[i].vah = sat2b[i].vah - sat1a[i].vah;


///////////////////////
//vag:

    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvag;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vag = fint[i];

////////////////////////
//vrg:
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvrg;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vrg = fint[i];

/////////////////////
//vap
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvap;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vap = fint[i];

//////////////////
//vrp
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvrp;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vrp = fint[i];


/////////////////////
//vaa
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvaa;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vaa = fint[i];

//////////////////
//vra
    for (i = 0; i < NDATA; i ++)
        fkn[i] = sat12[i].dvra;
    intydt(NDATA, DT, fkn, fint);
    for (i = 0; i < NDATA; i ++)
        sat12[i].vra = fint[i];



    for (i = 0; i < NDATA; i ++)
    {

//        if (sat12[i].error != 0) continue;

        fprintf (fpout_dbg, "%10d ", gnva[i].gps_time);                     //1: gps_time
        fprintf (fpout_dbg, "%25.15f ", sat1a[i].llr[0]);                   //2
        fprintf (fpout_dbg, "%25.15f ", sat1a[i].llr[1]);                   //3
        fprintf (fpout_dbg, "%25.12f ", sat1a[i].llr[2]);                   //4
        fprintf (fpout_dbg, "%25.15f ", sat2b[i].llr[0]);                   //5
        fprintf (fpout_dbg, "%25.15f ", sat2b[i].llr[1]);                   //6
        fprintf (fpout_dbg, "%25.12f ", sat2b[i].llr[2]);                   //7

        fprintf (fpout_dbg, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh);                     //8
        fprintf (fpout_dbg, "%20.12e ", sat12[i].glos - sat12[i].grfh);                     //9
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh - sat12[i].vmvfit);   //10
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vl1c - sat12[i].vrfh - sat12[i].vpwfit);   //11
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vl1cflt - sat12[i].vrfh);                  //12
        fprintf (fpout_dbg, "%20.12e ", sat12[i].glos - sat12[i].grfh - sat12[i].gmvfit);   //13
        fprintf (fpout_dbg, "%20.12e ", sat12[i].glos - sat12[i].grfh - sat12[i].gpwfit);   //14
        fprintf (fpout_dbg, "%20.12e ", sat12[i].glosflt - sat12[i].grfh);                  //15
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vcsr - sat12[i].vrfl);      //16
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vgfz - sat12[i].vrfl);      //17
        fprintf (fpout_dbg, "%20.12e ", sat12[i].vjpl - sat12[i].vrfl);      //18
        fprintf (fpout_dbg, "%20.12e ", sat12[i].gcsr - sat12[i].grfl);      //19
        fprintf (fpout_dbg, "%20.12e ", sat12[i].ggfz - sat12[i].grfl);      //20
        fprintf (fpout_dbg, "%20.12e ", sat12[i].gjpl - sat12[i].grfl);      //21


        fprintf (fpout_dbg, "%d\t", sat12[i].error);                        //22
        fprintf (fpout_dbg, "%20.12f ", sat12[i].vmvfitres);                //23
        fprintf (fpout_dbg, "%20.12f ", sat12[i].vpwfitres);                //24
        fprintf (fpout_dbg, "%20.12f ", sat12[i].gmvfitres);                //25
        fprintf (fpout_dbg, "%20.12f ", sat12[i].gpwfitres);                //26
        fprintf (fpout_dbg, "%20.12f ", sat12[i].vmvfitstd);                //27
        fprintf (fpout_dbg, "%20.12f ", sat12[i].vpwfitstd);                //28
        fprintf (fpout_dbg, "%20.12f ", sat12[i].gmvfitstd);                //29
        fprintf (fpout_dbg, "%20.12f ", sat12[i].gpwfitstd);                //30

        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp2);                      //34
        fprintf (fpout_sim, "%20.12f ", sat12[i].Tp1);                      //33
        fprintf (fpout_dbg, "%20.12f ", gnva[i].Tp);                        //31
        fprintf (fpout_dbg, "%20.12f ", gnvb[i].Tp);                        //32

/*

        fprintf (fpout_sim, "%20.12f\t", (sat1a[i].llr[1] + sat2b[i].llr[1])/2.0);  //2:
        fprintf (fpout_sim, "%20.12f\t", (sat1a[i].llr[0] + sat2b[i].llr[0])/2.0);  //3:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vcm);                             //4:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vc);                              //5:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vl1cm );                           //6:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitm );                           //7:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vl1c );                            //8:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfit );                            //9:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vb);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vbl);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vbh);   //13:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vk);   //14:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vr);   //15:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vaf);   //16:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrf);   //17:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].van);   //18:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrn);   //19:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vas);   //20:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrs);   //21:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vao);   //22:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vro);   //23:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vah);   //24:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrh);   //25:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vag);   //26:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrg);   //27:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vap);   //28:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrp);   //29:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vaa);   //28:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vra);   //29:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vab);   //30:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vrb);   //31:

        fprintf (fpout_sim, "%20.12f\t", sat12[i].vtb);   //32:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vn);   //33:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vs);   //34:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vo);   //35:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vh);   //36:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vp);  //37:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].va); //38:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vr12 ); //39:

        fprintf (fpout_sim, "%25.15f\t", sat1a[i].llr[0]);  //40:
        fprintf (fpout_sim, "%25.15f\t", sat1a[i].llr[1]);  //41:
        fprintf (fpout_sim, "%25.12f\t", sat1a[i].llr[2]);  //42:
        fprintf (fpout_sim, "%25.15f\t", sat2b[i].llr[0]);  //43:
        fprintf (fpout_sim, "%25.15f\t", sat2b[i].llr[1]);  //44:
        fprintf (fpout_sim, "%25.12f\t", sat2b[i].llr[2]);  //45:
        fprintf (fpout_sim, "%20.12e\t", sat12[i].vl1c - sat12[i].vpwfit);   //46:
        fprintf (fpout_sim, "%20.12e\t", sat12[i].vl1c - sat12[i].vmvfit);   //47:
        fprintf (fpout_sim, "%d\t", sat12[i].error);   //48:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfit);                       //49
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitres);                       //50
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vpwfitstd);                       //51
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfit);                       //52
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfitres);  //53:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].vmvfitstd);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].glos - sat12[i].gpwfit);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].glos - sat12[i].gmvfit);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfit);                       //49
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfitres);                       //50
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gpwfitstd);                       //51
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfit);                       //52
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfitres);  //53:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gmvfitstd);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gcsr - sat12[i].gbl);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].ggfz - sat12[i].gbl);  //54:
        fprintf (fpout_sim, "%20.12f\t", sat12[i].gjpl - sat12[i].gbl);  //54:
*/
        fprintf (fpout_dbg, "\n");
    }


    fflush(fpout_dbg);

    fclose(fpout_dbg);
}

    free (coefb);
    free (coefbl);
    free (coefrfl);
    free (coefrfh);
    free (coefjpl);
    free (coefcsr);
    free (coefgfz);
    free (kbrx);
    free (gnva);
    free (gnvb);
    free (ACA_EPH);
    free (SCA_EPH);
    free (ACB_EPH);
    free (SCB_EPH);
    free (AOD_EPH);
    free (info);
    free (sat1a);
    free (sat2b);
    free (sat12);
    free (fkn);
    free (fint);
    fclose(fpout_sim);
   
    ephem_close();  /* remove this line for use with solsys version 2 */

    printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");

    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 


