/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

  GRACE L1B to L1C & L2

  Version: 20 Dec 2011 
           fixed and stochastic constraint solution

  Copyright (c) 2011 Kun Shang (shang.34@osu.edu) All Right Reserved 

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _L1B2L1C_H_
    #include "l1b2l1c.h"
#endif


#define MAXLINE 500



double estacc()
{
    int ndt, bias, scale, narc, npar, nmax, a, p, i, n, m, l, nobs, ns, is;
    double *amax, *ai, *pi, *x, *y, *yfit, *atpa, *atpy;

    ndt = 5400;
    bias = 3;
    scale = 0;
    narc = (int)(86400.0/ndt);
    npar = (bias + scale) * 2;
    nmax = ndt / DT;
/*
    amax = (double *) calloc ( npar * nmax, sizeof(double));
    ai = (double *) calloc ( npar * nmax, sizeof(double));
    pi = (double *) calloc ( npar, sizeof(double));    
    x = (double *) calloc ( npar, sizeof(double));    
    y = (double *) calloc ( nmax, sizeof(double));
    yfit = (double *) calloc ( nmax, sizeof(double));
    atpy = (double *) calloc ( npar, sizeof(double));
    atpa = (double *) calloc ( npar * npar, sizeof(double));
*/

    for (a = 0; a < narc; a++)
    {
    
        amax = (double *) calloc ( npar * nmax, sizeof(double));
        ai = (double *) calloc ( npar * nmax, sizeof(double));
        pi = (double *) calloc ( npar, sizeof(double));    
        x = (double *) calloc ( npar, sizeof(double));    
        y = (double *) calloc ( nmax, sizeof(double));
        yfit = (double *) calloc ( nmax, sizeof(double));
        atpy = (double *) calloc ( npar, sizeof(double));
        atpa = (double *) calloc ( npar * npar, sizeof(double));

        for (p = 0; p < npar; p ++)
        {
            pi[p] = 0;
        }
        nobs = 0;
        for (i = 0; i < nmax; i ++)
        {
            l = a * nmax + i;
            pi[0] = pi[0] + sat1a[l].rv[0]; 
            pi[1] = pi[1] + sat1a[l].rv[1]; 
            pi[2] = pi[2] + sat1a[l].rv[2]; 
            pi[3] = pi[3] - sat2b[l].rv[0]; 
            pi[4] = pi[4] - sat2b[l].rv[1]; 
            pi[5] = pi[5] - sat2b[l].rv[2]; 

            if (scale != 0)
            {
                pi[6]  = pi[6]  + sat1a[l].rv[0] * sat1a[l].af[0]; 
                pi[7]  = pi[7]  + sat1a[l].rv[1] * sat1a[l].af[1]; 
                pi[8]  = pi[8]  + sat1a[l].rv[2] * sat1a[l].af[2]; 
                pi[9]  = pi[9]  - sat2b[l].rv[0] * sat2b[l].af[0]; 
                pi[10] = pi[10] - sat2b[l].rv[1] * sat2b[l].af[1]; 
                pi[11] = pi[11] - sat2b[l].rv[2] * sat2b[l].af[2]; 
            }
            
            for (p = 0; p < npar; p ++)
                amax [i * npar + p] = pi[p];
            if (sat12[l].error == 0)
            {
                for (p = 0; p < npar; p ++)
                {
                    ai[nobs * npar + p ] = amax [i * npar + p];
                }
                y[nobs] = sat12[l].vl1c;
                nobs ++;
            }
        }
            

        for (n = 0; n < npar; n++)
        {
            ns = n * npar;
            for (i = 0; i < nobs; i ++)
            {
                is = i * npar;
                for (m = 0; m < npar; m++)
                {
                    atpa[ns + m] +=  ai[is + n] * ai[is + m];
                }
                atpy[n] +=  ai[is + n] * y[i];
            }
        }

        brinv(atpa, npar);
        brmul(atpa, atpy, npar, npar, 1, x);

//        solvels_chol(atpa, npar, atpy, x, 0);
//        solvegaus(atpa, solvefor, atpy, ksi);


        brmul (amax, x, nmax, npar, 1, yfit);
        for (i = 0; i < nmax; i ++)
        {
            l = a * nmax + i;
            sat12[l].fitacc = yfit[i];
        }


        free (amax);
        free (ai);
        free (pi);
        free (x);
        free (y);
        free (yfit);
        free (atpy);
        free (atpa);
    }
    return 0;

}




double openopt (char *f_pt, int nmax)
{
    FILE *fp_grv;
    double cm1, sm1, cm2, sm2;
    int n, m, ic, is, l, ind;
    char string[200], name[20];

    if ((fp_grv = fopen (f_pt,"r")) == NULL)
    {
        printf ("Cannot open pole tides file?\n");
        exit (0);
    }



    while (1)
    {
        if (fgets (string, 200, fp_grv) == NULL) break;
        sscanf (string, "%d%d%lf%lf%lf%lf", &n, &m, &cm1, &sm1, &cm2, &sm2);
        if (n > nmax || n < 2 )
        {
            continue;
        }
        else if (m == 0)
        {
            OPTM1[n] = cm1;
            OPTM2[n] = cm2;
        }
        else
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            OPTM1[ic] = cm1;
            OPTM1[is] = sm1;
            OPTM2[ic] = cm2;
            OPTM2[is] = sm2;
        }
    }

    fclose(fp_grv);
    return 0;

}


double poletide (InfStruct *info, int nmax, double *coef)
{
    double t, m1, m2;
    int n, m, k, l, ind, ic, is;

    mtpole (info->mjd, info->xp, info->yp, &m1, &m2);


    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                coef[n] = OPTM1[n] * m1 + OPTM2[n] * m2;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
                coef[ic] = OPTM1[ic] * m1 + OPTM2[ic] * m2;
                coef[is] = OPTM1[is] * m1 + OPTM2[is] * m2;
            }
        }
    }

    n = 2; m = 1;
    l = nmax - m + 1;
    ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);

    ic = ind + n - m;
    is = ind + n - m + l;

    coef[ic] += -1.333e-9 * (m1 + 0.0115 * m2);
    coef[is] += -1.333e-9 * (m2 - 0.0115 * m1);

//    coef[ic] += -2.1778e-10   * (m1 -  0.01724 * m2);
//    coef[is] += -1.7232e-10   * (m2 -  0.03365 * m1);


    return 0;
}







double intydt(int num, double dt, double *y, double *val)
{
    int i;
    double cum;

    cum = 0;

    val[0] = 0;
    for (i = 1; i < num; i ++)
    {
        cum = cum + (y[i] + y[i-1]) * dt / 2.0;
        val[i] = cum;
    }

    return;
}        

    


























double reltivdiff (double vb1, double vb2, double *vrel)
{
    *vrel = vb1 * vb1 / C / C - vb2 * vb2 / C / C;
//    *vrel = vb1 * vb1  - vb2 * vb2;
    return 0;
}


/*

void reltivint ()
{
    int i, n;
    
    double x1[6], x2[6], v1[3], v2[3], a1[3], a2[3], f1[3], f2[3],
            vf, frelv1, frelv2, tjd[2], df12;


    vf = 0;
    for (i = 0; i < NDATA; i ++)
    {
        if (gnva[i].r == 0) continue;

        tjd[0] = info[i].jd0;
        tjd[1] = info[i].tt/86400.0;

        for (n = 0; n < 3; n++)
        {
            x1[n] = data1c[i].p1m[n] / AU;
            x1[n + 3] = data1c[i].v1m[n] / AU * 86400.0; 
            v1[n] = data1c[i].v1m[n];
            x2[n] = data1c[i].p2m[n] / AU;
            x2[n + 3] = data1c[i].v2m[n] / AU * 86400.0;
            v2[n] = data1c[i].v2m[n];
        }
        accel_pmiers (tjd, x1, a1, f1);
        accel_pmiers (tjd, x2, a2, f2);

        frelv1 = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]) * AU/86400.0/86400.0;
        frelv2 = (f2[0] * v2[0] + f2[1] * v2[1] + f2[2] * v2[2]) * AU/86400.0/86400.0;
//vf = integral (dvfdt) dt
            
        df12 = frelv2 - frelv1;

//        printf ("%15.12f\t%15.12f\t%15.12f\n", frelv1, frelv2, df12);
//        printf ("%e\t%e\t%e\n", frelv1, frelv2, df12);

        vf = vf + df12 * DT;
        data1c[i].vrel = vf;

    }


}


*/





double accgr (double *tjd, double *p, double *v, double *fgr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b,
        pxv[3], vxJ[3], ps[3], vs[3], xsc[6], rs, vsxps[3], vsxpsxv[3],
        unit[9], ppt[9], r5, r3, term1[3], term2[3], term3[3];
    int n;
    short int sun = 10;

    GME = GMA[0]; //m^3/s^2
    r = modvect(p);

    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2

    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;

    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);

    planet_ephemeris (tjd, 2, 10, ps, vs);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        vs[n] = vs[n] * AU / 86400.0;
    }

    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C / C / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C / C / r / r / r * (1 + gamma) *
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C / C / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
    }


    return 0;




}



















double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b, p[3], v[3],
        pxv[3], vxJ[3], ps[3], vs[3], rs, vsxps[3], vsxpsxv[3],
        term1[3], term2[3], term3[3];
    int n;

    GME = 398600.44150E+09; //m^3/s^2
    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2

    GME = GME * 86400.0 * 86400.0 / AU / AU / AU;
    GMS = GMS * 86400.0 * 86400.0 / AU / AU / AU;
    J = J * 86400.0 / AU / AU;

    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;
    p[0] = x[0]; p[1] = x[1]; p[2] = x[2];
    v[0] = x[3]; v[1] = x[4]; v[2] = x[5];

    r = modvect(p);
    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);

    planet_ephemeris (tjd, 2, 10, ps, vs);
    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C_AUDAY / C_AUDAY / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C_AUDAY / C_AUDAY / r / r / r * (1 + gamma) *
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C_AUDAY / C_AUDAY / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
//        printf ("%15.12f\t%15.12f\t%15.12f\n", term1[n],term2[n],term2[n]);
    }


    fnt[0] = - GME / (r*r*r) * p[0];
    fnt[1] = - GME / (r*r*r) * p[1];
    fnt[2] = - GME / (r*r*r) * p[2];

    return 0;
}




/*                               -*- Mode: C -*- 
 * Filename: dft.c
 * Copyright (C) Dan Noland 2003
 * Author: Dan Noland
 * Created: Thu Mar 11 03:19:19 2004
 *           By: Dan Noland
 * Last-Updated: Thu Mar 11 03:51:54 2004
 *     Update #: 12
 * Status: 
 */


#define TRUE 1
#define FALSE 0

// Take, eat, this is my code which is hacked for you...
// (In case you didn't grasp that one consider this public domain)
// nolandda Thu Mar 11 03:49:18 EST 2004

/*
  Discrete Fourier Transform
*/

int dft(long int length, double real_sample[], double imag_sample[])
{
  long int i, j;
  double arg;
  double cosarg,sinarg;
  double *temp_real=NULL,*temp_imag=NULL;

  temp_real = calloc(length, sizeof(double));
  temp_imag = calloc(length, sizeof(double));
  if (temp_real == NULL || temp_imag == NULL)
  {
    return(FALSE);
  }

  for(i=0; i<length; i+=1) 
  {
    temp_real[i] = 0;
    temp_imag[i] = 0;
    arg = -1.0 * TWOPI * (double)i / (double)length;
    for(j=0; j<length; j+=1) 
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
    }
  }

  /* Copy the data back */
  for (i=0; i<length; i+=1) 
  {
    real_sample[i] = temp_real[i];
    imag_sample[i] = temp_imag[i];
  }

  free(temp_real);
  free(temp_imag);
  return(TRUE);
}


/*
  Inverse Discrete Fourier Transform
*/

int inverse_dft(long int length, double real_sample[], double imag_sample[])
{
  long int i, j;
  double arg;
  double cosarg,sinarg;
  double *temp_real=NULL,*temp_imag=NULL;

  temp_real = calloc(length, sizeof(double));
  temp_imag = calloc(length, sizeof(double));
  if (temp_real == NULL || temp_imag == NULL)
  {
    return(FALSE);
  }

  for(i=0; i<length; i+=1) 
  {
    temp_real[i] = 0;
    temp_imag[i] = 0;
    arg = TWOPI * (double)i / (double)length;
    for(j=0; j<length; j+=1) 
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
    }
  }

  /* Copy the data back */
  for (i=0; i<length; i+=1) 
  {
    real_sample[i] = temp_real[i] / (double)length;
    imag_sample[i] = temp_imag[i] / (double)length;
  }

  free(temp_real);
  free(temp_imag);
  return(TRUE);
}



double daydft (double *obs, double *fft, int num,  double dt, double tc)
{
    double *xr, *xi, *yr, *yi, *win, fc, df, kc, th, fh, kh;
    int i, n, m, k;

    n = 1; k = 0;

    while(n<num)
    {
        n = n * 2;
        k = k + 1;
    }

//    n = num;

    xr = (double *) calloc (n, sizeof(double));
    xi = (double *) calloc (n, sizeof(double));
    win = (double *) calloc (n, sizeof(double));
    yr = (double *) calloc (n, sizeof(double));
    yi = (double *) calloc (n, sizeof(double));


    for (i = 0; i < n; i++)
    {
        if (i>=num)
        {
            xr[i] = 0;
            xi[i] = 0;
        }
        else
        {
            xr[i] = obs[i];
            xi[i] = 0;
        }
    }



//    dft(n, xr, xi);

    kkfft(xr,xi,n,k,yr,yi,0,1); 


/*       
    for (i = 0; i < num; i++)
    {
        if (i < 40 || (num - i) < 40)
        {
            xr[i] = 0;
            xi[i] = 0;
        }
    }
*/

    df = 1.0 / n / dt;

     
    fc = 1.0 / tc;
    kc = fc / df;

    th = 60; // seconds
    fh = 1.0 / th;
    kh = fh / df;

//    printf("fc = %f\t df = %f\t dt = %f\t kc = %f\t kh = %f\n", fc, df, dt, kc, kh);


//    kc = 120;

    m = 2;

    getwinhp(n, win, kc, m);

//    getwinbp(n, win, kc, kh, m);


    for (i = 0; i < n; i++)
    {
//        xr[i] = xr[i] * win[i];        xi[i] = xi[i] * win[i];
        yr[i] = yr[i] * win[i];        yi[i] = yi[i] * win[i];
    }

//    inverse_dft(n, xr, xi);

  

    kkfft(yr,yi,n,k,xr,xi,1,1); 


    for (i = 0; i < num; i++)
    {           
        fft[i] = xr[i];
//        fft[i] = sqrt(xr[i] * xr[i] + xi[i] * xi[i]);

    }

    free(xr);
    free(yr);
    free(xi);
    free(yi);
    free(win);


    return 0;
}



double getwinhp(int n, double *win, double kc, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n - i)/kc), 2*m));
        }

    }

    return 0;
}


double getwinbp(int n, double *win, double kc, double kh, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n - i)/kc), 2*m));
        }

    }

    for (i = 0; i < n; i++)
    {
        if (i < n / 2 )
        {
            win[i] = win[i] * 1/sqrt(1 + pow((i/kh), 2*m));
        }
        else if (i >= n / 2 )
        {
            win[i] = win[i] * 1/sqrt(1 + pow(((n - i)/kh), 2*m));
        }

    }


    return 0;
}

double getwinbp0(int n, double *win, double kc, int m)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i < n / 4 )
        {
            win[i] = 1 - 1/sqrt(1 + pow((i/kc), 2*m));
        }
        else if (i >= n / 4 && i < n / 2)
        {
            win[i] = 1 - 1/sqrt(1 + pow(((n/2 - i)/kc), 2*m));
        }
        else if (i >= n / 2 && i < 3 * n / 4)
        {
            win[i] = 1 - 1/sqrt(1 + pow((( i - n/2 )/kc), 2*m));
        }
        else if (i >= 3 * n / 4)
        {
            win[i] = 1 - 1/sqrt(1 + pow((( n - i)/kc), 2*m));
        }

    }

    return 0;
}



double dayfft (double *data, double *fft, int num,  double fs)
/*
    ÓÃÍ¾£º  ¶Ôk5Êý¾Ý×öFFT£¬Êä³öFFT×î´óÊ±µÄÆµÂÊÖµ
    ÊäÈë£º  double data[]           Ô­Ê¼Êý¾Ý
            int res                 ·Ö±æÂÊ£¨Hz£©½øÐÐFFTµÄÊý¾Ý³¤¶È£¨¿ÉÒÔ´óÓÚÔ­Ê¼Êý¾Ý³¤¶È£¬²¹Áã£©
            int debug               µ÷ÊÔ¿ª¹Ø£¨=1Êä³öFFTµÄ½á¹ûÎÄ¼þ£¬=0²»Êä³ö£©£¬
            double fln              Êä³öÏÔÊ¾µÄµÍÆµ
            double fhn              Êä³öÏÔÊ¾µÄµÍÆµ
            double fs               ²ÉÑùÂÊ
            int loop                debugÄ£Ê½ÏÂÊä³öµÄÎÄ¼þÃû
    Êä³ö£º  double *f0              FFT×î´óÊ±µÄÆµÂÊÖµ
    ·µ»Ø£º  =0                      Õý³£
*/

{
    int i,  label = 0; 
    double *t, *xr, *xi, *yr, *yi, fn, fln, fhn, df,dfn, f;
    int n = 1, k = 0;   //n ×ÊÁÏ³¤¶È
    double fl;
    double fh;
        
    while(n<num)
    {
        n = n * 2;
        k = k + 1;
    }
//  n = n / 2;
//  k = k - 1;

    xr = (double *) calloc (n, sizeof(double));
    xi = (double *) calloc (n, sizeof(double));

    

    for (i = 0; i < n; i++)
    {
        if (i>=num)
        {
            label = 1;
            xr[i] = 0;
            xi[i] = 0;
        }
        else
        {
            xr[i] = data[i];
            xi[i] = 0;
        }
    }

//  µÍÆµ£¬¸ßÆµ£¬
    fln = fl/fs;
    fhn = fh/fs;

    df = fs/n;
    dfn = df/fs;
    fn = (fhn - fln) / dfn;


    yr = (double *) calloc (n, sizeof(double));
    yi = (double *) calloc (n, sizeof(double));



    kkfft(xr,xi,n,k,yr,yi,0,1);     //xr,xiÊäÈëµÄÊµ²¿ºÍÐé²¿£¬·µ»ØÄ£ºÍ·ù½Ç£¨¶È£©, yr,yi·µ»Ø±ä»»ºóµÄÊµ²¿ºÍÐé²¿, 
                                    //n¸öµã, n=2^k, 0±íÊ¾¸µÁ¢Ò¶±ä»»£¨1ÎªÄæ±ä»»£©, 1±íÊ¾¼ÆËãÄ£ºÍ·ù½Ç£¨¶È£©
/*

    max = 0;
    for (f = fl; f < fh; f = f + df)
    {
        i = (int)(f/df);
        
        if (xr[i]>max)
        {
            max = xr[i];
            *f0 = f;
        }
    if (debug == 1)
        fprintf (fpout, "%e\t%e\n", f, xr[i]/sqrt(n));
    }
*/


    for (i = 0; i < n; i++)
    {
        if (i < 600 || (n - i) < 600)
        {
            yr[i] = 0;
            yi[i] = 0;
        }
    }

    kkfft(yr,yi,n,k,xr,xi,1,1);     //xr,xiÊäÈëµÄÊµ²¿ºÍÐé²¿£¬·µ»ØÄ£ºÍ·ù½Ç£¨¶È£©, yr,yi·µ»Ø±ä»»ºóµÄÊµ²¿ºÍÐé²¿, 
                                    //n¸öµã, n=2^k, 0±íÊ¾¸µÁ¢Ò¶±ä»»£¨1ÎªÄæ±ä»»£©, 1±íÊ¾¼ÆËãÄ£ºÍ·ù½Ç£¨¶È£©


    for (i = 0; i < num; i++)
    {           
        fft[i] = xr[i];
    }


    free(xr);
    free(yr);
    free(xi);
    free(yi);
    

    return 0;
}




void kkfft(double pr[], double pi[], int n, int k, double fr[], 
           double fi[], int l, int il)
  { int it,m,is,i,j,nv,l0;
    double p,q,s,vr,vi,poddr,poddi;
    for (it=0; it<=n-1; it++)
      { m=it; is=0;
        for (i=0; i<=k-1; i++)
          { j=m/2; is=2*is+(m-2*j); m=j;}
        fr[it]=pr[is]; fi[it]=pi[is];
      }
    pr[0]=1.0; pi[0]=0.0;
    p=6.283185306/(1.0*n);
    pr[1]=cos(p); pi[1]=-sin(p);
    if (l!=0) pi[1]=-pi[1];
    for (i=2; i<=n-1; i++)
      { p=pr[i-1]*pr[1]; q=pi[i-1]*pi[1];
        s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
        pr[i]=p-q; pi[i]=s-p-q;
      }
    for (it=0; it<=n-2; it=it+2)
      { vr=fr[it]; vi=fi[it];
        fr[it]=vr+fr[it+1]; fi[it]=vi+fi[it+1];
        fr[it+1]=vr-fr[it+1]; fi[it+1]=vi-fi[it+1];
      }
    m=n/2; nv=2;
    for (l0=k-2; l0>=0; l0--)
      { m=m/2; nv=2*nv;
        for (it=0; it<=(m-1)*nv; it=it+nv)
          for (j=0; j<=(nv/2)-1; j++)
            { p=pr[m*j]*fr[it+j+nv/2];
              q=pi[m*j]*fi[it+j+nv/2];
              s=pr[m*j]+pi[m*j];
              s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
              poddr=p-q; poddi=s-p-q;
              fr[it+j+nv/2]=fr[it+j]-poddr;
              fi[it+j+nv/2]=fi[it+j]-poddi;
              fr[it+j]=fr[it+j]+poddr;
              fi[it+j]=fi[it+j]+poddi;
            }
      }
    if (l!=0)
      for (i=0; i<=n-1; i++)
        { fr[i]=fr[i]/(1.0*n);
          fi[i]=fi[i]/(1.0*n);
        }
    if (il!=0)
      for (i=0; i<=n-1; i++)
        { pr[i]=sqrt(fr[i]*fr[i]+fi[i]*fi[i]);
          if (fabs(fr[i])<0.000001*fabs(fi[i]))
            { if ((fi[i]*fr[i])>0) pi[i]=90.0;
              else pi[i]=-90.0;
            }
          else
            pi[i]=atan(fi[i]/fr[i])*360.0/6.283185306;
        }
    return;
  }






double resfft(double *res, double *resflt, int num)
{
    int n, i, isign;
    double *data;

    n = 1;
    while(n<num)
    {
        n = n * 2;
    }

    data = (double *) malloc (n * sizeof(double));

    for (i = 0; i < n; i++)
    {
        if (i>num)
        {
            data[i] = 0;
        }
        else
        {
            data[i] = res[i];
        }
    }
//  if (label == 1)
//      printf("²¹Áã!!!\n");

    isign = 1;

    realft(data,n,isign);

    for (i = 0; i < num; i++)
    {           
        resflt[i] = data[i];
    }



}

void realft(double data[], int n, int isign)
{
    int i,i1,i2,i3,i4,n2p3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=3.141592653589793/(double) n;
    if (isign == 1) {
        c2 = -0.5;
        four1(data,n,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    n2p3=2*n+3;
    for (i=2;i<=n/2;i++) {
        i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n,-1);
    }
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], int nn, int isign)
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
        theta=6.28318530717959/(isign*mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

#undef SWAP


/***************************************************************************************/
/*                                                                                     */
/*      Functions for k5vssp(32) filter                                                */
/*                                                                                     */
/*      Version:    2008-7-3                                                          */
/*                                                                                     */
/*      Copyright (c) 2008 shangkun@shao.ac.cn All Right Reserved                      */
/*                                                                                     */
/***************************************************************************************/

/*

  Version:    2008-6-17 
  Version:    2008-6-18     ÐÞ¸Äy[]ÀÛ¼ÓÎ´ÇåÁãµÄÖÂÃü´íÎó 
  Version:    2008-7-3      Ôö¼ÓÊäÈë²ÎÊýorderÂË²¨½×ÊýÐèÍâ²¿Ö¸¶¨ 

*/

double k5filter (double x[], double y[], int size, double fln, 
                 double fhn, double fs, int order)
/*
    ÓÃÍ¾£º  ¶Ôk5Ò»ÃëÊý¾ÝÂË²¨
    ÊäÈë£º  double x[]              Ô­Ê¼Ò»ÃëÊý¾Ý
            int size                Êý¾Ý³¤¶È
            double fln              µÍÆµ
            double fhn              ¸ßÆµ
            double fs                   ²ÉÑùÂÊ£¨Èç¹ûÊÇÒ»ÃëÊý¾Ý=size£©
    Êä³ö£º  double y[]              ÂË²¨ºóµÄÊý¾Ý£¬Ã»ÓÐ³Ë·Å´óÒò×Ó
    ·µ»Ø£º  =0                      Õý³£
*/
{
    double h[1000] = {0}, fl, fh;
    int i=0;
    int n,k;

//  h = (double *) calloc ( order , sizeof(double));
    fl = fln/fs;
    fh = fhn/fs;

    firwin(order,fl,fh,h);
    for(n = 0 ; n < order/2 ; n++)
    {
        y[n] = 0;
        for(k = 0;k<order/2+n;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    for(n = order/2 ; n < size - order/2 ; n++)
    {
        y[n] = 0;
        for(k=0;k<order;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    for(n = size - order/2 ; n < size  ; n++)
    {
        y[n] = 0;
        for(k=order/2-size+n;k<order;k++)
        {
            y[n] = y[n] + x[n-k+order/2] * h[k];
        }
    }
    return 0;
}



void firwin(int n, double fln, double fhn,double h[])
{
    int i,n2,mid;
    double s,wc1,wc2,delay;
    double pis = 3.141592653589;

    if(0 == (n%2))
    {
        n2  = n/2-1;
        mid = 1;
    }
    else
    {
        n2=n/2;
        mid=0;
    }
    
    delay=n/2.0;
    
    wc1=2.0*pis*fln;
    wc2=2.0*pis*fhn;

    for(i=0;i<=n2;i++)
    {
        s=i-delay;
        h[i]=(sin(wc2*s)-sin(wc1*s))/(pis*s);
        h[i]=h[i]*hanningwin(n,i);
        h[n-i]=h[i];
    }
    
    if(mid==1)
        h[n/2]=(wc2-wc1)/pis;
}


double hanningwin(int n,int i)
{
    double w;
    double pis = 3.141592653589;
    w=0.5*(1.0-cos(2*i*pis/(n-1)));
    return(w);
}


double kaiserwin(int n,int i)
{
    double w;
    double pis = 3.141592653589;
    w=0.5*(1.0-cos(2*i*pis/(n-1)));
    return(w);
}


double kbrvel2pos (int i)
{
    int n, iter;
    double rou, rp[3], rv[3], rvp[3], rvpv[3], pos, vel, rpnew[3], ev[3],evt[3];


    for (n = 0; n < 3; n ++)
    {
        rp[n] = sat12[i].rp[n];
        rv[n] = sat12[i].rv[n];
    }

    iter = 0;
    while(1)
    {
        iter ++;            
        crsvect(rv, rp, rvp);        
        crsvect(rvp, rv, rvpv);
        pos = modvect(rp);
        vel = modvect(rv);

        if (kbrx[i].rate != 0)
            rou = sat12[i].rate * pos / vel;
        else
            rou = dotvect(rp, rv) / vel;

        for (n = 0; n < 3; n++)
        {
            evt[n] = rvpv[n] / modvect(rvpv);
            ev[n] = rv[n] / modvect(rv);
        }

        for (n = 0; n < 3; n ++)
        {
            rpnew[n] = rou * ev[n] + sqrt(pos * pos - rou * rou) * evt[n];
        }

//        printf ("iter = %d rx = %f ry = %f rz = %f dx = %f dy = %f dz = %f\n",
//            iter, rp[0], rp[1], rp[2], rp[0] - rpnew[0], rp[1] - rpnew[1], rp[2] - rpnew[2]);
            
        if (iter > 1)  break;

        for (n = 0; n < 3; n ++)
        {
            rp[n] = rpnew[n];
        }
    }
    
    
    for (n = 0; n < 3; n ++)
    {
        sat12[i].rp[n] = rpnew[n];
        sat2b[i].rp[n] = sat1a[i].rp[n] + sat12[i].rp[n];
    }
    return 0;

}




double kbrpos2vel (int i)
{
    int n, iter;
    double rou, vel, pos, rvnew[3], ep[3], ept[3], rp[3], rv[3], rpv[3], rpvp[3];

    rou = sat12[i].rate;
    for (n = 0; n < 3; n ++)
    {
        rp[n] = sat12[i].rp[n];
        rv[n] = sat12[i].rv[n];
    }

    iter = 0;
    while(1)
    {
        iter ++;            
        crsvect(rp, rv, rpv);        
        crsvect(rpv, rp, rpvp);
//        vel = modvect(rv);
        vel = sat12[i].vel;
        pos = modvect(rp);

        for (n = 0; n < 3; n++)
        {
            ept[n] = rpvp[n] / modvect(rpvp);
            ep[n] = rp[n] / modvect(rp);
        }

        for (n = 0; n < 3; n ++)
        {
            rvnew[n] = rou * ep[n] + sqrt(vel * vel - rou * rou) * ept[n];
        }                
//        printf ("iter = %d vx = %f vy = %f vz = %f dvx = %e dvy = %e dvz = %e\n",
//            iter, rv[0], rv[1], rv[2], rv[0] - rvnew[0], rv[1] - rvnew[1], rv[2] - rvnew[2]);
            
        if (iter > 1)  break;

        for (n = 0; n < 3; n ++)
        {
            rv[n] = rvnew[n];
        }
    }
    
    
    for (n = 0; n < 3; n ++)
    {
        sat12[i].rv[n] = rvnew[n];
        sat2b[i].rv[n] = sat1a[i].rv[n] + sat12[i].rv[n];
    }
    return 0;

}















/*

double alignrrr(int i)
{
    double p1[3], p2[3], v1[3], v2[3], range, rate, rou0, modev, v12[3], p12[3],
        pxv[3], pxvxp[3], vxp[3], vxpxv[3], ep[3], ept[3], ev[3], evt[3];
    int i_b, i_a, n, iter;



        iter = 0;
        while (1)
        {

            iter ++;
            crsvect(data1c[i].p12m, data1c[i].v12m, pxv);
            crsvect(pxv, data1c[i].p12m, pxvxp);
             
  
            for (n=0;n<3;n++)
            {
                ept[n] = pxvxp[n] / modvect(pxvxp);
                ep[n] = data1c[i].p12m[n] / modvect(data1c[i].p12m);
            }

//            data1c[i].vel12m = sqrt(data1c[i].ratem * data1c[i].ratem + data1c[i].pos12m * data1c[i].acclm  - dotvect(data1c[i].p12m, data1c[i].a12m) );
//            data1c[i].vel12m = sqrt(kbrx[i].rate * kbrx[i].rate + data1c[i].pos12m * kbrx[i].accl  - dotvect(data1c[i].p12m, data1c[i].a12m) );
            data1c[i].vel12m = modvect(data1c[i].v12m);
            for (n=0;n<3;n++)
            {
                data1c[i].v12m[n] = data1c[i].ratem * ep[n] + sqrt(data1c[i].vel12m * data1c[i].vel12m - data1c[i].ratem * data1c[i].ratem) * ept[n];
            }


            if (iter>2)
            {
                break;
            }
        }
    for (n = 0; n < 3; n ++ )
            data1c[i].v2m[n] =  data1c[i].v12m[n] + data1c[i].v1m[n];


    data1c[i].vel2m = modvect(data1c[i].v2m);

//        printf("%d\t%f\t%f\n", i,kbrx[i].rate, sqrt(data1c[i].pos12m * kbrx[i].accl  - dotvect(p12, data1c[i].a12)));



    return 0;

}

*/



/*


double fnlgdr(double t, int n, int m)
{
    double *pbar;

    pbar = (double *) calloc ( n - m + 1, sizeof(double));
    
    lgdr(t, n, m, pbar);

    return pbar[n-m];

}

double soliddiff (int num, double *llr1, double *llr2, 
            double *spd, double *acc12, double *dvdt)
{
    double gp1, gp2, *stcs, *pt1, *pt2, dv1, dv2, acc1[3], acc2[3];

    int nmax, id_perm;

    nmax = 4;
        
    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt2  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    id_perm = PERM;
//    stidecs(num, id_perm, stcs);
    stidecs_Anelastic(num, 1, stcs);
    
//    cs2pt (llr1, stcs, GMA[0], GMA[1], nmax, &gp1, pt1);
//    cs2pt (llr2, stcs, GMA[0], GMA[1], nmax, &gp2, pt2);
    
    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);
    cs2acc (num, llr2, stcs, GMA[0], GMA[1], nmax, &gp2, &dv2, acc2);
 

    acc12[0] = acc2[0] - acc1[0];
    acc12[1] = acc2[1] - acc1[1];
    acc12[2] = acc2[2] - acc1[2];
   
    *spd = gp2 - gp1;
    free(pt1);
    free(pt2);
    free(stcs);
    return 0;
}

*/







//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_Anelastic(InfStruct *info, int id_perm, double *stcs)
{

//    double gms2e  =  332946.048166;                   
//    double gmm2e  =  1/81.3005690699;
//    double gms2e  =  332946.0487185;                   
//    double gmm2e  =  1/81.3005538970823;

    double GMsun = 1.32712442076e20;
    double gms2e, gmm2e  = 0.0123000383;

//    double c20pt = -4.1736e-9;
//    double c20pt = -4.201e-9;
    double c20pt = C20PERM;

    double k20   =  0.29525;
    double k21   =  0.29470;
    double k22   =  0.29801;

    double REk20 =  0.30190;
    double REk21 =  0.29830;
    double REk22 =  0.30102;
    double IMk21 = -0.00144;
    double IMk22 = -0.00130;
    double k20pa = -0.00089;
    double k21pa = -0.00080;
    double k22pa = -0.00057;

    double k20p  = -0.00087;
    double k21p  = -0.00079;
    double k22p  = -0.00057;

    double k30   =  0.093;
    double k31   =  0.093;
    double k32   =  0.093;
    double k33   =  0.094;

    short int moon = 9, earth = 2, sun = 10, n;

    double ps[3], vs[3], pm[3], vm[3], tjd[2],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f, c21p, s21p, m1, m2;

    double GM, radius;

    GM = 398600.44180E+09;
    radius = 6378136.6;

    gms2e = GMsun/GM;

    tjd[0] = info->jd0;
    tjd[1] = info->tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (info->c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (info->c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = radius / llrm[2];
    rers = radius / llrs[2];

// Frequency Independent Terms

// C20
    c20 = REk20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


    s21 = - IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


// C22/S22
    c22 = + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    s22 = - IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20pa/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
    
// Frequency Dependent Terms

    if(id_perm==1) 
    {
        stcs[2]  = stcs[2] - c20pt; 
    }
    
    if (STIDE == 3)
        return 0;

    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
    stfrqdep(info->jdt, info->gmst, &c20f, &c21f, &s21f, &c22f, &s22f);


//    mtpole (info->mjd, info->xp, info->yp, &m1, &m2);
//    c21p = -1.333e-9 * (m1 + 0.0115 * m2);
//    s21p = -1.333e-9 * (m2 - 0.0115 * m1);
    c21p = 0; s21p = 0;

    stcs[2]  += c20f; 
    stcs[6]  += c21f + c21p; 
    stcs[10] += s21f + s21p; 
    stcs[13] += c22f; 
    stcs[16] += s22f; 

    return 0;

}









double stfrqdep(double jdt, double gmst, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f)
{
    double sets[71][8] = {
     0,5,5,5,6,5,  16.6e-12,  -6.7e-12, 
     0,5,5,5,7,5,  -0.1e-12,   0.1e-12, 
     0,5,6,5,5,4,  -1.2e-12,   0.8e-12, 
     0,5,7,5,5,5,  -5.5e-12,   4.3e-12, 
     0,5,7,5,6,5,   0.1e-12,  -0.1e-12, 
     0,5,8,5,5,4,  -0.3e-12,   0.2e-12, 
     0,6,3,6,5,5,  -0.3e-12,   0.7e-12, 
     0,6,5,4,4,5,   0.1e-12,  -0.2e-12, 
     0,6,5,4,5,5,  -1.2e-12,   3.7e-12, 
     0,6,5,4,6,5,   0.1e-12,  -0.2e-12, 
     0,6,5,6,5,5,   0.1e-12,  -0.2e-12, 
     0,7,3,5,5,5,   0.0e-12,   0.6e-12, 
     0,7,5,3,5,5,   0.0e-12,   0.3e-12, 
     0,7,5,5,5,5,   0.6e-12,   6.3e-12, 
     0,7,5,5,6,5,   0.2e-12,   2.6e-12, 
     0,7,5,5,7,5,   0.0e-12,   0.2e-12, 
     0,8,3,6,5,5,   0.1e-12,   0.2e-12, 
     0,8,5,4,5,5,   0.4e-12,   1.1e-12, 
     0,8,5,4,6,5,   0.2e-12,   0.5e-12, 
     0,9,3,5,5,5,   0.1e-12,   0.2e-12, 
     0,9,5,3,5,5,   0.1e-12,   0.1e-12, 
     1,2,5,7,5,5,  -0.1e-12,   0.0e-12, 
     1,2,7,5,5,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,4,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,5,5,  -0.7e-12,   0.1e-12, 
     1,3,7,4,5,5,  -0.1e-12,   0.0e-12, 
     1,4,5,5,4,5,  -1.3e-12,   0.1e-12, 
     1,4,5,5,5,5,  -6.8e-12,   0.6e-12, 
     1,4,7,5,5,5,   0.1e-12,   0.0e-12, 
     1,5,3,6,5,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,5,5,   0.4e-12,   0.0e-12, 
     1,5,5,6,5,5,   1.3e-12,  -0.1e-12, 
     1,5,5,6,6,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,5,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,6,5,   0.1e-12,   0.0e-12, 
     1,6,2,5,5,6,  -1.9e-12,   0.1e-12, 
     1,6,3,5,4,5,   0.5e-12,   0.0e-12, 
     1,6,3,5,5,5, -43.4e-12,   2.9e-12, 
     1,6,4,5,5,4,   0.6e-12,   0.0e-12, 
     1,6,4,5,5,6,   1.6e-12,  -0.1e-12, 
     1,6,5,3,4,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,3,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,4,5,  -8.8e-12,   0.5e-12, 
     1,6,5,5,5,5, 470.9e-12, -30.2e-12, 
     1,6,5,5,6,5,  68.1e-12,  -4.6e-12, 
     1,6,5,5,7,5,  -1.6e-12,   0.1e-12, 
     1,6,6,4,5,5,   0.1e-12,   0.0e-12, 
     1,6,6,5,4,4,  -0.1e-12,   0.0e-12, 
     1,6,6,5,5,4, -20.6e-12,  -0.3e-12, 
     1,6,6,5,5,6,   0.3e-12,   0.0e-12, 
     1,6,6,5,6,4,  -0.3e-12,   0.0e-12, 
     1,6,7,3,5,5,  -0.2e-12,   0.0e-12, 
     1,6,7,3,6,5,  -0.1e-12,   0.0e-12, 
     1,6,7,5,5,5,  -5.0e-12,   0.3e-12, 
     1,6,7,5,6,5,   0.2e-12,   0.0e-12, 
     1,6,8,5,5,4,  -0.2e-12,   0.0e-12, 
     1,7,3,6,5,5,  -0.5e-12,   0.0e-12, 
     1,7,3,6,6,5,  -0.1e-12,   0.0e-12, 
     1,7,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,7,5,4,5,5,  -2.1e-12,   0.1e-12, 
     1,7,5,4,6,5,  -0.4e-12,   0.0e-12, 
     1,8,3,5,5,5,  -0.2e-12,   0.0e-12, 
     1,8,5,3,5,5,  -0.1e-12,   0.0e-12, 
     1,8,5,5,5,5,  -0.6e-12,   0.0e-12, 
     1,8,5,5,6,5,  -0.4e-12,   0.0e-12, 
     1,8,5,5,7,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,5,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,6,5,  -0.1e-12,   0.0e-12, 
     2,4,5,6,5,5,  -0.3e-12,   0.0e-12, 
     2,5,5,5,5,5,  -1.2e-12,   0.0e-12  
    };

    double doodarg[6], ang, c20 = 0, c21 = 0, s21 = 0, c22 = 0, s22 = 0;
    int i, nsets = 71, argn[6], ncon = 1;

    for (i=0;i<nsets;i++)
    {
        argn[0] = (int)sets[i][0];
        argn[1] = (int)sets[i][1] - 5;
        argn[2] = (int)sets[i][2] - 5;
        argn[3] = (int)sets[i][3] - 5;
        argn[4] = (int)sets[i][4] - 5;
        argn[5] = (int)sets[i][5] - 5;
        
//        DOODSN(&info[num].jdt, &info[num].gmst, argn, &ncon, doodarg, &ang);
//        DOODSN(&jdt, &gmst, argn, &ncon, doodarg, &ang);
    
        arg2theta (jdt, gmst, argn, &ang);
        
        

//  C20 correction: Long period tidal constituent

        if(argn[0]==0) 
        {
           c20 = c20 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C21/S21 correction: Diurnal period tidal constituent

        if(argn[0]==1) 
        {
            c21 = c21 + sets[i][6]*sin(ang) + sets[i][7]*cos(ang);
            s21 = s21 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C22/S22 correction: Semi-diurnal period tidal constituent

        if(argn[0]==2) 
        {
            c22 = c22 + sets[i][6]*cos(ang);
            s22 = s22 - sets[i][6]*sin(ang);
        }
    }

    *c20f = c20;
    *c21f = c21;
    *s21f = s21;
    *c22f = c22;
    *s22f = s22;



    return 0;

}




double arg2theta (double jdt, double gmst, int n[], double *ang)
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


















/*



double stidecs_Anelastic(int num, int id_perm, double *stcs)
{

//    double gms2e  =  332946.048166;                   
//    double gmm2e  =  1/81.3005690699;
//    double gms2e  =  332946.0487185;                   
//    double gmm2e  =  1/81.3005538970823;

    double GMsun = 1.32712442076e20;
    double gms2e, gmm2e  = 0.0123000383;

//    double c20pt = -4.1736e-9;
    double c20pt = C20PERM;

    double k20   =  0.29525;
    double k21   =  0.29470;
    double k22   =  0.29801;

    double REk20 =  0.30190;
    double REk21 =  0.29830;
    double REk22 =  0.30102;
    double IMk21 = -0.00144;
    double IMk22 = -0.00130;
    double k20pa = -0.00089;
    double k21pa = -0.00080;
    double k22pa = -0.00057;

    double k20p  = -0.00087;
    double k21p  = -0.00079;
    double k22p  = -0.00057;

    double k30   =  0.093;
    double k31   =  0.093;
    double k32   =  0.093;
    double k33   =  0.094;

    short int moon = 9, earth = 2, sun = 10, n;

    double ps[3], vs[3], pm[3], vm[3], tjd[2],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f;

    double GM, radius;

//    GM = 398600.44150E+09;
//    radius = 6378136.3;

    GM = 398600.44150E+09;
    radius = 6378136.3;
    gms2e = GMsun/GM;

    tjd[0] = info[num].jd0;
    tjd[1] = info[num].tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (info[num].c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (info[num].c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = radius / llrm[2];
    rers = radius / llrs[2];

// Frequency Independent Terms

// C20
    c20 = REk20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


    s21 = - IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


// C22/S22
    c22 = + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    s22 = - IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20pa/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
    
// Frequency Dependent Terms

    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
//    stfrqdep(info->jdt, info->gmst, &c20f, &c21f, &s21f, &c22f, &s22f);


    stcs[2]  = c20 + c20f; 
    stcs[6]  = c21 + c21f; 
    stcs[10] = s21 + s21f; 
    stcs[13] = c22 + c22f; 
    stcs[16] = s22 + s22f; 

    if(id_perm==1) 
    {
        stcs[2]  = c20 + c20f - c20pt; 
    }
    return 0;

}


*/



/*

double getvrm (int num, double *coef, int nmax, double *llr1, double *llr2, double *dvdt)
{
    double gp1, gp2, *stcs, *pt1, *pt2, dv1, dv2, acc1[3], acc2[3];
    double p1e[3], p2e[3], p1[3], p2[3], p1e_a[3], p1e_b[3], p2e_a[3], p2e_b[3], 
            llr1_b[3], llr1_a[3], llr2_a[3], llr2_b[3];
        
    double *stcs_b, gp1_b, gp2_b, *stcs_a, gp1_a, gp2_a, dv1m, dv2m;
    int num_b, num_a;


   
    
    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }

    llh2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);

    llh2xyz (llr2, p2e);
    brmul(info[num].c_ei, p2e, 3, 3, 1, p2);
    brmul(info[num_b].c_ie, p2, 3, 3, 1, p2e_b);
    brmul(info[num_a].c_ie, p2, 3, 3, 1, p2e_a);
    xyz2llh(p2e_b, llr2_b);
    xyz2llh(p2e_a, llr2_a);

    cs2acc (num_b, llr1_b, coef, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
    cs2acc (num_b, llr2_b, coef, GMA[0], GMA[1], nmax, &gp2_b, &dv2, acc2);
    cs2acc (num_a, llr1_a, coef, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
    cs2acc (num_a, llr2_a, coef, GMA[0], GMA[1], nmax, &gp2_a, &dv2, acc2);

    dv1m = (gp1_a - gp1_b)/DT/2.0;
    dv2m = (gp2_a - gp2_b)/DT/2.0;

    *dvdt = dv2m - dv1m;
//    *dvdt = dv2 - dv1;

    return 0;
}



*/






double otidev12 (int num, int nmax, double *llr1, double *llr2)
{
    double gp1, *pt1, *stcs, dv1, acc1[3], p1e[3], p1[3], p1e_a[3], p1e_b[3], 
        llr1_b[3], llr1_a[3], *stcs_b, gp1_b, *stcs_a, gp1_a, rpv[3];
    int num_b, num_a, n;

    if (nmax < 2)
    {
        return 0;
    }

    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_b = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_a = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    otidecs_csr(&info[num], nmax, stcs);

    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &sat1a[num].vo, &dv1, sat1a[num].ao);
    cs2acc (num, llr2, stcs, GMA[0], GMA[1], nmax, &sat2b[num].vo, &dv1, sat2b[num].ao);


    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }

    otidecs_csr(&info[num_b], nmax, stcs_b);
    otidecs_csr(&info[num_a], nmax, stcs_a);

    llh2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);
    cs2pt (llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, pt1);
    cs2pt (llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, pt1);
    sat1a[num].dvto = (gp1_a - gp1_b)/DT/2.0;


    llh2xyz (llr2, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);
    cs2pt (llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, pt1);
    cs2pt (llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, pt1);
    sat2b[num].dvto = (gp1_a - gp1_b)/DT/2.0;




    sat12[num].vo = sat2b[num].vo - sat1a[num].vo;
    sat12[num].dvto = sat2b[num].dvto - sat1a[num].dvto;

    crsvect(sat1a[num].rp, sat1a[num].ao, rpv);
    sat1a[num].dvro = dotvect(rpv, info[num].wi);
    crsvect(sat2b[num].rp, sat2b[num].ao, rpv);
    sat2b[num].dvro = dotvect(rpv, info[num].wi);

    sat12[num].dvro = sat2b[num].dvro - sat1a[num].dvro;

    for (n=0;n<3;n++)
        sat12[num].ao[n] = sat2b[num].ao[n] - sat1a[num].ao[n];
        

    sat1a[num].dvao = dotvect(sat1a[num].ao, sat1a[num].rv);
    sat2b[num].dvao = dotvect(sat2b[num].ao, sat2b[num].rv);
    sat12[num].dvao = sat2b[num].dvao - sat1a[num].dvao;





    free(pt1);
    free(stcs);
    free(stcs_b);
    free(stcs_a);
    return 0;
}





double openaod (char *file_aod, int nmax, int mmax)
{
    FILE *fp_aod;
    double c00, s00, c06, s06, c12, s12, c18, s18, c24, s24;
    int n,m, l, ind, dim_x;
    char string[400], name[20];

    if ((fp_aod = fopen (file_aod,"r")) == NULL)
    {
        printf ("Cannot open AOD file?\n");
        exit (0);
    }


    dim_x = (nmax + 1) * (nmax + 1) + 1;

    AOD_EPH[0 * dim_x] = 0.0;
    AOD_EPH[1 * dim_x] = 0.25;
    AOD_EPH[2 * dim_x] = 0.5;
    AOD_EPH[3 * dim_x] = 0.75;
    AOD_EPH[4 * dim_x] = 1.0;



    AOD_EPH[0 * dim_x + 1] = 0;
    AOD_EPH[1 * dim_x + 1] = 0;
    AOD_EPH[2 * dim_x + 1] = 0;
    AOD_EPH[3 * dim_x + 1] = 0;
    AOD_EPH[4 * dim_x + 1] = 0;


    while (1)
    {
        if (fgets (string, 400, fp_aod) == NULL) break;
        sscanf (string, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &n, &m, &c00, &s00, &c06, &s06, &c12, &s12, &c18, &s18, &c24, &s24);

        if (n > nmax || n < 2 || m > mmax)   // permanently exclude degree 1 @7/24/2012
            continue;
        else if (m == 0)
        {
            AOD_EPH[0 * dim_x + 1 + n] = c00;
            AOD_EPH[1 * dim_x + 1 + n] = c06;
            AOD_EPH[2 * dim_x + 1 + n] = c12;
            AOD_EPH[3 * dim_x + 1 + n] = c18;
            AOD_EPH[4 * dim_x + 1 + n] = c24;
        }
        else
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            AOD_EPH[0 * dim_x + 1 + ind + n - m] = c00;
            AOD_EPH[1 * dim_x + 1 + ind + n - m] = c06;
            AOD_EPH[2 * dim_x + 1 + ind + n - m] = c12;
            AOD_EPH[3 * dim_x + 1 + ind + n - m] = c18;
            AOD_EPH[4 * dim_x + 1 + ind + n - m] = c24;
            AOD_EPH[0 * dim_x + 1 + ind + n - m + l] = s00;
            AOD_EPH[1 * dim_x + 1 + ind + n - m + l] = s06;
            AOD_EPH[2 * dim_x + 1 + ind + n - m + l] = s12;
            AOD_EPH[3 * dim_x + 1 + ind + n - m + l] = s18;
            AOD_EPH[4 * dim_x + 1 + ind + n - m + l] = s24;

        }
    }
    fclose(fp_aod);
    return 0;
}   















double lgr_order (double *y, int dim_y, int dim_x, double t, double *z, int order)
{
    int i, j, k, m, dim;
    double s;

    if (order < 1)
    {
        printf ("error: order < 1 !\n");
        exit(0);
    }
    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y))
        i = i + 1;
    k = i - order;
    if (k < 0)
        k = 0;
    m = i + order - 1;
    if (m > dim_y - 1)
        m = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
    }

    for (i = k; i <= m; i++)
    {
        s = 1.0;
        for (j = k; j <= m; j++)
        {
            if (j != i)
            {
                s = s * (t - y[j * dim_x]) / (y[i * dim_x] - y[j * dim_x]);
            }
        }
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + s * y[i * dim_x + dim + 1];
        }
    }
    return 0;
}
















double otidev (int num, int nmax, double *llr1, double *vs, double *dvts, double *as)
{
    double gp1, *stcs, dv1, acc1[3], p1e[3], p1[3], p1e_a[3], p1e_b[3], 
        llr1_b[3], llr1_a[3], *stcs_b, gp1_b, *stcs_a, gp1_a;
    int num_b, num_a;

//    nmax = NOMAX;
    
    if (nmax < 2)
    {
        *vs = 0; 
        *dvts = 0; 
        as[0] = 0; as[1] = 0; as[2] = 0;
        return 0;
    }

    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    otidecs_csr(&info[num], nmax, stcs);

    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);

    *vs = gp1;

    as[0] = acc1[0];
    as[1] = acc1[1];
    as[2] = acc1[2];
   

    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }
    stcs_b = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_a = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    otidecs_csr(&info[num_b], nmax, stcs_b);
    otidecs_csr(&info[num_a], nmax, stcs_a);

    llh2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);


    cs2acc (num_b, llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
    cs2acc (num_a, llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
    

    *dvts = (gp1_a - gp1_b)/DT/2.0;

    free(stcs);
    free(stcs_b);
    free(stcs_a);
    return 0;
}









double stidev (int num, double *llr1, double *vs, double *dvts, double *as)
{
    double gp1, *stcs, dv1, acc1[3], p1e[3], p1[3], p1e_a[3], p1e_b[3], 
        llr1_b[3], llr1_a[3], *stcs_b, gp1_b, *stcs_a, gp1_a;
    int num_b, num_a, nmax;

    nmax = 4;
        
    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

//    id_perm = PERM;
//    stidecs(num, 1, stcs);
//    stidecs_Anelastic(num, 1, stcs);
            
    stidecs_Anelastic(&info[num], 1, stcs);

    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);

    *vs = gp1;

    as[0] = acc1[0];
    as[1] = acc1[1];
    as[2] = acc1[2];
   

    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }
    stcs_b = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_a = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

//    stidecs(num_b, id_perm, stcs_b);
//    stidecs(num_a, id_perm, stcs_a);

//    stidecs_Anelastic(num_b, 1, stcs_b);
//    stidecs_Anelastic(num_a, 1, stcs_a);
    stidecs_Anelastic(&info[num_b], 1, stcs_b);
    stidecs_Anelastic(&info[num_a], 1, stcs_a);

    llh2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);


    cs2acc (num_b, llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
    cs2acc (num_a, llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
    

    *dvts = (gp1_a - gp1_b)/DT/2.0;

    free(stcs);
    free(stcs_b);
    free(stcs_a);
    return 0;
}













/*

double soliddiffdv (int num, double *llr1, double *llr2, 
            double *spd, double *acc12, double *dvdt)
{
    double gp1, gp2, *stcs, *pt1, *pt2, dv1, dv2, acc1[3], acc2[3];
    double p1e[3], p2e[3], p1[3], p2[3], p1e_a[3], p1e_b[3], p2e_a[3], p2e_b[3], 
            llr1_b[3], llr1_a[3], llr2_a[3], llr2_b[3];
        
    double *stcs_b, gp1_b, gp2_b, *stcs_a, gp1_a, gp2_a, dv1m, dv2m;
    int num_b, num_a;

    int nmax, id_perm;

    nmax = 4;
        
    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt2  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    id_perm = PERM;
    stidecs(num, 1, stcs);
   
    
    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }
    stcs_b = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs_a = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    stidecs(num_b, id_perm, stcs_b);
    stidecs(num_a, id_perm, stcs_a);


    llh2xyz (llr1, p1e);
    brmul(info[num].c_ei, p1e, 3, 3, 1, p1);
    brmul(info[num_b].c_ie, p1, 3, 3, 1, p1e_b);
    brmul(info[num_a].c_ie, p1, 3, 3, 1, p1e_a);
    xyz2llh(p1e_b, llr1_b);
    xyz2llh(p1e_a, llr1_a);

    llh2xyz (llr2, p2e);
    brmul(info[num].c_ei, p2e, 3, 3, 1, p2);
    brmul(info[num_b].c_ie, p2, 3, 3, 1, p2e_b);
    brmul(info[num_a].c_ie, p2, 3, 3, 1, p2e_a);
    xyz2llh(p2e_b, llr2_b);
    xyz2llh(p2e_a, llr2_a);



//        brmul(info[i].c_ie, data1c[i].p1, 3, 3, 1, p1e);  //from fixed acc to inertial acc
//        brmul(info[i].c_ie, data1c[i].p2, 3, 3, 1, p2e);  //from fixed acc to inertial acc

//        xyz2llh(p1e, data1c[i].llr1);
//        xyz2llh(p2e, data1c[i].llr2);


//    cs2acc (num_b, llr1, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
//    cs2acc (num_b, llr2, stcs_b, GMA[0], GMA[1], nmax, &gp2_b, &dv2, acc2);
//    cs2acc (num_a, llr1, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
//    cs2acc (num_a, llr2, stcs_a, GMA[0], GMA[1], nmax, &gp2_a, &dv2, acc2);

    cs2acc (num_b, llr1_b, stcs_b, GMA[0], GMA[1], nmax, &gp1_b, &dv1, acc1);
    cs2acc (num_b, llr2_b, stcs_b, GMA[0], GMA[1], nmax, &gp2_b, &dv2, acc2);
    cs2acc (num_a, llr1_a, stcs_a, GMA[0], GMA[1], nmax, &gp1_a, &dv1, acc1);
    cs2acc (num_a, llr2_a, stcs_a, GMA[0], GMA[1], nmax, &gp2_a, &dv2, acc2);

//    cs2pt (llr1, stcs, GMA[0], GMA[1], nmax, &gp1, pt1);
//    cs2pt (llr2, stcs, GMA[0], GMA[1], nmax, &gp2, pt2);
    
    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);
    cs2acc (num, llr2, stcs, GMA[0], GMA[1], nmax, &gp2, &dv2, acc2);
  

    acc12[0] = acc2[0] - acc1[0];
    acc12[1] = acc2[1] - acc1[1];
    acc12[2] = acc2[2] - acc1[2];
   
    *spd = gp2 - gp1;





    dv1m = (gp1_a - gp1_b)/DT/2.0;
    dv2m = (gp2_a - gp2_b)/DT/2.0;

    *dvdt = dv2m - dv1m;
//    *dvdt = dv2 - dv1;


    free(pt1);
    free(pt2);
    free(stcs);
    free(stcs_b);
    free(stcs_a);
    return 0;
}


*/





/*


//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs(int num, int id_perm, double *stcs)
{

    double gms2e  =  332946.048166;                   
    double gmm2e  =  1/81.3005690699;
//    double c20pt = -4.1736e-9;
    double c20pt = C20PERM;
    
    double k20   =  0.29525;
    double k21   =  0.29470;
    double k22   =  0.29801;
    double k20p  = -0.00087;
    double k21p  = -0.00079;
    double k22p  = -0.00057;
    double k30   =  0.093;
    double k31   =  0.093;
    double k32   =  0.093;
    double k33   =  0.094;

    short int moon = 9, earth = 2, sun = 10, n;

    double tjd[2], ps[3], vs[3], pm[3], vm[3],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f;


    tjd[0] = info[num].jd0;
    tjd[1] = info[num].tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (info[num].c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (info[num].c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = GMA[1] / llrm[2];
    rers = GMA[1] / llrs[2];

// Frequency Independent Terms

// C20
    c20 = k20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = k21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s21 = k21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C22/S22
    c22 = k22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s22 = k22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20p/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21p/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21p/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22p/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22p/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
    
// Frequency Dependent Terms
    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
//    stfrqdep(num, &c20f, &c21f, &s21f, &c22f, &s22f);


    stcs[2]  = c20 + c20f; 
    stcs[6]  = c21 + c21f; 
    stcs[10] = s21 + s21f; 
    stcs[13] = c22 + c22f; 
    stcs[16] = s22 + s22f; 

    if(id_perm==1) 
    {
        stcs[2]  = c20 + c20f - c20pt; 
    }
    return 0;

}


*/




/*
double stfrqdep(int num, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f)
{
    double sets[71][8] = {
     0,5,5,5,6,5,  16.6e-12,  -6.7e-12, 
     0,5,5,5,7,5,  -0.1e-12,   0.1e-12, 
     0,5,6,5,5,4,  -1.2e-12,   0.8e-12, 
     0,5,7,5,5,5,  -5.5e-12,   4.3e-12, 
     0,5,7,5,6,5,   0.1e-12,  -0.1e-12, 
     0,5,8,5,5,4,  -0.3e-12,   0.2e-12, 
     0,6,3,6,5,5,  -0.3e-12,   0.7e-12, 
     0,6,5,4,4,5,   0.1e-12,  -0.2e-12, 
     0,6,5,4,5,5,  -1.2e-12,   3.7e-12, 
     0,6,5,4,6,5,   0.1e-12,  -0.2e-12, 
     0,6,5,6,5,5,   0.1e-12,  -0.2e-12, 
     0,7,3,5,5,5,   0.0e-12,   0.6e-12, 
     0,7,5,3,5,5,   0.0e-12,   0.3e-12, 
     0,7,5,5,5,5,   0.6e-12,   6.3e-12, 
     0,7,5,5,6,5,   0.2e-12,   2.6e-12, 
     0,7,5,5,7,5,   0.0e-12,   0.2e-12, 
     0,8,3,6,5,5,   0.1e-12,   0.2e-12, 
     0,8,5,4,5,5,   0.4e-12,   1.1e-12, 
     0,8,5,4,6,5,   0.2e-12,   0.5e-12, 
     0,9,3,5,5,5,   0.1e-12,   0.2e-12, 
     0,9,5,3,5,5,   0.1e-12,   0.1e-12, 
     1,2,5,7,5,5,  -0.1e-12,   0.0e-12, 
     1,2,7,5,5,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,4,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,5,5,  -0.7e-12,   0.1e-12, 
     1,3,7,4,5,5,  -0.1e-12,   0.0e-12, 
     1,4,5,5,4,5,  -1.3e-12,   0.1e-12, 
     1,4,5,5,5,5,  -6.8e-12,   0.6e-12, 
     1,4,7,5,5,5,   0.1e-12,   0.0e-12, 
     1,5,3,6,5,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,5,5,   0.4e-12,   0.0e-12, 
     1,5,5,6,5,5,   1.3e-12,  -0.1e-12, 
     1,5,5,6,6,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,5,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,6,5,   0.1e-12,   0.0e-12, 
     1,6,2,5,5,6,  -1.9e-12,   0.1e-12, 
     1,6,3,5,4,5,   0.5e-12,   0.0e-12, 
     1,6,3,5,5,5, -43.4e-12,   2.9e-12, 
     1,6,4,5,5,4,   0.6e-12,   0.0e-12, 
     1,6,4,5,5,6,   1.6e-12,  -0.1e-12, 
     1,6,5,3,4,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,3,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,4,5,  -8.8e-12,   0.5e-12, 
     1,6,5,5,5,5, 470.9e-12, -30.2e-12, 
     1,6,5,5,6,5,  68.1e-12,  -4.6e-12, 
     1,6,5,5,7,5,  -1.6e-12,   0.1e-12, 
     1,6,6,4,5,5,   0.1e-12,   0.0e-12, 
     1,6,6,5,4,4,  -0.1e-12,   0.0e-12, 
     1,6,6,5,5,4, -20.6e-12,  -0.3e-12, 
     1,6,6,5,5,6,   0.3e-12,   0.0e-12, 
     1,6,6,5,6,4,  -0.3e-12,   0.0e-12, 
     1,6,7,3,5,5,  -0.2e-12,   0.0e-12, 
     1,6,7,3,6,5,  -0.1e-12,   0.0e-12, 
     1,6,7,5,5,5,  -5.0e-12,   0.3e-12, 
     1,6,7,5,6,5,   0.2e-12,   0.0e-12, 
     1,6,8,5,5,4,  -0.2e-12,   0.0e-12, 
     1,7,3,6,5,5,  -0.5e-12,   0.0e-12, 
     1,7,3,6,6,5,  -0.1e-12,   0.0e-12, 
     1,7,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,7,5,4,5,5,  -2.1e-12,   0.1e-12, 
     1,7,5,4,6,5,  -0.4e-12,   0.0e-12, 
     1,8,3,5,5,5,  -0.2e-12,   0.0e-12, 
     1,8,5,3,5,5,  -0.1e-12,   0.0e-12, 
     1,8,5,5,5,5,  -0.6e-12,   0.0e-12, 
     1,8,5,5,6,5,  -0.4e-12,   0.0e-12, 
     1,8,5,5,7,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,5,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,6,5,  -0.1e-12,   0.0e-12, 
     2,4,5,6,5,5,  -0.3e-12,   0.0e-12, 
     2,5,5,5,5,5,  -1.2e-12,   0.0e-12  
    };

    double doodarg[6], ang, c20 = 0, c21 = 0, s21 = 0, c22 = 0, s22 = 0;
    int i, nsets = 71, argn[6], ncon = 1;

    for (i=0;i<nsets;i++)
    {
        argn[0] = (int)sets[i][0];
        argn[1] = (int)sets[i][1] - 5;
        argn[2] = (int)sets[i][2] - 5;
        argn[3] = (int)sets[i][3] - 5;
        argn[4] = (int)sets[i][4] - 5;
        argn[5] = (int)sets[i][5] - 5;
        
        DOODSN(&info[num].jdt, &info[num].gmst, argn, &ncon, doodarg, &ang);

//  C20 correction: Long period tidal constituent

        if(argn[0]==0) 
        {
           c20 = c20 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C21/S21 correction: Diurnal period tidal constituent

        if(argn[0]==1) 
        {
            c21 = c21 + sets[i][6]*sin(ang) + sets[i][7]*cos(ang);
            s21 = s21 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C22/S22 correction: Semi-diurnal period tidal constituent

        if(argn[0]==2) 
        {
            c22 = c22 + sets[i][6]*cos(ang);
            s22 = s22 - sets[i][6]*sin(ang);
        }
    }

    *c20f = c20;
    *c21f = c21;
    *s21f = s21;
    *c22f = c22;
    *s22f = s22;



    return 0;

}

*/


/*
double oceandiff (int num, double *llr1, double *llr2, 
                    double *opd, double *acc12)
{
    double gp1, gp2, *stcs, *pt1, *pt2, dv1, dv2, acc1[3], acc2[3];

    int nmax, id_long;

    nmax = 100;

    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt2  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    id_long = 1;
    otidecs(num, nmax, id_long, stcs);

//    cs2pt (llr1, stcs, GMA[0], GMA[1], nmax, &gp1, pt1);
//    cs2pt (llr2, stcs, GMA[0], GMA[1], nmax, &gp2, pt2);

    cs2acc (num, llr1, stcs, GMA[0], GMA[1], nmax, &gp1, &dv1, acc1);
    cs2acc (num, llr2, stcs, GMA[0], GMA[1], nmax, &gp2, &dv2, acc2);
 

    acc12[0] = acc2[0] - acc1[0];
    acc12[1] = acc2[1] - acc1[1];
    acc12[2] = acc2[2] - acc1[2];
   

    
//    printf("%e\t%e\n", gp1, gp2);

    *opd = gp2 - gp1;
    free(pt1);
    free(pt2);
    free(stcs);
    return 0;

}
*/



double atmocdiff (int gps, double *p1, double *p2, double *apd, double *acc12)
{
    *apd = 0;
    return 0;
}

/*
int openotcs (char *infile, int *n)
{
    FILE *fp_ot;
    int i;
    char string[MAXLINE];

    if ((fp_ot = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open otide file?\n");
        exit (0);
    }

    fgets (string, MAXLINE, fp_ot);
    fgets (string, MAXLINE, fp_ot);
    fgets (string, MAXLINE, fp_ot);
    fgets (string, MAXLINE, fp_ot);

    i = 0;
    while (1)
    {
        if (fgets (string, MAXLINE, fp_ot) == NULL) break;
        sscanf (string, "%lf%s%d%d%lf%lf%lf%lf", 
            &otfes[i].ds, &otfes[i].name, &otfes[i].n, &otfes[i].m, &otfes[i].cp, &otfes[i].sp, &otfes[i].cm, &otfes[i].sm);    

        otfes[i].argn[0] = (int)(otfes[i].ds/100)%10;
        otfes[i].argn[1] = (int)(otfes[i].ds/10)%10 - 5;
        otfes[i].argn[2] = (int)(otfes[i].ds/1)%10 - 5;
        otfes[i].argn[3] = (int)(otfes[i].ds*10)%10 - 5;
        otfes[i].argn[4] = (int)(otfes[i].ds*100)%10 - 5;
        otfes[i].argn[5] = (int)(otfes[i].ds*1000)%10 - 5;
        i++;
    }

    (*n) = i;
    fclose(fp_ot);

    return 0;

}
*/


/*
double otidecs(int num, int nmax, int id_long, double *coef)
{
    double doodarg[6], ang, cp, sp, cm, sm;
    int i, ncon = 1, n,m, l, ind;

    for (i = 0; i < NFES; i++)
    {
        if (otfes[i].n > nmax)
        {
            continue;
        }

        n = otfes[i].n;
        m = otfes[i].m;
        cp = otfes[i].cp;
        sp = otfes[i].sp;
        cm = otfes[i].cm;
        sm = otfes[i].sm;

        DOODSN(&info[num].jdt, &info[num].gmst, otfes[i].argn, &ncon, doodarg, &ang);
//  ang=0;

        if (m == 0)
        {
            coef[n] = coef[n] + 1e-11 * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = coef[ind + n - m] + 1e-11 * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
            coef[ind + n - m + l] = coef[ind + n - m + l] + 1e-11 * ((sp-sm) * cos(ang) - (cp-cm)*sin(ang));
        }
    }



    return 0;

}


*/



int openotcs_csr (char *infile)
{
    FILE *fp_ot;
    int i, ncon, ntrm, nmax, mmax;
    char string[500];
    double Re, G, ge, rou, mass, pfcn, cfcn, fconst, knp[100];

    if ((fp_ot = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open otide file?\n");
        exit (0);
    }

    fgets (string, 500, fp_ot);
    fgets (string, 500, fp_ot);
    sscanf (string, "%3d%4d%4d%4d", &ncon, &ntrm, &nmax, &mmax);


//    ntrm = ntrm + 38;
    NFES = ntrm;
    otfes = (OTStruct *) calloc ( NFES, sizeof(OTStruct));

    fgets (string, 500, fp_ot);
    fgets (string, 500, fp_ot);
    sscanf (string, "%21lf%21lf%21lf%21lf%21lf", &Re, &rou, &mass, &pfcn, &cfcn);
    
//    fconst = 2 * TWOPI * rou * Re * Re / mass * 0.01;
    G = 6.67428e-11;
    ge = 9.7803278;

    fconst = 2 * TWOPI * rou * G / ge * 0.01;

    knp[0] = 0; knp[1] = 0;  

    i = 0;
    while (i<nmax)
    {
        if (fgets (string, 500, fp_ot) == NULL) break;
        sscanf (string, "%21lf%21lf%21lf%21lf%21lf%21lf", 
            &knp[i+1], &knp[i+2], &knp[i+3], &knp[i+4], &knp[i+5], &knp[i+6]);
        i = i + 6;
    }


    for (i = 0; i < ncon; i ++)
        if (fgets (string, 500, fp_ot) == NULL) break;

    for (i = 0; i < ntrm; i ++)
    {
        if (fgets (string, 500, fp_ot) == NULL) break;

//        printf ("i = %d\n", i);
        sscanf (string, "%*12c%7lf%c%c%c%c%2d%2d%lf%lf%lf%lf", 
            &otfes[i].ds, &otfes[i].name[0], &otfes[i].name[1], &otfes[i].name[2], &otfes[i].name[3], 
            &otfes[i].n, &otfes[i].m, &otfes[i].cp, &otfes[i].sp, &otfes[i].cm, &otfes[i].sm);
        otfes[i].name[4] = '\0';

        otfes[i].argn[0] = (int)(otfes[i].ds/100)%10;
        otfes[i].argn[1] = (int)(otfes[i].ds/10)%10 - 5;
        otfes[i].argn[2] = (int)(otfes[i].ds/1)%10 - 5;
        otfes[i].argn[3] = (int)(otfes[i].ds*10)%10 - 5;
        otfes[i].argn[4] = (int)(otfes[i].ds*100)%10 - 5;
        otfes[i].argn[5] = (int)(otfes[i].ds*1000)%10 - 5;
        otfes[i].fnm = fconst / normfct (otfes[i].n, otfes[i].m) 
            * (1.0 + knp[otfes[i].n]) / (2.0 * otfes[i].n + 1.0);
//            * (1.0 ) / (2.0 * otfes[i].n + 1.0);

    }

//    (*n) = i;
    fclose(fp_ot);
//    free (knp);

    return 0;

}


double normfct (int n, int m)
{
    int i;
    double nmm, npm;

    if (m == 0)
        return sqrt(2.0 * n + 1.0);
    else 
    {
        nmm = 1.0;
        npm = 1.0;
        for (i = 1; i <= n - m; i++)
            nmm = nmm * i;
        for (i = 1; i <= n + m; i++)
            npm = npm * i;
        
        return sqrt(2.0 * (2.0 * n + 1.0) * nmm / npm);
    }


}


double otidecs_csr(InfStruct *info, int nmax, double *coef)
{
    double doodarg[6], ang, cp, sp, cm, sm, fnm, m1, m2, 
        c10p, c11p, s11p, c20p, c21p, s21p;
    int i, ncon = 1, n,m, l, ind;

    for (i = 0; i < (nmax + 1) * (nmax + 1); i++)
    {
        coef[i] = 0;
    }

    for (i = 0; i < NFES; i++)
    {
        if (otfes[i].n > nmax) 
        {
            continue;
        }

        n = otfes[i].n;
        m = otfes[i].m;
        cp = otfes[i].cp;
        sp = otfes[i].sp;
        cm = otfes[i].cm;
        sm = otfes[i].sm;
        fnm = otfes[i].fnm;

        arg2theta (info->jdt, info->gmst, otfes[i].argn, &ang);
        
        if (m == 0)
        {
            coef[n] = coef[n] + fnm * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
        }
        else
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = coef[ind + n - m] + fnm * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
            coef[ind + n - m + l] = coef[ind + n - m + l] + fnm * ((sp-sm) * cos(ang) - (cp-cm)*sin(ang));
        }
    }

//    if (OTIDE == 1 || OTIDE == 3)
        return;

    mtpole (info->mjd, info->xp, info->yp, &m1, &m2);
    
//    c10p =  4.003249e-11 * (m1 + 1.603077 * m2); 
//    c11p =  6.061300e-11 * (m1 + 0.854717 * m2);
//    s11p =  1.225582e-10 * (m2 + 0.377549 * m1);
//    c20p = -3.544100e-12 * (m1 - 0.167507 * m2);  
    c21p = -2.1778e-10   * (m1 -  0.01724 * m2);
    s21p = -1.7232e-10   * (m2 -  0.03365 * m1);

    n = 2; m = 1;
    l = nmax - m + 1;
    ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
    coef[ind + n - m] = coef[ind + n - m] + c21p;
    coef[ind + n - m + l] = coef[ind + n - m + l] + s21p;

    return 0;

}





double mtpole (double mjd, double xp, double yp, double *m1, double *m2)
{
    double t, xpm, ypm, mjd2000 = 51544.0;
    int label;
    
    t = (mjd - mjd2000)/365.25;
    
    xpm = 0.055974 + t * 1.8243e-3
                   + t * t * 0.18413e-3
                   + t * t * t * 0.007024e-3;
    ypm = 0.346346 + t * 1.7896e-3
                   - t * t * 0.10729e-3
                   - t * t * t * 0.000908e-3;

    xpm = 0.054 + t * 0.00083;
    ypm = 0.357 + t * 0.00395;

    *m1 = xp - xpm;
    *m2 = - yp + ypm;

    return 0;
}









double otpole (double mjd, double xp, double yp, double *c21p, double *s21p)
{
    double t, m1, m2, xpm, ypm, mjd2000 = 51544.0;
    int label;
    
    t = (mjd - mjd2000)/365.25;
    
    xpm = 0.055974 + t * 1.8243e-3
                   + t * t * 0.18413e-3
                   + t * t * t * 0.007024e-3;
    ypm = 0.346346 + t * 1.7896e-3
                   - t * t * 0.10729e-3
                   - t * t * t * 0.000908e-3;

    xpm = 0.054 + t * 0.00083;
    ypm = 0.357 + t * 0.00395;

    m1 = xp - xpm;
    m2 = - yp + ypm;

//    printf ("m1 = %f\t m2 = %f\t xp = %f\t yp = %f\t xpm = %f\t ypm = %f\n",  
//            m1, m2, xp, yp, xpm, ypm);

    label = 1;
    *c21p = -2.1778e-10 * (m1 - label * 0.01724 * m2);
    *s21p = -1.7232e-10 * (m2 - label * 0.03365 * m1);


    return 0;
}








int openotcs_fes (char *infile)
{
    FILE *fp_ot;
    int i;
    char string[100];

    if ((fp_ot = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open otide file?\n");
        exit (0);
    }

    fgets (string, 100, fp_ot);
    fgets (string, 100, fp_ot);
    fgets (string, 100, fp_ot);
    fgets (string, 100, fp_ot);

    NFES = 59462;
    otfes = (OTStruct *) calloc ( NFES, sizeof(OTStruct));

    i = 0;
    while (1)
    {
        if (fgets (string, 100, fp_ot) == NULL) break;
        sscanf (string, "%lf%s%d%d%lf%lf%lf%lf", 
            &otfes[i].ds, &otfes[i].name, &otfes[i].n, &otfes[i].m, &otfes[i].cp, &otfes[i].sp, &otfes[i].cm, &otfes[i].sm);    

        otfes[i].argn[0] = (int)(otfes[i].ds/100)%10;
        otfes[i].argn[1] = (int)(otfes[i].ds/10)%10 - 5;
        otfes[i].argn[2] = (int)(otfes[i].ds/1)%10 - 5;
        otfes[i].argn[3] = (int)(otfes[i].ds*10)%10 - 5;
        otfes[i].argn[4] = (int)(otfes[i].ds*100)%10 - 5;
        otfes[i].argn[5] = (int)(otfes[i].ds*1000)%10 - 5;
        otfes[i].fnm = 1e-11;
        i++;
    }

//    (*n) = i;
    fclose(fp_ot);

    return 0;

}

































double getinfo(char *infile)
{
    double gmsth, jd0;


    FILE *fp_eop;
    double mjd0, mjdt, eopdat[6];
    int i, n, mjdi, lag;
    int eopflag = 0, noeop = 1, dim_eop;
    double *eop_eph, we[3], ux[3] = {1,0,0}, uy[3] = {0,1,0}, uz[3] = {0,0,1}, tx[3], ty[3], tz[3];
    char string[160];

    lag = 30; //day


    jd0 = GPS_S / 86400.0 + T0;
    mjd0 = jd0 - 2400000.5;

    if ((fp_eop = fopen (infile,"r")) == NULL)
    {   
        printf ("Cannot open eop file?\n");
        exit (0);
    }

    while (feof(fp_eop) == 0)
    {
        if (fgets (string, 160, fp_eop) ==NULL) break;
        sscanf (string, "%*d%*d%*d%d%", &mjdi);
        if (mjdi == (int)mjd0 - lag)
        {
            break;
        }
    }    
            
    dim_eop = lag + lag + lag;
    eop_eph  = (double *) calloc (dim_eop * 7, sizeof(double));
    
    for (i = 0; i < dim_eop; i++)
    {
        if (fgets (string, 160, fp_eop) ==NULL) break;
        sscanf (string, "%*d%*d%*d%lf%lf%lf%lf%lf%lf%lf", 
                &eop_eph[i * 7 + 0], &eop_eph[i * 7 + 1], &eop_eph[i * 7 + 2],
                &eop_eph[i * 7 + 3], &eop_eph[i * 7 + 4], &eop_eph[i * 7 + 5], 
                &eop_eph[i * 7 + 6]);
    }
    
    fclose(fp_eop);


    we[0] = 0; we[1] = 0; we[2] = ANGVEL;

    for (i = 0; i < NDATA; i++)
    {
        info[i].gps = sat12[i].t;

        info[i].jd0 = jd0;

        info[i].tt = info[i].gps - GPS_S + 19 + 32.184;     
        info[i].jdt = info[i].jd0 + info[i].tt / 86400.0;
        info[i].mjd = mjd0 + info[i].tt / 86400.0;
        
        info[i].leaps = getlps (info[i].jdt);

        info[i].utc = info[i].gps - GPS_S - (info[i].leaps - 19);


        mjdt = info[i].jdt - 2400000.5;

//  printf("mjdt = %d\n", DT);
//        lagrange (eop_eph, dim_eop, 7, mjdt, eopdat);
        lgr_order (eop_eph, dim_eop, 7, mjdt, eopdat, 1);
    
        for (n = 0; n < 6; n++)
        {
            if (eopdat[n] > 2.0)
            {
                printf ("eop err in lagrange! mjdt = %f\n", mjdt);
                exit(0);
            }
        }

        info[i].xp      = eopdat[0];
        info[i].yp      = eopdat[1];
        info[i].ut1_utc = eopdat[2];
        info[i].dx      = eopdat[4];
        info[i].dy      = eopdat[5];

//        geteop (utc, &xp, &yp, &ut1_utc, &dx, &dy);   

        info[i].deltat = 32.184 + info[i].leaps - info[i].ut1_utc;
        info[i].ut1 = info[i].utc + info[i].ut1_utc;

        sidereal_time (info[i].jd0, info[i].ut1/86400.0, info[i].deltat,0,1,ACCURACY, &gmsth);

        info[i].gmst = gmsth / 24 * 360.0 * DEG2RAD;


        cel_pole (info[i].jd0 + info[i].tt / 86400.0, 2, 
            info[i].dx * 1e3, info[i].dy * 1e3);

        cel2ter (info[i].jd0, info[i].ut1 / 86400.0, info[i].deltat, 1, ACCURACY, 0,
            info[i].xp, info[i].yp, ux, tx);
        cel2ter (info[i].jd0, info[i].ut1 / 86400.0, info[i].deltat, 1, ACCURACY, 0,
            info[i].xp, info[i].yp, uy, ty); 
        cel2ter (info[i].jd0, info[i].ut1 / 86400.0, info[i].deltat, 1, ACCURACY, 0,
            info[i].xp, info[i].yp, uz, tz);  

        for (n = 0; n < 3; n++)
        {
            info[i].c_ie[n*3] = tx[n];   //ie: i to e
            info[i].c_ie[n*3+1] = ty[n];
            info[i].c_ie[n*3+2] = tz[n];
        }
        
        mt(info[i].c_ie, 3, 3, info[i].c_ei);


        cip2gcrf (i, we, info[i].wi);



//    iau_pns(info[i].jdt, info[i].c_ei, 2);
//    mt(info[i].c_ei, 3, 3, info[i].c_ie);

//        printf ("%d\tjd0 = %.10f\t ut1 = %.10f\t utc = %.10f\n", i,info[i].jd0, info[i].ut1, info[i].utc );
//        for (n = 0; n < 9; n++)
//            printf ("%e\n", info[i].c_ie[n]);

    }

    return 0;
}


double nbodyv (double *tjd, double *rp, double *vn, double *dvtn, double *an)
{
    double tjd_s[2], tjd_e[2], pte, pts, pt,  acc[3], dt;
    dt = 5.0;

    nbodypt (tjd, rp, &pt, acc);

    *vn = pt;
    an[0] = acc[0];
    an[1] = acc[1];
    an[2] = acc[2];


    tjd_s[0] = tjd[0];
    tjd_s[1] = tjd[1] - dt / 86400.0;
    tjd_e[0] = tjd[0];
    tjd_e[1] = tjd[1] + dt / 86400.0;

    nbodypt (tjd_s, rp, &pts, acc);
    nbodypt (tjd_e, rp, &pte, acc);


    *dvtn = (pte - pts) / dt / 2.0;

    return 0;
}













/*

double nbodydiff (int num, double *p1, double *p2, double *npd, double *acc12, double *dvdt)
{
    double tjd[2], gp1, gp2, acc1[3], acc2[3];

    tjd[0] = info[num].jd0;
    tjd[1] = info[num].tt/86400.0;

    nbodypt (tjd, p1, &gp1, acc1);
    nbodypt (tjd, p2, &gp2, acc2);

    *npd = gp2 - gp1;
            
    acc12[0] = acc2[0] - acc1[0];
    acc12[1] = acc2[1] - acc1[1];
    acc12[2] = acc2[2] - acc1[2];

    return 0;
}


double nbodydiffdv (int num, double *p1, double *p2, double *npd, double *acc12, double *dvdt)
{
    double tjd[2], gp1, gp2, acc1[3], acc2[3];
    double tjd_b[2], gp1_b, gp2_b, tjd_a[2], gp1_a, gp2_a, dv1, dv2;
    int num_b, num_a;


    tjd[0] = info[num].jd0;
    tjd[1] = info[num].tt/86400.0;

    if (num != 0)
    {
        num_b = num - 1;
    }
    else
    {
        num_b = 0;
    }

    if (num != NDATA - 1)
    {
        num_a = num + 1;
    }
    else 
    {
        num_a = NDATA - 1;
    }


    nbodypt (tjd, p1, &gp1, acc1);
    nbodypt (tjd, p2, &gp2, acc2);

    *npd = gp2 - gp1;

    acc12[0] = acc2[0] - acc1[0];
    acc12[1] = acc2[1] - acc1[1];
    acc12[2] = acc2[2] - acc1[2];

    data1c[num].a1n[0] = acc1[0];
    data1c[num].a1n[1] = acc1[1];
    data1c[num].a1n[2] = acc1[2];

    data1c[num].a2n[0] = acc2[0];
    data1c[num].a2n[1] = acc2[1];
    data1c[num].a2n[2] = acc2[2];



    tjd_b[0] = info[num_b].jd0;
    tjd_b[1] = info[num_b].tt/86400.0;

    nbodypt (tjd_b, p1, &gp1_b, acc1);
    nbodypt (tjd_b, p2, &gp2_b, acc2);

    tjd_a[0] = info[num_a].jd0;
    tjd_a[1] = info[num_a].tt/86400.0;

    nbodypt (tjd_a, p1, &gp1_a, acc1);
    nbodypt (tjd_a, p2, &gp2_a, acc2);

    dv1 = (gp1_a - gp1_b)/DT/2.0;
    dv2 = (gp2_a - gp2_b)/DT/2.0;

    *dvdt = dv2 - dv1;
//    *dvdt = dv2;


    return 0;
}

*/


double nbodypt (double *tjd, double *pi, double *ptt, double *acc)
{
    double pj[3],vj[3], pij[3],f[3], rij, ri, rj, gm[11], 
        ptt1, ptt2, coszt, zt;
    short int earth, j, n;
    
    gm[0] =   2.203208082807623e+13;
    gm[1] =      3.248586038641429e+14;
    gm[2] =     398600.44180E+09;
    gm[3] =     4.28283719012840e+13;
    gm[4] =      1.267127698227696e+17;
    gm[5] =     3.794062664949063e+16;
    gm[6] =      5.794549096929744e+15;
    gm[7] =     6.836534169987595e+15;
    gm[8] =    9.816009029289940e+11;
    gm[9] =      4.902801056E+12;
    gm[10] =      1.32712442076e20;


    earth = 2;
    ptt1 = 0; ptt2 = 0; f[0] = 0; f[1] = 0; f[2] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (j == earth)
           continue;
        planet_ephemeris (tjd, j, earth, pj, vj);
        for (n = 0; n < 3; n++)
        {
            pj[n]  = pj[n] * AU;
            pij[n] = pj[n]  - pi[n];
        }
        rij= sqrt(pij[0] * pij[0] + pij[1] * pij[1] + pij[2] * pij[2]);
        ri = sqrt(pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2]);
        rj = sqrt(pj[0] * pj[0] + pj[1] * pj[1] + pj[2] * pj[2]);

        coszt = (pi[0] * pj[0] + pi[1] * pj[1] + pi[2] * pj[2]) / ri/ rj; 
        
        zt = acos(coszt);

        ptt1 = ptt1 + gm[j] * (1 / rij - 1 / rj - 1 / rj / rj * ri * coszt);
//        ptt1 = ptt1 + gm[j] * (1 / rij);

        ptt2 = ptt2 + gm[j] * ( 
            + pow(ri,2)/pow(rj,3) * ( 0.75*cos(2.0*zt) + 0.25)
            + pow(ri,3)/pow(rj,4) * ( 5.0/8.0*cos(3.0*zt) + 3.0/8.0*cos(zt) ) ); 
        for (n = 0; n < 3; n++)
            f[n] = f[n]
            + gm[j] / (rij * rij * rij) * pij[n] 
            - gm[j] / (rj * rj * rj) * pj[n];
        
    }
    

    *ptt = ptt1;
    
    for (n = 0; n < 3; n++)    
        acc[n] = f[n];


    return 0;
}






short int getlps (double jd)
{
/*
 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S
 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S
 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S
*/

    short int lps;

    if (jd >= 2456109.5)
        lps = 35;
    else if (jd >= 2454832.5)
        lps = 34;
    else if (jd >= 2453736.5)
        lps = 33;
    else if (jd >= 2451179.5)
        lps = 32;

    return lps;  
        
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* itrf2icrf ¨C µØ¹ÌÏµµ½¿Õ¹ÌÏµ×ª»»
* @param1: description of param1
* @param2: description of param2
* ¸Ä½ø: 
        
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double icrf2itrf(int num, double *vc, double *vt)
{

    cel_pole (info[num].jd0 + info[num].tt / 86400.0, 2, 
        info[num].dx * 1e3, info[num].dy * 1e3);

    cel2ter (info[num].jd0, info[num].ut1 / 86400.0, info[num].deltat, 1, ACCURACY, 0,
        info[num].xp, info[num].yp, vc, vt); /*--vc µ¥Î»m--*/    


    return 0;


}


void llh2xyz (double *llh, double *vt)
{
    double r, lon, lat;
    lat = llh[0] * DEG2RAD;
    lon = llh[1] * DEG2RAD;
    r = llh[2];
    vt[0] = r * cos(lat) * cos(lon);
    vt[1] = r * cos(lat) * sin(lon);
    vt[2] = r * sin(lat);

}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* xyz2llh ¨C xyz×ø±ê×ª»»µ½(Î³¶È, ¾­¶È, ¸ß¶È) (latitude, longitude, height)
* @param1: description of param1
* @param2: description of param2
* ¸Ä½ø: 
        1 Æ«µ¼Êý: ÏµÍ³²î, ÏµÍ³²î±äÂÊ
        2 µ¥³Ì¾àÀë...
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void xyz2llh (double *vt, double *llh)
{
    double r, cosphi, phi, costhe, sinthe;
    r = sqrt (vt[0] * vt[0] + vt[1] * vt[1] + vt[2] * vt[2]);

    cosphi = vt[2] / r;
    phi = acos(cosphi) ;
    costhe = vt[0] / r / sin(phi);
    sinthe = vt[1] / r / sin(phi);
    llh[2] = r;
    llh[1] = chosephase(sinthe, costhe) * RAD2DEG;
    llh[0] = 90.0 - phi * RAD2DEG;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* chosephase ¨C ÅÐ¶ÏÏóÏÞ
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double chosephase (double sinvalue, double cosvalue)
{
    double sv = sinvalue, cv = cosvalue;
    if (sv >= 0 && cv >= 0) 
        return (asin (sv));
    if (sv > 0 && cv < 0) 
        return (acos (cv));
    if (sv < 0 && cv < 0) 
        return ( - asin (sv) + TWOPI / 2.0);
    else 
        return (asin (sv) + TWOPI);
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lagrange (double *y, int dim_y, int dim_x, double t, double *z)
{ 
    int i, j, k, m, dim, order = 8;
    double s;
    
    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y)) 
        i = i + 1;
    k = i - order;
    if (k < 0) 
        k = 0;
    m = i + order - 1;
    if (m > dim_y - 1) 
        m = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
    }

    for (i = k; i <= m; i++)
    { 
        s = 1.0; 
        for (j = k; j <= m; j++)
        {
            if (j != i) 
            {
                s = s * (t - y[j * dim_x]) / (y[i * dim_x] - y[j * dim_x]);
            }
        }        
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + s * y[i * dim_x + dim + 1];
        }
    }
    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lagrangelow (double *y, int dim_y, int dim_x, double t, double *z)
{ 
    int i, j, k, m, dim, order = 2;
    double s;
    
    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y)) 
        i = i + 1;
    k = i - order;
    if (k < 0) 
        k = 0;
    m = i + order - 1;
    if (m > dim_y - 1) 
        m = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
    }

    for (i = k; i <= m; i++)
    { 
        s = 1.0; 
        for (j = k; j <= m; j++)
        {
            if (j != i) 
            {
                s = s * (t - y[j * dim_x]) / (y[i * dim_x] - y[j * dim_x]);
            }
        }        
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + s * y[i * dim_x + dim + 1];
        }
    }
    return 0;
}























/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double disse_a (double *v1, int i, double *ef1i, double *f1)
{
    double mat1[9], mtt1[9], qvec1[4], acc1[3], v1s[3];

    if (scaa[i].gps_time == 0 || scab[i].gps_time == 0)
    {
        *ef1i = 0;
        return 1;
    }

    qvec1[0] = scaa[i].quatangle;
    qvec1[1] = scaa[i].quaticoeff;
    qvec1[2] = scaa[i].quatjcoeff;
    qvec1[3] = scaa[i].quatkcoeff;
    acc1[0] = acca[i].lin_accl_x;
    acc1[1] = acca[i].lin_accl_y;
    acc1[2] = acca[i].lin_accl_z;

    quat2mat_i2s (qvec1, mat1);
//    brmul(mat1, acc1, 3,3,1, f1);
    mt(mat1, 3, 3, mtt1);
    brmul(mtt1, acc1, 3,3,1, f1);

    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]);

//    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]) * DT;

//    brmul(mat1, v1, 3, 3, 1, v1s);

//    *ef1i = (acc1[0] * v1s[0] + acc1[1] * v1s[1] + acc1[2] * v1s[2]) * DT;

    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/



double disse_b (double *v1, int i, double *ef1i, double *f1)
{
    double mat1[9], mtt1[9], qvec1[4], acc1[3], v1s[3];

    if (scaa[i].gps_time == 0 || scab[i].gps_time == 0)
    {
        *ef1i = 0;
        return 1;
    }

    qvec1[0] = scab[i].quatangle;
    qvec1[1] = scab[i].quaticoeff;
    qvec1[2] = scab[i].quatjcoeff;
    qvec1[3] = scab[i].quatkcoeff;
    acc1[0] = accb[i].lin_accl_x;
    acc1[1] = accb[i].lin_accl_y;
    acc1[2] = accb[i].lin_accl_z;

    quat2mat_i2s (qvec1, mat1);
//    brmul(mat1, acc1, 3,3,1, f1);
    mt(mat1, 3, 3, mtt1);
    brmul(mtt1, acc1, 3,3,1, f1);

    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]);

//    *ef1i = (f1[0] * v1[0] + f1[1] * v1[1] + f1[2] * v1[2]) * DT;
          
//    brmul(mat1, v1, 3, 3, 1, v1s);

//    *ef1i = (acc1[0] * v1s[0] + acc1[1] * v1s[1] + acc1[2] * v1s[2]) * DT;

    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


double accia (double *tjd, int label, double *acc)
{
    double tt, as[3], qvec[4], c_is[9], c_si[9];

    tt = tjd[1] * 86400.0;
    
    if (label == 1)
    {
        lagrangelow (ACA_EPH, DIM_ACA, 4, tt, as);
        lagrangelow (SCA_EPH, DIM_SCA, 5, tt, qvec);
    }

    if (label == 2)
    {
        lagrangelow (ACB_EPH, DIM_ACB, 4, tt, as);
        lagrangelow (SCB_EPH, DIM_SCB, 5, tt, qvec);
    }


    quat2mat_i2s (qvec, c_is);
    mt(c_is, 3, 3, c_si);


    brmul(c_si, as, 3,3,1, acc);

    return 0;

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double quat2mat_i2s (double *qvec, double *mat)
/*
 * SCA_1B data file containing edited quaternion 
 * rotating inertial frame to SRF (5-sec data interval)
 * Page 17 in
 * Algorithm Theoretical Basis Document for GRACE Level-1B Data Processing V1.2
 * Sien-Chong Wu Gerhard Kruizinga Willy Bertiger
 * May 9, 2006
*/
{
    double q[4];
    
    q[0]=qvec[0];
    q[1]=qvec[1]; 
    q[2]=qvec[2]; 
    q[3]=qvec[3]; 

    mat[0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    mat[1] = 2.0 * (q[1]*q[2] + q[3]*q[0]);
    mat[2] = 2.0 * (q[1]*q[3] - q[2]*q[0]);
    mat[3] = 2.0 * (q[1]*q[2] - q[3]*q[0]);
    mat[4] = q[0] * q[0] + q[2] * q[2] - q[1] * q[1] - q[3] * q[3];
    mat[5] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
    mat[6] = 2.0 * (q[1]*q[3] + q[2]*q[0]);
    mat[7] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
    mat[8] = q[0] * q[0] + q[3] * q[3] - q[1] * q[1] - q[2] * q[2];
 
    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/



double mean (double * array, double N)
{
    int i;
    double sum = 0 ;
    for (i = 0; i < N; i++)
        sum = sum + array [i];
    return sum/N;
} // function calculating mean


double std_dev (double * array, double N)
{  
    int i;
    double sum = 0;
    double STD_DEV = 0; // returning zero's

    for (i = 0; i < N; i++)
    {
        sum = sum + array [i];
        STD_DEV = STD_DEV + pow(array [i], 2);
    }
    return sqrt ((STD_DEV/N) - (pow(sum/N,2)));
} // function calculating standard deviation




double lsfl1c1d (double *xmax, double *x, double *y, int nmax, int n, 
    double tpsec, double *ymaxflt, int nply, int ncpr, int nplycpr)
{

    int i, j, m, err, k, cpr;

    double T, *a, *f, *fmax, *ft, *fpf, *fpy, *yflt, *res, stdres;

    k = nply + ncpr * nplycpr;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    yflt = (double *) calloc ( (int)(n) , sizeof(double));
    res = (double *) calloc ( (int)(n) , sizeof(double));
    fmax = (double *) calloc ( (int)(nmax*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));




    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < nply; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = 0; j < ncpr; j = j + 2)
        {
            for (m = 0; m < nplycpr; m ++)
            {
                f[i * k + nply + j * nplycpr + m ] = pow(x[i], m) * sin (TWOPI / T * x[i]);
                f[i * k + nply + j * nplycpr + nplycpr + m] = pow(x[i], m) *  cos (TWOPI / T * x[i]);
            }
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    for ( i = 0; i < nmax; i ++)
    {
        for (j = 0; j < nply; j ++)
        {
            fmax[i*k+j] = pow(xmax[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = 0; j < ncpr; j = j + 2)
        {
            for (m = 0; m < nplycpr; m ++)
            {
                fmax[i * k + nply + j * nplycpr + m ] = pow(xmax[i], m) * sin (TWOPI / T * xmax[i]);
                fmax[i * k + nply + j * nplycpr + nplycpr + m] = pow(xmax[i], m) *  cos (TWOPI / T * xmax[i]);
            }
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    mtrans (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);


    brmul (fmax, a, nmax, k, 1, ymaxflt);

    brmul (f, a, n, k, 1, yflt);
    for (j = 0; j < n; j ++)
        res[j] = yflt[j] - y[j];

    stdres = std_dev (res, n);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    free (res);
    free (yflt);
    free (fmax);



    return stdres;
}







void fitl1c_mv (int nply, int ncpr, int nplycpr, int ndt)
{
    int np, np2, is, ie, i, n, j, nlsf;
    double t0, x, *xlsf, *ylsf, *prd, tpmean;
    
    t0 = sat12[0].t;

    np = (int)(ndt / DT);
    np2 = (int) (np/2);

    xlsf = (double *) calloc ( np, sizeof(double));
    ylsf = (double *) calloc ( np, sizeof(double));
    prd = (double *) calloc ( np, sizeof(double));


    for (i = 0; i < NDATA; i ++)
    {
        t0 = sat12[i].t;

        x = (sat12[i].t - t0) /  86400.0;
        is = i - np2;
        ie = i + np2 - 1;
        if (is < 0)
        {
            is = 0;
            ie = np - 1;
        }
        if (ie > NDATA - 1)
        {
            ie = NDATA - 1;
            is = NDATA - np;
        }

        n = 0;
        for (j = is; j <= ie; j ++)
        {
            if (sat12[j].error == 9) continue;

            xlsf[n] = (sat12[j].t - t0) / 86400.0;
            ylsf[n] = sat12[j].vl1c;
            ylsf[n] = ylsf[n] - sat12[j].vref;  ///////////////!!!!!!!!!!!//////////
//            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            prd[n] = (gnva[j].Tp + gnvb[j].Tp) / 2.0;
            n ++;
        }

        if (n < np2)
        {
            sat12[i].vmvfit = 0;
            sat12[i].vmvfitres = 999;
            continue;
        }


        nlsf = n;
        tpmean = mean (prd, nlsf);
        sat12[i].Tp2 = tpmean;

//        printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);
        sat12[i].vmvfitres = lsfl1c1d (&x, xlsf, ylsf, 1, nlsf, tpmean, &sat12[i].vmvfit, nply, ncpr, nplycpr);
//        sat12[i].ystd = fabs(sat12[i].res - sat12[i].y - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
    
    }
    
    return;

}






void fitlos_mv (int nply, int ncpr, int nplycpr, int ndt)
{
    int np, np2, is, ie, i, n, j, nlsf;
    double t0, x, *xlsf, *ylsf, *prd, tpmean;
    
    t0 = sat12[0].t;

    np = (int)(ndt / DT);
    np2 = (int) (np/2);

    xlsf = (double *) calloc ( np, sizeof(double));
    ylsf = (double *) calloc ( np, sizeof(double));
    prd = (double *) calloc ( np, sizeof(double));


    for (i = 0; i < NDATA; i ++)
    {
        x = (sat12[i].t - t0) /  86400.0;
        is = i - np2;
        ie = i + np2 - 1;
        if (is < 0)
        {
            is = 0;
            ie = np - 1;
        }
        if (ie > NDATA - 1)
        {
            ie = NDATA - 1;
            is = NDATA - np;
        }

        n = 0;
        for (j = is; j <= ie; j ++)
        {
            if (sat12[j].error == 9) continue;

            xlsf[n] = (sat12[j].t - t0) / 86400.0;
            ylsf[n] = sat12[j].glos;
            ylsf[n] = ylsf[n] - sat12[j].gref;  ///////////////!!!!!!!!!!!//////////
//            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            prd[n] = (gnva[j].Tp + gnvb[j].Tp) / 2.0;
            n ++;
        }

        if (n < np2)
        {
            sat12[i].gmvfit = 0;
            sat12[i].gmvfitres = 999;
            continue;
        }

        nlsf = n;
        tpmean = mean (prd, nlsf);

        sat12[i].gmvfitres = lsfl1c1d (&x, xlsf, ylsf, 1, nlsf, tpmean, &sat12[i].gmvfit, nply, ncpr, nplycpr);
//        sat12[i].ystd = fabs(sat12[i].res - sat12[i].y - (sat12[i].vcsr + sat12[i].vgfz + sat12[i].vjpl - 3 * sat12[i].vbl) / 3.0);
    
    }
    
    return;

}











void fitl1c_pw (int order_poly, int order_cpr)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres, stdresm;

    nmax = 5400/DT;
//    nmax = 3600/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
//        t0 = sat12[i].t;

        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].vl1c;
            res[n] = res[n] - sat12[i].vref;  ///////////////!!!!!!!!!!!//////////
            resm[n] = sat12[i].vl1cm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].vpwfit = 0;  
//                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            tpmean = 5633.5;
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].vpwfit = fit[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitm = fitm[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitres = stdres;
                sat12[i - (nmax-1) + in].Tp1 = tpmean;
//                if (stdres > 0.003)
//                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}







void fitlos_pw (int order_poly, int order_cpr)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres, stdresm;

    nmax = 5400/DT;
//    nmax = 3600/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].glos;
            res[n] = res[n] - sat12[i].gref;  ///////////////!!!!!!!!!!!//////////
//            resm[n] = sat12[i].resm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].gpwfit = 0;  
//                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
//            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].gpwfit = fit[in];   // possible useful for GPS orbit
//                sat12[i - (nmax-1) + in].fitm = fitm[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].gpwfitres = stdres;
//                if (stdres > 0.003)
//                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}















double lsf_cpr_new (double *xmax, double *x, double *y, int nmax, int n, 
    double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *fmax, *ft, *fpf, *fpy, *yf, *res, stdres;

    k = order_poly + 1 + order_cpr * 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    yf = (double *) calloc ( (int)(n) , sizeof(double));
    res = (double *) calloc ( (int)(n) , sizeof(double));
    fmax = (double *) calloc ( (int)(nmax*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));




    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    for ( i = 0; i < nmax; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            fmax[i*k+j] = pow(xmax[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            fmax[i*k+j] = sin (TWOPI / T * xmax[i]);
            fmax[i*k+j + 1] = cos (TWOPI / T * xmax[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }



    mtrans (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);


    brmul (fmax, a, nmax, k, 1, yflt);

    brmul (f, a, n, k, 1, yf);
    for (j = 0; j < n; j ++)
        res[j] = yf[j] - y[j];

    stdres = std_dev (res, n);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    free (res);
    free (yf);
    free (fmax);



    return stdres;
}









void fitres_new_overlap (int order_poly, int order_cpr, int overlap)
{
    int i, n, in, nmax;
    double *xn, *xi, *res, *resm, *prd, *fit, *fitm, tpmean, t0, stdres, stdresm;

    nmax = 5400/DT;
            
    xn  = (double *) calloc ( nmax, sizeof(double));
    res = (double *) calloc ( nmax, sizeof(double));
    resm = (double *) calloc ( nmax, sizeof(double));
    prd = (double *) calloc ( nmax, sizeof(double));
    fit = (double *) calloc ( nmax, sizeof(double));
    fitm = (double *) calloc ( nmax, sizeof(double));
    xi  = (double *) calloc ( nmax, sizeof(double));
            
    n = 0;
    t0 = sat12[0].t;
    for (i = 0; i < NDATA; i ++)
    {
        if (sat12[i].error == 0)
        {
            xn[n] = (sat12[i].t - t0) / 86400.0;
            prd[n] = (gnva[i].Tp + gnvb[i].Tp) / 2.0;
            res[n] = sat12[i].vl1c;
            resm[n] = sat12[i].vl1cm;
            n++;
        }
        xi[(i)%nmax] = (sat12[i].t - t0) / 86400.0;

        if (i == 0) continue;
        if (((i+1)%(5400/DT)) == 0)   // possible useful for GPS orbit (orbit cut at 86399s)
//        if (((i)%(5400/DT)) == 0)
        {
            if (n < nmax / 2)
            {
                for (in = 0; in < nmax; in ++)
                {
                    sat12[i - (nmax-1) + in].vpwfit = 0;  
                    sat12[i - (nmax-1) + in].error = 3; 
                }
                t0 = sat12[i+1].t;
                n = 0;
                continue;
            }
            tpmean = mean (prd, n);
//            printf ("n = %d\t i = %d\t tpmean = %f\n", n, i, tpmean);

            stdres = lsf_cpr_new (xi, xn, res, nmax, n, tpmean, fit, order_poly, order_cpr);
            stdresm =lsf_cpr_new (xi, xn, resm, nmax, n, tpmean, fitm, order_poly, order_cpr);
//            printf ("n = %d\t i = %d\t tpmean = %f\t, stdres = %f\n", n, i, tpmean, stdres);
    
            for (in = 0; in < nmax; in ++)
            {
                sat12[i - (nmax-1) + in].vpwfit = fit[in];   // possible useful for GPS orbit
                sat12[i - (nmax-1) + in].vpwfitm = fitm[in];   // possible useful for GPS orbit
                if (stdres > 0.003)
                    sat12[i - (nmax-1) + in].error = 2;
//                sat12[i - nmax + in].fit = fit[in];
            }

            t0 = sat12[i+1].t;
            n = 0;
        }

    }
        
    free(xi);
    free(xn);
    free(res);
    free(resm);
    free(prd);
    free(fit);
    free(fitm);

}




















void fitres (int order_poly, int order_cpr, int offsetcyc)
{
    int i, is, ie, n, cnum, ofs, ofe, offset;
    double *xt, *res, *resm, *prd, *fit, *fitm, tpmean;


//    offset = 5400/DT;
    offset = offsetcyc;
    is = 0;
    for (i = 0; i < NDATA; i ++)
    {
        if (i == 0 || gnva[i].r == 0) continue;
//        if ( ((gnva[i+1].lat > gnva[i].lat) && (gnva[i-1].lat > gnva[i].lat)) 
        if ( ((i%(5400/DT)) == 0)
            || i == NDATA - 1)
        {
            ie = i + 1;
            cnum = ie - is;

            if (is - offset < 0)
            {
                ofs = 0;
                ofe = offset * 2;    
            }
            else if (ie + offset >=NDATA)
            {
                ofs = offset * 2;
                ofe = 0;
            }
            else 
            {
                ofs = offset;
                ofe = offset;
            }
            
            xt = (double *) calloc ( cnum + offset * 2, sizeof(double));
            res = (double *) calloc ( cnum + offset * 2, sizeof(double));
            prd = (double *) calloc ( cnum + offset * 2, sizeof(double));
            fit = (double *) calloc ( cnum + offset * 2, sizeof(double));
            
            for (n = 0; n < cnum + offset * 2; n ++)
            {
                xt[n] = (gnva[is - ofs + n].gps_time - gnva[is].gps_time) / 86400.0;
                prd[n] = gnva[is - ofs + n].Tp;
                res[n] = sat12[is - ofs + n].vl1c;
            }
            tpmean = mean (prd, cnum + offset * 2);
//            tpmean = 5633.13;
//            printf ("cnum = %d\t is = %d\t ie = %d\t tpmean = %f\n", cnum, is, ie, tpmean);

            lsf_cpr (xt, res, cnum + offset * 2, tpmean, fit, order_poly, order_cpr);
//            lsf_cpr_day (xt, res, cnum + offset * 2, tpmean, fit, order_poly, order_cpr);

            for (n = 0; n < cnum; n ++)
            {
                sat12[is + n].vpwfit = fit[n + ofs];
            }
            free(xt);
            free(res);
            free(prd);
            free(fit);
            is = ie;
        }
    }
        

}



int lsf_cpr_day ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *ft, *fpf, *fpy;

    k = order_poly + 1 + order_cpr * 2 + 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));


    y[0] = y[4];
    y[1] = y[4];
    y[2] = y[4];


    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k-2; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
        for (j = k-2; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI  * x[i]);
            f[i*k+j + 1] = cos (TWOPI  * x[i]);
        }

    }

/*
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", f[i*k+j]);
        }

        printf("%f\t%f\n", y[i], x[i]);

    }
*/


    mtrans (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);

/*
    for ( i = 0; i < k; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", fpf[i*k+j]);
        }

        printf("%f\n", fpy[i]);

    }

 
    for (j = 0; j < k; j ++)
    {
        printf("%f\t", a[j]);
    }
    printf("\n");

*/

    brmul (f, a, n, k, 1, yflt);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    return 0;
}














int lsf_cpr ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr)
{

    int i, j, err, k, cpr;

    double T, *a, *f, *ft, *fpf, *fpy;

    k = order_poly + 1 + order_cpr * 2;

    a = (double *) calloc ( (int)(k) , sizeof(double));
    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));


    y[0] = y[4];
    y[1] = y[4];
    y[2] = y[4];


    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < order_poly + 1; j ++)
        {
            f[i*k+j] = pow(x[i], j);
        }

        T = tpsec/86400.0;
        cpr = 1;
        for (j = order_poly + 1; j < k; j = j + 2)
        {
            f[i*k+j] = sin (TWOPI / T * x[i]);
            f[i*k+j + 1] = cos (TWOPI / T * x[i]);
            cpr ++;
//            T = T / cpr;
            T = tpsec/86400.0 / cpr;
        }
        
    }

/*
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", f[i*k+j]);
        }

        printf("%f\t%f\n", y[i], x[i]);

    }
*/


    mtrans (f, n, k, ft);
    brmul (ft, f, k, n, k, fpf);
    brmul (ft, y, k, n, 1, fpy);
    err = brinv(fpf,k);
    if (err == 0)
        return 1;
    brmul (fpf, fpy, k, k, 1, a);

/*
    for ( i = 0; i < k; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            printf("%f\t", fpf[i*k+j]);
        }

        printf("%f\n", fpy[i]);

    }

 
    for (j = 0; j < k; j ++)
    {
        printf("%f\t", a[j]);
    }
    printf("\n");

*/

    brmul (f, a, n, k, 1, yflt);

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    free (a);
    return 0;
}


int lsf_poly ( double *x, double *y, int n, double *a, int k)
{

    int i, j, err;

    double *f, *ft, *fpf, *fpy;

    f = (double *) calloc ( (int)(n*k) , sizeof(double));
    ft = (double *) calloc ( (int)(n*k) , sizeof(double));
    fpf = (double *) calloc ( (int)(k*k) , sizeof(double));
    fpy = (double *) calloc ( (int)(k) , sizeof(double));

    
    for ( i = 0; i < n; i ++)
    {
        for (j = 0; j < k; j ++)
        {
            f[i*k+j] = pow(x[i], j); 
        }
    }

    mtrans (f, n, k, ft);           
    brmul (ft, f, k, n, k, fpf);    
    brmul (ft, y, k, n, 1, fpy);    
    err = brinv(fpf,k); 
    if (err == 0)
        return 1;   
    brmul (fpf, fpy, k, k, 1, a);   //¾ØÕóÏà³Ë, ¾ÍÊÇ¼ÆËã(FTPyF)-1FTPyFÕâ¸ö¾ØÕóÁË, ¾ÍÊÇËã³öÁË¹«Ê½1, ·µ»Ø¾ØÕóa¾ÍÊÇ¹«Ê½1ÀïµÄc, Ò²¾ÍËã³öÀ´ÁË

    free (f);
    free (ft);
    free (fpf);
    free (fpy);
    return 0;
}


void mtrans(double a[], int m, int n, double b[])
{
    //a[m][n] -> b[n][m]
    int i,j;
    for (i=0; i<=m-1; i++)
        for (j=0; j<=n-1; j++)
            b[j*m + i]=a[i*n + j];
    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int calbias(char *infile_a, char *infile_b)
{

    double x0a, y0a, z0a, x0b, y0b, z0b;
    FILE *fpACC_a, *fpACC_b;
    int n_acca, n_accb, i, gps;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }


    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &x0a);
    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &y0a);
    fgets(line, 100, fpACC_a);
    sscanf (line, "%lf", &z0a);


    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &x0b);
    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &y0b);
    fgets(line, 100, fpACC_b);
    sscanf (line, "%lf", &z0b);


//    printf ("%e\n%e\n%e\n%e\n%e\n%e\n", x0a, y0a, z0a, x0b, y0b, z0b);


    for (i = 0; i < DIM_ACA; i++)
    {
        ACA_EPH[i * 4 + 1] = ACA_EPH[i * 4 + 1] + x0a;
        ACA_EPH[i * 4 + 2] = ACA_EPH[i * 4 + 2] + y0a;
        ACA_EPH[i * 4 + 3] = ACA_EPH[i * 4 + 3] + z0a;
    }
    

    for (i = 0; i < DIM_ACB; i++)
    {
        ACB_EPH[i * 4 + 1] = ACB_EPH[i * 4 + 1] + x0b;
        ACB_EPH[i * 4 + 2] = ACB_EPH[i * 4 + 2] + y0b;
        ACB_EPH[i * 4 + 3] = ACB_EPH[i * 4 + 3] + z0b;
    }
   
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/










int calbiaseph(char *infile_a, char *infile_b)
{

    double *MACC_EPBS, *MACC_EPSL, tt;
    FILE *fpACC_A, *fpACC_B;
    int n_acca, n_accb, i,k,m,n, lbs, lsl, MACC_ARBS, MACC_ARSL, 
        MACC_NOBS, MACC_NOSL, MACC_PRBS, MACC_PRSL, MACC;
//    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_A = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_B = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }

    if (MACC_BIAS != 0)
    {
        MACC_ARBS = (int)(86400 / MACC_DTBS);
        MACC_NOBS = 3 * MACC_BIAS;
        MACC_PRBS = MACC_NOBS * MACC_ARBS;
        MACC_EPBS = (double *) calloc (MACC_PRBS, sizeof(double));
    }

    if (MACC_SCAL != 0)
    {
        MACC_ARSL = (int)(86400 / MACC_DTSL);
        MACC_NOSL = 3 * MACC_SCAL;
        MACC_PRSL = MACC_NOSL * MACC_ARSL;
        MACC_EPSL = (double *) calloc (MACC_PRSL, sizeof(double));
    }

/*A*/   
 
    if (MACC_BIAS != 0)
    {
        for (i = 0; i < MACC_ARBS; i ++)
        {
            fscanf (fpACC_A, "%*s%*d");
            for (n = 0; n < MACC_NOBS; n++)
            {
                fscanf(fpACC_A, "%lf", &MACC_EPBS[i * MACC_NOBS + n]);
            }
        }
    }

    if (MACC_SCAL != 0)
    {
        for (i = 0; i < MACC_ARSL; i ++)
        {
            fscanf (fpACC_A, "%*s%*d");
            for (n = 0; n < MACC_NOSL; n++)
            {
                fscanf(fpACC_A, "%lf", &MACC_EPSL[i * MACC_NOSL + n]);
            }
        }
    }

    for (i = 0; i < DIM_ACA; i++)
    {
        tt = ACA_EPH[i * 4] / 86400.0;
//        tt = ( ACA_EPH[i * 4] - 51.184) / 86400.0;
    
        if (MACC_SCAL != 0)
        {
            lsl = (int)((ACA_EPH[i * 4] - 51.184)/MACC_DTSL);
            if (lsl < 0) lsl = 0;
            if (lsl > MACC_ARSL - 1) lsl = MACC_ARSL - 1;

            for (k = 0; k < 3; k ++)
            {
                ACA_EPH[i * 4 + k + 1] = ACA_EPH[i * 4 + k + 1] * MACC_EPSL[lsl * MACC_NOSL + k];
            }
        }   

        if (MACC_BIAS != 0)
        {
            lbs = (int)((ACA_EPH[i * 4] - 51.184)/MACC_DTBS);
            if (lbs < 0) lbs = 0;
            if (lbs > MACC_ARBS - 1) lbs = MACC_ARBS - 1;

            for (k = 0; k < 3; k ++)
            {
                for (m = 0; m < MACC_BIAS; m ++)
                {
                    ACA_EPH[i * 4 + k + 1] = ACA_EPH[i * 4 + k + 1] + MACC_EPBS[lbs * MACC_NOBS + k + 3 * m] * pow (tt, m);
                }
            }
        }
    }
    
/*B*/
    
    if (MACC_BIAS != 0)
    {
        for (i = 0; i < MACC_ARBS; i ++)
        {
            fscanf (fpACC_B, "%*s%*d");
            for (n = 0; n < MACC_NOBS; n++)
            {
                fscanf(fpACC_B, "%lf", &MACC_EPBS[i * MACC_NOBS + n]);
            }
        }
    }

    if (MACC_SCAL != 0)
    {
        for (i = 0; i < MACC_ARSL; i ++)
        {
            fscanf (fpACC_B, "%*s%*d");
            for (n = 0; n < MACC_NOSL; n++)
            {
                fscanf(fpACC_B, "%lf", &MACC_EPSL[i * MACC_NOSL + n]);
            }
        }
    }

    for (i = 0; i < DIM_ACB; i++)
    {
        tt = ACB_EPH[i * 4] / 86400.0;
    
        if (MACC_SCAL != 0)
        {
            lsl = (int)((ACB_EPH[i * 4] - 19 - 32.184)/MACC_DTSL);
            if (lsl < 0) lsl = 0;
            if (lsl > MACC_ARSL - 1) lsl = MACC_ARSL - 1;

            for (k = 0; k < 3; k ++)
            {
                ACB_EPH[i * 4 + k + 1] = ACB_EPH[i * 4 + k + 1] * MACC_EPSL[lsl * MACC_NOSL + k];
            }
        }   

        if (MACC_BIAS != 0)
        {
            lbs = (int)((ACB_EPH[i * 4] - 19 - 32.184)/MACC_DTBS);
            if (lbs < 0) lbs = 0;
            if (lbs > MACC_ARBS - 1) lbs = MACC_ARBS - 1;

            for (k = 0; k < 3; k ++)
            {
                for (m = 0; m < MACC_BIAS; m ++)
                {
                    ACB_EPH[i * 4 + k + 1] = ACB_EPH[i * 4 + k + 1] + MACC_EPBS[lbs * MACC_NOBS + k + 3 * m] * pow (tt, m);
                }
            }
        }
    }
    
   

    if (MACC_BIAS != 0)
        free (MACC_EPBS);
    if (MACC_SCAL != 0)
        free (MACC_EPSL);
    fclose (fpACC_A);
    fclose (fpACC_B);
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cal_acc_01(void)
{

    double scale_x_A, scale_y_A, scale_z_A, scale_x_B, scale_y_B, scale_z_B,
        bias_x_A, bias_y_A, bias_z_A, bias_x_B, bias_y_B, bias_z_B, 
        mjd, mjd0, mjdm, gps, gpss, mjds;
    int i;

    scale_x_A = 0.9595;
    scale_y_A = 0.9797;
    scale_z_A = 0.9485;
    scale_x_B = 0.9465;
    scale_y_B = 0.9842;
    scale_z_B = 0.9303;

    gpss = ACA_EPH[0] + GPS_S - 19 - 32.184;
    mjds = gpss / 86400.0 + T0 - 2400000.5;

//    if (mjds > 55562) // Jan. 1, 2011
    {
        scale_x_A = scale_x_A * 0.98;
        scale_z_A = scale_z_A * 0.98;
        scale_x_B = scale_x_B * 0.98;
        scale_z_B = scale_z_B * 0.98;
    }


    mjdm = 52705; //March 7, 2003

    for (i = 0; i < DIM_ACA; i++)
    {
        gps = ACA_EPH[i * 4] + GPS_S - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_A = - 1.106 
                       + 2.233e-4 * ( mjd - mjd0 ) 
                       + 2.5e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   27.042 
                       + 4.46e-3  * ( mjd - mjd0 ) 
                       + 1.1e-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5486 
                       - 1.139e-6 * ( mjd - mjd0 ) 
                       + 1.7e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_A = - 1.2095 
                       - 4.128e-5 * ( mjd - mjd0 ) 
                       + 9.7e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_A =   29.3370 
                       + 6.515e-4 * ( mjd - mjd0 ) 
                       - 3.9e-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_A = - 0.5606 
                       - 2.352e-6 * ( mjd - mjd0 ) 
                       + 3.8e-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }

        ACA_EPH[i * 4 + 1] = bias_x_A * 1e-6 + scale_x_A * ACA_EPH[i * 4 + 1];
        ACA_EPH[i * 4 + 2] = bias_y_A * 1e-6 + scale_y_A * ACA_EPH[i * 4 + 2];
        ACA_EPH[i * 4 + 3] = bias_z_A * 1e-6 + scale_z_A * ACA_EPH[i * 4 + 3];

    }
    

    for (i = 0; i < DIM_ACB; i++)
    {
        gps = ACA_EPH[i * 4] + GPS_S - 19 - 32.184;
        mjd = gps / 86400.0 + T0 - 2400000.5;
        if ( mjd < mjdm)
        {
            mjd0 = 52532;
            bias_x_B = - 0.5647 
                       - 7.788e-5 * ( mjd - mjd0 ) 
                       + 2.4E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   7.5101 
                       + 7.495E-3 * ( mjd - mjd0 ) 
                       - 9.6E-6   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.8602 
                       + 1.399E-4 * ( mjd - mjd0 ) 
                       + 2.5E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }
        else
        {
            mjd0 = 53736;
            bias_x_B = - 0.6049 
                       - 1.982E-5 * ( mjd - mjd0 ) 
                       + 3.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
            bias_y_B =   10.6860 
                       + 1.159E-3 * ( mjd - mjd0 ) 
                       - 4.3E-7   * ( mjd - mjd0 ) * ( mjd - mjd0 );                       
            bias_z_B = - 0.7901 
                       + 4.783E-5 * ( mjd - mjd0 ) 
                       - 6.5E-9   * ( mjd - mjd0 ) * ( mjd - mjd0 );
        }

        ACB_EPH[i * 4 + 1] = bias_x_B * 1e-6 + scale_x_B * ACB_EPH[i * 4 + 1];
        ACB_EPH[i * 4 + 2] = bias_y_B * 1e-6 + scale_y_B * ACB_EPH[i * 4 + 2];
        ACB_EPH[i * 4 + 3] = bias_z_B * 1e-6 + scale_z_B * ACB_EPH[i * 4 + 3];

    }
   
    return 0;

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





double modvect (double *v)
{
    return  sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


double dotvect (double *v1, double *v2)
{
    return  v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


void crsvect (double *v1, double *v2, double *v)
{
    v[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v[2] = v1[0] * v2[1] - v1[1] * v2[0]; 
}



void cip2gcrf (int i, double *vt, double *vc)
{
    double lps, tt, jd_tdb, v3[3], v4[3];

//    lps = getlps (jd + utc/86400.0);
//    tt = utc + (lps + 32.184);  
    jd_tdb = info[i].jdt;
    cel_pole (jd_tdb, 2, info[i].dx * 1e3, info[i].dy * 1e3);
    nutation (jd_tdb,-1,1,vt, v3);
    precession (jd_tdb,v3,T0, v4);
    frame_tie (v4,-1, vc);

}




void tod2gcrf (int i, double *vt, double *vc)
{
    double lps, tt, jd_tdb, v3[3], v4[3], tepn[9];
//    i = (int)(NDATA / 2);
        iau_pn (info[i].jdt, tepn, 3);        //IAU inertial to J2000 inertial
        brmul (tepn,vt,3,3,1,vc);

}

void gcrf2tod (int i, double *vc, double *vt)
{
    double lps, tt, jd_tdb, v3[3], v4[3], tepn[9], te[9];

        iau_pn (info[i].jdt, tepn, 3);        //IAU inertial to J2000 inertial

        mt(tepn, 3, 3, te);
        brmul (te,vc,3,3,1,vt);
}

void iau_pns (double jd, double *te, int cent)
{
    double tes[9] = {0}, tepn[9] ={0}, tb[9], utc;

    double vx[3] = {1,0,0}, vy[3] = {0,1,0}, vz[3] = {0,0,1}, te2[9];
    int i;

        cent = cent +1;
        iau_s (jd, tes, cent);          //IAU fixed to IAU inertial
        iau_pn (jd, tepn, cent);        //IAU inertial to J2000 inertial
        brmul (tepn,tes,3,3,3,te);
    return;
}

void iau_s (double jd, double *tes, int cent)
{
    double d, str, cosst, sinst;

    d = jd - 2451545.0;

    switch (cent)       //sun0, mercury1, ..., pluto9 
    {
    case 3 : str = 190.147 + 360.9856235 * d; break;
    }

    cosst = cos (str * DEG2RAD);
    sinst = sin (str * DEG2RAD);

    tes[0] = cosst;
    tes[1] = -sinst;
    tes[2] = 0;
    tes[3] = sinst;
    tes[4] = cosst;
    tes[5] = 0;
    tes[6] = 0;
    tes[7] = 0;
    tes[8] = 1;
    return;
}



void iau_pn (double jd, double *tes, int cent)
{
    double ra0, dec0, jcent, cr, sr, cd, sd;

    jcent = jd - 2451545.0;
    jcent = (jcent) / 36525.0;

    switch (cent)       //sun0, mercury1, ..., pluto9 
    {
    case 3 :
        ra0  = 0.00 - 0.641 * jcent;
        dec0 = 90.0 - 0.557 * jcent;
        break;
    }

    cr = cos (ra0 * DEG2RAD);
    sr = sin (ra0 * DEG2RAD);
    cd = cos (dec0 * DEG2RAD);
    sd = sin (dec0 * DEG2RAD);

    tes[0] = -sr;
    tes[1] = -cr * sd;
    tes[2] = cr * cd;
    tes[3] = cr;
    tes[4] = -sr * sd;
    tes[5] = sr * cd;
    tes[6] = 0;
    tes[7] = cd;
    tes[8] = sd;
    return;
}





/*

double potdiff_tod (int i, double *p1, double *v1, double *p2, double *v2, 
                double *acc12)
{
    double term1, term2, term1m, term2m, term4a, term4am, termacc;
    double range, rate, rou0, modev, pa1[3], pa2[3], va1, va2,
        pxv[3], pxvxp[3], vxp[3], vxpxv[3], ep[3], ept[3], ev[3], evt[3];
    int n, iter;

    double pv1[3], pv2[3], pv1m[3], pv2m[3], we[3], wi[3], vr1, vr2, vr1m, vr2m;
       


    term1 = - dotvect( data1c[i].v12, data1c[i].v12) / 2; 
    term2 = dotvect( data1c[i].v2, data1c[i].v12);

    term1m =    dotvect( data1c[i].v2m, data1c[i].v2m) / 2; 
    term2m =  - dotvect( data1c[i].v1m, data1c[i].v1m) / 2;


    term4a = -ANGVEL * (data1c[i].p12[0] * data1c[i].v2[1] - data1c[i].p2[1] * data1c[i].v12[0] 
        + data1c[i].p1[0] * data1c[i].v12[1] - data1c[i].p12[1] * data1c[i].v1[0]);
    term4am = +ANGVEL * (data1c[i].p1[0] * data1c[i].v1[1] - data1c[i].p1[1] * data1c[i].v1[0]) 
              -ANGVEL * (data1c[i].p2[0] * data1c[i].v2[1] - data1c[i].p2[1] * data1c[i].v2[0]);


        we[0] = 0; we[1] = 0; we[2] = - ANGVEL;
        tod2gcrf (info[i].jd0, info[(int)(NDATA/2)].utc, we, wi);

        crsvect(data1c[i].p1, data1c[i].v1, pv1);
        crsvect(data1c[i].p2, data1c[i].v2, pv2);

        vr1 = dotvect(pv1, wi);
        vr2 = dotvect(pv2, wi);


        crsvect(data1c[i].p1m, data1c[i].v1m, pv1m);
        crsvect(data1c[i].p2m, data1c[i].v2m, pv2m);

        vr1m = dotvect(pv1m, wi);
        vr2m = dotvect(pv2m, wi);


        crsvect(data1c[i].p1, data1c[i].a1n, pa1);
        crsvect(data1c[i].p2, data1c[i].a2n, pa2);

        va1 = dotvect(pa1, wi);
        va2 = dotvect(pa2, wi);
    data1c[i].dvra = va2-va1;

    data1c[i].vk1 = term1;
    data1c[i].vk1m = term1m;
    data1c[i].vk2 = term2;
    data1c[i].vk2m = term2m;
    data1c[i].vr = vr2 - vr1;
    data1c[i].vrm = vr2m - vr1m;

    if ( kbrx[i].gps_time == 0) 
    {
        data1c[i].vk1m = term1;
        data1c[i].vk2m = term2;
    }


    return 0;
}

*/


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/



double xyz2aei(double ele[6], double pos[3], double vel[3])
{
    double a, e, omega, i, w, E, M, r, v, h, HV[3], n, 
        GM, radius, Pz, Qz;
    
    GM = 398600.44150E+09;
    radius = 6378136.3;
    r = modvect (pos);
    v = modvect (vel);

    a = 1.0 / (2.0 / r - v * v / GM);

    if(a <= 0)
    {    
        ele[0]=0,ele[1]=0,ele[2]=0,ele[3]=0,ele[4]=0,ele[5]=0;
//        printf("error: a <= 0 !\n");
        return 0;
    }
    
    if(a <= radius)
    {    
//        printf("warning: a <= radius !\n");
    }
    
    n = sqrt ( GM / (a*a*a) );

    crsvect (pos, vel, HV);
    h = modvect(HV);

    e = sqrt (1.0 - h * h / GM / a);

    i = acos (HV[2] / h);     //unit: rad

    omega = chosephase (HV[0] / h / sin(i), - HV[1] / h / sin(i));   //unit: rad

    E = chosephase ( dotvect(pos, vel) / (a * a * n * e), (1.0 - r / a) / e);  //unit: rad

    M = E - e * sin(E);        //unit: rad

    Pz = (cos(E) / r * pos[2] - sin(E) / n / a * vel[2]);
    Qz = (sin(E) / r / sqrt(1.0-e*e) * pos[2] + (cos(E) - e) / n / a / sqrt(1.0-e*e) * vel[2]);

    w = chosephase ( Pz / sin(i), Qz /sin(i));  //unit: rad

    ele[0] = a;
    ele[1] = e;
    ele[2] = i * RAD2DEG;
    ele[3] = omega * RAD2DEG;
    ele[4] = w * RAD2DEG;
    ele[5] = M * RAD2DEG;
                          
    return TWOPI / n;                              
//    return 5400;                              
}
                




int readgnv (char *infile_a, char *infile_b, int *n_a, int *n_b)
{
    FILE *fpGNV_a, *fpGNV_b;
    int n_gnva, n_gnvb, i, gps_i;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpGNV_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        ////getch();
        exit (0);
    }
    if ( (fpGNV_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        ////getch();
        exit (0);
    }


    n_gnva = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpGNV_a) ==NULL) break;
        n_gnva ++;
        sscanf (line, "%d", &gps_i);
//  printf ("%d \t %d \t %d \n",  GPS_S, gps_i, DT); 
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%s%s%lf%lf%lf%lf%lf%lf%lf", 
                &gnva[i].gps_time, &gnva[i].GRACE_id, &gnva[i].coord_ref, 
                &gnva[i].xpos, &gnva[i].ypos, &gnva[i].zpos,
                &gnva[i].xvel, &gnva[i].yvel, &gnva[i].zvel,
                &gnva[i].r);

            gnva[i].pos[0] = gnva[i].xpos;
            gnva[i].pos[1] = gnva[i].ypos;            
            gnva[i].pos[2] = gnva[i].zpos;
            gnva[i].vel[0] = gnva[i].xvel;
            gnva[i].vel[1] = gnva[i].yvel;            
            gnva[i].vel[2] = gnva[i].zvel;
            
            gnva[i].Tp = xyz2aei(gnva[i].ele, gnva[i].pos, gnva[i].vel);

//            printf ("%f\n",  gnva[i].Tp); 
        }

    }

    n_gnvb = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpGNV_b) ==NULL) break;
        n_gnvb ++;
        sscanf (line, "%d", &gps_i);
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%s%s%lf%lf%lf%lf%lf%lf%lf", 
                &gnvb[i].gps_time, &gnvb[i].GRACE_id, &gnvb[i].coord_ref, 
                &gnvb[i].xpos, &gnvb[i].ypos, &gnvb[i].zpos,
                &gnvb[i].xvel, &gnvb[i].yvel, &gnvb[i].zvel,
                &gnvb[i].r);

            gnvb[i].pos[0] = gnvb[i].xpos;
            gnvb[i].pos[1] = gnvb[i].ypos;
            gnvb[i].pos[2] = gnvb[i].zpos;
            gnvb[i].vel[0] = gnvb[i].xvel;
            gnvb[i].vel[1] = gnvb[i].yvel;
            gnvb[i].vel[2] = gnvb[i].zvel;
            
            gnvb[i].Tp = xyz2aei(gnvb[i].ele, gnvb[i].pos, gnvb[i].vel); 
//            printf ("%f\n",  gnvb[i].Tp); 
        }

    }

    fclose(fpGNV_a);
    fclose(fpGNV_b);
    *n_a = n_gnva;
    *n_b = n_gnvb;

    return 0;

}




int readkbr (char *infile, int *n)
{
    FILE *fpKBR;
    int n_kbr, i, gps_i;
    char line[MAXLINE];    


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpKBR = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open fpKBR file!\n");
        ////getch();
        exit (0);
    }

    n_kbr = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpKBR) ==NULL) break;
        n_kbr ++;
        sscanf (line, "%d", &gps_i);
        if ((gps_i - GPS_S) % DT == 0 && (gps_i - GPS_S) / DT < NDATA)
        {
            i = (gps_i - GPS_S) / DT;
            sscanf (line, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d%d", 
                &kbrx[i].gps_time, &kbrx[i].biased_range, &kbrx[i].range_rate, 
                &kbrx[i].range_accl, &kbrx[i].iono_corr, 
                &kbrx[i].lighttime_corr, &kbrx[i].lighttime_rate, &kbrx[i].lighttime_accl,
                &kbrx[i].ant_centr_corr, &kbrx[i].ant_centr_rate, &kbrx[i].ant_centr_accl,
                &kbrx[i].K_A_SNR, &kbrx[i].Ka_A_SNR, &kbrx[i].K_B_SNR, &kbrx[i].Ka_B_SNR,
                &kbrx[i].qualflg);
        }
        kbrx[i].range = kbrx[i].biased_range + kbrx[i].lighttime_corr + kbrx[i].ant_centr_corr;
        kbrx[i].rate  = kbrx[i].range_rate   + kbrx[i].lighttime_rate + kbrx[i].ant_centr_rate;
        kbrx[i].accl  = kbrx[i].range_accl   + kbrx[i].lighttime_accl + kbrx[i].ant_centr_accl;

    }


    fclose(fpKBR);

    *n = n_kbr;
    return 0;

}



int readacc (char *infile_a, char *infile_b)
{
    FILE *fpACC_a, *fpACC_b;
    int n_acca, n_accb, i, gps;
    char line[MAXLINE];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpACC_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpACC_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }


    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpACC_a) ==NULL) break;
//        sscanf (line, "%d%*s%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf",
            &gps, &ACA_EPH[i * 4 + 1], &ACA_EPH[i * 4 + 2], &ACA_EPH[i * 4 + 3]);
        ACA_EPH[i * 4] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_ACA = i;

    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpACC_b) ==NULL) break;
//        sscanf (line, "%d%*s%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf",
            &gps, &ACB_EPH[i * 4 + 1], &ACB_EPH[i * 4 + 2], &ACB_EPH[i * 4 + 3]);
        ACB_EPH[i * 4] = gps - GPS_S + 19 + 32.184;
        i++;
    }


    DIM_ACB = i;

    fclose(fpACC_a);
    fclose(fpACC_b);


    return 0;

}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int readsca (char *infile_a, char *infile_b)
{
    FILE *fpSCA_a, *fpSCA_b;
    int n_scaa, n_scab, i, gps;
    char line[MAXLINE];    


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    if ( (fpSCA_a = fopen (infile_a,"r")) == NULL)
    {
        printf ("Cannot open fpin_a file!\n");
        //getch();
        exit (0);
    }
    if ( (fpSCA_b = fopen (infile_b,"r")) == NULL)
    {
        printf ("Cannot open fpin_b file!\n");
        //getch();
        exit (0);
    }

    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpSCA_a) ==NULL) break;
//        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf%lf",
            &gps, &SCA_EPH[i * 5 + 1], &SCA_EPH[i * 5 + 2],
            &SCA_EPH[i * 5 + 3], &SCA_EPH[i * 5 + 4]);
        SCA_EPH[i * 5] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_SCA = i;


    i = 0;
    while (1)
    {
        if (fgets(line,MAXLINE, fpSCA_b) ==NULL) break;
//        sscanf (line, "%d%*s%*d%lf%lf%lf%lf",
        sscanf (line, "%d%lf%lf%lf%lf",
            &gps, &SCB_EPH[i * 5 + 1], &SCB_EPH[i * 5 + 2],
            &SCB_EPH[i * 5 + 3], &SCB_EPH[i * 5 + 4]);
        SCB_EPH[i * 5] = gps - GPS_S + 19 + 32.184;
        i++;
    }

    DIM_SCB = i;


    fclose(fpSCA_a);
    fclose(fpSCA_b);


    return 0;
}







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brmul ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void brmul (double *a, double *b, int m,int n, int k,double *c)
{ 
    int i, j, l, u;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= k - 1; j++)
        {
            u = i * k + j; 
            c[u] = 0.0;
            for (l = 0; l <= n - 1; l++)
                c[u] = c[u] + a[i * n + l] * b[l * k + j];
        }
    }
    return;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* mt ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void mt (double *a, int m, int n, double *b)
{
    int i, j;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= n - 1; j++)
            b[j * m + i] = a[i * n + j];
    }
    return;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double opengrav (char file_grv[2][200], double *coef, double *gma, 
        int nmax, int mmax, int head)
{
    FILE *fp_grv;
    double value, c,s;
    int n,m, l, ind;
    char string[MAXLINE], name[20];

    if ((fp_grv = fopen (file_grv[0],"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
        exit (0);
    }

    if (head == 1)
    {
        while (1)
        {
            if (fgets (string, MAXLINE, fp_grv) == NULL) break;
            sscanf (string, "%s%lf", name, &value); 
            if (strcmp (name,"Gm") ==0) 
            {
                gma[0] = value;
            }
            if (strcmp (name,"RefDistance") ==0)
            {
                gma[1] = value;
            }
            if (strcmp (name,"BEGIN") ==0)  
                break;
        }
    }
    else 
    {
        gma[0] = 398600.44150E+09;
        gma[1] = 6378136.3;
    }

//    coef[0] = 1;
    coef[0] = 0;
    while (1)
    {
        if (fgets (string, MAXLINE, fp_grv) == NULL) break;
//        sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);    


        if (strlen(file_grv[1])==0)
        {
            sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);  
        }
        else 
        {
            sscanf (string, "%s", name);    
            if (strcmp (name,file_grv[1]) ==0)  
            {
                sscanf (string, "%*s%d%d%lf%lf", &n, &m, &c, &s);
//                printf ("n = %d m = %d c = %e s = %e\n", n, m, c, s);
            }
        }


//        if (n > nmax || n < 0)
        if (n > nmax || n < 2 || m > mmax)   // permanently exclude degree 1 @7/24/2012
            continue;
        else if (m == 0)
        {
            coef[n] = c;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = c;
            coef[ind + n - m + l] = s;
        }
    }

    fclose(fp_grv);
    return 0;
}



/*
double pointmass (int num, double *p1, double *p2, double gm, double a, double *ap12)
{
    double r1, r2, ap1[3], ap2[3];  
    int n;

       
    r1 = sqrt (p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);

    r2 = sqrt (p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);


    for (n = 0; n < 3; n ++)
    {
        ap1[n] = - gm / (r1*r1*r1) * p1[n];
        ap2[n] = - gm / (r2*r2*r2) * p2[n];
        ap12[n] = ap2[n] - ap1[n];
    }

    return 0;

}
*/



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2acc (int num, double *llr, double *cs, double gm, double a, int nmax, 
               double *v, double *dvdt, double *acc)
{
    int n, m, k, l, ind;

    double sinf, cosf, sinlon, coslon, sincolat, coscolat, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, accn[3], c_ei[9], c_en[9], c_in[9],
        *pt, *ptt, lat, lon, r, vi, dvdr, dvdcolat, dvdlon, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);
    sinlon = sin(lon * DEG2RAD);
    coslon = cos(lon * DEG2RAD);
    sincolat = cosf;
    coscolat = sinf;

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    ptt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    for (m = 0; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = coscolat; vi = 0; dvdlon = 0; dvdcolat = 0; dvdr = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
//        lgdr(t, nmax, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
                ptt[k] = aprn[n] * pbar1[k];
                vi = vi + pt[k] * cs[k];
            
//                if (n>=2)
                {
                    dvdr = dvdr + (n+1) * pt[k] * cs[k];
                    dvdcolat = dvdcolat + ptt[k] * cs[k];
                }
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
                ptt[ind + n - m] = aprn[n] * pbar1[k] * cosml[m];
                ptt[ind + n - m + l] = aprn[n] * pbar1[k] * sinml[m];
                vi = vi + pt[ind + n - m] * cs[ind + n - m];
                vi = vi + pt[ind + n - m + l] * cs[ind + n - m + l];

                dvdcolat = dvdcolat + ptt[ind + n - m] * cs[ind + n - m];
                dvdcolat = dvdcolat + ptt[ind + n - m + l] * cs[ind + n - m + l];

                dvdlon = dvdlon - m * pt[ind + n - m + l] * cs[ind + n - m];
                dvdlon = dvdlon + m * pt[ind + n - m] * cs[ind + n - m + l];          

                dvdr = dvdr + (n+1) * pt[ind + n - m] * cs[ind + n - m];
                dvdr = dvdr + (n+1) * pt[ind + n - m + l] * cs[ind + n - m + l];
            }
        }
    }

//    dvdcolat = - dvdcolat * sincolat; //tmd!!
    dvdcolat = dvdcolat;
    dvdlon = + dvdlon;
    dvdr = - dvdr / r;


    accn[0] = - dvdcolat / r;
    accn[1] = + dvdlon / r / sincolat;
    accn[2] = - dvdr;


    c_en[0] = - sinf * coslon;      //from fixed to up-east-north system: rmat
    c_en[1] = - sinlon;
    c_en[2] = - cosf * coslon;
    c_en[3] = - sinf * sinlon;
    c_en[4] = coslon;
    c_en[5] = - cosf * sinlon;
    c_en[6] = cosf;
    c_en[7] = 0;
    c_en[8] = - sinf;


//    mt(info[num].c_ie, 3, 3, info[num].c_ei);
    brmul (info[num].c_ei, c_en, 3, 3, 3, c_in);  //inertial to fixed matrix gmat = rmat*tbt

    brmul(c_in, accn, 3, 3, 1, acc);  //from fixed acc to inertial acc


    *v = vi;
//    *v = gm/r + gm/r * a/r * a/r * (sqrt(5.0) * (1.5 * t * t - 0.5)) * cs[2];
    *dvdt = - ANGVEL * dvdlon;
//    if (lon < 180) *dvdt = - ANGVEL * dvdlon - 3.08e-12 * dvdcolat;
//    if (lon >=180) *dvdt = - ANGVEL * dvdlon + 3.08e-12 * dvdcolat;

    free (pbar);
    free (pbar1);
    free (pbar2);
    free (pt);
    free (ptt);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2vdvdt (double *llr, double *cs, double gm, double a, int NMAX, double *v, double *dvdt, double *pt)
{
    int n, m, k, l, ind;

    double sinf, cosf, *cosml, *sinml, *aprn, *pbar, lat, lon, r, 
        vi, dvdti, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( NMAX + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( NMAX + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( NMAX + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( NMAX + 1, sizeof(double));

    for (m = 0; m <= NMAX; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= NMAX; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = sinf; vi = 0; dvdti = 0;
    for (m = 0; m <= NMAX; m ++)
    {
        l = NMAX - m + 1;
//        lgdr2(t, NMAX, m, pbar, pbar1, pbar2);
        lgdr(t, NMAX, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
                vi = vi + pt[k] * cs[k];
            }
            else
            {
                ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
                vi = vi + pt[ind + n - m] * cs[ind + n - m];
                vi = vi + pt[ind + n - m + l] * cs[ind + n - m + l];
                dvdti = dvdti + m * pt[ind + n - m + l] * cs[ind + n - m];
                dvdti = dvdti - m * pt[ind + n - m] * cs[ind + n - m + l];
            }
        }
    }



//! ATPA
//  gpti = 0;
//    for(k = 0; k < (NMAX + 1) * (NMAX + 1); k++)
//    {
//        gpti = gpti + pt[k] * cs[k];
//        pt[k] = pnmc[k];
//    }

    *v = vi;
    *dvdt = dvdti * ANGVEL;

    free (pbar);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cspt2gp (double *pt, double *cs, int NMAX, double *gpt)
{
    int k;
    double gpti;

    gpti = 0;
    for(k = 0; k < (NMAX + 1) * (NMAX + 1); k++)
    {
        gpti = gpti + pt[k] * cs[k];
    }

    *gpt = gpti; 

    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/














/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2pt (double *llr, double *cs, double gm, double a, int NMAX, double *gpt, double *pt)
{
    int n, m, k, l, ind;

    double sinf, cosf, *cosml, *sinml, *aprn, *pbar, lat, lon, r, 
        gpti, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( NMAX + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( NMAX + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( NMAX + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( NMAX + 1, sizeof(double));

    for (m = 0; m <= NMAX; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= NMAX; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = sinf;
    for (m = 0; m <= NMAX; m ++)
    {
        l = NMAX - m + 1;
//        lgdr2(t, NMAX, m, pbar, pbar1, pbar2);
        lgdr(t, NMAX, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
            }
            else
            {
                ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
            }
        }
    }



//! ATPA
    gpti = 0;
    for(k = 0; k < (NMAX + 1) * (NMAX + 1); k++)
    {
        gpti = gpti + pt[k] * cs[k];
//        printf ("k = %d pt = %e\t cs = %e\t gp = %e \n", k, pt[k], cs[k], pt[k] * cs[k]);
//        pt[k] = pnmc[k];
    }


    *gpt = gpti; 


//    *gpt = gm/r + gm/r * a/r * a/r * (sqrt(5.0) * (1.5 * t * t - 0.5)) * cs[2];
    
//    printf ("cs[2] = %e\n", cs[2]);


    free (pbar);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double lgdr(double t, int nmax, int m, double *pbar)

/*
! THIS CALCULATES THE FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M),
! MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
! THIS RETURNS ALL Pn,m, P'n,m, AND P''n,m (m=<n<=Nmax).
! THE RECURSION FORMULAR FOR THE FUNCTION ITSELF IS GIVEN IN JEKELI(1996).
! THE RECURSION FORMULAR FOR THE 1ST DERIVATIVE IS GIVEN IN TSCHERNING, ET AL(1983).
! THE FORMULAR FOR THE 2ND DERIVATIVE IS FROM THE ASSOCIATE LEGENDRE EQUATION.
! NOTE : EQUATIONS GIVEN IN TSCHERNING, ET AL(1983) HAVE ERRATA.
!
! S.C. Han, 1/24/01 (MODIFIED FOR CRAY T94 2/13/01)
!
*/
{
    int i;
//REAL*8 :: PBAR(NMAX-M+1),PBAR1(NMAX-M+1),PBAR2(NMAX-M+1),T,P00,P11,C,D
    double p00, p11, c, d;
//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION
//! Pm,m : JEKEIL (A.3c) & (A.3d) , P'm,m : TSCHERNING (7)

    p00 = 1.0; 
    p11 = sqrt (3.0*(1.0-t*t));
    if (m>=1)
    {
        pbar[0] = p11; 

        for (i = 2; i <= m; i++)
        {
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
    }

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0)) 
            / ((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

    return 0;
}

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2)

/*
! THIS CALCULATES THE FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M),
! MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
! THIS RETURNS ALL Pn,m, P'n,m, AND P''n,m (m=<n<=Nmax).
! THE RECURSION FORMULAR FOR THE FUNCTION ITSELF IS GIVEN IN JEKELI(1996).
! THE RECURSION FORMULAR FOR THE 1ST DERIVATIVE IS GIVEN IN TSCHERNING, ET AL(1983).
! THE FORMULAR FOR THE 2ND DERIVATIVE IS FROM THE ASSOCIATE LEGENDRE EQUATION.
! NOTE : EQUATIONS GIVEN IN TSCHERNING, ET AL(1983) HAVE ERRATA.
!
! S.C. Han, 1/24/01 (MODIFIED FOR CRAY T94 2/13/01)
!
*/
{
    int i;
//REAL*8 :: PBAR(NMAX-M+1),PBAR1(NMAX-M+1),PBAR2(NMAX-M+1),T,P00,P11,C,D
    double p00, p11, c, d;
//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION
//! Pm,m : JEKEIL (A.3c) & (A.3d) , P'm,m : TSCHERNING (7)

    p00 = 1.0; 
    p11 = sqrt (3.0*(1.0-t*t));
    if (m>=1)
    {
        pbar[0] = p11; 
        pbar1[0] = sqrt(3.0) * t;

        for (i = 2; i <= m; i++)
        {
            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]);
//            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]/(-sqrt(1.0-t*t)));
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
        pbar1[0]=0.0;
    }

//    ! Pm+1,m : JEKEIL (A.3b)

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

//  ! Pn,m (n>=m+2) : JEKEIL (A.3a)

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0))/((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

//  ! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 1ST DERIVATIVE
//  ! P'n,m (n>=m+1) : TSCHERNING (8)
    for (i=2; i<=nmax-m+1; i++)
    {
        c = 1.0/sqrt(1.0-t*t)*t*(m+i-1);
        d = 1.0/sqrt(1.0-t*t)*sqrt((((m+i-1)*(m+i-1)-m*m)*(2.0*(m+i-1)+1.0))/(2.0*(m+i-1)-1.0));
//!! found it different from TSCHERNING (8),dcl-2010-2-14
//!! Jianbin confirms code is correct, dcl-2010-2-15
//!!      D=1D0/SQRT(1D0-T**2)/SQRT((((M+I-1)**2-M**2)*(2D0*(M+I-1)+1D0))/(2D0*(M+I-1)-1D0))
        pbar1[i-1] = c * pbar[i-1] - d * pbar[i-2];
    }

//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 2ND DERIVATIVE
//! P''n,m (n>=m) : ASSOCIATE LEGENDRE EQUATION (2ND ORDER DIFFERENTIAL EQN.)

    for (i=1;i<=nmax-m+1;i++)
    {
        pbar2[i-1] = (-t/sqrt(1.0-t*t)) * pbar1[i-1] 
            - ((m+i-1)*(m+i)-m*m/(1.0-t*t)) * pbar[i-1];
    }
    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/









  int bssgj(double a[], int n)
  { int i,j,k,m;
    double w,g,*b;
    b=malloc(n*sizeof(double));
    for (k=0; k<=n-1; k++)
      { w=a[0];
        if (fabs(w)+1.0==1.0)
          { 
//            free(b); 
            printf("fail\n"); 
//            return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++)
          { g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
              a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
          }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
          a[(n-1)*n+i-1]=b[i];
      }
    for (i=0; i<=n-2; i++)
    for (j=i+1; j<=n-1; j++)
    free(b);
    return(2);
  }




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brinv ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int brinv (double *a,int n)
{ 
    int *is, *js, i, j, k, l, u, v;
    double d, p;

    is = (int *)malloc (n * sizeof(int));
    js = (int *)malloc (n * sizeof(int));
    for (k = 0; k <= n - 1; k++)
    { 
        d=0.0;

        for (i = k; i <= n - 1; i++)
        {
            for (j = k; j <= n - 1; j++)
            { 
                l = i * n + j; 
                p = fabs (a[l]);
                if (p > d) 
                { 
                    d = p; 
                    is[k] = i; 
                    js[k] = j;
                }
            }
        }

        if (d + 1.0 == 1.0)
        { 
//            free (is); 
//            free (js); 
//            printf ("err**not inv\n");
//            return(0);
        }

        if (is[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = is[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        if (js[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + js[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        l = k * n + k;
        a[l] = 1.0 / a[l];
        
        for (j = 0; j <= n - 1; j++)
        {
            if (j != k)
            { 
                u = k * n + j; 
                a[u] = a[u] * a[l];
            }
        }
        
        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
                for (j = 0; j <= n - 1; j++)
                    if (j != k)
                    { 
                        u = i * n + j;
                        a[u] = a[u] - a[i * n + k] * a[k * n + j];
                    }
        }

        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
            { 
                u = i * n + k; 
                a[u] = - a[u] * a[l];
            }
        }
    }
    
    for (k = n - 1; k >= 0; k--)
    { 
        if (js[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = js[k] * n + j;
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
        
        if (is[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + is[k];
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
    }
    

    free(is); 
    free(js);
    return(1);
}


