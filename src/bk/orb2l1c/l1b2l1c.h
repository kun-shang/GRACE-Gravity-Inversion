
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _L1B2L1C_H_
    #define _L1B2L1C_H_

    #ifndef __OMP_H__
        #include <omp.h>
    #endif

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

    #ifndef _NOVAS_
       #include "novas.h"
    #endif

    #ifndef _EPHMAN_
       #include "eph_manager.h"
    #endif

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    typedef struct GNV1C
    {
        int gps_time;
        char GRACE_id;
        char coord_ref;
        double xpos;
        double ypos;
        double zpos;
        double xvel;
        double yvel;
        double zvel;
        double pos[3];
        double vel[3];
        double lat;
        double lon;
        double r;
        double ele[6];
        double Tp;
    }GNV1C;



    typedef struct KBR1B
    {
        int gps_time;
        double biased_range;
        double range_rate;
        double range_accl;
        double iono_corr;
        double lighttime_corr;
        double lighttime_rate;
        double lighttime_accl;
        double ant_centr_corr;
        double ant_centr_rate;
        double ant_centr_accl;
        int K_A_SNR;
        int Ka_A_SNR;
        int K_B_SNR;
        int Ka_B_SNR;
        int qualflg;
        double range;
        double rate;
        double accl;
    }KBR1B;


    typedef struct ACC1B
    {
        int gps_time;
        char GRACE_id;
        double lin_accl_x;
        double lin_accl_y;
        double lin_accl_z;
        double ang_accl_x;
        double ang_accl_y;
        double ang_accl_z;
        double acl_x_res;
        double acl_y_res;
        double acl_z_res;
        int qualflg;
    }ACC1B;


    typedef struct SCA1B
    {
        int gps_time;
        char GRACE_id;
        int sca_id;
        double quatangle;
        double quaticoeff;
        double quatjcoeff;
        double quatkcoeff;
        double qual_rss;
        int qualflg;
    }SCA1B;



    typedef struct OTStruct
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
    }OTStruct;



    

    typedef struct InfStruct
    {
        int gps;
        int leaps;
        double jd0;
        double jdt;
        double mjd;        
        double gmst;
        double deltat;
        double utc;
        double ut1;
        double tt;
        double xp;
        double yp;
        double ut1_utc;
        double dx;
        double dy;
        double c_ie[9];
        double c_ei[9];
        double wi[3];
    }InfStruct;


    typedef struct DATAL1C
    {
        int error;  //0: normal; 1: no orbit; 2: no kbrr; 3: no kbra; 9: outlier
        double t;
        
        double vl1cflt;
        double glosflt;

        double vl1c;
        double fitacc;
        double vpwfit;
        double vpwfitstd;
        double vpwfitres;
        double vmvfit;
        double vmvfitstd;
        double vmvfitres;
        double vc;      // L1C potential difference 

        double vl1cm;
        double vcm;
        double vpwfitm;

        double vcsr;      // Official L2 (n = 0 ~ 60)
        double vjpl;      // Official L2 (n = 0 ~ 60)
        double vgfz;      // Official L2 (n = 0 ~ 60)
        double vggm;      // Official L2 (n = 0 ~ 60)
   
        double acsr[3];
        double ajpl[3];
        double agfz[3];
        double glos;
        double gpwfit;
        double gpwfitstd;
        double gpwfitres;
        double gmvfit;
        double gmvfitstd;
        double gmvfitres;
        double gcsr;
        double gjpl;
        double ggfz;
        
        double vcorr;
        double dvcorr;
        double vref;
        double gref;
//        double losg;
//        double losb;
//        double dvtg;

//        double vg; 
//        double dvg;
//        double ag[3];  


        double vk;      // vk = vk1 + vk2;

        double vr;      // rotation term
        double vr3;
        double vr12;
        double vrm;


        double vrfl;
        double arfl[3];   
        double grfl; 
        double vrfh;  
        double arfh[3];
        double grfh;
        double vb;  
        double ab[3];
        double gb;
        double vbl;  
        double abl[3];
        double gbl;
        double dvab;
        double vab; 
        double dvrb;
        double vrb;
        double dvtb;
        double vtb;        
    

        double af[3];  
        double dvaf;
        double vaf; 
        double dvrf;     
        double vrf;    

        double vn;
        double an[3];
        double dvan;
        double van;
        double dvrn;     
        double vrn;
        double dvtn;
        double vtn;

        double vs;
        double as[3];
        double dvas;
        double vas;
        double dvrs;     
        double vrs;
        double dvts;
        double vts;

        double vo;
        double ao[3];
        double dvao;
        double vao;
        double dvro;     
        double vro;
        double dvto;
        double vto;

        double vh;
        double ah[3];
        double dvah;
        double vah;
        double dvrh;     
        double vrh;
        double dvth;
        double vth;


        double va;
        double aa[3];
        double dvaa;
        double vaa;
        double dvra;     
        double vra;
        double dvta;
        double vta;

        double vp;
        double ap[3];
        double dvap;
        double vap;
        double dvrp;
        double vrp;
        double dvtp;
        double vtp;

        double ag[3];
        double dvag;
        double vag;
        double dvrg;
        double vrg;


        double rp[3];
        double llr[3];
        double pos;
        double rv[3];
        double vel;
        double ra[3];

        double wr[3];

        double range;
        double rate;
        double accl;
        double Tp1; 
        double Tp2; 
   
    }DATAL1C;


    double *ACA_EPH, *ACB_EPH, *SCA_EPH, *SCB_EPH, *AOD_EPH, *OT_EPH, *AT_EPH, *OPTM1, *OPTM2;
    int DIM_ACA, DIM_ACB, DIM_SCA, DIM_SCB;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    GNV1C *gnva, *gnvb;

    KBR1B *kbrx;

    ACC1B *acca, *accb;

    SCA1B *scaa, *scab;

    OTStruct *otfes;
    
    InfStruct *info;

    DATAL1C *sat1a, *sat2b, *sat12;

    int NDATA, DT, GPS_S, NFES, PERMT, STIDE, OTIDE, ATIDE, 
        MACC_BIAS, MACC_DTBS, MACC_SCAL, MACC_DTSL;
    double GMA[2], C20PERM;
    short int ACCURACY, INTERP;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
#ifdef _cplusplus
extern "C" 
{
#endif 

extern void DOODSN (double *jd_tdb, double *gmst, int *argn, int *ncon, 
        double *doodarg, double *ang);   
#define DOODSN doodsn_
#ifdef _cplusplus
}
#endif
*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    double estacc();
    void iau_pn (double jd, double *tes, int cent);
    void iau_pns (double jd, double *te, int cent);
    void iau_s (double jd, double *tes, int cent);
    double poletide (InfStruct *info, int nmax, double *coef);
    double openopt (char *f_pt, int nmax);

    void gcrf2tod (int i, double *vc, double *vt);

    void cip2gcrf (int i, double *vt, double *vc);
    double kbrvel2pos (int i);
    double kbrpos2vel (int i);

    double nbodyv (double *tjd, double *rp, double *vn, double *dvtn, double *an);

    double stidev (int num, double *llr1, double *vs, double *dvts, double *as);

    double otidev12 (int num, int nmax, double *llr1, double *llr2);
    double otidev (int num, int nmax, double *llr1, double *vs, double *dvts, double *as);

    double mtpole (double mjd, double xp, double yp, double *m1, double *m2);
    double normfct (int n, int m);
    double arg2theta (double jdt, double gmst, int argn[], double *ang);
//    double stpole (double mjd, double xp, double yp, double *c21p, double *s21p);
    double otpole (double mjd, double xp, double yp, double *c21p, double *s21p);

    double stfrqdep(double jdt, double gmst, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f);

    int openotcs_csr (char *infile);

    double otidecs_csr(InfStruct *info, int nmax, double *coef);

    int openotcs_fes (char *infile);


    double stidecs_Anelastic(InfStruct *info, int id_perm, double *stcs);


    double intydt(int num, double dt, double *y, double *val);

    void tod2gcrf (int i, double *vt, double *vc);


    double potdiff_tod (int i, double *p1, double *v1, double *p2, double *v2, 
                double *acc12);

//    double stidecs_Anelastic(int num, int id_perm, double *stcs);
//    double getvrm (int num, double *coef, int nmax, double *llr1, double *llr2, double *dvdt);
    double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr);

    void reltivint ();
    double reltivdiff (double vb1, double vb2, double *vrel);
    double mean (double * array, double N);
    double std_dev (double * array, double N);
    void fitres (int order_poly, int order_cpr, int offsetcyc);
//    void fitl1cdcl (int nply, int ncpr, int nplycpr, int ndt);
    double lsfl1c1d (double *xmax, double *x, double *y, int nmax, int n,
        double tpsec, double *ymaxflt, int nply, int ncpr, int nplycpr);

    void fitl1c_pw (int order_poly, int order_cpr);
    void fitlos_pw (int order_poly, int order_cpr);

    void fitl1c_mv (int nply, int ncpr, int nplycpr, int ndt);
    void fitlos_mv (int nply, int ncpr, int nplycpr, int ndt);
//    void fitres_new (int order_poly, int order_cpr);
    void fitres_new_overlap (int order_poly, int order_cpr, int overlap);

    double xyz2aei(double ele[6], double pos[3], double vel[3]);

    int lsf_cpr ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr);

    double lsf_cpr_new (double *xmax, double *x, double *y, int nmax, int n, 
        double tpsec, double *yflt, int order_poly, int order_cpr);

    int lsf_cpr_day ( double *x, double *y, int n, double tpsec, double *yflt, int order_poly, int order_cpr);

    double alignrrr(int i);
    double getwinbp(int n, double *win, double kc, double kh, int m);
    double getwinbp0(int n, double *win, double kc, int m);

    double getwinhp(int n, double *win, double kc, int m);

    int dft(long int length, double real_sample[], double imag_sample[]);
    int inverse_dft(long int length, double real_sample[], double imag_sample[]);
    double daydft (double *data, double *fft, int num,  double dt, double tc);

    double dayfft (double *data, double *fft, int num,  double fs);
    void kkfft(double pr[], double pi[], int n, int k, double fr[], 
		   double fi[], int l, int il);

    void realft(double data[], int n, int isign);
    void four1(double data[], int nn, int isign);

    double resfft(double *res, double *resflt, int num);

    double k5filter (double x[], double y[], int size, double fln, 
                 double fhn, double fs, int order);

    void firwin(int n, double fln, double fhn,double h[]);

    double hanningwin(int n,int i);

    double kaiserwin(int n,int i);

    double modvect (double *v);
    double dotvect (double *v1, double *v2);
    void crsvect (double *v1, double *v2, double *v);


//    double fnlgdr(double t, int n, int m);

    double getinfo(char *infile);

//    int openotcs (char *infile, int *n);

//    double otidecs (int num, int nmax, int id_long, double *stcs);

    double nbodypt (double *tjd, double *p, double *pt, double *a);
//    double nbodydiff (int num, double *p1, double *p2, double *npd, 
//        double *acc12, double *dvdt);
//    double nbodydiffdv (int num, double *p1, double *p2, double *npd, 
//        double *acc12, double *dvdt);

//    double stidecs (int num, int id_perm, double *stcs);
//    double stfrqdep (int num, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f);
//    double soliddiff (int num, double *llr1, double *llr2, 
//            double *spd, double *acc12, double *dvdt);
//    double soliddiffdv (int num, double *llr1, double *llr2, 
//            double *spd, double *acc12, double *dvdt);

//    double oceandiff (int num, double *p1, double *p2, double *opd, double *a);
    double atmocdiff (int num, double *p1, double *p2, double *apd, double *a);

    double cspt2gp (double *pt, double *cs, int NMAX, double *gpt);
    double cs2pt (double *llr, double *cs, double gm, double a, int NMAX, double *gpt, double *pt);
    double cs2vdvdt (double *llr, double *cs, double gm, double a, int NMAX, double *v, double *dvdt, double *pt);
    double cs2acc (int num, double *llr, double *cs, double gm, double a, int nmax, 
                   double *v, double *dvdt, double *acc);

//    double pointmass (int num, double *p1, double *p2, double gm, double a, double *ap12);

    double opengrav (char file_grv[2][200], double *coef, double *gma, 
        int nmax, int mmax, int head);

    int brinv (double *a,int n);

    double lgdr(double t, int nmax, int m, double *pbar);

    double lgdr2(double t, int nmax, int m, double *pbar, 
        double *pbar1, double *pbar2);

    double disse_a (double *v1, int i, double *ef1i, double *f1);
    double disse_b (double *v1, int i, double *ef1i, double *f1);


    int readkbr (char *infile, int *n);
    int readgnv (char *infile_a, char *infile_b, int *n_a, int *n_b);
    int readacc (char *infile_a, char *infile_b);
    int readsca (char *infile_a, char *infile_b);
    double cal_acc_01(void);
    double cal_acc_00(void);
    int calbias(char *infile_a, char *infile_b);
    int calbiaseph(char *infile_a, char *infile_b);
    double accia (double *tjd, int label, double *acc);
    double cal2_acc(void);
    double bias_acc(int order);
    int lsf_poly ( double *x, double *y, int n, double *a, int k);
    void mtrans(double a[], int m, int n, double b[]);
  
    int bssgj(double a[], int n);

    double quat2mat_i2s (double *qvec, double *mat);

    void brmul (double *a, double *b, int m,int n, int k,double *c);

    void mt (double *a, int m, int n, double *b);

    short int getlps (double jd);

    void xyz2llh (double *vt, double *llh);

    void llh2xyz (double *llh, double *vt);

    double lagrange (double *y, int dim_y, int dim_x, double t, double *z);
    double lgr_order (double *y, int dim_y, int dim_x, double t, double *z, int order);
    double openaod (char *file_aod, int nmax, int mmax);
    double lagrangelow (double *y, int dim_y, int dim_x, double t, double *z);
    double accgr (double *tjd, double *p, double *v, double *fgr);


    double icrf2itrf(int num, double *vc, double *vt);

    double chosephase (double sinvalue, double cosvalue);



#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

