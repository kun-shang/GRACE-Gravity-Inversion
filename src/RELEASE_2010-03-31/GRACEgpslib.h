/* $Id: GRACEgpslib.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

#define NAMELENGTH          8
#define NAMELENGTH_STRING   9

#define MAX_NSTATIONS     300
#define MAX_NSATELLITES    50
#define MAX_DATATYPE       15
#define MAX_PRN_NUMBER     32
#define MIN_PRN_NUMBER      1

#define CHIP              10.23e6
#define CHIP_CA           1.023e6
#define F1                (double)(154.0*CHIP)
#define F2                (double)(120.0*CHIP)
#define CSPEED            0.299792458e6
#define LAMBDA_1          (double)(CSPEED/F1)
#define LAMBDA_2          (double)(CSPEED/F2)
#define LAMBDA_CA         (double)(CSPEED/CHIP_CA)
#define GPS2UNIXTIME      946728000
#define TDPNAMELENGTH     16

/* grace specific timing issues */

#define USO_A_NOM_FREQ           4.832000e6
#define USO_A_RED_FREQ           4.832000e6
#define USO_B_NOM_FREQ           4.832099e6
#define USO_B_RED_FREQ           4.832099e6
#define K_FREQ_MULTIPLIER               5076
#define KA_FREQ_MULTIPLIER              6768
#define ONE_MEGAHERTZ_MULTIPLIER (4.0/19.0)

#define L1_NDX 0
#define L2_NDX 1
#define LC_NDX 2
#define P1_NDX 3
#define P2_NDX 4
#define PC_NDX 5
#define CA_NDX 6
#define DP_NDX 7
#define R2_NDX 8
#define AZ_NDX 9
#define EL_NDX 10

#define CA_RANGE_NDX  0
#define L1_RANGE_NDX  1
#define L2_RANGE_NDX  2
#define CA_PHASE_NDX  3
#define L1_PHASE_NDX  4
#define L2_PHASE_NDX  5
#define CA_SNR_NDX    6
#define L1_SNR_NDX    7
#define L2_SNR_NDX    8
#define CA_CHAN_NDX   9
#define L1_CHAN_NDX  10
#define L2_CHAN_NDX  11
#define K_PHASE_NDX  12
#define KA_PHASE_NDX 13
#define K_SNR_NDX    14
#define KA_SNR_NDX   15

#define L1_LOCK_NDX   0
#define L2_LOCK_NDX   1

typedef struct qm_obs_t {
        double RelTime;          /* relative time wrt to EpochTime    */
        int32_t   RcvId;            /* receiver index id                 */
        int32_t   TrnId;            /* transmitter index id              */
        int32_t   dtyp;             /* data type index                   */
        int32_t   mtyp;             /* measurement type index            */
                                 /*     rcv    trn                    */
                                 /* 1 = sat -> sat                    */
                                 /* 2 = sta -> sat                    */
                                 /* 3 = sat -> sta                    */
                                 /* 4 = sta -> sta                    */
        double qmbreak;          /* start time of phase bias          */
        float  sigma;            /* apriori sigma of observations     */
        int32_t   nobs;             /* number of observations            */
        double observations[MAX_DATATYPE]; /* array containing obs    */
        } qm_obs_t;

typedef struct qm_head_t {
        int32_t   Nstations;        /* number of stations in qm file     */
        int32_t   Nsatellites;      /* number of satellites in qm file   */
        int8_t  *StaName[MAX_NSTATIONS];  /* station names in qm file   */
        int8_t  *SatName[MAX_NSATELLITES];/* satellite names in qm file */
        double EpochTime;        /* epoch time of qmfile (sec2000)    */
        int32_t   nmax_dtype;       /* maximum number of data types      */                                   
        int32_t   nmax_obs;         /* maximum number of obs per record  */
        int32_t   n_meas;           /* number of measurements            */
        int32_t   sort_stat[5];     /* sort index array                  */
        } qm_head_t;

typedef struct Prn_t {
        int32_t flg;                 /* prn flag                          */ 
        int32_t flg_prev;            /* previous prn flag                 */ 
        double ttag;             /* Time tag                          */
        double ttag_last;        /* Previous Time Tag                 */
        double tu[3];            /* Pseudo Range (m) 1) CA 2) P1 3) P2*/
        double fz[3];            /* Phase (cycles)   1) CA 2) L1 3) L2*/
        double fz_prev[3];       /* Previous phase   1) CA 2) L1 3) L2*/
        int32_t snr[3];              /* SNR (V/V)        1) CA 2) L1 3) L2*/
        int32_t ss[3];               /* encoded SNR      1) CA 2) L1 3) L2*/
        int32_t nwrap;               /* nwrap counter (used for K and Ka) */
        } Prn_t;

typedef struct rnx_header_t{
        int8_t pgm[21];            /* program name                      */
        int8_t marker[61];         /* marker name                       */
        int8_t observer[21];       /* observer name                     */
        int8_t agency[41];         /* agency  name                      */
        int8_t recnum[21];         /* receiver number                   */
        int8_t rectype[21];        /* receiver type                     */
        int8_t recvers[21];        /* receiver version                  */
        int8_t antnum[21];         /* antenna number                    */
        int8_t anttype[21];        /* antenna type                      */
        double pos[3];           /* approximate position              */
        double ant[3];           /* antenna delta                     */
        double StartTime;        /* time of first observation         */
        double FinalTime;        /* time of last  observation         */
        double Interval;         /* data interval                     */
        } rnx_header_t;

typedef struct tdp_t {
        double time;             /* time tag (seconds)                */ 
        double apriori;          /* apriori value of tdp value        */ 
        double value;            /* value of tdp                      */
        double sigma;            /* sigma of value                    */
        int8_t name[17];           /* parameter name of value           */
        } tdp_t;

typedef struct pos_goa_t {
        int32_t   nfield;           /* number of fields in record        */ 
        int8_t   frame;            /* E or I = Earth fixed or Inertial  */ 
        int8_t   satname[NAMELENGTH_STRING];/* satellite name           */
        double time;             /* time tag (seconds)                */ 
        double xpos;             /* x coordinate (km)                 */ 
        double ypos;             /* y coordinate (km)                 */ 
        double zpos;             /* z coordinate (km)                 */ 
        double xvel;             /* x velocity   (km/sec)             */ 
        double yvel;             /* y velocity   (km/sec)             */ 
        double zvel;             /* z velocity   (km/sec)             */ 
        double xpos_sig;         /* sigma x coordinate (km)           */ 
        double ypos_sig;         /* sigma y coordinate (km)           */ 
        double zpos_sig;         /* sigma z coordinate (km)           */ 
        double xvel_sig;         /* sigma x velocity   (km/sec)       */ 
        double yvel_sig;         /* sigma y velocity   (km/sec)       */ 
        double zvel_sig;         /* sigma z velocity   (km/sec)       */ 
        } pos_goa_t;

void    ReadRecByte(FILE* src, int32_t *RecBytes);
void    WriteRecByte(FILE* dst, int32_t RecBytes);
boolean ReadQmHeader(FILE *src, qm_head_t *qmhead, int32_t verbose_flag);
boolean WriteQmHeader(FILE *src, qm_head_t *qmhead);
boolean ReadQmRecord(FILE* src, qm_head_t *qmhead, qm_obs_t *qmobs);
boolean WriteQmRecord(FILE* src, qm_head_t *qmhead, qm_obs_t *qmobs);
void    ComputeLC_PC(double f1, double f2, double L1, double L2, 
                     double L1sigma, double L2sigma, double *LC,
                     double *LCsigma);
void    ComputeLW(double f1, double f2, double L1, double L2, double *LW);
void    ComputePW(double f1, double f2, double P1, double P2, double *PW);
int32_t    prn2gpsh(int32_t prn, double T2000);
int32_t    gps2prnh(int32_t gps, double T2000);
int32_t    WriteRnxHeader(FILE *dst,int32_t do_snrs,int32_t do_ca, rnx_header_t *rnxhead);
int32_t    PurgeRnxBlock(FILE *dst, Prn_t *pP,int32_t Npoints, int32_t do_snrs,
                      int32_t do_ca);
int32_t    ReadTdp(FILE* src, tdp_t *TdpRec);
int32_t    WriteTdp(FILE* dst, tdp_t *TdpRec);
int32_t     get_words(int8_t *str, int8_t **words, int32_t nwords_max);
int32_t    ReadGoaPosRecord(FILE* src, pos_goa_t *posgoa);
int32_t    WriteGoaPosRecord(FILE* dst, pos_goa_t *posgoa);
int32_t     lagrange_int( double x, double *y, int32_t ntab, double *xt, double *yt, 
                      int32_t ndeg, int32_t compute_deriv, double *yd);
void    WideCode (double f1, double f2, double PW, double LW, double *Code);
void    WideIon (double f1, double f2, double L1, double L2, double *Ion);
int32_t    eci2quaternion(int32_t T2000_int, double T2000_frac ,quaternion *Q,int32_t norder, FILE *eci);
int32_t eci2grace_quaternion(int32_t T2000_int, double T2000_frac , double clkerrA,
                          double clkerrB,
                          quaternion *Q_A, quaternion *Q_B, int32_t norder_int,
                          FILE *eciA, FILE *eciB);
void LoadECIOrbitFile(double **pXorb, double **pYorb, double **pZorb,
                   double **pXdot, double **pYdot, double **pZdot,
                   double **pTimeOrb, double *pTimeZero, int32_t *pNorb, FILE *eci);
int32_t OrbitEciInt(int32_t seconds, double sec_frac, double xyz[3], double llh[3],
             double xyzdot[3], double ae, double flat, int32_t norder, FILE *eci);
int32_t OrbitEciInt_second(int32_t seconds, double sec_frac, double xyz[3], double llh[3],
             double xyzdot[3], double ae, double flat, int32_t norder, FILE *eci);
void xyz2llh(double xyz[3] , double ae, double flat, double llh[3]);
void llh2xyz(double xyz[3] , double ae, double flat, double llh[3]);
double VectorNorm(double *x, int32_t n);
double DotProduct(double *x, double *y, int32_t n );
void CrossProduct(double *x, double *y, int32_t n, double *cross);

