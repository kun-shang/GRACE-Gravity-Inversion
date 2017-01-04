#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Version Information: $Id: Cal_inc.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

typedef double Real;

typedef struct {
  int32_t n;
  double *x;
} Vec;

typedef struct {
  int32_t nr, nc;
  double **x;
} Mat;

typedef struct {
  int8_t *name;
  FILE *fp;
} File_t;

typedef struct {
  int32_t n;
  double ep;
  double dt;
  double *p[3];
  double *v[3];
  int8_t *name;
} Eci_t;

#define MAXLINE 2048
#define MAXDEG 16
#define MAXN 2000
#define MAXDIM 2000
#define NFIT 9

#define elsif else if
#define streq(A,B) (!strcmp(A,B))

#define loop(I,N) for(I=0;I<N;I++)
#define loop3(I) loop(I,3)
#define loop33(I,J) loop3(I) loop3(J)

#define rnd(A) floor(A + 0.5)
#define ACC (1.0e-12)


/* global variables */

#ifdef MAIN
#define EXT
#define EXTINIT(A,B) A = B
#else
#define EXT extern
#define EXTINIT(A,B) extern A
#endif

EXT Vec OldSoln, OldSig, Soln2, Sig2;

EXT Real PI, D2R, R2D, a_K, a_Ka;
EXT Real xhat[3], yhat[3], zhat[3];
EXT Vec doff;
EXT Real xdum[16];
EXT int32_t idum[16];
EXT int8_t cdum[16];

EXT int8_t *pname;
EXT int8_t *GRACE_id;
EXT File_t File_log, File_acc, File_mag, File_sca, File_kbr, File_eci, File_eci2,
  File_cg, File_old, File_new, File_soln, File_param, File_oldsoln, File_soln2;
EXTINIT ( FILE *fp_log, stdout );

EXT Eci_t eci, eci2;

EXTINIT ( int8_t *xyz, "xyz" );

EXTINIT ( Real solx, 0 );  /* meters */
EXTINIT ( Real soly, 0 );  /* meters */
EXTINIT ( Real solz, 0 );  /* meters */

EXTINIT ( Real sig_ap, 1e3 );  /* meters */
EXTINIT ( Real sig_tight, 1e-9 );  /* meters */
EXT Real t1, t2, period;
EXT int32_t dummy_set;
EXTINIT ( int32_t t1_set, 0 );
EXTINIT ( int32_t t2_set, 0 );
EXTINIT ( int32_t period_set, 0 );
EXTINIT ( int32_t do_real, 0 );
EXTINIT ( int32_t do_sim, 0 );
EXTINIT ( int32_t do_hola, 0 );
EXTINIT ( int32_t GRACE_id_set, 0 );
EXTINIT ( int32_t do_write, 0 );

EXT Real Mass, Lhorn;
EXT Mat MoI, iMoI;
EXT Mat Normal_AA, Cov;
EXT Vec Normal_Ay, Soln, Sig, sigA, wtA;
EXT Mat Amat, iAmat;
EXT Vec Yvec, Bvec;
EXT int8_t Line[1024];

EXTINIT ( Real Xfac, 1.0 );
EXTINIT ( Real Yfac, 1.0 );
EXTINIT ( Real Zfac, 1.0 );

EXTINIT ( int32_t acc_n, 0 );
EXTINIT ( int32_t mag_n, 0 );
EXTINIT ( int32_t sca_n, 0 );
EXTINIT ( int32_t kbr_n, 0 );
EXTINIT ( int32_t pos_n, 0 );
EXT Real acc_t[MAXN], mag_t[MAXN], sca_t[MAXN], kbr_t[MAXN], pos_t[MAXN], kbr_X[MAXN];
EXT Mat acc_A, mag_al, sca_Q, pos_R;

/* prototypes */

/* from Cal_CG.c */
int32_t read_normal_eqns();
int32_t write_normal_eqns();

/* from cmdline.c */
int32_t Get_Real_Arg(Real *p, int32_t *flg);
int32_t Get_String_Arg(int8_t **p, int32_t *flg);
int32_t File_Open_Arg(File_t *f, const int8_t *mode);
int32_t No_Arg();
int32_t read_cmdline(int32_t argc, int8_t **argv);
int32_t do_usage();

/* from rx.c */
int32_t bomb(int8_t *error_text);
int32_t mat_inv(Mat a, Mat ai);
int32_t mat_trans(Mat a);
int32_t my_ludcmp(Real **a, int32_t n, int32_t *indx, Real *d);
int32_t my_lubksb(Real **a, int32_t n, int32_t *indx, Real *b);
Real pythag(Real a, Real b);
int32_t tred2(Real **a, int32_t n, Real *d, Real *e);
int32_t tqli(Real *d, Real *e, int32_t n, Real **z);
Real fsign(Real a, Real b);

/* from utils.c */
Vec Vec_alloc(int32_t n);
Mat Mat_alloc(int32_t nrows, int32_t ncols);
int32_t fit_pd(int32_t n, Real tep, Real *ts, Real *xs, Real *sigs, Real pd, int32_t *flags, Real *as, Real *amp, Real *res);
int32_t set_v(Vec vn, Vec v);
int32_t m_v(Vec mv, Mat m, Vec v);
int32_t v_m(Vec vm, Vec v, Mat m);
int32_t m_m(Mat mab, Mat ma, Mat mb);
int32_t rot_x(Mat mat, Real th);
int32_t rot_y(Mat mat, Real th);
int32_t rot_z(Mat mat, Real th);
int32_t rot_gen(Mat mat, Real th, Vec ax, Vec u, Vec v);
Real atan3(Real x, Real y);
Real theta(Real x, Real y);
Real Vnorm(Vec vhat, Vec v);
Real Vdot(Vec a, Vec b);
int32_t cross(Vec axb, Vec a, Vec b);
int32_t diff(Vec ab, Vec a, Vec b);
Real costh(Vec a, Vec b);
int32_t ahola(int32_t n);
int32_t neg(Vec a);
int32_t getq(Vec q, Vec x, Vec y, Vec z, Vec u, Vec v, Vec w);
Real getrand();

