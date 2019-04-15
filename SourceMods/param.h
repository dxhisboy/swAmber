#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define ToGBFloat(x) (x)
//#define ToGBFloat(x) SNGL(x)
//#define GBFloat float
#define GBFloat double
#define PI 3.1415926535897932384626433832795
#define C_GBL 310
#define C_NPTRA 301

#define CACHE_HBIT 5
#define CACHE_SBIT 3
#define CACHE_LINESIZE (1 << CACHE_SBIT)
#define CACHE_LINECNT (1 << CACHE_HBIT)
#define CACHE_MMASK (CACHE_LINESIZE - 1)
#define CACHE_LMASK (CACHE_LINECNT - 1)

#define W_CACHE_SIZE_T 5
#define W_CACHE_SIZE_P 3
#define W_C_S_T W_CACHE_SIZE_T 
#define W_C_S_P W_CACHE_SIZE_P
#define T_CPE_SIZE (1 << W_C_S_T)
#define P_CPE_SIZE (1 << W_C_S_P)
#define T_C_S  T_CPE_SIZE
#define P_C_S  P_CPE_SIZE
#define TU_SIZE 1
#define ZNUM 2048

typedef struct listdata_rec_c
{
  int offset;
  int cnt;
}listdata_rec_c;

typedef struct cit_tbl_rec_c
{
    int img_lo;
    int img_hi;
}cit_tbl_rec_c;


typedef struct pme_mask_param
{
  int img_lo;
  int img_hi;
  int my_bkt_lo; 
  int my_bkt_hi;
  int cit_tbl_x_dim;
  int cit_tbl_y_dim;
  int cit_tbl_z_dim;
  int ndbkt;
  int *dbkt;
  listdata_rec_c *atm_maskdata;
  int *atm_mask;
  listdata_rec_c *img_maskdata;
  int *img_mask;
  int *gbl_atm_img_map;
  int *gbl_img_atm_map;
  cit_tbl_rec_c *cit_save;
  int ncit_save;
  int myrank;
  //int atm_cnt;

}pme_mask_param;

typedef struct atm_lst_rec_c 
{
  int idx;
  int nxt;
}atm_lst_rec_c;


typedef struct cit_param_c
{
  double * fraction;
  int *ifrc_idx;
  int *f_idx;
  int atm_cnt;
  atm_lst_rec_c *atm_lst;
  int *atm_sorting_key;
  int *tail_idx_tbl;
  int *crd_idx_cnt;
  int *crd_idx_cnt_two;
  int *crd_idx_start;
  int *crd_idx_lst_tbl;
  int cit_tbl_x_dim;
  int cit_tbl_y_dim;
  int cit_tbl_z_dim;
  int numflat;
  int nxt_atm_lst_idx;
  double scale_fac_x;
  double scale_fac_y;
  double scale_fac_z;
  double sub_fac_x;       
  double sub_fac_y;        
  double sub_fac_z;
  int *hilbert3d4;
  int *hilbert3d8;
}cit_param_c;

typedef struct memset_param
{
  int cnt;
  double *array;
}memset_param;

typedef struct dihed_rec_c
{
  int atm_i;
  int atm_j;
  int atm_k;
  int atm_l;
  int parm_idx;
}dihed_rec_c;
typedef struct dihed_crd_c
{
  double crd_i[3];
  double crd_j[3];
  double crd_k[3];
  double crd_l[3];
  int atmi_map;
  int atmj_map;
  int atmk_map;
  int atml_map;
}dihed_crd_c;

//typedef struct dihed_crd_c
//{
//  double crd_i0, crd_i1, crd_i2;
//  double crd_j0, crd_j1, crd_j2;
//  double crd_k0, crd_k1, crd_k2;
//  double crd_l0, crd_l1, crd_l2;
//}dihed_crd_c;

typedef struct dihed_frc_c
{
  double frc_i[3];
  double frc_j[3];
  double frc_k[3];
  double frc_l[3];
}dihed_frc_c;

typedef struct dihed_atm_c
{
  double x[3];
  int atm_map;
  int padding;
}dihed_atm_c;

typedef struct dfrc_packed_c
{
  double frc[3];
  int atm_idx;
  int padding;
}dfrc_packed_c;


typedef struct dihed_param_c
{
  dihed_rec_c *dihed;
  dihed_crd_c *dihed_crd;
  double *x;
  double *frc;
  double *gbl_pk;
  double *gbl_pn;
  double *bar_cont;
  double *bar_lambda;
  double *gbl_gamc;
  double *gbl_gams;
  double *ti_signs;
  double *ti_ene;
  double *ti_ene_delta;
  double *ti_weights;
  int *gbl_ipn;
  int *ti_lst;
  int *ti_sc_lst;
  int *gbl_atm_owner_map;
  int dihed_cnt;
  int atm_cnt;
  int ti_mode;
  int ti_region;
  int ti_lscale;
  int mytaskid;
  int si_dihedral_ene;
  int ifmbar_lcl;
  int do_mbar;
  int bar_states;
  int si_dvdl;
  double *ep;//It's a variable.
  dihed_atm_c *atm_in;
  //double *mpe_frc;
  dihed_frc_c *dfrc;

  dfrc_packed_c *atm_force;/////
  int atm_force_len;
  int atm_force_step;
  int *pos_mpe;
  int myrank;

  //double *mpe_epl;
  GBFloat *gmul;
  GBFloat tm24;
  GBFloat tm06;
  GBFloat tenm3;
  GBFloat zero;
  GBFloat one;

}dihed_param_c;

//////////sqrt////////////
#define inv_sqrt(x, r) {                        \
    double y = x;                               \
    double xhalf = 0.5 * y;                     \
    long i = *(long*)(&y);                      \
    i = 0x5fe6ec85e7de30daLL - (i >> 1);        \
    y = *(double*)(&i);                         \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    r = y;                                      \
  }

