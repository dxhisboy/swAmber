#include <simd.h>
typedef struct gbl_accum_vars{
  double eedvir, eed, evdw, ehb;
  double vxx, vxy, vxz, vyy, vyz, vzz;
} gbl_accum_vars_t;

typedef struct gbl_accum_vars_v4{
  doublev4 eedvir, eed, evdw, ehb;
  doublev4 vxx, vxy, vxz, vyy, vyz, vzz;
} gbl_accum_vars_v4_t;

typedef struct ico_tbl {
  double *cn1, *cn2, *asol, *bsol;
  int *ico;
  int ntypes, padding;
} ico_tbl_t;

typedef struct packed_ico_tbl {
  //ic > 0: cn1 = p1, cn2 = p2, icd = 2
  //ic < 0: asol = p1, bsol = p2, icd = -2
  double *p1, *p2, *icd;
  int ntypes, padding;
} packed_ico_tbl_t;

typedef struct a2i_target {
  double bkt_tran[3];
  int img_i, jlo, jhi, padding;
} a2i_target_t;

typedef struct a2b_target {
  double bkt_tran[3];
  double *i_frc;
  double *i_crd;
  double i_qterm;
  int i_iac;
  double (*j_frc)[3];
  double (*j_crd)[3];
  double *j_qterm;
  int *j_iac;
  int j_cnt;  
} a2b_target_t;

typedef struct v2b_target {
  doublev4 i_qterm;
  doublev4 iunmask;
  doublev4 *ij_unmask;
  doublev4 bkt_tran[3];
  int iimg[4];
  doublev4 *i_frc;
  doublev4 *i_crd;
  doublev4 *i_p1, *i_p2, *i_icd;
  double (*j_frc)[3];
  doublev4 (*j_frc_v4)[3];
  double (*j_crd)[3];
  double *j_qterm;
  int *j_iac;
  int j_cnt, jlo;  
} v2b_target_t;


typedef struct share_vars {
  double max_nb_cut2, es_cut2, p12, p6, dxdr;
  double cut, cut6, cut3;
  double cutinv, cut2inv, cut3inv, cut6inv;
  double fswitch2, fswitch3, fswitch6;
  double invfswitch6cut6, invfswitch3cut3;
  int ifunc, padding;
} share_vars_t;

typedef struct share_vars_v4 {
  doublev4 max_nb_cut2, es_cut2, p12, p6, dxdr;
  doublev4 cut, cut6, cut3;
  doublev4 cutinv, cut2inv, cut3inv, cut6inv;
  doublev4 fswitch2, fswitch3, fswitch6;
  doublev4 invfswitch6cut6, invfswitch3cut3;
  int ifunc, padding;
} share_vars_v4_t;

typedef struct img_vars {
  double (*img_frc)[3];
  double (*img_frc_rep)[3];
  double (*img_crd)[3];
  double *img_qterm;
  int *img_iac, *img_atm_map;
  double tranvec[18][3];
  int my_img_lo, my_img_hi, my_bkt_lo, my_bkt_hi, atm_cnt, padding;
} img_vars_t;

typedef struct cit_tbl_rec {
  int img_lo, img_hi;
} cit_tbl_rec_t;
typedef struct neigh_cit_entry{
  int img_lo, img_hi, itran, cit_no, first, irep;
} neigh_cit_entry_t;
typedef struct cit_vars {
  cit_tbl_rec_t *cits;
  neigh_cit_entry_t (*neigh_cits)[64];
  double (*bbox)[8];
  int *bkt_upd_cnt;
  int cit_tbl_x_dim, cit_tbl_y_dim, cit_tbl_z_dim, cit_bkt_delta;
  int dbkt[64][3];
  int ndbkt, padding;
} cit_vars_t;

typedef struct nb_energy_arg {
  img_vars_t img_vars;
  cit_vars_t cit_vars;
  share_vars_t share_vars;
  ico_tbl_t ico_tbl;
  packed_ico_tbl_t packed_ico;
  gbl_accum_vars_t *accum;
  int mytaskid, padding;
} nb_energy_arg_t;


typedef struct pairs_calc_a2i_arg {
  img_vars_t img_vars;
  share_vars_t share_vars;
  a2i_target_t a2i;
  ico_tbl_t ico_tbl;
  packed_ico_tbl_t packed_ico;
  gbl_accum_vars_t *accum;
  int mytaskid, padding;
} pairs_calc_a2i_arg_t;

typedef struct pairs_calc_a2b_arg {
  //img_vars_t img_vars;
  share_vars_t share_vars;
  a2b_target_t a2b;
  ico_tbl_t ico_tbl;
  packed_ico_tbl_t packed_ico;
  gbl_accum_vars_t *accum;
  int mytaskid, padding;
} pairs_calc_a2b_arg_t;

typedef struct pairs_calc_v2b_arg {
  img_vars_t img_vars;
  share_vars_v4_t share_vars;
  v2b_target_t v2b;
  ico_tbl_t ico_tbl;
  gbl_accum_vars_v4_t *accum;
  int mytaskid, padding;
} pairs_calc_v2b_arg_t;
