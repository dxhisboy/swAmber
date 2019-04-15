#include <simd.h>
#include <math.h>
#include <stdio.h>
#define RSQRT_PI 0.56418958354775628
#define MAX_BKT_SIZE 64
#define MAX_TYPE_SQ 1024
#define MAX_TYPE 32
#define dbg -197372

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

typedef struct v2l_target {
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
  int j_nlst, *j_lst;
} v2l_target_t;

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

typedef struct pairs_calc_v2l_arg {
  img_vars_t img_vars;
  share_vars_v4_t share_vars;
  v2l_target_t v2l;
  ico_tbl_t ico_tbl;
  gbl_accum_vars_v4_t *accum;
  int mytaskid, padding;
} pairs_calc_v2l_arg_t;

void cit_wrap_tran(int *cx, int *cy, int *cz, int *tran, int xlim, int ylim, int zlim){
  *tran = 9 + 3 + 1;
  if (*cx < 0){
    *cx += xlim;
    *tran -= 1;
  }
  if (*cx >= xlim){
    *cx -= xlim;
    *tran += 1;
  }
  if (*cy < 0){
    *cy += ylim;
    *tran -= 3;
  }
  if (*cy >= ylim){
    *cy -= ylim;
    *tran += 3;
  }
  if (*cz < 0){
    *cz += zlim;
    *tran -= 9;
  }
  if (*cz >= zlim){
    *cz -= zlim;
    *tran += 9;
  }
}

void pack_ico_tbl(ico_tbl_t *ico, packed_ico_tbl_t *packed_ico){
  int ntypes = ico->ntypes;
  packed_ico->ntypes = ntypes;
  int iac, jac;
  for (iac = 0; iac < ntypes; iac ++){
    for (jac = 0; jac < ntypes; jac ++){
      int ic = ico->ico[iac * ntypes + jac];
      if (ic != 0) {
#ifdef HAS_10_12
        packed_ico->icd[iac * ntypes + jac] = ic > 0 ? 2.0 : -2.0;
        packed_ico->p1[iac * ntypes + jac] = ic > 0 ? ico->cn1[ic - 1] : ico->asol[-ic - 1];
        packed_ico->p2[iac * ntypes + jac] = ic > 0 ? ico->cn2[ic - 1] : ico->bsol[-ic - 1];
#else
        packed_ico->icd[iac * ntypes + jac] = 2.0;
        packed_ico->p1[iac * ntypes + jac] = ico->cn1[ic - 1];
        packed_ico->p2[iac * ntypes + jac] = ico->cn2[ic - 1];
#endif
      } else {
        packed_ico->icd[iac * ntypes + jac] = 0.0;
      }
    }
  }
}

void scan_bbox(int j_cnt, double (*j_crd)[3], double *bbox_c, double *bbox_h){
  //lwpf_start(BBOX);
  int img_j;
  double bbox_lo[3], bbox_hi[3];
  bbox_lo[0] = 1e19;
  bbox_lo[1] = 1e19;
  bbox_lo[2] = 1e19;
  bbox_hi[0] = -1e19;
  bbox_hi[1] = -1e19;
  bbox_hi[2] = -1e19;

  for (img_j = 0; img_j < j_cnt; img_j ++){
    if (j_crd[img_j][0] < bbox_lo[0]) bbox_lo[0] = j_crd[img_j][0];
    if (j_crd[img_j][1] < bbox_lo[1]) bbox_lo[1] = j_crd[img_j][1];
    if (j_crd[img_j][2] < bbox_lo[2]) bbox_lo[2] = j_crd[img_j][2];
    if (j_crd[img_j][0] > bbox_hi[0]) bbox_hi[0] = j_crd[img_j][0];
    if (j_crd[img_j][1] > bbox_hi[1]) bbox_hi[1] = j_crd[img_j][1];
    if (j_crd[img_j][2] > bbox_hi[2]) bbox_hi[2] = j_crd[img_j][2];
  }
  bbox_c[0] = (bbox_lo[0] + bbox_hi[0]) * 0.5;
  bbox_c[1] = (bbox_lo[1] + bbox_hi[1]) * 0.5;
  bbox_c[2] = (bbox_lo[2] + bbox_hi[2]) * 0.5;

  bbox_h[0] = bbox_hi[0] - bbox_c[0];
  bbox_h[1] = bbox_hi[1] - bbox_c[1];
  bbox_h[2] = bbox_hi[2] - bbox_c[2];
  //lwpf_stop(BBOX);
}

#ifdef MPE
#include <athread.h>
extern void slave_get_nb_energy_enterance(nb_energy_arg_t *);
extern void slave_get_nb_energy_enterance2(nb_energy_arg_t *);
extern void slave_get_nb_energy_reduction_excl(nb_energy_arg_t *);
static int r = 0;
#define LWPF_UNITS U(NB_ENERGY)
#include <lwpf2/lwpf2.h>
#include <gptl.h>
perf_config_t conf;

double icd[MAX_TYPE_SQ], p1[MAX_TYPE_SQ], p2[MAX_TYPE_SQ];
gbl_accum_vars_t accum[64], *accum_in_save;

void get_nb_energy_main_loop_spawn_(nb_energy_arg_t *arg){
  if (!athread_idle()) athread_init();
  conf.pcrc = PCRC_ALL;
  conf.pcr0 = PC0_N_INST;
  conf.pcr1 = PC1_T_DATA_REL_BLOCK;
  conf.pcr2 = PC2_T_PIPE_MISS;
  lwpf_init(&conf);
  arg->packed_ico.icd = icd;
  arg->packed_ico.p1 = p1;
  arg->packed_ico.p2 = p2;
  pack_ico_tbl(&(arg->ico_tbl), &(arg->packed_ico));

  accum_in_save = arg->accum;
  memset(accum, 0, sizeof(gbl_accum_vars_t) * 64);
  arg->accum = accum;
  int icit;
  int my_bkt_lo = arg->img_vars.my_bkt_lo;
  int my_bkt_hi = arg->img_vars.my_bkt_hi;
  double (*bbox)[8] = arg->cit_vars.bbox;
  double (*img_crd)[3] = arg->img_vars.img_crd;
  int ncit = arg->cit_vars.cit_tbl_x_dim * arg->cit_vars.cit_tbl_y_dim * arg->cit_vars.cit_tbl_z_dim;
  //printf("%p\n", bbox);
  /* GPTLstart("bbox scan"); */
  /* for (icit = 0; icit < ncit; icit ++){ */
  /*   if (arg->cit_vars.bkt_upd_cnt[icit] > 0) { */
  /*     int lo = arg->cit_vars.cits[icit].img_lo - 1; */
  /*     int hi = arg->cit_vars.cits[icit].img_hi; */
  /*     scan_bbox(hi - lo, img_crd + lo, bbox[icit], bbox[icit] + 4); */
  /*   } */
  /* } */
  /* GPTLstop("bbox scan"); */
  GPTLstart("main loop cpe");
  athread_spawn(get_nb_energy_enterance, arg);
}

void get_nb_energy_main_loop_join_(nb_energy_arg_t *arg){
  athread_join();
  int    atm_cnt       = arg->img_vars.atm_cnt;
  int watch = atm_cnt + 197372;
  double (*img_frc_rep)[3] = arg->img_vars.img_frc_rep;
  //if (arg->mytaskid == 1) printf("after main: %f %f %f\n", img_frc_rep[watch][0], img_frc_rep[watch][1], img_frc_rep[watch][2]);

  GPTLstop("main loop cpe");
  /* athread_spawn(get_nb_energy_enterance2, arg); */
  /* athread_join(); */
  //double (*img_frc_rep)[3] = arg->img_vars.img_frc_rep;
  //if (arg->mytaskid == 1) printf("after reduction: %f %f %f\n", img_frc_rep[watch][0], img_frc_rep[watch][1], img_frc_rep[watch][2]);

  GPTLstart("get_nb_energy_reduction");
  athread_spawn(get_nb_energy_reduction_excl, arg);
  if (arg->mytaskid == 0) {
    lwpf_report_summary(stdout, &conf);
    //lwpf_report_detail(stdout, &conf);
  }
  r++;
  int i;
  for (i = 0; i < 64; i ++){
    accum_in_save->vxx += accum[i].vxx;
    accum_in_save->vxy += accum[i].vxy;
    accum_in_save->vxz += accum[i].vxz;
    accum_in_save->vyy += accum[i].vyy;
    accum_in_save->vyz += accum[i].vyz;
    accum_in_save->vzz += accum[i].vzz;
    accum_in_save->eedvir += accum[i].eedvir;
    accum_in_save->eed += accum[i].eed;
    accum_in_save->evdw += accum[i].evdw;
    accum_in_save->ehb += accum[i].ehb;
  }
  athread_join();
  GPTLstop("get_nb_energy_reduction");
}

#endif

#ifdef CPE
#include <slave.h>
#include <dma.h>
#include <assert.h>
//#define J_UPD_SCALAR
#define vmatch(r, a, b) \
  asm("vmatch %1, %2, %0" : "=&r"(r) : "r"(a), "r"(b))
#define simd_vselgt(a, b, c) simd_vselle(a, c, b)
#define simd_vselge(a, b, c) simd_vsellt(a, c, b)
#define simd_vselne(a, b, c) simd_vseleq(a, c, b)

#define simd_bcastf(x) simd_vshff((doublev4)(x), (doublev4)(x), 0)
//#define simd_vlde(v4, fptr) asm volatile("vlde %0, 0(%s)" : "=r"());
#define align_ptr(x) (void*)(((long)(x)) + 31 & ~31)
#define transpose4x4_2x4(in0, in1, in2, in3, ot0, ot1) {                \
  doublev4 o0 = simd_vshff(in1,in0,0x44);                             \
  doublev4 o2 = simd_vshff(in3,in2,0x44);                             \
  ot0 = simd_vshff(o2,o0,0x88);                                       \
  ot1 = simd_vshff(o2,o0,0xDD);                                       \
  }

#define transpose4x4(in0, in1, in2, in3, ot0, ot1, ot2, ot3) {  \
  doublev4 o0 = simd_vshff(in1,in0,0x44);                     \
  doublev4 o1 = simd_vshff(in1,in0,0xEE);                     \
  doublev4 o2 = simd_vshff(in3,in2,0x44);                     \
  doublev4 o3 = simd_vshff(in3,in2,0xEE);                     \
  ot0 = simd_vshff(o2,o0,0x88);                               \
  ot1 = simd_vshff(o2,o0,0xDD);                               \
  ot2 = simd_vshff(o3,o1,0x88);                               \
  ot3 = simd_vshff(o3,o1,0xDD);                               \
  }
#define vshuffd_rc(a, b, c, d) (d | (c << 2) | (b << 4) | (a << 6))
#define simd_vsumd(x) {                                 \
  x += simd_vshff(x, x, vshuffd_rc(2, 3, 0, 1));      \
  x += simd_vshff(x, x, vshuffd_rc(1, 0, 3, 2));      \
  }

#define LWPF_UNIT U(NB_ENERGY)
//#define LWPF_KERNELS K(ALL) K(MAIN) K(GET_J) K(GET_I) K(UPD_FRC) K(COMP) K(V2B) K(EE) K(VDW) K(SWITCH) K(CALC) K(REDUCE) K(INVSQRT) K(DELR) K(LOOP_J) K(EXT) K(SCAN) K(BBOX)
//#define LWPF_KERNELS K(ALL) K(MAIN) K(GET_J) K(GET_I) K(UPD_FRC) K(COMP) K(V2B) K(EE) K(VDW) K(SWITCH) K(CALC) K(REDUCE) K(INVSQRT) K(DELR) K(LOOP_J) K(EXT) K(SCAN) K(BBOX)
//#define LWPF_KERNELS K(ALL) K(BBOX) K(PREF) K(DIST2) K(UPD_FRC) K(LDD_I) K(V2B) K(VDW) K(EE) K(REDUCE) K(LOOP_J)//K(MAIN) K(GET_J) K(GET_I)  K(COMP)  K(EE)  K(SWITCH) K(CALC)  K(INVSQRT) K(DELR) K(LOOP_J) K(EXT) K(SCAN) K(BBOX)
#define LWPF_KERNELS K(ALL) K(BBOX) K(MAIN) K(REDUCE_CPE) K(REDUCE_EXCL) K(J_CHK) K(J_LOOP)
#include <lwpf2/lwpf2.h>

#include <dma_macros.h>
#include <pairs_calc_wrapper.h>

void get_task_range(int my_bkt_hi, int my_bkt_lo, int *pe_bkt_hi, int *pe_bkt_lo){
  int ncit = my_bkt_hi - my_bkt_lo + 1;
  int ncit_rem = ncit & 63;
  int ncit_pe = (ncit + 63) >> 6;
  int s1 = ncit_rem * ncit_pe;
  int lo, hi;
  if (ncit_rem == 0){
    lo = _MYID * ncit_pe + my_bkt_lo;
    hi = lo + ncit_pe;
  } else {
    if (_MYID < ncit_rem){
      lo = _MYID * ncit_pe + my_bkt_lo;
      hi = lo + ncit_pe;
    } else {
      lo = (_MYID - ncit_rem) * (ncit_pe - 1) + s1 + my_bkt_lo;
      hi = lo + ncit_pe - 1;
    }
  }
  *pe_bkt_hi = hi;
  *pe_bkt_lo = lo;
}
void parameters_prefetch(int ntypes, double *p1, double *p2, double *icd, doublev4 *p1_v4, doublev4 *p2_v4, doublev4 *icd_v4, int iac0, int iac1, int iac2, int iac3) {
  //lwpf_start(PREF);
  int jac;
  for (jac = 0; jac < ntypes; jac ++){
    icd_v4[jac] = simd_vinsf0((doublev4)icd[(iac0 - 1) * MAX_TYPE + jac], icd_v4[jac]);
    icd_v4[jac] = simd_vinsf1((doublev4)icd[(iac1 - 1) * MAX_TYPE + jac], icd_v4[jac]);
    icd_v4[jac] = simd_vinsf2((doublev4)icd[(iac2 - 1) * MAX_TYPE + jac], icd_v4[jac]);
    icd_v4[jac] = simd_vinsf3((doublev4)icd[(iac3 - 1) * MAX_TYPE + jac], icd_v4[jac]);

    p1_v4[jac] = simd_vinsf0((doublev4)p1[(iac0 - 1) * MAX_TYPE + jac], p1_v4[jac]);
    p1_v4[jac] = simd_vinsf1((doublev4)p1[(iac1 - 1) * MAX_TYPE + jac], p1_v4[jac]);
    p1_v4[jac] = simd_vinsf2((doublev4)p1[(iac2 - 1) * MAX_TYPE + jac], p1_v4[jac]);
    p1_v4[jac] = simd_vinsf3((doublev4)p1[(iac3 - 1) * MAX_TYPE + jac], p1_v4[jac]);

    p2_v4[jac] = simd_vinsf0((doublev4)p2[(iac0 - 1) * MAX_TYPE + jac], p2_v4[jac]);
    p2_v4[jac] = simd_vinsf1((doublev4)p2[(iac1 - 1) * MAX_TYPE + jac], p2_v4[jac]);
    p2_v4[jac] = simd_vinsf2((doublev4)p2[(iac2 - 1) * MAX_TYPE + jac], p2_v4[jac]);
    p2_v4[jac] = simd_vinsf3((doublev4)p2[(iac3 - 1) * MAX_TYPE + jac], p2_v4[jac]);
  }
  //lwpf_stop(PREF);
}

void parameters_prefetch_v4(int ntypes, double *p1, double *p2, double *icd, doublev4 *p1_v4, doublev4 *p2_v4, doublev4 *icd_v4, int iac0, int iac1, int iac2, int iac3) {
  int jac;
  //lwpf_start(PREF);
  for (jac = 0; jac < ntypes; jac += 4){
    doublev4 i0_icd_v4, i1_icd_v4, i2_icd_v4, i3_icd_v4;
    simd_load(i0_icd_v4, icd + (iac0 - 1) * MAX_TYPE + jac);
    simd_load(i1_icd_v4, icd + (iac1 - 1) * MAX_TYPE + jac);
    simd_load(i2_icd_v4, icd + (iac2 - 1) * MAX_TYPE + jac);
    simd_load(i3_icd_v4, icd + (iac3 - 1) * MAX_TYPE + jac);
    transpose4x4(i0_icd_v4, i1_icd_v4, i2_icd_v4, i3_icd_v4, icd_v4[jac], icd_v4[jac + 1], icd_v4[jac + 2], icd_v4[jac + 3]);
  }
  for (jac = 0; jac < ntypes; jac += 4){
    doublev4 i0_p1_v4, i1_p1_v4, i2_p1_v4, i3_p1_v4;
    simd_load(i0_p1_v4, p1 + (iac0 - 1) * MAX_TYPE + jac);
    simd_load(i1_p1_v4, p1 + (iac1 - 1) * MAX_TYPE + jac);
    simd_load(i2_p1_v4, p1 + (iac2 - 1) * MAX_TYPE + jac);
    simd_load(i3_p1_v4, p1 + (iac3 - 1) * MAX_TYPE + jac);
    transpose4x4(i0_p1_v4, i1_p1_v4, i2_p1_v4, i3_p1_v4, p1_v4[jac], p1_v4[jac + 1], p1_v4[jac + 2], p1_v4[jac + 3]);
  }
  for (jac = 0; jac < ntypes; jac += 4){
    doublev4 i0_p2_v4, i1_p2_v4, i2_p2_v4, i3_p2_v4;
    simd_load(i0_p2_v4, p2 + (iac0 - 1) * MAX_TYPE + jac);
    simd_load(i1_p2_v4, p2 + (iac1 - 1) * MAX_TYPE + jac);
    simd_load(i2_p2_v4, p2 + (iac2 - 1) * MAX_TYPE + jac);
    simd_load(i3_p2_v4, p2 + (iac3 - 1) * MAX_TYPE + jac);
    transpose4x4(i0_p2_v4, i1_p2_v4, i2_p2_v4, i3_p2_v4, p2_v4[jac], p2_v4[jac + 1], p2_v4[jac + 2], p2_v4[jac + 3]);
  }
  //lwpf_stop(PREF);
}

void parameters_prefetch32(int ntypes, double *p, doublev4 *p_v4, int iac0, int iac1, int iac2, int iac3) {
  int jac;
  //lwpf_start(PREF);
  for (jac = 0; jac < ntypes; jac += 16){
    doublev4 i0_p_v4_1, i1_p_v4_1, i2_p_v4_1, i3_p_v4_1;
    doublev4 i0_p_v4_2, i1_p_v4_2, i2_p_v4_2, i3_p_v4_2;
    simd_load(i0_p_v4_1, p + (iac0 - 1) * MAX_TYPE + jac);
    simd_load(i1_p_v4_1, p + (iac1 - 1) * MAX_TYPE + jac);
    simd_load(i2_p_v4_1, p + (iac2 - 1) * MAX_TYPE + jac);
    simd_load(i3_p_v4_1, p + (iac3 - 1) * MAX_TYPE + jac);
    simd_load(i0_p_v4_2, p + (iac0 - 1) * MAX_TYPE + jac + 4);
    simd_load(i1_p_v4_2, p + (iac1 - 1) * MAX_TYPE + jac + 4);
    simd_load(i2_p_v4_2, p + (iac2 - 1) * MAX_TYPE + jac + 4);
    simd_load(i3_p_v4_2, p + (iac3 - 1) * MAX_TYPE + jac + 4);
    transpose4x4(i0_p_v4_1, i1_p_v4_1, i2_p_v4_1, i3_p_v4_1, p_v4[jac], p_v4[jac + 1], p_v4[jac + 2], p_v4[jac + 3]);
    transpose4x4(i0_p_v4_2, i1_p_v4_2, i2_p_v4_2, i3_p_v4_2, p_v4[jac + 4], p_v4[jac + 5], p_v4[jac + 6], p_v4[jac + 7]);
    simd_load(i0_p_v4_1, p + (iac0 - 1) * MAX_TYPE + jac + 8);
    simd_load(i1_p_v4_1, p + (iac1 - 1) * MAX_TYPE + jac + 8);
    simd_load(i2_p_v4_1, p + (iac2 - 1) * MAX_TYPE + jac + 8);
    simd_load(i3_p_v4_1, p + (iac3 - 1) * MAX_TYPE + jac + 8);
    simd_load(i0_p_v4_2, p + (iac0 - 1) * MAX_TYPE + jac + 12);
    simd_load(i1_p_v4_2, p + (iac1 - 1) * MAX_TYPE + jac + 12);
    simd_load(i2_p_v4_2, p + (iac2 - 1) * MAX_TYPE + jac + 12);
    simd_load(i3_p_v4_2, p + (iac3 - 1) * MAX_TYPE + jac + 12);
    transpose4x4(i0_p_v4_1, i1_p_v4_1, i2_p_v4_1, i3_p_v4_1, p_v4[jac + 8], p_v4[jac + 9], p_v4[jac + 10], p_v4[jac + 11]);
    transpose4x4(i0_p_v4_2, i1_p_v4_2, i2_p_v4_2, i3_p_v4_2, p_v4[jac + 12], p_v4[jac + 13], p_v4[jac + 14], p_v4[jac + 15]);
  }
  //lwpf_stop(PREF);
}
extern void parameters_prefetch_32(double *p, doublev4 *p_v4, int iac0, int iac1, int iac2, int iac3);
//                                        a0              a1       a2        a3        a4         a5
inline double dist2_crd_box(double *crd, double *box_lo, double *box_hi, double *tranvec){
  double box_center[3];
  box_center[0] = (box_lo[0] + box_hi[0]) * 0.5;
  box_center[1] = (box_lo[1] + box_hi[1]) * 0.5;
  box_center[2] = (box_lo[2] + box_hi[2]) * 0.5;
  double v[3];
  v[0] = box_center[0] + tranvec[0] - crd[0];
  v[1] = box_center[1] + tranvec[1] - crd[1];
  v[2] = box_center[2] + tranvec[2] - crd[2];
  asm("fcpys $31, %3, %0\n\t"
      "fcpys $31, %4, %1\n\t"
      "fcpys $31, %5, %2\n\t"
      : "=r"(v[0]), "=r"(v[1]), "=r"(v[2])
      : "r"(v[0]), "r"(v[1]), "r"(v[2]));
  double h[3];
  h[0] = box_hi[0] - box_center[0];
  h[1] = box_hi[1] - box_center[1];
  h[2] = box_hi[2] - box_center[2];
  double u[3];
  u[0] = v[0] - h[0];
  u[1] = v[1] - h[1];
  u[2] = v[2] - h[2];
  if (u[0] < 0) u[0] = 0;
  if (u[1] < 0) u[1] = 0;
  if (u[2] < 0) u[2] = 0;
  return u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
}

double dist2_crd_box_new(double *crd, double *box_c, double *box_h){
  //lwpf_start(DIST2);
  double v[3];
  v[0] = box_c[0] - crd[0];
  v[1] = box_c[1] - crd[1];
  v[2] = box_c[2] - crd[2];
  asm("fcpys $31, %3, %0\n\t"
      "fcpys $31, %4, %1\n\t"
      "fcpys $31, %5, %2\n\t"
      : "=r"(v[0]), "=r"(v[1]), "=r"(v[2])
      : "r"(v[0]), "r"(v[1]), "r"(v[2]));
  double u[3];
  u[0] = v[0] - box_h[0];
  u[1] = v[1] - box_h[1];
  u[2] = v[2] - box_h[2];
  asm ("fsellt %3, $31, %3, %0\n\t"
       "fsellt %4, $31, %4, %1\n\t"
       "fsellt %5, $31, %5, %2\n\t"
       : "=r"(u[0]), "=r"(u[1]), "=r"(u[2])
       : "r"(u[0]), "r"(u[1]), "r"(u[2]));
  /* if (u[0] < 0) u[0] = 0; */
  /* if (u[1] < 0) u[1] = 0; */
  /* if (u[2] < 0) u[2] = 0; */
  double ret = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  //lwpf_stop(DIST2);

  return ret;
}

void copy_share_vars_to_vec(share_vars_t *in, share_vars_v4_t *out){
  out->max_nb_cut2     = in->max_nb_cut2    ;
  out->es_cut2         = in->es_cut2        ;
  out->p12             = in->p12            ;
  out->p6              = in->p6             ;
  out->dxdr            = in->dxdr           ;
  out->cut             = in->cut            ;
  out->cut6            = in->cut6           ;
  out->cut3            = in->cut3           ;
  out->cutinv          = in->cutinv         ;
  out->cut2inv         = in->cut2inv        ;
  out->cut3inv         = in->cut3inv        ;
  out->cut6inv         = in->cut6inv        ;
  out->fswitch2        = in->fswitch2       ;
  out->fswitch3        = in->fswitch3       ;
  out->fswitch6        = in->fswitch6       ;
  out->invfswitch6cut6 = in->invfswitch6cut6;
  out->invfswitch3cut3 = in->invfswitch3cut3;
  out->ifunc           = in->ifunc;
}
#define ICIT_STEP 16
void get_nb_energy_bbox_scan(nb_energy_arg_t *arg){
  dma_init();
  double (*gbl_bbox)[8] = arg->cit_vars.bbox;
  double (*img_crd)[3] = arg->img_vars.img_crd;
  int ncit = arg->cit_vars.cit_tbl_x_dim * arg->cit_vars.cit_tbl_y_dim * arg->cit_vars.cit_tbl_z_dim;
  int bkt_upd_cnt[ICIT_STEP];
  cit_tbl_rec_t cits[ICIT_STEP];
  double bbox[ICIT_STEP][8];
  double j_crd[MAX_BKT_SIZE][3];

  int ist, ied, icit;
  for (ist = _MYID * ICIT_STEP; ist < ncit; ist += ICIT_STEP * 64){
    ied = ist + ICIT_STEP;
    if (ied > ncit) ied = ncit;
    pe_get(arg->cit_vars.bkt_upd_cnt + ist, bkt_upd_cnt, sizeof(int) * (ied - ist));
    pe_get(arg->cit_vars.cits        + ist, cits       , sizeof(cit_tbl_rec_t) * (ied - ist));
    dma_syn();
    for (icit = ist; icit < ied; icit ++){
      if (bkt_upd_cnt[icit - ist] > 0){
        int lo = cits[icit - ist].img_lo - 1;
        int hi = cits[icit - ist].img_hi;
        if (hi - lo > MAX_BKT_SIZE) printf("icit=%d hi-lo=%d excceed %d", icit, hi - lo, MAX_BKT_SIZE);
        if (lo < hi){
          pe_get(img_crd + lo, j_crd, sizeof(double) * 3 * (hi - lo));
          dma_syn();
          scan_bbox(hi - lo, j_crd, bbox[icit - ist], bbox[icit - ist] + 4);
        }
      }
    }
    pe_put(gbl_bbox + ist, bbox, sizeof(double) * 8 * (ied - ist));
    dma_syn();
  }
}

void get_nb_energy_main_loop(nb_energy_arg_t *arg){
  /* lwpf_enter(NB_ENERGY); */
  /* lwpf_start(ALL); */
  dma_init();

  double (*img_frc)[3] = arg->img_vars.img_frc  ;
  double (*img_crd)[3] = arg->img_vars.img_crd  ;
  double *img_qterm    = arg->img_vars.img_qterm;
  int    *img_iac      = arg->img_vars.img_iac  ;
  int    my_bkt_lo     = arg->img_vars.my_bkt_lo;
  int    my_bkt_hi     = arg->img_vars.my_bkt_hi;
  int    my_img_lo     = arg->img_vars.my_img_lo - 1;
  int    my_img_hi     = arg->img_vars.my_img_hi;
  int    atm_cnt       = arg->img_vars.atm_cnt;
  double (*img_frc_rep)[3] = arg->img_vars.img_frc_rep;
  cit_tbl_rec_t *cits = arg->cit_vars.cits;
  double (*bbox)[8] = arg->cit_vars.bbox;
  neigh_cit_entry_t (*gbl_neigh_cits)[64] = arg->cit_vars.neigh_cits;
  int cit_tbl_x_dim   = arg->cit_vars.cit_tbl_x_dim;
  int cit_tbl_y_dim   = arg->cit_vars.cit_tbl_y_dim;
  int cit_tbl_z_dim   = arg->cit_vars.cit_tbl_z_dim;
  int ndbkt           = arg->cit_vars.ndbkt;
  int idbkt, idim;

  int ntypes   = arg->packed_ico.ntypes;
  int ico[MAX_TYPE_SQ];
  double p1[MAX_TYPE_SQ], p2[MAX_TYPE_SQ], icd[MAX_TYPE_SQ];
  doublev4 i_p1[MAX_TYPE], i_p2[MAX_TYPE], i_icd[MAX_TYPE];
  pe_get(arg->packed_ico.p1 , p1 , ntypes * ntypes * sizeof(double));
  pe_get(arg->packed_ico.p2 , p2 , ntypes * ntypes * sizeof(double));
  pe_get(arg->packed_ico.icd, icd, ntypes * ntypes * sizeof(double));
  dma_syn();
  int iac, jac;
  for (iac = ntypes - 1; iac >= 0; iac --){
    for (jac = ntypes - 1; jac >= 0; jac --){
      p1[iac * MAX_TYPE + jac] = p1[iac * ntypes + jac];
      p2[iac * MAX_TYPE + jac] = p2[iac * ntypes + jac];
      icd[iac * MAX_TYPE + jac] = icd[iac * ntypes + jac];
    }
  }
  pairs_calc_v2l_arg_t pairs_v2l_arg;
  gbl_accum_vars_v4_t accum;
  memset(&accum, 0, sizeof(gbl_accum_vars_v4_t));

  copy_share_vars_to_vec(&(arg->share_vars), &(pairs_v2l_arg.share_vars));
  pairs_v2l_arg.img_vars = arg->img_vars;
  pairs_v2l_arg.accum    = &accum;
  pairs_v2l_arg.mytaskid = arg->mytaskid;
  pairs_v2l_arg.v2l.i_p1  = i_p1;
  pairs_v2l_arg.v2l.i_p2  = i_p2;
  pairs_v2l_arg.v2l.i_icd = i_icd;

  double tranvec[18][3];
  int itran;
  for (itran = 0; itran < 18; itran ++){
    tranvec[itran][0] = arg->img_vars.tranvec[itran][0];
    tranvec[itran][1] = arg->img_vars.tranvec[itran][1];
    tranvec[itran][2] = arg->img_vars.tranvec[itran][2];
  }

  void (*pairs_v2l_calc_func)(pairs_calc_v2l_arg_t *) = pairs_calc_v2l_funcs[arg->share_vars.ifunc];
  int icit;
  double j_crd[MAX_BKT_SIZE][3], j_qterm[MAX_BKT_SIZE], j_frc[MAX_BKT_SIZE][3];
  doublev4 j_frc_v4[MAX_BKT_SIZE][3];
  doublev4 i_vec[MAX_BKT_SIZE];
  doublev4 ij_unmask[MAX_BKT_SIZE];
  int j_lst[MAX_BKT_SIZE], j_nlst;
  int img_j;
  for (img_j = 0; img_j < MAX_BKT_SIZE; img_j ++) {
    ij_unmask[img_j] = 2.0;
  }
  int j_iac[MAX_BKT_SIZE];
  double i_frc[MAX_BKT_SIZE][3], i_crd[MAX_BKT_SIZE][3], i_qterm[MAX_BKT_SIZE];
  long i_inrange[MAX_BKT_SIZE];
  int i_iac[MAX_BKT_SIZE];
  //lwpf_start(MAIN);
  int pe_bkt_lo, pe_bkt_hi;
  /* if (_MYID) return; */
  /* for (_MYID = 0; _MYID < 64; _MYID ++){ */
  get_task_range(my_bkt_hi, my_bkt_lo, &pe_bkt_hi, &pe_bkt_lo);
  doublev4 i_unmask[5];
  i_unmask[0] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);
  i_unmask[1] = simd_set_doublev4(2.0, 0.0, 0.0, 0.0);
  i_unmask[2] = simd_set_doublev4(2.0, 2.0, 0.0, 0.0);
  i_unmask[3] = simd_set_doublev4(2.0, 2.0, 2.0, 0.0);
  i_unmask[4] = simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
  double max_nb_cut2 = arg->share_vars.max_nb_cut2;
  //int cnt_calc = 0, cnt_calc_group = 0;
  for (icit = pe_bkt_lo; icit < pe_bkt_hi; icit ++){
    //for (icit = my_bkt_lo + _MYID; icit <= my_bkt_hi; icit += 64){
    neigh_cit_entry_t neigh_cits[64];
    pe_get(gbl_neigh_cits[icit - my_bkt_lo], neigh_cits, sizeof(neigh_cit_entry_t) * ndbkt);
    dma_syn();
    int cit_lo = neigh_cits[ndbkt - 1].img_lo - 1;
    int cit_hi = neigh_cits[ndbkt - 1].img_hi;

    //assert(i_cnt <= MAX_BKT_SIZE);
    if (my_img_lo > cit_lo) cit_lo = my_img_lo;
    if (my_img_hi < cit_hi) cit_hi = my_img_hi;
    int i_cnt = cit_hi - cit_lo;
    //lwpf_start(GET_I);
    if (i_cnt <= 0) continue;
    doublev4 ibox_c_v4, ibox_h_v4;
    double ibox_c[4], ibox_h[4];
    pe_get(img_crd   + cit_lo, i_crd  , i_cnt * 3 * sizeof(double));
    pe_get(img_qterm + cit_lo, i_qterm, i_cnt * 1 * sizeof(double));
    pe_get(img_iac   + cit_lo, i_iac  , i_cnt * 1 * sizeof(int)   );
    simd_load(ibox_c_v4, bbox[icit]);
    simd_load(ibox_h_v4, bbox[icit] + 4);

    ibox_c[0] = simd_vextf0(ibox_c_v4);
    ibox_c[1] = simd_vextf1(ibox_c_v4);
    ibox_c[2] = simd_vextf2(ibox_c_v4);
    
    ibox_h[0] = simd_vextf0(ibox_h_v4);
    ibox_h[1] = simd_vextf1(ibox_h_v4);
    ibox_h[2] = simd_vextf2(ibox_h_v4);

    memset(i_frc, 0, i_cnt * 3 * sizeof(double));
    dma_syn();

    int img_i;
    for (img_i = 0; img_i < i_cnt; img_i ++){
      i_vec[img_i] = simd_set_doublev4(i_crd[img_i][0], i_crd[img_i][1], i_crd[img_i][2], i_qterm[img_i]);
    }
    //lwpf_stop(GET_I);

    //int idbkt;
    for (idbkt = 0; idbkt < ndbkt; idbkt ++){
      int jhi = neigh_cits[idbkt].img_hi;
      int jlo = neigh_cits[idbkt].img_lo - 1;
      //if (cit_lo == 0 && jlo <= 400566 && jhi > 400566) printf("%f %f %f\n", img_crd[400565][0], img_crd[400565][1], img_crd[400565][2]);
      int irep = neigh_cits[idbkt].irep - 1;
      int first = neigh_cits[idbkt].first;
      int itran = neigh_cits[idbkt].itran;
      int jcit = neigh_cits[idbkt].cit_no;
      int j_cnt = jhi - jlo;
      doublev4 jbox_c_v4, jbox_h_v4;
      double jbox_c[4], jbox_h[4];
      if (j_cnt <= 0) continue;

      //lwpf_start(GET_J);
      pe_get(img_crd   + jlo, j_crd  , j_cnt * 3 * sizeof(double));
      pe_get(img_qterm + jlo, j_qterm, j_cnt * 1 * sizeof(double));
      pe_get(img_iac   + jlo, j_iac  , j_cnt * 1 * sizeof(int)   );
      if (first) {
        memset(j_frc, 0, j_cnt * 3 * sizeof(double));
      } else {
        pe_get(img_frc_rep + jlo + irep * atm_cnt, j_frc  , j_cnt * 3 * sizeof(double));
      }
      simd_load(jbox_c_v4, bbox[jcit]);
      simd_load(jbox_h_v4, bbox[jcit] + 4);
      /* simd_store(jbox_c_v4, jbox_c); */
      /* simd_store(jbox_h_v4, jbox_h); */
      jbox_c[0] = simd_vextf0(jbox_c_v4);
      jbox_c[1] = simd_vextf1(jbox_c_v4);
      jbox_c[2] = simd_vextf2(jbox_c_v4);

      jbox_h[0] = simd_vextf0(jbox_h_v4);
      jbox_h[1] = simd_vextf1(jbox_h_v4);
      jbox_h[2] = simd_vextf2(jbox_h_v4);

      //add a fake node after all atoms
      i_iac[cit_hi - cit_lo] = 1;
      //if (itran < 0 || itran > 20) printf("%d %d %d %d %d\n", i_cnt, _MYID, arg->mytaskid, itran, jcit);
      pairs_v2l_arg.v2l.bkt_tran[0] = tranvec[itran][0];
      pairs_v2l_arg.v2l.bkt_tran[1] = tranvec[itran][1];
      pairs_v2l_arg.v2l.bkt_tran[2] = tranvec[itran][2];
      dma_syn();
      j_nlst = 0;

      /* lwpf_start(J_CHK); */
      for (img_j = 0; img_j < j_cnt; img_j ++){
        double j_crd_tran[3];
        j_crd_tran[0] = j_crd[img_j][0] + tranvec[itran][0];
        j_crd_tran[1] = j_crd[img_j][1] + tranvec[itran][1];
        j_crd_tran[2] = j_crd[img_j][2] + tranvec[itran][2];
        if (dist2_crd_box_new(j_crd_tran, ibox_c, ibox_h) < max_nb_cut2){
          j_lst[j_nlst ++] = img_j;
        }
      }
      lwpf_stop(J_CHK);
      pairs_v2l_arg.v2l.j_lst = j_lst;
      pairs_v2l_arg.v2l.j_nlst = j_nlst;
      /* if (arg->mytaskid == 0 && _MYID == 0){ */
      /*   printf("%d %d\n", j_cnt, j_nlst); */
      /* } */
      if (j_nlst > 0){
        /* int o = 197372 - jlo; */
        /* if (arg->mytaskid == 1 && jlo <= 197372 && irep == 1 && jhi > 197372) printf("read: %d %d %d %f %f %f\n", _MYID, idbkt, first, j_frc[o][0], j_frc[o][1], j_frc[o][2]); */

        //scan_bbox(j_cnt, j_crd, bbox_c, bbox_h);
        jbox_c[0] = jbox_c[0] + tranvec[itran][0];
        jbox_c[1] = jbox_c[1] + tranvec[itran][1];
        jbox_c[2] = jbox_c[2] + tranvec[itran][2];
        //lwpf_stop(SCAN);
        //lwpf_stop(GET_J);
        //Last bucket, update myself.

#ifndef J_UPD_SCALAR      
        for (img_j = 0; img_j < j_cnt; img_j ++) {
          j_frc_v4[img_j][0] = 0.0;
          j_frc_v4[img_j][1] = 0.0;
          j_frc_v4[img_j][2] = 0.0;
        }
#endif

        if (idbkt == ndbkt - 1){
          for (img_i = 0; img_i < i_cnt; img_i ++){
            j_frc[img_i + cit_lo - jlo][0] += i_frc[img_i][0];
            j_frc[img_i + cit_lo - jlo][1] += i_frc[img_i][1];
            j_frc[img_i + cit_lo - jlo][2] += i_frc[img_i][2];
          }
        }
        img_i = cit_lo;
        while (img_i < cit_hi) {        
          doublev4 i_frc_v4[3];
          doublev4 i_crd_v4[3];
          doublev4 i_qterm_v4;
          int ni = 0;
          int iarr[4];
          //lwpf_start(V2B);
          while (ni < 4 && img_i < cit_hi){
            double dist2 = dist2_crd_box_new(i_crd[img_i - cit_lo], jbox_c, jbox_h);
            if (dist2 < max_nb_cut2) {
              iarr[ni ++] = img_i - cit_lo;
            }
            img_i ++;
          }
          if (ni == 0) {
            //lwpf_stop(V2B);
            continue;
          }
          doublev4 iunmask = i_unmask[ni];
          //int nibak = ni;
          while (ni < 4) {
            iarr[ni] = cit_hi - cit_lo; //iarr[ni - 1];
            ni ++;
          }
          //lwpf_start(PREF);
          parameters_prefetch_32(p1, i_p1, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          parameters_prefetch_32(p2, i_p2, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          parameters_prefetch_32(icd, i_icd, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          //lwpf_stop(PREF);

          //lwpf_start(LDD_I);
          i_frc_v4[0] = 0.0;
          i_frc_v4[1] = 0.0;
          i_frc_v4[2] = 0.0;
          transpose4x4(i_vec[iarr[0]], i_vec[iarr[1]], i_vec[iarr[2]], i_vec[iarr[3]], i_crd_v4[0], i_crd_v4[1], i_crd_v4[2], i_qterm_v4);
          /* pairs_v2b_arg.v2b.iimg[0] = 1 << 27; */
          /* pairs_v2b_arg.v2b.iimg[1] = 1 << 27; */
          /* pairs_v2b_arg.v2b.iimg[2] = 1 << 27; */
          /* pairs_v2b_arg.v2b.iimg[3] = 1 << 27; */
          pairs_v2l_arg.v2l.iimg[0] = iarr[0] + cit_lo;
          pairs_v2l_arg.v2l.iimg[1] = iarr[1] + cit_lo;
          pairs_v2l_arg.v2l.iimg[2] = iarr[2] + cit_lo;
          pairs_v2l_arg.v2l.iimg[3] = iarr[3] + cit_lo;

          pairs_v2l_arg.v2l.i_frc   = i_frc_v4; //i_frc[img_i - cit_lo];
          pairs_v2l_arg.v2l.i_crd   = i_crd_v4; //i_crd[img_i - cit_lo];
          pairs_v2l_arg.v2l.i_qterm = i_qterm_v4; //i_qterm[img_i - cit_lo];
          pairs_v2l_arg.v2l.iunmask = iunmask; //i_qterm[img_i - cit_lo];
          pairs_v2l_arg.v2l.j_frc   = j_frc;
          pairs_v2l_arg.v2l.j_frc_v4 = j_frc_v4;
          pairs_v2l_arg.v2l.j_crd   = j_crd;
          pairs_v2l_arg.v2l.j_qterm = j_qterm;
          pairs_v2l_arg.v2l.j_iac   = j_iac;

          if (idbkt == ndbkt - 1) {
            jhi = iarr[ni - 1] + cit_lo;
            pairs_v2l_arg.v2l.iimg[0] = iarr[0] + cit_lo;
            pairs_v2l_arg.v2l.iimg[1] = iarr[1] + cit_lo;
            pairs_v2l_arg.v2l.iimg[2] = iarr[2] + cit_lo;
            pairs_v2l_arg.v2l.iimg[3] = iarr[3] + cit_lo;
          }
          pairs_v2l_arg.v2l.jlo   = jlo;
          pairs_v2l_arg.v2l.j_cnt = jhi - jlo;
          int j_cnt_eval = jhi - jlo;
          if (idbkt == ndbkt - 1){
            doublev4 jlti[5];
            doublev4 c0_v4 = 0.0;
            jlti[0] = simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
            jlti[1] = simd_set_doublev4(0.0, 2.0, 2.0, 2.0);
            jlti[2] = simd_set_doublev4(0.0, 0.0, 2.0, 2.0);
            jlti[3] = simd_set_doublev4(0.0, 0.0, 0.0, 2.0);
            jlti[4] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);
            int jlti_ptr = 0;
            for (img_j = iarr[0] + cit_lo - jlo; img_j < j_cnt; img_j ++){
              while (jlti_ptr < 4 && img_j + jlo >= iarr[jlti_ptr] + cit_lo) jlti_ptr ++;
              ij_unmask[img_j] = jlti[jlti_ptr];
            }
          }
          pairs_v2l_arg.v2l.ij_unmask = ij_unmask;
          //lwpf_stop(LDD_I);
          //lwpf_start(V2B);
          lwpf_start(J_LOOP);
          pairs_v2l_calc_func(&pairs_v2l_arg);
          lwpf_stop(J_LOOP);
          //lwpf_stop(V2B);
          //lwpf_start(UPD_FRC);
          if (idbkt == ndbkt - 1){
            for (img_j = iarr[0] + cit_lo - jlo; img_j < j_cnt; img_j ++){
              ij_unmask[img_j] = 2.0;
            }

            int ijoff = cit_lo - jlo;
            j_frc[iarr[0] + ijoff][0] += simd_vextf0(i_frc_v4[0]);
            j_frc[iarr[1] + ijoff][0] += simd_vextf1(i_frc_v4[0]);
            j_frc[iarr[2] + ijoff][0] += simd_vextf2(i_frc_v4[0]);
            j_frc[iarr[3] + ijoff][0] += simd_vextf3(i_frc_v4[0]);

            j_frc[iarr[0] + ijoff][1] += simd_vextf0(i_frc_v4[1]);
            j_frc[iarr[1] + ijoff][1] += simd_vextf1(i_frc_v4[1]);
            j_frc[iarr[2] + ijoff][1] += simd_vextf2(i_frc_v4[1]);
            j_frc[iarr[3] + ijoff][1] += simd_vextf3(i_frc_v4[1]);

            j_frc[iarr[0] + ijoff][2] += simd_vextf0(i_frc_v4[2]);
            j_frc[iarr[1] + ijoff][2] += simd_vextf1(i_frc_v4[2]);
            j_frc[iarr[2] + ijoff][2] += simd_vextf2(i_frc_v4[2]);
            j_frc[iarr[3] + ijoff][2] += simd_vextf3(i_frc_v4[2]);
          } else {
            i_frc[iarr[0]][0] += simd_vextf0(i_frc_v4[0]);
            i_frc[iarr[1]][0] += simd_vextf1(i_frc_v4[0]);
            i_frc[iarr[2]][0] += simd_vextf2(i_frc_v4[0]);
            i_frc[iarr[3]][0] += simd_vextf3(i_frc_v4[0]);

            i_frc[iarr[0]][1] += simd_vextf0(i_frc_v4[1]);
            i_frc[iarr[1]][1] += simd_vextf1(i_frc_v4[1]);
            i_frc[iarr[2]][1] += simd_vextf2(i_frc_v4[1]);
            i_frc[iarr[3]][1] += simd_vextf3(i_frc_v4[1]);

            i_frc[iarr[0]][2] += simd_vextf0(i_frc_v4[2]);
            i_frc[iarr[1]][2] += simd_vextf1(i_frc_v4[2]);
            i_frc[iarr[2]][2] += simd_vextf2(i_frc_v4[2]);
            i_frc[iarr[3]][2] += simd_vextf3(i_frc_v4[2]);
          }
          //lwpf_stop(V2B);
          //lwpf_stop(UPD_FRC);
        }

        //lwpf_stop(COMP);
#ifndef J_UPD_SCALAR
        for (img_j = 0; img_j < j_cnt; img_j += 4) {
          doublev4 frc_v4[4][3];
          doublev4 frc[3];
          transpose4x4(j_frc_v4[img_j][0], j_frc_v4[img_j + 1][0], j_frc_v4[img_j + 2][0], j_frc_v4[img_j + 3][0], frc_v4[0][0], frc_v4[1][0], frc_v4[2][0], frc_v4[3][0]);
          transpose4x4(j_frc_v4[img_j][1], j_frc_v4[img_j + 1][1], j_frc_v4[img_j + 2][1], j_frc_v4[img_j + 3][1], frc_v4[0][1], frc_v4[1][1], frc_v4[2][1], frc_v4[3][1]);
          transpose4x4(j_frc_v4[img_j][2], j_frc_v4[img_j + 1][2], j_frc_v4[img_j + 2][2], j_frc_v4[img_j + 3][2], frc_v4[0][2], frc_v4[1][2], frc_v4[2][2], frc_v4[3][2]);

          frc[0] = frc_v4[0][0] + frc_v4[1][0] + frc_v4[2][0] + frc_v4[3][0];
          frc[1] = frc_v4[0][1] + frc_v4[1][1] + frc_v4[2][1] + frc_v4[3][1];
          frc[2] = frc_v4[0][2] + frc_v4[1][2] + frc_v4[2][2] + frc_v4[3][2];
          j_frc[img_j + 0][0] += simd_vextf0(frc[0]);
          j_frc[img_j + 1][0] += simd_vextf1(frc[0]);
          j_frc[img_j + 2][0] += simd_vextf2(frc[0]);
          j_frc[img_j + 3][0] += simd_vextf3(frc[0]);


          j_frc[img_j + 0][1] += simd_vextf0(frc[1]);
          j_frc[img_j + 1][1] += simd_vextf1(frc[1]);
          j_frc[img_j + 2][1] += simd_vextf2(frc[1]);
          j_frc[img_j + 3][1] += simd_vextf3(frc[1]);



          j_frc[img_j + 0][2] += simd_vextf0(frc[2]);
          j_frc[img_j + 1][2] += simd_vextf1(frc[2]);
          j_frc[img_j + 2][2] += simd_vextf2(frc[2]);
          j_frc[img_j + 3][2] += simd_vextf3(frc[2]);

        }
#endif
      }
      //if (arg->mytaskid == 1 && jlo <= 197372 && irep == 1 && jhi > 197372) printf("write: %d %f %f %f\n", jcit, j_frc[o][0], j_frc[o][1], j_frc[o][2]);
      pe_put(img_frc_rep + jlo + irep * atm_cnt, j_frc, j_cnt * 3 * sizeof(double));
      dma_syn();
      //}

      //MAYBE SYNC HERE!!!
    }
  }
  gbl_accum_vars_t accum_scalar;
  simd_vsumd(accum.vxx   );
  simd_vsumd(accum.vxy   );
  simd_vsumd(accum.vxz   );
  simd_vsumd(accum.vyy   );
  simd_vsumd(accum.vyz   );
  simd_vsumd(accum.vzz   );
  simd_vsumd(accum.eedvir);
  simd_vsumd(accum.eed   );
  simd_vsumd(accum.evdw  );
  simd_vsumd(accum.ehb   );

  accum_scalar.vxx    = accum.vxx   ;
  accum_scalar.vxy    = accum.vxy   ;
  accum_scalar.vxz    = accum.vxz   ;
  accum_scalar.vyy    = accum.vyy   ;
  accum_scalar.vyz    = accum.vyz   ;
  accum_scalar.vzz    = accum.vzz   ;
  accum_scalar.eedvir = accum.eedvir;
  accum_scalar.eed    = accum.eed   ;
  accum_scalar.evdw   = accum.evdw  ;
  accum_scalar.ehb    = accum.ehb   ;

  //lwpf_stop(MAIN);
  pe_put(arg->accum + _MYID, &accum_scalar, sizeof(gbl_accum_vars_t));
  dma_syn();
  /* lwpf_stop(ALL); */
  /* lwpf_exit(NB_ENERGY); */
  /* } */
  /*   _MYID = 0; */
}

void get_nb_energy_reduction(nb_energy_arg_t *arg){
  dma_init();
  double (*img_crd)[3] = arg->img_vars.img_crd;
  int ncit = arg->cit_vars.cit_tbl_x_dim * arg->cit_vars.cit_tbl_y_dim * arg->cit_vars.cit_tbl_z_dim;
  int bkt_upd_cnt[ICIT_STEP];
  int atm_cnt = arg->img_vars.atm_cnt;
  cit_tbl_rec_t cits[ICIT_STEP];
  double (*img_frc_rep)[3] = arg->img_vars.img_frc_rep;
  double j_frc[MAX_BKT_SIZE][3], j_frc_rep[MAX_BKT_SIZE][3];
  int ist, ied, icit;
  for (ist = _MYID * ICIT_STEP; ist < ncit; ist += ICIT_STEP * 64){
    ied = ist + ICIT_STEP;
    if (ied > ncit) ied = ncit;
    pe_get(arg->cit_vars.bkt_upd_cnt + ist, bkt_upd_cnt, sizeof(int) * (ied - ist));
    pe_get(arg->cit_vars.cits        + ist, cits       , sizeof(cit_tbl_rec_t) * (ied - ist));
    dma_syn();
    for (icit = ist; icit < ied; icit ++){
      int nrep = bkt_upd_cnt[icit - ist];
      int lo = cits[icit - ist].img_lo - 1;
      int hi = cits[icit - ist].img_hi;
      int cnt = hi - lo;
      if (nrep > 1 && cnt > 0){
        pe_get(img_frc_rep + lo, j_frc, sizeof(double) * 3 * (hi - lo));
        int irep;
        for (irep = 1; irep < nrep; irep ++){
          pe_get(img_frc_rep + lo + irep * atm_cnt, j_frc_rep, sizeof(double) * 3 * (hi - lo));
          dma_syn();

          int img_j;
          for (img_j = 0; img_j < cnt; img_j ++){
            j_frc[img_j][0] += j_frc_rep[img_j][0];
            j_frc[img_j][1] += j_frc_rep[img_j][1];
            j_frc[img_j][2] += j_frc_rep[img_j][2];
          }
        }
        pe_put(img_frc_rep + lo, j_frc, sizeof(double) * 3 * (hi - lo));
        dma_syn();
      }
    }
  }
}

void get_nb_energy_reduction_excl(nb_energy_arg_t *gl_arg){
  dma_init();
  nb_energy_arg_t arg;
  pe_get(gl_arg, &arg, sizeof(nb_energy_arg_t));
  dma_syn();
  double (*img_crd)[3] = arg.img_vars.img_crd;
  int ncit = arg.cit_vars.cit_tbl_x_dim * arg.cit_vars.cit_tbl_y_dim * arg.cit_vars.cit_tbl_z_dim;
  int bkt_upd_cnt[ICIT_STEP];
  int atm_cnt = arg.img_vars.atm_cnt;
  cit_tbl_rec_t cits[ICIT_STEP];
  double (*img_frc_rep)[3] = arg.img_vars.img_frc_rep;
  double (*img_frc)[3] = arg.img_vars.img_frc;
  double j_frc[MAX_BKT_SIZE][3], j_frc_rep[MAX_BKT_SIZE][3];
  int ist, ied, icit;
  for (ist = _MYID * ICIT_STEP; ist < ncit; ist += ICIT_STEP * 64){
    ied = ist + ICIT_STEP;
    if (ied > ncit) ied = ncit;
    pe_get(arg.cit_vars.bkt_upd_cnt + ist, bkt_upd_cnt, sizeof(int) * (ied - ist));
    pe_get(arg.cit_vars.cits        + ist, cits       , sizeof(cit_tbl_rec_t) * (ied - ist));
    dma_syn();
    for (icit = ist; icit < ied; icit ++){
      int nrep = bkt_upd_cnt[icit - ist];
      if (nrep > 0){
        int lo = cits[icit - ist].img_lo - 1;
        int hi = cits[icit - ist].img_hi;
        if (lo < hi){
          int cnt = hi - lo;
          pe_get(img_frc + lo, j_frc, sizeof(double) * 3 * (hi - lo));
          pe_get(img_frc_rep + lo, j_frc_rep, sizeof(double) * 3 * (hi - lo));
          dma_syn();
          int img_j;
          for (img_j = 0; img_j < cnt; img_j ++){
            j_frc[img_j][0] += j_frc_rep[img_j][0];
            j_frc[img_j][1] += j_frc_rep[img_j][1];
            j_frc[img_j][2] += j_frc_rep[img_j][2];
          }
          pe_put(img_frc + lo, j_frc, sizeof(double) * 3 * (hi - lo));
          dma_syn();
        }
      }/*  else { */
      /*   int lo = cits[icit - ist].img_lo - 1; */
      /*   int hi = cits[icit - ist].img_hi; */
      /*   int cnt = hi - lo; */
      /*   int img_j; */
      /*   for (img_j = 0; img_j < cnt; img_j ++){ */
      /*     j_frc[img_j][0] = 0; */
      /*     j_frc[img_j][1] = 0; */
      /*     j_frc[img_j][2] = 0; */
      /*   } */
      /*   pe_put(img_frc_rep + lo + atm_cnt, j_frc, sizeof(double) * 3 * (hi - lo)); */
      /*   dma_syn(); */
      /* } */
    }
  }
}


void get_nb_energy_enterance(nb_energy_arg_t *gl_arg){
  lwpf_enter(NB_ENERGY);
  lwpf_start(ALL);
  dma_init();
  nb_energy_arg_t arg;
  pe_get(gl_arg, &arg, sizeof(nb_energy_arg_t));
  dma_syn();
  lwpf_start(BBOX);
  get_nb_energy_bbox_scan(&arg);
  lwpf_stop(BBOX);
  athread_syn(ARRAY_SCOPE, 0xffff);
  lwpf_start(MAIN);
  get_nb_energy_main_loop(&arg);
  lwpf_stop(MAIN);
  //int atm_cnt = arg.img_vars.atm_cnt;
  //int watch = atm_cnt + 197372;
  //double (*img_frc_rep)[3] = arg.img_vars.img_frc_rep;
  //i (arg.mytaskid == 1 && !_MYID) printf("after main: %f %f %f\n", img_frc_rep[watch][0], img_frc_rep[watch][1], img_frc_rep[watch][2]);
  athread_syn(ARRAY_SCOPE, 0xffff);
  lwpf_start(REDUCE_CPE);
  get_nb_energy_reduction(&arg);
  lwpf_stop(REDUCE_CPE);
  //if (arg.mytaskid == 1 && !_MYID) printf("after reduce: %f %f %f\n", img_frc_rep[watch][0], img_frc_rep[watch][1], img_frc_rep[watch][2]);

  lwpf_stop(ALL);
  lwpf_exit(NB_ENERGY);
}

void get_nb_energy_enterance2(nb_energy_arg_t *gl_arg){
  lwpf_enter(NB_ENERGY);
  lwpf_start(ALL);
  dma_init();
  nb_energy_arg_t arg;
  pe_get(gl_arg, &arg, sizeof(nb_energy_arg_t));
  dma_syn();
  lwpf_start(REDUCE_CPE);
  get_nb_energy_reduction(&arg);
  lwpf_stop(REDUCE_CPE);

  lwpf_stop(ALL);
  lwpf_exit(NB_ENERGY);
}

#endif
