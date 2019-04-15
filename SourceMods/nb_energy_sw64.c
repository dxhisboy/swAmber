#include <simd.h>
#include <math.h>
#include <stdio.h>
#define RSQRT_PI 0.56418958354775628
#define MAX_BKT_SIZE 64
#define MAX_TYPE_SQ 1024
#define MAX_TYPE 32
#define dbg -197372

#include <nb_energy_struct.h>

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
  conf.pcr0 = PC0_N_ICACHE_ACCESS;
  conf.pcr1 = PC1_N_ICACHE_MISS;
  conf.pcr2 = PC2_N_ICACHE_MISS;
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
#define J_UPD_SCALAR
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
#define LWPF_KERNELS K(ALL) K(BBOX) K(MAIN) K(REDUCE_CPE) K(REDUCE_EXCL) K(V2B) K(PREF)
#include <lwpf2/lwpf2.h>
#define DMA_FAST
#include <dma_macros.h>
__thread_local int idbkt_dbg;
//#include <pairs_calc_wrapper.h>
#define SPARSE_CITS 0x18fff18c7fffffffL
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
  if (u[0] < 0) u[0] = 0;
  if (u[1] < 0) u[1] = 0;
  if (u[2] < 0) u[2] = 0;
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

#include <nb_energy_gen.h>

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
  //get_nb_energy_main_loop(&arg);
  get_nb_energy_funcs[arg.share_vars.ifunc](&arg);
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
