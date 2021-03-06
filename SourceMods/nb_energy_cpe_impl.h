#include <suffix.h>
#include <pairs_calc_v2b.h>
static double CAT(dist2_crd_box_new, SUFFIX)(double *crd, double *box_c, double *box_h){
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

static void CAT(copy_share_vars_to_vec, SUFFIX)(share_vars_t *in, share_vars_v4_t *out){
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

static void CAT(cpe_vec_zfill_dbl, SUFFIX)(double *arr, int ndbl){
  int i;
  for (i = 0; i < ndbl; i += 16){
    asm ("vstd $31, 0(%0)\n\t"
         "vstd $31, 32(%0)\n\t"
         "vstd $31, 64(%0)\n\t"
         "vstd $31, 96(%0)\n\t"
         ::"r"(arr + i) : "memory");
  }
}

static void CAT(cpe_vec_add_dbl, SUFFIX)(double *dst, double *src, int ndbl){
  int i;
  float t0, t1, t2, t3, t4, t5, t6, t7;
  for (i = 0; i < ndbl; i += 16){
    asm ("vldd %2,  0(%0)\n\t"
         "vldd %3, 32(%0)\n\t"
         "vldd %4, 64(%0)\n\t"
         "vldd %5, 96(%0)\n\t"
         "vldd %6,  0(%1)\n\t"
         "vldd %7, 32(%1)\n\t"
         "vldd %8, 64(%1)\n\t"
         "vldd %9, 96(%1)\n\t"
         "vaddd %2, %6, %2\n\t"
         "vaddd %3, %7, %3\n\t"
         "vaddd %4, %8, %4\n\t"
         "vaddd %5, %9, %5\n\t"
         "vstd %2,  0(%0)\n\t"
         "vstd %3, 32(%0)\n\t"
         "vstd %4, 64(%0)\n\t"
         "vstd %5, 96(%0)\n\t"
         :: "r"(dst + i), "r"(src + i)
          , "r"(t0), "r"(t1), "r"(t2), "r"(t3), "r"(t4), "r"(t5), "r"(t6), "r"(t7) : "memory");
  }
}

static int CAT(pack_jatoms, SUFFIX)(int j_cnt, double *ibox_c, double *ibox_h, double (*j_crd)[3], double *j_qterm, int *j_iac, int *j_lst, double max_nb_cut2){
  //lwpf_start(J_CHK);
  /* if (SPARSE_CITS & (1L << idbkt)){ */
  int img_j, j_nlst = 0;
  for (img_j = 0; img_j < j_cnt; img_j ++){
    if (CAT(dist2_crd_box_new, SUFFIX)(j_crd[img_j], ibox_c, ibox_h) < max_nb_cut2){
      j_lst[j_nlst] = img_j;
      j_crd[j_nlst][0] = j_crd[img_j][0];
      j_crd[j_nlst][1] = j_crd[img_j][1];
      j_crd[j_nlst][2] = j_crd[img_j][2];
      j_qterm[j_nlst] = j_qterm[img_j];
      j_iac[j_nlst] = j_iac[img_j];
      j_nlst ++;
    }
  }

  //lwpf_stop(J_CHK);
  return j_nlst;
}
void CAT(parameters_prefetch8, SUFFIX)(int ntypes, double *p, doublev4 *p_v4, int iac0, int iac1, int iac2, int iac3) {
  int jac;
  //lwpf_start(PREF);
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
  //lwpf_stop(PREF);
}
void CAT(get_nb_energy, SUFFIX)(nb_energy_arg_t *arg)
{
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
  doublev4 jlti[5];
  jlti[0] = simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
  jlti[1] = simd_set_doublev4(0.0, 2.0, 2.0, 2.0);
  jlti[2] = simd_set_doublev4(0.0, 0.0, 2.0, 2.0);
  jlti[3] = simd_set_doublev4(0.0, 0.0, 0.0, 2.0);
  jlti[4] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);
  
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
  pairs_calc_v2b_arg_t pairs_v2b_arg;
  gbl_accum_vars_v4_t accum;
  //memset(&accum, 0, sizeof(gbl_accum_vars_v4_t));
  doublev4 c0_v4 = 0.0, c2_v4 = 2.0;
  accum.eedvir = c0_v4;
  accum.eed    = c0_v4;
  accum.evdw   = c0_v4;
  accum.ehb    = c0_v4;
  accum.vxx    = c0_v4;
  accum.vxy    = c0_v4;
  accum.vxz    = c0_v4;
  accum.vyy    = c0_v4;
  accum.vyz    = c0_v4;
  accum.vzz    = c0_v4;

  CAT(copy_share_vars_to_vec, SUFFIX)(&(arg->share_vars), &(pairs_v2b_arg.share_vars));
  pairs_v2b_arg.img_vars = arg->img_vars;
  pairs_v2b_arg.accum    = &accum;
  pairs_v2b_arg.mytaskid = arg->mytaskid;

  double tranvec[18][3];
  int itran;
  for (itran = 0; itran < 18; itran ++){
    tranvec[itran][0] = arg->img_vars.tranvec[itran][0];
    tranvec[itran][1] = arg->img_vars.tranvec[itran][1];
    tranvec[itran][2] = arg->img_vars.tranvec[itran][2];
  }

  //void (*pairs_v2b_calc_func)(pairs_calc_v2b_arg_t *) = pairs_calc_v2b_funcs[arg->share_vars.ifunc];
  int icit;
  double j_crd[MAX_BKT_SIZE][3], j_qterm[MAX_BKT_SIZE], j_frc[MAX_BKT_SIZE][3], j_frc_acc[MAX_BKT_SIZE][3];
  //doublev4 j_frc_v4[MAX_BKT_SIZE][3];
  doublev4 i_vec[MAX_BKT_SIZE];
  doublev4 ij_unmask[MAX_BKT_SIZE];
  int j_lst[MAX_BKT_SIZE], j_nlst;
  int img_j;
  for (img_j = 0; img_j < MAX_BKT_SIZE; img_j ++) {
    ij_unmask[img_j] = 2.0;
  }

  pairs_v2b_arg.v2b.ij_unmask = ij_unmask;
  int j_iac[MAX_BKT_SIZE];
  double i_frc[MAX_BKT_SIZE][3], i_crd[MAX_BKT_SIZE][3], i_qterm[MAX_BKT_SIZE];
  long i_inrange[MAX_BKT_SIZE];
  int i_iac[MAX_BKT_SIZE];
  //lwpf_start(MAIN);
  int pe_bkt_lo, pe_bkt_hi;

  get_task_range(my_bkt_hi, my_bkt_lo, &pe_bkt_hi, &pe_bkt_lo);
  doublev4 i_unmask[5];
  i_unmask[0] = simd_set_doublev4(0.0, 0.0, 0.0, 0.0);
  i_unmask[1] = simd_set_doublev4(2.0, 0.0, 0.0, 0.0);
  i_unmask[2] = simd_set_doublev4(2.0, 2.0, 0.0, 0.0);
  i_unmask[3] = simd_set_doublev4(2.0, 2.0, 2.0, 0.0);
  i_unmask[4] = simd_set_doublev4(2.0, 2.0, 2.0, 2.0);
  double max_nb_cut2 = arg->share_vars.max_nb_cut2;
  long need_pack_mask = 0x47710477fdc77fL;

  doublev4 i_crd_v4[3], i_frc_v4[3];
  pairs_v2b_arg.v2b.i_p1  = i_p1;
  pairs_v2b_arg.v2b.i_p2  = i_p2;
  pairs_v2b_arg.v2b.i_icd = i_icd;
  pairs_v2b_arg.v2b.i_frc   = i_frc_v4;
  pairs_v2b_arg.v2b.i_crd   = i_crd_v4;

  pairs_v2b_arg.v2b.j_frc   = j_frc;
  pairs_v2b_arg.v2b.j_crd   = j_crd;
  pairs_v2b_arg.v2b.j_qterm = j_qterm;
  pairs_v2b_arg.v2b.j_iac   = j_iac;
  //pairs_v2b_arg.v2b.jlo   = jlo;

  //int cnt_calc = 0, cnt_calc_group = 0;
  for (icit = my_bkt_lo + _MYID; icit <= my_bkt_hi; icit += 64){
    //for (icit = pe_bkt_lo; icit < pe_bkt_hi; icit ++){
    //for (icit = my_bkt_lo + _MYID; icit <= my_bkt_hi; icit += 64){
    neigh_cit_entry_t neigh_cits[64];
    pe_get(gbl_neigh_cits[icit - my_bkt_lo], neigh_cits, sizeof(neigh_cit_entry_t) * ndbkt);
    dma_syn();
    int cit_lo = neigh_cits[ndbkt - 1].img_lo - 1;
    int cit_hi = neigh_cits[ndbkt - 1].img_hi;
    int ilo = cit_lo;
    int ihi = cit_hi;
    //assert(i_cnt <= MAX_BKT_SIZE);
    if (my_img_lo > cit_lo) ilo = my_img_lo;
    if (my_img_hi < cit_hi) ihi = my_img_hi;
    int i_cnt = ihi - ilo;
    ilo = ilo - cit_lo;
    ihi = ihi - cit_lo;
    int cit_cnt = cit_hi - cit_lo;
    //lwpf_start(GET_I);
    if (i_cnt <= 0) continue;
    doublev4 ibox_c_v4, ibox_h_v4;
    double ibox_c[4], ibox_h[4];
    pe_get(img_crd   + cit_lo, i_crd  , cit_cnt * 3 * sizeof(double));
    pe_get(img_qterm + cit_lo, i_qterm, cit_cnt * 1 * sizeof(double));
    pe_get(img_iac   + cit_lo, i_iac  , cit_cnt * 1 * sizeof(int)   );
    simd_load(ibox_c_v4, bbox[icit]);
    simd_load(ibox_h_v4, bbox[icit] + 4);

    ibox_c[0] = simd_vextf0(ibox_c_v4);
    ibox_c[1] = simd_vextf1(ibox_c_v4);
    ibox_c[2] = simd_vextf2(ibox_c_v4);
    
    ibox_h[0] = simd_vextf0(ibox_h_v4);
    ibox_h[1] = simd_vextf1(ibox_h_v4);
    ibox_h[2] = simd_vextf2(ibox_h_v4);

    int img_i;
    CAT(cpe_vec_zfill_dbl, SUFFIX)(i_frc[0], cit_cnt * 3);
    dma_syn();

    for (img_i = ilo; img_i < ihi; img_i ++){
      i_vec[img_i] = simd_set_doublev4(i_crd[img_i][0], i_crd[img_i][1], i_crd[img_i][2], i_qterm[img_i]);
    }
    //lwpf_stop(GET_I);

    //int idbkt;

    for (idbkt = 0; idbkt < ndbkt; idbkt ++){
      //if ((SPARSE_CITS & 1L << idbkt)) continue;
      int jhi = neigh_cits[idbkt].img_hi;
      int jlo = neigh_cits[idbkt].img_lo - 1;
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
      CAT(cpe_vec_zfill_dbl, SUFFIX)(j_frc[0], j_cnt * 3);
      /* for (img_j = 0; img_j < j_cnt; img_j ++){ */
      /*   j_frc[img_j][0] = 0; */
      /*   j_frc[img_j][1] = 0; */
      /*   j_frc[img_j][2] = 0; */
      /* } */

      if (first) {
        CAT(cpe_vec_zfill_dbl, SUFFIX)(j_frc_acc[0], j_cnt * 3);
        /* for (img_j = 0; img_j < j_cnt; img_j ++){ */
        /*   j_frc_acc[img_j][0] = 0; */
        /*   j_frc_acc[img_j][1] = 0; */
        /*   j_frc_acc[img_j][2] = 0; */
        /* } */
      } else {
        pe_get(img_frc_rep + jlo + irep * atm_cnt, j_frc_acc , j_cnt * 3 * sizeof(double));
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
      i_iac[cit_hi - cit_lo + 0] = 1;
      i_iac[cit_hi - cit_lo + 1] = 1;
      i_iac[cit_hi - cit_lo + 2] = 1;
      i_iac[cit_hi - cit_lo + 3] = 1;

      pairs_v2b_arg.v2b.bkt_tran[0] = tranvec[itran][0];
      pairs_v2b_arg.v2b.bkt_tran[1] = tranvec[itran][1];
      pairs_v2b_arg.v2b.bkt_tran[2] = tranvec[itran][2];

      dma_syn();


      //scan_bbox(j_cnt, j_crd, bbox_c, bbox_h);
      jbox_c[0] = jbox_c[0] + tranvec[itran][0];
      jbox_c[1] = jbox_c[1] + tranvec[itran][1];
      jbox_c[2] = jbox_c[2] + tranvec[itran][2];
      //lwpf_stop(SCAN);
      //lwpf_stop(GET_J);
      //Last bucket, update myself.

      if (idbkt == ndbkt - 1){
        CAT(cpe_vec_add_dbl, SUFFIX)(j_frc[0], i_frc[0], cit_cnt * 3);
      }
      img_i = ilo;
      while (img_i < ihi) {
        /* doublev4 i_frc_v4[3]; */
        /* doublev4 i_crd_v4[3]; */
        doublev4 i_qterm_v4;
        int ni = 0;
        //intv8 iv8 = cit_hi - cit_lo;
        int iarr[8];
        //simd_store(iv8, iarr);
        iarr[4] = MAX_BKT_SIZE;
        //lwpf_start(V2B);
        if (SPARSE_CITS & 1L << idbkt) {
          //if () {
          while (ni < 4 && img_i < ihi){
            double dist2 = CAT(dist2_crd_box_new, SUFFIX)(i_crd[img_i], jbox_c, jbox_h);
            if (dist2 < max_nb_cut2) {
              iarr[ni ++] = img_i;
            }
            img_i ++;
          } 
        } else {
          iarr[0] = img_i + 0;
          iarr[1] = img_i + 1;
          iarr[2] = img_i + 2;
          iarr[3] = img_i + 3;
          ni = ihi - img_i;
          if (ni > 4) ni = 4;
          img_i += 4;
        }
        if (ni == 0) {
          continue;
        }
        doublev4 iunmask = i_unmask[ni];

        while (ni < 4) {
          iarr[ni] = cit_hi - cit_lo; //iarr[ni - 1];
          ni ++;
        }
        //lwpf_start(PREF);

        if (ntypes <= 8){
          CAT(parameters_prefetch8, SUFFIX)(ntypes, p1, i_p1, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          CAT(parameters_prefetch8, SUFFIX)(ntypes, p2, i_p2, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          CAT(parameters_prefetch8, SUFFIX)(ntypes, icd, i_icd, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
        } else {
          parameters_prefetch_32(p1, i_p1, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          parameters_prefetch_32(p2, i_p2, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
          parameters_prefetch_32(icd, i_icd, i_iac[iarr[0]], i_iac[iarr[1]], i_iac[iarr[2]], i_iac[iarr[3]]);
        }
        //lwpf_stop(PREF);

        //lwpf_start(LDD_I);
        i_frc_v4[0] = 0.0;
        i_frc_v4[1] = 0.0;
        i_frc_v4[2] = 0.0;
        transpose4x4(i_vec[iarr[0]], i_vec[iarr[1]], i_vec[iarr[2]], i_vec[iarr[3]], i_crd_v4[0], i_crd_v4[1], i_crd_v4[2], i_qterm_v4);
        pairs_v2b_arg.v2b.iimg[0] = iarr[0] + cit_lo;
        pairs_v2b_arg.v2b.iimg[1] = iarr[1] + cit_lo;
        pairs_v2b_arg.v2b.iimg[2] = iarr[2] + cit_lo;
        pairs_v2b_arg.v2b.iimg[3] = iarr[3] + cit_lo;


        int j_cnt_eval = jhi - jlo;

        /* if (SPARSE_CITS & 1L << idbkt) { */
        /*   j_cnt_eval = j_nlst; */
        /* } */
        if (idbkt == ndbkt - 1){
          jhi = iarr[3] + cit_lo;
          j_cnt_eval = jhi - jlo;
          int jlti_ptr = 0;
          for (img_j = iarr[0]; img_j < j_cnt_eval; img_j ++){
            while (img_j + jlo >= iarr[jlti_ptr] + cit_lo) jlti_ptr ++;
            ij_unmask[img_j] = jlti[jlti_ptr];
          }
        }
        pairs_v2b_arg.v2b.iunmask = iunmask; //i_qterm[img_i - cit_lo];
        pairs_v2b_arg.v2b.j_cnt = j_cnt_eval;
        pairs_v2b_arg.v2b.jlo = jlo;
        pairs_v2b_arg.v2b.i_qterm = i_qterm_v4;
        //lwpf_stop(LDD_I);
        //lwpf_start(V2B);
        //pairs_v2b_calc_func(&pairs_v2b_arg);
        CAT(pairs_calc_v2b, SUFFIX)(&pairs_v2b_arg);
        //lwpf_stop(V2B);
        //lwpf_start(UPD_FRC);
        /* double (*i_frc_wb)[3] = i_frc; */
        if (idbkt == ndbkt - 1){
          for (img_j = iarr[0]; img_j < j_cnt_eval; img_j ++){
            ij_unmask[img_j] = c2_v4;
          }
          /*   i_frc_wb = j_frc; */
          /* } */

          /* i_frc_wb[iarr[0]][0] += simd_vextf0(i_frc_v4[0]); */
          /* i_frc_wb[iarr[1]][0] += simd_vextf1(i_frc_v4[0]); */
          /* i_frc_wb[iarr[2]][0] += simd_vextf2(i_frc_v4[0]); */
          /* i_frc_wb[iarr[3]][0] += simd_vextf3(i_frc_v4[0]); */

          /* i_frc_wb[iarr[0]][1] += simd_vextf0(i_frc_v4[1]); */
          /* i_frc_wb[iarr[1]][1] += simd_vextf1(i_frc_v4[1]); */
          /* i_frc_wb[iarr[2]][1] += simd_vextf2(i_frc_v4[1]); */
          /* i_frc_wb[iarr[3]][1] += simd_vextf3(i_frc_v4[1]); */

          /* i_frc_wb[iarr[0]][2] += simd_vextf0(i_frc_v4[2]); */
          /* i_frc_wb[iarr[1]][2] += simd_vextf1(i_frc_v4[2]); */
          /* i_frc_wb[iarr[2]][2] += simd_vextf2(i_frc_v4[2]); */
          /* i_frc_wb[iarr[3]][2] += simd_vextf3(i_frc_v4[2]); */

          j_frc[iarr[0]][0] += simd_vextf0(i_frc_v4[0]);
          j_frc[iarr[1]][0] += simd_vextf1(i_frc_v4[0]);
          j_frc[iarr[2]][0] += simd_vextf2(i_frc_v4[0]);
          j_frc[iarr[3]][0] += simd_vextf3(i_frc_v4[0]);

          j_frc[iarr[0]][1] += simd_vextf0(i_frc_v4[1]);
          j_frc[iarr[1]][1] += simd_vextf1(i_frc_v4[1]);
          j_frc[iarr[2]][1] += simd_vextf2(i_frc_v4[1]);
          j_frc[iarr[3]][1] += simd_vextf3(i_frc_v4[1]);

          j_frc[iarr[0]][2] += simd_vextf0(i_frc_v4[2]);
          j_frc[iarr[1]][2] += simd_vextf1(i_frc_v4[2]);
          j_frc[iarr[2]][2] += simd_vextf2(i_frc_v4[2]);
          j_frc[iarr[3]][2] += simd_vextf3(i_frc_v4[2]);
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
      }
      /* if (SPARSE_CITS & (1L << idbkt)){ */
      /*   int j; */
      /*   for (j = 0; j < j_nlst; j ++){ */
      /*   } */
      /* } */
      /* for (img_j = 0; img_j < j_cnt; img_j ++){ */
      /*   j_frc_acc[img_j][0] += j_frc[img_j][0]; */
      /*   j_frc_acc[img_j][1] += j_frc[img_j][1]; */
      /*   j_frc_acc[img_j][2] += j_frc[img_j][2]; */
      /* } */
      CAT(cpe_vec_add_dbl, SUFFIX)(j_frc_acc[0], j_frc[0], j_cnt * 3);

      pe_put(img_frc_rep + jlo + irep * atm_cnt, j_frc_acc, j_cnt * 3 * sizeof(double));
      dma_syn();
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

  pe_put(arg->accum + _MYID, &accum_scalar, sizeof(gbl_accum_vars_t));
  dma_syn();
}

#undef SUFFIX
