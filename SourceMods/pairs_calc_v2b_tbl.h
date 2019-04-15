extern void calc_sw_v4(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx);
extern void simd_calc_sw(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx);
#define RSQRT_PI 0.56418958354775628
#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_v2b_tbl_efv_2cut_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_efv_2cut(pairs_calc_v2b_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_v2b_tbl_efv_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_efv(pairs_calc_v2b_arg_t *arg)
#endif /*FSWITCH*/
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_v2b_tbl_fv_2cut_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_fv_2cut(pairs_calc_v2b_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_v2b_tbl_fv_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_fv(pairs_calc_v2b_arg_t *arg)
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_v2b_tbl_f_2cut_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_f_2cut(pairs_calc_v2b_arg_t *arg)
#endif /* FSWITCH */
#else //NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_v2b_tbl_f_fs(pairs_calc_v2b_arg_t *arg)
#else
void pairs_calc_v2b_tbl_f(pairs_calc_v2b_arg_t *arg)
#endif //FSWITCH
#endif //CALC_NOVEC_2CUT
#endif //CALC_F
{
  doublev4 c1_v4 = 1, c0_v4 = 0, c6_v4 = 6, c12_v4 = 12, c2_v4 = 2;
  //lwpf_start(V2B);
  int j_cnt   = arg->v2b.j_cnt;
  int jlo     = arg->v2b.jlo;
  doublev4 *i_frc = arg->v2b.i_frc;

  doublev4 *i_ic   = arg->v2b.i_icd ;
  doublev4 *i_p1   = arg->v2b.i_p1  ;
  doublev4 *i_p2   = arg->v2b.i_p2  ;
  doublev4 iunmask = arg->v2b.iunmask;
  int *iimg        = arg->v2b.iimg;
  //lwpf_start(EXT);
  doublev4 jlti[5];
  doublev4 dix = arg->v2b.bkt_tran[0] - arg->v2b.i_crd[0];
  doublev4 diy = arg->v2b.bkt_tran[1] - arg->v2b.i_crd[1];
  doublev4 diz = arg->v2b.bkt_tran[2] - arg->v2b.i_crd[2];

  jlti[0] = simd_vselne(iunmask, simd_set_doublev4(2.0, 2.0, 2.0, 2.0), c0_v4);
  jlti[1] = simd_vselne(iunmask, simd_set_doublev4(0.0, 2.0, 2.0, 2.0), c0_v4);
  jlti[2] = simd_vselne(iunmask, simd_set_doublev4(0.0, 0.0, 2.0, 2.0), c0_v4);
  jlti[3] = simd_vselne(iunmask, simd_set_doublev4(0.0, 0.0, 0.0, 2.0), c0_v4);
  jlti[4] = simd_vselne(iunmask, simd_set_doublev4(0.0, 0.0, 0.0, 0.0), c0_v4);
  doublev4 (*j_frc_v4)[3] = arg->v2b.j_frc_v4;
  double (*j_frc)[3] = arg->v2b.j_frc;
  double (*j_crd)[3] = arg->v2b.j_crd  ;
  double *j_qterm    = arg->v2b.j_qterm;
  int    *j_iac      = arg->v2b.j_iac  ;

  doublev4 max_nb_cut2     = arg->share_vars.max_nb_cut2    ;
  doublev4 es_cut2         = arg->share_vars.es_cut2        ;
  doublev4 p12             = arg->share_vars.p12            ;
  doublev4 p6              = arg->share_vars.p6             ;
  doublev4 dxdr            = arg->share_vars.dxdr           ;
  doublev4 cut3inv         = arg->share_vars.cut3inv        ;
  doublev4 cut6inv         = arg->share_vars.cut6inv        ;
  doublev4 fswitch2        = arg->share_vars.fswitch2       ;
  doublev4 invfswitch6cut6 = arg->share_vars.invfswitch6cut6;
  doublev4 invfswitch3cut3 = arg->share_vars.invfswitch3cut3;
  doublev4 lowest_efs_delr2 = arg->share_vars.lowest_efs_delr2;
  doublev4 efs_tbl_dens    = arg->share_vars.efs_tbl_dens;
  doublev4 del_efs         = arg->share_vars.del_efs;
  double *efs_tbl          = (double*)arg->share_vars.efs_tbl;
  //lwpf_stop(EXT);
  /* doublev4 eedvir = 0; */
  /* doublev4 eed   = 0; */
  /* doublev4 evdw  = 0; */
  /* doublev4 ehb   = 0; */

  /* doublev4 vxx   = 0; */
  /* doublev4 vxy   = 0; */
  /* doublev4 vxz   = 0; */
  /* doublev4 vyy   = 0; */
  /* doublev4 vyz   = 0; */
  /* doublev4 vzz   = 0; */
  gbl_accum_vars_v4_t *accum = arg->accum;
  doublev4 dumx = 0, dumy = 0, dumz = 0;
  doublev4 cgi = arg->v2b.i_qterm;

  int jlti_ptr = 0;
  int img_j;
  /* if (arg->mytaskid == 0 && _MYID == 0){ */
  /*   long sp; */
  /*   asm("addl $30, $31, %0\n\t" : "=r"(sp)); */
  /*   printf("%d\n", sp); */
  /* } */
  //lwpf_start(LOOP_J);
  for (img_j = 0; img_j < j_cnt; img_j ++){
    //lwpf_start(DELR);
    while (img_j + jlo >= iimg[jlti_ptr]) jlti_ptr ++;
    doublev4 unmask = jlti[jlti_ptr];
    doublev4 j_crd_v4[3];
    j_crd_v4[0] = j_crd[img_j][0];
    j_crd_v4[1] = j_crd[img_j][1];
    j_crd_v4[2] = j_crd[img_j][2];

    doublev4 delx_v4 = j_crd_v4[0] + dix;
    doublev4 dely_v4 = j_crd_v4[1] + diy;
    doublev4 delz_v4 = j_crd_v4[2] + diz;
    doublev4 delr2_v4 = delx_v4 * delx_v4 + dely_v4 * dely_v4 + delz_v4 * delz_v4;
    doublev4 df_v4 = 0;
    //lwpf_stop(DELR);
    //lwpf_start(CALC);
    doublev4 lt_max_nb_cut = simd_vfcmplt(delr2_v4, max_nb_cut2);
    lt_max_nb_cut = simd_vselne(unmask, lt_max_nb_cut, c0_v4);
    //lt_max_nb_cut = simd_vselne(iunmask, lt_max_nb_cut, c0_v4);
    int has_lt_max_nb_cut;
    vmatch(has_lt_max_nb_cut, lt_max_nb_cut, 0x40000000);
    if (has_lt_max_nb_cut) {
      //lwpf_start(INVSQRT);
      doublev4 delr_v4 = simd_vsqrtd(delr2_v4);
      doublev4 delrinv_v4 = c1_v4 / delr_v4;
      doublev4 delr2inv_v4 = delrinv_v4 * delrinv_v4;
      //lwpf_stop(INVSQRT);
      //lwpf_start(EE);
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      doublev4 lt_es_cut = simd_vfcmplt(delr2_v4, max_nb_cut2);
      int has_lt_es_cut;
      vmatch(has_lt_es_cut, lt_es_cut, 0x40000000);
      if (has_lt_es_cut){
#endif //BUILD_PAIRS_CALC_NOVEC_2CUT
        doublev4 cgj_v4;
        cgj_v4 = j_qterm[img_j];
        doublev4 cgi_cgj_v4 = simd_vselne(lt_max_nb_cut, cgi * cgj_v4, c0_v4);
        doublev4 lt_lowest_efs = simd_vfcmplt(delr2_v4, lowest_efs_delr2);
        int has_lt_lowest_efs;
        vmatch(has_lt_lowest_efs, lt_lowest_efs, 0x40000000);
        doublev4 b0_v4, b1_v4, dfee_v4;
        if (has_lt_lowest_efs) {
          doublev4 x_v4 = dxdr * delr_v4;
          doublev4 sw_v4, d_sw_dx_v4;
          //lwpf_start(SWITCH);
          //simd_calc_sw_inline_exp(x_v4, &sw_v4, &d_sw_dx_v4);
          simd_calc_sw_inline_exp(x_v4, &sw_v4, &d_sw_dx_v4);
          //lwpf_stop(SWITCH);
          b0_v4 = cgi_cgj_v4 * delrinv_v4 * sw_v4;
          b1_v4 = b0_v4 - cgi_cgj_v4 * d_sw_dx_v4 * dxdr;
          dfee_v4 = b1_v4 * delr2inv_v4;
        } else {
          long t0, t1, t2, t3, f0, f1, f2, f3;
          doublev4 e_tbl[4], f_tbl[4];
          doublev4 idx_v4 = simd_vselne(lt_max_nb_cut, delr2_v4 * efs_tbl_dens, c0_v4);
          doublev4 idx_floor_v4;
          asm ("vextf %9, 0, %0\n\t"
               "vextf %9, 1, %1\n\t"
               "vextf %9, 2, %2\n\t"
               "vextf %9, 3, %3\n\t"
               "fcvtdlr %0, 5, %0\n\t"
               "fcvtdlr %1, 5, %1\n\t"
               "fcvtdlr %2, 5, %2\n\t"
               "fcvtdlr %3, 5, %3\n\t"
               "fcvtld %0, %4\n\t"
               "fcvtld %1, %5\n\t"
               "fcvtld %2, %6\n\t"
               "fcvtld %3, %7\n\t"
               "vinsf %4, %8, 0, %8\n\t"
               "vinsf %5, %8, 1, %8\n\t"
               "vinsf %6, %8, 2, %8\n\t"
               "vinsf %7, %8, 3, %8\n\t"
               : "=&r"(t0), "=&r"(t1), "=&r"(t2), "=&r"(t3),
                 "=&r"(f0), "=&r"(f1), "=&r"(f2), "=&r"(f3),
                 "=&r"(idx_floor_v4)
               : "r"(idx_v4));
          doublev4 du = delr2_v4 - idx_floor_v4 * del_efs;
          doublev4 du2 = du * du;
          doublev4 du3 = du * du2;
          simd_load(f_tbl[0], efs_tbl + (t0 << 3) + 4);
          simd_load(f_tbl[1], efs_tbl + (t1 << 3) + 4);
          simd_load(f_tbl[2], efs_tbl + (t2 << 3) + 4);
          simd_load(f_tbl[3], efs_tbl + (t3 << 3) + 4);
          transpose4x4(f_tbl[0], f_tbl[1], f_tbl[2], f_tbl[3], f_tbl[0], f_tbl[1], f_tbl[2], f_tbl[3]);
          dfee_v4 = cgi_cgj_v4 * (f_tbl[0] + du * f_tbl[1] + du2 * f_tbl[2] + du3 * f_tbl[3]);
          /* doublev4 gt = simd_vfcmpgt(dfee_v4_tmp, dfee_v4 * 1.1); */
          /* doublev4 lt = simd_vfcmplt(dfee_v4_tmp, dfee_v4 * 0.9); */
          /* int has_gt, has_lt; */
          /* vmatch(has_gt, gt, 0x40000000); */
          /* vmatch(has_lt, lt, 0x40000000); */

#ifdef NEED_VIR
          b1_v4 = dfee_v4 * delr2_v4;
#endif
#ifdef NEED_ENE
          simd_load(e_tbl[0], efs_tbl + (t0 << 3));
          simd_load(e_tbl[1], efs_tbl + (t1 << 3));
          simd_load(e_tbl[2], efs_tbl + (t2 << 3));
          simd_load(e_tbl[3], efs_tbl + (t3 << 3));
          transpose4x4(e_tbl[0], e_tbl[1], e_tbl[2], e_tbl[3], e_tbl[0], e_tbl[1], e_tbl[2], e_tbl[3]);
          b0_v4 = cgi_cgj_v4 * (e_tbl[0] + du * e_tbl[1] + du2 * e_tbl[2] + du3 * e_tbl[3]);
#endif
        }
#ifdef NEED_ENE
        accum->eed += simd_vselne(lt_max_nb_cut, b0_v4, c0_v4);
#endif //NEED_ENE
#ifdef NEED_VIR
        accum->eedvir += simd_vselne(lt_max_nb_cut, -b1_v4, c0_v4);
#endif //NEED_VIR
        //doublev4 dfee_v4 = simd_vselne(lt_max_nb_cut, b1_v4 * delr2inv_v4, c0_v4);
        dfee_v4 = simd_vselne(lt_max_nb_cut, dfee_v4, c0_v4);
        doublev4 uno = simd_vfcmpun(dfee_v4, dfee_v4);
        int has_uno;
        vmatch(has_uno, uno, 0x40000000);
        if (has_uno){
          printf("%d %d %d %d %d\n", img_j + jlo, iimg[0], iimg[1], iimg[2], has_lt_lowest_efs);
          /* simd_print_doublev4(lt_max_nb_cut); */
          /* simd_print_doublev4(idx_v4); */
          /* simd_print_doublev4(idx_floor_v4); */
          /* simd_print_doublev4(delr2_v4); */
          /* simd_print_doublev4(du); */
          /* simd_print_doublev4(f_tbl[0]); */
          /* simd_print_doublev4(f_tbl[1]); */
          /* simd_print_doublev4(f_tbl[2]); */
          /* simd_print_doublev4(f_tbl[3]); */
          /* simd_print_doublev4(dfee_v4); */
          /* simd_print_doublev4(cgi_cgj_v4 * (f_tbl[0] + du * f_tbl[1] + du2 * f_tbl[2] + du3 * f_tbl[3])); */
          /* printf("%d %d %d %d\n", t0, t1, t2, t3); */
        }

        df_v4 += dfee_v4;
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      }
#endif //BUILD_PAIRS_CALC_NOVEC_2CUT
      //lwpf_stop(EE);
      //lwpf_start(VDW);
      int jac = j_iac[img_j] - 1;
      doublev4 ic_v4 = simd_vselne(lt_max_nb_cut, i_ic[jac], c0_v4);
      int has_ic_gt, has_ic_lt;
      vmatch(has_ic_gt, ic_v4, 0x40000000);
      vmatch(has_ic_lt, ic_v4, 0xc0000000);
      int has_ic_nz = has_ic_gt || has_ic_lt;
      if (has_ic_nz){
#ifdef HAS_10_12
        if (has_ic_gt){
          doublev4 vdw_cond = simd_vfcmpeq(c2_v4, ic_v4);
#else //HAS_10_12
          doublev4 vdw_cond = ic_v4; //simd_vselne(ic_v4, lt_max_nb_cut, c0_v4);
#endif
          doublev4 delr6inv_v4 = delr2inv_v4 * delr2inv_v4 * delr2inv_v4;
          doublev4 delr12inv_v4 = delr6inv_v4 * delr6inv_v4;
#ifdef FSWITCH
          doublev4 delr3inv_v4 = delr2inv_v4 * delrinv_v4;
          doublev4 delr6inv_cut6inv_v4 = delr6inv_v4 - cut6inv;
          doublev4 delr3inv_cut3inv_v4 = delr3inv_v4 - cut3inv;

          doublev4 sw_cond_v4 = delr2_v4 - fswitch2;
          doublev4 f12_sw_v4  = i_p1[jac] * p12 * delr6inv_cut6inv_v4 * delr6inv_cut6inv_v4;
          doublev4 f6_sw_v4   = i_p2[jac] * p6 * delr3inv_cut3inv_v4 * delr3inv_cut3inv_v4;
          doublev4 df12_sw_v4 = -c12_v4 * p12 * delr2inv_v4 * delr6inv_v4 * delr6inv_cut6inv_v4;
          doublev4 df6_sw_v4  = -c6_v4 * p6 * delr3inv_v4 * delr2inv_v4 * delr3inv_cut3inv_v4;

          doublev4 f12_nosw_v4  = i_p1[jac] * (delr12inv_v4 - invfswitch6cut6);
          doublev4 f6_nosw_v4   = i_p2[jac] * (delr6inv_v4 - invfswitch3cut3);
          doublev4 df12_nosw_v4 = -c12_v4 * delr2inv_v4 * delr12inv_v4;
          doublev4 df6_nosw_v4  = -c6_v4 * delr2inv_v4 * delr6inv_v4;

          doublev4 f12_v4 = simd_vselgt(sw_cond_v4, f12_sw_v4, f12_nosw_v4);
          doublev4 f6_v4  = simd_vselgt(sw_cond_v4, f6_sw_v4 , f6_nosw_v4 );
          doublev4 df12_v4 = simd_vselgt(sw_cond_v4, df12_sw_v4, df12_nosw_v4);
          doublev4 df6_v4  = simd_vselgt(sw_cond_v4, df6_sw_v4 , df6_nosw_v4 );

#ifdef NEED_ENE
          doublev4 e_vdw_v4 = simd_vselne(vdw_cond, f12_v4 - f6_v4, c0_v4);
          accum->evdw += e_vdw_v4;
#endif //NEED_ENE
          doublev4 df_vdw_v4 = simd_vselne(vdw_cond, df6_v4 * i_p2[jac] - df12_v4 * i_p1[jac], c0_v4);
          df_v4   += df_vdw_v4;
#else //FSWITCH
          doublev4 f6_v4 = i_p2[jac] * delr6inv_v4;
          doublev4 f12_v4 = i_p1[jac] * delr12inv_v4;
#ifdef NEED_ENE
          doublev4 e_vdw_v4 = simd_vselne(vdw_cond, f12_v4 - f6_v4, c0_v4);
          accum->evdw += e_vdw_v4;
#endif //NEED_ENE
          doublev4 df_vdw_v4 = simd_vselne(vdw_cond, (c12_v4 * f12_v4 - c6_v4 * f6_v4) * delr2inv_v4, c0_v4);
          df_v4 += df_vdw_v4; //(12.0 * f12 - 6.0 * f6) * delr2inv;
#endif //FSWITCH
#ifdef HAS_10_12
        }
        if (has_ic_lt) {
          doublev4 delr10inv_v4 = delr2inv_v4 * delr2inv_v4 * delr2inv_v4 * delr2inv_v4 * delr2inv_v4;
          doublev4 delr12inv_v4 = delr10inv_v4 * delr2inv_v4;
          doublev4 f10_v4 = i_p2[jac] * delr10inv_v4;
          doublev4 f12_v4 = i_p1[jac] * delr12inv_v4;
          doublev4 hb_cond = simd_vfcmpeq(-c2_v4, ic_v4);
#ifdef NEED_ENE
          accum->ehb += simd_vsellne(hb_cond, f12_v4 - f10_v4);
#endif //NEED_ENE
          doublev4 df_hb_v4 = simd_vselne(hb_cond, (c12_v4 * f12_v4 - c10_v4 * f10_v4) * delr2inv_v4, c0_v4);
          df_v4 += df_hb_v4;
        }
#endif //HAS_10_12
      }

      //lwpf_stop(VDW);
      //lwpf_start(REDUCE);
      doublev4 dfx_v4 = delx_v4 * df_v4;
      doublev4 dfy_v4 = dely_v4 * df_v4;
      doublev4 dfz_v4 = delz_v4 * df_v4;
#ifdef NEED_VIR
      accum->vxx -= delx_v4 * dfx_v4;
      accum->vxy -= delx_v4 * dfy_v4;
      accum->vxz -= delx_v4 * dfz_v4;
      accum->vyy -= dely_v4 * dfy_v4;
      accum->vyz -= dely_v4 * dfz_v4;
      accum->vzz -= delz_v4 * dfz_v4;
#endif // NEED_VIR
      dumx += dfx_v4;
      dumy += dfy_v4;
      dumz += dfz_v4;

      simd_vsumd(dfx_v4);
      simd_vsumd(dfy_v4);
      simd_vsumd(dfz_v4);
      j_frc[img_j][0] += dfx_v4;
      j_frc[img_j][1] += dfy_v4;
      j_frc[img_j][2] += dfz_v4;
      //lwpf_stop(REDUCE);
    }
    //lwpf_stop(CALC);
  }
  //lwpf_stop(LOOP_J);
  i_frc[0] -= dumx;
  i_frc[1] -= dumy;
  i_frc[2] -= dumz;

  //lwpf_stop(V2B);
}
#undef NEED_ENE
#undef NEED_VIR
