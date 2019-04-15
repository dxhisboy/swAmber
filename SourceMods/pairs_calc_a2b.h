extern double erfc_wiki(double);
#define RSQRT_PI 0.56418958354775628
#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2b_efv_2cut_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_efv_2cut(pairs_calc_a2b_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_a2b_efv_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_efv(pairs_calc_a2b_arg_t *arg)
#endif /*FSWITCH*/
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2b_fv_2cut_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_fv_2cut(pairs_calc_a2b_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_a2b_fv_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_fv(pairs_calc_a2b_arg_t *arg)
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2b_f_2cut_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_f_2cut(pairs_calc_a2b_arg_t *arg)
#endif /* FSWITCH */
#else //NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2b_f_fs(pairs_calc_a2b_arg_t *arg)
#else
void pairs_calc_a2b_f(pairs_calc_a2b_arg_t *arg)
#endif //FSWITCH
#endif //CALC_NOVEC_2CUT
#endif //CALC_F
{
  //lwpf_start(A2B);
  int j_cnt   = arg->a2b.j_cnt;

  double dix = arg->a2b.bkt_tran[0] - arg->a2b.i_crd[0];
  double diy = arg->a2b.bkt_tran[1] - arg->a2b.i_crd[1];
  double diz = arg->a2b.bkt_tran[2] - arg->a2b.i_crd[2];
  double *i_frc = arg->a2b.i_frc;
  /* double bkt_tran[3]; */
  /* bkt_tran[0] = arg->a2b.bkt_tran[0]; */
  /* bkt_tran[1] = arg->a2b.bkt_tran[1]; */
  /* bkt_tran[2] = arg->a2b.bkt_tran[2]; */

  /* int ntypes   = arg->ico_tbl.ntypes; */
  /* int *ico     = arg->ico_tbl.ico   ; */
  /* double *cn1  = arg->ico_tbl.cn1   ; */
  /* double *cn2  = arg->ico_tbl.cn2   ; */
  /* double *asol = arg->ico_tbl.asol  ; */
  /* double *bsol = arg->ico_tbl.bsol  ; */

  int    ntypes = arg->packed_ico.ntypes;
  double *p1    = arg->packed_ico.p1    ;
  double *p2    = arg->packed_ico.p2    ;
  long   *icd   = (long*)arg->packed_ico.icd;
  double (*j_frc)[3] = arg->a2b.j_frc  ;
  double (*j_crd)[3] = arg->a2b.j_crd  ;
  double *j_qterm    = arg->a2b.j_qterm;
  int    *j_iac      = arg->a2b.j_iac  ;

  double max_nb_cut2     = arg->share_vars.max_nb_cut2    ;
  double es_cut2         = arg->share_vars.es_cut2        ;
  double p12             = arg->share_vars.p12            ;
  double p6              = arg->share_vars.p6             ;
  double dxdr            = arg->share_vars.dxdr           ;
  double cut             = arg->share_vars.cut            ;
  double cut6            = arg->share_vars.cut6           ;
  double cut3            = arg->share_vars.cut3           ;
  double cutinv          = arg->share_vars.cutinv         ;
  double cut2inv         = arg->share_vars.cut2inv        ;
  double cut3inv         = arg->share_vars.cut3inv        ;
  double cut6inv         = arg->share_vars.cut6inv        ;
  double fswitch2        = arg->share_vars.fswitch2       ;
  double fswitch3        = arg->share_vars.fswitch3       ;
  double fswitch6        = arg->share_vars.fswitch6       ;
  double invfswitch6cut6 = arg->share_vars.invfswitch6cut6;
  double invfswitch3cut3 = arg->share_vars.invfswitch3cut3;

  double eedvir = 0;
  double eed   = 0;
  double evdw  = 0;
  double ehb   = 0;

  double vxx   = 0;
  double vxy   = 0;
  double vxz   = 0;
  double vyy   = 0;
  double vyz   = 0;
  double vzz   = 0;

  double dumx = 0, dumy = 0, dumz = 0;
  //if (arg->mytaskid == 0) printf("%d %d\n", img_i, jhi);
  //if (j_cnt > 100) printf("%d\n", j_cnt);
  double cgi = arg->a2b.i_qterm;
  long   *i_icd = icd + (arg->a2b.i_iac - 1) * ntypes;
  double *i_p1  = p1  + (arg->a2b.i_iac - 1) * ntypes;
  double *i_p2  = p2  + (arg->a2b.i_iac - 1) * ntypes;
  int img_j;
  //lwpf_start(LOOP_J);
  for (img_j = 0; img_j < j_cnt; img_j ++){
    double delx = j_crd[img_j][0] + dix;
    double dely = j_crd[img_j][1] + diy;
    double delz = j_crd[img_j][2] + diz;
    double delr2 = delx * delx + dely * dely + delz * delz;
    double df = 0;
    //lwpf_start(CALC);
    if (delr2 < max_nb_cut2) {
      double delr = sqrt(delr2);
      double delrinv = 1 / delr;
      double delr2inv = delrinv * delrinv;
      //lwpf_start(EE);
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 < es_cut2){
#endif //BUILD_PAIRS_CALC_NOVEC_2CUT
        double cgi_cgj = cgi * j_qterm[img_j];
        double x = dxdr * delr;
        //lwpf_start(SWITCH);
#ifdef MPE
        double sw;
        derfcfun_(&x, &sw);
#endif
/* #ifdef CPE */
/*         double sw = erfc_wiki(x); */
/* #endif */
/*         double d_sw_dx = -(2*exp(-x*x)) * RSQRT_PI; */
        double sw, d_sw_dx;
        calc_sw(x, &sw, &d_sw_dx);
        //lwpf_stop(SWITCH);
        double b0 = cgi_cgj * delrinv * sw;
        double b1 = b0 - cgi_cgj * d_sw_dx * dxdr;
#ifdef NEED_ENE
        eed += b0;
#endif //NEED_ENE
#ifdef NEED_VIR
        eedvir += -b1;
#endif //NEED_VIR
        df += b1 * delr2inv;
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      }
#endif //BUILD_PAIRS_CALC_NOVEC_2CUT
      //lwpf_stop(EE);
      //lwpf_start(VDW);
      int jac = j_iac[img_j] - 1;
      long icd = i_icd[jac];
      if (icd != 0){
#ifdef HAS_10_12
        //if (ic > -1){
        if (icd > 0){
#endif //HAS_10_12
          double delr6inv = delr2inv * delr2inv * delr2inv;
          double delr12inv = delr6inv * delr6inv;
          double f12, f6;
#ifdef FSWITCH
          double delr3inv = delr2inv * delrinv;
          double df12, df6;
          if (delr2 > fswitch2){
            double delr6inv_cut6inv = delr6inv - cut6inv;
            double delr3inv_cut3inv = delr3inv - cut3inv;
            //if (img_i == 999) puts("sw");
            f12  = i_p1[jac] * p12 * delr6inv_cut6inv * delr6inv_cut6inv;
            f6   = i_p2[jac] * p6  * delr3inv_cut3inv * delr3inv_cut3inv;
            df12 = -12.0 * p12 * delr2inv * delr6inv * delr6inv_cut6inv;
            df6  =  -6.0 * p6  * delr3inv * delr2inv * delr3inv_cut3inv;
          } else {
            //if (img_i == 999) puts("nosw");
            f12  = i_p1[jac] * (delr12inv - invfswitch6cut6);
            f6   = i_p2[jac] * (delr6inv - invfswitch3cut3);
            df12 = -12.0 * delr2inv * delr12inv;
            df6  = -6.0 * delr2inv * delr6inv;
          }
#ifdef NEED_ENE
          evdw += f12 - f6;
#endif //NEED_ENE
          df    = df - i_p1[jac] * df12 + i_p2[jac] * df6;
#else //FSWITCH
          f6  = i_p2[jac] * delr6inv;
          f12 = i_p1[jac] * delr12inv;
#ifdef NEED_ENE
          evdw += f12 - f6;
#endif //NEED_ENE
          df += (12.0 * f12 - 6.0 * f6) * delr2inv;
#endif //FSWITCH
#ifdef HAS_10_12
        } else {
          ic = -1 - ic;
          double delr10inv = delr2inv * delr2inv * delr2inv * delr2inv * delr2inv;
          double delr12inv = delr10inv * delr2inv;
          double f10 = i_p1[jac] * delr10inv;
          double f12 = i_p2[jac] * delr12inv;
#ifdef NEED_ENE
          ehb += f12 - f10;
#endif //NEED_ENE
          df += (12.0 * f12 - 10.0 * f10) * delr2inv;
        }
#endif //HAS_10_12
      }
      //lwpf_stop(VDW);
      double dfx = delx * df;
      double dfy = dely * df;
      double dfz = delz * df;
      //if (img_i == 999) printf("full: %e %e %e %e\n", df, delx, dely, delz);
#ifdef NEED_VIR
      vxx -= delx * dfx;
      vxy -= delx * dfy;
      vxz -= delx * dfz;
      vyy -= dely * dfy;
      vyz -= dely * dfz;
      vzz -= delz * dfz;
#endif // NEED_VIR
      dumx += dfx;
      dumy += dfy;
      dumz += dfz;

      j_frc[img_j][0] += dfx;
      j_frc[img_j][1] += dfy;
      j_frc[img_j][2] += dfz;
    }
    //lwpf_stop(CALC);
  }
  //lwpf_stop(LOOP_J);
  i_frc[0] -= dumx;
  i_frc[1] -= dumy;
  i_frc[2] -= dumz;

  arg->accum->vxx    += vxx;
  arg->accum->vxy    += vxy;
  arg->accum->vxz    += vxz;
  arg->accum->vyy    += vyy;
  arg->accum->vyz    += vyz;
  arg->accum->vzz    += vzz;
  arg->accum->eedvir += eedvir;
  arg->accum->eed    += eed;
  arg->accum->evdw   += evdw;
  arg->accum->ehb    += ehb;
  //lwpf_stop(A2B);
}
#undef NEED_ENE
#undef NEED_VIR
