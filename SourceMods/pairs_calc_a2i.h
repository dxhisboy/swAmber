#define RSQRT_PI 0.56418958354775628
#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2i_efv_2cut_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_efv_2cut(pairs_calc_a2i_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_a2i_efv_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_efv(pairs_calc_a2i_arg_t *arg)
#endif /*FSWITCH*/
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2i_fv_2cut_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_fv_2cut(pairs_calc_a2i_arg_t *arg)
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
void pairs_calc_a2i_fv_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_fv(pairs_calc_a2i_arg_t *arg)
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2i_f_2cut_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_f_2cut(pairs_calc_a2i_arg_t *arg)
#endif /* FSWITCH */
#else //NOVEC_2CUT
#ifdef FSWITCH
void pairs_calc_a2i_f_fs(pairs_calc_a2i_arg_t *arg)
#else
void pairs_calc_a2i_f(pairs_calc_a2i_arg_t *arg)
#endif //FSWITCH
#endif //CALC_NOVEC_2CUT
#endif //CALC_F
{
  int img_i = arg->a2i.img_i - 1;
  int jlo   = arg->a2i.jlo - 1;
  int jhi   = arg->a2i.jhi;

  double bkt_tran[3];
  bkt_tran[0] = arg->a2i.bkt_tran[0];
  bkt_tran[1] = arg->a2i.bkt_tran[1];
  bkt_tran[2] = arg->a2i.bkt_tran[2];

  int ntypes   = arg->ico_tbl.ntypes;
  int *ico     = arg->ico_tbl.ico   ;
  double *cn1  = arg->ico_tbl.cn1   ;
  double *cn2  = arg->ico_tbl.cn2   ;
  double *asol = arg->ico_tbl.asol  ;
  double *bsol = arg->ico_tbl.bsol  ;

  double (*img_frc)[3] = arg->img_vars.img_frc  ;
  double (*img_crd)[3] = arg->img_vars.img_crd  ;
  double *img_qterm    = arg->img_vars.img_qterm;
  int    *img_iac      = arg->img_vars.img_iac  ;

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
  double dix = bkt_tran[0] - img_crd[img_i][0];
  double diy = bkt_tran[1] - img_crd[img_i][1];
  double diz = bkt_tran[2] - img_crd[img_i][2];
  //if (arg->mytaskid == 0) printf("%d %d\n", img_i, jhi);
  if (jhi - jlo > 100) printf("%d %d %d\n", img_i, jlo, jhi);
  double cgi = img_qterm[img_i];
  int iaci = ntypes * (img_iac[img_i] - 1);
  int nerf = 0;
  int img_j;
  for (img_j = jlo; img_j < jhi; img_j ++){
    double delx = img_crd[img_j][0] + dix;
    double dely = img_crd[img_j][1] + diy;
    double delz = img_crd[img_j][2] + diz;
    double delr2 = delx * delx + dely * dely + delz * delz;
    double df = 0;
    if (delr2 < max_nb_cut2) {
      double delr = sqrt(delr2);
      double delrinv = 1 / delr;
      double delr2inv = delrinv * delrinv;
      int ic = ico[iaci + img_iac[img_j] - 1] - 1;
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 < es_cut2){
#endif //BUILD_PAIRS_CALC_NOVEC_2CUT
        double cgi_cgj = cgi * img_qterm[img_j];
        double x = dxdr * delr;
        double sw;
        derfcfun_(&x, &sw);
        double d_sw_dx = -(2*exp(-x*x)) * RSQRT_PI;
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
      //if (img_i == 999) printf("ee: %d %e %e %e %e %d\n", img_j, df, img_crd[img_j][0], img_crd[img_j][1], img_crd[img_j][2], ic);
      if (ic != -1){
#ifdef HAS_10_12
        if (ic > -1){
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
            f12  = cn1[ic] * p12 * delr6inv_cut6inv * delr6inv_cut6inv;
            f6   = cn2[ic] * p6  * delr3inv_cut3inv * delr3inv_cut3inv;
            df12 = -12.0 * p12 * delr2inv * delr6inv * delr6inv_cut6inv;
            df6  =  -6.0 * p6  * delr3inv * delr2inv * delr3inv_cut3inv;
          } else {
            //if (img_i == 999) puts("nosw");
            f12  = cn1[ic] * (delr12inv - invfswitch6cut6);
            f6   = cn2[ic] * (delr6inv - invfswitch3cut3);
            df12 = -12.0 * delr2inv * delr12inv;
            df6  = -6.0 * delr2inv * delr6inv;
          }
#ifdef NEED_ENE
          evdw += f12 - f6;
#endif //NEED_ENE
          df    = df - cn1[ic] * df12 + cn2[ic] * df6;
#else //FSWITCH
          f6  = cn2[ic] * delr6inv;
          f12 = cn1[ic] * delr12inv;
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
          double f10 = bsol[ic] * delr10inv;
          double f12 = asol[ic] * delr12inv;
#ifdef NEED_ENE
          ehb += f12 - f10;
#endif //NEED_ENE
          df += (12.0 * f12 - 10.0 * f10) * delr2inv;
        }
#endif //HAS_10_12
      }
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

      img_frc[img_j][0] += dfx;
      img_frc[img_j][1] += dfy;
      img_frc[img_j][2] += dfz;
    }
  }
  img_frc[img_i][0] -= dumx;
  img_frc[img_i][1] -= dumy;
  img_frc[img_i][2] -= dumz;

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
}
#undef NEED_ENE
#undef NEED_VIR
