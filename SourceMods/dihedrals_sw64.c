#include <param.h>

#ifdef MPE
#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include <mpi.h>
#include <math.h>

extern SLAVE_FUN(get_dihedrals_energy_c_para)();
extern SLAVE_FUN(dihedrals_pack_c_para)();

extern SLAVE_FUN(pack_atom_slave)();
extern SLAVE_FUN(add_frc_slave)();
extern SLAVE_FUN(memset_slave)();
extern SLAVE_FUN(memset_cpe)();
static int r = 0;

#define LWPF_UNITS U(BOND_PACK)
#include <lwpf2/lwpf2.h>

GBFloat sign(GBFloat x, GBFloat y)
{
  if(y > 0)
    return fabs(x);
  else if(y < 0)
    return -fabs(x);
  else
    return 0;
}

void ti_update_ene(double pot_ene, 
                  int ene_type, 
                  int region,
                  double *ti_ene,
                  double *ti_ene_delta,
                  double *ti_signs,
                  double *ti_weights,
                  int ifmbar_lcl,
                  int do_mbar,
                  int bar_states,
                  double *bar_cont,
                  double *bar_lambda,
                  int si_dvdl)
{
  int i;
  double dvdl_temp;

  dvdl_temp = pot_ene * ti_signs[region];
  ti_ene[si_dvdl * 2 + 0] = ti_ene[si_dvdl * 2 + 0] + dvdl_temp;
  ti_ene_delta[ene_type] = ti_ene_delta[ene_type] + dvdl_temp;

  if (ifmbar_lcl != 0 && do_mbar)
  {
    for(i = 1; i <=  bar_states; i++)
    {
      bar_cont[i] = bar_cont[i] 
      + pot_ene * (bar_lambda[i * 2 + region-1] 
      - ti_weights[region]);
    }
  }
  return;
}



void get_dihedrals_energy_c_(int *dihed_cnt_, 
                            dihed_rec_c *dihed_,
                            /*dihed_crd_c *dihed_crd_,*/
                            double *x_,
                            double *frc_,
                            double *ep_,
                            double *gbl_pk_,
                            double *gbl_pn_,
                            int *gbl_ipn_,
                            double *gbl_gamc_,
                            double *gbl_gams_,
                            int *ti_mode_,
                            int *ti_region_,
                            int *ti_lscale_,
                            double *ti_signs_,
                            int *ti_lst_,
                            int *ti_sc_lst_,
                            double *ti_ene_,
                            double *ti_ene_delta_,
                            double *ti_weights_,
                            int *gbl_atm_owner_map_,
                            int *mytaskid_,
                            int *si_dihedral_ene_,
                            int *ifmbar_lcl_,
                            int *do_mbar_,
                            int *bar_states_,
                            double *bar_cont_,
                            double *bar_lambda_,
                            int *si_dvdl_,
                            int *atm_cnt_
                            )
{
  int dihed_cnt           = *dihed_cnt_;
  int atm_cnt             = *atm_cnt_;
  double *ep              = ep_;
  int do_mbar             = *do_mbar_;
  int bar_states          = *bar_states_;
  int ifmbar_lcl          = *ifmbar_lcl_;
  int mytaskid            = *mytaskid_;
  int ti_mode             = *ti_mode_;
  int ti_region           = *ti_region_;
  int ti_lscale           = *ti_lscale_;
  int si_dihedral_ene     = *si_dihedral_ene_;
  int si_dvdl             = *si_dvdl_;

  dihed_rec_c *dihed      = dihed_ - 1;
  //dihed_crd_c *dihed_crd  = dihed_crd_ - 1;
  double *x               = x_ - 3;
  double *frc             = frc_ - 3;
  double *gbl_pk          = gbl_pk_ - 1;
  double *gbl_pn          = gbl_pn_ - 1;
  double *gbl_gamc        = gbl_gamc_ -1;
  double *gbl_gams        = gbl_gams_ -1;
  double *ti_signs        = ti_signs_ -1;
  double *ti_ene          = ti_ene_ -2;
  double *ti_ene_delta    = ti_ene_delta_ -1;
  double *ti_weights      = ti_weights_ - 1;
  double *bar_cont        = bar_cont_ -1;
  double *bar_lambda      = bar_lambda_ -2;
  
  int *gbl_ipn            = gbl_ipn_ - 1;
  int *ti_lst             = ti_lst_ -3;
  int *ti_sc_lst          = ti_sc_lst_ -1;
  int *gbl_atm_owner_map  = gbl_atm_owner_map_ -1;

#ifdef SPDP
  float gmul[10] = {0.0, 2.0, 0.0, 4.0, 0.0, 6.0, 0.0, 8.0, 0.0, 10.0};
  float tm24   = 1.0E-18;
  float tm06   = 1.0E-06;
  float tenm3  = 1.0E-03;
  float zero   = 0.0;
  float one    = 1.0;
#else
  double gmul[10] = {0.0, 2.0, 0.0, 4.0, 0.0, 6.0, 0.0, 8.0, 0.0, 10.0};
  double tm24   = 1.E-18;
  double tm06   = 1.E-06;
  double tenm3  = 1.E-03;
  double zero   = 0.0;
  double one    = 1.0;
#endif

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  double time1,time2,time3,time4;
  double time5,time6,time7,time8;
 
  time1 = MPI_Wtime(); 
  int i, j;
  ///////packing///////
  dihed_atm_c *atm_in;
  atm_in = (dihed_atm_c *)malloc(sizeof(dihed_atm_c)*(atm_cnt + 64));
  dihed_frc_c *dfrc;
  dfrc = (dihed_frc_c *)malloc(sizeof(dihed_frc_c) * (dihed_cnt+10));
  
  dihed_crd_c *dihed_crd;
  dihed_crd = (dihed_crd_c *)malloc(sizeof(dihed_crd_c) * (dihed_cnt + 64));

  dihed_param_c pm;
  pm.dihed_cnt          = dihed_cnt         ;
  pm.atm_cnt            = atm_cnt           ;
  pm.ep                 = ep                ;
  pm.do_mbar            = do_mbar           ;
  pm.bar_states         = bar_states        ;
  pm.ifmbar_lcl         = ifmbar_lcl        ;
  pm.mytaskid           = mytaskid          ;
  pm.ti_mode            = ti_mode           ;
  pm.ti_region          = ti_region         ;
  pm.ti_lscale          = ti_lscale         ;
  pm.si_dihedral_ene    = si_dihedral_ene   ;
  pm.si_dvdl            = si_dvdl           ;
                           
  pm.dihed              = dihed             ;
  pm.dihed_crd          = dihed_crd             ;//
  pm.dfrc               = dfrc             ;//////
  pm.x                  = x                 ;
  pm.frc                = frc               ;
  pm.gbl_pk             = gbl_pk            ;
  pm.gbl_pn             = gbl_pn            ;
  pm.gbl_gamc           = gbl_gamc          ;
  pm.gbl_gams           = gbl_gams          ;
  pm.ti_signs           = ti_signs          ;
  pm.ti_ene             = ti_ene            ;
  pm.ti_ene_delta       = ti_ene_delta      ;
  pm.ti_weights         = ti_weights        ;
  pm.bar_cont           = bar_cont          ;
  pm.bar_lambda         = bar_lambda        ;
  
  pm.gbl_ipn            = gbl_ipn           ;
  pm.gbl_atm_owner_map  = gbl_atm_owner_map ;
  pm.atm_in             = atm_in            ;//packing
  pm.gmul               = gmul              ;
  pm.tm24               = tm24              ;
  pm.tm06               = tm06              ;
  pm.tenm3              = tenm3             ;
  pm.zero               = zero              ;
  pm.one                = one               ;

  //int dd;
  //for(dd = 1; dd <= dihed_cnt; dd++)
  //{
  //  int ai = dihed[dd].atm_i;
  //  int aj = dihed[dd].atm_j;
  //  int ak = abs(dihed[dd].atm_k);
  //  int al = abs(dihed[dd].atm_l);

  //  dihed_crd[dd].crd_i[0] = x[ai * 3 + 0];
  //  dihed_crd[dd].crd_i[1] = x[ai * 3 + 1];
  //  dihed_crd[dd].crd_i[2] = x[ai * 3 + 2];
  //  dihed_crd[dd].atmi_map = gbl_atm_owner_map[ai];
  //  
  //  dihed_crd[dd].crd_j[0] = x[aj * 3 + 0];
  //  dihed_crd[dd].crd_j[1] = x[aj * 3 + 1];
  //  dihed_crd[dd].crd_j[2] = x[aj * 3 + 2];
  //  dihed_crd[dd].atmj_map = gbl_atm_owner_map[aj];
  //
  //  dihed_crd[dd].crd_k[0] = x[ak * 3 + 0];
  //  dihed_crd[dd].crd_k[1] = x[ak * 3 + 1];
  //  dihed_crd[dd].crd_k[2] = x[ak * 3 + 2];
  //  dihed_crd[dd].atmk_map = gbl_atm_owner_map[ak];
  //  
  //  dihed_crd[dd].crd_l[0] = x[al * 3 + 0];
  //  dihed_crd[dd].crd_l[1] = x[al * 3 + 1];
  //  dihed_crd[dd].crd_l[2] = x[al * 3 + 2];
  //  dihed_crd[dd].atml_map = gbl_atm_owner_map[al];
  //}

  //time2 = MPI_Wtime();

#define SUNWAY
#ifdef SUNWAY

    
  if(athread_idle() == 0)
    athread_init();
  
  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  //conf.pcr2 = PC2_N_DMA_REQ;
  //lwpf_init(&conf);
  //time3 = MPI_Wtime();

  athread_spawn(dihedrals_pack_c_para, &pm);
  athread_join();

  time4 = MPI_Wtime();
  athread_spawn(get_dihedrals_energy_c_para, &pm);
  athread_join();
  
  time5 = MPI_Wtime();
  //if (myrank == 0) {
  //  lwpf_report_summary(stdout, &conf);
  //}
  //r++;
  ////////compute for frc////////
  for(i = 1; i <= dihed_cnt; i++)
  {
    frc[dihed[i].atm_i * 3 + 0] += pm.dfrc[i].frc_i[0];
    frc[dihed[i].atm_i * 3 + 1] += pm.dfrc[i].frc_i[1];
    frc[dihed[i].atm_i * 3 + 2] += pm.dfrc[i].frc_i[2];

    frc[dihed[i].atm_j * 3 + 0] += pm.dfrc[i].frc_j[0];
    frc[dihed[i].atm_j * 3 + 1] += pm.dfrc[i].frc_j[1];
    frc[dihed[i].atm_j * 3 + 2] += pm.dfrc[i].frc_j[2];

    frc[abs(dihed[i].atm_k) * 3 + 0] += pm.dfrc[i].frc_k[0];
    frc[abs(dihed[i].atm_k) * 3 + 1] += pm.dfrc[i].frc_k[1];
    frc[abs(dihed[i].atm_k) * 3 + 2] += pm.dfrc[i].frc_k[2];

    frc[abs(dihed[i].atm_l) * 3 + 0] += pm.dfrc[i].frc_l[0];
    frc[abs(dihed[i].atm_l) * 3 + 1] += pm.dfrc[i].frc_l[1];
    frc[abs(dihed[i].atm_l) * 3 + 2] += pm.dfrc[i].frc_l[2];
  }

  time6 = MPI_Wtime();
  

  //time out
  //if(myrank == 0)
  //  {
  //  printf("total: %.10f, crd_pack = %.10f, calculate_time = %.10f,write_frc = %.10f\n",time6-time1,time4-time1, time5-time4, time6-time5);
  //  }

  free(atm_in);
  free(dfrc);
  free(dihed_crd);


#else

  // Local variabl:
  GBFloat ap;
  GBFloat cosnp;
  GBFloat cphi;
  GBFloat ct, ct0;
  GBFloat dc1, dc2, dc3, dc4, dc5, dc6;
  GBFloat df;
  GBFloat dr1, dr2, dr3, dr4, dr5, dr6;
  GBFloat drx, dry, drz;
  GBFloat dums;
  GBFloat dx, dy, dz;
  GBFloat epl;
  double epw;
  GBFloat f1, f2;
  GBFloat fxi, fyi, fzi;
  GBFloat fxj, fyj, fzj;
  GBFloat fxk, fyk, fzk;
  GBFloat fxl, fyl, fzl;
  GBFloat g;
  GBFloat gx, gy, gz;
  int ic, ic0;
  int inc;
  int k, kt, l, lt;
  int jn, j_out,jmax;
  GBFloat s;
  GBFloat sinnp;
  GBFloat sphi;
  GBFloat xa, ya, za;
  GBFloat xij, yij, zij;
  GBFloat xkj, ykj, zkj;
  GBFloat xkl, ykl, zkl;
  GBFloat z1, z2;
  GBFloat z11, z12, z22;
  
  epl = zero;

  //for(j_out = 1; j_out <= dihed_cnt; j_out += 64)
  for(j_out = 1; j_out <= dihed_cnt; j_out ++)
  {
    //jmax = j_out + 63;

    //if(jmax > dihed_cnt) jmax = dihed_cnt;

    //for(jn = j_out; jn <= jmax; jn++)
    {
      jn  = j_out;
      i   = dihed[jn].atm_i;
      j   = dihed[jn].atm_j;
      kt  = dihed[jn].atm_k;
      lt  = dihed[jn].atm_l;
      k   = abs(kt);
      l   = abs(lt);

      // Calculation of ij, kj, kl vectors:
      //xij = ToGBFloat(x[i * 3 + 0] - x[j* 3 + 0]);
      //yij = ToGBFloat(x[i * 3 + 1] - x[j* 3 + 1]);
      //zij = ToGBFloat(x[i * 3 + 2] - x[j* 3 + 2]);
      //xkj = ToGBFloat(x[k * 3 + 0] - x[j* 3 + 0]);
      //ykj = ToGBFloat(x[k * 3 + 1] - x[j* 3 + 1]);
      //zkj = ToGBFloat(x[k * 3 + 2] - x[j* 3 + 2]);
      //xkl = ToGBFloat(x[k * 3 + 0] - x[l* 3 + 0]);
      //ykl = ToGBFloat(x[k * 3 + 1] - x[l* 3 + 1]);
      //zkl = ToGBFloat(x[k * 3 + 2] - x[l* 3 + 2]);
      double xi0, xi1, xi2; 
      double xj0, xj1, xj2; 
      double xk0, xk1, xk2; 
      double xl0, xl1, xl2; 
      xi0 = dihed_crd[jn].crd_i[0];
      xi1 = dihed_crd[jn].crd_i[1];
      xi2 = dihed_crd[jn].crd_i[2];
      xj0 = dihed_crd[jn].crd_j[0];
      xj1 = dihed_crd[jn].crd_j[1];
      xj2 = dihed_crd[jn].crd_j[2];
      xk0 = dihed_crd[jn].crd_k[0];
      xk1 = dihed_crd[jn].crd_k[1];
      xk2 = dihed_crd[jn].crd_k[2];
      xl0 = dihed_crd[jn].crd_l[0];
      xl1 = dihed_crd[jn].crd_l[1];
      xl2 = dihed_crd[jn].crd_l[2];

      xij = ToGBFloat(xi0 - xj0);
      yij = ToGBFloat(xi1 - xj1);
      zij = ToGBFloat(xi2 - xj2);
      xkj = ToGBFloat(xk0 - xj0);
      ykj = ToGBFloat(xk1 - xj1);
      zkj = ToGBFloat(xk2 - xj2);
      xkl = ToGBFloat(xk0 - xl0);
      ykl = ToGBFloat(xk1 - xl1);
      zkl = ToGBFloat(xk2 - xl2);

      //if(xi0 != atm_in[i].x[0] || xi1 != atm_in[i].x[1] || xi2 != atm_in[i].x[2] ||
      //   xj0 != atm_in[j].x[0] || xj1 != atm_in[j].x[1] || xj2 != atm_in[j].x[2] ||
      //   xk0 != atm_in[k].x[0] || xk1 != atm_in[k].x[1] || xk2 != atm_in[k].x[2] ||
      //   xl0 != atm_in[l].x[0] || xl1 != atm_in[l].x[1] || xl2 != atm_in[l].x[2] 
      //)
      //{
      //  printf("%lf, %lf, %lf; %lf, %lf, %lf\n", xi0, xi1, xi2, atm_in[i].x[0], atm_in[i].x[1], atm_in[i].x[2]);
      //  printf("%lf, %lf, %lf; %lf, %lf, %lf\n", xj0, xj1, xj2, atm_in[j].x[0], atm_in[j].x[1], atm_in[j].x[2]);
      //  printf("%lf, %lf, %lf; %lf, %lf, %lf\n", xk0, xk1, xk2, atm_in[k].x[0], atm_in[k].x[1], atm_in[k].x[2]);
      //  printf("%lf, %lf, %lf; %lf, %lf, %lf\n\n\n\n", xl0, xl1, xl2, atm_in[l].x[0], atm_in[l].x[1], atm_in[l].x[2]);

      //}
      //xij = ToGBFloat(atm_in[i].x[0] - atm_in[j].x[0]);
      //yij = ToGBFloat(atm_in[i].x[1] - atm_in[j].x[1]);
      //zij = ToGBFloat(atm_in[i].x[2] - atm_in[j].x[2]);
      //xkj = ToGBFloat(atm_in[k].x[0] - atm_in[j].x[0]);
      //ykj = ToGBFloat(atm_in[k].x[1] - atm_in[j].x[1]);
      //zkj = ToGBFloat(atm_in[k].x[2] - atm_in[j].x[2]);
      //xkl = ToGBFloat(atm_in[k].x[0] - atm_in[l].x[0]);
      //ykl = ToGBFloat(atm_in[k].x[1] - atm_in[l].x[1]);
      //zkl = ToGBFloat(atm_in[k].x[2] - atm_in[l].x[2]);

      // Get the normal vector:
      dx = yij * zkj - zij * ykj;
      dy = zij * xkj - xij * zkj;
      dz = xij * ykj - yij * xkj;
      gx = zkj * ykl - ykj * zkl;
      gy = xkj * zkl - zkj * xkl;
      gz = ykj * xkl - xkj * ykl;

      fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24);
      fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24);
      ct  = dx * gx + dy * gy + dz * gz;

      // Branch if linear dihedral:
      if (tenm3 <= fxi) z1 = one / fxi;
      else z1 = zero;

      if (tenm3 <= fyi) z2 = one / fyi;
      else z2 = zero;
      
      z12 = z1 * z2;

      if (z12 != zero) fzi = one;
      else fzi = zero;

      s = xkj * (dz * gy - dy * gz) 
        + ykj * (dx * gz - dz * gx)
        + zkj * (dy * gx - dx * gy);

      ap = PI - sign(acos(max(-one,min(one,ct*z12))), s);

      cphi = cos(ap);
      sphi = sin(ap);

      // Calculate the energy and the derivatives with respect to cosphi:
      ic    = dihed[jn].parm_idx;
      inc   = gbl_ipn[ic];
      ct0   = gbl_pn[ic] * ap;
      cosnp = cos(ct0);
      sinnp = sin(ct0);
      epw   = (gbl_pk[ic] + cosnp*gbl_gamc[ic] + sinnp*gbl_gams[ic]) * fzi;

      dums = sphi + sign(tm24, sphi);
      if (tm06 > fabs(dums))
      {
        df = fzi * gbl_gamc[ic] * (gbl_pn[ic] - 
              gmul[inc-1] + gmul[inc-1] * cphi);
      }
      else
      {
        df = fzi * gbl_pn[ic] * (gbl_gamc[ic] * sinnp - 
             gbl_gams[ic] * cosnp)/dums;
      }
      
      //if (ti_mode != 0) 
      //{
      //  ti_region = 0;
      //  if (ti_lst[i * 3 + 0] == 1 || ti_lst[j * 3 + 0] == 1 ||
      //      ti_lst[k * 3 + 0] == 1 || ti_lst[l * 3 + 0] == 1)
      //  {
      //    ti_region = 1;
      //  }
      //  else if(ti_lst[i * 3 + 1] == 1 || ti_lst[j * 3 + 1] == 1 ||
      //          ti_lst[k * 3 + 1] == 1 || ti_lst[l * 3 + 1] == 1)
      //  {
      //    ti_region = 2;
      //  }

      //  if (ti_region > 0)
      //  {
      //    if (ti_mode != 1)
      //    {
      //      ti_lscale = (ti_sc_lst[i] < 2 && ti_sc_lst[j] < 2 &&
      //                   ti_sc_lst[k] < 2 && ti_sc_lst[l] < 2);
      //    }
      //    //if (gbl_atm_owner_map[i] == mytaskid)
      //    if (atm_in[i].atm_map == mytaskid)
      //    {
      //       if (ti_mode == 1 || ti_lscale)
      //       {//!linear scaling
      //         ti_update_ene(epw, 
      //                       si_dihedral_ene, 
      //                       ti_region,
      //                       ti_ene,
      //                       ti_ene_delta,
      //                       ti_signs,
      //                       ti_weights,
      //                       ifmbar_lcl,
      //                       do_mbar,
      //                       bar_states,
      //                       bar_cont,
      //                       bar_lambda,
      //                       si_dvdl);

      //         epw = epw * ti_weights[ti_region];
      //       }
      //       else
      //       {
      //         ti_ene[si_dihedral_ene * 2 + ti_region - 1] = 
      //            ti_ene[si_dihedral_ene * 2 + ti_region - 1] + epw;
      //         epw = 0.0;
      //       }
      //    }
      //    if(ti_mode == 1 || ti_lscale) 
      //      df = df * ti_weights[ti_region];
      //  }//if
      //}//if-ti_mode

      // Now do torsional first derivatives:
      // Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:
      z11 = z1 * z1;
      z12 = z1 * z2;
      z22 = z2 * z2;
      dc1 = -gx * z12 - cphi * dx * z11;
      dc2 = -gy * z12 - cphi * dy * z11;
      dc3 = -gz * z12 - cphi * dz * z11;
      dc4 =  dx * z12 + cphi * gx * z22;
      dc5 =  dy * z12 + cphi * gy * z22;
      dc6 =  dz * z12 + cphi * gz * z22;

      // Update the first derivative array:
      dr1 = df * ( dc3 * ykj - dc2 * zkj);
      dr2 = df * ( dc1 * zkj - dc3 * xkj);
      dr3 = df * ( dc2 * xkj - dc1 * ykj);
      dr4 = df * ( dc6 * ykj - dc5 * zkj);
      dr5 = df * ( dc4 * zkj - dc6 * xkj);
      dr6 = df * ( dc5 * xkj - dc4 * ykj);
      drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl);
      dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl);
      drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl);
      fxi = -dr1;
      fyi = -dr2;
      fzi = -dr3;
      fxj = -drx + dr1;
      fyj = -dry + dr2;
      fzj = -drz + dr3;
      fxk = +drx + dr4;
      fyk = +dry + dr5;
      fzk = +drz + dr6;
      fxl = -dr4;
      fyl = -dr5;
      fzl = -dr6;

      //if (gbl_atm_owner_map[i] == mytaskid)
      if (atm_in[i].atm_map == mytaskid)
      {
         epl = epl + epw;
         frc[i * 3 + 0] = frc[i * 3 + 0] + fxi;
         frc[i * 3 + 1] = frc[i * 3 + 1] + fyi;
         frc[i * 3 + 2] = frc[i * 3 + 2] + fzi;
      }

      //if (gbl_atm_owner_map[j] == mytaskid)
      if (atm_in[j].atm_map == mytaskid)
      {
         frc[j * 3 + 0] = frc[j * 3 + 0] + fxj;
         frc[j * 3 + 1] = frc[j * 3 + 1] + fyj;
         frc[j * 3 + 2] = frc[j * 3 + 2] + fzj;
      }  

      //if(gbl_atm_owner_map[k] == mytaskid)
      if(atm_in[k].atm_map == mytaskid)
      {
         frc[k *3 + 0] = frc[k * 3 + 0] + fxk;
         frc[k *3 + 1] = frc[k * 3 + 1] + fyk;
         frc[k *3 + 2] = frc[k * 3 + 2] + fzk;
      }

      //if(gbl_atm_owner_map[l] == mytaskid)
      if(atm_in[l].atm_map == mytaskid)
      {
         frc[l * 3 + 0] = frc[l * 3 + 0] + fxl;
         frc[l * 3 + 1] = frc[l * 3 + 1] + fyl;
         frc[l * 3 + 2] = frc[l * 3 + 2] + fzl;
      }
    }//do
  }//do

  *ep = epl;
  
  free(atm_in);
  free(dfrc);
  //free(mpe_frc);
  free(dihed_crd);
  
#endif

  return;
  
}

#endif



#ifdef CPE
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <dma_macros.h>

#define LWPF_UNIT U(BOND_PACK)
#define LWPF_KERNELS K(ALL) K(GET) K(FOR) K(CACHE) K(PUT) 
#include <lwpf2/lwpf2.h>
#include <poly_math.h>

#define DMA_FAST

#define ISTEP 32
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

void memset_cpe(memset_param *pa)
{
  //printf("memset_cpe\n");
  dma_init();
  int i, j, k;
  int ist, ied, isz;
  doublev4 zerov4;
  double *tmp = &zerov4;
  tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0;
  double szero[ZNUM];
  memset_param la;
  pe_get(pa, &la, sizeof(memset_param));
  dma_syn();

  double *array = la.array;
  int n = la.cnt;
  for(i = 0; i < ZNUM; i+=4)
  {
    simd_store(zerov4, &szero[i]);
  }
  
  for(ist = _MYID*ZNUM; ist < n; ist += ZNUM*64)
  {
    ied = ist + ZNUM;
    if(ied > n)
      ied = n;
    isz = ied - ist;

    pe_put(array+ist, szero, sizeof(double)*isz);
    dma_syn();

  }
}

void memset_slave(memset_param *pa)
{
  int i, j, k;
  int ist, ied, isz;
  doublev4 zerov4;
  zerov4 = 0;
  double szero[ZNUM];

  memset_param la;
  dma_init();
  pe_get(pa, &la, sizeof(memset_param));
  dma_syn();

  int n = la.cnt;
  double *array = la.array;
  for(i = 0; i < ZNUM; i++)
  {
    szero[i] = 0;
  }
  
  for(ist = _MYID*ZNUM; ist < n; ist += ZNUM*64)
  {
    ied = ist + ZNUM;
    if(ied > n)
      ied = n;
    isz = ied - ist;

    pe_put(array+ist, szero, sizeof(double)*isz);
    dma_syn();

  }
}

static void* miss_get( int *c_t, void *cache_data, void *data, int data_size, int id)
{
  int t_s = T_CPE_SIZE;
  int p_s = P_CPE_SIZE;
  int tu_size = TU_SIZE;
  int p_ps = id / t_s;
  int cache_p_ps = p_ps % p_s;
  int cache_t_ps = id - p_ps * t_s;
  char *C_D = (char*)(cache_data);
  char *D = (char*)(data);
  int s = data_size / sizeof(char);
  int ps = cache_p_ps * (tu_size + 1);
  int get_cache_ps = c_t[ps + tu_size];
  int i,j,k;
  for(i = 0;i < tu_size;i++)
  {
    if(c_t[ps + i] == p_ps)
    {
        return (void*)(&C_D[(((cache_p_ps * tu_size + i) * t_s) + cache_t_ps) * s]);
    }
  }

  {
    volatile unsigned long get_reply, put_reply;
    get_reply = 0;
    athread_get(PE_MODE, D + (p_ps * t_s) * s, 
                          C_D + ((cache_p_ps * tu_size + get_cache_ps ) * t_s) * s, 
                          data_size * t_s, (void*)&get_reply,0,0,0);
    while(get_reply != 1);
    c_t[ps + get_cache_ps] = p_ps;
  }
  c_t[ps + tu_size] = (get_cache_ps + 1) % tu_size;
  return (void*)(&C_D[(((cache_p_ps * tu_size + get_cache_ps) * t_s) + cache_t_ps) * s]);
}

void pack_atom_slave(dihed_param_c *pm)
{
  dma_init();
  dihed_param_c lp;
  pe_get(pm,&lp,sizeof(dihed_param_c));
  dma_syn();

  int          atm_cnt           = lp.atm_cnt;
  //dihed_atm_c *atm_in            = lp.atm_in;
  //double      *x                 = lp.x;
  //int         *gbl_atm_owner_map = lp.gbl_atm_owner_map;
 
  dihed_atm_c    atm_in[ISTEP];
  double         x[ISTEP * 3];
  int            gbl_atm_owner_map[ISTEP];

  int i,j,k,ii;
  int ist,ied,isz;
  for(ist = 1 + _MYID * ISTEP; ist < atm_cnt+1; ist += ISTEP * 64) 
  {
    ied = ist + ISTEP;
    if(ied > atm_cnt + 1)
      ied = atm_cnt + 1;
    isz = ied - ist;
    
    pe_get(lp.x + ist * 3, x , sizeof(double)*isz*3);
    pe_get(lp.gbl_atm_owner_map + ist , gbl_atm_owner_map , sizeof(int)*isz);
    dma_syn();
    
    for(i = ist; i < ied; i++ )
    {
      ii = i - ist;  
      //ii = i;
      atm_in[ii].x[0]  = x[ii * 3 + 0];
      atm_in[ii].x[1]  = x[ii * 3 + 1];
      atm_in[ii].x[2]  = x[ii * 3 + 2];
      atm_in[ii].atm_map = gbl_atm_owner_map[ii];
    }
    pe_put(lp.atm_in + ist,atm_in,sizeof(dihed_atm_c)*isz);
    dma_syn();
  }
}

#define ADDSTEP 1024
//void add_frc_slave(dihed_param_c *pm)
//{
//  dihed_param_c lp;
//  dma_init();
//  pe_get(pm,&lp,sizeof(dihed_param_c));
//  dma_syn();
//  
//  int     atm_cnt    = lp.atm_cnt;
//  
//  int i,j,k,ii;
//  int ist,ied,isz;
//  double frc[ADDSTEP * 3];
//  double mpe_frc[ADDSTEP * 3];
//
//  for(ist = 1 + ADDSTEP * _MYID; ist <= atm_cnt; ist += ADDSTEP * 64)
//  {
//    ied = ist + ADDSTEP;
//    if(ied > atm_cnt + 1)
//      ied = atm_cnt + 1;
//    isz = ied - ist;
//    
//    pe_get(lp.frc + 3 * ist, frc,  isz * sizeof(double) * 3);
//    dma_syn();
//    
//    for(j = 0; j < 64; j++)
//    {
//      pe_get(lp.mpe_frc + 3 * (j * (atm_cnt + 512) + ist), 
//                 mpe_frc, isz * sizeof(double) * 3);
//      dma_syn();
//      for(i = 0; i < isz; i++)
//      {
//        frc[i * 3 + 0] += mpe_frc[i * 3 + 0]; 
//        frc[i * 3 + 1] += mpe_frc[i * 3 + 1];
//        frc[i * 3 + 2] += mpe_frc[i * 3 + 2];
//      }
//    }
//    pe_put( lp.frc + 3 * ist,frc, isz * sizeof(double) * 3);
//    dma_syn();
//  }
//}

void reg_reduce_inplace_doublev4(doublev4 *arr, int len)
{
  int i, j, n;
  doublev4 tmp;
  for (i = 1; i < 8; i += i)
  {
    if ((_ROW & i) == i)
    {
      for (j = 0; j < len; j ++)
        asm("putc %0, %1": : "r"(arr[j]), "r"(_ROW ^ i));
    }
    if ((_ROW & i) == 0)
    {
      for (j = 0; j < len; j ++)
      {
        asm("getc %0\n\t"
            "vaddd %0, %1, %1\n\t"
            : "=r"(tmp), "+r"(arr[j]));
      }
    }
    athread_syn(COL_SCOPE, 0xff);
  }

  athread_syn(ARRAY_SCOPE, 0xffff);
  if (_ROW == 0)
  {
    for (i = 1; i < 8; i += i)
    {
      if ((_COL & i) == i)
      {
        for (j = 0; j < len; j ++)
          asm("putr %0, %1": : "r"(arr[j]), "r"(_COL ^ i));
      }
      if ((_COL & i) == 0)
      {
        for (j = 0; j < len; j ++)
        {
          asm("getr %0\n\t" 
              "vaddd %0, %1, %1\n\t"
              : "=r"(tmp), "+r"(arr[j]));
        }
      }
    }
    athread_syn(ROW_SCOPE, 0xff);
  }

}

double fabs_sw(double a)
{
  if(a < 0)
    return -a;
  return a;
}

GBFloat sign_sw(GBFloat x, GBFloat y)
{
  if(y > 0)
    return fabs_sw(x);
  else if(y < 0)
    return -fabs_sw(x);
  else
    return 0;
}

void write_cache_sw(int i,
                    int *tag_cache, 
                    double frc_cache[][CACHE_LINESIZE][3],
                    double *frc,
                    double fx, double fy, double fz, 
                    int *tot_cmp, int *miss_cache)
{
  dma_init();
  double *pfrc;
  if(tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] == -1) 
  {
     pe_get(frc + (i & ~CACHE_MMASK)*3, 
            frc_cache[(i >> CACHE_SBIT) & CACHE_LMASK],
            sizeof(double)*3*CACHE_LINESIZE);
     dma_syn();
     tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] 
     = i & ~CACHE_MMASK;//i >> CACHE_SBIT;
     
     *miss_cache = (*miss_cache) + 1;
  }
  else if(tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] != 
        (i & ~CACHE_MMASK)) //j >> JCACHE_SBIT) 
  {
    pe_put(frc + tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK]*3, 
            frc_cache[(i >> CACHE_SBIT) & CACHE_LMASK],
            sizeof(double)*3*CACHE_LINESIZE);
    pe_get(frc + (i & ~CACHE_MMASK)*3, 
            frc_cache[(i >> CACHE_SBIT) & CACHE_LMASK],
            sizeof(double)*3*CACHE_LINESIZE);
    dma_syn();
    tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] 
     = i & ~CACHE_MMASK;

     *miss_cache = (*miss_cache) + 1;
  }
  //else if(tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] 
  //       == (i & ~CACHE_SBIT))

  *tot_cmp = (*tot_cmp) + 1;

   pfrc = frc_cache[(i >> CACHE_SBIT) & CACHE_LMASK][i & CACHE_MMASK];
   pfrc[0] += fx;
   pfrc[1] += fy;
   pfrc[2] += fz;

}

#define DSTEP 32
#define PACK_HBIT 10
#define PACK_SBIT 3
#define PACK_LINESIZE (1 << PACK_SBIT)
#define PACK_LINECNT (1 << (PACK_HBIT - PACK_SBIT))
#define PACK_MMASK (PACK_LINESIZE - 1)
#define PACK_LMASK (PACK_LINECNT - 1)
void read_cache_sw( int ai,
                    double x_cache[][PACK_LINESIZE][3], 
                    double *x,
                    int map_cache[][PACK_LINESIZE], 
                    int *gbl_atm_owner_map,
                    int *tag_cache,
                    double *xci,
                    int *mci)
{
  dma_init();
  //cache
  if(tag_cache[(ai >> PACK_SBIT) & PACK_LMASK] != ai >> PACK_SBIT)
  {
    pe_get(x + (ai & ~PACK_MMASK) * 3, x_cache[(ai >> PACK_SBIT) & PACK_LMASK], 
          sizeof(double) * 3 * PACK_LINESIZE);
    pe_get(gbl_atm_owner_map + (ai & ~PACK_MMASK), 
          map_cache[(ai >> PACK_SBIT) & PACK_LMASK], 
          sizeof(int) * PACK_LINESIZE);
    dma_syn();
    tag_cache[(ai >> PACK_SBIT) & PACK_LMASK] = ai >> PACK_SBIT;
  }
  double *tmpx =   x_cache[(ai >> PACK_SBIT) & PACK_LMASK][ai & PACK_MMASK];
  xci[0] = tmpx[0];
  xci[1] = tmpx[1];
  xci[2] = tmpx[2];

  *mci = map_cache[(ai >> PACK_SBIT) & PACK_LMASK][ai & PACK_MMASK];
}
void dihedrals_pack_c_para(dihed_param_c *pm)
{
  //lwpf_enter(BOND_PACK)
  //lwpf_start(ALL);

  dma_init();
  dihed_param_c lp;
  pe_get(pm, &lp, sizeof(dihed_param_c));
  dma_syn();

  dihed_rec_c *dihed      = lp.dihed;
  dihed_crd_c *dihed_crd  = lp.dihed_crd;
  int dihed_cnt           = lp.dihed_cnt;
  double *x               = lp.x;
  double (*x3)[3]         = lp.x;
  int *gbl_atm_owner_map  = lp.gbl_atm_owner_map;
  
  dihed_rec_c l_dihed[DSTEP];
  dihed_crd_c l_dihed_crd[DSTEP];

  double x_cache[PACK_LINECNT][PACK_LINESIZE][3];
  int  map_cache[PACK_LINECNT][PACK_LINESIZE];
  int  tag_cache[PACK_LINECNT];
  double *xci, *xcj, *xck, *xcl;
  double tmp_xci[3], tmp_xcj[3], tmp_xck[3], tmp_xcl[3];
  int mci, mcj, mck, mcl;
  int i;
  
  for (i = 0; i < PACK_LINECNT; i ++)
    tag_cache[i] = -1;

  int dd;
  int dst, ded, dsz, doff;
  for(dst = 1 + _MYID * DSTEP; dst < dihed_cnt + 1; dst += DSTEP * 64)
  {
    ded = dst + DSTEP;
    if(ded > dihed_cnt + 1)
      ded = dihed_cnt + 1;
    dsz = ded - dst;
    
    //lwpf_start(GET);
    pe_get(dihed + dst, l_dihed, sizeof(dihed_rec_c) * dsz);
    dma_syn();
    //lwpf_stop(GET);

    //lwpf_start(FOR);
    for(dd = dst; dd < ded; dd++)
    {
      doff = dd - dst;

      int ai = l_dihed[doff].atm_i;
      int aj = l_dihed[doff].atm_j;
      int ak = abs(l_dihed[doff].atm_k);
      int al = abs(l_dihed[doff].atm_l);

      //lwpf_start(CACHE);
      ////cache
      read_cache_sw(ai, x_cache, x, 
                    map_cache, gbl_atm_owner_map,
                    tag_cache, tmp_xci, &mci);
      read_cache_sw(aj, x_cache, x, 
                    map_cache, gbl_atm_owner_map,
                    tag_cache, tmp_xcj, &mcj);
      read_cache_sw(ak, x_cache, x, 
                    map_cache, gbl_atm_owner_map,
                    tag_cache, tmp_xck, &mck);
      read_cache_sw(al, x_cache, x, 
                    map_cache, gbl_atm_owner_map,
                    tag_cache, tmp_xcl, &mcl);

      //lwpf_stop(CACHE);
      l_dihed_crd[doff].crd_i[0] = tmp_xci[0];
      l_dihed_crd[doff].crd_i[1] = tmp_xci[1];
      l_dihed_crd[doff].crd_i[2] = tmp_xci[2];
      l_dihed_crd[doff].atmi_map = mci;
      l_dihed_crd[doff].crd_j[0] = tmp_xcj[0];
      l_dihed_crd[doff].crd_j[1] = tmp_xcj[1];
      l_dihed_crd[doff].crd_j[2] = tmp_xcj[2];
      l_dihed_crd[doff].atmj_map = mcj;
      l_dihed_crd[doff].crd_k[0] = tmp_xck[0];
      l_dihed_crd[doff].crd_k[1] = tmp_xck[1];
      l_dihed_crd[doff].crd_k[2] = tmp_xck[2];
      l_dihed_crd[doff].atmk_map = mck;
      l_dihed_crd[doff].crd_l[0] = tmp_xcl[0];
      l_dihed_crd[doff].crd_l[1] = tmp_xcl[1];
      l_dihed_crd[doff].crd_l[2] = tmp_xcl[2];
      l_dihed_crd[doff].atml_map = mcl;
    }//for-dd
    //lwpf_stop(FOR);

    //lwpf_start(PUT);
    pe_put(dihed_crd + dst, l_dihed_crd, sizeof(dihed_crd_c) * dsz);
    dma_syn();
    //lwpf_stop(PUT);
  }//for-dst
  
  //lwpf_stop(ALL);
  //lwpf_exit(BOND_PACK);

}

void get_dihedrals_energy_c_para(dihed_param_c *pm)
{
  //lwpf_enter(BOND)
  //lwpf_start(ALL);
  
  dma_init();
  dihed_param_c lp;
  pe_get(pm, &lp, sizeof(dihed_param_c));
  dma_syn();

  int dihed_cnt           = lp.dihed_cnt;
  int atm_cnt             = lp.atm_cnt;
  double *ep              = lp.ep;/////!!!!!*
  int mytaskid            = lp.mytaskid;

  dihed_rec_c *dihed      = lp.dihed;
  dihed_crd_c *dihed_crd  = lp.dihed_crd;
  dihed_frc_c *dfrc       = lp.dfrc;
  double *x               = lp.x;
  double *gbl_pk          = lp.gbl_pk;
  double *gbl_pn          = lp.gbl_pn;
  double *gbl_gamc        = lp.gbl_gamc;
  double *gbl_gams        = lp.gbl_gams;
  double *ti_signs        = lp.ti_signs;
  double *ti_ene          = lp.ti_ene;
  double *ti_ene_delta    = lp.ti_ene_delta;
  double *ti_weights      = lp.ti_weights;
  double *bar_cont        = lp.bar_cont;
  double *bar_lambda      = lp.bar_lambda;
  
  int *gbl_ipn            = lp.gbl_ipn;
  int *gbl_atm_owner_map  = lp.gbl_atm_owner_map;
  //dihed_atm_c *atm_in     = lp.atm_in;//packing
  GBFloat *gmul           = lp.gmul;
  GBFloat tm24            = lp.tm24 ;
  GBFloat tm06            = lp.tm06 ;
  GBFloat tenm3           = lp.tenm3;
  GBFloat zero            = lp.zero ;
  GBFloat one             = lp.one  ;

  // Local variable:
  GBFloat ap;
  GBFloat cosnp;
  GBFloat cphi;
  GBFloat ct, ct0;
  GBFloat dc1, dc2, dc3, dc4, dc5, dc6;
  GBFloat df;
  GBFloat dr1, dr2, dr3, dr4, dr5, dr6;
  GBFloat drx, dry, drz;
  GBFloat dums;
  GBFloat dx, dy, dz;
  GBFloat epl;
  //double epl;
  double epw;
  GBFloat f1, f2;
  GBFloat fxi, fyi, fzi;
  GBFloat fxj, fyj, fzj;
  GBFloat fxk, fyk, fzk;
  GBFloat fxl, fyl, fzl;
  GBFloat g;
  GBFloat gx, gy, gz;
  int ic, ic0;
  int inc;
  int i, j, k, kt, l, lt;
  int jn, j_out,jmax;
  GBFloat s;
  GBFloat sinnp;
  GBFloat sphi;
  GBFloat xa, ya, za;
  GBFloat xij, yij, zij;
  GBFloat xkj, ykj, zkj;
  GBFloat xkl, ykl, zkl;
  GBFloat z1, z2;
  GBFloat z11, z12, z22;
  
  GBFloat l_gmul[10];
  int l_gbl_ipn[C_GBL];            
  double l_gbl_pn[C_GBL], l_gbl_pk[C_GBL];
  double l_gbl_gamc[C_GBL], l_gbl_gams[C_GBL];
  dihed_rec_c l_dihed[ISTEP];
  dihed_crd_c l_dihed_crd[ISTEP];
  dihed_frc_c l_dfrc[ISTEP];
  doublev4 eplv4[2];
  double *vep = (double *)(void *)eplv4;
  
  //nptra = 300
  pe_get(gmul,      l_gmul,      sizeof(GBFloat)*10);
  pe_get(gbl_ipn,   l_gbl_ipn,   sizeof(int)*C_NPTRA);
  pe_get(gbl_pn,    l_gbl_pn,    sizeof(double)*C_NPTRA);
  pe_get(gbl_pk,    l_gbl_pk,    sizeof(double)*C_NPTRA);
  pe_get(gbl_gamc,  l_gbl_gamc,  sizeof(double)*C_NPTRA);
  pe_get(gbl_gams,  l_gbl_gams,  sizeof(double)*C_NPTRA);
  dma_syn();

  vep[0] = vep[1] = vep[2] = vep[3] = 0;
  vep[4] = vep[5] = vep[6] = vep[7] = 0;
  epl = zero;
  vep[0] = zero;

  int ist, ied, isz;
  int joff;

  for(ist = 1 + ISTEP * _MYID; ist < dihed_cnt + 1; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > dihed_cnt + 1)
      ied = dihed_cnt + 1;
    isz = ied - ist;
    
    pe_get(dihed+ist,     l_dihed,     sizeof(dihed_rec_c)*isz);
    pe_get(dihed_crd+ist, l_dihed_crd, sizeof(dihed_crd_c)*isz);
    dma_syn();

    for(j_out = ist; j_out < ied; j_out++)
    {
      jn  = j_out;
      joff = j_out- ist;
      i   = l_dihed[joff].atm_i;
      j   = l_dihed[joff].atm_j;
      kt  = l_dihed[joff].atm_k;
      lt  = l_dihed[joff].atm_l;
      k   = abs(kt);
      l   = abs(lt);
      
      //lwpf_start(CACHE_READ);
     
      //double xi0, xi1, xi2; 
      //double xj0, xj1, xj2; 
      //double xk0, xk1, xk2; 
      //double xl0, xl1, xl2; 
      //xi0 = l_dihed_crd[joff].crd_i[0];
      //xi1 = l_dihed_crd[joff].crd_i[1];
      //xi2 = l_dihed_crd[joff].crd_i[2];
      //xj0 = l_dihed_crd[joff].crd_j[0];
      //xj1 = l_dihed_crd[joff].crd_j[1];
      //xj2 = l_dihed_crd[joff].crd_j[2];
      //xk0 = l_dihed_crd[joff].crd_k[0];
      //xk1 = l_dihed_crd[joff].crd_k[1];
      //xk2 = l_dihed_crd[joff].crd_k[2];
      //xl0 = l_dihed_crd[joff].crd_l[0];
      //xl1 = l_dihed_crd[joff].crd_l[1];
      //xl2 = l_dihed_crd[joff].crd_l[2];

      //xij = ToGBFloat(xi0 - xj0);
      //yij = ToGBFloat(xi1 - xj1);
      //zij = ToGBFloat(xi2 - xj2);
      //xkj = ToGBFloat(xk0 - xj0);
      //ykj = ToGBFloat(xk1 - xj1);
      //zkj = ToGBFloat(xk2 - xj2);
      //xkl = ToGBFloat(xk0 - xl0);
      //ykl = ToGBFloat(xk1 - xl1);
      //zkl = ToGBFloat(xk2 - xl2);
      
      xij = ToGBFloat(l_dihed_crd[joff].crd_i[0] - l_dihed_crd[joff].crd_j[0]);
      yij = ToGBFloat(l_dihed_crd[joff].crd_i[1] - l_dihed_crd[joff].crd_j[1]);
      zij = ToGBFloat(l_dihed_crd[joff].crd_i[2] - l_dihed_crd[joff].crd_j[2]);
      xkj = ToGBFloat(l_dihed_crd[joff].crd_k[0] - l_dihed_crd[joff].crd_j[0]);
      ykj = ToGBFloat(l_dihed_crd[joff].crd_k[1] - l_dihed_crd[joff].crd_j[1]);
      zkj = ToGBFloat(l_dihed_crd[joff].crd_k[2] - l_dihed_crd[joff].crd_j[2]);
      xkl = ToGBFloat(l_dihed_crd[joff].crd_k[0] - l_dihed_crd[joff].crd_l[0]);
      ykl = ToGBFloat(l_dihed_crd[joff].crd_k[1] - l_dihed_crd[joff].crd_l[1]);
      zkl = ToGBFloat(l_dihed_crd[joff].crd_k[2] - l_dihed_crd[joff].crd_l[2]);


      dx = yij * zkj - zij * ykj;
      dy = zij * xkj - xij * zkj;
      dz = xij * ykj - yij * xkj;
      gx = zkj * ykl - ykj * zkl;
      gy = xkj * zkl - zkj * xkl;
      gz = ykj * xkl - xkj * ykl;
      
      double dxyz = dx * dx + dy * dy + dz * dz + tm24;
      double dtmp;
      inv_sqrt(dxyz, dtmp);
      fxi = dxyz * dtmp;

      double gxyz = gx * gx + gy * gy + gz * gz + tm24;
      double gtmp;
      inv_sqrt(gxyz, gtmp);
      fyi = gxyz * gtmp;

      //fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24);
      //fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24);


      ct  = dx * gx + dy * gy + dz * gz;

      if (tenm3 <= fxi) z1 = one / fxi;
      else z1 = zero;

      if (tenm3 <= fyi) z2 = one / fyi;
      else z2 = zero;
      
      z12 = z1 * z2;

      if (z12 != zero) fzi = one;
      else fzi = zero;

      s = xkj * (dz * gy - dy * gz) 
        + ykj * (dx * gz - dz * gx)
        + zkj * (dy * gx - dx * gy);
      
      //lwpf_start(ONE);
      ap = PI - sign_sw(acos(max(-one,min(one,ct*z12))), s);

      //cphi = p_cosnpi_pid(ap);//0.5;//cos(ap);
      //sphi = p_sinnpi_pid(ap);//0.5;//cos(ap);
      cphi = p_cosd(ap);
      sphi = p_sind(ap);

      ic    = l_dihed[joff].parm_idx;
      inc   = l_gbl_ipn[ic];
     
      ct0   = l_gbl_pn[ic] * ap;
      
      //cosnp = p_cosnpi_pid(ct0);//0.5;//cos(ap);
      //sinnp = p_sinnpi_pid(ct0);//0.5;//cos(ap);
      cosnp = p_cosd(ct0);
      sinnp = p_sind(ct0);
      //cosnp = cos(ct0);
      //sinnp = sin(ct0);

      epw   = (l_gbl_pk[ic] + cosnp*l_gbl_gamc[ic] + sinnp*l_gbl_gams[ic]) * fzi;

      dums = sphi + sign_sw(tm24, sphi);
      
      if (tm06 > fabs_sw(dums))
      {
        df = fzi * l_gbl_gamc[ic] * (l_gbl_pn[ic] - 
              l_gmul[inc-1] + l_gmul[inc-1] * cphi);
      }
      else
      {
        df = fzi * l_gbl_pn[ic] * (l_gbl_gamc[ic] * sinnp - 
             l_gbl_gams[ic] * cosnp)/dums;
      }
      
      //lwpf_stop(ONE);
      
      z11 = z1 * z1;
      z12 = z1 * z2;
      z22 = z2 * z2;
      dc1 = -gx * z12 - cphi * dx * z11;
      dc2 = -gy * z12 - cphi * dy * z11;
      dc3 = -gz * z12 - cphi * dz * z11;
      dc4 =  dx * z12 + cphi * gx * z22;
      dc5 =  dy * z12 + cphi * gy * z22;
      dc6 =  dz * z12 + cphi * gz * z22;

      dr1 = df * ( dc3 * ykj - dc2 * zkj);
      dr2 = df * ( dc1 * zkj - dc3 * xkj);
      dr3 = df * ( dc2 * xkj - dc1 * ykj);
      dr4 = df * ( dc6 * ykj - dc5 * zkj);
      dr5 = df * ( dc4 * zkj - dc6 * xkj);
      dr6 = df * ( dc5 * xkj - dc4 * ykj);
      drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl);
      dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl);
      drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl);
      fxi = -dr1;
      fyi = -dr2;
      fzi = -dr3;
      fxj = -drx + dr1;
      fyj = -dry + dr2;
      fzj = -drz + dr3;
      fxk = +drx + dr4;
      fyk = +dry + dr5;
      fzk = +drz + dr6;
      fxl = -dr4;
      fyl = -dr5;
      fzl = -dr6;

//CACHE_TEST
      //lwpf_start(CACHE_WRITE);
      l_dfrc[joff].frc_i[0] = 0;
      l_dfrc[joff].frc_i[1] = 0;
      l_dfrc[joff].frc_i[2] = 0;
      l_dfrc[joff].frc_j[0] = 0;
      l_dfrc[joff].frc_j[1] = 0;
      l_dfrc[joff].frc_j[2] = 0;
      l_dfrc[joff].frc_k[0] = 0;
      l_dfrc[joff].frc_k[1] = 0;
      l_dfrc[joff].frc_k[2] = 0;
      l_dfrc[joff].frc_l[0] = 0;
      l_dfrc[joff].frc_l[1] = 0;
      l_dfrc[joff].frc_l[2] = 0;



      if (l_dihed_crd[joff].atmi_map == mytaskid)
      {
        vep[0] += epw;

        l_dfrc[joff].frc_i[0] = fxi;
        l_dfrc[joff].frc_i[1] = fyi;
        l_dfrc[joff].frc_i[2] = fzi;

      }

      if (l_dihed_crd[joff].atmj_map == mytaskid)
      {
        l_dfrc[joff].frc_j[0] = fxj;
        l_dfrc[joff].frc_j[1] = fyj;
        l_dfrc[joff].frc_j[2] = fzj;

      }  

      if (l_dihed_crd[joff].atmk_map == mytaskid)
      {
        l_dfrc[joff].frc_k[0] = fxk;
        l_dfrc[joff].frc_k[1] = fyk;
        l_dfrc[joff].frc_k[2] = fzk;

      }

      if (l_dihed_crd[joff].atml_map == mytaskid)
      {
        l_dfrc[joff].frc_l[0] = fxl;
        l_dfrc[joff].frc_l[1] = fyl;
        l_dfrc[joff].frc_l[2] = fzl;

      }
      //lwpf_stop(CACHE_WRITE);
    }//do
    pe_put(dfrc+ist, l_dfrc, sizeof(dihed_frc_c)*isz);
    dma_syn();

  }//do
  reg_reduce_inplace_doublev4(eplv4, 1);
  if(_MYID == 0)
  {
    pe_put(ep, &vep[0], sizeof(double));
    dma_syn();
  }
  //lwpf_stop(OTHER);
  //lwpf_stop(ALL);
  //lwpf_exit(BOND);
  return;
}


#endif
