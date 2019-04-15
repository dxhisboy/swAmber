#include <params.h> 
#ifdef MPE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#define D_R_M (4.6566128752457969E-10)
#define LWPF_NOCOLOR
#define LWPF_UNITS U(COMP)
#include <lwpf2/lwpf2.h>


extern SLAVE_FUN(lang_set_sunway_slave)();
void randn_(double *mu_, double *sigma_, double *result)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;       
  double mu = *mu_;
  double sigma = *sigma_;

  if (call == 1)
  {
    call = !call;
    *result = (mu + sigma * (double) X2);
    return;
  }

  do
  {
    U1 = -1 + ((double) rand () / RAND_MAX) * 2;
    U2 = -1 + ((double) rand () / RAND_MAX) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  call = !call;
  *result = (mu + sigma * (double) X1);
}

double randn(double mu, double sigma )
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;       

  if (call == 1)
  {
    call = !call;
    return (mu + sigma * (double) X2);
  }


  do
  {
    U1 = -1 + ((double) rand () * D_R_M) * 2;
    U2 = -1 + ((double) rand () * D_R_M) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  call = !call;
  return (mu + sigma * (double) X1);
}


void lang_set_( double *vel, double *frc, double *mass, double *mass_inv, 
                double *dtx_, double *c_explic_, double  *fln1_, 
                double *wfac_, double *c_implic_, double *sdfac_,
                int *gbl_atm_owner_map,
                int *atm_cnt_, int *no_ntt3_sync_, int *mytaskid_)
{
  double dtx = *dtx_;
  double c_explic = *c_explic_;
  double fln1 = *fln1_;
  double wfac = *wfac_;
  double c_implic = *c_implic_;
  double sdfac = *sdfac_;

  int atm_cnt = *atm_cnt_;
  int no_ntt3_sync = *no_ntt3_sync_;
  int mytaskid = *mytaskid_;
  
  int j;
  lang_set_param param;

  param.dtx               = dtx               ;  
  param.c_explic          = c_explic          ;
  param.fln1              = fln1              ;
  param.wfac              = wfac              ;
  param.c_implic          = c_implic          ;
  param.sdfac             = sdfac             ;
  param.atm_cnt           = atm_cnt           ;
  param.no_ntt3_sync      = no_ntt3_sync      ;
  param.mytaskid          = mytaskid          ;
  param.vel               = vel               ;
  param.frc               = frc               ;
  param.mass              = mass              ;
  param.mass_inv          = mass_inv          ;
  param.gbl_atm_owner_map = gbl_atm_owner_map ;

#define SUNWAY
#ifdef SUNWAY
  if(athread_idle() == 0)
    athread_init();
 
#ifdef L_W
  
    perf_config_t conf;
    conf.pcrc = PCRC_ALL;
    conf.pcr0 = PC0_CNT_INST;
    conf.pcr1 = PC1_CNT_DIV_SQRT;
    conf.pcr2 = PC2_CNT_GLD;
    lwpf_init(&conf);
  
#endif


  athread_spawn(lang_set_sunway_slave, &param);
  athread_join();


#ifdef L_W
  if(mytaskid == 0)
  {
    lwpf_report_summary(stdout, &conf);
  }
#endif



#else

  //int cnt = 0;
  //if(mytaskid == 0)
  //{
  //  for(j = 1; j < atm_cnt; j++)
  //  {
  //    if(j)
  //      if(gbl_atm_owner_map[j - 1] != gbl_atm_owner_map[j])
  //      {
  //        printf("the gbl_map = %d cnt = %d j = %d \n", gbl_atm_owner_map[j - 1], cnt, j);
  //        cnt = 0;
  //      }
  //    cnt++;
  //  }
  //}
  
  if (no_ntt3_sync)
  {
    for( j = 0; j < atm_cnt; j++)
    {
      if(gbl_atm_owner_map[j] == mytaskid)
      {
        double wfac = mass_inv[j] * dtx;
        double aamass = mass[j];
        double rsd = sdfac * sqrt(aamass);
        double fln1 =  randn(0.0, rsd); 
        double fln2 =  randn(0.0, rsd);
        double fln3 =  randn(0.0, rsd);
        vel[0 + 3 * j] = (vel[0 + 3 * j] * c_explic + (frc[0 + 3 * j] + fln1) * wfac) * c_implic; 
        vel[1 + 3 * j] = (vel[1 + 3 * j] * c_explic + (frc[1 + 3 * j] + fln2) * wfac) * c_implic;
        vel[2 + 3 * j] = (vel[2 + 3 * j] * c_explic + (frc[2 + 3 * j] + fln3) * wfac) * c_implic;
      }
    }
  }
  else
  {
    for( j = 1; j < atm_cnt; j++)
    {
      if(gbl_atm_owner_map[j] == mytaskid)
      {
        double wfac = mass_inv[j] * dtx;
        double aamass = mass[j];
        double rsd = sdfac * sqrt(aamass);
        double fln1 =  randn(0.0, rsd); 
        double fln2 =  randn(0.0, rsd);
        double fln3 =  randn(0.0, rsd);
        vel[0 + 3 * j] = (vel[0 + 3 * j] * c_explic + (frc[0 + 3 * j] + fln1) * wfac) * c_implic; 
        vel[1 + 3 * j] = (vel[1 + 3 * j] * c_explic + (frc[1 + 3 * j] + fln2) * wfac) * c_implic;
        vel[2 + 3 * j] = (vel[2 + 3 * j] * c_explic + (frc[2 + 3 * j] + fln3) * wfac) * c_implic;
      }
      else
      {
        double rsd = 1.0;
        randn(0.0, rsd); 
        randn(0.0, rsd);
        randn(0.0, rsd);
      }
    }
  }
#endif
}
#endif


#ifdef CPE

#include <math.h>
#include <stdlib.h>
#define D_R_M (3.0518509475997192E-5)
#define LWPF_UNIT U(COMP)
#define LWPF_KERNELS K(ATHREAD_SYN) K(A_P) K(A_G) K(ALL) K(FRC_CACHE) K(T1)  K(T2)  K(T3)  K(T4)
#include <lwpf2/lwpf2.h>

__thread_local static double X1, X2;
__thread_local static int call = 0;
__thread_local static int _T_ = 1;
__thread_local static unsigned long NEXT;

int slave_rand()
{
  NEXT = NEXT * 1103515245 +  12345;
  return (unsigned int)(NEXT / 65536)%32768;
}

double randn(double mu, double sigma )
{
  double U1, U2, W, mult;
  if (call == 1)
  {
    call = !call;
    return (mu + sigma * (double) X2);
  }

  do
  {

#ifdef L_W
  lwpf_start(T1);
#endif
    U1 = -1 + ((double) slave_rand () * D_R_M) * 2;
    U2 = -1 + ((double) slave_rand () * D_R_M) * 2;
   
#ifdef L_W
  lwpf_stop(T1);
#endif

    W = U1 * U1 + U2 * U2;
  }while (W >= 1 || W == 0);
   
  mult = sqrt ((-2 * log (W)) / W);
   
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;
  return (mu + sigma * (double) X1);
}


#define READ_SIZE 128
void lang_set_sunway_slave(lang_set_param *param)
{
  
  pe_init();
  lang_set_param lp;
  pe_get(param, &lp, sizeof(lang_set_param));
  pe_syn();
  int j;
  
  int jst, jed;
  int j_size;
 
  double  dtx               = lp.dtx               ;  
  double  c_explic          = lp.c_explic          ;
  double  fln1              = lp.fln1              ;
  double  wfac              = lp.wfac              ;
  double  c_implic          = lp.c_implic          ;
  double  sdfac             = lp.sdfac             ;
  int     atm_cnt           = lp.atm_cnt           ;
  int     no_ntt3_sync      = lp.no_ntt3_sync      ;
  int     mytaskid          = lp.mytaskid          ;

  call = 0;
  if(_T_)
    NEXT = rand();
  _T_ = 0;
#define TEST
#ifndef TEST
  double *vel               = lp.vel               ;
  double *frc               = lp.frc               ;
  double *mass              = lp.mass              ;
  double *mass_inv          = lp.mass_inv          ;
  int    *gbl_atm_owner_map = lp.gbl_atm_owner_map ;
  
  j_size = atm_cnt;
  if(_MYID)
    return;
#else

#ifdef L_W
  lwpf_start(ALL);
#endif


  double vel[3 * READ_SIZE], frc[3 * READ_SIZE];
  double mass[READ_SIZE], mass_inv[READ_SIZE];
  int gbl_atm_owner_map[READ_SIZE];
#endif
  if (no_ntt3_sync)
  {
    
#ifdef TEST
    for( jst = _MYID * READ_SIZE; jst < atm_cnt; jst += 64 * READ_SIZE)
    {
      jed = jst + READ_SIZE;
      if( jed > atm_cnt)
        jed = atm_cnt;
        j_size = jed - jst;
 
#ifdef L_W
  lwpf_start(A_G);
#endif 
      pe_get( lp.gbl_atm_owner_map + jst, gbl_atm_owner_map, sizeof(int) * j_size);
      pe_syn();
     
      int cnt = 0;
      for( j = 0; j <  j_size; j++)
        if(gbl_atm_owner_map[j] == mytaskid)
          cnt++;

#ifdef L_W
  lwpf_stop(A_G);
#endif 
      if(!cnt)
        continue;
 
#ifdef L_W
  lwpf_start(A_G);
#endif     
      pe_get( lp.vel + 3 * jst, vel, sizeof(double) * 3 * j_size);
      pe_get( lp.frc + 3 * jst, frc, sizeof(double) * 3 * j_size);
      pe_get( lp.mass + jst, mass, sizeof(double) * j_size);
      pe_get( lp.mass_inv + jst, mass_inv, sizeof(double) * j_size);
      pe_syn();

#ifdef L_W
  lwpf_stop(A_G);
  lwpf_start(T3);
#endif

#endif

      for(j = 0; j < j_size; j++ )
      {
        if(gbl_atm_owner_map[j] == mytaskid)
        {
          double wfac = mass_inv[j] * dtx;
          double aamass = mass[j];
          double rsd = sdfac * sqrt(aamass);
#ifdef L_W
  lwpf_start(T2);
#endif
          double fln1 =  randn(0.0, rsd); 
          double fln2 =  randn(0.0, rsd);
          double fln3 =  randn(0.0, rsd);
#ifdef L_W
  lwpf_stop(T2);
#endif
          vel[0 + 3 * j] = 
            (vel[0 + 3 * j] * c_explic + (frc[0 + 3 * j] + fln1) * wfac) * c_implic; 
          vel[1 + 3 * j] = 
            (vel[1 + 3 * j] * c_explic + (frc[1 + 3 * j] + fln2) * wfac) * c_implic;
          vel[2 + 3 * j] = 
            (vel[2 + 3 * j] * c_explic + (frc[2 + 3 * j] + fln3) * wfac) * c_implic;
        }
      }
#ifdef TEST

#ifdef L_W
  lwpf_stop(T3);
  lwpf_start(A_P);
#endif

      pe_put( lp.vel + 3 * jst, vel, sizeof(double) * 3 * j_size);
      pe_syn();


#ifdef L_W
  lwpf_stop(A_P);
#endif


    }
#endif

  }
  else
  {
#ifdef TEST
    for( jst = _MYID * READ_SIZE; jst < atm_cnt; jst += 64 * READ_SIZE)
    {
      jed = jst + READ_SIZE;
      if( jed > atm_cnt)
        jed = atm_cnt;
        j_size = jed - jst;

#ifdef L_W
  lwpf_start(A_G);
  
  lwpf_start(T2);
#endif

      pe_get( lp.vel + 3 * jst, vel, sizeof(double) * 3 * j_size);
      pe_get( lp.frc + 3 * jst, frc, sizeof(double) * 3 * j_size);
      pe_get( lp.mass + jst, mass, sizeof(double) * j_size);
      pe_get( lp.mass_inv + jst, mass_inv, sizeof(double) * j_size);
      pe_get( lp.gbl_atm_owner_map + jst, gbl_atm_owner_map, sizeof(int) * j_size);
      pe_syn();

#ifdef L_W
  lwpf_stop(A_G);
  
  lwpf_stop(T2);
#endif
      int cnt = 0;


#endif

      for(j = 0; j < j_size; j++ )
      {
        if(gbl_atm_owner_map[j] == mytaskid)
        {
          cnt++;
          double wfac = mass_inv[j] * dtx;
          double aamass = mass[j];
          double rsd = sdfac * sqrt(aamass);
          double fln1 =  randn(0.0, rsd); 
          double fln2 =  randn(0.0, rsd);
          double fln3 =  randn(0.0, rsd);
          vel[0 + 3 * j] = 
            (vel[0 + 3 * j] * c_explic + (frc[0 + 3 * j] + fln1) * wfac) * c_implic; 
          vel[1 + 3 * j] = 
            (vel[1 + 3 * j] * c_explic + (frc[1 + 3 * j] + fln2) * wfac) * c_implic;
          vel[2 + 3 * j] = 
            (vel[2 + 3 * j] * c_explic + (frc[2 + 3 * j] + fln3) * wfac) * c_implic;
        }
        else
        {
          double rsd = 1.0;
          randn(0.0, rsd); 
          randn(0.0, rsd);
          randn(0.0, rsd);
        }
      }

#ifdef TEST

#ifdef L_W
  lwpf_start(A_P);
#endif
      if(cnt)
        pe_put( lp.vel + 3 * jst, vel, sizeof(double) * 3 * j_size);
      pe_syn();
#ifdef L_W
  lwpf_stop(A_P);
#endif

    }
#endif
  }

#ifdef L_W
  lwpf_stop(ALL);
#endif

}
#endif






