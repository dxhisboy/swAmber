#include <params.h>

#ifdef MPE
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

extern SLAVE_FUN(pbc_CPE)(void *);

void get_fract_crds_sunway_(int    *atm_cnt_t,
                            double *fraction,
                            double *crd,
                            int    *is_orthog_t,
                            double *recip)
{
 
 Pbc param;
 param.atm_cnt         = *atm_cnt_t;
 param.fraction        = fraction;
 param.crd             = crd;
 param.is_orthog       = *is_orthog_t;
 param.recip           = recip;
 int i,j,k;

 int atm_cnt   = *atm_cnt_t;
 int is_orthog = *is_orthog_t;

 double recip_11,recip_22,recip_33;
 double f1,f2,f3;

 //double fraction_test[atm_cnt * 3];
 //param.fraction = fraction_test;

#define SUNWAY
#ifdef SUNWAY
  
  //printf("get_fract_cpe\n");
  
  if(athread_idle()==0) 
    athread_init();

  athread_spawn(pbc_CPE,&param);
  athread_join();

#else
 
 //printf("get_fract_mpe\n");
 if (is_orthog != 0)
 {
   recip_11 = recip[0];
   recip_22 = recip[4];
   recip_33 = recip[8];

   for( i = 0; i < atm_cnt; i++)
   {
     f1 = recip_11 * crd[i*3+0]; 
     f2 = recip_22 * crd[i*3+1]; 
     f3 = recip_33 * crd[i*3+2]; 
     //fraction[i*3+0] = f1 - (double)((int)(f1+0.5)) + 0.5;
     //fraction[i*3+1] = f2 - (double)((int)(f2+0.5)) + 0.5;
     //fraction[i*3+2] = f3 - (double)((int)(f3+0.5)) + 0.5;
   
     fraction[i*3+0] = f1 -round(f1) + 0.5;
     fraction[i*3+1] = f2 -round(f2) + 0.5;
     fraction[i*3+2] = f3 -round(f3) + 0.5;
     
     if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
     if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
     if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
     if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
     if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
     if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;
          
   }

 }else{

   for( i = 0; i < atm_cnt; i++)
   {
     f1 = crd[i*3+0] * recip[0] + crd[i*3+1] * recip[1] + 
          crd[i*3+2] * recip[2];                    
                                                   
     f2 = crd[i*3+0] * recip[3] + crd[i*3+1] * recip[4] + 
          crd[i*3+2] * recip[5];                    
                                                   
     f3 = crd[i*3+0] * recip[6] + crd[i*3+1] * recip[7] + 
          crd[i*3+2] * recip[8];
                
     //fraction[i*3+0] = f1 - (double)(int)(f1+0.5) + 0.5;
     //fraction[i*3+1] = f2 - (double)(int)(f2+0.5) + 0.5;
     //fraction[i*3+2] = f3 - (double)(int)(f3+0.5) + 0.5;
     
     fraction[i*3+0] = f1 -round(f1) + 0.5;
     fraction[i*3+1] = f2 -round(f2) + 0.5;
     fraction[i*3+2] = f3 -round(f3) + 0.5;
     
     if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
     if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
     if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
     if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
     if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
     if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;

   }

 }

 for( i = 0; i < atm_cnt; i++)
 {  
   if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
   if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
   if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
   if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
   if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
   if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;
 }

#endif 

}

#endif


#ifdef CPE

#include <stdio.h>
#include <stdlib.h>
#define ISTEP  32
#define SNSTEP 128

void pbc_CPE(Pbc *param)
{
  
  Pbc lp;
  pe_init();
  pe_get(param,&lp,sizeof(Pbc));
  pe_syn();

  int atm_cnt = lp.atm_cnt;
  
  //double *fraction = lp.fraction;
  //double *crd      = lp.crd;
  //double *recip    = lp.recip;
  int    is_orthog = lp.is_orthog;
 
  double recip_11,recip_22,recip_33;
  double f1,f2,f3;
  
  int i,j,k,ii;
  int ist,ied,isz;

  double recip[9];
  double crd[ISTEP*3];
  double fraction[ISTEP*3];

  pe_get(lp.recip,recip,sizeof(double)*9);
  pe_syn();
  
  if (is_orthog != 0)
  {
    recip_11 = recip[0];
    recip_22 = recip[4];
    recip_33 = recip[8];
    
    for(ist = _MYID * ISTEP;ist < atm_cnt;ist += ISTEP * 64)
    { 
    
      ied = ist + ISTEP;
      if(ied > atm_cnt)
        ied = atm_cnt;
      isz = ied - ist;

      pe_get(lp.crd + ist * 3,crd,sizeof(double)*isz*3);
      pe_syn();
    
      for( i = 0; i < isz; i++)
      {
        f1 = recip_11 * crd[i*3+0]; 
        f2 = recip_22 * crd[i*3+1]; 
        f3 = recip_33 * crd[i*3+2]; 
        //fraction[i*3+0] = f1 - (double)((int)(f1+0.5)) + 0.5;
        //fraction[i*3+1] = f2 - (double)((int)(f2+0.5)) + 0.5;
        //fraction[i*3+2] = f3 - (double)((int)(f3+0.5)) + 0.5;
      
        fraction[i*3+0] = f1 -round(f1) + 0.5;
        fraction[i*3+1] = f2 -round(f2) + 0.5;
        fraction[i*3+2] = f3 -round(f3) + 0.5;
        
        if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
        if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
        if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
        if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
        if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
        if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;
       
       }
       pe_put(lp.fraction + ist * 3,fraction,sizeof(double)*isz*3);
       pe_syn();
    }
  }else{

    for(ist = _MYID * ISTEP;ist < atm_cnt;ist += ISTEP * 64)
    { 
      ied = ist + ISTEP;
      if(ied > atm_cnt)
        ied = atm_cnt;
      isz = ied - ist;
      
      pe_get(lp.crd + ist * 3,crd,sizeof(double)*isz*3);
      pe_syn();
      
      for( i = 0; i < isz; i++)
      {
        f1 = crd[i*3+0] * recip[0] + crd[i*3+1] * recip[1] + 
             crd[i*3+2] * recip[2];                    
                                                      
        f2 = crd[i*3+0] * recip[3] + crd[i*3+1] * recip[4] + 
             crd[i*3+2] * recip[5];                    
                                                      
        f3 = crd[i*3+0] * recip[6] + crd[i*3+1] * recip[7] + 
             crd[i*3+2] * recip[8];
                   
        //fraction[i*3+0] = f1 - (double)(int)(f1+0.5) + 0.5;
        //fraction[i*3+1] = f2 - (double)(int)(f2+0.5) + 0.5;
        //fraction[i*3+2] = f3 - (double)(int)(f3+0.5) + 0.5;
        
        fraction[i*3+0] = f1 -round(f1) + 0.5;
        fraction[i*3+1] = f2 -round(f2) + 0.5;
        fraction[i*3+2] = f3 -round(f3) + 0.5;
      
        if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
        if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
        if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
        if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
        if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
        if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;


      }
      pe_put(lp.fraction + ist * 3,fraction,sizeof(double)*isz*3);
      pe_syn();
    }
  }

  //for(ist = _MYID * ISTEP;ist < atm_cnt;ist += ISTEP * 64)
  //{ 
  //  ied = ist + ISTEP;
  //  if(ied > atm_cnt)
  //    ied = atm_cnt;
  //  isz = ied - ist;
  //  for( i = ist; i < ied; i++)
  //  {  
  //    if (fraction[i*3+0] <  0) fraction[i*3+0] = fraction[i*3+0] + 1;
  //    if (fraction[i*3+0] >= 1) fraction[i*3+0] = fraction[i*3+0] - 1;
  //    if (fraction[i*3+1] <  0) fraction[i*3+1] = fraction[i*3+1] + 1;
  //    if (fraction[i*3+1] >= 1) fraction[i*3+1] = fraction[i*3+1] - 1;
  //    if (fraction[i*3+2] <  0) fraction[i*3+2] = fraction[i*3+2] + 1;
  //    if (fraction[i*3+2] >= 1) fraction[i*3+2] = fraction[i*3+2] - 1;
  //  }
  //}
}


#endif
