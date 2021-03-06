#ifndef PARAM_H
#define PARAM_H
#ifdef MPE
#include <athread.h>
#endif

#ifdef CPE
#include <slave.h>
#include <simd.h>

#define pe_init() volatile int reply = 0; int pe_cnt = 0;
#define pe_get(mem, ldm, size) {athread_get(PE_MODE, (mem), (ldm), (size), (void*)&reply, 0, 0, 0); pe_cnt ++;}
#define pe_put(mem, ldm, size) {athread_put(PE_MODE, (ldm), (mem), (size), (void*)&reply, 0, 0); pe_cnt ++;}
#define pe_syn() {while (reply != pe_cnt); asm volatile("memb");}
//#define pe_syn() {asm volatile("memb");}

#endif

typedef struct lang_set_param{
  double dtx, c_explic;
  double fln1, wfac;
  double c_implic, sdfac;
  int atm_cnt, no_ntt3_sync;
  int mytaskid;

  double *vel, *frc, *mass, *mass_inv;
  int *gbl_atm_owner_map;
}lang_set_param;
  
typedef struct Pbc{
  int atm_cnt;
  double *fraction;
  double *crd;
  int is_orthog;
  double *recip;
}Pbc;

#ifdef CPE
void reg_reduce_inplace_doublev4(doublev4 *arr, int len)
{
  //if(_MYID == 0)
  //printf("1-------------\n");
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
  //if(_MYID == 0)
  //  printf("2-------------\n");

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
  //if(_MYID == 0)
  //  printf("3-------------\n");
}
#endif

#endif
