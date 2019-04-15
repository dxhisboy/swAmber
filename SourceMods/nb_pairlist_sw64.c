typedef struct cit_tbl_rec {
  int img_lo, img_hi;
} cit_tbl_rec_t;
typedef struct adjust_imgcrds_arg{
  double (*img_crd)[3], (*saved_imgcrd)[3], (*crd)[3], (*saved_crd)[3];
  int *img_atm_map;
  char *used_img_map;
  cit_tbl_rec_t *cit;
  double box_del[3];
  int ncits, ntp, mytaskid, padding;
} adjust_imgcrds_arg_t;

#ifdef MPE
#include <athread.h>
extern void slave_adjust_imgcrds_para(adjust_imgcrds_arg_t *gl_arg);
void adjust_imgcrds_c_(double (*img_crd)[3], int *img_atm_map,
                       int *ncits_ptr, char *used_img_map, cit_tbl_rec_t *cit,
                       double (*saved_imgcrd)[3], double (*crd)[3], double (*saved_crd)[3],
                       int *ntp_ptr, double *box_del, int *mytaskid_ptr){
  int ncits = *ncits_ptr;
  int ntp = *ntp_ptr;
  int mytaskid = *mytaskid_ptr;
  adjust_imgcrds_arg_t arg;
  arg.ncits        = ncits;
  arg.ntp          = ntp;
  arg.mytaskid     = mytaskid;
  arg.img_crd      = img_crd;
  arg.saved_imgcrd = saved_imgcrd;
  arg.crd          = crd;
  arg.saved_crd    = saved_crd;
  arg.img_atm_map  = img_atm_map;
  arg.used_img_map = used_img_map;
  arg.cit          = cit;
  arg.box_del[0]   = box_del[0];
  arg.box_del[1]   = box_del[1];
  arg.box_del[2]   = box_del[2];

  athread_spawn(adjust_imgcrds_para, &arg);
  athread_join();
  return;
  int icit;
  if (!ntp){
    for (icit = 0; icit < ncits; icit ++){
      if (used_img_map[icit]){
        int img_i;
        int cit_lo = cit[icit].img_lo;
        int cit_hi = cit[icit].img_hi;
        for (img_i = cit_lo - 1; img_i < cit_hi; img_i ++){
          int atm_i = img_atm_map[img_i] - 1;
          img_crd[img_i][0] = crd[atm_i][0] + saved_imgcrd[img_i][0] - saved_crd[atm_i][0];
          img_crd[img_i][1] = crd[atm_i][1] + saved_imgcrd[img_i][1] - saved_crd[atm_i][1];
          img_crd[img_i][2] = crd[atm_i][2] + saved_imgcrd[img_i][2] - saved_crd[atm_i][2];
        }
      }
    }
  } else {
    for (icit = 0; icit < ncits; icit ++){
      if (used_img_map[icit]){
        int img_i;
        int cit_lo = cit[icit].img_lo;
        int cit_hi = cit[icit].img_hi;
        for (img_i = cit_lo - 1; img_i < cit_hi; img_i ++){
          int atm_i = img_atm_map[img_i] - 1;
          img_crd[img_i][0] = crd[atm_i][0] + (saved_imgcrd[img_i][0] - saved_crd[atm_i][0]) * box_del[0];
          img_crd[img_i][1] = crd[atm_i][1] + (saved_imgcrd[img_i][1] - saved_crd[atm_i][1]) * box_del[1];
          img_crd[img_i][2] = crd[atm_i][2] + (saved_imgcrd[img_i][2] - saved_crd[atm_i][2]) * box_del[2];
        }
      }
    }
  }
}
#endif

#ifdef CPE
#include <slave.h>
#include <dma.h>
#include <simd.h>
#include <dma_macros.h>

#define MAX_BKT_SIZE 64
#define CIT_STEP 8
#define CACHE_SBIT 3
#define CACHE_HBIT 9
#define CACHE_LINESIZE (1 << CACHE_SBIT)
#define CACHE_LINECNT  (1 << CACHE_HBIT - CACHE_SBIT)
#define CACHE_MMASK    (CACHE_LINESIZE - 1)
#define CACHE_LMASK    (CACHE_LINECNT - 1)
void adjust_imgcrds_para(adjust_imgcrds_arg_t *gl_arg){
  //if (_MYID) return;
  dma_init();
  adjust_imgcrds_arg_t arg;
  pe_get(gl_arg, &arg, sizeof(adjust_imgcrds_arg_t));
  dma_syn();

  int ncits = arg.ncits;
  int ntp   = arg.ntp;
  int mytaskid = arg.mytaskid;
  cit_tbl_rec_t *gbl_cit        = arg.cit;
  char *gbl_used_img_map        = arg.used_img_map;
  int *gbl_img_atm_map          = arg.img_atm_map;
  double (*gbl_img_crd)[3]      = arg.img_crd;
  double (*gbl_saved_imgcrd)[3] = arg.saved_imgcrd;
  double (*gbl_crd)[3]          = arg.crd;
  double (*gbl_saved_crd)[3]    = arg.saved_crd;

  double box_del[3];
  box_del[0]   = arg.box_del[0];
  box_del[1]   = arg.box_del[1];
  box_del[2]   = arg.box_del[2];
  int atm_cachetag[CACHE_LINECNT];
  double crd_cache[CACHE_LINECNT][CACHE_LINESIZE][3], saved_crd_cache[CACHE_LINECNT][CACHE_LINESIZE][3];
  double img_crd[MAX_BKT_SIZE][3], saved_imgcrd[MAX_BKT_SIZE][3];
  int img_atm_map[MAX_BKT_SIZE];
  char used_img_map_base[MAX_BKT_SIZE], *used_img_map;
  cit_tbl_rec_t cit[CIT_STEP];
  int i, icit, cit_st, cit_ed;
  for (i = 0; i < CACHE_LINECNT; i ++){
    atm_cachetag[i] = -1;
  }
  if (!ntp){
    box_del[0] = 1;
    box_del[1] = 1;
    box_del[2] = 1;
  }
  for (cit_st = _MYID * CIT_STEP; cit_st < ncits; cit_st += 64 * CIT_STEP){
    cit_ed = cit_st + CIT_STEP;
    if (cit_ed > ncits) cit_ed = ncits;
    int cit_sz = cit_ed - cit_st;
    /* if (mytaskid == 47 && cit_st == 64){ */
    /*   printf("%d %d\n", cit_st, cit_ed); */
    /* } */

    pe_get(gbl_cit          + cit_st, cit         , cit_sz * sizeof(cit_tbl_rec_t));
    dma_syn();
    for (icit = 0; icit < cit_sz; icit ++){
      int img_i;
      int cit_lo = cit[icit].img_lo - 1;
      int cit_hi = cit[icit].img_hi;
      /* if (mytaskid == 47 && icit + cit_st == 71){ */
      /*   printf("%d %d\n", cit_lo, cit_hi); */
      /* } */

      if (cit_lo >= cit_hi) continue;

      int img_map_base = cit_lo & ~3;
      int img_map_size = (cit_hi - img_map_base + 3) & ~3;
      pe_get(gbl_used_img_map + img_map_base, used_img_map_base, sizeof(char) * img_map_size);
      used_img_map = used_img_map_base + (cit_lo & 3);
      dma_syn();

      int upd_cnt = 0;
      for (img_i = 0; img_i < cit_hi - cit_lo; img_i ++){
        upd_cnt += used_img_map[img_i];
      }
      if (!upd_cnt) continue;
      pe_get(gbl_saved_imgcrd + cit_lo, saved_imgcrd, sizeof(double) * 3 * (cit_hi - cit_lo));
      pe_get(gbl_img_crd      + cit_lo, img_crd     , sizeof(double) * 3 * (cit_hi - cit_lo));
      pe_get(gbl_img_atm_map  + cit_lo, img_atm_map , sizeof(int)        * (cit_hi - cit_lo));
      dma_syn();
      /* if (mytaskid == 5 && _MYID == 0) printf("%d\n", gbl_used_img_map[81542]); */

      for (img_i = cit_lo; img_i < cit_hi; img_i ++){
        int img_off = img_i - cit_lo;
        /* if (mytaskid == 62 && img_i == 976040)  */
        /*   printf("!!!: %d %d %d %d %d %d\n", used_img_map[img_off], cit_lo, cit_hi, img_map_base, img_map_size, icit + cit_st); */

        if (!used_img_map[img_off]) continue;
        upd_cnt ++;
        int atm_i = img_atm_map[img_off] - 1;
        int atm_tag  = atm_i >> CACHE_SBIT;
        int atm_line = atm_tag & CACHE_LMASK;
        int atm_off  = atm_i & CACHE_MMASK;
        if (atm_cachetag[atm_line] != atm_tag){
          pe_get(gbl_crd       + (atm_i & ~CACHE_MMASK), crd_cache      [atm_line], sizeof(double) * 3 * CACHE_LINESIZE);
          pe_get(gbl_saved_crd + (atm_i & ~CACHE_MMASK), saved_crd_cache[atm_line], sizeof(double) * 3 * CACHE_LINESIZE);
          dma_syn();
        }
        /* if (mytaskid == 62 && img_i == 976040) { */
        /*   printf("%d %d %d\n", atm_i, atm_tag, atm_line, atm_off); */
        /*   printf("%f %f %f\n", crd_cache[atm_line][atm_off][0], crd_cache[atm_line][atm_off][1], crd_cache[atm_line][atm_off][2]); */
        /*   printf("%f %f %f\n", saved_crd_cache[atm_line][atm_off][0], saved_crd_cache[atm_line][atm_off][1], saved_crd_cache[atm_line][atm_off][2]); */
        /*   printf("%f %f %f\n", saved_imgcrd[img_off][0], saved_imgcrd[img_off][1], saved_imgcrd[img_off][2]); */
        /* } */

        img_crd[img_off][0] = crd_cache[atm_line][atm_off][0] + (saved_imgcrd[img_off][0] - saved_crd_cache[atm_line][atm_off][0]) * box_del[0];
        img_crd[img_off][1] = crd_cache[atm_line][atm_off][1] + (saved_imgcrd[img_off][1] - saved_crd_cache[atm_line][atm_off][1]) * box_del[1];
        img_crd[img_off][2] = crd_cache[atm_line][atm_off][2] + (saved_imgcrd[img_off][2] - saved_crd_cache[atm_line][atm_off][2]) * box_del[2];

      }
      pe_put(gbl_img_crd + cit_lo, img_crd, sizeof(double) * 3 * (cit_hi - cit_lo));
    }
  }

}
#endif
