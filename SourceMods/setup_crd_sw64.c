
typedef struct atm_sort_rec
{
  int key;
  int atm;
}atm_sort_rec;

typedef struct setup_atm_lst_param
{
  int atm_cnt;
  double *fraction;
  int cit_tbl_x_dim;
  int cit_tbl_y_dim;
  int cit_tbl_z_dim;
  atm_sort_rec *atm_lst;
  int *hilbert3d4;
  
}setup_atm_lst_param;

typedef struct cit_tbl_rec
{
  int img_lo;
  int img_hi;

}cit_tbl_rec;

typedef struct setup_map_param
{
  int atm_cnt;
  int *atm_img_map;
  int *img_atm_map;
  atm_sort_rec *atm_lst;
  cit_tbl_rec *crd_idx_tbl;
  int *icit;
  int *map;

}setup_map_param;


#ifdef MPE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <athread.h>
#include <mpi.h>


extern SLAVE_FUN(setup_atm_lst_c_para)(void *);
extern SLAVE_FUN(setup_map_c_para)(void *);

void setup_map_c_(int *atm_cnt_,
                  int *atm_img_map_,
                  int *img_atm_map_,
                  atm_sort_rec *atm_lst_,
                  cit_tbl_rec *crd_idx_tbl
                  )
{  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  int atm_cnt = *atm_cnt_;
  int *atm_img_map = atm_img_map_ - 1;
  int *img_atm_map = img_atm_map_ - 1;
  atm_sort_rec *atm_lst = atm_lst_ - 1;
  int i,j,k,atm_id,citno;

  int icit = -1;

  //int *map;
  //map = (int *)malloc(sizeof(int) * (atm_cnt + 64));
  //if(map == NULL) printf("map: Memory Failed\n");

  setup_map_param param;
  param.atm_cnt     = atm_cnt;
  param.atm_img_map = atm_img_map;
  param.img_atm_map = img_atm_map;
  param.atm_lst     = atm_lst;
  param.crd_idx_tbl = crd_idx_tbl;
  param.icit        = &icit;
  //param.map         = map;

#define SUNWAY
#ifdef SUNWAY
  double time1,time2,time3;

  time1 = MPI_Wtime();
  for(atm_id = 1; atm_id <= atm_cnt; atm_id++)
  {
     //atm_img_map[atm_lst[atm_id].atm] = atm_id;
     //img_atm_map[atm_id] = atm_lst[atm_id].atm;
     citno = atm_lst[atm_id].key >> 12;
     if (icit != citno)
     {
        if (icit >= 0) 
          crd_idx_tbl[icit].img_hi = atm_id - 1;
        icit = citno;
        crd_idx_tbl[icit].img_lo = atm_id;
     }
  }
  crd_idx_tbl[icit].img_hi = atm_cnt;
 
  time2 = MPI_Wtime();

  if(athread_idle()==0) 
    athread_init();

  athread_spawn(setup_map_c_para,&param);
  athread_join();

  time3 = MPI_Wtime();

  //if(myrank == 0)
  //  printf("all time = %lf,one = %lf,two = %lf\n",time3-time1,time2-time1,time3-time2);
  
  //crd_idx_tbl[icit].img_hi = atm_cnt;
  
#else

  for(atm_id = 1; atm_id <= atm_cnt; atm_id++)
  {
     atm_img_map[atm_lst[atm_id].atm] = atm_id;
     img_atm_map[atm_id] = atm_lst[atm_id].atm;
     citno = atm_lst[atm_id].key >> 12;
     if (icit != citno)
     {
        if (icit >= 0) 
          crd_idx_tbl[icit].img_hi = atm_id - 1;
        icit = citno;
        crd_idx_tbl[icit].img_lo = atm_id;
     }
  }
  crd_idx_tbl[icit].img_hi = atm_cnt;

#endif

}

void setup_atm_lst_c_(int *atm_cnt_,
                  int *cit_tbl_x_dim_,
                  int *cit_tbl_y_dim_,
                  int *cit_tbl_z_dim_,
                  double *fraction,
                  int *hilbert3d4,
                  atm_sort_rec *atm_lst)
{
  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  int atm_cnt = *atm_cnt_;
  int cit_tbl_x_dim = *cit_tbl_x_dim_;
  int cit_tbl_y_dim = *cit_tbl_y_dim_;
  int cit_tbl_z_dim = *cit_tbl_z_dim_;
  
  int subcit_dim  = 4;
  int subcit_mask = 3;
  
  double scale_fac_x = (double)cit_tbl_x_dim;
  double scale_fac_y = (double)cit_tbl_y_dim;
  double scale_fac_z = (double)cit_tbl_z_dim;
  
  double sub_fac_x = (double)(cit_tbl_x_dim * subcit_dim);
  double sub_fac_y = (double)(cit_tbl_y_dim * subcit_dim);
  double sub_fac_z = (double)(cit_tbl_z_dim * subcit_dim);

  setup_atm_lst_param param;
  param.atm_cnt    = atm_cnt;
  param.fraction   = fraction;
  param.hilbert3d4 = hilbert3d4;
  param.atm_lst    = atm_lst;
  param.cit_tbl_x_dim = cit_tbl_x_dim;
  param.cit_tbl_y_dim = cit_tbl_y_dim;
  param.cit_tbl_z_dim = cit_tbl_z_dim;

#define SUNWAY
#ifdef SUNWAY

  if(athread_idle()==0) 
    athread_init();

  athread_spawn(setup_atm_lst_c_para,&param);
  athread_join();

#else

  int x_idx,y_idx,z_idx;
  int x_sub,y_sub,z_sub;
  int flat_sub,flat_idx;
 
  int i,j,k,atm_id;
  for(atm_id = 0; atm_id < atm_cnt; atm_id++)
  {
     x_idx = (int)(fraction[atm_id*3+0] * scale_fac_x);
     y_idx = (int)(fraction[atm_id*3+1] * scale_fac_y);
     z_idx = (int)(fraction[atm_id*3+2] * scale_fac_z);
     x_sub = ((int)(fraction[atm_id*3+0] * sub_fac_x)) & subcit_mask;
     y_sub = ((int)(fraction[atm_id*3+1] * sub_fac_y)) & subcit_mask;
     z_sub = ((int)(fraction[atm_id*3+2] * sub_fac_z)) & subcit_mask;
     flat_sub = (z_sub * subcit_dim + y_sub) * subcit_dim + x_sub;
     flat_idx = (z_idx * cit_tbl_y_dim + y_idx) * cit_tbl_x_dim + x_idx;
     atm_lst[atm_id].key = (flat_idx << 12) | hilbert3d4[flat_sub];
     atm_lst[atm_id].atm = atm_id + 1;
  }
#endif

}

#endif


#ifdef CPE
#include<slave.h>
#include<simd.h>
#include<dma_macros.h>

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

#define CACHE_HBIT 5
#define CACHE_SBIT 3
#define CACHE_LINESIZE (1 << CACHE_SBIT)
#define CACHE_LINECNT (1 << CACHE_HBIT)
#define CACHE_MMASK (CACHE_LINESIZE - 1)
#define CACHE_LMASK (CACHE_LINECNT - 1)

void write_cache_sw(int i,
                    int *tag_cache, 
                    double frc_cache[][CACHE_LINESIZE][3],
                    double *frc,
                    double fx, double fy, double fz, 
                    int *tot_cmp, int *miss_cache)
{
  dma_init();
  double *pfrc;
  if(tag_cache[(i >> CACHE_SBIT) & CACHE_LMASK] == -1)// get line number 
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

void write_cache(int i,
                    int *tag_cache, 
                    int atm_cache[][CACHE_LINESIZE],//take care of here
                    int *atm_img_map,
                    int val)
{
  dma_init();

  int cache_tag   = i >> CACHE_SBIT;
  int cache_line  = cache_tag & CACHE_LMASK;
  int cache_head  = i & ~CACHE_MMASK;
  int cache_off   = i & CACHE_MMASK;
  
  if(tag_cache[cache_line] == -1)// get line number 
  {
     pe_get(atm_img_map + cache_head,//get from the NO1 of the line 
            atm_cache[cache_line],
            sizeof(int) * CACHE_LINESIZE);
     dma_syn();
     tag_cache[cache_line] = cache_head;//mark the number of the head of the line
     
  }
  else if(tag_cache[cache_line] != cache_head )
  {
    pe_put(atm_img_map + tag_cache[cache_line], atm_cache[cache_line],
            sizeof(int) * CACHE_LINESIZE);
    pe_get(atm_img_map + cache_head, atm_cache[cache_line],
            sizeof(int) * CACHE_LINESIZE);
    dma_syn();
    tag_cache[cache_line] = cache_head;

  }
  atm_cache[cache_line][cache_off] = val;

}

void read_cache(int i,
                int *tag_cache,
                int cache[][CACHE_LINESIZE],
                int *read
                )
{
  
  dma_init();
  int cache_tag  = i >> CACHE_SBIT;
  int cache_line = i >> CACHE_SBIT & CACHE_LMASK;
  int cache_head = i & ~CACHE_MMASK;
  int cache_off  = i & CACHE_MMASK;

  if(tag_cache[cache_line] != cache_head)
  {
    pe_get(read + cache_head,cache[cache_line],sizeof(int)*CACHE_LINESIZE);
    dma_syn();
    tag_cache[cache_line] = cache_head;
  }
  return cache[cache_line][cache_off];
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


void setup_atm_lst_c_para(setup_atm_lst_param *pm)
{
  dma_init();
  setup_atm_lst_param lp;
  pe_get(pm,&lp,sizeof(setup_atm_lst_param));
  dma_syn();

  int atm_cnt          = lp.atm_cnt;
  int cit_tbl_x_dim    = lp.cit_tbl_x_dim;
  int cit_tbl_y_dim    = lp.cit_tbl_y_dim;
  int cit_tbl_z_dim    = lp.cit_tbl_z_dim;

  int subcit_dim  = 4;
  int subcit_mask = 3;
  
  double scale_fac_x = (double)cit_tbl_x_dim;
  double scale_fac_y = (double)cit_tbl_y_dim;
  double scale_fac_z = (double)cit_tbl_z_dim;
  
  double sub_fac_x = (double)(cit_tbl_x_dim * subcit_dim);
  double sub_fac_y = (double)(cit_tbl_y_dim * subcit_dim);
  double sub_fac_z = (double)(cit_tbl_z_dim * subcit_dim);
  
  int x_idx,y_idx,z_idx;
  int x_sub,y_sub,z_sub;
  int flat_sub,flat_idx;
  int i,j,atm_id;
  
  //double *fraction      = lp.fraction;
  //atm_sort_rec *atm_lst = lp.atm_lst;

  int hilbert3d4[64];
  pe_get(lp.hilbert3d4,hilbert3d4,sizeof(int)*64);
  dma_syn();
  
  double fraction[ISTEP * 3];
  atm_sort_rec atm_lst[ISTEP];

  int ist,ied,isz;
  
  for(ist = _MYID * ISTEP;ist < atm_cnt; ist += ISTEP *  64 )
  { 
    ied = ist + ISTEP;
    if(ied > atm_cnt)
      ied = atm_cnt;
    isz = ied - ist;
   
    pe_get(lp.fraction + ist * 3,fraction,sizeof(double)*isz*3);
    dma_syn();

    for(i = ist; i < ied; i++)
    //for(atm_id = ist; atm_id < ied; atm_id++)
    {
       atm_id = i - ist;
       x_idx = (int)(fraction[atm_id*3+0] * scale_fac_x);
       y_idx = (int)(fraction[atm_id*3+1] * scale_fac_y);
       z_idx = (int)(fraction[atm_id*3+2] * scale_fac_z);
       x_sub = ((int)(fraction[atm_id*3+0] * sub_fac_x)) & subcit_mask;
       y_sub = ((int)(fraction[atm_id*3+1] * sub_fac_y)) & subcit_mask;
       z_sub = ((int)(fraction[atm_id*3+2] * sub_fac_z)) & subcit_mask;
       flat_sub = (z_sub * subcit_dim + y_sub) * subcit_dim + x_sub;
       flat_idx = (z_idx * cit_tbl_y_dim + y_idx) * cit_tbl_x_dim + x_idx;
       atm_lst[atm_id].key = (flat_idx << 12) | hilbert3d4[flat_sub];
       atm_lst[atm_id].atm = i + 1;
    }
    pe_put(lp.atm_lst + ist,atm_lst,sizeof(atm_sort_rec)*isz);
    dma_syn();
  }

}

void setup_map_c_para(setup_map_param *pm)
{
  dma_init();
  setup_map_param lp;
  pe_get(pm,&lp,sizeof(setup_map_param));
  dma_syn();

  int atm_cnt               = lp.atm_cnt;
  int *atm_img_map          = lp.atm_img_map;
  //int *img_atm_map          = lp.img_atm_map;
  //atm_sort_rec *atm_lst     = lp.atm_lst;
  
  int img_atm_map[ISTEP];
  atm_sort_rec atm_lst[ISTEP];//yes

  int i,j,k,atm_id,citno,ii;
  
  int tag_cache[CACHE_LINECNT];
  int atm_cache[CACHE_LINECNT][CACHE_LINESIZE];

  for(i = 0;i < CACHE_LINECNT; i++)
      tag_cache[i] = -1;

  int icit = -1;
  
  int ist,ied,isz;

  for(ist = _MYID * ISTEP + 1;ist < atm_cnt + 1; ist += ISTEP *  64 )
  { 
    ied = ist + ISTEP;
    if(ied > atm_cnt + 1)
      ied = atm_cnt + 1;
    isz = ied - ist;
   
    pe_get(lp.atm_lst + ist,atm_lst,sizeof(atm_sort_rec)*isz);
    dma_syn();

    for(i = ist; i < ied; i++)
    {
          
     atm_id = i - ist;
     ii = atm_lst[atm_id].atm;
     
     atm_img_map[ii] = i;
     
     //write_cache(ii,tag_cache,atm_cache,atm_img_map,i);
     
     img_atm_map[atm_id] = ii;
     
    }
    pe_put(lp.img_atm_map + ist,img_atm_map,sizeof(int)*isz);
    dma_syn();

  }

  //for(i = 0;i < CACHE_LINECNT; i++)
  //{
  //  if(tag_cache[i] != -1)
  //  {  
  //    pe_put(atm_img_map + tag_cache[i], atm_cache[i],sizeof(int)*CACHE_LINESIZE);
  //    dma_syn();
  //  }
  //}

 
}

#endif














