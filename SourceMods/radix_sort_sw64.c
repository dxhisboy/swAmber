#define RADIX_BITS 8
#define RADIX_MASK ((1 << RADIX_BITS) - 1)
#define NBKTS (1 << RADIX_BITS)
typedef struct kvpair{
  int k, v;
} kvpair_t;


typedef struct radix_distrib_param
{
  kvpair_t *data;
  kvpair_t *bkt_data_stor;
  int *head; 
  int *cnt;
  int n;
  int shift;
  int myrank;
}radix_distrib_param;

typedef struct radix_gather_param
{
  kvpair_t *data; 
  kvpair_t **bkt_data;
  int *bkt_head;
  int *bkt_head_distrib;
  int myrank;
}radix_gather_param;

#ifdef MPE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<athread.h>
#include<mpi.h>
extern SLAVE_FUN(distrib_para)();
extern SLAVE_FUN(gather_para)();
void radix_sort_distrib(int shift, 
                        kvpair_t *data, 
                        kvpair_t *bkt, 
                        int *head, 
                        int *cnt, 
                        int st, 
                        int ed)
{

  int i;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //if(myrank == 0) 
  //  printf("distrib-----on-----mpe\n");

 ///////////////////
  kvpair_t tmp_cache[NBKTS][15];
  int cache_cnt[NBKTS];
  for(i = 0; i < NBKTS; i++)
  {
    cache_cnt[i] = 0;
  }

  for (i = 0; i < NBKTS; i ++)
  {
    cnt[i] = 0;
  }
  for (i = st; i < ed; i ++)
  {
    int bktid = (data[i].k >> shift) & RADIX_MASK;
    cnt[bktid] ++;
  }
  int accu_cnt = 0;
  for (i = 0; i < NBKTS; i ++)
  {
    head[i] = accu_cnt;
    accu_cnt += cnt[i];
    cnt[i] = 0;
  }
  int ii;
  head[NBKTS] = accu_cnt;

  for (i = st; i < ed; i ++)
  {
      int bktid = (data[i].k >> shift) & RADIX_MASK;
      int index = head[bktid] + cnt[bktid];
      bkt[index] = data[i];
      cnt[bktid]++;
      
  }//for-i

//  #define DSTEP 64
//  //for (i = st; i < ed; i ++)
//  for (i = st; i < ed; i+= DSTEP)
//  {
//    int ist, ied, isz, ii;
//    ist = i;
//    ied = i + DSTEP;
//    if(ied  > ed)
//      ied = ed;
//    isz = ied - ist;
//    
//    for(ii = ist; ii < ied; ii++)
//    {
//      //int bktid = (data[i].k >> shift) & RADIX_MASK;
//      int bktid = (data[ii].k >> shift) & RADIX_MASK;
//      int index = head[bktid] + cnt[bktid];
//      
//      //tmp_cache[bktid][cache_cnt[bktid]] = data[i];
//      tmp_cache[bktid][cache_cnt[bktid]] = data[ii];
//      cache_cnt[bktid]++;
//      if(cache_cnt[bktid] == 15)
//      {
//        int b;
//        for(b = 0; b < 15; b++)
//        {
//          bkt[index + b] = tmp_cache[bktid][b];
//        }
//        cache_cnt[bktid] = 0;
//        cnt[bktid] +=15;
//      }
//    }//for-ii
//  }//for-i
//
//  int bid; 
//  for(bid = 0; bid < NBKTS; bid++)
//  {
//    if(cache_cnt[bid] != 0)
//    {
//      for(ii = 0; ii < cache_cnt[bid]; ii++)
//      {
//        bkt[head[bid] + cnt[bid] + ii] = tmp_cache[bid][ii];
//      }
//      cnt[bid] += cache_cnt[bid];
//    }
//  }

  //for (i = st; i < ed; i ++)
  //{
  //  int bktid = (data[i].k >> shift) & RADIX_MASK;
  //  int index = head[bktid] + cnt[bktid];
  //  
  //  tmp_cache[bktid][cache_cnt[bktid]] = data[i];
  //  cache_cnt[bktid]++;
  //  if(cache_cnt[bktid] == 15)
  //  {
  //    for(ii = 0; ii < 15; ii++)
  //    {
  //      bkt[index + ii] = tmp_cache[bktid][ii];
  //    }
  //    cache_cnt[bktid] = 0;
  //    cnt[bktid] +=15;
  //  }
  //}
  //int bid; 
  //for(bid = 0; bid < NBKTS; bid++)
  //{
  //  if(cache_cnt[bid] != 0)
  //  {
  //    for(ii = 0; ii < cache_cnt[bid]; ii++)
  //    {
  //      bkt[head[bid] + cnt[bid] + ii] = tmp_cache[bid][ii];
  //    }
  //    cnt[bid] += cache_cnt[bid];
  //  }
  //}




}
void gather_bkt_data(kvpair_t *data, 
                      kvpair_t *bkt_data[64], 
                      int *bkt_head, 
                      int (*bkt_head_distrib)[NBKTS + 1], 
                      int st, 
                      int ed)
{
  int i, peid;
  for (i = st; i < ed; i ++)
  {
    int ndat_bkt = 0;
    for (peid = 0; peid < 64; peid ++)
    {
      int ndat_pebkt = bkt_head_distrib[peid][i + 1] - bkt_head_distrib[peid][i];
      memcpy(data + bkt_head[i] + ndat_bkt, bkt_data[peid] + bkt_head_distrib[peid][i], ndat_pebkt * sizeof(kvpair_t));
      ndat_bkt += ndat_pebkt;
    }
  }
}

//void gather_bkt_data(kvpair_t *data, 
//                      kvpair_t **bkt_data, 
//                      int *bkt_head, 
//                      int *bkt_head_distrib, 
//                      int st, 
//                      int ed)
//{
//  int i, peid;
//  
//  int myrank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  //if(myrank == 0) 
//  //  printf("gather-----on-----mpe\n");
//
//  ///////(1)/////////
//  //for (i = st; i < ed; i ++)
//  //{
//  //  //int ndat_bkt = 0;
//  //  for (peid = 0; peid < 64; peid ++)
//  //  {
//  //    int ndat_pebkt = bkt_head_distrib[peid * (NBKTS+1) + i + 1] 
//  //                  - bkt_head_distrib[peid * (NBKTS+1) + i];//< 3000
//  //    //memcpy(data + bkt_head[i] + ndat_bkt, 
//  //    memcpy(data + bkt_head[i] + nn[i-st], 
//  //          bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i], 
//  //          ndat_pebkt * sizeof(kvpair_t));
//  //    nn[i-st] += ndat_pebkt;
//  //    //ndat_bkt += ndat_pebkt;
//  //  }
//  //}
//
//  ///////(2)/////////
//  //for (i = st; i < ed; i ++)
//  //{
//  //  int ndat_bkt = 0;
//  //  for (peid = 0; peid < 64; peid ++)
//  //  {
//  //    int ndat_pebkt = bkt_head_distrib[peid * (NBKTS+1) + i + 1] 
//  //                  - bkt_head_distrib[peid * (NBKTS+1) + i];//< 3000
//  //    int j;
//  //    for(j = 0; j < ndat_pebkt; j++)
//  //    {
//  //      data[bkt_head[i] + j] = 
//  //        *(bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i] + j); 
//  //    }
//  //    //memcpy(data + bkt_head[i] + ndat_bkt, 
//  //    //      bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i], 
//  //    //      ndat_pebkt * sizeof(kvpair_t));
//  //    ndat_bkt += ndat_pebkt;
//  //  }
//  //}
//
//  ///////(2)/////////
//  for (i = st; i < ed; i ++)
//  {
//    int ndat_bkt = 0;
//    for (peid = 0; peid < 64; peid ++)
//    {
//      int ndat_pebkt = bkt_head_distrib[peid * (NBKTS+1) + i + 1] 
//                    - bkt_head_distrib[peid * (NBKTS+1) + i];//< 3000
//      int j, jst, jed, jsz;
//      //for(j = 0; j < ndat_pebkt; j++)
//      #define JSTEP 3
//      for(j = 0; j < ndat_pebkt; j+=JSTEP)
//      {
//        jst = j;
//        jed = jst + JSTEP;
//        if(jed > ndat_pebkt)
//          jed = ndat_pebkt;
//        jsz = jed - jst;
//        
//        data[bkt_head[i] + ndat_bkt + jst] = 
//          (bkt_data[peid])[bkt_head_distrib[peid * (NBKTS+1) + i] + jst]; 
//
//        data[bkt_head[i] + ndat_bkt + jst+1] = 
//          (bkt_data[peid])[bkt_head_distrib[peid * (NBKTS+1) + i] + jst+1]; 
//
//        data[bkt_head[i] + ndat_bkt + jst+2] = 
//          (bkt_data[peid])[bkt_head_distrib[peid * (NBKTS+1) + i] + jst+2]; 
//          //*(bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i] + j); 
//      }
//      //memcpy(data + bkt_head[i] + ndat_bkt, 
//      //      bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i], 
//      //      ndat_pebkt * sizeof(kvpair_t));
//
//      ndat_bkt += ndat_pebkt;
//    }
//  }
//
//
//  ///////(3)/////////
////  int j; 
////  for (peid = 0; peid < 64; peid ++)
////  {
////    for (i = st; i < ed; i ++)
////    {
////      int ndat_pebkt = bkt_head_distrib[peid * (NBKTS+1) + i + 1] 
////                    - bkt_head_distrib[peid * (NBKTS+1) + i];//< 3000
////      for(j = 0; j < ndat_pebkt; j++)
////      {
////        data[bkt_head[i] + nn[i-st] + j]=
////        bkt_data[peid][bkt_head_distrib[peid * (NBKTS+1) + i] + j]
////
////      }
////      nn[i-st] += ndat_pebkt;
////      //if(ndat_pebkt > 0)
////      //{
////      //  memcpy(data + bkt_head[i] + nn[i-st], 
////      //        bkt_data[peid] + bkt_head_distrib[peid * (NBKTS+1) + i], 
////      //        ndat_pebkt * sizeof(kvpair_t));
////      //  nn[i-st] += ndat_pebkt;
////      //}
////    }
////  }
//  
//}

void radix_sort_c_(kvpair_t *data, int *n_ptr)
{
  int n = *n_ptr;
  int key_max;
  int i, myid;
  int bkt_cnt_distrib[64][NBKTS];
  int bkt_head_distrib[64][NBKTS + 1];
  int bkt_cnt[NBKTS], bkt_head[NBKTS + 1];
  kvpair_t *bkt_data_stor = malloc(64 * (n + 512) * sizeof(kvpair_t));
  kvpair_t *bkt_data[64];
  for (myid = 0; myid < 64; myid ++)
  {
    bkt_data[myid] = bkt_data_stor + (n + 512) * myid;
  }
  key_max = 0;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  double time0, time1, time2, time3;
  double time4, time5, time6, time7;
  double time8, time9, time10, time11;
  time0 = MPI_Wtime();
  for (i = 0; i < n; i ++)
    if (data[i].k > key_max) key_max = data[i].k;
  int shift;
  
  time1 = MPI_Wtime();
  if(athread_idle() == 0)
    athread_init();

  radix_distrib_param pm;
  for (shift = 0; key_max >> shift; shift += RADIX_BITS)
  {
    
    time3 = MPI_Wtime();
    pm.n              = n;
    pm.shift          = shift;
    pm.data           = data;
    pm.bkt_data_stor  = bkt_data_stor; 
    pm.head           = &bkt_head_distrib[0][0];
    pm.cnt            = &bkt_cnt_distrib[0][0];
    pm.myrank         = myrank;

#define DISTRIB
#ifdef DISTRIB  
    athread_spawn(distrib_para, &pm);
    athread_join();

#else
    for (myid = 0; myid < 64; myid ++)
    {
      int ndat_pe = (n + 63) >> 6;
      int dat_st = myid * ndat_pe;
      int dat_ed = dat_st + ndat_pe;
      if (dat_ed > n) dat_ed = n;
      //radix_sort_distrib(shift, 
      //                  data, 
      //                  bkt_data_stor + (n + 512) * myid, 
      //                  bkt_head_distrib[myid], 
      //                  bkt_cnt_distrib[myid], 
      //                  dat_st, dat_ed);
      
      radix_sort_distrib(shift, 
                        data, 
                        bkt_data_stor + (n + 512) * myid, 
                        &bkt_head_distrib[0][0] +myid * (NBKTS + 1), 
                        &bkt_cnt_distrib[0][0]  +myid * NBKTS, 
                        dat_st, dat_ed);
    }
#endif

    time4 = MPI_Wtime();

    int accu_cnt = 0;
    for (i = 0; i < NBKTS; i ++)
    {
      bkt_cnt[i] = 0;
      bkt_head[i] = accu_cnt;
      for (myid = 0; myid < 64; myid ++)
      {
        bkt_cnt[i] += bkt_cnt_distrib[myid][i];
      }
      accu_cnt += bkt_cnt[i];
    }
    bkt_head[NBKTS] = accu_cnt;
    
    time5 = MPI_Wtime();
    
    radix_gather_param pm2;
    pm2.data              = data;
    pm2.bkt_data          = bkt_data;
    pm2.bkt_head          = bkt_head; 
    pm2.bkt_head_distrib  = &bkt_head_distrib[0][0];
    pm2.myrank            = myrank;
#define GATHER
#ifdef GATHER
    athread_spawn(gather_para, &pm2);
    athread_join();

#else
    for (myid = 0; myid < 64; myid ++) 
    {
      int nbkt_pe = (NBKTS + 63) >> 6;
      int bkt_st = myid * nbkt_pe;
      int bkt_ed = bkt_st + nbkt_pe;
      if (bkt_ed > NBKTS) bkt_ed = NBKTS;
      gather_bkt_data(data, 
                      bkt_data, 
                      bkt_head, 
                      &bkt_head_distrib[0][0], 
                      bkt_st, 
                      bkt_ed);

      //gather_bkt_data(data, 
      //                bkt_data, 
      //                bkt_head, 
      //                bkt_head_distrib, 
      //                bkt_st, 
      //                bkt_ed);
    }
#endif

    time6 = MPI_Wtime();
    //if(myrank == 0)
    //{
    //  printf("big-for = %.10f, %.10f, %.10f, %.10f\n", 
    //  time6-time3, time4-time3, time5-time4, time6-time5);
    //}
  }
  time2 = MPI_Wtime();
  //if(myrank == 0)
  //{
  //  printf("all = %.10f, %.10f, %.10f\n", time2-time0, time1-time0, time2-time1);
  //}
  free(bkt_data_stor);
}
#endif



#ifdef CPE
#include<slave.h>
#include<simd.h>
#include<dma_macros.h>
#define DSTEP 128
#define C_SZ  20
void distrib_para(radix_distrib_param *pm)
{
  dma_init();
  radix_distrib_param lp;
  pe_get(pm, &lp, sizeof(radix_distrib_param));
  dma_syn();

  int n           = lp.n;
  int shift       = lp.shift;
  kvpair_t *data  = lp.data;
  kvpair_t *bkt   = lp.bkt_data_stor + _MYID * (n + 512);
  int *head       = lp.head  + _MYID * (NBKTS + 1);
  int *cnt        = lp.cnt   + _MYID * NBKTS;
  
  int l_cnt[NBKTS], l_head[NBKTS+1];
  kvpair_t l_data[DSTEP];
  kvpair_t l_bkt[DSTEP];

  int ndat_pe = (n + 63) >> 6;
  int dat_st  = _MYID * ndat_pe;
  int dat_ed  = dat_st + ndat_pe;
  if (dat_ed > n) dat_ed = n;
  
  int i, j;
  intv8 zerov8;
  int *tmp = &zerov8;
  tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0;
  tmp[4] = tmp[5] = tmp[6] = tmp[7] = 0;
  
  kvpair_t tmp_cache[NBKTS][C_SZ];
  int cache_cnt[NBKTS];

  for(i = 0; i < NBKTS; i+= 8)
  {
    simd_store(zerov8, &l_cnt[i]);
    simd_store(zerov8, &cache_cnt[i]);
  }

  //for (i = dat_st; i < dat_ed; i ++)
  for (i = dat_st; i < dat_ed; i+=DSTEP)
  {
    int ist, ied, isz, ii;
    ist = i;
    ied = i + DSTEP;
    if(ied  > dat_ed)
      ied = dat_ed;
    isz = ied - ist;
    pe_get(data + ist, l_data, sizeof(kvpair_t)*isz);
    dma_syn();
    for(ii = ist; ii < ied; ii++)
    {
      int ioff = ii - ist;
      //int bktid = (data[i].k >> shift) & RADIX_MASK;
      int bktid = (l_data[ioff].k >> shift) & RADIX_MASK;
      l_cnt[bktid] ++;
    }
  }

  int accu_cnt = 0;
  for (i = 0; i < NBKTS; i ++)
  {
    l_head[i] = accu_cnt;
    accu_cnt += l_cnt[i];
    l_cnt[i] = 0;
  }
  l_head[NBKTS] = accu_cnt;

  //for (i = dat_st; i < dat_ed; i ++)
  for (i = dat_st; i < dat_ed; i+=DSTEP)
  {
    int ist, ied, isz, ii;
    ist = i;
    ied = i + DSTEP;
    if(ied  > dat_ed)
      ied = dat_ed;
    isz = ied - ist;
    pe_get(data + ist, l_data, sizeof(kvpair_t)*isz);
    dma_syn();

    for(ii = ist; ii < ied; ii++)
    {
      int ioff = ii - ist;
      int bktid = (l_data[ioff].k >> shift) & RADIX_MASK;
      int index = l_head[bktid] + l_cnt[bktid];
      
      //bkt[index] = l_data[ioff];
      //l_cnt[bktid] ++;

      tmp_cache[bktid][cache_cnt[bktid]] = l_data[ioff];
      cache_cnt[bktid]++;
      if(cache_cnt[bktid] == C_SZ)
      {
        pe_put(bkt+index, tmp_cache[bktid], sizeof(kvpair_t) * C_SZ);
        dma_syn();
        cache_cnt[bktid] = 0;
        l_cnt[bktid]+=C_SZ;
      }
    }//for-ii
  }//for-i

  int bid;
  for(bid = 0; bid < NBKTS; bid++)
  {
    if(cache_cnt[bid] != 0)
    {
      pe_put(bkt+(l_head[bid] + l_cnt[bid]), 
                  tmp_cache[bid],
                  sizeof(kvpair_t) * cache_cnt[bid]);
      dma_syn();
      l_cnt[bid] += cache_cnt[bid];
    }
  }

  pe_put(cnt, l_cnt, sizeof(int)*NBKTS);
  pe_put(head, l_head, sizeof(int)*(NBKTS+1));
  dma_syn();

}

#define JSTEP 128
#define PSTEP 4
void gather_para(radix_gather_param *pm)
{
  dma_init();
  radix_gather_param lp;
  pe_get(pm, &lp, sizeof(radix_gather_param));
  dma_syn();

  kvpair_t *data        = lp.data; 
  kvpair_t **bkt_data   = lp.bkt_data;
  int *bkt_head         = lp.bkt_head;
  int *bkt_head_distrib = lp.bkt_head_distrib;

  int nbkt_pe = (NBKTS + 63) >> 6;
  int bkt_st = _MYID* nbkt_pe;
  int bkt_ed = bkt_st + nbkt_pe;
  if (bkt_ed > NBKTS) bkt_ed = NBKTS;
  
  kvpair_t l_data[JSTEP];
  int l_bkt_head[NBKTS+1];
  kvpair_t *l_bkt_data[64];

  pe_get(bkt_head, l_bkt_head, sizeof(int) * (NBKTS+1));
  pe_get(bkt_data, l_bkt_data, sizeof(kvpair_t *) * 64);
  dma_syn();
  
  int i, peid, j;
  //for (i = bkt_st; i < bkt_ed; i ++)
  //{
  //  int ndat_bkt = 0;
  //  for (peid = 0; peid < 64; peid ++)
  //  {
  //    int ndat_pebkt = bkt_head_distrib[peid * (NBKTS+1) + i + 1] 
  //                  - bkt_head_distrib[peid * (NBKTS+1) + i];//< 3000
  //    int j, jst, jed, jsz;
  //    for(j = 0; j < ndat_pebkt; j+=JSTEP)
  //    {
  //      jst = j;
  //      jed = jst + JSTEP;
  //      if(jed > ndat_pebkt)
  //        jed = ndat_pebkt;
  //      jsz = jed - jst;
  //      if(jsz > 0)
  //      {
  //        pe_get((l_bkt_data[peid]) + (bkt_head_distrib[peid*(NBKTS+1)+i] + jst),
  //                l_data,
  //                sizeof(kvpair_t) * jsz);
  //        dma_syn();

  //        pe_put(data + (l_bkt_head[i] + ndat_bkt + jst),
  //                l_data,
  //                sizeof(kvpair_t) * jsz);
  //        dma_syn();
  //      }
  //    }//for-j
  //    ndat_bkt += ndat_pebkt;
  //  }//for-peid
  //}//for-i
 
  int ndat_cnt[NBKTS];
  int l_bkt_head_distrib[PSTEP * (NBKTS +1)];
  intv8 zerov8;
  int *tmp = &zerov8;
  tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0;
  tmp[4] = tmp[5] = tmp[6] = tmp[7] = 0;
  for(i = 0; i < NBKTS; i+= 8)
  {
    simd_store(zerov8, &ndat_cnt[i]);
  }


  //for (peid = 0; peid < 64; peid ++)
  for (peid = 0; peid < 64; peid += PSTEP)
  {
    int pp, pst, ped, psz, poff;
    pst = peid; 
    ped = pst + PSTEP;
    if(ped > 64) ped = 64;
    psz = ped - pst;
    pe_get(bkt_head_distrib + pst*(NBKTS+1), l_bkt_head_distrib, sizeof(int) * psz *(NBKTS+1));
    dma_syn();

    for(pp = pst; pp < ped; pp++)
    {
      poff = pp - pst;
      //int ndat_bkt = 0;
      for (i = bkt_st; i < bkt_ed; i ++)
      {
        int ndat_pebkt = l_bkt_head_distrib[poff * (NBKTS+1) + i + 1] 
                       - l_bkt_head_distrib[poff * (NBKTS+1) + i];
        int j, jst, jed, jsz;
        for(j = 0; j < ndat_pebkt; j+=JSTEP)
        {
          jst = j;
          jed = jst + JSTEP;
          if(jed > ndat_pebkt)
            jed = ndat_pebkt;
          jsz = jed - jst;
          if(jsz > 0)
          {
            pe_get((l_bkt_data[pp]) + (l_bkt_head_distrib[poff*(NBKTS+1)+i] + jst),
                    l_data,
                    sizeof(kvpair_t) * jsz);
            dma_syn();

            //pe_put(data + (l_bkt_head[i] + ndat_bkt + jst),
            pe_put(data + (l_bkt_head[i] + ndat_cnt[i-bkt_st] + jst),
                    l_data,
                    sizeof(kvpair_t) * jsz);
            dma_syn();
          }
        }//for-j
        //ndat_bkt += ndat_pebkt;
        ndat_cnt[i - bkt_st] += ndat_pebkt;
      }//for-peid
    }//for-pp
  }//for-i

}

#endif
