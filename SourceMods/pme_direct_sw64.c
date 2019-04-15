#ifdef MPE
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<athread.h>
#define LWPF_UNITS U(PME_SETUP_MASK)
#include <lwpf2/lwpf2.h>

#include<param.h>
extern SLAVE_FUN(pme_mask_setup_c_pre)();
extern SLAVE_FUN(pme_mask_setup_c_para)();
void wrap_c(int x, int xlim, int *ox_, int *tx_)
{
  int ox = *ox_;
  int tx = *tx_;
  ox = x;
  tx = 0;
  if(ox < 0)
  {
     ox = ox + xlim;
     tx = -1;
  }
  if(ox >= xlim)
  {
     ox = ox - xlim;
     tx = 1;
  }
  *ox_ = ox;
  *tx_ = tx;
}

void cit_wrap_c(int x, int y, int z, int *ox, int *oy, int *oz, int *tx, int *ty, int *tz, int xlim, int ylim, int zlim)
{
  wrap_c(x, xlim, ox, tx);
  wrap_c(y, ylim, oy, ty);
  wrap_c(z, zlim, oz, tz);
}

 void pme_mask_setup_c_(int *img_lo_, 
                    int *img_hi_, 
                    int *my_bkt_lo_,
                    int *my_bkt_hi_,
                    int *cit_tbl_x_dim_,
                    int *cit_tbl_y_dim_,
                    int *cit_tbl_z_dim_,
                    int *ndbkt_,
                    int *dbkt_,
                    listdata_rec_c *atm_maskdata_, 
                    int *atm_mask_, 
                    listdata_rec_c *img_maskdata_, 
                    int *img_mask_,
                    int *gbl_atm_img_map_, 
                    int *gbl_img_atm_map_,  
                    cit_tbl_rec_c *cit_save_,         
                    int *ncit_save_)
                    
 {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  int img_lo                    = *img_lo_;
  int img_hi                    = *img_hi_;
  int my_bkt_lo                 = *my_bkt_lo_;
  int my_bkt_hi                 = *my_bkt_hi_;
  int cit_tbl_x_dim             = *cit_tbl_x_dim_;
  int cit_tbl_y_dim             = *cit_tbl_y_dim_;
  int cit_tbl_z_dim             = *cit_tbl_z_dim_;
  int ndbkt                     = *ndbkt_;
  int *dbkt                     = dbkt_         - 3;
  listdata_rec_c *atm_maskdata  = atm_maskdata_ - 1;
  int *atm_mask                 = atm_mask_     - 1;
  listdata_rec_c *img_maskdata  = img_maskdata_ - 1;
  int *img_mask                 = img_mask_     - 1;
  int *gbl_atm_img_map          = gbl_atm_img_map_ - 1; 
  int *gbl_img_atm_map          = gbl_img_atm_map_ - 1;
  cit_tbl_rec_c *cit_save       = cit_save_;         
  int ncit_save                 = *ncit_save_;
  //int atm_cnt                   = *atm_cnt_;
  
  pme_mask_param pm;
  pm.img_lo                             = img_lo;
  pm.img_hi                             = img_hi;
  pm.my_bkt_lo                          = my_bkt_lo; 
  pm.my_bkt_hi                          = my_bkt_hi;
  pm.cit_tbl_x_dim                      = cit_tbl_x_dim;
  pm.cit_tbl_y_dim                      = cit_tbl_y_dim;
  pm.cit_tbl_z_dim                      = cit_tbl_z_dim;
  pm.ndbkt                              = ndbkt;
  pm.dbkt                               = dbkt;
  pm.atm_maskdata                       = atm_maskdata;
  pm.atm_mask                           = atm_mask;
  pm.img_maskdata                       = img_maskdata;
  pm.img_mask                           = img_mask;
  pm.gbl_atm_img_map                    = gbl_atm_img_map;
  pm.gbl_img_atm_map                    = gbl_img_atm_map;
  pm.cit_save                           = cit_save;
  pm.ncit_save                          = ncit_save;
  pm.myrank                             = myrank;
  //pm.atm_cnt                            = atm_cnt;
  
  int  idbkt, itran;
  int  dx, dy, dz, xhi, yhi;
  int  mx, my, mz, nx, ny, nz, tx, ty, tz, cx, cy, cz;
  int  icit, cit_lo, cit_hi, jcit;
  int  img_i, img_j, jhi, jlo;
  int  atm_i, i, img_off, atm_off;
  int  img_mask_cnt;
  int  mask27 = 0x07ffffff; 

  img_mask_cnt = 0;

  double time0, time1, time2, time3;
  time0 = MPI_Wtime();
  
  if(athread_idle() == 0)
    athread_init();
  
  for(img_i = img_lo; img_i <= img_hi; img_i ++)
  {
     atm_i = gbl_img_atm_map[img_i];

     img_maskdata[img_i].offset = img_mask_cnt;
     img_maskdata[img_i].cnt    = atm_maskdata[atm_i].cnt;

     img_mask_cnt = img_mask_cnt + atm_maskdata[atm_i].cnt;

     atm_off = atm_maskdata[atm_i].offset;
     img_off = img_maskdata[img_i].offset;
     for(i = 1; i <= atm_maskdata[atm_i].cnt; i++)
     {
        img_mask[img_off+i] = 
        ((unsigned int)(18 << 27))| gbl_atm_img_map[atm_mask[atm_off+i]];
     }
  }//for-img_i

////#define MASK_PRE
//#ifdef MASK_PRE
//  athread_spawn(pme_mask_setup_c_pre, &pm);
//  athread_join();
//#else
//  for(img_i = img_lo; img_i <= img_hi; img_i ++)
//  {
//     atm_i = gbl_img_atm_map[img_i];
//
//     //img_maskdata[img_i].offset = img_mask_cnt;
//     //img_maskdata[img_i].cnt    = atm_maskdata[atm_i].cnt;
//
//     //img_mask_cnt = img_mask_cnt + atm_maskdata[atm_i].cnt;
//
//     atm_off = atm_maskdata[atm_i].offset;
//     img_off = img_maskdata[img_i].offset;
//     if(myrank == 0)
//     printf("img_off = %d, cnt = %d\n", img_off,
//            atm_maskdata[atm_i].cnt);
//
//     for(i = 1; i <= atm_maskdata[atm_i].cnt; i++)
//     {
//        img_mask[img_off+i] = 
//        ((unsigned int)(18 << 27))| gbl_atm_img_map[atm_mask[atm_off+i]];
//     }
//  }//for-img_i
//#endif

  time1 = MPI_Wtime();
     
  //perf_config_t conf;
  //conf.pcrc = PCRC_ALL;
  //conf.pcr0 = PC0_CYCLE;
  //conf.pcr1 = PC1_CYCLE;
  //conf.pcr2 = PC2_N_DMA_REQ;
  ////conf.pcr2 = PC2_CNT_GLD;
  //lwpf_init(&conf);

#define PME_MASK
#ifdef PME_MASK
  
  athread_spawn(pme_mask_setup_c_para, &pm);
  athread_join();

  //int r = 0; 
  //if(myrank == 0) 
  //{
  //  lwpf_report_summary(stdout, &conf);
  //  //lwpf_report_detail(stdout, &conf);
  //}

#else

  for(icit = my_bkt_lo; icit <=  my_bkt_hi; icit++)
  {
    my = icit / cit_tbl_x_dim;
    mz = my / cit_tbl_y_dim;
    mx = icit - my * cit_tbl_x_dim;
    my = my - mz * cit_tbl_y_dim;

    cit_lo = cit_save[icit].img_lo;
    cit_hi = cit_save[icit].img_hi;

    cit_lo = max(img_lo,cit_lo);
    cit_hi = min(img_hi,cit_hi);
    //version(1)
    //for(idbkt = 1; idbkt <= ndbkt; idbkt++)
    //{
    //   dx = dbkt[idbkt * 3 + 0];
    //   dy = dbkt[idbkt * 3 + 1];
    //   dz = dbkt[idbkt * 3 + 2];
    //   cx = mx + dx;
    //   cy = my + dy;
    //   cz = mz + dz;
    //   cit_wrap_c(cx, cy, cz, 
    //               &nx, &ny, &nz, 
    //               &tx, &ty, &tz,
    //               cit_tbl_x_dim, 
    //               cit_tbl_y_dim, 
    //               cit_tbl_z_dim);//
    //   itran = (tz + 1) * 9 + (ty + 1) * 3 + tx + 1;
    //   jcit = (nz * cit_tbl_y_dim + ny) * cit_tbl_x_dim + nx;
    //   jlo = cit_save[jcit].img_lo;
    //   jhi = cit_save[jcit].img_hi;

    //   for(img_i = cit_lo; img_i <= cit_hi; img_i++)
    //   {
    //      if (idbkt == ndbkt) jhi = img_i - 1;
    //      img_off = img_maskdata[img_i].offset;

    //      for(i = img_off + 1; i <= img_off + img_maskdata[img_i].cnt; i++)
    //      {
    //         img_j = img_mask[i] &  mask27;
    //         if (img_j >= jlo &&  img_j <= jhi)
    //         {
    //            img_mask[i] = (((unsigned int )itran) << 27) |  img_j;
    //         }
    //      }
    //   }//for-img_i
    //}//for-idbkt


    //version(2)
    int dbkt_lo[64], dbkt_hi[64], itran[64];
    for(idbkt = 1; idbkt <= ndbkt; idbkt++)
    {
       dx = dbkt[idbkt * 3 + 0];
       dy = dbkt[idbkt * 3 + 1];
       dz = dbkt[idbkt * 3 + 2];
       cx = mx + dx;
       cy = my + dy;
       cz = mz + dz;
       cit_wrap_c(cx, cy, cz, 
                   &nx, &ny, &nz, 
                   &tx, &ty, &tz,
                   cit_tbl_x_dim, 
                   cit_tbl_y_dim, 
                   cit_tbl_z_dim);//
       itran[idbkt] = (tz + 1) * 9 + (ty + 1) * 3 + tx + 1;
       jcit = (nz * cit_tbl_y_dim + ny) * cit_tbl_x_dim + nx;
       dbkt_lo[idbkt] = cit_save[jcit].img_lo;
       dbkt_hi[idbkt] = cit_save[jcit].img_hi;
    }//for-idbkt

    for(img_i = cit_lo; img_i <= cit_hi; img_i++)
    {
       //if (idbkt == ndbkt) jhi = img_i - 1;
       img_off = img_maskdata[img_i].offset;

       for(i = img_off + 1; i <= img_off + img_maskdata[img_i].cnt; i++)
       {
          img_j = img_mask[i] &  mask27;
          for(idbkt = 1; idbkt <= ndbkt; idbkt++)
          {
            if(idbkt == ndbkt) dbkt_hi[idbkt] = img_i - 1;
            if (img_j >= dbkt_lo[idbkt] &&  img_j <= dbkt_hi[idbkt])
            {
               img_mask[i] = (((unsigned int )itran[idbkt]) << 27) |  img_j;
            }
          }
       }
    }//for-img_i
  }//for-icit
  #endif
  time2 = MPI_Wtime();
  //if(myrank == 0)
  //{
  //  printf("total = %.10f, %.10f, %.10f\n", 
  //        time2-time0, time1-time0, time2-time1);
  //}

}
#endif



#ifdef CPE
#define DMA_FAST
#include<slave.h>
#include<dma_macros.h>
//#define MSTEP 64
#define MSTEP 32
#define C_P 64
//#define CSTEP 64
#define CSTEP 4
#define C_M 200
#define C_D 200

#define STEP 64

#define cCACHE_HBIT 11
#define cCACHE_SBIT 4
#define cCACHE_LINESIZE (1 << cCACHE_SBIT)
#define cCACHE_LINECNT (1 << (cCACHE_HBIT - cCACHE_SBIT))
#define cCACHE_MMASK (cCACHE_LINESIZE - 1)
#define cCACHE_LMASK (cCACHE_LINECNT - 1)

#define LWPF_UNIT U(PME_SETUP_MASK)
#define LWPF_KERNELS K(ALL)  K(MAIN) K(ICIT) K(BEFORE_IMG_I) K(CACHE) K(IMG_I)  K(MASK) K(GET) K(PUT)
#include <lwpf2/lwpf2.h>
#include<param.h>

static inline void wrap_c_cpe(int x, int xlim, int *ox_, int *tx_)
{
  int ox = *ox_;
  int tx = *tx_;
  ox = x;
  tx = 0;
  if(ox < 0)
  {
     ox = ox + xlim;
     tx = -1;
  }
  if(ox >= xlim)
  {
     ox = ox - xlim;
     tx = 1;
  }
  *ox_ = ox;
  *tx_ = tx;
}

static inline void cit_wrap_c_cpe(int x, int y, int z, int *ox, int *oy, int *oz, int *tx, int *ty, int *tz, int xlim, int ylim, int zlim)
{
  wrap_c_cpe(x, xlim, ox, tx);
  wrap_c_cpe(y, ylim, oy, ty);
  wrap_c_cpe(z, zlim, oz, tz);
}

void pme_mask_setup_c_pre(pme_mask_param *pm)//unfinished
{

  dma_init();
  pme_mask_param lp;
  pe_get(pm, &lp, sizeof(pme_mask_param));
  dma_syn();

  int img_lo                      = lp.img_lo;
  int img_hi                      = lp.img_hi;
  listdata_rec_c *atm_maskdata    = lp.atm_maskdata;
  int *atm_mask                   = lp.atm_mask;
  listdata_rec_c *img_maskdata    = lp.img_maskdata;
  int *img_mask                   = lp.img_mask;
  int *gbl_atm_img_map            = lp.gbl_atm_img_map;
  int *gbl_img_atm_map            = lp.gbl_img_atm_map;
  int myrank                      = lp.myrank;

  int  idbkt, itran;
  int  img_i, img_j, jhi, jlo;
  int  atm_i, i, img_off, atm_off;
  int  img_mask_cnt;
  int  mask27 = 0x07ffffff; 
  img_mask_cnt = 0;

  int mst, med, msz, moff;
  //listdata_rec_c l_atm_maskdata[C_P];
  int l_atm_mask[C_P];
  listdata_rec_c l_img_maskdata[MSTEP];
  int l_img_mask[C_P];
  //int l_gbl_atm_img_map[C_P];
  int l_gbl_img_atm_map[MSTEP];
  unsigned int mask1827 = (unsigned int)(18 << 27);
  
  for(mst = img_lo + _MYID * MSTEP; mst < img_hi +1; mst += MSTEP * 64)
  {
    med = mst + MSTEP;
    if(med > img_hi + 1)
      med = img_hi + 1;
    msz = med - mst;

    pe_get(gbl_img_atm_map + mst, l_gbl_img_atm_map, 
            sizeof(int) * msz);
    pe_get(img_maskdata + mst, l_img_maskdata, 
            sizeof(listdata_rec_c) * msz);
    dma_syn();

    for(img_i = mst; img_i < med; img_i ++)
    {
      int moff    = img_i - mst;
      atm_i       = l_gbl_img_atm_map[moff];
      atm_off     = atm_maskdata[atm_i].offset;
      int a_cnt   = atm_maskdata[atm_i].cnt;
      img_off     = l_img_maskdata[moff].offset;
      
      if(a_cnt > 0)
      {
        pe_get(img_mask + img_off + 1, l_img_mask, sizeof(int) * a_cnt);
        pe_get(atm_mask + atm_off + 1, l_atm_mask, sizeof(int) * a_cnt);
        dma_syn();

        for(i = 1; i <= a_cnt; i++)
        {
           l_img_mask[i-1] = mask1827 | gbl_atm_img_map[l_atm_mask[i-1]];
        }

        pe_put(img_mask + img_off + 1, l_img_mask, sizeof(int) * a_cnt);
        dma_syn();
      }

    }//for-img_i
  }//for-mst

}
void pme_mask_setup_c_para(pme_mask_param *pm)
{
  //lwpf_enter(PME_SETUP_MASK);
  //lwpf_start(ALL);

  dma_init();
  pme_mask_param lp;
  pe_get(pm, &lp, sizeof(pme_mask_param));
  dma_syn();

  int img_lo                      = lp.img_lo;
  int img_hi                      = lp.img_hi;
  int my_bkt_lo                   = lp.my_bkt_lo; 
  int my_bkt_hi                   = lp.my_bkt_hi;
  int cit_tbl_x_dim               = lp.cit_tbl_x_dim;
  int cit_tbl_y_dim               = lp.cit_tbl_y_dim;
  int cit_tbl_z_dim               = lp.cit_tbl_z_dim;
  int ndbkt                       = lp.ndbkt;
  int *dbkt                       = lp.dbkt;
  listdata_rec_c *atm_maskdata    = lp.atm_maskdata;
  int *atm_mask                   = lp.atm_mask;
  listdata_rec_c *img_maskdata    = lp.img_maskdata;
  int *img_mask                   = lp.img_mask;
  int *gbl_atm_img_map            = lp.gbl_atm_img_map;
  int *gbl_img_atm_map            = lp.gbl_img_atm_map;
  cit_tbl_rec_c *cit_save         = lp.cit_save;
  int ncit_save                   = lp.ncit_save;
  int myrank                      = lp.myrank;

  int  idbkt;//, itran;
  int  dx, dy, dz, xhi, yhi;
  int  mx, my, mz, nx, ny, nz, tx, ty, tz, cx, cy, cz;
  int  icit, cit_lo, cit_hi, jcit;
  int  img_i, img_j, jhi, jlo;
  int  atm_i, i, img_off, atm_off;
  int  img_mask_cnt;
  int  mask27 = 0x07ffffff; 
  img_mask_cnt = 0;

  int mst, med, msz, moff;
  int l_gbl_img_atm_map[MSTEP];
  
  int cst, ced, csz, coff;
  cit_tbl_rec_c l_cit_save[CSTEP];
  cit_tbl_rec_c cj_cache[cCACHE_LINECNT][cCACHE_LINESIZE];
  int tag_cache[cCACHE_LINECNT];

  listdata_rec_c l_img_maskdata[STEP];
  int l_img_mask[C_M];
  int l_dbkt[C_D];
  pe_get(dbkt, l_dbkt, sizeof(int)*(ndbkt+1) * 3);
  dma_syn();
  
  for(i = 0; i < cCACHE_LINECNT; i++)
    tag_cache[i] = -1;
  //lwpf_start(MAIN);
  for(cst = my_bkt_lo + _MYID * CSTEP; cst < my_bkt_hi +1; cst += CSTEP * 64)
  {
    ced = cst + CSTEP;
    if(ced > my_bkt_hi + 1)
      ced = my_bkt_hi + 1;
    csz = ced - cst;
    
    pe_get(cit_save + cst, l_cit_save, sizeof(cit_tbl_rec_c) * csz);
    dma_syn();

    //lwpf_start(ICIT);
    for(icit = cst; icit < ced; icit++)
    {
      coff = icit - cst;
      my = icit / cit_tbl_x_dim;
      mz = my / cit_tbl_y_dim;
      mx = icit - my * cit_tbl_x_dim;
      my = my - mz * cit_tbl_y_dim;

      cit_lo = l_cit_save[coff].img_lo;
      cit_hi = l_cit_save[coff].img_hi;

      cit_lo = max(img_lo,cit_lo);
      cit_hi = min(img_hi,cit_hi);

      int dbkt_lo[64], dbkt_hi[64], itran[64];
      //lwpf_start(BEFORE_IMG_I);
      for(idbkt = 1; idbkt <= ndbkt; idbkt++)
      {
          dx = l_dbkt[idbkt * 3 + 0];
          dy = l_dbkt[idbkt * 3 + 1];
          dz = l_dbkt[idbkt * 3 + 2];
          cx = mx + dx;
          cy = my + dy;
          cz = mz + dz;
          cit_wrap_c_cpe(cx, cy, cz, 
                     &nx, &ny, &nz, 
                     &tx, &ty, &tz,
                     cit_tbl_x_dim, 
                     cit_tbl_y_dim, 
                     cit_tbl_z_dim);//
          //itran = (tz + 1) * 9 + (ty + 1) * 3 + tx + 1;
          itran[idbkt] = (tz + 1) * 9 + (ty + 1) * 3 + tx + 1;
          jcit = (nz * cit_tbl_y_dim + ny) * cit_tbl_x_dim + nx;

          if (jcit < 0 || jcit >= ncit_save) printf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n", icit, idbkt, jcit, dx, dy, dz, nx, ny, nz);
          //lwpf_start(CACHE);
          if (tag_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK] 
              != jcit >> cCACHE_SBIT)
          {
            pe_get(cit_save + (jcit & ~cCACHE_MMASK), 
                  cj_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK], 
                  sizeof(cit_tbl_rec_c) * cCACHE_LINESIZE);
            dma_syn();
            tag_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK] 
                    = jcit >> cCACHE_SBIT;
          }
          cit_tbl_rec_c pcj = 
          cj_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK][jcit & cCACHE_MMASK];
          //lwpf_stop(CACHE);

          //jlo = pcj.img_lo;
          //jhi = pcj.img_hi;
          dbkt_lo[idbkt] = pcj.img_lo;
          dbkt_hi[idbkt] = pcj.img_hi;
      }//for_idbkt

      //if(_MYID == 46)
      //  if(cit_hi - cit_lo + 1 > 64)
      //  printf("cit_hi - cit_lo = %d\n", cit_hi-cit_lo+1);

    if(cit_hi >= cit_lo)
    {
      pe_get(img_maskdata + cit_lo, l_img_maskdata, 
             sizeof(listdata_rec_c)*(cit_hi - cit_lo + 1));
      dma_syn();
    }

      //lwpf_stop(BEFORE_IMG_I);

      //////////no-partition/////////////////////////
      //lwpf_start(IMG_I);
      for(img_i = cit_lo; img_i <= cit_hi; img_i++)
      {
        int ioff = img_i - cit_lo;
        //if (idbkt == ndbkt) jhi = img_i - 1;
        img_off      = l_img_maskdata[ioff].offset;
        int img_cnt  = l_img_maskdata[ioff].cnt;
       
       ////lwpf_start(GET);
       if(img_cnt > 0)
       {
        pe_get(img_mask + img_off + 1, l_img_mask, 
              sizeof(int) * img_cnt);
        dma_syn();
       }
       //lwpf_stop(GET);

       //lwpf_start(MASK);
       int put_cnt = 0;
       for(i = img_off + 1; i <= img_off + img_cnt; i++)
       {
         int offset = i - img_off -1;
         img_j = l_img_mask[offset] &  mask27;
         for(idbkt = 1; idbkt <= ndbkt; idbkt++)
         {
           if(idbkt == ndbkt) dbkt_hi[idbkt] = img_i - 1;
           if (img_j >= dbkt_lo[idbkt] &&  img_j <= dbkt_hi[idbkt])
           {
              l_img_mask[offset] = 
               (((unsigned int )itran[idbkt]) << 27) |  img_j;
              put_cnt++;
           }
         }
       }//for-i
       ////lwpf_stop(MASK);

       ////lwpf_start(PUT);
       if(img_cnt > 0 && put_cnt > 0)
       {
         pe_put(img_mask + img_off + 1, l_img_mask, 
               sizeof(int) * img_cnt);
         dma_syn();
       }
       //lwpf_stop(PUT);


//      /////////////subpartition//////////
//       if(img_cnt > 0)
//       {
//          int base_pos = img_off + 1;
//          int icnt, icsz, icst, iced;
//          for(icnt = base_pos; icnt < base_pos + img_cnt; icnt += C_P)
//          {
//            icst = icnt;
//            iced = icst + C_P;
//            if(iced > base_pos + img_cnt)
//              iced = base_pos + img_cnt;
//            icsz = iced - icst;
//
//            pe_get(img_mask + icnt, l_img_mask, 
//                  sizeof(int) * icsz);
//            dma_syn();
//            
//            int put_cnt = 0;
//            for(i = icnt; i < icnt + icsz; i++)
//            //for(i = img_off + 1; i <= img_off + img_cnt; i++)
//            {
//              int offset = i - icnt;
//              img_j = l_img_mask[offset] &  mask27;
//              for(idbkt = 1; idbkt <= ndbkt; idbkt++)
//              {
//                if(idbkt == ndbkt) dbkt_hi[idbkt] = img_i - 1;
//                if (img_j >= dbkt_lo[idbkt] &&  img_j <= dbkt_hi[idbkt])
//                {
//                   l_img_mask[offset] = 
//                    (((unsigned int )itran[idbkt]) << 27) |  img_j;
//                   put_cnt++;
//                }
//              }
//            }//for-i
//            
//            if(put_cnt > 0)
//            {
//              pe_put(img_mask + icnt, l_img_mask, 
//                    sizeof(int) * icsz);
//              dma_syn();
//            }
//
//
//          }//for-icnt
//      }//if-img_cnt




      }//for-img_i

      //lwpf_stop(IMG_I);
    }//for-icit

    //lwpf_stop(ICIT);

  }//for-cst
  //lwpf_stop(MAIN);
  //lwpf_stop(ALL);
  //lwpf_exit(PME_SETUP_MASK);

}

void pme_mask_setup_c_para_test(pme_mask_param *pm)
{
  //lwpf_enter(PME_SETUP_MASK);
  lwpf_start(ALL);

  dma_init();
  pme_mask_param lp;
  pe_get(pm, &lp, sizeof(pme_mask_param));
  dma_syn();

  int img_lo                      = lp.img_lo;
  int img_hi                      = lp.img_hi;
  int my_bkt_lo                   = lp.my_bkt_lo; 
  int my_bkt_hi                   = lp.my_bkt_hi;
  int cit_tbl_x_dim               = lp.cit_tbl_x_dim;
  int cit_tbl_y_dim               = lp.cit_tbl_y_dim;
  int cit_tbl_z_dim               = lp.cit_tbl_z_dim;
  int ndbkt                       = lp.ndbkt;
  int *dbkt                       = lp.dbkt;
  listdata_rec_c *atm_maskdata    = lp.atm_maskdata;
  int *atm_mask                   = lp.atm_mask;
  listdata_rec_c *img_maskdata    = lp.img_maskdata;
  int *img_mask                   = lp.img_mask;
  int *gbl_atm_img_map            = lp.gbl_atm_img_map;
  int *gbl_img_atm_map            = lp.gbl_img_atm_map;
  cit_tbl_rec_c *cit_save         = lp.cit_save;
  int ncit_save                   = lp.ncit_save;
  int myrank                      = lp.myrank;

  int  idbkt, itran;
  int  dx, dy, dz, xhi, yhi;
  int  mx, my, mz, nx, ny, nz, tx, ty, tz, cx, cy, cz;
  int  icit, cit_lo, cit_hi, jcit;
  int  img_i, img_j, jhi, jlo;
  int  atm_i, i, img_off, atm_off;
  int  img_mask_cnt;
  int  mask27 = 0x07ffffff; 
  img_mask_cnt = 0;

  int mst, med, msz, moff;
  int l_gbl_img_atm_map[MSTEP];
  
  int cst, ced, csz, coff;
  cit_tbl_rec_c l_cit_save[CSTEP];
  cit_tbl_rec_c cj_cache[cCACHE_LINECNT][cCACHE_LINESIZE];
  int tag_cache[cCACHE_LINECNT];

  listdata_rec_c l_img_maskdata[CSTEP];
  int l_img_mask[C_M];
  int l_dbkt[C_D];
  pe_get(dbkt, l_dbkt, sizeof(int)*(ndbkt+1) * 3);
  dma_syn();
  
  for(i = 0; i < cCACHE_LINECNT; i++)
    tag_cache[i] = -1;
  lwpf_start(MAIN);
  for(cst = my_bkt_lo + _MYID * CSTEP; cst < my_bkt_hi +1; cst += CSTEP * 64)
  {
    ced = cst + CSTEP;
    if(ced > my_bkt_hi + 1)
      ced = my_bkt_hi + 1;
    csz = ced - cst;
    
    pe_get(cit_save + cst, l_cit_save, sizeof(cit_tbl_rec_c) * csz);
    dma_syn();

    lwpf_start(ICIT);
    for(icit = cst; icit < ced; icit++)
    {
      coff = icit - cst;
      my = icit / cit_tbl_x_dim;
      mz = my / cit_tbl_y_dim;
      mx = icit - my * cit_tbl_x_dim;
      my = my - mz * cit_tbl_y_dim;

      cit_lo = l_cit_save[coff].img_lo;
      cit_hi = l_cit_save[coff].img_hi;

      cit_lo = max(img_lo,cit_lo);
      cit_hi = min(img_hi,cit_hi);

      for(idbkt = 1; idbkt <= ndbkt; idbkt++)
      {
          lwpf_start(BEFORE_IMG_I);
          dx = l_dbkt[idbkt * 3 + 0];
          dy = l_dbkt[idbkt * 3 + 1];
          dz = l_dbkt[idbkt * 3 + 2];
          cx = mx + dx;
          cy = my + dy;
          cz = mz + dz;
          cit_wrap_c_cpe(cx, cy, cz, 
                     &nx, &ny, &nz, 
                     &tx, &ty, &tz,
                     cit_tbl_x_dim, 
                     cit_tbl_y_dim, 
                     cit_tbl_z_dim);//
          itran = (tz + 1) * 9 + (ty + 1) * 3 + tx + 1;
          jcit = (nz * cit_tbl_y_dim + ny) * cit_tbl_x_dim + nx;

          lwpf_start(CACHE);
          if (tag_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK] 
              != jcit >> cCACHE_SBIT)
          {
            pe_get(cit_save + (jcit & ~cCACHE_MMASK), 
                  cj_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK], 
                  sizeof(cit_tbl_rec_c) * cCACHE_LINESIZE);
            dma_syn();
            tag_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK] 
                    = jcit >> cCACHE_SBIT;
          }
          cit_tbl_rec_c pcj = 
          cj_cache[(jcit >> cCACHE_SBIT) & cCACHE_LMASK][jcit & cCACHE_MMASK];
          lwpf_stop(CACHE);

          jlo = pcj.img_lo;
          jhi = pcj.img_hi;

         pe_get(img_maskdata + cit_lo, l_img_maskdata, 
                sizeof(listdata_rec_c)*(cit_hi - cit_lo + 1));
         dma_syn();

         lwpf_stop(BEFORE_IMG_I);

         lwpf_start(IMG_I);
         for(img_i = cit_lo; img_i <= cit_hi; img_i++)
         {
           int ioff = img_i - cit_lo;
           if (idbkt == ndbkt) jhi = img_i - 1;
           img_off      = l_img_maskdata[ioff].offset;
           int img_cnt  = l_img_maskdata[ioff].cnt;
          
          lwpf_start(GET);
          if(img_cnt > 0)
          {
           pe_get(img_mask + img_off + 1, l_img_mask, 
                 sizeof(int) * img_cnt);
           dma_syn();
          }
          lwpf_stop(GET);

          lwpf_start(MASK);
          int put_cnt = 0;
          for(i = img_off + 1; i <= img_off + img_cnt; i++)
          {
            int offset = i - img_off -1;
            img_j = l_img_mask[offset] &  mask27;
            if (img_j >= jlo &&  img_j <= jhi)
            {
               l_img_mask[offset] = 
                (((unsigned int )itran) << 27) |  img_j;
               put_cnt++;
            }
          }
          lwpf_stop(MASK);

          lwpf_start(PUT);
          if(img_cnt > 0 && put_cnt > 0)
          {
            pe_put(img_mask + img_off + 1, l_img_mask, 
                  sizeof(int) * img_cnt);
            dma_syn();
          }
          lwpf_stop(PUT);
         }//for-img_i
          lwpf_stop(IMG_I);
      }//for-idbkt
    }//for-icit
    lwpf_stop(ICIT);
  }//for-cst
  lwpf_stop(MAIN);
  lwpf_stop(ALL);
  lwpf_exit(PME_SETUP_MASK);

}
#endif
