#include <math.h>
#include <simd.h>
double derfcfun(double x){
  double absx, c, p, q, nonexperfc, erf, erfc;
  if (x < 0)
    absx = -x;
  else
    absx = x;
  if (x > 26) {
    erfc = 0;
  } else if (x < -5.5) {
    erfc = 2.0e0;
  } else if (absx <= 0.5) {
    c = x * x;
    p=((-0.356098437018154e-1*c+0.699638348861914e1)*c + 0.219792616182942e2)*c+0.242667955230532e3;
    q=((c+0.150827976304078e2)*c+0.911649054045149e2)*c + 0.215058875869861e3;
    erf = x*p/q;
    erfc = 1.0-erf;
  } else if (absx < 4.0) {
    c=absx;
    p=((((((-0.136864857382717e-6*c+0.564195517478974e0)*c+0.721175825088309e1)*c+0.431622272220567e2)*c+
         0.152989285046940e3)*c+0.339320816734344e3)*c+0.451918953711873e3)*c+0.300459261020162e3;
    q=((((((c+0.127827273196294e2)*c+0.770001529352295e2)*c + 0.277585444743988e3)*c+0.638980264465631e3)*c+ 
        0.931354094850610e3)*c+0.790950925327898e3)*c + 0.300459260956983e3;
    if ( x > 0.e0 ) {
      nonexperfc = p/q;
    } else {
      nonexperfc = 2.e0*exp(x*x) - p/q;
    }
    erfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.e0) erfc = 2.0 - erfc;

  } else {
    c=1.0/(x*x);
    p=(((0.223192459734185e-1*c+0.278661308609648e0)*c+0.226956593539687e0)*c+0.494730910623251e-1)*c+
      0.299610707703542e-2;
    q=(((c+0.198733201817135e1)*c+0.105167510706793e1)*c+0.191308926107830e0)*c+0.106209230528468e-1;
    c=(-c*p/q + 0.564189583547756e0)/absx;
    if( x > 0.e0 ) {
      nonexperfc = c;
    } else {
      nonexperfc = 2.e0*exp(x*x) - c;
    }
    erfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.e0) erfc = 2.0- erfc;
  }
  return erfc;
}

double erfc_wiki(double x){
  if (x >= 0) {
    double p = 0.3275911;
    double t = 1 / (1 + p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = exp(-x*x);
    return polypart * exppart;
  } else {
    double p = 0.3275911;
    double t = 1 / (1 - p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = exp(-x*x);
    return 2.0 + polypart * exppart;
  }
}

double fast_exp(double x){
  double ln2inv = 1.4426950408889634;
  double ln2 = 0.6931471805599453;
  double a1 = -0.9999999995;
  double a2 = 0.4999999206;
  double a3 = -0.1666653019;
  double a4 = 0.0416573475;
  double a5 = -0.0083013598;
  double a6 = 0.0013298820;
  double a7 = -0.0001413161;
  //x = a * ln2 - b                                                                                                                                                                                                                                                                       
  double a = ceil(x * ln2inv);
  double b = a * ln2 - x;
  //printf("%e %e %d\n", a, b, (long)a);
  double polypart = (((((((a7 * b + a6) * b + a5) * b + a4) * b + a3) * b) + a2) * b + a1) * b + 1.0;
  double e2part;
  *(long*)&e2part = (((long)a + 0x3ffUL) << 52L);
  return e2part * polypart;
}

#define RSQRT_PI 0.56418958354775628
//-(2*exp(-x*x)) * RSQRT_PI;
//#ifdef MPE
void calc_sw(double x, double *sw, double *d_sw_dx){
  if (x >= 0) {
    double p = 0.3275911;
    double t = 1 / (1 + p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = fast_exp(-x*x);
    *d_sw_dx = -(2 * exppart) * RSQRT_PI;
    *sw = polypart * exppart;
  } else {
    double p = 0.3275911;
    double t = 1 / (1 - p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = fast_exp(-x*x);
    *d_sw_dx = -(2 * exppart) * RSQRT_PI;
    *sw = 2.0 - polypart * exppart;
  }  
}

void calc_sw_(double *x_ptr, double *sw, double *d_sw_dx){
  double x = *x_ptr;
  if (x >= 0) {
    double p = 0.3275911;
    double t = 1 / (1 + p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = fast_exp(-x*x);
    *d_sw_dx = -(2 * exppart) * RSQRT_PI;
    *sw = polypart * exppart;
  } else {
    double p = 0.3275911;
    double t = 1 / (1 - p * x);
    double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
    double polypart = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
    double exppart = fast_exp(-x*x);
    *d_sw_dx = -(2 * exppart) * RSQRT_PI;
    *sw = 2.0 - polypart * exppart;
  }  
}

//#endif
#ifdef CPE
#include <slave.h>
extern __thread_local int debug;

void calc_sw_v4(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx){
  double *sw4 = (double*)sw;
  double *dsw4 = (double*)d_sw_dx;
  calc_sw(simd_vextf0(x), sw4 + 0, dsw4 + 0);
  calc_sw(simd_vextf1(x), sw4 + 1, dsw4 + 1);
  calc_sw(simd_vextf2(x), sw4 + 2, dsw4 + 2);
  calc_sw(simd_vextf3(x), sw4 + 3, dsw4 + 3);
}

doublev4 simd_fast_exp(doublev4 x){
  doublev4 ln2inv = 1.4426950408889634;
  doublev4 ln2 = 0.6931471805599453;
  doublev4 a1 = -0.9999999995;
  doublev4 a2 = 0.4999999206;
  doublev4 a3 = -0.1666653019;
  doublev4 a4 = 0.0416573475;
  doublev4 a5 = -0.0083013598;
  doublev4 a6 = 0.0013298820;
  doublev4 a7 = -0.0001413161;
  
  doublev4 xln2inv = x * ln2inv;
  doublev4 a_f;
  int256 a_i;
  float t0, t1, t2, t3;
  asm ("vextf %2, 0, %3\n\t"
       "vextf %2, 1, %4\n\t"
       "vextf %2, 2, %5\n\t"
       "vextf %2, 3, %6\n\t"
       "fcvtdlr %3, 7, %3\n\t"
       "fcvtdlr %4, 7, %4\n\t"
       "fcvtdlr %5, 7, %5\n\t"
       "fcvtdlr %6, 7, %6\n\t"
       "vinsf %3, $31, 0, %1\n\t"
       "vinsf %4,  %1, 1, %1\n\t"
       "vinsf %5,  %1, 2, %1\n\t"
       "vinsf %6,  %1, 3, %1\n\t"
       "fcvtld %3, %3\n\t"
       "fcvtld %4, %4\n\t"
       "fcvtld %5, %5\n\t"
       "fcvtld %6, %6\n\t"
       "vinsf %3, %0, 0, %0\n\t"
       "vinsf %4, %0, 1, %0\n\t"
       "vinsf %5, %0, 2, %0\n\t"
       "vinsf %6, %0, 3, %0\n\t"
       : "=&r"(a_f), "=&r"(a_i)
       : "r"(xln2inv), "r"(t0), "r"(t1), "r"(t2), "r"(t3));
  doublev4 b = a_f * ln2 - x;
  //doublev4 polypart = (((((((a7 * b + a6) * b + a5) * b + a4) * b + a3) * b) + a2) * b + a1) * b + 1.0;
  doublev4 p = a7 * b + a6;
  p = p * b + a5;
  p = p * b + a4;
  p = p * b + a3;
  p = p * b + a2;
  p = p * b + a1;
  p = p * b + 1.0;
  doublev4 e2part;
  int256 e2parti;
  int shfw_mask = 0x67452301;
  asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
       "vaddl %0, $31, %1\n\t"
       : "=&r"(e2part) : "r"(a_i), "r"(shfw_mask));
  return e2part * p;
}

doublev4 simd_very_fast_exp(doublev4 x){
  doublev4 ln2inv = 1.4426950408889634;
  doublev4 ln2 = 0.6931471805599453;
  /* doublev4 a1 = -0.9999999995; */
  /* doublev4 a2 = 0.4999999206; */
  /* doublev4 a3 = -0.1666653019; */
  /* doublev4 a4 = 0.0416573475; */
  /* doublev4 a5 = -0.0083013598; */
  /* doublev4 a6 = 0.0013298820; */
  /* doublev4 a7 = -0.0001413161; */
  doublev4 a1=-0.9998684;
  doublev4 a2=0.4982926;
  doublev4 a3=-0.1595332;
  doublev4 a4=0.0293641;
  doublev4 xln2inv = x * ln2inv;
  doublev4 a_f;
  int256 a_i;
  float t0, t1, t2, t3;
  asm ("vextf %2, 0, %3\n\t"
       "vextf %2, 1, %4\n\t"
       "vextf %2, 2, %5\n\t"
       "vextf %2, 3, %6\n\t"
       "fcvtdlr %3, 7, %3\n\t"
       "fcvtdlr %4, 7, %4\n\t"
       "fcvtdlr %5, 7, %5\n\t"
       "fcvtdlr %6, 7, %6\n\t"
       "vinsf %3, $31, 0, %1\n\t"
       "vinsf %4,  %1, 1, %1\n\t"
       "vinsf %5,  %1, 2, %1\n\t"
       "vinsf %6,  %1, 3, %1\n\t"
       "fcvtld %3, %3\n\t"
       "fcvtld %4, %4\n\t"
       "fcvtld %5, %5\n\t"
       "fcvtld %6, %6\n\t"
       "vinsf %3, %0, 0, %0\n\t"
       "vinsf %4, %0, 1, %0\n\t"
       "vinsf %5, %0, 2, %0\n\t"
       "vinsf %6, %0, 3, %0\n\t"
       : "=&r"(a_f), "=&r"(a_i)
       : "r"(xln2inv), "r"(t0), "r"(t1), "r"(t2), "r"(t3));
  doublev4 b = a_f * ln2 - x;
  //doublev4 polypart = (((((((a7 * b + a6) * b + a5) * b + a4) * b + a3) * b) + a2) * b + a1) * b + 1.0;
  /* doublev4 p = a7 * b + a6; */
  /* p = p * b + a5; */
  /* p = p * b + a4; */
  /* p = p * b + a3; */
  /* p = p * b + a2; */
  /* p = p * b + a1; */
  /* p = p * b + 1.0; */
  doublev4 p = a4 * b + a3;
  p = p * b + a2;
  p = p * b + a1;
  p = p * b + 1.0;
  doublev4 e2part;
  int256 e2parti;
  int shfw_mask = 0x67452301;
  asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
       "vaddl %0, $31, %1\n\t"
       : "=&r"(e2part) : "r"(a_i), "r"(shfw_mask));
  /* if (debug){ */
  /*   puts("detached"); */
  /*   simd_print_doublev4(e2part * p); */
  /*   simd_print_doublev4(e2part); */
  /*   simd_print_doublev4(p); */
  /* } */
  return e2part * p;
}

void simd_calc_sw(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx){
  doublev4 p = 0.3275911;
  doublev4 absx;
  asm ("vcpys $31, %1, %0\n\t"
       : "=r"(absx) : "r"(x));
  doublev4 c1 = 1.0;
  doublev4 t = c1 / (c1 + p * absx);
  doublev4 a1 = 0.254829592;
  doublev4 a2 = -0.284496736;
  doublev4 a3 = 1.421413741;
  doublev4 a4 = -1.453152027;
  doublev4 a5 = 1.061405429;

  //doublev4 p = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
  doublev4 pp = a5 * t + a4;
  pp = pp * t + a3;
  pp = pp * t + a2;
  pp = pp * t + a1;
  pp = pp * t;
  doublev4 exppart = simd_very_fast_exp(-x*x);
  *d_sw_dx = -(2 * exppart) * RSQRT_PI;
  doublev4 erfabs = c1 - pp * exppart;
  *sw = c1 - simd_vcpys(x, erfabs);
}


doublev4 simd_very_fast_exp_trans(doublev4 negxsq){
  doublev4 ln2inv = 1.4426950408889634;
  doublev4 ln2 = 0.6931471805599453;
  doublev4 exp_a1=-0.9998684;
  doublev4 exp_a2=0.4982926;
  doublev4 exp_a3=-0.1595332;
  doublev4 exp_a4=0.0293641;
  doublev4 xln2inv = negxsq * ln2inv;
  doublev4 nln2_f;
  int256 nln2_i;
  float t0, t1, t2, t3;
  asm ("vextf %2, 0, %3\n\t"
       "vextf %2, 1, %4\n\t"
       "vextf %2, 2, %5\n\t"
       "vextf %2, 3, %6\n\t"
       "fcvtdlr %3, 7, %3\n\t"
       "fcvtdlr %4, 7, %4\n\t"
       "fcvtdlr %5, 7, %5\n\t"
       "fcvtdlr %6, 7, %6\n\t"
       "vinsf %3, $31, 0, %1\n\t"
       "vinsf %4,  %1, 1, %1\n\t"
       "vinsf %5,  %1, 2, %1\n\t"
       "vinsf %6,  %1, 3, %1\n\t"
       "fcvtld %3, %3\n\t"
       "fcvtld %4, %4\n\t"
       "fcvtld %5, %5\n\t"
       "fcvtld %6, %6\n\t"
       "vinsf %3, %0, 0, %0\n\t"
       "vinsf %4, %0, 1, %0\n\t"
       "vinsf %5, %0, 2, %0\n\t"
       "vinsf %6, %0, 3, %0\n\t"
       : "=&r"(nln2_f), "=&r"(nln2_i)
       : "r"(xln2inv), "r"(t0), "r"(t1), "r"(t2), "r"(t3));
  doublev4 exp_b = nln2_f * ln2 - negxsq;
  doublev4 exp_poly = exp_a4 * exp_b + exp_a3;
  exp_poly = exp_poly * exp_b + exp_a2;
  exp_poly = exp_poly * exp_b + exp_a1;
  exp_poly = exp_poly * exp_b + 1.0;
  doublev4 exp_p2n;
  int exp_shfw_mask = 0x67452301;
  asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
       "vaddl %0, $31, %1\n\t"
       : "=&r"(exp_p2n) : "r"(nln2_i), "r"(exp_shfw_mask));
  return exp_p2n * exp_poly;
}

__thread_local double sw_tbl[] = {
  /* 0*/1.0, // 1.0
  /* 1*/2.0, // 2.0
  /* 2*/0.3275911, //erf_p
  /* 3*/0.254829592, //erf_a1
  /* 4*/-0.284496736, //erf_a2
  /* 5*/1.421413741, //erf_a3
  /* 6*/-1.453152027, //erf_a4
  /* 7*/1.061405429, //erf_a5
  /* 8*/1.4426950408889634, //ln2inv
  /* 9*/0.6931471805599453, //ln2
  /*10*/-0.9998684, //exp_a1
  /*11*/0.4982926, //exp_a2
  /*12*/-0.1595332, //exp_a3
  /*13*/0.0293641 //exp_a4
};
void simd_calc_sw_inline_exp(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx) {
  /* doublev4 c1 = 1.0; */
  /* doublev4 c2 = 2.0; */
  /* doublev4 erf_p = 0.3275911; */
  /* doublev4 erf_a1 = 0.254829592; */
  /* doublev4 erf_a2 = -0.284496736; */
  /* doublev4 erf_a3 = 1.421413741; */
  /* doublev4 erf_a4 = -1.453152027; */
  /* doublev4 erf_a5 = 1.061405429; */
  /* doublev4 ln2inv = 1.4426950408889634; */
  /* doublev4 ln2 = 0.6931471805599453; */
  /* doublev4 exp_a1=-0.9998684; */
  /* doublev4 exp_a2=0.4982926; */
  /* doublev4 exp_a3=-0.1595332; */
  /* doublev4 exp_a4=0.0293641; */
  doublev4 negxsq = -x * x;
  doublev4 absx;
  asm ("vcpys $31, %1, %0\n\t"
       : "=r"(absx) : "r"(x));
  doublev4 xln2inv = negxsq * sw_tbl[8];

  doublev4 nln2_f;
  int256 nln2_i;
  float t0, t1, t2, t3;
  asm ("vextf %2, 0, %3\n\t"
       "vextf %2, 1, %4\n\t"
       "vextf %2, 2, %5\n\t"
       "vextf %2, 3, %6\n\t"
       "fcvtdlr %3, 7, %3\n\t"
       "fcvtdlr %4, 7, %4\n\t"
       "fcvtdlr %5, 7, %5\n\t"
       "fcvtdlr %6, 7, %6\n\t"
       "vinsf %3, $31, 0, %1\n\t"
       "vinsf %4,  %1, 1, %1\n\t"
       "vinsf %5,  %1, 2, %1\n\t"
       "vinsf %6,  %1, 3, %1\n\t"
       "fcvtld %3, %3\n\t"
       "fcvtld %4, %4\n\t"
       "fcvtld %5, %5\n\t"
       "fcvtld %6, %6\n\t"
       "vinsf %3, %0, 0, %0\n\t"
       "vinsf %4, %0, 1, %0\n\t"
       "vinsf %5, %0, 2, %0\n\t"
       "vinsf %6, %0, 3, %0\n\t"
       : "=&r"(nln2_f), "=&r"(nln2_i)
       : "r"(xln2inv), "r"(t0), "r"(t1), "r"(t2), "r"(t3));
  doublev4 exp_p2n;
  int exp_shfw_mask = 0x67452301;
  asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
       "vaddl %0, $31, %1\n\t"
       : "=&r"(exp_p2n) : "r"(nln2_i), "r"(exp_shfw_mask));

  doublev4 erf_t = sw_tbl[0] / (sw_tbl[0] + sw_tbl[3] * absx);
  doublev4 exp_b = nln2_f * sw_tbl[9] - negxsq;
  doublev4 erf_poly = sw_tbl[7] * erf_t + sw_tbl[6];
  doublev4 exp_poly = sw_tbl[13] * exp_b + sw_tbl[12];
  erf_poly = erf_poly * erf_t + sw_tbl[5];
  exp_poly = exp_poly * exp_b + sw_tbl[11];
  erf_poly = erf_poly * erf_t + sw_tbl[4];
  exp_poly = exp_poly * exp_b + sw_tbl[10];
  erf_poly = erf_poly * erf_t + sw_tbl[3];
  exp_poly = exp_poly * exp_b + sw_tbl[0];
  erf_poly = erf_poly * erf_t;

  doublev4 expnxsq = exp_p2n * exp_poly;
  *d_sw_dx = -(sw_tbl[1] * expnxsq) * RSQRT_PI;
  doublev4 erfabs = sw_tbl[0] - erf_poly * expnxsq;
  *sw = sw_tbl[0] - simd_vcpys(x, erfabs);

}

void simd_fast_calc_sw(doublev4 x, doublev4 *sw, doublev4 *d_sw_dx) {
  doublev4 c1 = 1.0;
  doublev4 c2 = 2.0;
  doublev4 erf_p = 0.47047;
  doublev4 erf_a1 = 0.3480242;
  doublev4 erf_a2 = -0.0958798;
  doublev4 erf_a3 = 0.7478556;
  doublev4 ln2inv = 1.4426950408889634;
  doublev4 ln2 = 0.6931471805599453;
  doublev4 exp_a1=-0.9998684;
  doublev4 exp_a2=0.4982926;
  doublev4 exp_a3=-0.1595332;
  doublev4 exp_a4=0.0293641;
  doublev4 negxsq = -x * x;
  doublev4 absx;
  asm ("vcpys $31, %1, %0\n\t"
       : "=r"(absx) : "r"(x));
  doublev4 xln2inv = negxsq * ln2inv;

  doublev4 nln2_f;
  int256 nln2_i;
  float t0, t1, t2, t3;
  asm ("vextf %2, 0, %3\n\t"
       "vextf %2, 1, %4\n\t"
       "vextf %2, 2, %5\n\t"
       "vextf %2, 3, %6\n\t"
       "fcvtdlr %3, 7, %3\n\t"
       "fcvtdlr %4, 7, %4\n\t"
       "fcvtdlr %5, 7, %5\n\t"
       "fcvtdlr %6, 7, %6\n\t"
       "vinsf %3, $31, 0, %1\n\t"
       "vinsf %4,  %1, 1, %1\n\t"
       "vinsf %5,  %1, 2, %1\n\t"
       "vinsf %6,  %1, 3, %1\n\t"
       "fcvtld %3, %3\n\t"
       "fcvtld %4, %4\n\t"
       "fcvtld %5, %5\n\t"
       "fcvtld %6, %6\n\t"
       "vinsf %3, %0, 0, %0\n\t"
       "vinsf %4, %0, 1, %0\n\t"
       "vinsf %5, %0, 2, %0\n\t"
       "vinsf %6, %0, 3, %0\n\t"
       : "=&r"(nln2_f), "=&r"(nln2_i)
       : "r"(xln2inv), "r"(t0), "r"(t1), "r"(t2), "r"(t3));
  doublev4 exp_p2n;
  int exp_shfw_mask = 0x67452301;
  asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
       "vaddl %0, $31, %1\n\t"
       : "=&r"(exp_p2n) : "r"(nln2_i), "r"(exp_shfw_mask));

  doublev4 erf_t_divisor = c1 + erf_p * absx;
  doublev4 erf_t;// = c1 / (c1 + erf_p * absx);
  asm ("vdivs %1, %2, %0" : "=r"(erf_t): "r"(c1), "r"(erf_t_divisor));
  doublev4 exp_b = nln2_f * ln2 - negxsq;
  doublev4 erf_poly = erf_a3 * erf_t + erf_a2;
  erf_poly = erf_poly * erf_t + erf_a1;
  erf_poly = erf_poly * erf_t;

  doublev4 exp_poly = exp_a4 * exp_b + exp_a3;
  exp_poly = exp_poly * exp_b + exp_a2;
  exp_poly = exp_poly * exp_b + exp_a1;
  exp_poly = exp_poly * exp_b + c1;


  doublev4 expnxsq = exp_p2n * exp_poly;
  *d_sw_dx = -(c2 * expnxsq) * RSQRT_PI;
  doublev4 erfabs = c1 - erf_poly * expnxsq;
  *sw = c1 - simd_vcpys(x, erfabs);
}

#endif
