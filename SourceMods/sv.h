#define SVF 2
typedef doublev4 sv_t[SVF];

#define sv_sv_bop(r, op, a, b)                  \
  r[0] = a[0] op b[0];                          \
  r[1] = a[1] op b[1];                          \

#define sv_v_bop(r, op, a, b)                   \
  r[0] = a[0] op b;                             \
  r[1] = a[1] op b;                             \
  
#define v_sv_bop(r, op, a, b)                   \
  r[0] = a op b[0];                             \
  r[1] = a op b[1];                             \

#define sv_sv_sv_fma(r, a, b, c)                \
  r[0] = a[0] * b[0] + c[0];                    \
  r[1] = a[1] * b[1] + c[1];                    \

#define v_sv_sv_fma(r, a, b, c)                 \
  r[0] = a * b[0] + c[0];                       \
  r[1] = a * b[1] + c[1];                       \

#define sv_sv_v_fma(r, a, b, c)                 \
  r[0] = a[0] * b[0] + c;                       \
  r[1] = a[1] * b[1] + c;                       \

#define sv_sv_add(r, a, b) sv_sv_bop(r, + , a, b)
#define sv_sv_sub(r, a, b) sv_sv_bop(r, - , a, b)
#define sv_sv_div(r, a, b) sv_sv_bop(r, / , a, b)
#define sv_sv_mul(r, a, b) sv_sv_bop(r, * , a, b)

#define v_sv_add(r, a, b) v_sv_bop(r, + , a, b)
#define v_sv_sub(r, a, b) v_sv_bop(r, - , a, b)
#define v_sv_div(r, a, b) v_sv_bop(r, / , a, b)
#define v_sv_mul(r, a, b) v_sv_bop(r, * , a, b)

#define sv_v_add(r, a, b) sv_v_bop(r, + , a, b)
#define sv_v_sub(r, a, b) sv_v_bop(r, - , a, b)
#define sv_v_div(r, a, b) sv_v_bop(r, / , a, b)
#define sv_v_mul(r, a, b) sv_v_bop(r, * , a, b)

#define sv_sqrt(r, a)                           \
  r[0] = simd_vsqrtd(a[0]);                     \
  r[1] = simd_vsqrtd(a[1]);                     \
