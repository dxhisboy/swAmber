#define CAT(x, y) __CAT__(x, y)
#define __CAT__(x, y) x ## y
#define STR(x) __STR__(x)
#define __STR__(x) #x

/* #define PREFIX pairs_calc_a2i */
/* #include <pairs_calc_gen.h> */
/* #undef PREFIX */

/* #define PREFIX pairs_calc_a2b */
/* #include <pairs_calc_gen.h> */
/* #undef PREFIX */

#define PREFIX pairs_calc_v2b
#include <pairs_calc_gen.h>
#undef PREFIX

/* #define PREFIX pairs_calc_v2l */
/* #include <pairs_calc_gen.h> */
/* #undef PREFIX */

/* typedef pairs_calc_v2b_arg_t pairs_calc_v2b_tbl_arg_t; */
/* #define PREFIX pairs_calc_v2b_tbl */
/* #include <pairs_calc_gen.h> */
/* #undef PREFIX */

#undef CAT
#undef STR

#if 0
//#include <pairs_calc_a2i.h>
//#include <pairs_calc_a2b.h>
#include <pairs_calc_v2b.h>
//#include <pairs_calc_v2l.h>
//#include <pairs_calc_v2b_tbl.h>
#endif
