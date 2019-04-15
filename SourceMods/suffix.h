#define CAT(x, y) __CAT__(x, y)
#define __CAT__(x, y) x ## y
#define STR(x) __STR__(x)
#define __STR__(x) #x

#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
#define SUFFIX _efv_2cut_fs
#else
#define SUFFIX _efv_2cut
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
#define SUFFIX _efv_fs
#else
#define SUFFIX _efv
#endif /*FSWITCH*/
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
#define SUFFIX _fv_2cut_fs
#else
#define SUFFIX _fv_2cut
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
#define SUFFIX _fv_fs
#else
#define SUFFIX _fv
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
#define SUFFIX _f_2cut_fs
#else
#define SUFFIX _f_2cut
#endif /* FSWITCH */
#else //NOVEC_2CUT
#ifdef FSWITCH
#define SUFFIX _f_fs
#else
#define SUFFIX _f
#endif //FSWITCH
#endif //CALC_NOVEC_2CUT
#endif //CALC_F

