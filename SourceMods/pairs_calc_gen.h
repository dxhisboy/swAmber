#define HDR <PREFIX.h>
#define PAIRS_CALC_EFV PAIRS_CALC_EFV
#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include HDR
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include HDR
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include HDR
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include HDR
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include HDR
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include HDR
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define FSWITCH
#include HDR
#undef FSWITCH
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define FSWITCH
#include HDR
#undef FSWITCH
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define FSWITCH
#include HDR
#undef FSWITCH
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define FSWITCH
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include HDR
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef FSWITCH
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define FSWITCH
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include HDR
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef FSWITCH
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#define FSWITCH
#include HDR
#undef FSWITCH
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

static void (*CAT(PREFIX, _funcs[16]))(CAT(PREFIX, _arg_t) *arg) = {
  /* E, V, C, S */
  /* 0  0  0  0 */ CAT(PREFIX, _f_2cut),
  /* 0  0  0  1 */ CAT(PREFIX, _f_2cut_fs),
  /* 0  0  1  0 */ CAT(PREFIX, _f),
  /* 0  0  1  1 */ CAT(PREFIX, _f_fs),
  /* 0  1  0  0 */ CAT(PREFIX, _fv_2cut),
  /* 0  1  0  1 */ CAT(PREFIX, _fv_2cut_fs),
  /* 0  1  1  0 */ CAT(PREFIX, _fv),
  /* 0  1  1  1 */ CAT(PREFIX, _fv_fs),
  /* 1  0  0  0 */ CAT(PREFIX, _efv_2cut),
  /* 1  0  0  1 */ CAT(PREFIX, _efv_2cut_fs),
  /* 1  0  1  0 */ CAT(PREFIX, _efv),
  /* 1  0  1  1 */ CAT(PREFIX, _efv_fs),
  /* 1  1  0  0 */ CAT(PREFIX, _efv_2cut),
  /* 1  1  0  1 */ CAT(PREFIX, _efv_2cut_fs),
  /* 1  1  1  0 */ CAT(PREFIX, _efv),
  /* 1  1  1  1 */ CAT(PREFIX, _efv_fs)
};
void CAT(PREFIX, _efv_2cut_fs_c_)(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _efv_2cut_fs_c_)))));
void CAT(PREFIX, _efv_2cut_c_   )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _efv_2cut_c_   )))));
void CAT(PREFIX, _efv_fs_c_     )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _efv_fs_c_     )))));
void CAT(PREFIX, _efv_c_        )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _efv_c_        )))));
void CAT(PREFIX, _fv_2cut_fs_c_ )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _fv_2cut_fs_c_ )))));
void CAT(PREFIX, _fv_2cut_c_    )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _fv_2cut_c_    )))));
void CAT(PREFIX, _fv_fs_c_      )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _fv_fs_c_      )))));
void CAT(PREFIX, _fv_c_         )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _fv_c_         )))));
void CAT(PREFIX, _f_2cut_fs_c_  )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _f_2cut_fs_c_  )))));
void CAT(PREFIX, _f_2cut_c_     )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _f_2cut_c_     )))));
void CAT(PREFIX, _f_fs_c_       )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _f_fs_c_       )))));
void CAT(PREFIX, _f_c_          )(CAT(PREFIX, _arg_t) *) __attribute__((alias(STR(CAT(PREFIX, _f_c_          )))));

#undef HDR
