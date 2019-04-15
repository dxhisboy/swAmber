#define PAIRS_CALC_EFV PAIRS_CALC_EFV
#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include <nb_energy_cpe_impl.h>
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include <nb_energy_cpe_impl.h>
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define FSWITCH
#include <nb_energy_cpe_impl.h>
#undef FSWITCH
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define FSWITCH
#include <nb_energy_cpe_impl.h>
#undef FSWITCH
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define FSWITCH
#include <nb_energy_cpe_impl.h>
#undef FSWITCH
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define FSWITCH
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef FSWITCH
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define FSWITCH
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include <nb_energy_cpe_impl.h>
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef FSWITCH
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#define FSWITCH
#include <nb_energy_cpe_impl.h>
#undef FSWITCH
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F
static void (*get_nb_energy_funcs[16])(nb_energy_arg_t *arg) = {
  /* E, V, C, S */
  /* 0  0  0  0 */ get_nb_energy_f_2cut,
  /* 0  0  0  1 */ get_nb_energy_f_2cut_fs,
  /* 0  0  1  0 */ get_nb_energy_f,
  /* 0  0  1  1 */ get_nb_energy_f_fs,
  /* 0  1  0  0 */ get_nb_energy_fv_2cut,
  /* 0  1  0  1 */ get_nb_energy_fv_2cut_fs,
  /* 0  1  1  0 */ get_nb_energy_fv,
  /* 0  1  1  1 */ get_nb_energy_fv_fs,
  /* 1  0  0  0 */ get_nb_energy_efv_2cut,
  /* 1  0  0  1 */ get_nb_energy_efv_2cut_fs,
  /* 1  0  1  0 */ get_nb_energy_efv,
  /* 1  0  1  1 */ get_nb_energy_efv_fs,
  /* 1  1  0  0 */ get_nb_energy_efv_2cut,
  /* 1  1  0  1 */ get_nb_energy_efv_2cut_fs,
  /* 1  1  1  0 */ get_nb_energy_efv,
  /* 1  1  1  1 */ get_nb_energy_efv_fs
};
