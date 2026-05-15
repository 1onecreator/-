
#ifndef PUBLIC_PARAMS_H
#define PUBLIC_PARAMS_H

#include <stdint.h>
#include "config_params.h"

#ifdef __cplusplus
extern "C" {
#endif


extern uint32_t A[N_ROWS][M_COLS];
extern int32_t  T_A[M_COLS][M_COLS];
extern const uint32_t B_mat[N_ROWS][M_COLS];
extern const uint32_t U_vec[N_ROWS];

// 每个字段 i 的 A_i（不随 t 变）
extern uint32_t A_I[K_FIELDS][N_ROWS][M_COLS];
extern uint32_t Ao[N_ROWS][M_COLS];
extern int32_t  T_Ao[M_COLS][M_COLS];
extern uint32_t Uu_inv_cache[M_COLS][M_COLS];
extern uint32_t Uo_inv_cache[M_COLS][M_COLS];
extern double SIGMA_SEQ[T_PERIOD];
extern double DELTA_SEQ[T_PERIOD];
extern int          Uo_inv_cache_t;
extern int          Uu_inv_cache_t;


#ifdef __cplusplus
}
#endif
#endif 
