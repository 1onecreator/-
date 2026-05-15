
#ifndef PUBLIC_INIT_H
#define PUBLIC_INIT_H

#include <stdbool.h>
#include "config_params.h"

// 计算序列：每个 t 取 max(sigma_pre, sigma_left) 作为 SIGMA_SEQ[t]；
// DELTA_SEQ[t] 取 DELTA_R（若未定义则用 sqrt(n ln q)*OMEGA_M）
bool Public_compute_sequences(double sigma_seq[T_PERIOD],
                              double delta_seq[T_PERIOD]);

// 仅将 SIGMA_SEQ 与 DELTA_SEQ 写入/更新到 public_params.h
bool Public_write_sequences(const char* path,
                            const double sigma_seq[T_PERIOD],
                            const double delta_seq[T_PERIOD]);

#endif 
