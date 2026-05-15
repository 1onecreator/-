
#ifndef TRAPGEN_H
#define TRAPGEN_H

#include <stdint.h>
#include <stdbool.h>
#include "config_params.h"

// 生成 A (n×m) 与 m×m 的格基 T_A，满足 A·T_A ≡ 0 (mod q)
bool TrapGen_generate(uint32_t A[N_ROWS][M_COLS],
                      int32_t  T_A[M_COLS][M_COLS]);

// 验证：检查 A·T_A ≡ 0 (mod q)
bool TrapGen_verify(const uint32_t A[N_ROWS][M_COLS],
                    const int32_t  T_A[M_COLS][M_COLS]);

// 写入 public_params.h：覆盖旧的 A / T_A / T_A_SMALL，并写入新 A 与 m×m 的 T_A
bool TrapGen_write_public_params(const char* path,
                                 const uint32_t A[N_ROWS][M_COLS],
                                 const int32_t  T_A[M_COLS][M_COLS]);

#endif 
