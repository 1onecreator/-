
#ifndef SAMPLE_PRE_H
#define SAMPLE_PRE_H

#include <stdint.h>
#include <stdbool.h>
#include "config_params.h"

// 生成 e \in Z_q^m，使 A·e = u (mod q)
// 说明：使用 T_A 的前 (m-n) 列（短基）进行噪声采样；e 为 [0,q) 表示
bool SamplePre_run(const uint32_t A[N_ROWS][M_COLS],
                   const int32_t  T_A[M_COLS][M_COLS],
                   const uint32_t u[N_ROWS],
                   double sigma,
                   uint64_t seed,  
                   uint32_t e[M_COLS]);

// 验证：检查 A·e ≡ u (mod q)
bool SamplePre_verify(const uint32_t A[N_ROWS][M_COLS],
                      const uint32_t e[M_COLS],
                      const uint32_t u[N_ROWS]);

#endif 
