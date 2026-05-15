
#ifndef GEN_U_H
#define GEN_U_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "config_params.h"

// 生成均匀随机 u ∈ Z_q^n，写入数组 u[N_ROWS]
bool U_generate(uint32_t u[N_ROWS]);

// 将 u 写入 ./public_params.h

bool U_write_public_params(const char* path, const uint32_t u[N_ROWS]);

#endif 
