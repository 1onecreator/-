
#ifndef INDEX_AND_U_H
#define INDEX_AND_U_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "config_params.h"


#ifndef K_FIELDS
#define K_FIELDS KEYWORD_FIELD_COUNT
#endif


#if __STDC_VERSION__ >= 201112L
_Static_assert(K_FIELDS > 0, "KEYWORD_FIELD_COUNT must be > 0");
#endif

// 生成：A_all[k][n][m] 与 B[n][m]（元素均匀分布于 Z_q）
bool IndexAndU_generate(uint32_t A_all[K_FIELDS][N_ROWS][M_COLS],
                        uint32_t B[N_ROWS][M_COLS]);

/* 写入 ./public_params.h
 * - 若文件不存在：创建并写入（含 include guard）
 * - 若文件存在：先查找并替换任何旧的
 *       static const uint32_t A_I[...]
 *       static const uint32_t B_mat[...]
 *   两个块，再把新的 A_I/B_mat 插到 #endif 之前
 */
bool IndexAndU_write_public_params(const char* path,
                                   const uint32_t A_all[K_FIELDS][N_ROWS][M_COLS],
                                   const uint32_t B[N_ROWS][M_COLS]);

#endif 
