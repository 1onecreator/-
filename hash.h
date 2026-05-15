#ifndef HASH_H
#define HASH_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "config_params.h"

/* 约定
 * - 全部结果按行主序写入二维数组（out[r][c]）
 * - 全部元素均匀采样于 Z_q，且结果不为 0（矩阵“非全零”、向量“至少一项非零”）
 * - H1 输出 m×m 且在 Z_q 上可逆
 */

/* H1 : Z_q^{n×m} × N -> Z_q^{m×m}   （可逆）
 * 参数：
 *   A : 仅用于吸收熵形成域分隔（不昂贵地 hash 全矩阵）
 *   t : 域分隔（时间段等）
 * 输出：
 *   out[m][m]  （可逆）
 * 返回 true 恒成立
 */

bool H1_make(const uint32_t A[N_ROWS][M_COLS],
             uint32_t t,
             uint32_t out[M_COLS][M_COLS]);


void H2_make(const uint8_t* msg, size_t msg_len,
             uint32_t t,
             uint32_t out[N_ROWS][N_ROWS]);


void H3_make(const uint32_t* mat_row_major, size_t m_rows, size_t cols,
             const uint8_t* beta_bits, size_t l_bits,
             uint32_t out[N_ROWS]);

#endif 
