
#ifndef VERIFY_H
#define VERIFY_H

#include <stdbool.h>
#include "Trapdoor.h"      
#include "config_params.h" 

/*
 * @brief 密文记录的数据结构。
 * @note 此结构现在是公共接口的一部分，以便调用者可以填充它并传递给 verify_search_from_data。
 */
typedef struct {
    int t, record_id, filled;
    uint32_t C_ut[M_COLS][L_COLS];
    uint32_t c_beta[L_COLS];
    uint32_t C_i[K_LEVEL][M_COLS][L_COLS];
    uint8_t  have_Ci[K_LEVEL];
    uint32_t s_t[M_COLS];
    uint32_t h_t[N_ROWS];
} cipher_rec_t;


/* ---------------- 公共验证函数 ---------------- */

/**
 * @brief 从文件加载密文和陷门，并验证是否匹配。
 * @param cipher_path   指向密文文件的路径 (例如 "cipher/cipher_1.txt")。
 * @param trapdoor_path 指向陷门文件的路径。
 * @return true 表示命中，false 表示未命中或发生错误。
 */
bool verify_search(const char* cipher_path, const char* trapdoor_path);

/**
 * @brief 直接使用内存中的数据结构验证密文和陷门是否匹配。
 * @param C      指向已填充的密文记录结构体的指针。
 * @param TW     指向已填充的陷门记录结构体的指针。
 * @param Ao_t   指向与时间 t 对应的 Ao_t 公共矩阵的指针。
 * @return true 表示命中，false 表示未命中。
 */
bool verify_search_from_data(const cipher_rec_t* C, 
                             const trapdoor_rec_t* TW, 
                             const uint32_t Ao_t[N_ROWS][M_COLS]);

#endif 
