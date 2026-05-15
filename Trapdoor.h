
#ifndef TRAPDOOR_H
#define TRAPDOOR_H

#include <stdint.h>
#include <stdbool.h>
#include "config_params.h"

/* 陷门结构：包含 K_LEVEL 个 e_i（每个 2M 维） */
typedef struct {
    int t;                              
    unsigned kprime;                    
    uint8_t  present[K_LEVEL];          
    uint32_t e[K_LEVEL][2 * M_COLS];    
    uint8_t nonce[K_LEVEL][16];
} TrapdoorTW;


typedef TrapdoorTW trapdoor_rec_t;

/**
 * 生成陷门（算法 3.4）
 * @param t         时间段
 * @param fields    字段索引数组（int 类型，取值 0..K_LEVEL-1）
 * @param values    对应的字段值
 * @param kprime    关键字个数
 * @param out       输出陷门结构
 * @return          成功返回 true，失败返回 false
 */
bool Trapdoor_make(int t,
                   const int* fields,
                   const char* const* values,
                   int kprime,
                   TrapdoorTW* out);

/**
 * 保存陷门到文件
 */
bool Trapdoor_save(const TrapdoorTW* tw, const char* path);

/**
 * 从文件读取陷门
 */
bool load_trapdoor_for_t(const char* path, trapdoor_rec_t* out);

/**
 * 从 temp.txt 加载 A_ut / T_ut（用于 SampleLeft）
 */
int load_Aut_from_temp(int t, const char* path, uint32_t A[N_ROWS][M_COLS], int32_t T_A[M_COLS][M_COLS]);

#endif // TRAPDOOR_H
