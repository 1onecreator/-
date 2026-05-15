
#ifndef GLOBAL_SAMPLING_H
#define GLOBAL_SAMPLING_H

#include <stdint.h>
#include <stdbool.h>
#include "config_params.h"

// 全局采样空间初始化
bool global_sampling_init(void);

// 全局采样空间清理
void global_sampling_cleanup(void);

// 获取预构建的采样空间
bool get_global_gauss_X(double sigma, void** gauss_table);
bool get_global_gauss_x(double sigma, void** gauss_table);
bool get_global_sampleleft_gauss(double sigma, void** gauss_table);

// 设置全局随机数种子（用于采样）
void global_rng_seed(uint64_t seed);
void global_sampleleft_rng_seed(uint64_t seed);

// 全局采样函数
int32_t global_sample_gauss(double sigma);
int32_t global_sampleleft_dgauss_sample(void* gauss_table);
uint32_t global_sample_mod_q(void);

#endif 
