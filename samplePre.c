
#include "samplePre.h"
#include "global_sampling.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config_params.h"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

// ===== 轻量 PRNG + 高斯采样 =====
static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static inline uint64_t xorshift64star(void){
    uint64_t x=rng_state; x^=x>>12; x^=x<<25; x^=x>>27; rng_state=x;
    return x*0x2545F4914F6CDD1DULL;
}
static inline void rng_seed(uint64_t s){
    rng_state = s ? s : 0x9e3779b97f4a7c15ULL;
    for(int i=0;i<8;++i) (void)xorshift64star();
}
static int32_t sample_gauss(double sigma){
    double u1=((xorshift64star()>>11)&0x1fffff)/(double)(1u<<21);
    double u2=((xorshift64star()>>11)&0x1fffff)/(double)(1u<<21);
    if(u1<=0.0) u1=1e-12;
    double z = sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
    double v=z*sigma, cap=8.0*sigma;
    if(v>cap) v=cap; if(v<-cap) v=-cap;
    long long r=llround(v);
    if(r> 2147483647LL) r= 2147483647LL;
    if(r<-2147483647LL) r=-2147483647LL;
    return (int32_t)r;
}

// ===== 模运算 =====
static inline uint32_t mod_q_i64(long long x){
    long long r=x%(long long)Q_MODULUS; if(r<0) r+=Q_MODULUS; return (uint32_t)r;
}
static inline uint32_t mod_q_i128(__int128 x){
    long long r=(long long)(x%(__int128)Q_MODULUS); if(r<0) r+=Q_MODULUS; return (uint32_t)r;
}

// ===== 矩阵求逆（高斯消元，模 q）=====
static int matrix_inv_modq(const uint32_t A[N_ROWS][N_ROWS], uint32_t inv[N_ROWS][N_ROWS]){
    
    uint32_t aug[N_ROWS][2*N_ROWS];
    for(int i=0;i<N_ROWS;++i){
        for(int j=0;j<N_ROWS;++j) aug[i][j] = A[i][j];
        for(int j=0;j<N_ROWS;++j) aug[i][N_ROWS+j] = (i==j) ? 1u : 0u;
    }
    
    
    for(int i=0;i<N_ROWS;++i){
        
        int pivot = -1;
        for(int k=i;k<N_ROWS;++k){
            if(aug[k][i] != 0){ pivot=k; break; }
        }
        if(pivot == -1) return 0; 
        
        
        if(pivot != i){
            for(int j=0;j<2*N_ROWS;++j){
                uint32_t tmp = aug[i][j]; aug[i][j] = aug[pivot][j]; aug[pivot][j] = tmp;
            }
        }
        
        
        uint32_t inv_pivot = 1;
        for(uint32_t k=1;k<Q_MODULUS;++k){
            if((k * aug[i][i]) % Q_MODULUS == 1){ inv_pivot = k; break; }
        }
        for(int j=0;j<2*N_ROWS;++j) aug[i][j] = (aug[i][j] * inv_pivot) % Q_MODULUS;
        
        
        for(int k=0;k<N_ROWS;++k){
            if(k==i) continue;
            uint32_t factor = aug[k][i];
            for(int j=0;j<2*N_ROWS;++j){
                aug[k][j] = mod_q_i64((long long)aug[k][j] - (long long)factor * aug[i][j]);
            }
        }
    }
    
    
    for(int i=0;i<N_ROWS;++i){
        for(int j=0;j<N_ROWS;++j) inv[i][j] = aug[i][N_ROWS+j];
    }
    return 1;
}

bool SamplePre_run(const uint32_t A[N_ROWS][M_COLS],
                   const int32_t  T_A[M_COLS][M_COLS],
                   const uint32_t u[N_ROWS],
                   double sigma,
                   uint64_t seed,  
                   uint32_t e[M_COLS]){
    global_rng_seed(seed);  

    // 1) A1^{-1}
    uint32_t A1[N_ROWS][N_ROWS];
    for(int i=0;i<N_ROWS;++i) for(int j=0;j<N_ROWS;++j) A1[i][j] = A[i][j];
    uint32_t A1_inv[N_ROWS][N_ROWS];
    if(!matrix_inv_modq(A1, A1_inv)) return false;

    // 2) 采样噪声 e2 ~ D_{Z^s,σ}
    uint32_t e2[M_COLS - N_ROWS];
    for(int j=0;j<M_COLS-N_ROWS;++j) e2[j] = (uint32_t)mod_q_i64(global_sample_gauss(sigma));

    // 3) 计算 e1 = A1^{-1} (u - A2 e2)
    uint32_t A2e2[N_ROWS];
    for(int i=0;i<N_ROWS;++i){
        __int128 acc=0;
        for(int j=0;j<M_COLS-N_ROWS;++j) acc += (__int128)A[i][N_ROWS+j] * e2[j];
        A2e2[i] = mod_q_i128(acc);
    }
    uint32_t target[N_ROWS];
    for(int i=0;i<N_ROWS;++i) target[i] = mod_q_i64((long long)u[i] - (long long)A2e2[i]);
    
    uint32_t e1[N_ROWS];
    for(int i=0;i<N_ROWS;++i){
        __int128 acc=0;
        for(int j=0;j<N_ROWS;++j) acc += (__int128)A1_inv[i][j] * target[j];
        e1[i] = mod_q_i128(acc);
    }

    // 4) 拼接 e = [e1 || e2]
    for(int j=0;j<N_ROWS;++j) e[j] = e1[j];
    for(int j=0;j<M_COLS-N_ROWS;++j) e[N_ROWS+j] = e2[j];

    return true;
}

// 验证：检查 A·e ≡ u (mod q)
bool SamplePre_verify(const uint32_t A[N_ROWS][M_COLS],
                      const uint32_t e[M_COLS],
                      const uint32_t u[N_ROWS]){
    for(int i=0;i<N_ROWS;++i){
        __int128 acc=0;
        for(int j=0;j<M_COLS;++j) acc += (__int128)A[i][j] * e[j];
        if(mod_q_i128(acc) != (u[i]%Q_MODULUS)) return false;
    }
    return true;
}
