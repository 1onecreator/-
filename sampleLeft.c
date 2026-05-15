
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sampleLeft.h"
#include "samplePre.h"
#include "global_sampling.h"
#include "config_params.h"

/* ---------------- 离散高斯采样 ---------------- */
static uint64_t sampleleft_rng_state = 0x9e3779b97f4a7c15ULL;
static inline uint64_t sampleleft_xorshift64star(void){
    uint64_t x=sampleleft_rng_state; x^=x>>12; x^=x<<25; x^=x>>27; sampleleft_rng_state=x;
    return x*0x2545F4914F6CDD1DULL;
}
static inline void sampleleft_rng_seed(uint64_t s){
    sampleleft_rng_state = s ? s : 0x9e3779b97f4a7c15ULL;
    for(int i=0;i<8;++i) (void)sampleleft_xorshift64star();
}
static inline double sampleleft_uniform01(void){
    uint64_t r=sampleleft_xorshift64star(); return (double)((r>>11)*(1.0/9007199254740992.0));
}
typedef struct { double sigma; int K; double *cdf; int ready; } sampleleft_dgauss_t;
static void sampleleft_dgauss_build(double sigma, sampleleft_dgauss_t* T){
    if (T->ready && fabs(T->sigma - sigma) < 1e-12) return;
    if (T->cdf) { free(T->cdf); T->cdf=NULL; }
    int K = (int)ceil(12.0 * sigma); if (K<8) K=8;
    T->cdf = (double*)malloc((size_t)(K+1)*sizeof(double));
    T->sigma=sigma; T->K=K;
    double mass=0.0;
    for (int k=0;k<=K;k++){ double w=exp(-(double)k*(double)k/(2.0*sigma*sigma)); T->cdf[k]=w; mass+=w; }
    double acc=0.0; for(int k=0;k<=K;k++){ acc += T->cdf[k]/mass; T->cdf[k]=acc; }
    if (T->cdf[K]<1.0) T->cdf[K]=1.0;
    T->ready=1;
}
static inline int sampleleft_dgauss_sample(sampleleft_dgauss_t* T){
    double u = sampleleft_uniform01();
    int lo=0,hi=T->K;
    while (lo<hi){ int mid=(lo+hi)/2; if(u<=T->cdf[mid]) hi=mid; else lo=mid+1; }
    int k=lo; if(k==0) return 0;
    return (sampleleft_xorshift64star()&1ULL) ? k : -k;
}
static sampleleft_dgauss_t GAUSS_TAB = {0};

/* ---------------- 模运算 ---------------- */
static inline long long mod_q_i64(long long x){
    long long r=x%(long long)Q_MODULUS; if(r<0) r+=Q_MODULUS; return r;
}
static inline uint32_t mod_q_i128(__int128 x){
    long long r=(long long)(x%(__int128)Q_MODULUS); if(r<0) r+=Q_MODULUS; return (uint32_t)r;
}

bool SampleLeft_run(const uint32_t A[N_ROWS][M_COLS],
                    const int32_t  T_A[M_COLS][M_COLS],
                    const uint32_t G[N_ROWS][M_COLS],
                    const uint32_t u[N_ROWS],
                    double sigma,
                    uint64_t seed,  
                    uint32_t e[2*M_COLS]){
    
    uint64_t seed_eG = seed;
    uint64_t seed_eA = seed ^ 0xDEADBEEFCAFEBABEULL;
    
    // 1) 采样 e_G：使用全局采样空间
    global_sampleleft_rng_seed(seed_eG);
    void* gauss_table = NULL;
    if (!get_global_sampleleft_gauss(sigma, &gauss_table)) return false;
    
    uint32_t eG[M_COLS];
    for(int j=0;j<M_COLS;++j) {
        int sample = global_sampleleft_dgauss_sample(gauss_table);
        eG[j] = (uint32_t)mod_q_i64((long long)sample);
    }

    // 2) 计算 target = u - G·e_G (mod q)
    uint32_t target[N_ROWS];
    for(int i=0;i<N_ROWS;++i){
        __int128 acc=0;
        for(int j=0;j<M_COLS;++j) acc += (__int128)G[i][j]*eG[j];
        uint32_t Ge = mod_q_i128(acc);
        target[i]=mod_q_i64((long long)u[i] - (long long)Ge);
    }

    // 3) 用 SamplePre 采样 e_A：A·e_A ≡ target (mod q)
    uint32_t eA[M_COLS];
    if(!SamplePre_run(A, T_A, target, sigma, seed_eA, eA)) return false;

    // 4) 拼接 e = [e_A || e_G]
    for(int j=0;j<M_COLS;++j) e[j]=eA[j];
    for(int j=0;j<M_COLS;++j) e[M_COLS+j]=eG[j];

    return true;
}


bool SampleLeft_verify(const uint32_t A[N_ROWS][M_COLS],
                       const uint32_t G[N_ROWS][M_COLS],
                       const uint32_t e[2*M_COLS],
                       const uint32_t u[N_ROWS]){
    for (int r=0;r<N_ROWS;r++){
        __int128 acc=0;
        for (int j=0;j<M_COLS;j++) acc += (__int128)A[r][j]*e[j];
        for (int j=0;j<M_COLS;j++) acc += (__int128)G[r][j]*e[M_COLS+j];
        uint32_t lhs = mod_q_i128(acc);
        if (lhs != u[r]) return false;
    }
    return true;
}
