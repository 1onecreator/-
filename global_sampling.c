
#define _GNU_SOURCE 
#include "global_sampling.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

// ==================== 数据结构定义 ====================


typedef struct { double sigma; int K; double *cdf; int ready; } dgauss_tab_t;


typedef struct { double sigma; int K; double *cdf; int ready; } sampleleft_dgauss_t;


typedef struct {
    double sigma;
    void* gauss_table;
    bool is_sampleleft;  
} sampling_cache_entry_t;


typedef struct {
    
    sampling_cache_entry_t* cipher_gauss_X_cache;
    sampling_cache_entry_t* cipher_gauss_x_cache;
    sampling_cache_entry_t* sampleleft_gauss_cache;
    int cipher_X_cache_size;
    int cipher_x_cache_size;
    int sampleleft_cache_size;
    
    
    uint64_t global_rng_state;
    uint64_t global_sampleleft_rng_state;
    
    
    pthread_mutex_t cache_mutex;
    
    
    bool initialized;
} global_sampling_manager_t;

static global_sampling_manager_t g_sampling_manager = {0};

// ==================== 随机数生成器 ====================


static uint64_t xorshift64star(void) {
    uint64_t x = g_sampling_manager.global_rng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    g_sampling_manager.global_rng_state = x;
    return x * 0x2545F4914F6CDD1DULL;
}


static uint64_t sampleleft_xorshift64star(void) {
    uint64_t x = g_sampling_manager.global_sampleleft_rng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    g_sampling_manager.global_sampleleft_rng_state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

// ==================== 高斯分布构建函数 ====================


static void dgauss_build(double sigma, dgauss_tab_t* T) {
    if (T->ready && fabs(T->sigma - sigma) < 1e-12) return;
    if (T->cdf) { free(T->cdf); T->cdf = NULL; }
    
    int K = (int)ceil(12.0 * sigma);
    if (K < 8) K = 8;
    T->cdf = (double*)malloc((size_t)(K + 1) * sizeof(double));
    T->sigma = sigma;
    T->K = K;
    
    double sum = 0.0;
    for (int i = 0; i <= K; i++) {
        double x = (double)i;
        sum += exp(-(x * x) / (2.0 * sigma * sigma));
        T->cdf[i] = sum;
    }
    for (int i = 0; i <= K; i++) {
        T->cdf[i] /= sum;
    }
    T->ready = 1;
}


static void sampleleft_dgauss_build(double sigma, sampleleft_dgauss_t* T) {
    if (T->ready && fabs(T->sigma - sigma) < 1e-12) return;
    if (T->cdf) { free(T->cdf); T->cdf = NULL; }
    
    int K = (int)ceil(12.0 * sigma);
    if (K < 8) K = 8;
    T->cdf = (double*)malloc((size_t)(K + 1) * sizeof(double));
    T->sigma = sigma;
    T->K = K;
    
    double sum = 0.0;
    for (int i = 0; i <= K; i++) {
        double x = (double)i;
        sum += exp(-(x * x) / (2.0 * sigma * sigma));
        T->cdf[i] = sum;
    }
    for (int i = 0; i <= K; i++) {
        T->cdf[i] /= sum;
    }
    T->ready = 1;
}

// ==================== 采样函数 ====================


static int32_t sample_gauss(double sigma) {
    double u1 = ((xorshift64star() >> 11) & 0x1fffff) / (double)(1u << 21);
    double u2 = ((xorshift64star() >> 11) & 0x1fffff) / (double)(1u << 21);
    if (u1 <= 0.0) u1 = 1e-12;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    double v = z * sigma, cap = 8.0 * sigma;
    if (v > cap) v = cap;
    if (v < -cap) v = -cap;
    long long r = llround(v);
    return (int32_t)r;
}


static int32_t dgauss_sample_sym(dgauss_tab_t* T) {
    if (!T || !T->ready) return 0;
    
    double u = ((xorshift64star() >> 11) & 0x1fffff) / (double)(1u << 21);
    if (u < 0.5) {
        u = 2.0 * u;
        for (int i = 0; i <= T->K; i++) {
            if (u <= T->cdf[i]) return (int32_t)i;
        }
        return (int32_t)T->K;
    } else {
        u = 2.0 * (1.0 - u);
        for (int i = 0; i <= T->K; i++) {
            if (u <= T->cdf[i]) return -(int32_t)i;
        }
        return -(int32_t)T->K;
    }
}


static int32_t sampleleft_dgauss_sample(sampleleft_dgauss_t* T) {
    if (!T || !T->ready) return 0;
    
    double u = ((sampleleft_xorshift64star() >> 11) & 0x1fffff) / (double)(1u << 21);
    if (u < 0.5) {
        u = 2.0 * u;
        for (int i = 0; i <= T->K; i++) {
            if (u <= T->cdf[i]) return (int32_t)i;
        }
        return (int32_t)T->K;
    } else {
        u = 2.0 * (1.0 - u);
        for (int i = 0; i <= T->K; i++) {
            if (u <= T->cdf[i]) return -(int32_t)i;
        }
        return -(int32_t)T->K;
    }
}


static uint32_t sample_mod_q(void) {
    const uint64_t qmod = (uint64_t)Q_MODULUS;
    const uint64_t RANGE = 1ull << 32;
    const uint64_t limit = (RANGE / qmod) * qmod;
    for (;;) {
        uint32_t x = (uint32_t)(xorshift64star() & 0xffffffffu);
        if ((uint64_t)x < limit) return (uint32_t)((uint64_t)x % qmod);
    }
}

// ==================== 缓存管理 ====================


static sampling_cache_entry_t* find_cache_entry(sampling_cache_entry_t* cache, int cache_size, double sigma, bool is_sampleleft) {
    for (int i = 0; i < cache_size; i++) {
        if (fabs(cache[i].sigma - sigma) < 1e-12 && cache[i].is_sampleleft == is_sampleleft) {
            return &cache[i];
        }
    }
    return NULL;
}


static bool add_cache_entry(sampling_cache_entry_t** cache, int* cache_size, double sigma, void* gauss_table, bool is_sampleleft) {
    *cache = (sampling_cache_entry_t*)realloc(*cache, (*cache_size + 1) * sizeof(sampling_cache_entry_t));
    if (!*cache) return false;
    
    (*cache)[*cache_size].sigma = sigma;
    (*cache)[*cache_size].gauss_table = gauss_table;
    (*cache)[*cache_size].is_sampleleft = is_sampleleft;
    (*cache_size)++;
    return true;
}

// ==================== 公共API ====================

bool global_sampling_init(void) {
    if (g_sampling_manager.initialized) return true;
    
    
    if (pthread_mutex_init(&g_sampling_manager.cache_mutex, NULL) != 0) {
        fprintf(stderr, "[global_sampling] Failed to initialize mutex\n");
        return false;
    }
    
    
    g_sampling_manager.global_rng_state = 0x9e3779b97f4a7c15ULL;
    g_sampling_manager.global_sampleleft_rng_state = 0x9e3779b97f4a7c15ULL;
    
    
    for (int i = 0; i < 8; i++) {
        (void)xorshift64star();
        (void)sampleleft_xorshift64star();
    }
    
    // 预构建常用的高斯分布表
    printf("[global_sampling] Initializing sampling spaces...\n");
    
    // 预构建常用的sigma值
    double common_sigmas[] = {1.0, 2.0, 3.0, 5.0, 10.0};  
    int num_sigmas = sizeof(common_sigmas) / sizeof(common_sigmas[0]);
    
    for (int i = 0; i < num_sigmas; i++) {
        double sigma = common_sigmas[i];
        
        
        dgauss_tab_t* gauss_X = (dgauss_tab_t*)malloc(sizeof(dgauss_tab_t));
        memset(gauss_X, 0, sizeof(dgauss_tab_t));
        dgauss_build(sigma, gauss_X);
        add_cache_entry(&g_sampling_manager.cipher_gauss_X_cache, 
                       &g_sampling_manager.cipher_X_cache_size, sigma, gauss_X, false);
        
        
        dgauss_tab_t* gauss_x = (dgauss_tab_t*)malloc(sizeof(dgauss_tab_t));
        memset(gauss_x, 0, sizeof(dgauss_tab_t));
        dgauss_build(sigma, gauss_x);
        add_cache_entry(&g_sampling_manager.cipher_gauss_x_cache, 
                       &g_sampling_manager.cipher_x_cache_size, sigma, gauss_x, false);
        
        
        sampleleft_dgauss_t* sampleleft_gauss = (sampleleft_dgauss_t*)malloc(sizeof(sampleleft_dgauss_t));
        memset(sampleleft_gauss, 0, sizeof(sampleleft_dgauss_t));
        sampleleft_dgauss_build(sigma, sampleleft_gauss);
        add_cache_entry(&g_sampling_manager.sampleleft_gauss_cache, 
                       &g_sampling_manager.sampleleft_cache_size, sigma, sampleleft_gauss, true);
    }
    
    g_sampling_manager.initialized = true;
    printf("[global_sampling] Sampling spaces initialized successfully\n");
    printf("[global_sampling] Cached %d GAUSS_X tables, %d GAUSS_x tables, %d sampleleft tables\n",
           g_sampling_manager.cipher_X_cache_size,
           g_sampling_manager.cipher_x_cache_size,
           g_sampling_manager.sampleleft_cache_size);
    
    return true;
}

void global_sampling_cleanup(void) {
    if (!g_sampling_manager.initialized) return;
    
    pthread_mutex_lock(&g_sampling_manager.cache_mutex);
    
    
    for (int i = 0; i < g_sampling_manager.cipher_X_cache_size; i++) {
        dgauss_tab_t* table = (dgauss_tab_t*)g_sampling_manager.cipher_gauss_X_cache[i].gauss_table;
        if (table && table->cdf) free(table->cdf);
        free(table);
    }
    free(g_sampling_manager.cipher_gauss_X_cache);
    
    
    for (int i = 0; i < g_sampling_manager.cipher_x_cache_size; i++) {
        dgauss_tab_t* table = (dgauss_tab_t*)g_sampling_manager.cipher_gauss_x_cache[i].gauss_table;
        if (table && table->cdf) free(table->cdf);
        free(table);
    }
    free(g_sampling_manager.cipher_gauss_x_cache);
    
    
    for (int i = 0; i < g_sampling_manager.sampleleft_cache_size; i++) {
        sampleleft_dgauss_t* table = (sampleleft_dgauss_t*)g_sampling_manager.sampleleft_gauss_cache[i].gauss_table;
        if (table && table->cdf) free(table->cdf);
        free(table);
    }
    free(g_sampling_manager.sampleleft_gauss_cache);
    
    pthread_mutex_unlock(&g_sampling_manager.cache_mutex);
    pthread_mutex_destroy(&g_sampling_manager.cache_mutex);
    
    memset(&g_sampling_manager, 0, sizeof(g_sampling_manager));
    printf("[global_sampling] Sampling spaces cleaned up\n");
}

bool get_global_gauss_X(double sigma, void** gauss_table) {
    if (!g_sampling_manager.initialized) return false;
    
    pthread_mutex_lock(&g_sampling_manager.cache_mutex);
    sampling_cache_entry_t* entry = find_cache_entry(g_sampling_manager.cipher_gauss_X_cache, 
                                                    g_sampling_manager.cipher_X_cache_size, sigma, false);
    
    if (!entry) {
        
        dgauss_tab_t* new_table = (dgauss_tab_t*)malloc(sizeof(dgauss_tab_t));
        memset(new_table, 0, sizeof(dgauss_tab_t));
        dgauss_build(sigma, new_table);
        add_cache_entry(&g_sampling_manager.cipher_gauss_X_cache, 
                       &g_sampling_manager.cipher_X_cache_size, sigma, new_table, false);
        *gauss_table = new_table;
    } else {
        *gauss_table = entry->gauss_table;
    }
    
    pthread_mutex_unlock(&g_sampling_manager.cache_mutex);
    return true;
}

bool get_global_gauss_x(double sigma, void** gauss_table) {
    if (!g_sampling_manager.initialized) return false;
    
    pthread_mutex_lock(&g_sampling_manager.cache_mutex);
    sampling_cache_entry_t* entry = find_cache_entry(g_sampling_manager.cipher_gauss_x_cache, 
                                                    g_sampling_manager.cipher_x_cache_size, sigma, false);
    
    if (!entry) {
        
        dgauss_tab_t* new_table = (dgauss_tab_t*)malloc(sizeof(dgauss_tab_t));
        memset(new_table, 0, sizeof(dgauss_tab_t));
        dgauss_build(sigma, new_table);
        add_cache_entry(&g_sampling_manager.cipher_gauss_x_cache, 
                       &g_sampling_manager.cipher_x_cache_size, sigma, new_table, false);
        *gauss_table = new_table;
    } else {
        *gauss_table = entry->gauss_table;
    }
    
    pthread_mutex_unlock(&g_sampling_manager.cache_mutex);
    return true;
}

bool get_global_sampleleft_gauss(double sigma, void** gauss_table) {
    if (!g_sampling_manager.initialized) return false;
    
    pthread_mutex_lock(&g_sampling_manager.cache_mutex);
    sampling_cache_entry_t* entry = find_cache_entry(g_sampling_manager.sampleleft_gauss_cache, 
                                                    g_sampling_manager.sampleleft_cache_size, sigma, true);
    
    if (!entry) {
        
        sampleleft_dgauss_t* new_table = (sampleleft_dgauss_t*)malloc(sizeof(sampleleft_dgauss_t));
        memset(new_table, 0, sizeof(sampleleft_dgauss_t));
        sampleleft_dgauss_build(sigma, new_table);
        add_cache_entry(&g_sampling_manager.sampleleft_gauss_cache, 
                       &g_sampling_manager.sampleleft_cache_size, sigma, new_table, true);
        *gauss_table = new_table;
    } else {
        *gauss_table = entry->gauss_table;
    }
    
    pthread_mutex_unlock(&g_sampling_manager.cache_mutex);
    return true;
}

void global_rng_seed(uint64_t seed) {
    g_sampling_manager.global_rng_state = seed ? seed : 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 8; i++) (void)xorshift64star();
}

void global_sampleleft_rng_seed(uint64_t seed) {
    g_sampling_manager.global_sampleleft_rng_state = seed ? seed : 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 8; i++) (void)sampleleft_xorshift64star();
}

int32_t global_sample_gauss(double sigma) {
    return sample_gauss(sigma);
}

int32_t global_sampleleft_dgauss_sample(void* gauss_table) {
    return sampleleft_dgauss_sample((sampleleft_dgauss_t*)gauss_table);
}

uint32_t global_sample_mod_q(void) {
    return sample_mod_q();
}
