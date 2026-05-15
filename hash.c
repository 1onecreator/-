#ifdef q
#undef q
#endif
#ifdef n
#undef n
#endif
#ifdef m
#undef m
#endif

#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#include "config_params.h"
#include "hash.h"

/* ------------ 基本工具 ------------ */

static inline uint32_t mod_q_u32(uint64_t x) {
    uint32_t r = (uint32_t)(x % (uint64_t)Q_MODULUS);
    return r;
}


static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}


static uint64_t bytes_to_seed64(const uint8_t* p, size_t n, uint64_t domain) {
    uint64_t h = 1469598103934665603ULL ^ (domain * 0x9e3779b97f4a7c15ULL);
    for (size_t i = 0; i < n; ++i) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ULL;
        h ^= h >> 23;
        h *= 0x2127599bf4325c37ULL;
        h ^= h >> 47;
    }
    return h ^ (n * 0x9e3779b97f4a7c15ULL);
}


static uint32_t prg_u32(uint64_t* state) {
    *state = splitmix64(*state ? *state : 0x6a09e667f3bcc909ULL);
    return (uint32_t)(*state & 0xffffffffu);
}


static uint32_t sample_mod_q(uint64_t* st) {
    const uint64_t q_mod = (uint64_t)Q_MODULUS;
    const uint64_t R = 1ull << 32;                 
    const uint64_t limit = (R / q_mod) * q_mod;    
    for (;;) {
        uint32_t x = prg_u32(st);
        if ((uint64_t)x < limit) {
            return (uint32_t)((uint64_t)x % q_mod); 
        }
        
    }
}


static void put_u32(uint8_t* tag, size_t* pos, size_t cap, uint32_t x) {
    if (*pos + 4 > cap) return;
    tag[(*pos)++] = (uint8_t)(x);
    tag[(*pos)++] = (uint8_t)(x >> 8);
    tag[(*pos)++] = (uint8_t)(x >> 16);
    tag[(*pos)++] = (uint8_t)(x >> 24);
}

/* ------------ H1: 可逆 m×m ------------
 * 构造：随机置换 P、单位下三角 L、上三角 U（对角随机∈Z_q^*），
 * 最后 M = P·L·U (mod q)。恒可逆。 */

bool H1_make(const uint32_t A[N_ROWS][M_COLS],
             uint32_t t,
             uint32_t out[M_COLS][M_COLS]) {
    uint8_t tag[64]; size_t pos = 0;
    const int sr = (N_ROWS < 8 ? N_ROWS : 8);
    const int sc = (M_COLS < 16 ? M_COLS : 16);
    put_u32(tag, &pos, sizeof(tag), (uint32_t)N_ROWS);
    put_u32(tag, &pos, sizeof(tag), (uint32_t)M_COLS);
    put_u32(tag, &pos, sizeof(tag), t);
    for (int r = 0; r < sr; ++r)
        for (int c = 0; c < sc; ++c)
            put_u32(tag, &pos, sizeof(tag), A[r][c]);
    uint64_t st = bytes_to_seed64(tag, pos, 0x484131ULL); 

    int perm[M_COLS];
    for (int i = 0; i < M_COLS; ++i) perm[i] = i;
    for (int i = M_COLS - 1; i > 0; --i) {
        uint32_t r = prg_u32(&st);
        int j = (int)(r % (uint32_t)(i + 1));
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }

    uint32_t L[M_COLS][M_COLS], U[M_COLS][M_COLS];
    for (int i = 0; i < M_COLS; ++i)
        for (int j = 0; j < M_COLS; ++j)
            L[i][j] = U[i][j] = 0u;

    for (int i = 0; i < M_COLS; ++i) {
        L[i][i] = 1u;
        for (int j = 0; j < i; ++j)
            L[i][j] = sample_mod_q(&st);
    }
    for (int i = 0; i < M_COLS; ++i) {
        uint32_t d = 0; while (d == 0) d = sample_mod_q(&st);
        U[i][i] = d;
        for (int j = i + 1; j < M_COLS; ++j)
            U[i][j] = sample_mod_q(&st);
    }

    uint32_t M1[M_COLS][M_COLS];
    for (int i = 0; i < M_COLS; ++i)
        for (int j = 0; j < M_COLS; ++j) {
            uint64_t acc = 0;
            for (int k = 0; k < M_COLS; ++k)
                acc += (uint64_t)L[i][k] * (uint64_t)U[k][j];
            M1[i][j] = mod_q_u32(acc);
        }

    for (int i = 0; i < M_COLS; ++i) {
        int pi = perm[i];
        for (int j = 0; j < M_COLS; ++j)
            out[i][j] = M1[pi][j];
    }

    bool all_zero = true;
    for (int i = 0; i < M_COLS && all_zero; ++i)
        for (int j = 0; j < M_COLS && all_zero; ++j)
            if (out[i][j] != 0u) all_zero = false;
    if (all_zero) out[0][0] = 1u;
    return true;
}

/* ------------ H2: 均匀 n×n，第二参数显式为 t ∈ N ------------ */
void H2_make(const uint8_t* msg, size_t msg_len,
             uint32_t t,
             uint32_t out[N_ROWS][N_ROWS]) {
    
    uint8_t tag[52]; size_t pos = 0;
    put_u32(tag, &pos, sizeof(tag), (uint32_t)msg_len);
    put_u32(tag, &pos, sizeof(tag), t);

    size_t tlim = msg_len < 12 ? msg_len : 12;
    for (size_t u = 0; u < tlim; ++u) tag[pos++] = msg[u];
    for (size_t u = 0; u < tlim; ++u) tag[pos++] = msg[msg_len - 1 - u];

    uint64_t st = bytes_to_seed64(tag, pos, 0x484132ULL); 

    bool any_nonzero = false;
    for (int r = 0; r < N_ROWS; ++r) {
        for (int c = 0; c < N_ROWS; ++c) {
            uint32_t v = sample_mod_q(&st);
            out[r][c] = v;
            if (v != 0u) any_nonzero = true;
        }
    }
    if (!any_nonzero) out[0][0] = 1u;
}

/* ------------ H3: 均匀 n 维，非全零 ------------ */
void H3_make(const uint32_t* mat_row_major, size_t m_rows, size_t cols,
             const uint8_t* beta_bytes, size_t beta_len,  
             uint32_t out[N_ROWS]) {

    uint8_t tag[96]; size_t pos = 0;
    put_u32(tag, &pos, sizeof(tag), (uint32_t)m_rows);
    put_u32(tag, &pos, sizeof(tag), (uint32_t)cols);
    put_u32(tag, &pos, sizeof(tag), (uint32_t)beta_len);

    
    const size_t total = m_rows * cols;
    for (size_t t = 0; t < total; ++t) {
        put_u32(tag, &pos, sizeof(tag), mat_row_major[t]);
        if (pos + 4 > sizeof(tag)) {                
            uint64_t st = bytes_to_seed64(tag, pos, 0x484133ULL);
            
            pos = 0;
            put_u32(tag, &pos, sizeof(tag), (uint32_t)st);
            put_u32(tag, &pos, sizeof(tag), (uint32_t)(st >> 32));
        }
    }

    
    for (size_t i = 0; i < beta_len; ++i) {
        tag[pos++] = beta_bytes[i];
        if (pos >= sizeof(tag)) {                   
            uint64_t st = bytes_to_seed64(tag, pos, 0x484133ULL);
            pos = 0;
            put_u32(tag, &pos, sizeof(tag), (uint32_t)st);
            put_u32(tag, &pos, sizeof(tag), (uint32_t)(st >> 32));
        }
    }

    
    uint64_t st = bytes_to_seed64(tag, pos, 0x484133ULL); 

    
    bool any_nonzero = false;
    for (int i = 0; i < N_ROWS; ++i) {
        uint32_t v = sample_mod_q(&st);
        out[i] = v;
        if (v != 0u) any_nonzero = true;
    }
    if (!any_nonzero) out[0] = 1u;
}
