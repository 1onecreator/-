


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>

#define _SC_NPROCESSORS_ONLN 0
long sysconf(int name) {
    if (name == _SC_NPROCESSORS_ONLN) {
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        return sysinfo.dwNumberOfProcessors;
    }
    return 1;
}
#else
#include <unistd.h>
#endif

#include "config_params.h"
#include "public_params.h"
#include "hash.h"
#include "verify.h"     
#include "Trapdoor.h"

#define TEMP_PATH_DEFAULT "temp.txt"


#define MALLOC_OR_GOTO(ptr, type, size, label) \
    do { \
        ptr = (type)malloc(size); \
        if (!ptr) { \
            fprintf(stderr, "[verify] ERROR: Memory allocation failed for %s.\n", #ptr); \
            goto label; \
        } \
    } while (0)

/* =========================== 内部工具 =========================== */
static inline long long mod_q_i64(long long x){
    long long r=x%(long long)Q_MODULUS; if(r<0) r+=Q_MODULUS; return r;
}
static inline uint32_t mod_q_i128(__int128 x){
    long long r=(long long)(x%(__int128)Q_MODULUS); if(r<0) r+=Q_MODULUS; return (uint32_t)r;
}
static inline int32_t center_q(uint32_t x){
    uint32_t h = Q_MODULUS >> 1;
    return (x <= h) ? (int32_t)x : (int32_t)x - (int32_t)Q_MODULUS;
}
static inline uint32_t dist_to_half(uint32_t x){
    uint32_t h = Q_MODULUS >> 1;
    int32_t d = center_q((x + Q_MODULUS - h) % Q_MODULUS);
    return (d >= 0) ? (uint32_t)d : (uint32_t)(-d);
}
static int starts_with(const char* s, const char* p){
    for(; *p; ++s, ++p) if (*s != *p) return 0;
    return 1;
}
static int read_uint32_array_multiline(FILE* fp, const char* first_line, uint32_t* out, int max_count){
    int count = 0;
    char line[16384];
    const char* p = strchr(first_line, '{');
    if (p) {
        p++;
        while (*p && count < max_count) {
            while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
            if (*p == '}') return count;
            if (!*p) break;
            unsigned long val; char* endp;
            val = strtoul(p, &endp, 10);
            if (endp != p) { out[count++] = (uint32_t)val; p = endp; } else { break; }
            while (*p == ',' || *p == ' ' || *p == '\t') p++;
        }
        if (strchr(first_line, '}')) return count;
    }
    while (fgets(line, sizeof(line), fp) && count < max_count) {
        p = line;
        while (*p && count < max_count) {
            while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
            if (*p == '}') return count;
            if (!*p) break;
            unsigned long val; char* endp;
            val = strtoul(p, &endp, 10);
            if (endp != p) { out[count++] = (uint32_t)val; p = endp; } else { break; }
            while (*p == ',' || *p == ' ' || *p == '\t') p++;
        }
        if (strchr(line, '}')) return count;
    }
    return count;
}
static int parse_u32_row(const char* line, uint32_t* out, int count){
    const char* p = line;
    for(int i=0; i<count; ++i){
        while(*p && (*p==' '||*p=='\t'||*p=='\n'||*p=='\r')) ++p;
        if(!*p) return i;
        unsigned long x; if(sscanf(p,"%lu",&x)!=1) return i;
        out[i]=(uint32_t)x;
        while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') ++p;
    }
    return count;
}

/* =========================== 文件加载辅助函数 =========================== */
static int load_block_by_t_pair(const char* path, int t, const char* pair_tag,
                                const char* name, int rows, int cols,
                                uint32_t* out_u32)
{
    FILE* fp=fopen(path,"r"); if(!fp) return 0;
    char line[8192]; int in_target=0;
    char tkey[64]; snprintf(tkey,sizeof(tkey),"t = %d", t);
    char pkey[64]; snprintf(pkey,sizeof(pkey),"Pair = %s", pair_tag);
    char nkey[64]; snprintf(nkey,sizeof(nkey),"%st%d (", name, t);
    while (fgets(line,sizeof(line),fp)){
        if (!in_target){
            if (strstr(line,tkey) && strstr(line,pkey)) in_target=1;
            continue;
        }else{
            if (starts_with(line, nkey)){
                for (int r=0;r<rows;r++){
                    if (!fgets(line,sizeof(line),fp)){ fclose(fp); return 0; }
                    if (parse_u32_row(line, out_u32+r*cols, cols)!=cols){ fclose(fp); return 0; }
                }
                fclose(fp); return 1;
            }
        }
    }
    fclose(fp); return 0;
}
static int load_Ao_from_temp(int t, const char* temp_path, uint32_t Ao[N_ROWS][M_COLS]){
    return load_block_by_t_pair(temp_path, t, "Ao/T_Ao", "Ao", N_ROWS, M_COLS, (uint32_t*)Ao);
}
static int load_latest_cipher_for_t(int t, const char* cipher_path, cipher_rec_t* rec){
    FILE* fp = fopen(cipher_path, "r");
    if (!fp) return 0;
    memset(rec, 0, sizeof(*rec));
    rec->t = t;
    char line[8192];
    int found_t = 0, found_record = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (starts_with(line, "t=")) {
            int read_t, read_record_id;
            if (sscanf(line, "t=%d record_id=%d", &read_t, &read_record_id) == 2) {
                if (read_t == t) { found_t = 1; found_record = 1; rec->record_id = read_record_id; }
                else { found_t = 0; found_record = 0; }
            } else if (sscanf(line, "t=%d", &read_t) == 1) {
                if (read_t == t) { found_t = 1; found_record = 0; } else { found_t = 0; found_record = 0; }
            }
        } else if (found_t && !found_record && starts_with(line, "record_id=")) {
            sscanf(line, "record_id=%d", &rec->record_id); found_record = 1;
        } else if (found_t && found_record && starts_with(line, "C_ut=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_ut, M_COLS * L_COLS) != M_COLS * L_COLS) { fclose(fp); return 0; }
        } else if (found_t && found_record && starts_with(line, "c_beta=")) {
            if (read_uint32_array_multiline(fp, line, rec->c_beta, L_COLS) != L_COLS) { fclose(fp); return 0; }
        } else if (found_t && found_record && starts_with(line, "C_1=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[0], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[0] = 1;
        } else if (found_t && found_record && starts_with(line, "C_2=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[1], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[1] = 1;
        } else if (found_t && found_record && starts_with(line, "C_3=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[2], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[2] = 1;
        } else if (found_t && found_record && starts_with(line, "C_4=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[3], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[3] = 1;
        } else if (found_t && found_record && starts_with(line, "C_5=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[4], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[4] = 1;
        } else if (found_t && found_record && starts_with(line, "C_6=")) {
            if (read_uint32_array_multiline(fp, line, (uint32_t*)rec->C_i[5], M_COLS * L_COLS) == M_COLS * L_COLS) rec->have_Ci[5] = 1;
        } else if (found_t && found_record && starts_with(line, "s_t=")) {
            if (read_uint32_array_multiline(fp, line, rec->s_t, M_COLS) != M_COLS) { fclose(fp); return 0; }
        } else if (found_t && found_record && starts_with(line, "h_t=")) {
            if (read_uint32_array_multiline(fp, line, rec->h_t, N_ROWS) != N_ROWS) { fclose(fp); return 0; }
        }
    }
    fclose(fp);
    rec->filled=1; return 1;
}

/* ---------------- β' 和 H3 工具函数 ---------------- */
static void make_beta_closest(const uint32_t* c_beta, const uint32_t* z,
                              uint8_t* beta_bits_out ){
    memset(beta_bits_out, 0, (L_COLS+7)/8);
    uint32_t q_half = Q_MODULUS >> 1;
    uint32_t q_quarter = Q_MODULUS >> 2;
    for (int c=0;c<L_COLS;c++){
        uint32_t beta_prime = mod_q_i64((long long)c_beta[c] - (long long)z[c]);
        uint32_t diff = (beta_prime >= q_half) ? (beta_prime - q_half) : (q_half - beta_prime);
        if (diff < q_quarter) {
            beta_bits_out[c>>3] |= (uint8_t)(1u << (c & 7));
        }
    }
}
static void run_h3_all_ci(const cipher_rec_t* C,
                          const uint8_t* beta_bits, uint32_t* out_h3){
    size_t tot = (size_t)(1 + K_LEVEL) * (size_t)M_COLS * (size_t)L_COLS;
    uint32_t* features = (uint32_t*)malloc(sizeof(uint32_t)*tot);
    if (!features) return;
    size_t p = 0;
    for (int c=0;c<L_COLS;c++) for (int r=0;r<M_COLS;r++) features[p++] = C->C_ut[r][c];
    for (int i=0;i<K_LEVEL;i++) for (int c=0;c<L_COLS;c++) for (int r=0;r<M_COLS;r++) features[p++] = C->C_i[i][r][c];
    fprintf(stderr,"[diag] H3 input: C_ut + all %d C_i, total features=%zu\n", K_LEVEL, tot);
    H3_make(features, 1u, (uint32_t)tot, beta_bits, (uint32_t)((L_COLS+7)/8), out_h3);
    free(features);
}

/* =========================== 并行计算 z 的工作函数 =========================== */
typedef struct {
    int thread_id;
    int num_threads;
    const trapdoor_rec_t *TW;
    const cipher_rec_t *C;
    uint32_t *output_z;
} CalcZArgs;

void* calculate_z_worker(void* arg) {
    CalcZArgs* args = (CalcZArgs*)arg;
    const trapdoor_rec_t* TW = args->TW;
    const cipher_rec_t* C = args->C;
    int start_c = (args->thread_id * L_COLS) / args->num_threads;
    int end_c = ((args->thread_id + 1) * L_COLS) / args->num_threads;
    for (int c = start_c; c < end_c; c++) {
        __int128 total_z_c = 0;
        for (int i = 0; i < K_LEVEL; i++) {
            if (TW->present[i]) {
                const uint32_t* eA = &TW->e[i][0];
                const uint32_t* eG = &TW->e[i][M_COLS];
                __int128 a1 = 0, a2 = 0;
                for (int r = 0; r < M_COLS; r++) a1 += (__int128)C->C_ut[r][c] * eA[r];
                for (int r = 0; r < M_COLS; r++) a2 += (__int128)C->C_i[i][r][c] * eG[r];
                total_z_c += a1 + a2;
            }
        }
        args->output_z[c] = mod_q_i128(total_z_c);
    }
    return NULL;
}

/* =========================== 核心验证逻辑 =========================== */
static bool verify_core(const trapdoor_rec_t* TW, const cipher_rec_t* C, const uint32_t Ao_unused[N_ROWS][M_COLS])
{
    (void)Ao_unused; 

    /* 0) 确定线程数（用于并行计算 z） */
    long num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    int num_threads = (num_cores > 1) ? (int)num_cores : 2;
    fprintf(stderr, "[verify_core] Using %d threads for parallel computation (z only).\n", num_threads);
    pthread_t threads[num_threads];

    /* 1) 检查所有 K_LEVEL 个 C_i 是否都存在 */
    for (int i=0;i<K_LEVEL;i++){
        if (!C->have_Ci[i]){ 
            fprintf(stderr,"[verify_core] missing C_%d in cipher (expected all %d C_i)\n", i+1, K_LEVEL); 
            return false; 
        }
    }

    /* 2) 并行计算 z = T_{β'}(C) */
    uint32_t z[L_COLS];
    {
        fprintf(stderr,"[verify_core] Calculating z vector in parallel...\n");
        CalcZArgs z_args[num_threads];
        for (int i = 0; i < num_threads; i++) {
            z_args[i] = (CalcZArgs){ .thread_id = i, .num_threads = num_threads, .TW = TW, .C = C, .output_z = z };
            pthread_create(&threads[i], NULL, calculate_z_worker, &z_args[i]);
        }
        for (int i = 0; i < num_threads; i++) {
            pthread_join(threads[i], NULL);
        }
    }
    
    /* 3) 由 (c_beta - z) 进行最近星座判决得到 β' 的位串 */
    uint8_t beta_bits[(L_COLS+7)/8];
    make_beta_closest(C->c_beta, z, beta_bits);

    /* 4) 计算 H3(C, β') 作为 LHS */
    uint32_t h3[N_ROWS];
    run_h3_all_ci(C, beta_bits, h3);

    /* 5) 最终与密文中的 h_t 直接比较 */
    for (int r=0;r<N_ROWS;r++){
        if (h3[r] != C->h_t[r]){
            fprintf(stderr,"[verify_core] ❌ 未命中：H3(C,β') != h_t(cipher)\n");
            return false;
        }
    }

    fprintf(stderr,"[verify_core] ✅ 命中成功！H3(C,β') == h_t(cipher)\n");
    return true;
}

/* =========================== 公共接口实现 =========================== */
/**
 * 接口1：从文件加载数据并验证
 */
bool verify_search(const char* cipher_path, const char* trapdoor_path) {
    fprintf(stderr, "[verify] Starting verification from files...\n");

    bool result = false; 
    trapdoor_rec_t* TW = NULL;
    cipher_rec_t* C = NULL;
    uint32_t (*Ao)[M_COLS] = NULL; 

    MALLOC_OR_GOTO(TW, trapdoor_rec_t*, sizeof(trapdoor_rec_t), cleanup);
    MALLOC_OR_GOTO(C,  cipher_rec_t*,  sizeof(cipher_rec_t),   cleanup);
    
    MALLOC_OR_GOTO(Ao, uint32_t (*)[M_COLS], sizeof(uint32_t[N_ROWS][M_COLS]), cleanup);

    
    if(!load_trapdoor_for_t(trapdoor_path, TW)){ 
        fprintf(stderr,"[verify] Error: Failed to load trapdoor from '%s'.\n", trapdoor_path); 
        goto cleanup;
    }
    int t = TW->t;
    if(t <= 0 || t > T_PERIOD){ 
        fprintf(stderr,"[verify] Error: Invalid time t=%d found in trapdoor.\n", t); 
        goto cleanup;
    }

    
    
    
    

    
    if(!load_latest_cipher_for_t(t, cipher_path, C)){ 
        fprintf(stderr,"[verify] Error: Failed to load ciphertext for t=%d from '%s'.\n", t, cipher_path); 
        goto cleanup; 
    }
    
    
    result = verify_core(TW, C, Ao);

cleanup:
    free(TW);
    free(C);
    free(Ao);
    return result;
}

/**
 * 接口2：直接从内存数据结构验证
 */
bool verify_search_from_data(const cipher_rec_t* C, 
                             const trapdoor_rec_t* TW, 
                             const uint32_t Ao_t_unused[N_ROWS][M_COLS])
{
    fprintf(stderr, "[verify] Starting verification from data structures (no self-check)...\n");

    if (C == NULL || TW == NULL) {
        fprintf(stderr, "[verify] Error: One or more data structure pointers are NULL.\n");
        return false;
    }
    if (C->t != TW->t) {
        fprintf(stderr, "[verify] Error: Mismatch in time t between ciphertext (%d) and trapdoor (%d).\n", C->t, TW->t);
        return false;
    }

    
    (void)Ao_t_unused;

    return verify_core(TW, C, NULL);
}
