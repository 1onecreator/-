#define _POSIX_C_SOURCE 200809L 
#include "cipher.h"
#include "global_sampling.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>

#include <sys/stat.h> 
#ifdef _WIN32
#include <direct.h> 
#include <windows.h> 
#else
#include <sys/stat.h> 
#include <sys/types.h>
#include <dirent.h> 
#include <unistd.h> 
#endif

#include "generateKW.h"
#include "public_params.h"   
#include "hash.h"
#include "samplePre.h"
#include "config_params.h"


#ifndef CIPHER_VERBOSE_MODE
#define CIPHER_VERBOSE_MODE 1
#endif

#if CIPHER_VERBOSE_MODE
#define LOG_CIPHER(...) fprintf(stderr, __VA_ARGS__)
#else
#define LOG_CIPHER(...) do {} while (0)
#endif

#define MALLOC_OR_GOTO(ptr, type, size, label) \
    do { \
        ptr = (type)malloc(size); \
        if (!ptr) { \
            LOG_CIPHER("[cipher] ERROR: Memory allocation failed for %s in %s.\n", #ptr, __func__); \
            goto label; \
        } \
    } while (0)

#define CIPHER_DIR "cipher"


#ifndef DIAG_H2
#define DIAG_H2 1
#endif
#ifndef DIAG_AIT
#define DIAG_AIT 1
#endif
#ifndef DIAG_LOAD
#define DIAG_LOAD 1
#endif

#define SIGMA_X_DEFAULT 0.1
#define SIGMA_x_DEFAULT 0.1
#define GAUSS_TAIL 12.0

// 用于实时记录密文数量。-1 表示尚未初始化。
static int cipher_num = -1;

/* ---------------- 模运算 ---------------- */
static inline uint32_t mod_q_i64(long long x){
    long long r = x % (long long)Q_MODULUS;
    if (r < 0) r += Q_MODULUS;
    return (uint32_t)r;
}
static inline uint32_t mod_q_i128(__int128 x){
    long long r = (long long)(x % (__int128)Q_MODULUS);
    if (r < 0) r += Q_MODULUS;
    return (uint32_t)r;
}

/* ---------------- 轻量 PRNG ---------------- */
static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static inline uint64_t xorshift64star(void){
    uint64_t x=rng_state; x^=x>>12; x^=x<<25; x^=x>>27; rng_state=x;
    return x*0x2545F4914F6CDD1DULL;
}
static inline void rng_seed(uint64_t s){
    rng_state = s ? s : 0x9e3779b97f4a7c15ULL;
    for(int i=0;i<8;++i) (void)xorshift64star();
}

/* ---------------- 线代 ---------------- */
static void nm_mul_mm_mod(const uint32_t A_nm[N_ROWS][M_COLS],
                          const uint32_t R_mm[M_COLS][M_COLS],
                          uint32_t Out_nm[N_ROWS][M_COLS]){
    for (int i=0;i<N_ROWS;i++){
        for (int j=0;j<M_COLS;j++){
            __int128 acc=0;
            for (int k=0;k<M_COLS;k++) acc += (__int128)A_nm[i][k]*R_mm[k][j];
            Out_nm[i][j] = mod_q_i128(acc);
        }
    }
}
static void mtimesS_add(const uint32_t A_nm[N_ROWS][M_COLS],
                        const uint32_t S_nl[N_ROWS][L_COLS],
                        uint32_t Out_ml[M_COLS][L_COLS]){
    for (int j=0;j<M_COLS;j++){
        for (int c=0;c<L_COLS;c++){
            __int128 acc=0;
            for (int r=0;r<N_ROWS;r++) acc += (__int128)A_nm[r][j]*S_nl[r][c];
            Out_ml[j][c] = mod_q_i128(acc);
        }
    }
}
static void nn_mul_nm_mod(const uint32_t H_nn[N_ROWS][N_ROWS],
                          const uint32_t Bnm[N_ROWS][M_COLS],
                          uint32_t Out_nm[N_ROWS][M_COLS]){
    for (int i=0;i<N_ROWS;i++){
        for (int j=0;j<M_COLS;j++){
            __int128 acc=0;
            for (int k=0;k<N_ROWS;k++) acc += (__int128)H_nn[i][k]*Bnm[k][j];
            Out_nm[i][j] = mod_q_i128(acc);
        }
    }
}
static void RT_timesX_add_signed(const int8_t R_mm[M_COLS][M_COLS],
                                 const int32_t Xs_ml[M_COLS][L_COLS],
                                 int64_t Acc_ml[M_COLS][L_COLS]){
    for (int j=0;j<M_COLS;j++){
        for (int c=0;c<L_COLS;c++){
            long long acc = (long long)Acc_ml[j][c];
            for (int r=0;r<M_COLS;r++) acc += (long long)R_mm[r][j]*(long long)Xs_ml[r][c];
            Acc_ml[j][c] = acc;
        }
    }
}

/* ---------------- 规范化 kw ---------------- */
static void normalize_kw(const char* field, const char* value, char* out, size_t cap){
    size_t w=0;
    for(const char* p=field; *p && w+1<cap; ++p){
        unsigned char ch=(unsigned char)*p;
        if((ch>='A'&&ch<='Z')||(ch>='a'&&ch<='z')||(ch>='0'&&ch<='9')) out[w++]=(char)tolower(ch);
    }
    if(w+1<cap) out[w++]=':';
    for(const char* p=value; *p && w+1<cap; ++p){
        unsigned char ch=(unsigned char)*p;
        if((ch>='A'&&ch<='Z')||(ch>='a'&&ch<='z')||(ch>='0'&&ch<='9')) out[w++]=(char)tolower(ch);
        else if(ch==' ') { if(w+1<cap) out[w++]=' '; }
    }
    out[(w<cap)?w:(cap-1)]='\0';
}

/* ---------------- 离散高斯 ---------------- */
typedef struct { double sigma; int K; double *cdf; int ready; } dgauss_tab_t;

static void dgauss_build(double sigma, dgauss_tab_t* T){
    if (T->ready && fabs(T->sigma - sigma) < 1e-12) return;
    if (T->cdf) { free(T->cdf); T->cdf=NULL; }
    int K = (int)ceil(GAUSS_TAIL * sigma);
    if (K < 8) K = 8;
    T->cdf = (double*)malloc((size_t)(K+1)*sizeof(double));
    T->sigma = sigma; T->K = K;
    double mass = 0.0;
    for (int k=0;k<=K;k++){ double w=exp(-(double)k*(double)k/(2.0*sigma*sigma)); T->cdf[k]=w; mass+=w; }
    double acc=0.0; for(int k=0;k<=K;k++){ acc += T->cdf[k]/mass; T->cdf[k]=acc; }
    if (T->cdf[K] < 1.0) T->cdf[K] = 1.0;
    T->ready = 1;
}
static uint64_t xorshift64star_for_gauss(void){ return xorshift64star(); }
static inline double rng_uniform01_gauss(void){
    uint64_t r = xorshift64star_for_gauss(); return (double)((r >> 11) * (1.0/9007199254740992.0));
}
static inline int dgauss_sample_sym(dgauss_tab_t* T){
    double u = rng_uniform01_gauss();
    int lo=0, hi=T->K;
    while (lo<hi){ int mid=(lo+hi)/2; if (u<=T->cdf[mid]) hi=mid; else lo=mid+1; }
    int k = lo; if (k==0) return 0;
    return (xorshift64star() & 1ULL) ? k : -k;
}
static dgauss_tab_t GAUSS_X = {0};
static dgauss_tab_t GAUSS_x = {0};

/* ===================================================================== */
/*                       文件与目录操作工具                                */
/* ===================================================================== */
static void ensure_cipher_directory_exists(void) {
    struct stat st = {0};
    if (stat(CIPHER_DIR, &st) == -1) {
        LOG_CIPHER("[cipher] Directory '%s' not found. Creating...\n", CIPHER_DIR);
        #ifdef _WIN32
            if (_mkdir(CIPHER_DIR) != 0 && errno != EEXIST) {
        #else
            if (mkdir(CIPHER_DIR, 0755) != 0 && errno != EEXIST) {
        #endif
            LOG_CIPHER("[cipher] ERROR: Failed to create directory '%s'.\n", CIPHER_DIR);
            #if CIPHER_VERBOSE_MODE
            perror("[cipher] mkdir");
            #endif
        }
    }
}

/* ===================================================================== */
/*                       temp.txt 解析工具                                */
/* ===================================================================== */
static int starts_with(const char* s, const char* prefix){
    for(; *prefix; ++s, ++prefix){ if(*s != *prefix) return 0; } return 1;
}
static int parse_u32_row(const char* line, uint32_t* out, int cols){
    int n=0; const char* p=line;
    while (*p && n<cols){
        while (*p && isspace((unsigned char)*p)) ++p; if(!*p) break;
        char* endp=NULL; unsigned long v=strtoul(p,&endp,10); if (endp==p) break;
        out[n++] = (uint32_t)(v % Q_MODULUS); p=endp;
    }
    return n;
}
static int parse_i32_row(const char* line, int32_t* out, int cols){
    int n=0; const char* p=line;
    while (*p && n<cols){
        while (*p && isspace((unsigned char)*p)) ++p; if(!*p) break;
        char* endp=NULL; long v=strtol(p,&endp,10); if (endp==p) break;
        out[n++] = (int32_t)v; p=endp;
    }
    return n;
}
static int load_block_by_t_pair(const char* path, int t, const char* pair_tag,
                                const char* name, int rows, int cols,
                                uint32_t* out_u32, int is_i32)
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
                    if (is_i32){
                        if (parse_i32_row(line, ((int32_t*)out_u32)+r*cols, cols)!=cols){ fclose(fp); return 0; }
                    }else{
                        if (parse_u32_row(line, out_u32+r*cols, cols)!=cols){ fclose(fp); return 0; }
                    }
                }
                fclose(fp); return 1;
            }
        }
    }
    fclose(fp); return 0;
}
static int load_from_temp_A4(uint32_t A_out[N_ROWS][M_COLS], int t){
    return load_block_by_t_pair("temp.txt",t,"A/T_A","A",N_ROWS,M_COLS,&A_out[0][0],0);
}
static int load_from_temp_Ao4(uint32_t Ao_out[N_ROWS][M_COLS], int t){
    return load_block_by_t_pair("temp.txt",t,"Ao/T_Ao","Ao",N_ROWS,M_COLS,&Ao_out[0][0],0);
}
static int load_from_temp_TAo4(int32_t T_out[M_COLS][M_COLS], int t){
    return load_block_by_t_pair("temp.txt",t,"Ao/T_Ao","T_Ao",M_COLS,M_COLS,(uint32_t*)&T_out[0][0],1);
}
static int load_from_temp_UuInv4(uint32_t Uinv[M_COLS][M_COLS], int t){
    (void)t;
    FILE* fp=fopen("temp.txt","r"); if(!fp) return 0;
    char hdr[64]; 
    snprintf(hdr,sizeof(hdr),"Uu_inv_cache%d (", t);
    char line[8192];
    while (fgets(line,sizeof(line),fp)){
        if (strncmp(line, hdr, (int)strlen(hdr))==0){
            for (int r=0;r<M_COLS;r++){
                if (!fgets(line,sizeof(line),fp)){ fclose(fp); return 0; }
                int n=0; const char* p=line;
                while (*p && n<M_COLS){
                    while (*p && isspace((unsigned char)*p)) ++p; if(!*p) break;
                    char* endp=NULL; unsigned long v=strtoul(p,&endp,10);
                    if (endp==p) break;
                    Uinv[r][n++] = (uint32_t)(v % Q_MODULUS);
                    p=endp;
                }
                if (n!=M_COLS){ fclose(fp); return 0; }
            }
            fclose(fp); return 1;
        }
    }
    fclose(fp); return 0;
}

/* ===================================================================== */
/*                               主实现                                   */
/* ===================================================================== */

int Cipher_count_record(void) {
    LOG_CIPHER("[cipher] Counting existing records in '%s/' directory...\n", CIPHER_DIR);
    ensure_cipher_directory_exists();
    int total_records = 0;

#ifdef _WIN32
    
    char search_path[256];
    snprintf(search_path, sizeof(search_path), "%s\\cipher_*.txt", CIPHER_DIR);
    
    WIN32_FIND_DATA fd;
    HANDLE hFind = FindFirstFile(search_path, &fd);
    if (hFind == INVALID_HANDLE_VALUE) {
        LOG_CIPHER("[cipher] No cipher files found. Record count is 0.\n");
        cipher_num = 0;
        return 0;
    }
    
    do {
        char file_path[256];
        snprintf(file_path, sizeof(file_path), "%s\\%s", CIPHER_DIR, fd.cFileName);
        FILE* fp = fopen(file_path, "r");
        if (fp) {
            char line[512];
            while (fgets(line, sizeof(line), fp)) {
                
                if (strstr(line, "----")) {
                    total_records++;
                }
            }
            fclose(fp);
        } else {
            LOG_CIPHER("[cipher] WARNING: Could not open %s to count records.\n", file_path);
        }
    } while (FindNextFile(hFind, &fd) != 0);
    FindClose(hFind);

#else
    
    DIR *d = opendir(CIPHER_DIR);
    if (!d) {
        LOG_CIPHER("[cipher] ERROR: Cannot open directory '%s'. Record count is 0.\n", CIPHER_DIR);
        cipher_num = 0;
        return 0;
    }

    struct dirent *dir;
    while ((dir = readdir(d)) != NULL) {
        if (strncmp(dir->d_name, "cipher_", 7) == 0 && strstr(dir->d_name, ".txt")) {
            char file_path[256];
            snprintf(file_path, sizeof(file_path), "%s/%s", CIPHER_DIR, dir->d_name);
            FILE* fp = fopen(file_path, "r");
            if (fp) {
                char line[512];
                while(fgets(line, sizeof(line), fp)) {
                    if (strstr(line, "----")) {
                        total_records++;
                    }
                }
                fclose(fp);
            } else {
                 LOG_CIPHER("[cipher] WARNING: Could not open %s to count records.\n", file_path);
            }
        }
    }
    closedir(d);
#endif

    cipher_num = total_records;
    LOG_CIPHER("[cipher] Total records found: %d\n", cipher_num);
    return cipher_num;
}

int Cipher_get_record_count(void) {
    if (cipher_num < 0) {
        return Cipher_count_record();
    }
    return cipher_num;
}


bool Cipher_encrypt_with_fields(int t,
                                const kw_field_t* fields,
                                const char* const* values,
                                int kprime)
{
    if (t <= 0 || !fields || !values || kprime < 0 || kprime > K_LEVEL) return false;

    if (cipher_num < 0) {
        Cipher_count_record();
    }

    int record_id = cipher_num + 1;

    LOG_CIPHER("\n========== CIPHER ENCRYPTION START ==========\n");
    LOG_CIPHER("[cipher] t=%d, record_id=%d, kprime=%d\n", t, record_id, kprime);

    bool success = false;
    uint32_t (*A_loaded)[M_COLS] = NULL;
    uint32_t (*Ao_loaded)[M_COLS] = NULL;
    int32_t  (*T_Ao_loaded)[M_COLS] = NULL;
    uint32_t (*Uinv)[M_COLS] = NULL;
    uint32_t (*A_it)[N_ROWS][M_COLS] = NULL;
    uint32_t (*S)[L_COLS] = NULL;
    int32_t  (*Xs)[L_COLS] = NULL;
    uint32_t (*Cu_tmp)[L_COLS] = NULL;
    int64_t  (*Cu_acc)[L_COLS] = NULL;
    uint32_t (*C_ut)[L_COLS] = NULL;
    uint32_t (*C_i)[M_COLS][L_COLS] = NULL;
    uint32_t (*Sum)[M_COLS] = NULL;
    uint32_t (*H_i_mat)[N_ROWS] = NULL;
    uint32_t (*HB)[M_COLS] = NULL;
    uint32_t (*tmp_ci)[L_COLS] = NULL;
    int64_t  (*Ci_acc)[L_COLS] = NULL;
    int8_t   (*R_mm)[M_COLS] = NULL;
    uint32_t* features = NULL;

    MALLOC_OR_GOTO(A_loaded,    uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(Ao_loaded,   uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(T_Ao_loaded, int32_t  (*)[M_COLS],         sizeof(int32_t[M_COLS][M_COLS]),           cleanup);
    MALLOC_OR_GOTO(Uinv,        uint32_t (*)[M_COLS],         sizeof(uint32_t[M_COLS][M_COLS]),          cleanup);
    MALLOC_OR_GOTO(A_it,        uint32_t (*)[N_ROWS][M_COLS], sizeof(uint32_t[K_LEVEL][N_ROWS][M_COLS]), cleanup);
    MALLOC_OR_GOTO(S,           uint32_t (*)[L_COLS],         sizeof(uint32_t[N_ROWS][L_COLS]),          cleanup);
    MALLOC_OR_GOTO(Xs,          int32_t  (*)[L_COLS],         sizeof(int32_t[M_COLS][L_COLS]),           cleanup);
    MALLOC_OR_GOTO(Cu_tmp,      uint32_t (*)[L_COLS],         sizeof(uint32_t[M_COLS][L_COLS]),          cleanup);
    MALLOC_OR_GOTO(Cu_acc,      int64_t  (*)[L_COLS],         sizeof(int64_t[M_COLS][L_COLS]),           cleanup);
    MALLOC_OR_GOTO(C_ut,        uint32_t (*)[L_COLS],         sizeof(uint32_t[M_COLS][L_COLS]),          cleanup);
    MALLOC_OR_GOTO(C_i,         uint32_t (*)[M_COLS][L_COLS], sizeof(uint32_t[K_LEVEL][M_COLS][L_COLS]), cleanup);
    MALLOC_OR_GOTO(Sum,         uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(H_i_mat,     uint32_t (*)[N_ROWS],         sizeof(uint32_t[N_ROWS][N_ROWS]),          cleanup);
    MALLOC_OR_GOTO(HB,          uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(tmp_ci,      uint32_t (*)[L_COLS],         sizeof(uint32_t[M_COLS][L_COLS]),          cleanup);
    MALLOC_OR_GOTO(Ci_acc,      int64_t  (*)[L_COLS],         sizeof(int64_t[M_COLS][L_COLS]),           cleanup);
    MALLOC_OR_GOTO(R_mm,        int8_t   (*)[M_COLS],         sizeof(int8_t[M_COLS][M_COLS]),            cleanup);

    if (!load_from_temp_A4(A_loaded, t))   { LOG_CIPHER("[cipher] load At%d failed\n", t);   goto cleanup; }
    if (!load_from_temp_Ao4(Ao_loaded, t)) { LOG_CIPHER("[cipher] load Aot%d failed\n", t);  goto cleanup; }
    if (!load_from_temp_TAo4(T_Ao_loaded,t)){ LOG_CIPHER("[cipher] load T_Aot%d failed\n", t);goto cleanup; }
    if (!load_from_temp_UuInv4(Uinv, t))   { LOG_CIPHER("[cipher] load Uu_inv_cache failed\n"); goto cleanup; }

    #if (DIAG_LOAD && CIPHER_VERBOSE_MODE)
    LOG_CIPHER("[cipher] loaded At%d/Aot%d/T_Aot%d/Uu_inv_cache for t=%d\n", t, t, t, t);
    #endif

    uint8_t present[K_LEVEL]; memset(present, 0, sizeof(present));
    uint8_t beta_bits[(L_COLS+7)/8]; memset(beta_bits, 0, sizeof(beta_bits));
    kw_params_t P = { .l_bits = L_COLS, .k_hashes = 3 };

    for (int j=0;j<kprime;j++){
        int idx = (int)fields[j];
        if (idx>=0 && idx<K_LEVEL) present[idx]=1;
        uint8_t one[(L_COLS+7)/8];
        kw_encode_field(&P, fields[j], values[j], one);
        for (size_t b=0;b<sizeof(beta_bits);++b) beta_bits[b] |= one[b];
    }

    for (int i=0;i<K_LEVEL;i++) {
        nm_mul_mm_mod(A_I[i], Uinv, A_it[i]);
    }

    uint64_t main_seed = ((uint64_t)t<<32) ^ (uint64_t)record_id ^ 0xC1BEEF9EULL;
    global_rng_seed(main_seed);
    
    // 使用全局采样空间
    void* gauss_X_table = NULL;
    void* gauss_x_table = NULL;
    if (!get_global_gauss_X(SIGMA_X_DEFAULT, &gauss_X_table) || 
        !get_global_gauss_x(SIGMA_x_DEFAULT, &gauss_x_table)) {
        LOG_CIPHER("[cipher] Failed to get global sampling tables\n");
        return false;
    }

    for (int r=0;r<N_ROWS;r++)
        for (int c=0;c<L_COLS;c++)
            S[r][c] = global_sample_mod_q();

    for (int r=0;r<M_COLS;r++)
        for (int c=0;c<L_COLS;c++)
            Xs[r][c] = global_sampleleft_dgauss_sample(gauss_X_table);

    memset(Cu_tmp, 0, sizeof(Cu_tmp));
    mtimesS_add(A_loaded, S, Cu_tmp);
    for (int r=0;r<M_COLS;r++)
        for (int c=0;c<L_COLS;c++)
            Cu_acc[r][c] = (int64_t)Cu_tmp[r][c] + (int64_t)Xs[r][c];
    for (int r=0;r<M_COLS;r++)
        for (int c=0;c<L_COLS;c++)
            C_ut[r][c] = mod_q_i64(Cu_acc[r][c]);

    uint32_t c_beta[L_COLS];
    for (int c=0;c<L_COLS;c++){
        __int128 acc=0;
        for (int r=0;r<N_ROWS;r++) acc += (__int128)U_vec[r]*S[r][c];
        uint32_t val = mod_q_i128(acc);
        int32_t x_small = global_sampleleft_dgauss_sample(gauss_x_table);
        long long v = (long long)val + (long long)x_small;
        uint32_t bump = (beta_bits[c>>3] & (1u<<(c&7))) ? (Q_MODULUS/2u) : 0u;
        c_beta[c] = mod_q_i64(v + (long long)bump);
    }

    memset(C_i, 0, sizeof(C_i));
    const kw_field_meta_t* field_meta = kw_get_fields();
    for (int i=0;i<K_LEVEL;i++){
        uint32_t Sum[N_ROWS][M_COLS];
        if (present[i]) {
            const char* fname = field_meta[i].key;
            const char* vstr = "None";
            for (int j=0;j<kprime;j++){ if ((int)fields[j]==i){ vstr=values[j]; break; } }
            char kwbuf[1024]; normalize_kw(fname, vstr, kwbuf, sizeof(kwbuf));
            uint8_t inb[1200];
            size_t wl = strnlen(kwbuf, sizeof(kwbuf));
            memcpy(inb, kwbuf, wl);
            inb[wl+0]=(uint8_t)(t & 0xFF); inb[wl+1]=(uint8_t)((t>>8)&0xFF);
            inb[wl+2]=(uint8_t)((t>>16)&0xFF); inb[wl+3]=(uint8_t)((t>>24)&0xFF);
            uint32_t H_i_mat[N_ROWS][N_ROWS];
            H2_make(inb, (uint32_t)(wl+4), (uint32_t)(i+1), H_i_mat);
            uint32_t HB[N_ROWS][M_COLS];
            nn_mul_nm_mod(H_i_mat, B_mat, HB);
            for (int r=0;r<N_ROWS;r++)
                for (int c=0;c<M_COLS;c++)
                    Sum[r][c] = mod_q_i64((long long)A_it[i][r][c] + (long long)HB[r][c]);
        } else {
            for (int r=0;r<N_ROWS;r++) for (int c=0;c<M_COLS;c++) Sum[r][c] = A_it[i][r][c];
        }
        uint32_t tmp[M_COLS][L_COLS]; memset(tmp, 0, sizeof(tmp));
        mtimesS_add(Sum, S, tmp);
        int64_t Ci_acc[M_COLS][L_COLS];
        for (int r=0;r<M_COLS;r++) for (int c=0;c<L_COLS;c++) Ci_acc[r][c] = (int64_t)tmp[r][c];
        int8_t R_mm[M_COLS][M_COLS];
        for (int r=0;r<M_COLS;r++)
            for (int c=0;c<M_COLS;c++){
                uint64_t rand = xorshift64star() % 100;
                if (rand < 70) R_mm[r][c] = 0; else if (rand < 85) R_mm[r][c] = 1; else R_mm[r][c] = -1;
            }
        RT_timesX_add_signed(R_mm, Xs, Ci_acc);
        for (int r=0;r<M_COLS;r++) for (int c=0;c<L_COLS;c++) C_i[i][r][c] = mod_q_i64(Ci_acc[r][c]);
    }

    size_t tot = (size_t)(1 + K_LEVEL) * (size_t)M_COLS * (size_t)L_COLS;
    features = (uint32_t*)malloc(sizeof(uint32_t)*tot);
    if (!features) {
    LOG_CIPHER("[cipher] ERROR: Memory allocation failed for features.\n");
    goto cleanup; 
    }
    
    size_t p = 0;
    for (int c=0;c<L_COLS;c++) for (int r=0;r<M_COLS;r++) features[p++] = C_ut[r][c];
    for (int idx=0; idx<K_LEVEL; idx++) for (int c=0;c<L_COLS;c++) for (int r=0;r<M_COLS;r++) features[p++] = C_i[idx][r][c];

    uint32_t h_t[N_ROWS];
    H3_make(features, 1, (uint32_t)tot, beta_bits, (uint32_t)sizeof(beta_bits), h_t);

    uint32_t s_t[M_COLS];
    uint64_t seed = ((uint64_t)t<<32) ^ (uint64_t)record_id ^ 0xC1BEEF9EULL;
    if (!SamplePre_run(Ao_loaded, T_Ao_loaded, h_t, SIGMA_SEQ[t-1], seed, s_t)) return false;

    ensure_cipher_directory_exists();
    char out_path[256];
    snprintf(out_path, sizeof(out_path), "%s/cipher_%d_%d.txt", CIPHER_DIR, t, kprime);

    FILE* fp = fopen(out_path, "a");
    if (!fp) {
        LOG_CIPHER("[cipher] ERROR: Failed to open output file %s\n", out_path);
        return false;
    }
    LOG_CIPHER("[cipher] Appending ciphertext to %s\n", out_path);

    fprintf(fp, "t=%d record_id=%d\n", t, record_id);
    fprintf(fp, "c_beta={");
    for (int i=0;i<L_COLS;i++){ fprintf(fp, "%u", c_beta[i]); if (i+1<L_COLS) fputc(',', fp); }
    fputs("}\n", fp);
    fprintf(fp, "C_ut={\n");
    for (int r=0;r<M_COLS;r++){
        for (int c=0;c<L_COLS;c++){ fprintf(fp, "%u", C_ut[r][c]); if (c+1<L_COLS) fputc(',', fp); }
        fputc('\n', fp);
    }
    fputs("}\n", fp);
    for (int idx=0; idx<K_LEVEL; idx++){
        fprintf(fp, "C_%d={\n", idx+1);
        for (int r=0;r<M_COLS;r++){
            for (int c=0;c<L_COLS;c++){ fprintf(fp, "%u", C_i[idx][r][c]); if (c+1<L_COLS) fputc(',', fp); }
            fputc('\n', fp);
        }
        fputs("}\n", fp);
    }
    fprintf(fp, "h_t={");
    for (int i=0;i<N_ROWS;i++){ fprintf(fp, "%u", h_t[i]); if (i+1<N_ROWS) fputc(',', fp); }
    fputs("}\n", fp);
    fprintf(fp, "s_t={");
    for (int i=0;i<M_COLS;i++){ fprintf(fp, "%u", s_t[i]); if (i+1<M_COLS) fputc(',', fp); }
    fputs("}\n----\n", fp);

    fclose(fp);
    
    success = true;

    cipher_num++;
    
    LOG_CIPHER("[cipher] Encryption completed successfully for t=%d, record_id=%d. New record count: %d\n", t, record_id, cipher_num);
    LOG_CIPHER("========== CIPHER ENCRYPTION END ==========\n\n");
    
    cleanup:
        free(A_loaded);
        free(Ao_loaded);
        free(T_Ao_loaded);
        free(Uinv);
        free(A_it);
        free(S);
        free(Xs);
        free(Cu_tmp);
        free(Cu_acc);
        free(C_ut);
        free(C_i);
        free(Sum);
        free(H_i_mat);
        free(HB);
        free(tmp_ci);
        free(Ci_acc);
        free(R_mm);
        free(features);
    
    return success;
}



void Cipher_clear_all_files(void) {
    LOG_CIPHER("[cipher] Clearing contents of all cipher text files from '%s/' directory...\n", CIPHER_DIR);
    ensure_cipher_directory_exists();

#ifdef _WIN32

    char search_path[256];
    snprintf(search_path, sizeof(search_path), "%s\\cipher_*.txt", CIPHER_DIR);
    
    WIN32_FIND_DATA fd;
    HANDLE hFind = FindFirstFile(search_path, &fd);
    if (hFind == INVALID_HANDLE_VALUE) {
        LOG_CIPHER("[cipher] No cipher files found to clear.\n");
        cipher_num = 0; 
        return;
    }
    
    int files_cleared = 0;
    do {
        char file_path[256];
        snprintf(file_path, sizeof(file_path), "%s\\%s", CIPHER_DIR, fd.cFileName);
        FILE* fp = fopen(file_path, "w"); 
        if (fp) {
            fclose(fp);
            LOG_CIPHER("[cipher] Cleared contents of: %s\n", file_path);
            files_cleared++;
        } else {
            LOG_CIPHER("[cipher] ERROR: Failed to open and clear %s\n", file_path);
        }
    } while (FindNextFile(hFind, &fd) != 0);
    
    FindClose(hFind);
    if (files_cleared > 0) {
        LOG_CIPHER("[cipher] Total files cleared: %d\n", files_cleared);
    }

#else
    
    DIR *d = opendir(CIPHER_DIR);
    if (!d) {
        LOG_CIPHER("[cipher] ERROR: Cannot open directory '%s'.\n", CIPHER_DIR);
        cipher_num = 0; 
        return;
    }

    struct dirent *dir;
    int files_cleared = 0;
    while ((dir = readdir(d)) != NULL) {
        if (strncmp(dir->d_name, "cipher_", 7) == 0 && strstr(dir->d_name, ".txt")) {
            char file_path[256];
            snprintf(file_path, sizeof(file_path), "%s/%s", CIPHER_DIR, dir->d_name);
            FILE* fp = fopen(file_path, "w"); 
            if (fp) {
                fclose(fp);
                LOG_CIPHER("[cipher] Cleared contents of: %s\n", file_path);
                files_cleared++;
            } else {
                LOG_CIPHER("[cipher] ERROR: Failed to open and clear %s\n", file_path);
            }
        }
    }
    closedir(d);
    
    if (files_cleared == 0) {
        LOG_CIPHER("[cipher] No cipher files found to clear.\n");
    } else {
        LOG_CIPHER("[cipher] Total files cleared: %d\n", files_cleared);
    }
#endif

    cipher_num = 0;
    LOG_CIPHER("[cipher] Internal record count has been reset to 0.\n");
    LOG_CIPHER("[cipher] Clearing process finished.\n");
}
