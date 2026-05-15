



#define _POSIX_C_SOURCE 200809L 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>

#include "Trapdoor.h"
#include "config_params.h"
#include "public_params.h"   
#include "hash.h"            
#include "sampleLeft.h"      
#include "generateKW.h"      


#ifdef _WIN32
  #include <windows.h>
  #include <bcrypt.h>
  #pragma comment(lib, "bcrypt.lib")
#else
  #include <fcntl.h>
  #include <unistd.h>
#endif

/* ===== 内存分配辅助宏  ===== */
#define MALLOC_OR_GOTO(ptr, type, size, label) \
    do { \
        ptr = (type)malloc(size); \
        if (!ptr) { \
            fprintf(stderr, "[trapdoor] ERROR: Memory allocation failed for %s.\n", #ptr); \
            goto label; \
        } \
    } while (0)

/* ---------------- 内部工具 ---------------- */
static inline long long mod_q_i64(long long x){
    long long r=x%(long long)Q_MODULUS; if(r<0) r+=Q_MODULUS; return r;
}
static inline uint32_t mod_q_i128(__int128 x){
    long long r=(long long)(x%(__int128)Q_MODULUS); if(r<0) r+=Q_MODULUS; return (uint32_t)r;
}

/* ---------------- 安全随机数：secure_random_bytes ---------------- */
static int secure_random_bytes(void* buf, size_t len){
#ifdef _WIN32
    NTSTATUS st = BCryptGenRandom(NULL, (PUCHAR)buf, (ULONG)len, BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    return (st == 0) ? 1 : 0;
#else
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) return 0;
    uint8_t* p = (uint8_t*)buf;
    size_t got = 0;
    while (got < len){
        ssize_t r = read(fd, p + got, len - got);
        if (r <= 0){ close(fd); return 0; }
        got += (size_t)r;
    }
    close(fd);
    return 1;
#endif
}

/* ---------------- 64-bit FNV-1a 哈希（作 seed 派生） ---------------- */
static inline uint64_t fnv1a64_update(uint64_t h, const void* data, size_t n){
    const uint8_t* p = (const uint8_t*)data;
    while (n--){
        h ^= *p++;
        h *= 1099511628211ull; 
    }
    return h;
}
static inline uint64_t derive_seed_from_ctx(int t, int i, const uint8_t* kw, size_t kwlen,
                                            const uint8_t nonce[16]){
    uint64_t h = 1469598103934665603ull; 
    h = fnv1a64_update(h, &t, sizeof(t));
    h = fnv1a64_update(h, &i, sizeof(i));
    h = fnv1a64_update(h, kw, kwlen);
    h = fnv1a64_update(h, nonce, 16);
    
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

/* ---------------- 关键字标准化 ---------------- */
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

/* ---------------- 从 temp.txt 读取辅助 ---------------- */
static int starts_with(const char* s, const char* p){
    for(; *p; ++s, ++p) if (*s != *p) return 0;
    return 1;
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
static int parse_i32_row(const char* line, int32_t* out, int count){
    const char* p = line;
    for(int i=0; i<count; ++i){
        while(*p && (*p==' '||*p=='\t'||*p=='\n' ||*p=='\r')) ++p;
        if(!*p) return i;
        long x; if(sscanf(p,"%ld",&x)!=1) return i;
        out[i]=(int32_t)x;
        while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') ++p;
    }
    return count;
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

/* 加载 A_ut 和 T_ut（用于 SampleLeft）- 读取 At{t} 和 T_At{t} */
int load_Aut_from_temp(int t, const char* path, uint32_t A[N_ROWS][M_COLS], int32_t T_A[M_COLS][M_COLS]){
    if (!load_block_by_t_pair(path, t, "A/T_A", "A",   N_ROWS, M_COLS, (uint32_t*)A, 0)) return 0;
    if (!load_block_by_t_pair(path, t, "A/T_A", "T_A", M_COLS, M_COLS, (uint32_t*)T_A, 1)) return 0;
    return 1;
}

/* 加载 Uu_inv_cache4 并计算 A_t = A_I · U_inv */
static int load_At_from_temp(int t, const char* path, uint32_t A_t[K_LEVEL][N_ROWS][M_COLS]){
    
    
    uint32_t (*U_inv)[M_COLS] = malloc(sizeof(uint32_t[M_COLS][M_COLS]));
    if (!U_inv) {
        fprintf(stderr, "[trapdoor] Memory allocation failed for U_inv.\n");
        return 0;
    }

    FILE* fp=fopen(path,"r");
    if(!fp) {
        free(U_inv);
        return 0;
    }

    char hdr[64]; 
    snprintf(hdr,sizeof(hdr),"Uu_inv_cache%d (", t);
    char line[16384]; 
    int found = 0;
    while (fgets(line,sizeof(line),fp)){
        if (strncmp(line, hdr, (int)strlen(hdr))==0){
            for (int r=0;r<M_COLS;r++){
                if (!fgets(line,sizeof(line),fp)){ fclose(fp); return 0; }
                int n=0; const char* p=line;
                while (*p && n<M_COLS){
                    while (*p && (*p==' '||*p=='\t'||*p=='\n'||*p=='\r')) ++p;
                    if(!*p) break;
                    char* endp=NULL; unsigned long v=strtoul(p,&endp,10);
                    if (endp==p) break;
                    U_inv[r][n++] = (uint32_t)(v % Q_MODULUS);
                    p=endp;
                }
                if (n!=M_COLS){ fclose(fp); return 0; }
            }
            found = 1;
            break;
        }
    }
    fclose(fp);
    if (!found) {
        free(U_inv);
        return 0;
    }

    for (int i=0; i<K_LEVEL; i++){
        for (int r=0; r<N_ROWS; r++){
            for (int c=0; c<M_COLS; c++){
                __int128 acc = 0;
                for (int k=0; k<M_COLS; k++){
                    acc += (__int128)A_I[i][r][k] * U_inv[k][c];
                }
                A_t[i][r][c] = mod_q_i128(acc);
            }
        }
        uint64_t fph=0;
        for(int r=0;r<N_ROWS;r++)
            for(int c=0;c<M_COLS;c++)
                fph = fph*11400714819323198485ull ^ A_t[i][r][c];
        fprintf(stderr,"[Ait-trapdoor] i=%d t=%d fp=%llu\n", i+1, t, (unsigned long long)fph);
    }

    
    free(U_inv);
    return 1;
}

/* ---------------- 保存陷门到文件 ---------------- */
bool Trapdoor_save(const TrapdoorTW* tw, const char* path){
    if (!tw || !path) return false;
    FILE* fp = fopen(path, "w");
    if (!fp) return false;

    fprintf(fp, "t=%d\n", tw->t);
    for (int i = 0; i < K_LEVEL; i++) {
        fprintf(fp, "i=%d present=%d\n", i+1, tw->present[i]);
        
        fprintf(fp, "e_%d={", i+1);
        for (int j = 0; j < 2 * M_COLS; j++) {
            if (j > 0) fprintf(fp, ",");
            fprintf(fp, "%u", tw->e[i][j]);
        }
        fprintf(fp, "}\n");
        
        fprintf(fp, "nonce_%d=", i+1);
        for (int k=0; k<16; k++) fprintf(fp, "%02x", tw->nonce[i][k]);
        fputc('\n', fp);
    }
    fclose(fp);
    return true;
}

/* ---------------- 从文件读取陷门 ---------------- */
bool load_trapdoor_for_t(const char* path, trapdoor_rec_t* out){
    if (!path || !out) return false;
    FILE* fp = fopen(path, "r");
    if (!fp) return false;

    memset(out, 0, sizeof(*out));
    char line[8192];
    int found_t = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "t=", 2) == 0) {
            sscanf(line, "t=%d", &out->t);
            found_t = 1;
        } else if (found_t && strncmp(line, "i=", 2) == 0) {
            int i, present;
            if (sscanf(line, "i=%d present=%d", &i, &present) == 2 && i >= 1 && i <= K_LEVEL) {
                out->present[i-1] = (uint8_t)present;
            }
        } else if (found_t && strncmp(line, "e_", 2) == 0) {
            int i;
            if (sscanf(line, "e_%d={", &i) == 1 && i >= 1 && i <= K_LEVEL) {
                char* start = strchr(line, '{');
                if (start) {
                    start++;
                    char* end = strchr(start, '}');
                    if (end) *end = '\0';
                    int count = 0; char* token = strtok(start, ",");
                    while (token && count < 2 * M_COLS) {
                        sscanf(token, "%u", &out->e[i-1][count]);
                        count++;
                        token = strtok(NULL, ",");
                    }
                }
            }
        }
        
    }
    fclose(fp);
    return found_t;
}

/* ---------------- 主实现 ---------------- */
bool Trapdoor_make(int t,
                   const int* fields,
                   const char* const* values,
                   int kprime,
                   TrapdoorTW* out)
{
    if (!out || t<=0) return false;
    memset(out, 0, sizeof(*out));
    out->t      = t;
    out->kprime = (unsigned)kprime;

    
    bool success = false;
    uint32_t (*A_ut)[M_COLS] = NULL;
    int32_t  (*T_ut)[M_COLS] = NULL;
    uint32_t (*A_t)[N_ROWS][M_COLS] = NULL;
    uint32_t (*H_i)[N_ROWS] = NULL;
    uint32_t (*HB)[M_COLS] = NULL;
    uint32_t (*G_i)[M_COLS] = NULL;

    
    MALLOC_OR_GOTO(A_ut,    uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(T_ut,    int32_t  (*)[M_COLS],         sizeof(int32_t[M_COLS][M_COLS]),           cleanup);
    MALLOC_OR_GOTO(A_t,     uint32_t (*)[N_ROWS][M_COLS], sizeof(uint32_t[K_LEVEL][N_ROWS][M_COLS]), cleanup);
    MALLOC_OR_GOTO(H_i,     uint32_t (*)[N_ROWS],         sizeof(uint32_t[N_ROWS][N_ROWS]),          cleanup);
    MALLOC_OR_GOTO(HB,      uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);
    MALLOC_OR_GOTO(G_i,     uint32_t (*)[M_COLS],         sizeof(uint32_t[N_ROWS][M_COLS]),         cleanup);


    
    
    
    if (!load_Aut_from_temp(t, "temp.txt", A_ut, T_ut)) {
        fprintf(stderr, "[trapdoor] load A_ut/T_ut failed\n");
        goto cleanup;
    }

    
    
    if (!load_At_from_temp(t, "temp.txt", A_t)) {
        fprintf(stderr, "[trapdoor] load A_t failed\n");
        goto cleanup;
    }

    
    const kw_field_meta_t* field_meta = kw_get_fields();

    
    uint32_t u_parts[K_LEVEL][N_ROWS];
    memset(u_parts, 0, sizeof(u_parts));
    if (kprime > 0) {
        int present_indices[K_LEVEL]; int present_count = 0;
        for (int j = 0; j < kprime; j++) {
            if (fields[j] >= 0 && fields[j] < K_LEVEL) {
                present_indices[present_count++] = fields[j];
            }
        }
        if (present_count > 0) {
            for (int r = 0; r < N_ROWS; r++) {
                uint32_t base = U_vec[r] / present_count;
                uint32_t rem  = U_vec[r] % present_count;
                for (int j = 0; j < present_count; j++) {
                    u_parts[present_indices[j]][r] = base + (j < rem ? 1 : 0);
                }
            }
        }
    }

    
    for (int i = 0; i < K_LEVEL; i++) {
        out->present[i] = 0;

        int found = 0, field_idx = -1;
        for (int j = 0; j < kprime; j++) {
            if (fields[j] == i) { found = 1; field_idx = j; break; }
        }
        if (!found) {
            memset(out->e[i], 0, sizeof(out->e[i]));
            memset(out->nonce[i], 0, sizeof(out->nonce[i])); 
            continue;
        }

        out->present[i] = 1;

        
        const char* fname = field_meta[i].key;
        const char* vstr  = values[field_idx];
        char kwbuf[1024];
        normalize_kw(fname, vstr, kwbuf, sizeof(kwbuf));

        uint8_t inb[1200];
        size_t wl = strnlen(kwbuf, sizeof(kwbuf));
        memcpy(inb, kwbuf, wl);
        inb[wl+0] = (uint8_t)(t & 0xFF);
        inb[wl+1] = (uint8_t)((t>>8) & 0xFF);
        inb[wl+2] = (uint8_t)((t>>16)& 0xFF);
        inb[wl+3] = (uint8_t)((t>>24)& 0xFF);

        
        
        fprintf(stderr,"[H2-trapdoor] i=%d t=%d dom=%u len=%u kwbuf=\"%s\"\n",
                i+1, t, (unsigned)(i+1), (unsigned)(wl+4), kwbuf);
        H2_make(inb, (uint32_t)(wl+4), (uint32_t)(i+1), H_i);

        
        for (int r = 0; r < N_ROWS; r++) {
            for (int c = 0; c < M_COLS; c++) {
                __int128 acc = 0;
                for (int k = 0; k < N_ROWS; k++) acc += (__int128)H_i[r][k] * B_mat[k][c];
                HB[r][c] = mod_q_i128(acc);
            }
        }

        
        for (int r = 0; r < N_ROWS; r++) {
            for (int c = 0; c < M_COLS; c++) {
                G_i[r][c] = mod_q_i64((long long)A_t[i][r][c] + (long long)HB[r][c]);
            }
        }

        
        if (!secure_random_bytes(out->nonce[i], 16)){
            fprintf(stderr, "[trapdoor] secure_random_bytes failed on i=%d\n", i+1);
            goto cleanup;
        }
        uint64_t seed = derive_seed_from_ctx(t, i, (const uint8_t*)kwbuf, wl, out->nonce[i]);

        if (!SampleLeft_run(A_ut, T_ut, G_i, u_parts[i], SIGMA_S, seed, out->e[i])){
            fprintf(stderr, "[trapdoor] SampleLeft_run failed on i=%d\n", i+1);
            goto cleanup;
        }
        if (!SampleLeft_verify(A_ut, G_i, out->e[i], u_parts[i])) {
            fprintf(stderr, "[trapdoor] SampleLeft_verify failed on i=%d\n", i+1);
            goto cleanup;
        }

        fprintf(stderr,"[trapdoor] e[%d] generated: eA[0..3]=%u,%u,%u,%u eG[0..3]=%u,%u,%u,%u\n",
                i, out->e[i][0], out->e[i][1], out->e[i][2], out->e[i][3],
                out->e[i][M_COLS], out->e[i][M_COLS+1], out->e[i][M_COLS+2], out->e[i][M_COLS+3]);
        
        
    }

    success = true; 

    
cleanup:
    free(A_ut);
    free(T_ut);
    free(A_t);
    free(H_i);
    free(HB);
    free(G_i);
    
    return success;
}
