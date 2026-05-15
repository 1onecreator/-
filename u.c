
#include "u.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef q
#undef q
#endif
#ifdef n
#undef n
#endif
#ifdef m
#undef m
#endif

/* ===== 轻量 PRG ===== */
static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static uint32_t prg_u32(uint64_t* st) {
    *st = splitmix64(*st ? *st : 0x6a09e667f3bcc909ULL);
    return (uint32_t)(*st & 0xffffffffu);
}

/* —— 均匀采样到 [0, q) —— */
static uint32_t sample_mod_q(uint64_t* st) {
    const uint64_t qmod  = (uint64_t)Q_MODULUS;     
    const uint64_t RANGE = 1ull << 32;              
    const uint64_t limit = (RANGE / qmod) * qmod;    
    for (;;) {
        uint32_t x = prg_u32(st);
        if ((uint64_t)x < limit) return (uint32_t)((uint64_t)x % qmod);
    }
}

/* ---------------- 公开 API：生成 u ---------------- */
bool U_generate(uint32_t u[N_ROWS]) {
    
    uint64_t st = 0x1234d00df00dbeefULL
                ^ ((uint64_t)N_ROWS << 32)
                ^ (uint64_t)Q_MODULUS;

    bool any_nonzero = false;
    for (int i = 0; i < N_ROWS; ++i) {
        u[i] = sample_mod_q(&st);
        if (u[i] != 0u) any_nonzero = true;
    }
    
    if (!any_nonzero) u[0] = 1u;
    return true;
}

/* ================= 文件工具 ================= */

static char* read_all(const char* path, size_t* out_len) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return NULL;
    fseek(fp, 0, SEEK_END);
    long n = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (n < 0) { fclose(fp); return NULL; }
    char* buf = (char*)malloc((size_t)n + 1);
    if (!buf) { fclose(fp); return NULL; }
    size_t rd = fread(buf, 1, (size_t)n, fp);
    fclose(fp);
    buf[rd] = '\0';
    if (out_len) *out_len = rd;
    return buf;
}

static int find_last_endif(const char* s) {
    const char* last = NULL;
    const char* p = s;
    while ((p = strstr(p, "#endif")) != NULL) {
        last = p; p += 6;
    }
    return last ? (int)(last - s) : -1;
}


static void remove_block(char** text, size_t* len, const char* token) {
    char* s = *text;
    char* pos = strstr(s, token);
    if (!pos) return;
    char* end = strstr(pos, "};");
    if (!end) return;
    end += 2; 

    size_t left_len  = (size_t)(pos - s);
    size_t right_len = *len - (size_t)(end - s);

    char* out = (char*)malloc(left_len + right_len + 1);
    if (!out) return;

    memcpy(out, s, left_len);
    memcpy(out + left_len, end, right_len);
    out[left_len + right_len] = '\0';

    free(*text);
    *text = out;
    *len  = left_len + right_len;
}

static void fprint_u(FILE* fp, const uint32_t u[N_ROWS]) {
    fprintf(fp, "const uint32_t U_vec[N_ROWS] = {");
    for (int i = 0; i < N_ROWS; ++i) {
        fprintf(fp, "%u%s", u[i], (i + 1 == N_ROWS) ? "" : ",");
    }
    fprintf(fp, "};\n\n");
}

static bool create_new_header(const char* path, const uint32_t u[N_ROWS]) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return false;
    fprintf(fp, "// public_params.h — auto-generated (U_vec)\n");
    fprintf(fp, "#ifndef PUBLIC_PARAMS_H\n#define PUBLIC_PARAMS_H\n\n");
    fprintf(fp, "#include <stdint.h>\n#include \"config_params.h\"\n\n");
    fprintf(fp, "/* 本段由 U_generate 生成：U_vec[N_ROWS] ∈ Z_q^n */\n\n");
    fprint_u(fp, u);
    fprintf(fp, "#endif // PUBLIC_PARAMS_H\n");
    fclose(fp);
    return true;
}

/* ---------------- 公开 API ---------------- */
bool U_write_public_params(const char* path, const uint32_t u[N_ROWS]) {
    size_t len = 0;
    char* old = read_all(path, &len);
    if (!old) {
        return create_new_header(path, u);
    }

    
    remove_block(&old, &len, "const uint32_t U_vec[");

    
    int pos = find_last_endif(old);
    if (pos < 0) {
        
        FILE* fp = fopen(path, "wb");
        if (!fp) { free(old); return false; }
        fwrite(old, 1, len, fp);
        fprintf(fp, "\n/* === Appended by U_generate === */\n");
        fprint_u(fp, u);
        fclose(fp);
        free(old);
        return true;
    }

    
    FILE* fp = fopen(path, "wb");
    if (!fp) { free(old); return false; }
    fwrite(old, 1, (size_t)pos, fp);
    fprintf(fp, "\n/* === Appended by U_generate === */\n");
    fprint_u(fp, u);
    fwrite(old + (size_t)pos, 1, len - (size_t)pos, fp);
    fclose(fp);
    free(old);
    return true;
}
