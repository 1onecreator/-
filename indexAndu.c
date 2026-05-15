
#include "indexAndu.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}
static void prg_fill_u32(uint32_t* out, size_t count, uint64_t seed) {
    uint64_t s = seed ? seed : 0x6a09e667f3bcc909ULL;
    for (size_t i = 0; i < count; ++i) {
        s = splitmix64(s);
        out[i] = (uint32_t)(s & 0xffffffffu);
    }
}

/* ===== 生成 A_i 和 B（均匀模 q） ===== */
bool IndexAndU_generate(uint32_t A_all[K_FIELDS][N_ROWS][M_COLS],
                        uint32_t B[N_ROWS][M_COLS]) {
    if (K_FIELDS <= 0) return false;

    uint64_t seed_base = 0x1234cafebabedeadULL
                       ^ ((uint64_t)K_FIELDS << 48)
                       ^ ((uint64_t)N_ROWS  << 32)
                       ^ ((uint64_t)M_COLS  << 16)
                       ^ (uint64_t)Q_MODULUS;

    
    for (int idx = 0; idx < K_FIELDS; ++idx) {
        uint64_t seed = splitmix64(seed_base ^ (uint64_t)(0x1000 + idx));
        uint32_t buf[N_ROWS * M_COLS];
        prg_fill_u32(buf, (size_t)N_ROWS * (size_t)M_COLS, seed);
        size_t p = 0;
        for (int r = 0; r < N_ROWS; ++r) {
            for (int c = 0; c < M_COLS; ++c) {
                A_all[idx][r][c] = (uint32_t)(buf[p++] % Q_MODULUS);
            }
        }
    }

    
    {
        uint64_t seed = splitmix64(seed_base ^ 0xB00B1E5ULL);
        uint32_t buf[N_ROWS * M_COLS];
        prg_fill_u32(buf, (size_t)N_ROWS * (size_t)M_COLS, seed);
        size_t p = 0;
        for (int r = 0; r < N_ROWS; ++r) {
            for (int c = 0; c < M_COLS; ++c) {
                B[r][c] = (uint32_t)(buf[p++] % Q_MODULUS);
            }
        }
    }
    return true;
}

/* ===== 文件工具 ===== */
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


static void remove_block(char** text, size_t* len, const char* start_token) {
    char* s = *text;
    char* pos = strstr(s, start_token);
    if (!pos) return;

    char* end = strstr(pos, "};");
    if (!end) return;
    end += 2; 

    size_t left_len = (size_t)(pos - s);
    size_t right_len = *len - (size_t)(end - s);

    char* out = (char*)malloc(left_len + right_len + 1);
    if (!out) return;

    memcpy(out, s, left_len);
    memcpy(out + left_len, end, right_len);
    out[left_len + right_len] = '\0';

    free(*text);
    *text = out;
    *len = left_len + right_len;
}

static void fprint_matrix_B(FILE* fp, const uint32_t B[N_ROWS][M_COLS]) {
    fprintf(fp, "const uint32_t B_mat[N_ROWS][M_COLS] = {\n");
    for (int i = 0; i < N_ROWS; ++i) {
        fprintf(fp, "  {");
        for (int j = 0; j < M_COLS; ++j) {
            fprintf(fp, "%u%s", B[i][j], (j + 1 == M_COLS) ? "" : ",");
        }
        fprintf(fp, "}%s\n", (i + 1 == N_ROWS) ? "" : ",");
    }
    fprintf(fp, "};\n\n");
}
static void fprint_tensor_Ai(FILE* fp, const uint32_t A_all[K_FIELDS][N_ROWS][M_COLS]) {
    fprintf(fp, "const uint32_t A_I[K_FIELDS][N_ROWS][M_COLS] = {\n");
    for (int k = 0; k < K_FIELDS; ++k) {
        fprintf(fp, "  {\n");
        for (int i = 0; i < N_ROWS; ++i) {
            fprintf(fp, "    {");
            for (int j = 0; j < M_COLS; ++j) {
                fprintf(fp, "%u%s", A_all[k][i][j], (j + 1 == M_COLS) ? "" : ",");
            }
            fprintf(fp, "}%s\n", (i + 1 == N_ROWS) ? "" : ",");
        }
        fprintf(fp, "  }%s\n", (k + 1 == K_FIELDS) ? "" : ",");
    }
    fprintf(fp, "};\n\n");
}

static bool create_new_header(const char* path,
                              const uint32_t A_all[K_FIELDS][N_ROWS][M_COLS],
                              const uint32_t B[N_ROWS][M_COLS]) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return false;
    fprintf(fp, "/* 本段由 IndexAndU 生成：\n"
                " * - A_I[K_FIELDS][N_ROWS][M_COLS]  (K_FIELDS = KEYWORD_FIELD_COUNT)\n"
                " * - B_mat[N_ROWS][M_COLS]\n"
                " */\n\n");
    fprint_tensor_Ai(fp, A_all);
    fprint_matrix_B(fp, B);
    fclose(fp);
    return true;
}

bool IndexAndU_write_public_params(const char* path,
                                   const uint32_t A_all[K_FIELDS][N_ROWS][M_COLS],
                                   const uint32_t B[N_ROWS][M_COLS]) {
    size_t len = 0;
    char* old = read_all(path, &len);
    if (!old) {
        return create_new_header(path, A_all, B);
    }

    
    remove_block(&old, &len, "const uint32_t A_I[");
    remove_block(&old, &len, "const uint32_t B_mat[");

    
    int pos = find_last_endif(old);
    if (pos < 0) {
        
        FILE* fp = fopen(path, "wb");
        if (!fp) { free(old); return false; }
        fwrite(old, 1, len, fp);
        fprintf(fp, "\n/* === Appended by IndexAndU === */\n");
        fprint_tensor_Ai(fp, A_all);
        fprint_matrix_B(fp, B);
        fprintf(fp, "\n#endif // PUBLIC_PARAMS_H\n");
        fclose(fp);
        free(old);
        return true;
    }

    
    FILE* fp = fopen(path, "wb");
    if (!fp) { free(old); return false; }

    fwrite(old, 1, (size_t)pos, fp);
    fprintf(fp, "\n/* === Appended by IndexAndU === */\n");
    fprint_tensor_Ai(fp, A_all);
    fprint_matrix_B(fp, B);
    fwrite(old + (size_t)pos, 1, len - (size_t)pos, fp);

    fclose(fp);
    free(old);
    return true;
}
