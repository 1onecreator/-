
#include "public.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#ifdef q
#undef q
#endif
#ifdef n
#undef n
#endif
#ifdef m
#undef m
#endif
#ifdef t
#undef t
#endif


#include "public_params.h"

/* ===== 计算序列：SIGMA_SEQ / DELTA_SEQ =====
 * σ_pre  ≥ ||T|| * OMEGA_M * sqrt(log m)
 * σ_left ≥ ||T|| * OMEGA_M * sqrt(log (2m))  (m1=m)
 * SIGMA_SEQ[t] = max(σ_pre, σ_left)
 * DELTA_SEQ[t] = DELTA_R（若无定义则用 sqrt(n ln q)*OMEGA_M）
 */
bool Public_compute_sequences(double sigma_seq[T_PERIOD],
                              double delta_seq[T_PERIOD]) {
    // 计算 ||T||_∞（行和无穷范数）
    double norm_T = 0.0;
    for (int i = 0; i < M_COLS; ++i) {
        long long s = 0;
        for (int j = 0; j < M_COLS; ++j) {
            long long v = T_A[i][j];            
            s += (v >= 0 ? v : -v);
        }
        if ((double)s > norm_T) norm_T = (double)s;
    }
    if (norm_T < 1.0) norm_T = 1.0;

    const double omega  = OMEGA_M;
    const double log_m  = log((double)M_COLS);
    const double log_2m = log((double)(2 * M_COLS));

    const double sigma_pre   = norm_T * omega * sqrt(log_m  > 1e-9 ? log_m  : 1.0);
    const double sigma_left  = norm_T * omega * sqrt(log_2m > 1e-9 ? log_2m : 1.0);
    const double sigma_final = (sigma_pre > sigma_left) ? sigma_pre : sigma_left;

    for (int tt = 0; tt < T_PERIOD; ++tt) {      
        sigma_seq[tt] = sigma_final;
#ifdef DELTA_R
        delta_seq[tt] = (double)DELTA_R;
#else
        delta_seq[tt] = sqrt((double)N_ROWS * log((double)Q_MODULUS)) * omega;
#endif
    }
    return true;
}

/* ================= 文件写入工具（只处理 SIGMA_SEQ / DELTA_SEQ） ================= */

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
    const char* last = NULL; const char* p = s;
    while ((p = strstr(p, "#endif")) != NULL) { last = p; p += 6; }
    return last ? (int)(last - s) : -1;
}
static void remove_block(char** text, size_t* len, const char* token) {
    char* s = *text;
    for (;;) {
        char* pos = strstr(s, token);
        if (!pos) break;
        char* end = strstr(pos, "};");
        if (!end) break;
        end += 2;
        size_t left  = (size_t)(pos - s);
        size_t right = *len - (size_t)(end - s);
        char* out = (char*)malloc(left + right + 1);
        if (!out) break;
        memcpy(out, s, left);
        memcpy(out + left, end, right);
        out[left + right] = '\0';
        free(*text);
        *text = s = out;
        *len  = left + right;
    }
}


static void fprint_sigma_seq(FILE* fp, const double sigma_seq[T_PERIOD]) {
    fprintf(fp, "static const double SIGMA_SEQ[T_PERIOD] = {");
    for (int i = 0; i < T_PERIOD; ++i) {
        fprintf(fp, "%.10f%s", sigma_seq[i], (i + 1 == T_PERIOD) ? "" : ",");
    }
    fprintf(fp, "};\n\n");
}
static void fprint_delta_seq(FILE* fp, const double delta_seq[T_PERIOD]) {
    fprintf(fp, "static const double DELTA_SEQ[T_PERIOD] = {");
    for (int i = 0; i < T_PERIOD; ++i) {
        fprintf(fp, "%.10f%s", delta_seq[i], (i + 1 == T_PERIOD) ? "" : ",");
    }
    fprintf(fp, "};\n\n");
}


static bool create_new_header(const char* path,
                              const double sigma_seq[T_PERIOD],
                              const double delta_seq[T_PERIOD]) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return false;
    fprintf(fp, "// public_params.h — auto-generated (SIGMA/DELTA sequences)\n");
    fprintf(fp, "#ifndef PUBLIC_PARAMS_H\n#define PUBLIC_PARAMS_H\n\n");
    fprintf(fp, "#include <stdint.h>\n#include \"config_params.h\"\n\n");
    fprint_sigma_seq(fp, sigma_seq);
    fprint_delta_seq(fp, delta_seq);
    fprintf(fp, "#endif // PUBLIC_PARAMS_H\n");
    fclose(fp);
    return true;
}

bool Public_write_sequences(const char* path,
                            const double sigma_seq[T_PERIOD],
                            const double delta_seq[T_PERIOD]) {
    size_t len = 0;
    char* old = read_all(path, &len);
    if (!old) {
        return create_new_header(path, sigma_seq, delta_seq);
    }

    
    remove_block(&old, &len, "static const double SIGMA_SEQ");
    remove_block(&old, &len, "static const double DELTA_SEQ");

    int pos = find_last_endif(old);
    FILE* fp = fopen(path, "wb");
    if (!fp) { free(old); return false; }

    if (pos < 0) {
        fwrite(old, 1, len, fp);
        fprintf(fp, "\n/* === Appended by Public_write_sequences === */\n");
        fprint_sigma_seq(fp, sigma_seq);
        fprint_delta_seq(fp, delta_seq);
        fprintf(fp, "\n#endif // PUBLIC_PARAMS_H\n");
    } else {
        fwrite(old, 1, (size_t)pos, fp);
        fprintf(fp, "\n/* === Appended by Public_write_sequences === */\n");
        fprint_sigma_seq(fp, sigma_seq);
        fprint_delta_seq(fp, delta_seq);
        fwrite(old + (size_t)pos, 1, len - (size_t)pos, fp);
    }
    fclose(fp);
    free(old);
    return true;
}
