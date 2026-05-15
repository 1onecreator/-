
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "config_params.h"
#include "public_params.h"
#include "hash.h"      
#include "NewBasis.h"  


#ifndef START_I
#define START_I 0
#endif

/* ===== 小工具 ===== */
static inline uint32_t mod_q_i128(__int128 x){
    long long r = (long long)(x % (__int128)Q_MODULUS);
    if (r < 0) r += Q_MODULUS;
    return (uint32_t)r;
}
static inline int64_t floormod_i64(int64_t x, int64_t q){
    int64_t r = x % q; return (r < 0) ? (r + q) : r;
}

/* ===== 线代：u32 (mod q) ===== */
static void mm_mul_mod_u32(const uint32_t A[M_COLS][M_COLS],
                           const uint32_t B[M_COLS][M_COLS],
                           uint32_t C[M_COLS][M_COLS]){
    for (int i=0;i<M_COLS;++i){
        for (int j=0;j<M_COLS;++j){
            __int128 acc = 0;
            for (int k=0;k<M_COLS;++k) acc += (__int128)A[i][k] * B[k][j];
            C[i][j] = mod_q_i128(acc);
        }
    }
}
static void nm_mul_mm_mod_u32(const uint32_t A_nm[N_ROWS][M_COLS],
                              const uint32_t R_mm[M_COLS][M_COLS],
                              uint32_t Out_nm[N_ROWS][M_COLS]){
    for (int i=0;i<N_ROWS;++i){
        for (int j=0;j<M_COLS;++j){
            __int128 acc = 0;
            for (int k=0;k<M_COLS;++k) acc += (__int128)A_nm[i][k] * R_mm[k][j];
            Out_nm[i][j] = mod_q_i128(acc);
        }
    }
}

/* ===== 把 H1 输出整形成上三角单位块：R=[I X;0 I], R^{-1}=[I -X;0 I] ===== */
static void make_block_unit_upper_from_hash(const uint32_t raw[M_COLS][M_COLS],
                                            uint32_t R[M_COLS][M_COLS],
                                            uint32_t Rinv[M_COLS][M_COLS]){
    const int n = N_ROWS;
    const int m = M_COLS;
    const int k = m - n;

    for (int i=0;i<m;++i) for (int j=0;j<m;++j){ R[i][j]=0; Rinv[i][j]=0; }
    for (int i=0;i<n;++i){ R[i][i]=1; Rinv[i][i]=1; }
    for (int i=0;i<k;++i){ R[n+i][n+i]=1; Rinv[n+i][n+i]=1; }

    for (int i=0;i<n;++i){
        for (int j=0;j<k;++j){
            uint32_t v = raw[i][n+j] % 5u;          
            R[i][n+j]    = v;
            Rinv[i][n+j] = v ? (uint32_t)(Q_MODULUS - v) : 0u;
        }
    }
}

/* ===== 文本写出：带 t 的区块名 ===== */
static void write_matrix_u32(FILE* fp, const char* name, int rows, int cols, const uint32_t *mat){
    fprintf(fp, "%s (%dx%d):\n", name, rows, cols);
    for (int i=0;i<rows;++i){
        for (int j=0;j<cols;++j){
            fprintf(fp, "%u%s", mat[i*cols + j], (j+1==cols) ? "" : " ");
        }
        fputc('\n', fp);
    }
}
static void write_matrix_i32(FILE* fp, const char* name, int rows, int cols, const int32_t *mat){
    fprintf(fp, "%s (%dx%d):\n", name, rows, cols);
    for (int i=0;i<rows;++i){
        for (int j=0;j<cols;++j){
            fprintf(fp, "%d%s", mat[i*cols + j], (j+1==cols) ? "" : " ");
        }
        fputc('\n', fp);
    }
}

/* 结果落盘（只追加）：块名随 t 改为
   - A/T_A:      At<t>,   T_At<t>,   Uu_inv_cache<t>
   - Ao/T_Ao:    Aot<t>,  T_Aot<t>,  Uo_inv_cache<t>
*/
static void dump_result_to_temp_txt(int t, const char* tag,
                                    const uint32_t A_upd[N_ROWS][M_COLS],
                                    const int32_t  T_upd[M_COLS][M_COLS],
                                    const uint32_t Uinv_cache[M_COLS][M_COLS])
{
    FILE* fp = fopen("temp.txt", "a"); 
    if (!fp){ perror("[renew] fopen temp.txt"); return; }
    fprintf(fp, "==============================\n");
    fprintf(fp, "t = %d  |  Pair = %s\n", t, tag);
    fprintf(fp, "Q = %u, N_ROWS = %d, M_COLS = %d\n\n", (unsigned)Q_MODULUS, N_ROWS, M_COLS);

    char nameA[64], nameTA[64], nameU[64];

    if (strcmp(tag, "A/T_A") == 0){
        snprintf(nameA,  sizeof(nameA),  "At%d", t);
        snprintf(nameTA, sizeof(nameTA), "T_At%d", t);
        snprintf(nameU,  sizeof(nameU),  "Uu_inv_cache%d", t);
    }else{ 
        snprintf(nameA,  sizeof(nameA),  "Aot%d", t);
        snprintf(nameTA, sizeof(nameTA), "T_Aot%d", t);
        snprintf(nameU,  sizeof(nameU),  "Uo_inv_cache%d", t);
    }

    write_matrix_u32(fp, nameA,  N_ROWS, M_COLS, &A_upd[0][0]);
    write_matrix_i32(fp, nameTA, M_COLS, M_COLS, &T_upd[0][0]);
    write_matrix_u32(fp, nameU,  M_COLS, M_COLS, &Uinv_cache[0][0]);

    fputc('\n', fp);
    fclose(fp);
}

/* ===== u32<->i64 的拷贝 ===== */
static void copy_u32_to_i64_mat(const uint32_t* src, int rows, int cols, int64_t* dst){
    for (int i=0;i<rows*cols;++i) dst[i] = (int64_t)src[i];
}
static void copy_i64_to_u32_mat_modq(const int64_t* src, int rows, int cols, uint32_t* dst, int64_t q){
    for (int i=0;i<rows*cols;++i){
        int64_t z = floormod_i64(src[i], q);
        dst[i] = (uint32_t)z;
    }
}

/* ===== 核心：构造 U 与 Uinv_chain，并用 NewBasisDel_C 计算 (A', T_A')；最后落盘 =====
   - U = ∏_{s=i}^{t-1} R_s
   - Uinv_chain = ∏_{s=i}^{t-1} R_s^{-1}   （按 s 递增次序右乘，保证 A' = A * Uinv_chain）
   - NewBasisDel_C(A, U, T_A, ...) 内部执行 B=A*U^{-1} 与 T'=U*T_A（并中心化/LLL）
*/
static int do_update_and_dump(const uint32_t A_in[N_ROWS][M_COLS],
                              const int32_t  T_in[M_COLS][M_COLS],
                              int t,
                              const char* tag_for_dump)
{
    /* (1) 链式累计 U 与 Uinv_chain，并同步演进 A_cur 方便 H1(A|s) */
    uint32_t U[M_COLS][M_COLS] = {0};
    uint32_t Uinv_chain[M_COLS][M_COLS] = {0};
    for (int i=0;i<M_COLS;++i){ U[i][i] = 1; Uinv_chain[i][i] = 1; }

    uint32_t A_cur[N_ROWS][M_COLS];
    memcpy(A_cur, A_in, sizeof(A_cur));

    for (int s = START_I; s < t; ++s){
        uint32_t raw[M_COLS][M_COLS];
        H1_make(A_cur, (uint32_t)s, raw);

        uint32_t R[M_COLS][M_COLS], Rinv[M_COLS][M_COLS];
        make_block_unit_upper_from_hash(raw, R, Rinv);

        
        uint32_t U_new[M_COLS][M_COLS];
        mm_mul_mod_u32(R, U, U_new);
        memcpy(U, U_new, sizeof(U));

        
        uint32_t Uinv_new[M_COLS][M_COLS];
        mm_mul_mod_u32(Uinv_chain, Rinv, Uinv_new);
        memcpy(Uinv_chain, Uinv_new, sizeof(Uinv_chain));

        
        uint32_t A_new[N_ROWS][M_COLS];
        nm_mul_mm_mod_u32(A_cur, Rinv, A_new);
        memcpy(A_cur, A_new, sizeof(A_cur));
    }

    /* (2) 调用 NewBasisDel_C：得到 A' 与 T_A' （i64 工作区） */
    int64_t *A_i64  = (int64_t*)malloc(sizeof(int64_t)*(size_t)N_ROWS*M_COLS);
    int64_t *U_i64  = (int64_t*)malloc(sizeof(int64_t)*(size_t)M_COLS*M_COLS);
    int64_t *T_i64  = (int64_t*)malloc(sizeof(int64_t)*(size_t)M_COLS*M_COLS);
    int64_t *B_i64  = (int64_t*)malloc(sizeof(int64_t)*(size_t)N_ROWS*M_COLS);
    int64_t *TB_i64 = (int64_t*)malloc(sizeof(int64_t)*(size_t)M_COLS*M_COLS);
    if(!A_i64||!U_i64||!T_i64||!B_i64||!TB_i64){
        free(A_i64); free(U_i64); free(T_i64); free(B_i64); free(TB_i64);
        return -100;
    }

    copy_u32_to_i64_mat(&A_in[0][0], N_ROWS, M_COLS, A_i64);
    copy_u32_to_i64_mat(&U[0][0],    M_COLS, M_COLS, U_i64);
    
    for (int r=0;r<M_COLS;++r)
        for (int c=0;c<M_COLS;++c)
            T_i64[(size_t)r*M_COLS + c] = (int64_t)T_in[r][c];

    int rc = NewBasisDel_C(A_i64, U_i64, T_i64,
                           N_ROWS, M_COLS, (int64_t)Q_MODULUS,
                           0, B_i64, TB_i64);
    if (rc != 0){
        fprintf(stderr, "[renew] NewBasisDel_C failed: rc=%d\n", rc);
        free(A_i64); free(U_i64); free(T_i64); free(B_i64); free(TB_i64);
        return rc;
    }

    /* (3) 转回并落盘：At*/ /*/T_At*/ /* 与 U*_inv_cache<t>（即 Uinv_chain） */
    uint32_t A_out[N_ROWS][M_COLS];
    int32_t  T_out[M_COLS][M_COLS];
    copy_i64_to_u32_mat_modq(B_i64,  N_ROWS, M_COLS, &A_out[0][0], (int64_t)Q_MODULUS);
    for (int r=0;r<M_COLS;++r)
        for (int c=0;c<M_COLS;++c)
            T_out[r][c] = (int32_t)TB_i64[(size_t)r*M_COLS + c];

    dump_result_to_temp_txt(t, tag_for_dump, A_out, T_out, Uinv_chain);

    free(A_i64); free(U_i64); free(T_i64); free(B_i64); free(TB_i64);
    return 0;
}

/* ===== 对外接口：仅计算并把结果写入 temp.txt；不写任何内存缓存 ===== */
int renew_update_A_TA_for_t(int t){
    return do_update_and_dump(A,  T_A,  t, "A/T_A");     
}
int renew_update_Ao_TAo_for_t(int t){
    return do_update_and_dump(Ao, T_Ao, t, "Ao/T_Ao");   
}
