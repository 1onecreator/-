
#include "TrapGen.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>



/* ---------------- 小工具：模运算 ---------------- */
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
static inline uint32_t rng_u32_modq(void){
    return (uint32_t)(xorshift64star() % (uint64_t)Q_MODULUS);
}
static inline int rng_small_trits(void){
    
    uint64_t r = xorshift64star();
    int v = (int)(r % 3ULL);
    return (v==0)?-1: (v==1?0:1);
}

/* ---------------- 线性代数：模 q 消元/求逆 ---------------- */
static long long egcd_ll(long long a, long long b, long long* x, long long* y){
    if (b == 0){ *x = 1; *y = 0; return a; }
    long long x1,y1; long long g = egcd_ll(b, a%b, &x1, &y1);
    *x = y1; *y = x1 - (a/b)*y1; return g;
}
static uint32_t inv_mod_q(uint32_t a){
    long long x,y; (void)egcd_ll((long long)a, (long long)Q_MODULUS, &x, &y);
    long long r = x % (long long)Q_MODULUS; if (r<0) r += Q_MODULUS;
    return (uint32_t)r;
}

/* 高斯–约旦：判断可逆并给出逆（模 q） */
static bool invert_mod_q(const uint32_t A1[N_ROWS][N_ROWS],
                         uint32_t A1_inv[N_ROWS][N_ROWS]){
    uint32_t aug[N_ROWS][2*N_ROWS];
    for(int i=0;i<N_ROWS;++i){
        for(int j=0;j<N_ROWS;++j){ aug[i][j] = A1[i][j] % Q_MODULUS; }
        for(int j=0;j<N_ROWS;++j){ aug[i][N_ROWS+j] = (i==j)?1u:0u; }
    }
    for(int col=0; col<N_ROWS; ++col){
        int piv=-1;
        for(int r=col; r<N_ROWS; ++r){ if(aug[r][col]!=0u){piv=r; break;} }
        if(piv<0) return false;
        if(piv!=col){
            for(int j=0;j<2*N_ROWS;++j){ uint32_t t=aug[col][j]; aug[col][j]=aug[piv][j]; aug[piv][j]=t; }
        }
        uint32_t inv = inv_mod_q(aug[col][col]);
        for(int j=0;j<2*N_ROWS;++j)
            aug[col][j] = (uint32_t)(((__int128)aug[col][j]*inv) % Q_MODULUS);
        for(int r=0;r<N_ROWS;++r) if(r!=col){
            uint32_t f = aug[r][col];
            if(!f) continue;
            for(int j=0;j<2*N_ROWS;++j){
                __int128 val = (__int128)aug[r][j] - (__int128)f*aug[col][j];
                aug[r][j] = mod_q_i128(val);
            }
        }
    }
    for(int i=0;i<N_ROWS;++i)
        for(int j=0;j<N_ROWS;++j)
            A1_inv[i][j] = aug[i][N_ROWS+j];
    return true;
}


bool TrapGen_generate(uint32_t A[N_ROWS][M_COLS],
                      int32_t  T_A[M_COLS][M_COLS]){
    if(!A || !T_A) return false;
    rng_seed(0); 

    /* 1) 随机可逆 A1 */
    uint32_t A1[N_ROWS][N_ROWS];
    uint32_t A1_inv[N_ROWS][N_ROWS];
    int tries=0;
    do{
        for(int i=0;i<N_ROWS;++i)
            for(int j=0;j<N_ROWS;++j)
                A1[i][j] = rng_u32_modq();
        tries++;
        if(tries>1000){ fprintf(stderr,"[TrapGen] fail to find invertible A1.\n"); return false; }
    }while(!invert_mod_q(A1, A1_inv));

    /* 2) 小整数 R (n × (m-n)) */
    const int K = M_COLS - N_ROWS;  
    int32_t R[N_ROWS][ (M_COLS> N_ROWS)?(M_COLS - N_ROWS):1 ];
    for(int i=0;i<N_ROWS;++i)
        for(int j=0;j<K;++j)
            R[i][j] = rng_small_trits(); 

    /* 3) A2 = A1 * R (mod q)  —— 尺寸 n×(m-n) */
    uint32_t A2[N_ROWS][ (M_COLS> N_ROWS)?(M_COLS - N_ROWS):1 ];
    for(int i=0;i<N_ROWS;++i){
        for(int j=0;j<K;++j){
            __int128 acc=0;
            for(int k=0;k<N_ROWS;++k) acc += (__int128)A1[i][k] * (long long)R[k][j];
            A2[i][j] = mod_q_i128(acc);
        }
    }

    /* 4) 拼 A = [A1 | A2] */
    for(int i=0;i<N_ROWS;++i){
        for(int j=0;j<N_ROWS;++j) A[i][j] = A1[i][j];
        for(int j=0;j<K;++j)      A[i][N_ROWS + j] = A2[i][j];
    }

    /* 5) 构造 m×m 的 T_A
     *    前 (m-n) 列是 [-R; I] （短），其余 n 列直接放 q*e（任意位置都行）
     */
    
    for(int r=0;r<M_COLS;++r)
        for(int c=0;c<M_COLS;++c)
            T_A[r][c] = 0;

    
    for(int c=0; c<K; ++c){
        
        for(int r=0; r<N_ROWS; ++r) T_A[r][c] = -R[r][c];
        
        for(int r=0; r<K; ++r)      T_A[N_ROWS + r][c] = (r==c)? 1 : 0;
    }
    
    for(int j=0; j<N_ROWS; ++j){
        int col = K + j;              
        int row = j % M_COLS;         
        T_A[row][col] = (int32_t)Q_MODULUS;
    }

    return true;
}

/* ---------------- TrapGen_verify ---------------- */
bool TrapGen_verify(const uint32_t A[N_ROWS][M_COLS],
                    const int32_t  T[M_COLS][M_COLS]){
    if(!A || !T) return false;
    for(int i=0;i<N_ROWS;++i){
        for(int c=0;c<M_COLS;++c){
            __int128 acc=0;
            for(int k=0;k<M_COLS;++k) acc += (__int128)A[i][k]*(long long)T[k][c];
            if(mod_q_i128(acc) != 0u){
                
                return false;
            }
        }
    }
    return true;
}

/* ---------------- TrapGen_write_public_params ----------------
 * 直接生成一个新的 public_params.h 文件（包含 A 与 T_A 的“定义”）
 * 注意：这是“写死常量到头文件”的方案；若你更倾向 .c/.h 分离，可自行改写为写入 .c 文件。
 */
static void fprint_mat_u32(FILE* fp, const char* name,
                           const uint32_t* data, int rows, int cols){
    fprintf(fp, "const uint32_t %s[%d][%d] = {\n", name, rows, cols);
    for(int i=0;i<rows;++i){
        fputs("  {", fp);
        for(int j=0;j<cols;++j){
            fprintf(fp, "%u", data[i*cols + j]);
            if(j+1<cols) fputc(',', fp);
        }
        fputs("}", fp);
        if(i+1<rows) fputc(',', fp);
        fputc('\n', fp);
    }
    fputs("};\n", fp);
}
static void fprint_mat_i32(FILE* fp, const char* name,
                           const int32_t* data, int rows, int cols){
    fprintf(fp, "const int32_t %s[%d][%d] = {\n", name, rows, cols);
    for(int i=0;i<rows;++i){
        fputs("  {", fp);
        for(int j=0;j<cols;++j){
            fprintf(fp, "%d", data[i*cols + j]);
            if(j+1<cols) fputc(',', fp);
        }
        fputs("}", fp);
        if(i+1<rows) fputc(',', fp);
        fputc('\n', fp);
    }
    fputs("};\n", fp);
}

bool TrapGen_write_public_params(const char* path,
                                 const uint32_t A_mat[N_ROWS][M_COLS],
                                 const int32_t  T_mat[M_COLS][M_COLS]){
    if(!path || !A_mat || !T_mat) return false;
    FILE* fp = fopen(path, "w");
    if(!fp) return false;

    fprintf(fp, "// public_params.h — generated by TrapGen\n");
    
    
    fprint_mat_u32(fp, "A",    &A_mat[0][0], N_ROWS, M_COLS);
    fprint_mat_i32(fp, "T_A",  &T_mat[0][0], M_COLS, M_COLS);

    
    

    fclose(fp);
    return true;
}
