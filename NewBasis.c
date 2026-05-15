
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "NewBasis.h"

/* ---------- 基础：安全取模 / 扩展欧几里得 / 逆元 ---------- */
static inline int64_t floormod_i64(int64_t x, int64_t q){
    int64_t r = x % q; return (r < 0) ? (r + q) : r;
}
static int64_t egcd_i64(int64_t a, int64_t b, int64_t* x, int64_t* y){
    if (b == 0){ *x = 1; *y = 0; return a; }
    int64_t x1,y1; int64_t g = egcd_i64(b, a % b, &x1, &y1);
    *x = y1; *y = x1 - (a / b) * y1; return g;
}
static int64_t modinv_i64(int64_t a, int64_t q){
    a = floormod_i64(a, q);
    int64_t x,y; int64_t g = egcd_i64(a, q, &x, &y);
    if (g != 1 && g != -1) return 0;
    int64_t inv = x % q; if (inv < 0) inv += q; return inv;
}
static inline int64_t center_coeff(int64_t v, int64_t q){
    v = floormod_i64(v, q);
    int64_t half = q/2;
    if (v > half) v -= q;
    return v;
}

/* ---------- 线代：矩阵乘法（mod q）与求逆（Z_q） ---------- */
static void mat_mul_mod_q(const int64_t* A, const int64_t* B,
                          int r, int k, int c, int64_t q,
                          int64_t* Out){
    for (int i=0;i<r;++i){
        for (int j=0;j<c;++j){
            __int128 acc = 0;
            for (int t=0;t<k;++t){
                int64_t a = floormod_i64(A[(size_t)i*k + t], q);
                int64_t b = floormod_i64(B[(size_t)t*c + j], q);
                acc += (__int128)a * b;
                
                if ((t & 15)==15){ acc %= q; }
            }
            int64_t z = (int64_t)(acc % q);
            if (z < 0) z += q;
            Out[(size_t)i*c + j] = z;
        }
    }
}

static int mat_inv_mod_q(const int64_t* R, int m, int64_t q, int64_t* Rinv_out){
    int64_t* A = (int64_t*)malloc(sizeof(int64_t)*(size_t)m*m);
    int64_t* I = (int64_t*)malloc(sizeof(int64_t)*(size_t)m*m);
    if (!A || !I) { free(A); free(I); return -2; }
    for (int i=0;i<m*m;++i) A[i] = floormod_i64(R[i], q);
    for (int i=0;i<m*m;++i) I[i] = 0;
    for (int i=0;i<m;  ++i) I[(size_t)i*m + i] = 1;

    int row = 0;
    for (int col=0; col<m && row<m; ++col){
        int pivot = -1;
        for (int r=row; r<m; ++r){
            if (A[(size_t)r*m + col] % q != 0){ pivot = r; break; }
        }
        if (pivot < 0){
            
            free(A); free(I); return -1;
        }
        if (pivot != row){
            for (int c=0;c<m;++c){
                int64_t ta = A[(size_t)row*m + c]; A[(size_t)row*m + c] = A[(size_t)pivot*m + c]; A[(size_t)pivot*m + c] = ta;
                int64_t ti = I[(size_t)row*m + c]; I[(size_t)row*m + c] = I[(size_t)pivot*m + c]; I[(size_t)pivot*m + c] = ti;
            }
        }
        int64_t pv  = floormod_i64(A[(size_t)row*m + col], q);
        int64_t inv = modinv_i64(pv, q); if (!inv){ free(A); free(I); return -1; }

        for (int c=0;c<m;++c){
            A[(size_t)row*m + c] = floormod_i64(A[(size_t)row*m + c] * inv, q);
            I[(size_t)row*m + c] = floormod_i64(I[(size_t)row*m + c] * inv, q);
        }
        for (int r = 0; r < m; ++r) if (r != row) {
    int64_t f = floormod_i64(A[(size_t)r*m + col], q);
    if (!f) continue;
    for (int c = 0; c < m; ++c) {
        A[(size_t)r*m + c] =
            floormod_i64(A[(size_t)r*m + c] - (__int128)f * A[(size_t)row*m + c], q);
        I[(size_t)r*m + c] =
            floormod_i64(I[(size_t)r*m + c] - (__int128)f * I[(size_t)row*m + c], q);
    }
}
        ++row;
    }
    
    for (int i=0;i<m*m;++i) Rinv_out[i] = I[i];
    free(A); free(I);
    return 0;
}


#ifndef NEWBASIS_ENABLE_LLL
#define NEWBASIS_ENABLE_LLL 1   
#endif

#if NEWBASIS_ENABLE_LLL
static void lll_reduction_int(const int64_t* Tin, int m, double delta, int64_t* Tout){
    
    int64_t* B = (int64_t*)malloc(sizeof(int64_t)*(size_t)m*m);
    memcpy(B, Tin, sizeof(int64_t)*(size_t)m*m);
    double* bcol =(double*)malloc(sizeof(double)*(size_t)m*m);
    double* bstar=(double*)malloc(sizeof(double)*(size_t)m*m);
    double* mu   =(double*)calloc((size_t)m*m, sizeof(double));
    double* n2   =(double*)malloc(sizeof(double)*(size_t)m);

    int k=1, loops=0, maxloops=4000;   
    while (k<m && loops++<maxloops){
        
        for(int j=0;j<m;++j) for(int r=0;r<m;++r) bcol[(size_t)j*m+r]=(double)B[(size_t)r*m+j];
        
        for(int i=0;i<m;++i){
            for(int r=0;r<m;++r) bstar[(size_t)i*m+r]=bcol[(size_t)i*m+r];
            for(int j=0;j<i;++j){
                double denom = n2[j] > 1e-18 ? n2[j] : 1e-18;
                double dot=0.0; for(int r=0;r<m;++r) dot += bcol[(size_t)i*m+r]*bstar[(size_t)j*m+r];
                double muij = dot/denom; mu[(size_t)i*m+j]=muij;
                for(int r=0;r<m;++r) bstar[(size_t)i*m+r]-=muij*bstar[(size_t)j*m+r];
            }
            double s=0.0; for(int r=0;r<m;++r){ double v=bstar[(size_t)i*m+r]; s+=v*v; }
            n2[i]=(s>1e-18)?s:1e-18;
        }
        
        for(int j=k-1;j>=0;--j){
            long qround=(long)llround(mu[(size_t)k*m+j]);
            if(qround) for(int r=0;r<m;++r) B[(size_t)r*m+k]-=qround*B[(size_t)r*m+j];
        }
        
        double lhs=n2[k];
        double mu_k_k1=mu[(size_t)k*m+(k-1)];
        double rhs=(delta - mu_k_k1*mu_k_k1)*n2[k-1];
        if(lhs >= rhs - 1e-12) ++k;
        else{
            for(int r=0;r<m;++r){
                int64_t tmp=B[(size_t)r*m+k];
                B[(size_t)r*m+k]=B[(size_t)r*m+(k-1)];
                B[(size_t)r*m+(k-1)]=tmp;
            }
            if(k>1) --k;
        }
    }
    memcpy(Tout,B,sizeof(int64_t)*(size_t)m*m);
    free(B); free(bcol); free(bstar); free(mu); free(n2);
}
#endif

/* ---------- 公开接口 ---------- */
int NewBasisDel_C(const int64_t* A, const int64_t* R, const int64_t* T_A,
                  int n, int m, int64_t q, int do_lll,
                  int64_t* B_out, int64_t* TB_out)
{
    if(!A || !R || !T_A || !B_out || !TB_out || n<=0 || m<=0 || q<=1) return -1000;

    /* (0) 验证 A*T_A ≡ 0 (mod q) */
    int64_t* AT = (int64_t*)malloc(sizeof(int64_t)*(size_t)n*m);
    if(!AT) return -1001;
    mat_mul_mod_q(A, T_A, n, m, m, q, AT);
    for(int i=0;i<n*m;++i){ if(AT[i] != 0){ free(AT); return -1; } }
    free(AT);

    /* (1) 计算 R^{-1} (mod q) */
    int64_t* Rinv = (int64_t*)malloc(sizeof(int64_t)*(size_t)m*m);
    if(!Rinv) return -1002;
    if(mat_inv_mod_q(R, m, q, Rinv) != 0){ free(Rinv); return -4; }

    /* (2) B = A * R^{-1} (mod q) */
    mat_mul_mod_q(A, Rinv, n, m, m, q, B_out);

    /* (3) TB = R * T_A */
    mat_mul_mod_q(R, T_A, m, m, m, q, TB_out);
    for(int i=0;i<m*m;++i) TB_out[i] = center_coeff(TB_out[i], q);

#if NEWBASIS_ENABLE_LLL
    if(do_lll){
        int64_t* Ttmp = (int64_t*)malloc(sizeof(int64_t)*(size_t)m*m);
        if(!Ttmp){ free(Rinv); return -1003; }
        lll_reduction_int(TB_out, m, 0.99, Ttmp);
        memcpy(TB_out, Ttmp, sizeof(int64_t)*(size_t)m*m);
        free(Ttmp);
    }
#else
    (void)do_lll; 
#endif

    /* (4) 一致性检查：B*TB ≡ 0 (mod q) */
    int64_t* BT = (int64_t*)malloc(sizeof(int64_t)*(size_t)n*m);
    if(!BT){ free(Rinv); return -1004; }
    mat_mul_mod_q(B_out, TB_out, n, m, m, q, BT);
    for(int i=0;i<n*m;++i){ if(BT[i] != 0){ free(BT); free(Rinv); return -5; } }
    free(BT);
    free(Rinv);
    return 0;
}


void NewBasis_SampleUpperR(int m, int64_t q, uint64_t seed, int64_t* R_out){
    (void)seed; 
    for(int i=0;i<m*m;++i) R_out[i]=0;
    for(int i=0;i<m;++i) R_out[(size_t)i*m+i]=1;
    
    for(int i=0;i<m;++i) for(int j=i+1;j<m;++j)
        R_out[(size_t)i*m + j] = (int64_t)((i + 3*j + 7) % 5) % q;
}
