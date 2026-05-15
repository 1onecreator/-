#ifndef NEWBASIS_H
#define NEWBASIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* NewBasisDel_C
 * 输入：
 *   A(n×m), R(m×m), T_A(m×m), 模数 q，开关 do_lll（非 0 则执行 LLL）
 * 输出：
 *   B_out(n×m) = A * R^{-1} (mod q)
 *   TB_out(m×m) = R * T_A （内部按 mod q 做中心化存储，之后再可选 LLL）
 * 返回：0 成功；负数表示错误。
 */
int NewBasisDel_C(const int64_t* A, const int64_t* R, const int64_t* T_A,
                  int n, int m, int64_t q, int do_lll,
                  int64_t* B_out, int64_t* TB_out);

/* 可选工具：采样一个上三角（对角 1）且在 Z_q 上可逆的 R（m×m）
 * 便于外部构造 NewBasisDel 所需的 R。
 */
void NewBasis_SampleUpperR(int m, int64_t q, uint64_t seed, int64_t* R_out);

#ifdef __cplusplus
}
#endif

#endif /* NEWBASIS_H */
