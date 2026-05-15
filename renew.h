#ifndef RENEW_H
#define RENEW_H

#ifdef __cplusplus
extern "C" {
#endif

/* 仅更新 A 与 T_A（写回 public_params.h 中的 A 与 T_A） */
int renew_update_A_TA_for_t(int t);

/* 仅更新 Ao 与 T_Ao（写回 public_params.h 中的 Ao 与 T_Ao） */
int renew_update_Ao_TAo_for_t(int t);

#ifdef __cplusplus
}
#endif
#endif 
