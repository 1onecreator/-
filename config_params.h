#ifndef CONFIG_PARAMS_H
#define CONFIG_PARAMS_H

#include <stdint.h>
#include <stddef.h>

/* ========= 基本参数 ========= */
enum {
    L_BITS      = 128,     // l
    Q_MODULUS   = 65537,   // q
    B_BITS      = 17,      // b
    K_LEVEL     = 6,       // k
    M_COLS      = 512,     // m
    N_ROWS      = 256,     // n
    T_PERIOD    = 7        // t
};

/* beta 向量（长度 4） */
#define BETA_LEN 4
static const uint8_t BETA[BETA_LEN] = {1, 0, 1, 0};

/* ========= 高斯 / 其他实数参数 ========= */
#define OMEGA_M     9.42               // omigam
#define DELTA_R     354.9187213444001  // delta_R = sqrt(n * ln(q)) * OMEGA_M
#define DELTA_BIG   2.12e15            // delta
#define SIGMA       3.0               // sigma
#define SIGMA_S     3.0                // sigma_s

/* ========= 医疗关键字字段 ========= */
typedef enum {
    FIELD_DOCTOR_NAME   = 1,
    FIELD_PATIENT_AGE   = 2,
    FIELD_PATIENT_GENDER= 3,
    FIELD_DEPARTMENT    = 4,
    FIELD_DIAGNOSIS     = 5,
    FIELD_TREATMENT     = 6
} KeywordField;

#define KEYWORD_FIELD_COUNT 6
#ifndef K_FIELDS
#define K_FIELDS KEYWORD_FIELD_COUNT
#endif

typedef struct {
    KeywordField id;
    const char*  name;
} FieldEntry;

static const FieldEntry MEDICAL_KEYWORD_FIELDS[KEYWORD_FIELD_COUNT] = {
    { FIELD_DOCTOR_NAME,    "doctor_name"    },
    { FIELD_PATIENT_AGE,    "patient_age"    },
    { FIELD_PATIENT_GENDER, "patient_gender" },
    { FIELD_DEPARTMENT,     "department"     },
    { FIELD_DIAGNOSIS,      "diagnosis"      },
    { FIELD_TREATMENT,      "treatment"      }
};

/* ========= 额外常量 ========= */
#define L_COLS 256   

#endif 
