#include "generateKW.h"
#include <string.h>
#include <ctype.h>

static const kw_field_meta_t g_fields[KWF_COUNT] = {
    { "doctor_name",    "Doctor name"    },
    { "patient_age",    "Patient age"    },
    { "patient_gender", "Patient gender" },
    { "department",     "Department"     },
    { "diagnosis",      "Diagnosis"      },
    { "treatment",      "Treatment"      }
};

const kw_field_meta_t* kw_get_fields(void) { return g_fields; }

static inline uint64_t fnv1a64(const void* data, size_t len, uint64_t seed)
{
    const uint8_t* p = (const uint8_t*)data;
    uint64_t h = 1469598103934665603ULL ^ (seed * 0x9e3779b97f4a7c15ULL);
    for (size_t i = 0; i < len; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static inline void set_bit(uint8_t* bits, size_t pos)
{
    bits[pos >> 3] |= (uint8_t)(1u << (pos & 7));
}

static void hash_token_into_bits(const kw_params_t* P,
                                 const char* token, size_t tlen,
                                 uint8_t* out_bits)
{
    if (tlen == 0) return;

    static const uint64_t seeds[5] = {
        0xA5A5A5A55AA55AA5ULL, 0x9E3779B97F4A7C15ULL,
        0x94D049BB133111EBULL, 0xBF58476D1CE4E5B9ULL,
        0x2545F4914F6CDD1DULL
    };
    const int kh = (P && P->k_hashes > 0 && P->k_hashes <= 5) ? P->k_hashes : 3;
    const size_t L = (P && P->l_bits > 0) ? (size_t)P->l_bits : 256u;

    for (int i = 0; i < kh; ++i) {
        uint64_t h = fnv1a64(token, tlen, seeds[i]);
        size_t    p = (size_t)(h % L);
        set_bit(out_bits, p);
    }
}

/* 规范化并分词 */
static void tokenize_and_hash(const kw_params_t* P,
                              const char* s,
                              uint8_t* out_bits)
{
    char buf[1024];
    size_t n = strlen(s);
    if (n >= sizeof(buf)) n = sizeof(buf) - 1;

    size_t i = 0;
    while (i < n) {
        while (i < n && !isalnum((unsigned char)s[i])) ++i;
        size_t st = i;
        while (i < n && isalnum((unsigned char)s[i])) ++i;
        size_t ed = i;
        if (ed > st) {
            size_t tlen = ed - st;
            if (tlen > sizeof(buf)-1) tlen = sizeof(buf)-1;
            for (size_t k = 0; k < tlen; ++k)
                buf[k] = (char)tolower((unsigned char)s[st + k]);
            buf[tlen] = '\0';
            hash_token_into_bits(P, buf, tlen, out_bits);
        }
    }
}

void kw_encode_field(const kw_params_t* P,
                     kw_field_t field_idx,          
                     const char* value,
                     uint8_t* out_bits)
{
    const size_t L   = (P && P->l_bits > 0) ? (size_t)P->l_bits : 256u;
    const char*  val = (value && value[0]) ? value : "None";
    memset(out_bits, 0, (L + 7) / 8);

    if (field_idx < 0 || field_idx >= KWF_COUNT) return;

    
    char mixed[1200];
    const char* key = g_fields[field_idx].key;
    size_t klen = strlen(key);
    size_t vlen = strlen(val);
    if (klen + 1 + vlen >= sizeof(mixed)) vlen = sizeof(mixed) - klen - 2;

    memcpy(mixed, key, klen);
    mixed[klen] = ':';
    memcpy(mixed + klen + 1, val, vlen);
    mixed[klen + 1 + vlen] = '\0';

    tokenize_and_hash(P, mixed, out_bits);
}

size_t kw_bits_to_hex(const uint8_t* bits, size_t l_bits, char* hex_out)
{
    static const char* HEX = "0123456789abcdef";
    const size_t nB = (l_bits + 7u) / 8u;
    size_t w = 0;
    for (size_t i = 0; i < nB; ++i) {
        uint8_t b = bits[i];
        hex_out[w++] = HEX[(b >> 4) & 0xF];
        hex_out[w++] = HEX[b & 0xF];
    }
    hex_out[w] = '\0';
    return w;
}

void kw_make_beta_ones(uint8_t* beta_bits, size_t l_bits)
{
    const size_t nB = (l_bits + 7u) / 8u;
    memset(beta_bits, 0xFF, nB);
    size_t extra = nB * 8u - l_bits;
    if (extra) {
        uint8_t mask = (uint8_t)(0xFFu >> extra);
        beta_bits[nB - 1] = mask;
    }
}
