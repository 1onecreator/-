#ifndef GENERATE_KW_H
#define GENERATE_KW_H

#include <stddef.h>
#include <stdint.h>


#include "config_params.h"


enum { KWF_COUNT = KEYWORD_FIELD_COUNT };
typedef int kw_field_t;  


typedef struct {
    int l_bits;
    int k_hashes;
} kw_params_t;


typedef struct {
    const char* key;    
    const char* title; 
} kw_field_meta_t;


const kw_field_meta_t* kw_get_fields(void);


void kw_encode_field(const kw_params_t* P,
                     kw_field_t field_idx,        
                     const char* value,           
                     uint8_t* out_bits);          


size_t kw_bits_to_hex(const uint8_t* bits, size_t l_bits, char* hex_out);


void kw_make_beta_ones(uint8_t* beta_bits, size_t l_bits);

#endif 
