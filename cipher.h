#ifndef CIPHER_H
#define CIPHER_H

#include <stdbool.h>
#include "generateKW.h" 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 统计所有 cipher_t.txt 文件中的密文总数。
 *
 * 该函数会遍历 CIPHER_DIR 目录，读取所有 'cipher_*.txt' 文件，
 * 并统计记录分隔符("----")的数量。结果会保存在一个内部静态变量中，
 * 用于决定下一个 record_id。
 * 该函数会在加密时自动调用（如果需要），也可以手动调用以初始化或刷新计数。
 *
 * @return 查找到的密文总数。
 */
int Cipher_count_record(void);

/**
 * @brief 获取当前已加密的密文数量。
 * 
 * 如果计数器尚未初始化，此函数会先调用 Cipher_count_record()。
 *
 * @return 当前内存中的密文数量。
 */
int Cipher_get_record_count(void);

/**
 * @brief 清空 cipher 目录下所有密文文件的内容。
 *
 * 此函数同时会将内部的密文计数器重置为 0。
 */
void Cipher_clear_all_files(void);

/**
 * @brief 对一组关键字字段进行加密。
 *
 * record_id 现在由内部自动管理，其值会自动设为 (当前密文总数 + 1)。
 * 当加密成功并写入文件后，内部密文计数器会自动加一。
 *
 * @param t 时间参数 t。
 * @param fields 指向 kw_field_t 类型数组的指针，表示关键字字段。
 * @param values 指向字符串数组的指针，表示对应的关键字值。
 * @param kprime 本次加密的关键字字段数量。
 * @return 成功返回 true，失败返回 false。
 */
bool Cipher_encrypt_with_fields(int t,
                                const kw_field_t* fields,
                                const char* const* values,
                                int kprime);

#ifdef __cplusplus
}
#endif

#endif