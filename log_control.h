#ifndef LOG_CONTROL_H
#define LOG_CONTROL_H

#include <stdio.h>

/**
 * @brief 全局日志详细模式开关
 *
 * 设置为 1: 输出详细的调试和过程信息。
 * 设置为 0: 只输出关键的状态信息、最终结果和错误。
 */
#define GLOBAL_VERBOSE_MODE 0

/*
 * @brief 详细日志宏
 *
 * 当 GLOBAL_VERBOSE_MODE 为 1 时，此宏会打印消息到 stderr。
 * 当 GLOBAL_VERBOSE_MODE 为 0 时，此宏不执行任何操作。
 */
#if GLOBAL_VERBOSE_MODE
    #define LOG_DETAIL(...) fprintf(stderr, __VA_ARGS__)
#else
    #define LOG_DETAIL(...) do {} while (0)
#endif

#endif 
