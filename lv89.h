#ifndef LV89_ALN_H
#define LV89_ALN_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Find the edit distance between two sequences
 *
 * @param tl         target sequence length
 * @param ts         target sequence
 * @param ql         query sequence length
 * @param qs         query sequence
 * @param is_ext     extension alignment if true (stop when reaching the end of either query or target)
 * @param score      (out) edit distance
 * @param t_endl     (out) length of target in the alignment
 * @param q_endl     (out) length of query in the alignment
 * @param n_cigar    (out) number of cigar operations
 * @param mem        temporary memory of lv_ed_bufsize(tl, ql) bytes
 *
 * @return CIGAR in the htslib packing
 */
uint32_t *lv_ed_unified(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t *score, int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar);

uint64_t *lv_ed_segment(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t step, int32_t *score, int32_t *n_slice);
int32_t lv_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_global, uint8_t *mem);
int32_t lv_ed_semi(int32_t tl, const char *ts, int32_t ql, const char *qs, uint8_t *mem, int32_t *t_end, int32_t *q_end);
int32_t lv_ed_bufsize(int32_t tl, int32_t ql);

#ifdef __cplusplus
}
#endif

#endif
