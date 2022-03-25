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
 * @param bw         bandwidth
 * @param step       step size; only when n_cigar is not NULL
 * @param score      (out) edit distance
 * @param t_endl     (out) length of target in the alignment
 * @param q_endl     (out) length of query in the alignment
 * @param n_cigar    (in/out) number of cigar operations; NULL if don't need CIGAR
 *
 * @return CIGAR in the htslib packing
 */
uint32_t *lv_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t bw, int32_t step, int32_t *score, int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar);

uint32_t *lv_ed_basic(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t bw, int32_t *score, int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar);
uint64_t *lv_ed_segment(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t bw, int32_t step, int32_t *score, int32_t *n_slice);

#ifdef __cplusplus
}
#endif

#endif
