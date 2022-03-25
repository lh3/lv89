#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "lv89.h"

typedef struct {
	int32_t d;
	int32_t k, p;
} wf_diag_t;

typedef struct {
	int32_t n;
	wf_diag_t *a;
} wf_tb1_t;

typedef struct {
	int32_t m, n;
	wf_tb1_t *a;
} wf_tb_t;

static void wf_tb_add(wf_tb_t *tb, int32_t n, const wf_diag_t *a)
{
	wf_tb1_t *p;
	if (tb->n == tb->m) {
		tb->m += (tb->m>>1) + 4;
		tb->a = (wf_tb1_t*)realloc(tb->a, tb->m * sizeof(*tb->a));
	}
	p = &tb->a[tb->n++];
	p->n = n;
	p->a = (wf_diag_t*)malloc(n * sizeof(wf_diag_t));
	memcpy(p->a, a, n * sizeof(wf_diag_t));
}

static int32_t wf_step(int32_t is_ext, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag_t *a, int32_t *t_end, int32_t *q_end, int32_t *p_end)
{
	int32_t j, st = 0, en = n + 2;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
	*t_end = *q_end = -1;
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, max_k;
		const char *ts_, *qs_;
		uint64_t cmp = 0;
		if (k >= tl || k + p->d >= ql) continue;
		max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		ts_ = ts + 1;
		qs_ = qs + p->d + 1;
		while (k + 7 < max_k) {
			uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
			uint64_t y = *(uint64_t*)(qs_ + k);
			cmp = x ^ y;
			if (cmp == 0) k += 8;
			else break;
		}
		if (cmp)
			k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
		else if (k + 7 >= max_k)
			while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
				++k;
		if (k + p->d == ql - 1 || k == tl - 1) {
			if (is_ext || (k + p->d == ql - 1 && k == tl - 1)) {
				*t_end = k, *q_end = k + p->d, *p_end = p->p;
				return -1;
			}
		}
		p->k = k;
	}

	// wfa_next
	b[0].d = a[0].d - 1;
	b[0].p = a[0].p;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].p =  n == 1 || a[0].k > a[1].k? a[0].p : a[1].p;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k, p = a[j-1].p;
		p = k > a[j+1].k + 1? p : a[j+1].p;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		p = k > a[j].k + 1? p : a[j].p;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k, b[j+1].p = p;
	}
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].p = a[n-2].k > a[n-1].k + 1? a[n-2].p : a[n-1].p;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].p = a[n-1].p;
	b[n+1].k = a[n-1].k;

	if (b[0].d <= -tl) ++st;
	if (b[n+1].d >= ql) --en;
	memcpy(a, &b[st], (en - st) * sizeof(*a));
	return en - st;
}

static uint64_t *wf_traceback(wf_tb_t *tb, int32_t t_end, int32_t q_end, int32_t p_end, int32_t *n_slice)
{
	uint64_t *slice;
	int32_t i, k = 0, p = p_end;
	slice = (uint64_t*)calloc(tb->n + 2, sizeof(uint64_t));
	slice[k++] = (uint64_t)(t_end+1) << 32 | (q_end+1);
	i = tb->n - 1;
	while (p >= 0) {
		wf_diag_t *t = &tb->a[i].a[p];
		slice[k++] = (uint64_t)t->k<<32 | (t->d + t->k);
		p = t->p, --i;
	}
	assert(i == -1);
	slice[k++] = 0;
	assert(k == tb->n + 2);
	for (i = 0; i < k>>1; ++i) { // reverse
		uint64_t tmp = slice[i];
		slice[i] = slice[k-i-1], slice[k-i-1] = tmp;
	}
	*n_slice = k;
	return slice;
}

uint64_t *lv_ed_segment(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t step, int32_t *score, int32_t *n_slice)
{
	int32_t s = 0, n = 1, t_end = -1, q_end = -1, p_end = -1, i;
	wf_diag_t *a;
	uint64_t *slice = 0;
	wf_tb_t tb = {0,0,0};
	assert(tl > 0 && ql > 0);
	a = (wf_diag_t*)malloc(2 * (tl + ql + 2) * sizeof(*a)); // without CIGAR, this would be all the memory needed
	a[0].d = 0, a[0].k = -1, a[0].p = -1;
	while (1) {
		n = wf_step(is_ext, tl, ts, ql, qs, n, a, &t_end, &q_end, &p_end);
		if (n < 0) break;
		++s;
		if (s % step == 0) {
			wf_tb_add(&tb, n, a);
			for (i = 0; i < n; ++i) a[i].p = i;
		}
	}
	free(a);
	slice = wf_traceback(&tb, t_end, q_end, p_end, n_slice);
	for (i = 0; i < tb.n; ++i) free(tb.a[i].a);
	free(tb.a);
	*score = s;
	return slice;
}
