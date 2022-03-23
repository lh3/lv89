#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "lv89.h"
#include "edlib.h"
#ifdef _USE_WFA2
#include "wavefront/wavefront_align.h"
#endif
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, s, is_global = 0, use_edlib = 0, use_wfa = 0, report_cigar = 0, is_semi = 0;
	int32_t n_cigar;
	uint32_t *cigar = 0;

	while ((c = ketopt(&o, argc, argv, 1, "gewcsw", 0)) >= 0) {
		if (c == 'g') is_global = 1;
		else if (c == 's') is_semi = 1;
		else if (c == 'e') use_edlib = 1;
		else if (c == 'w') use_wfa = 1;
		else if (c == 'c') report_cigar = 1;
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: ed-test [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -g    count gaps at the end of the target\n");
		fprintf(stderr, "  -s    stop when reaching the end of either query or target\n");
		fprintf(stderr, "  -c    report CIGAR (only working with -s; not supporting -e or -w)\n");
		fprintf(stderr, "  -e    use edlib (incompatible with -s or -w)\n");
#ifdef _USE_WFA2
		fprintf(stderr, "  -w    use WFA2 (implying -g; incompatible with -s or -e)\n");
#endif
		return 1;
	}
	assert(!use_edlib || !use_wfa);

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

	if (use_edlib) {
		EdlibAlignResult rst;
		fprintf(stderr, "Using edlib...\n");
		rst = edlibAlign(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l,
				edlibNewAlignConfig(-1, is_global? EDLIB_MODE_NW : EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));
		s = rst.editDistance;
#ifdef _USE_WFA2
	} else if (use_wfa) {
		fprintf(stderr, "Using WFA2-lib...\n");
		wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
		attributes.distance_metric = edit;
		attributes.alignment_scope = compute_score;
		attributes.heuristic.strategy = wf_heuristic_none;
		wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
		wavefront_align(wf_aligner, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		s = wf_aligner->align_status.score;
		wavefront_aligner_delete(wf_aligner);
#endif
	} else {
		fprintf(stderr, "Using lv89...\n");
		if (report_cigar) {
			if (is_semi)
				cigar = lv_ed_semi_cigar(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, &s, &n_cigar);
			else {
				fprintf(stderr, "ERROR: not implemented\n");
				abort();
			}
		} else {
			uint8_t *mem = (uint8_t*)malloc(lv_ed_bufsize(ks1->seq.l, ks2->seq.l));
			if (is_semi)
				s = lv_ed_semi(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, mem);
			else
				s = lv_ed(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, is_global, mem);
			free(mem);
		}
	}

	if (!report_cigar) {
		printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, s);
	} else {
		int32_t i;
		printf("%s\t%s\t%d\t", ks1->name.s, ks2->name.s, s);
		for (i = 0; i < n_cigar; ++i)
			printf("%d%c", cigar[i]>>4, "MID"[cigar[i]&0xf]);
		putchar('\n');
	}
	free(cigar);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
