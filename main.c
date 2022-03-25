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
	int c, s, step = 0, is_ext = 0, use_edlib = 0, use_wfa = 0, report_cigar = 0, mem_mode = 2, bw = 0;
	int32_t n_cigar, t_endl, q_endl;
	uint32_t *cigar = 0;
	char *cigar_str = 0;

	while ((c = ketopt(&o, argc, argv, 1, "ewcxs:m:b:", 0)) >= 0) {
		if (c == 'x') is_ext = 1;
		else if (c == 's') step = atoi(o.arg);
		else if (c == 'e') use_edlib = 1;
		else if (c == 'w') use_wfa = 1;
		else if (c == 'c') report_cigar = 1;
		else if (c == 'b') bw = atoi(o.arg);
		else if (c == 'm') mem_mode = atoi(o.arg);
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: ed-test [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -x      extension mode\n");
		fprintf(stderr, "  -c      report CIGAR (implying -u; not supporting -e or -w)\n");
		fprintf(stderr, "  -b INT  band width (<=0 to disable) [%d]\n", bw);
		fprintf(stderr, "  -e      use edlib (not supporting -x)\n");
#ifdef _USE_WFA2
		fprintf(stderr, "  -w      use WFA2 (not supporting -x)\n");
		fprintf(stderr, "  -m INT  memory mode: 1=low, 2=med and 3=high (for -w only) [%d]\n", mem_mode);
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

	t_endl = ks1->seq.l, q_endl = ks2->seq.l;

	if (use_edlib) {
		EdlibAlignResult rst;
		fprintf(stderr, "Using edlib...\n");
		rst = edlibAlign(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l,
				edlibNewAlignConfig(-1, is_ext? EDLIB_MODE_SHW : EDLIB_MODE_NW, report_cigar? EDLIB_TASK_PATH : EDLIB_TASK_DISTANCE, NULL, 0));
		s = rst.editDistance;
		if (report_cigar)
			cigar_str = edlibAlignmentToCigar(rst.alignment, rst.alignmentLength, EDLIB_CIGAR_EXTENDED);
		edlibFreeAlignResult(rst);
#ifdef _USE_WFA2
	} else if (use_wfa) {
		fprintf(stderr, "Using WFA2-lib...\n");
		wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
		attributes.distance_metric = edit;
		attributes.alignment_scope = report_cigar? compute_alignment : compute_score;
		attributes.memory_mode = mem_mode <= 1? wavefront_memory_low : mem_mode == 2? wavefront_memory_med : wavefront_memory_high;
		attributes.heuristic.strategy = wf_heuristic_none;
		wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
		wavefront_align(wf_aligner, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		s = wf_aligner->align_status.score;
		wavefront_aligner_delete(wf_aligner);
#endif
	} else {
		fprintf(stderr, "Using lv89...\n");
		cigar = lv_ed(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, is_ext, bw, step, &s, &t_endl, &q_endl, report_cigar? &n_cigar : 0);
	}

	printf("%s\t%ld\t0\t%d\t+\t%s\t%ld\t0\t%d\t%d", ks1->name.s, ks1->seq.l, t_endl, ks2->name.s, ks2->seq.l, q_endl, s);
	if (report_cigar && !use_wfa) {
		int32_t i, ed = 0;
		putchar('\t');
		if (cigar_str) {
			fputs(cigar_str, stdout);
		} else {
			for (i = 0; i < n_cigar; ++i) {
				printf("%d%c", cigar[i]>>4, "MIDNSHP=XB"[cigar[i]&0xf]);
				if ((cigar[i]&0xf) != 7) ed += cigar[i]>>4;
			}
			assert(ed == s);
		}
		free(cigar);
	}
	putchar('\n');

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
