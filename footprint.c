#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"
#include "klib/kvec.h"
#include <float.h>
#include <time.h>
#include <math.h>

/*
 * footprint.c
 * 
 * Jeremy Wang
 * 20160803
 *
 * Identify and score putative DNAse-seq footprints
 *
 * 1. Take genome-wide pileup and initial set of putative binding sites
 *    found ex. by tag density, binding site PWM scoring, etc.
 * 2. Compute height kernel (analogous to PWM) across entire footprint window (size given)
 *    at putative binding sites
 * 3. Score kernel (like a PWM) across entire genome
 * 4. Create new set of putative binding sites using some threshold
 *
 * Repeat 2-4 several times (until set converges)
 *
 * Other considerations:
 * - what if there are multiple different kernels?
 *   how do we identify and separate them?
 * - are positions really independent of each other?
*/

typedef struct read {
  uint32_t pos;
  float weight;
} read;

typedef kvec_t(read) readList;

// kernel defines a shaped PWM of a given size
typedef struct kernel {
  uint32_t width;
  double* vals;
} kernel;

double gauss(double x) {
  return pow(M_E, -(x * x)/2) / PI2;
}


int main(int argc, char *argv[]) {

  // seed random
  srand(time(NULL));

  if(argc < 6) {
    printf("Usage: footprint <BAM> <motif_bed> <footprint_width> <assay {d,a,c,f}> <output BED>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *motif_file = argv[2];
  int footprint_width = atoi(argv[3]);
  char assay = argv[4][0];
  char *output_bed = argv[5];

  if(assay != 'd' && assay != 'c' && assay != 'a' && assay != 'f') {
    fprintf(stderr, "Unknown assay type '%c', should be one of 'd' (DNase), 'c' (ChIP), 'f' (FAIRE), or 'a' (ATAC)\n", assay);
    return -1;
  }


  // initialize BAM file
  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int i, absent, ret_val, ct;
  khint_t bin; // hash bin (result of kh_put)

  bam = sam_open(bam_file, "rb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", bam_file);
    return -1;
  }
  header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return -1;
  }

  // load the index, just to get the total # seq without going through the whole thing
  hts_idx_t *idx;
  // searches for bam_file + .bai, etc
  if ((idx = sam_index_load2(bam, bam_file, NULL)) == 0) {
    fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
    return 1;
  }

  double genome_size = 0;
  uint64_t mapped, unmapped, total_mapped = 0;

  printf("%i chromosomes in BAM file.\n", header->n_targets);

  float **coverage = malloc(sizeof(float*) * header->n_targets); // ~12Gb

  for (i = 0; i < header->n_targets; i++) {
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);

    genome_size += header->target_len[i];

    coverage[i] = (float*)calloc(header->target_len[i], sizeof(float));
    printf("chrom %i\n", i);

    hts_idx_get_stat(idx, i, &mapped, &unmapped);
    total_mapped += mapped;
  }

  printf("Genome size: %f\n", genome_size);
  printf("Mapped reads: %llu\n", total_mapped);

  int read_len = -1;

  /*
   * Build simple arrays per chromosome of tag/coverage counts
   */
  aln = bam_init1();
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    int32_t tid = aln->core.tid;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1; // inclusive

    if(read_len == -1) {
      read_len = endpos - pos + 1;
    }

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    uint8_t *tag = bam_aux_get(aln, "XW");
    float weight = 1.0;
    if (tag != 0) { // tag exists
      weight = bam_aux2f(tag);
      if (weight == 0.0) { // some datasets have this tag set before I touch them, most (maybe all) are zero
        weight = 1.0;
      }
    }

    int offset = 0; // offset of cut site from read start
    if(assay == 'd') {
      offset = 0;
    } else if (assay == 'a') {
      offset = 4; // sometimes 5 is used on the fw strand and 4 on rv, we'll use 4 all the time
    }
    // faire and chip both use the entire read

    if(assay == 'd' || assay == 'a') {
      if (bam_is_rev(aln)) {
        coverage[tid][endpos - offset] += weight;
      } else {
        coverage[tid][pos + offset] += weight;
      }
    } else if(assay == 'c' || assay == 'f') {
      for(i = pos; i <= endpos; i++) {
        coverage[tid][i] += weight;
      }
    }
  }

  // clean up BAM file
  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing BAM file.\n");
  }

  FILE *bed_fout;
  // output bed format may be nonstandard: "chrom startPos endPos avgDensityScore"
  bed_fout = fopen(output_bed, "w");
  
  int32_t tid;
  float score;
  int start;
  for(tid = 0; tid < header->n_targets; tid++) {
    r = 0; // start at the first read on this chromosome kv_A(reads[tid], r)
    printf("Chrom %s, %i bp...\n", header->target_name[tid], header->target_len[tid]);

    for(i = 0; i < header->target_len[tid]; i++) {
      if(i % 1000000 == 0) {
        printf("At pos %i\n", i);
      }

      score = density_by_coverage(i, distr, b, w, coverage[tid], header->target_len[tid]);

      fprintf(bed_fout, "%s\t%i\t%i\t%f\n", header->target_name[tid], start, i, score); // end position is EXclusive
    }
    break;
  }

  ret_val = fclose(bed_fout);

  bam_hdr_destroy(header);

  return 0;
}
