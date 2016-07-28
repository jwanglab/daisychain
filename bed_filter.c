#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "klib/kvec.h"
#include <htslib/sam.h>
#include <htslib/khash.h>

/*
 * bedfilter.c
 * 
 * Jeremy Wang
 * 20160602
 *
 * Filter locations present in an input BED file
 * from an input BAM file
 * output a new, filtered BAM file
 *
 * Both BAM and BED file must be sorted
 * within each chromosome, but chromosomes
 * may be in any order
*/

// creates string:int hash
KHASH_MAP_INIT_STR(tidMap, int);

typedef struct {
  int start;
  int end;
} range;

typedef kvec_t(range) rangeVec;

// int:vec hash
KHASH_MAP_INIT_INT(rangeMap, rangeVec);

int main(int argc, char *argv[]) {

  if(argc < 4) {
    printf("Usage: bedfilter <BAM> <excluded regions BED> <output BAM>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *exclude_bed = argv[2];
  char *output_bam = argv[3];


  /*
   * Open BAM and read header
   * to create chromosome name -> tid map
   */

  khash_t(tidMap) *tidmap = kh_init(tidMap);
  khash_t(rangeMap) *rangemap = kh_init(rangeMap);

  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int i, absent, ret_val;
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

  int tid;
  for (tid = 0; tid < header->n_targets; tid++) {
    // add chrom: tid to tidmap (a tidMap)
    bin = kh_put(tidMap, tidmap, strdup(header->target_name[tid]), &absent);
    kh_val(tidmap, bin) = tid;

    // add tid: vec(range) to rangemap (a rangeMap)
    bin = kh_put(rangeMap, rangemap, tid, &absent);
    kv_init(kh_val(rangemap, bin));
  }
  // array of integers for the indices into the corresponding tid range array
  int *tid_range_index = calloc(header->n_targets, sizeof(int));


  /*
   * Load and hash regions in exclusion BED file
   */

  FILE* fp;
  fp = fopen(exclude_bed, "r");
  char* line = NULL;
  size_t len = 0;
  ssize_t chr_read;

  // all of these must be allocated or sscanf will segfault
  char chrom[100], remainder[1000];
  int start = 0, end = 0;

  // read BED file line by line
  while ((chr_read = getline(&line, &len, fp)) != -1) {

    // chrom  start  end  remainder...
    // chr1	29232	29246	...
    ret_val = sscanf(line, "%s\t%d\t%d\t%s", chrom, &start, &end, remainder);
    if (ret_val == -1) {
      fprintf(stderr, "Error scanning line '%s'\n", line);
      return 1;
    }

    // look up tid by chromosome
    bin = kh_get(tidMap, tidmap, chrom);
    if(bin == kh_end(tidmap)) { // key is missing
      // this happens with chrY sometimes
      // means there are no reads aligned to it anyway, so skip it
      continue;
    }
    tid = kh_value(tidmap, bin);

    // add range to the correct tid vector
    bin = kh_get(rangeMap, rangemap, tid);
    range r;
    r.start = start;
    r.end = end;
    kv_push(range, kh_value(rangemap, bin), r); // append
  }


  /*
   * Loop through BAM file, write only reads which
   * are aligned AND
   * do not land in exclusion file
   * - exclusion ranges are strand-agnostic
   */

  aln = bam_init1();

  samFile *out_bam = sam_open(output_bam, "wb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\" for writing\n", output_bam);
    return -1;
  }
  ret_val = sam_hdr_write(out_bam, header); // copy header

  uint32_t kept = 0;
  uint32_t removed = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;

    bin = kh_get(rangeMap, rangemap, tid);
    rangeVec ranges = kh_value(rangemap, bin);

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    while(pos > kv_A(ranges, tid_range_index[tid]).end) {
      tid_range_index[tid]++;
    }

    // overlaps
    if(endpos >= kv_A(ranges, tid_range_index[tid]).start) {
      removed++;
      continue;
    }

    ret_val = sam_write1(out_bam, header, aln);

    kept++;
  }
  printf("%i reads kept.\n", kept);
  printf("%i reads removed.\n", removed);

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(out_bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing output BAM.\n");
    return -1;
  }

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  return 0;
}

