#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"
#include "htslib/htslib/ksort.h"
#include "htslib/htslib/kseq.h"

/*
 * expander.c
 * 
 * Take a BAM file with XW weight tags and write a BAM file with a proportional
 * duplication of reads, removing the XW tag
 *
 * Jeremy Wang
 * 20160513 (Friday)
*/


int main(int argc, char *argv[]) {

  if(argc < 4) {
    printf("Usage: expander <BAM> <weight multiplier> <output BAM>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  int multiplier = atoi(argv[2]);
  char *output_bam = argv[3];

  samFile *bam = sam_open(bam_file, "rb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\" for reading\n", bam_file);
    return -1;
  }
  bam_hdr_t *header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return -1;
  }
  bam1_t *aln = bam_init1();

  samFile *out_bam = sam_open(output_bam, "wb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\" for writing\n", output_bam);
    return -1;
  }

  int ret_val = sam_hdr_write(out_bam, header); // copy header

  float weight = 1.0;
  uint8_t *tag;
  int i;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    tag = bam_aux_get(aln, "XW");
    if (tag != 0) { // tag exists
      weight = bam_aux2f(tag);
      bam_aux_del(aln, tag);
    } else {
      // not XW tag
      weight = 1.0;
    }
    for(i = 0; i < weight * multiplier; i++) {
      ret_val = sam_write1(out_bam, header, aln);
    }
  }

  bam_destroy1(aln);
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

