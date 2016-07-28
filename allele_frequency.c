#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/ksort.h>
#include <htslib/kseq.h>

/*
 * allele_frequencies.c
 * 
 * Jeremy Wang
 * 20160502
 *
 * We use the fast 2-bit encoding (ascii >> 1 & 3):
 * Aa   0
 * Cc   1
 * GgNn 3
 * Tt   2
 * (val ^ 2) gives reverse-complement
 * - of course, Ns are indistinguishable from Gs, and look like Cs on the reverse strand
*/

// warn if the coverage at a single locus is higher than this
uint32_t COVERAGE_WARNING_LIMIT = 200;

int32_t MARGIN = 10;

// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(FILE*, fileread);

// creates string:[array of uint8] hash
KHASH_MAP_INIT_STR(refSeq, uint8_t*);

// creates string:uint32 hash
KHASH_MAP_INIT_STR(refLen, uint32_t);


int main(int argc, char *argv[]) {

  if(argc < 3) {
    printf("Usage: allele_frequencies <BAM> <reference FASTA>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *ref_fasta = argv[2];


  // load FASTA file
  // convert into fw and rv arrays of 2-byte hexamer values

  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);

  FILE* fp;
  kseq_t* seq;
  int l, i, absent;

  fp = fopen(ref_fasta, "r");
  seq = kseq_init(fp);
  printf("Reading fasta file: %s\n", ref_fasta);

  khint_t bin; // hash bin (result of kh_put)

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    printf("Reading %s (%i bp).\n", seq->name.s, l);

    // pack 4 alleles into a byte (2 bits each)
    int n_bytes = (l/4) + (l%4 == 0 ? 0 : 1);
    uint8_t* bytes = malloc(sizeof(uint8_t) * n_bytes);
    for(i = 0; i < n_bytes; i++) bytes[i] = 0; // init

    uint8_t allele;
    for(i = 0; i < l; i++) {
      allele = (seq->seq.s[i] >> 1) & 3;
      bytes[i/4] += (allele << (2 * (3 - (i % 4))));
    }

    // seq array
    bin = kh_put(refSeq, ref, strdup(seq->name.s), &absent);
    kh_val(ref, bin) = bytes;

    // sequence length
    bin = kh_put(refLen, rlen, strdup(seq->name.s), &absent);
    kh_val(rlen, bin) = l;
  }

  kseq_destroy(seq);
  fclose(fp);


  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int ret_val;

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
  // construct array from reference information so that we can look it up with read.tid
  uint8_t **ref_array = malloc(sizeof(uint8_t*) * header->n_targets);
  int *rlen_array = malloc(sizeof(uint32_t) * header->n_targets); // has to be a plain int because that's what kseq gives out
  printf("BAM targets:\n");
  for (i = 0; i < header->n_targets; i++) {
    printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);
    bin = kh_get(refSeq, ref, header->target_name[i]);
    ref_array[i] = kh_value(ref, bin);
    bin = kh_get(refLen, rlen, header->target_name[i]);
    int fa_len = kh_value(rlen, bin);
    if (fa_len != header->target_len[i]) { // target_len is a uint32_t
      printf("Reference fasta length (%i) and BAM length (%u) of %s do not agree.\n", fa_len, header->target_len[i], header->target_name[i]);
    }
    rlen_array[i] = fa_len;
  }

  aln = bam_init1();

  int32_t read_len = 0;
  double* allele_counts;
  allele_counts = calloc(4 * 10000, sizeof(double));
  uint8_t al;

  // keep track of how many consecutive reads are at this position
  // to try to detect alignment or PCR artifacts
  // - this will not work if the BAM file is unsorted
  int32_t current_tid = -1;
  int32_t current_pos = -1;
  uint32_t current_count = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;
    /*
    printf("tid: %i, qlen: %i, pos: %i, endpos: %i, ", tid, qlen, pos, endpos);
    printf(bam_is_rev(aln) ? "rv" : "fw");
    printf("\n");
    */

    if (ct == 0) {
      read_len = qlen;
      /*
      allele_counts = malloc(sizeof(double) * 4 * (read_len + MARGIN*2));
      for (i = 0; i < 4 * (read_len + MARGIN*2); i++) { // init
        allele_counts[i] = 0.0;
      }
      */
    } else if (qlen != read_len) {
      // reads are not all the same length...
    }

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    // warn about too-high coverage, possible artifacts
    if(tid != current_tid || pos != current_pos) {
      if(current_count >= COVERAGE_WARNING_LIMIT) {
        printf("WARNING: There are %i reads aligned to %i:%i. These may be alignment or PCR artifacts.\n", current_count, current_tid, current_pos);
      }
      current_tid = tid;
      current_pos = pos;
      current_count = 0;
    }
    current_count++;
    // actually only evaluate reads up to the limit
    if(current_count >= COVERAGE_WARNING_LIMIT) {
      continue;
    }

    float weight = 1.0;
    uint8_t *tag = bam_aux_get(aln, "XW");
    if (tag != 0) { // tag exists
      weight = bam_aux2f(tag);
      if (weight == 0.0) { // some datasets have this tag set before I touch them, most (maybe all) are zero
        weight = 1.0;
      }
    }
    if (!bam_is_rev(aln)) {
      int bg_range_start = (pos >= 10 ? -10 : -1*pos);
      int bg_range_end = (endpos + 10 < rlen_array[tid] ? (endpos - pos) + 10 : rlen_array[tid] - pos - 1);
      for (i = bg_range_start; i < bg_range_end + 1; i++) {
        al = (ref_array[tid][(pos + i) / 4] >> (2 * (3 - ((pos + i) % 4)))) & 3;
        allele_counts[(i + 10) * 4 + al] += weight;
      }
    } else {
      int bg_range_start = (endpos + 10 < rlen_array[tid] ? -10 : endpos - (rlen_array[tid] - 1));
      int bg_range_end = (pos >= 10 ? (endpos - pos) + 10 : endpos);
      for (i = bg_range_start; i < bg_range_end + 1; i++) {
        al = (ref_array[tid][(endpos - i) / 4] >> (2 * (3 - ((endpos - i) % 4)))) & 3;
        allele_counts[(i + 10) * 4 + (al ^ 2)] += weight; // rv
      }
    }

    ct++;
    /*
    if(ct >= 10) {
      break;
    }
    */
  }

  for (i = 0; i < read_len + MARGIN*2; i++) {
    printf("%i\t%f\t%f\t%f\t%f\n", (i-10), allele_counts[i*4], allele_counts[i*4+1], allele_counts[i*4+2], allele_counts[i*4+3]);
  }

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  return 0;
}

