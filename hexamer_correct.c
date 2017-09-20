#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"
#include "htslib/htslib/ksort.h"
#include "htslib/htslib/kseq.h"

/*
 * hc.c
 * 
 * HC: Hexamer bias correction
 * 1. Compute observed hexamer frequencies
 * 2. Compute background (expected) hexamer frequencies
 * 3. Write a new BAM file with reads w/XW(f) tag weight values
 *
 * Jeremy Wang
 * 20160429
 *
 * We use the fast 2-bit encoding (ascii >> 1 & 3):
 * Aa   0
 * Cc   1
 * GgNn 3
 * Tt   2
 * (val ^ 2) gives reverse-complement
 * - of course, Ns are indistinguishable from Gs, and look like Cs on the reverse strand
*/

typedef struct hexWeight{
  uint16_t hex;
  float weight;
} hexWeight;

// define sort of hexWeight by weight
#define weight_lt(a, b) ((a).weight < (b).weight)
KSORT_INIT(by_weight, hexWeight, weight_lt)

// fix my endianness issues, should also work on big endian architectures, but have not tested
const int ONE = 1;
#define is_bigendian() ( (*(char*)&ONE) == 0 )
#define reverse_int(a) ( ((a & 255) << 24) + (((a >> 8) & 255) << 16) + (((a >> 16) & 255) << 8) + ((a >> 24) & 255) )

// convert 2-bit array (a), ref sequence (b), of alleles to uint32_t encoded hexamer at position (c)
#define tohex(a, b, c) ( ( (is_bigendian() ? (*((uint32_t*)(a[b]+((c)/4)))) : reverse_int(*((uint32_t*)(a[b]+((c)/4))))) >> (2*(10-((c)%4))) ) & 4095 )
// reverse complement 12-bit hexamer order (2730 === 1010 10101010)
#define rc(a) ( (((a & 3) << 10) + (((a >> 2) & 3) << 8) + (((a >> 4) & 3) << 6) + (((a >> 6) & 3) << 4) + (((a >> 8) & 3) << 2) + ((a >> 10) & 3)) ^ (2730) )

const int32_t BACKGROUND_MARGIN = 20;

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

char* hex2str(uint16_t hex) {
  char* s = malloc(sizeof(char) * 7);
  s[6] = '\0';
  int i;
  for(i = 0; i < 6; i++) {
    switch((hex >> ((5-i)*2)) & 3) {
      case 0:
        s[i] = 'A';
        break;
      case 1:
        s[i] = 'C';
        break;
      case 2:
        s[i] = 'T';
        break;
      case 3:
        s[i] = 'G';
    }
  }
  return s;
}


char* u32bin(uint32_t d) {
  int i;
  char *c = malloc(sizeof(char) * 33);
  c[32] = '\0';
  for(i = 0; i < 32; i++) {
    c[i] = (d & (2<<(31-i))) ? '1' : '0';
  }
  return c;
}
char* u16bin(uint16_t d) {
  int i;
  char *c = malloc(sizeof(char) * 17);
  c[16] = '\0';
  for(i = 0; i < 16; i++) {
    c[i] = (d & (2<<(15-i))) ? '1' : '0';
  }
  return c;
}
char* u8bin(uint8_t d) {
  int i;
  char *c = malloc(sizeof(char) * 9);
  c[8] = '\0';
  for(i = 0; i < 8; i++) {
    c[i] = (d & (2<<(7-i))) ? '1' : '0';
  }
  return c;
}


int main(int argc, char *argv[]) {

  if(argc < 4) {
    printf("Usage: hc <BAM> <reference FASTA> <output BAM>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *ref_fasta = argv[2];
  char *output_bam = argv[3];


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


  // load BAM file
  // compute observed and background hexamer frequencies


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

  // initialize read and background hexamer counts
  uint32_t hex_counts[4096];
  uint32_t background_counts[4096];
  for (i = 0; i < 4096; i++) {
    hex_counts[i] = 0;
    background_counts[i] = 0;
  }

  int32_t read_len = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;

    if (ct == 0) {
      read_len = qlen;
    } else if (qlen != read_len) {
      // reads are not all the same length...
    }

    if (aln->core.flag & 4) { // unmapped
      continue;
    }
    uint32_t hex, al;
    if (!bam_is_rev(aln)) {
      hex = tohex(ref_array, tid, pos-3); // hexamer is centered about cut site
      //printf("chrom %i, pos %i, offset: %i, ref bytes: %u, refbin: %s, bin: %u, hexbin: %s\n", tid, pos, pos%4, *((uint32_t*)(ref_array[tid]+(pos/4))), u32bin(*((uint32_t*)(ref_array[tid]+(pos/4)))), hex, u16bin((uint16_t)hex));
      hex_counts[hex]++;
      // background centered around cut site
      int bg_range_start = (pos >= BACKGROUND_MARGIN ? -1 * BACKGROUND_MARGIN : -1*pos);
      int bg_range_end = (pos + BACKGROUND_MARGIN < rlen_array[tid] ? BACKGROUND_MARGIN : rlen_array[tid] - pos - 1);
      for (i = bg_range_start; i < bg_range_end + 1; i++) {
        hex = tohex(ref_array, tid, pos + i);
        background_counts[hex]++;
      }
    } else {
      hex = tohex(ref_array, tid, endpos-2); // hexamer is centered about cut site
      hex = rc(hex);
      //printf("RC chrom %i, pos %i, ref bytes: %u, hex: %u, bin: %s\n", tid, pos, *((uint32_t*)(ref_array[tid]+(pos/4))), hex, u32bin(hex));
      hex_counts[hex]++;
      // background
      int bg_range_start = (endpos + BACKGROUND_MARGIN + 1 < rlen_array[tid] ? -1 * (BACKGROUND_MARGIN + 1) : endpos - (rlen_array[tid] - 1));
      int bg_range_end = (endpos >= BACKGROUND_MARGIN ? BACKGROUND_MARGIN : endpos);
      for (i = bg_range_start; i < bg_range_end + 1; i++) {
        hex = tohex(ref_array, tid, endpos - i);
        hex = rc(hex);
        background_counts[hex]++;
      }
    }

    ct++;
    /*
    if (ct > 10) {
      break;
    }
    */
  }

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }


  // compute hexamer weights
  float weights[4096];
  uint64_t hex_total = 0;
  uint64_t background_total = 0;
  for (i = 0; i < 4096; i++) {
    hex_total = hex_total + hex_counts[i];
    background_total = background_total + background_counts[i];
  }

  hexWeight hex_weights[4096];
  for (i = 0; i < 4096; i++) {
    weights[i] = ((float)background_counts[i] / background_total) / ((float)hex_counts[i] / hex_total);
    hex_weights[i].hex = i;
    hex_weights[i].weight = weights[i];
  }

  // sort by weight in place
  ks_mergesort(by_weight, 4096, hex_weights, 0);

  printf("Lowest weighted hexamers:\n");
  for(i = 0; i < 10; i++) {
    printf("%s: %f\n", hex2str(hex_weights[i].hex), hex_weights[i].weight);
  }
  printf("Highest weighted hexamers:\n");
  for(i = 4085; i < 4096; i++) {
    printf("%s: %f\n", hex2str(hex_weights[i].hex), hex_weights[i].weight);
  }



  // read BAM file again
  // apply weights and write to new BAM file

  bam = sam_open(bam_file, "rb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\" for reading\n", bam_file);
    return -1;
  }
  header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return -1;
  }

  samFile *out_bam = sam_open(output_bam, "wb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\" for writing\n", output_bam);
    return -1;
  }
  ret_val = sam_hdr_write(out_bam, header); // copy header

  // read input, assign weight, write output
  float weight;
  float unmapped_weight = 1.0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {
    int32_t tid = aln->core.tid;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;

    // in the off chance that they're already using the XW tag, just delete it :/
    uint8_t *tag = bam_aux_get(aln, "XW");
    if (tag != 0) {
      bam_aux_del(aln, tag);
    }

    // htslib supports 'f' (4-byte) and 'd' (8-byte) types
    // internally, aux_type2size('f') would return 4, hopefully sizeof(float) does the same
    // use bam_aux2f(bam_aux_get(aln, "XW")) to decode
    if (aln->core.flag & 4) { // unmapped, weight will be 1
      bam_aux_append(aln, "XW", 'f', sizeof(float), (uint8_t*)&unmapped_weight);
    } else {
      if (!bam_is_rev(aln)) { // fw
        bam_aux_append(aln, "XW", 'f', sizeof(float), (uint8_t*)&weights[tohex(ref_array, tid, pos-3)]);
      } else { // rv
        uint32_t hex = tohex(ref_array, tid, endpos-2);
        hex = rc(hex);
        bam_aux_append(aln, "XW", 'f', sizeof(float), (uint8_t*)&weights[hex]);
      }
    }
    ret_val = sam_write1(out_bam, header, aln);
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

