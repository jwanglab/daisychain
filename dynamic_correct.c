#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>

/*
 * dc.c
 * 
 * Dynamic (k-mer) bias correction
 * 1. Compute observed k-mer frequencies
 * 2. Compute background (expected) k-mer frequencies
 * 3. Write a new BAM file with reads w/XW(f) tag weight values
 *
 * Jeremy Wang
 * 20160525
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

// fix my endianness issues, should also work on big endian architectures, but have not tested
const int ONE = 1;
#define is_bigendian() ( (*(char*)&ONE) == 0 )
#define reverse_int(a) ( ((a & 255) << 24) + (((a >> 8) & 255) << 16) + (((a >> 16) & 255) << 8) + ((a >> 24) & 255) )

// convert 2-bit array (a), ref sequence (b), of alleles to uint32_t encoded hexamer at position (c)
#define tohex(a, b, c) ( ( (is_bigendian() ? (*((uint32_t*)(a[b]+((c)/4)))) : reverse_int(*((uint32_t*)(a[b]+((c)/4))))) >> (2*(10-((c)%4))) ) & 4095 )
// reverse complement 12-bit hexamer order (2730 === 1010 10101010)
#define rc(a) ( (((a & 3) << 10) + (((a >> 2) & 3) << 8) + (((a >> 4) & 3) << 6) + (((a >> 6) & 3) << 4) + (((a >> 8) & 3) << 2) + ((a >> 10) & 3)) ^ (2730) )

const int32_t BACKGROUND_MARGIN = 20;
const int32_t READ_MARGIN = 10;

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


float compute_weight(uint8_t **ref_array, int tid, int pos, float **weights, int32_t read_len, char rev) {
  float weight = 1;
  int i;
  for (i = -1 * READ_MARGIN; i < read_len + READ_MARGIN; i++) {
    if (!rev) {
      weight *= weights[i + READ_MARGIN][tohex(ref_array, tid, pos + i)];
    } else {
      weight *= weights[i + READ_MARGIN][rc(tohex(ref_array, tid, pos - i))];
    }
  }
  return weight;
}


int main(int argc, char *argv[]) {

  if(argc < 4) {
    printf("Usage: dc <BAM> <reference FASTA> <output BAM>\n");
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
  int l, i, j, absent;

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


  /*
   * Load BAM file
   * Compute observed and background hexamer frequencies
   */

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
  aln = bam_init1();

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

  // get the first read to get read length
  ret_val = sam_read1(bam, header, aln);
  int32_t read_len = aln->core.l_qseq;
  size_t bias_window_size = read_len + READ_MARGIN * 2;

  printf("Read length: %i\n", read_len);
  printf("Evaluating window of size: %i\n", bias_window_size);


  /*
   * Main loop to iteratively compute and apply weights
   * - until there are no severely biased kmers
   */

  float **weights = malloc(sizeof(float*) * bias_window_size);
  for (j = 0; j < bias_window_size; j++) {
    weights[j] = malloc(sizeof(float) * 4096);
    for (i = 0; i < 4096; i++) {
      weights[j][i] = 1.0;
    }
  }

  double background_counts[4096]; // only one background array
  for (i = 0; i < 4096; i++) {
    background_counts[i] = 0;
  }
  double baseline_allele_counts[4];
  uint64_t background_total = 0;

  int iter = 0; // loop counter

  while(1) {

    if (iter > 0) {
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
      aln = bam_init1();
      ret_val = sam_read1(bam, header, aln);
    }

    // initialize read hexamer counts
    double **hex_counts; // dynamically allocated array (one for each position in read) of 4096-arrays (one for each hexamer)
    hex_counts = malloc(sizeof(double*) * bias_window_size);
    for (i = 0; i < bias_window_size; i++) {
      hex_counts[i] = calloc(4096, sizeof(double)); // calloc allocates and sets to 0
    }

    // keep track of individual allele counts
    double** allele_counts = malloc(sizeof(double*) * bias_window_size);
    for (i = 0; i < 1000; i++) {
      allele_counts[i] = calloc(4, sizeof(double));
    }
    uint8_t al;

    uint32_t ct = 1;

    // keep track of how many consecutive reads are at this position
    // to try to detect alignment or PCR artifacts
    // - this will not work if the BAM file is unsorted
    int32_t current_tid = -1;
    int32_t current_pos = -1;
    uint32_t current_count = 0;

    printf("Looping through BAM file (iteration %i).\n", iter);
    
    // loop through BAM file
    do {
      int32_t tid = aln->core.tid;
      int32_t qlen = aln->core.l_qseq;
      int32_t pos = aln->core.pos;
      int32_t endpos = bam_endpos(aln) - 1;

      // read is not the same length
      if (qlen != read_len) {
        fprintf(stderr, "Read at %i:%i is %ibp, should be %ibp\n", tid, pos, qlen, read_len);
        continue;
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

      uint32_t hex, al;
      if (!bam_is_rev(aln)) {
        float weight = compute_weight(ref_array, tid, pos, weights, read_len, 0);

        // this loop could be more efficient if it kept the hex state and not recompute it every time
        int range_start = (pos >= READ_MARGIN ? pos - READ_MARGIN : 0);
        int range_end = (pos + read_len + READ_MARGIN - 1 < rlen_array[tid] ? pos + read_len + READ_MARGIN - 1 : rlen_array[tid]);
        for (i = range_start; i <= range_end; i++) {
          //printf("pos %i\n", i);
          hex = tohex(ref_array, tid, i);
          hex_counts[i - range_start][hex] += weight;

          //printf("chrom %i, pos %i, offset: %i, ref bytes: %u, refbin: %s, bin: %u, hexbin: %s\n", tid, pos, pos%4, *((uint32_t*)(ref_array[tid]+(pos/4))), u32bin(*((uint32_t*)(ref_array[tid]+(pos/4)))), hex, u16bin((uint16_t)hex));

          if (iter == 0) {
            // all go into a single background count
            background_counts[hex] += weight;
          }

          // per-locus allele counts
          al = (ref_array[tid][i / 4] >> (2 * (3 - (i % 4)))) & 3;
          allele_counts[i - range_start][al] += weight;
        }
      } else {
        float weight = compute_weight(ref_array, tid, endpos-5, weights, read_len, 1);

        // this loop could be more efficient if it kept the hex state and not recompute and RC it every time
        int range_start = (endpos - read_len + 1 >= READ_MARGIN ? endpos - read_len + 1 - READ_MARGIN : 0);
        int range_end = (endpos + READ_MARGIN < rlen_array[tid] ? endpos + READ_MARGIN : rlen_array[tid] - 1);
        for (i = range_start; i <= range_end; i++) {
          //printf("range_start %i, current(i) %i, range_end %i\n", range_start, i, range_end);

          // the canonical position relative to the bias window if this were on the forward strand
          // i.e. the position in the hexamer and allele frequency array
          int strand_pos = range_end - i;

          hex = tohex(ref_array, tid, i - 5);
          hex = rc(hex);
          hex_counts[strand_pos][hex] += weight;

          //printf("RC chrom %i, pos %i, ref bytes: %u, hex: %u, bin: %s\n", tid, pos, *((uint32_t*)(ref_array[tid]+(pos/4))), hex, u32bin(hex));

          if (iter == 0) {
            background_counts[hex] += weight;
          }

          al = (ref_array[tid][i / 4] >> (2 * (3 - (i % 4)))) & 3;
          allele_counts[strand_pos][al ^ 2] += weight; // rv
        }
      }

      ct++;
      /*
      if (ct > 10) {
        break;
      }
      */
    } while ((ret_val = sam_read1(bam, header, aln)) >= 0);


    /*
     * Clean up this pass through the BAM file
     */

    ret_val = sam_close(bam);
    if (ret_val < 0) {
      fprintf(stderr, "Error closing input BAM.\n");
      return -1;
    }
    bam_hdr_destroy(header);


    /*
     * Compute hexamer weights
     */
    uint64_t *hex_total = calloc(bias_window_size, sizeof(uint64_t));
    for (i = 0; i < 4096; i++) {
      for (j = 0; j < bias_window_size; j++) {
        hex_total[j] = hex_total[j] + hex_counts[j][i];
      }
      if(iter == 0) {
        background_total = background_total + background_counts[i];
      }
    }


    /*
     * compute allelic variance per hexamer window
     */
    // baseline will be computed the first time through then stay the same
    if (iter == 0) {
      for (i = 0; i < 4; i++) baseline_allele_counts[i] = 0;
      for (i = 0; i < 4096; i++) {
        baseline_allele_counts[(i >> 10) & 3] += background_counts[i];
      }
      for (i = 0; i < 4; i++) baseline_allele_counts[i] /= background_total;

      printf("Baseline\n");
      printf("A: %f\n", baseline_allele_counts[0]);
      printf("C: %f\n", baseline_allele_counts[1]);
      printf("G: %f\n", baseline_allele_counts[3]);
      printf("T: %f\n", baseline_allele_counts[2]);
    }

    for (i = 0; i < bias_window_size; i++) {
      double tot_alleles = allele_counts[i][0] + allele_counts[i][1] + allele_counts[i][2] + allele_counts[i][3];
      for (j = 0; j < 4; j++) {
        allele_counts[i][j] /= tot_alleles;
      }
    }

    double *avg_variance = calloc(bias_window_size, sizeof(double));
    for (i = 0; i < bias_window_size; i++) {
      int last = bias_window_size - i >= 6 ? 6 : bias_window_size - i;
      for (j = 0; j < last; j++) {
        int l;
        for (l = 0; l < 4; l++) {
          //printf("At pos %i + %i allele %i, obs %f, baseline %f\n", i, j, l, allele_counts[i + j][l], baseline_allele_counts[l]);
          double diff = allele_counts[i + j][l] - baseline_allele_counts[l];
          avg_variance[i] += (diff < 0 ? -1 * diff : diff) / baseline_allele_counts[l];
        }
      }
      avg_variance[i] /= 4 * last;

      //printf("Avg variance at locus %i: %f\n", i, avg_variance[i]);
    }


    /*
     * Find tile of maximum variance
     */

    int max_variance_tile = 0;
    for (i = 1; i < bias_window_size; i++)
      if(avg_variance[i] > avg_variance[max_variance_tile])
        max_variance_tile = i;

    printf("Most biased tile is %i, compounding the reweighting for this tile\n", max_variance_tile);

    // apply the adjustment for the most biased tile
    for (i = 0; i < 4096; i++) {
      weights[max_variance_tile][i] *= ((float)background_counts[i] / background_total) / ((float)hex_counts[max_variance_tile][i] / hex_total[max_variance_tile]);
    }

    // do I need to free the second dimension first?
    //free(hex_counts);
    //free(allele_counts);

    iter++;
    if(iter >= 4) {
      break;
    }
  } // end of main loop


  /*
   * Compute and write final weights to output BAM file
   */

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
  aln = bam_init1();

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
        weight = compute_weight(ref_array, tid, pos, weights, read_len, 0);
      } else { // rv
        weight = compute_weight(ref_array, tid, endpos-5, weights, read_len, 1);
      }
      bam_aux_append(aln, "XW", 'f', sizeof(float), (uint8_t*)&weight);
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

