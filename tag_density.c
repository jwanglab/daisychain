#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <htslib/sam.h>
#include <htslib/khash.h>

/*
 * tagdensity.c
 * 
 * Jeremy Wang
 * 20160519
 *
 * Compute tag density scores
 * a la footprinting
 * given a BAM file and motif BED file
*/

// creates string:int hash
KHASH_MAP_INIT_STR(tidMap, int);

// creates string:uint32 hash
KHASH_MAP_INIT_STR(refLen, uint32_t);

int32_t MARGIN_5P = 99; // bp in 5' direction from the center of a motif
int32_t MARGIN_3P = 100; // bp in 3' direction from the center of a motif

int main(int argc, char *argv[]) {

  if(argc < 4) {
    printf("Usage: tagdensity <BAM> <motif BED> <output CSV>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *motif_bed = argv[2];
  char *output_csv = argv[3];

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

  khash_t(tidMap) *rmap = kh_init(tidMap);
  khash_t(refLen) *rlen = kh_init(refLen);

  for (i = 0; i < header->n_targets; i++) {
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);
    bin = kh_put(tidMap, rmap, strdup(header->target_name[i]), &absent);
    kh_val(rmap, bin) = i;
    bin = kh_put(refLen, rlen, strdup(header->target_name[i]), &absent);
    kh_val(rlen, bin) = header->target_len[i];
  }

  // API scaffolding largely from htslib/test/test_view.c
  hts_idx_t *idx;
  // searches for bam_file + .bai, etc
  if ((idx = sam_index_load2(bam, bam_file, NULL)) == 0) {
    fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
    return 1;
  }

  aln = bam_init1();

  hts_itr_t *iter;

  FILE* fp;
  fp = fopen(motif_bed, "r");
  char* line = NULL;
  size_t len = 0;
  ssize_t chr_read;

  FILE* fout;
  fout = fopen(output_csv, "w");

  // all of these must be allocated or sscanf will segfault
  char chrom[100], motif[100];
  char strand = '*';
  int start = 0, end = 0;
  float score = 0, pval = 0;

  // read BED file line by line
  while ((chr_read = getline(&line, &len, fp)) != -1) {

    // chrom  start  end  motif  score?  strand  pval?
    // chr1	29232	29246	CTCF_disc8_CTCF	12.2347	-	2.25e-05
    ret_val = sscanf(line, "%s\t%d\t%d\t%s\t%f\t%c\t%f", chrom, &start, &end, motif, &score, &strand, &pval);
    if (ret_val == -1) {
      fprintf(stderr, "Error scanning line '%s'\n", line);
      return 1;
    }

    // look up tid by chromosome
    bin = kh_get(tidMap, rmap, chrom);
    int tid = kh_value(rmap, bin);

    int range_start = start + (end - start) / 2 - MARGIN_5P;
    int range_end = start + (end - start) / 2 + MARGIN_3P;
    //printf("Query range: %d: %d-%d\n", tid, range_start, range_end);
    if ((iter = sam_itr_queryi(idx, tid, range_start, range_end)) == 0) {
        fprintf(stderr, "[E::%s] fail to fetch region %i:%i-%i\n", __func__, tid, range_start, range_end);
        return 1;
    }

    float tag_count = 0;
    while ((ret_val = sam_itr_next(bam, iter, aln)) >= 0) {
      int32_t pos = aln->core.pos;
      int32_t endpos = bam_endpos(aln) - 1;
      /*
      printf("tid: %i, qlen: %i, pos: %i, endpos: %i, ", tid, qlen, pos, endpos);
      printf(bam_is_rev(aln) ? "rv" : "fw");
      printf("\n");
      */
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

      // all we assume is that the fetched reads *overlap* the target window
      if ((bam_is_rev(aln) && endpos <= range_end) || (!bam_is_rev(aln) && pos >= range_start)) {
        tag_count += weight;
      }
    }

    fprintf(fout, "%s,%d,%d,%c,%d\n", chrom, start, end, strand, (int)(tag_count + 1));

    ct++;
    /*
    if(ct >= 10) {
      break;
    }
    */
  }

  ret_val = fclose(fout);

  hts_itr_destroy(iter);

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  return 0;
}

