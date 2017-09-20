#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"
#include "htslib/htslib/ksort.h"
#include "klib/kvec.h"
#include <float.h>
#include <time.h>
#include <math.h>

/*
 * feature_density.c
 * 
 * Jeremy Wang
 * 20160711
 *
 * Compute feature density estimate, modeled off F-seq,
 * but not necessarily producing identical results
 * - we focus on DNase-seq and DNase hypersensitivity,
 *   but this feature density estimation can be used
 *   for assays including ChIP, FAIRE, and ATAC
*/

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#ifndef M_E
  #define M_E 2.718281828459045
#endif

double PI2 = sqrt(M_PI * 2);

typedef struct read {
  uint32_t pos;
  float weight;
} read;

typedef kvec_t(read) readList;

typedef struct norm {
  double mean;
  double std;
} norm;

// define sort of reads by pos
#define pos_lt(a, b) ((a).pos < (b).pos)
KSORT_INIT(by_pos, read, pos_lt)

// creates string:int hash
KHASH_MAP_INIT_STR(tidMap, int);

// creates string:uint32 hash
KHASH_MAP_INIT_STR(refLen, uint32_t);

double gauss(double x) {
  return pow(M_E, -(x * x)/2) / PI2;
}

/*
 * Precompute gaussian distribution across the given window
 * b: bandwidth
 * w: width
 * n: total number of reads
 */
double* precompute(int b, int w, int n) {
  double* distr = calloc(w+1, sizeof(double));
  for(int i = 0; i <= w; i++) {
    distr[i] = gauss(i / (double)b) * 20000000 / n;
  }
  return distr;
}

/*
 * x: position to compute PDF
 * b: feature bandwidth
 * w: compute window size
 * read_weights: *sorted* list of (position, read weight)
 * r: read index to start
 * n: length of above
 * read_len: length of reads, only applies to ChIP, all reads must be the same length
 */
double density_by_read_list(int x, double* distr, int b, int w, read* read_weights, int r, int n, char assay, int read_len) {

  double p = 0;
  int start = (x >= w ? x - w : 0);
  int end = (x + w <= read_weights[n-1].pos ? x + w : read_weights[n-1].pos); // inclusive

  int i = r;

  // collect reads in the target window, accumulating a gaussian pdf
  if(assay == 'd' || assay == 'a') {
    // have to give ourselves a margin of read_len before quitting because fw/rv strand start/end positions are out of order
    while(i < n && read_weights[i].pos + read_len <= end) {
      if(read_weights[i].pos < start || read_weights[i].pos > end) { // not in range
        continue;
      }
      if(abs(read_weights[i].pos - x) >= w+1) {
        printf("computing density at pos %i, window %i - %i\n", x, start, end);
        printf("last read was at pos %i\n", read_weights[i-1].pos);
        printf("adding read %i at pos %i, weight %f\n", i, read_weights[i].pos, read_weights[i].weight);
        printf("out of range of distr at %i\n", abs(read_weights[i].pos - x));
      }
      p += distr[abs(read_weights[i].pos - x)] * read_weights[i].weight;
      i++;
    }
  } else if(assay == 'c' || assay == 'f') {
    // these positions are still in order because they're all start on fw strand
    int pos;
    while(i < n && read_weights[i].pos <= end) {
      // this will be slower than others by a factor of the read length!
      for(pos = (read_weights[i].pos < x-w ? x-w : read_weights[i].pos); pos < read_weights[i].pos + read_len && pos <= x + w; pos++) {
        if(abs(pos - x) >= w+1) {
          printf("out of range of distr at %i", abs(pos - x));
        }
        p += distr[abs(pos - x)] * read_weights[i].weight;
      }
      i++;
    }
  }

  p /= b;
  return p;
}

/*
 * x: position to compute PDF
 * b: feature bandwidth
 * w: compute window size
 * coverage: chrom-sized array of float
 * n: length of coverage array (chrom length)
 */
double density_by_coverage(int x, double* distr, int b, int w, float *coverage, int n) {
  double p = 0;
  int start = (x >= w ? x - w : 0);
  int end = (x + w <= n-1 ? x + w : n-1); // inclusive

  int i;
  for(i = start; i <= end; i++) {
    if(coverage[i] > 0) {
      p += distr[abs(i-x)];
    }
  }

  p /= b;
  return p;
}


int compute_window(int b) { // b: feature bandwidth
  int i = 1400;
      
  double bw = (double)b;
  while(1) {
    double x = ++i / bw;
    if(gauss(x) < FLT_MIN) break; // FLT_MIN is probably 1e-37, from <float.h>
  }
  return (i-1);
}

norm norm_distr(double *arr, int n) {
  norm d;

  double tot = 0.0;
  for(int i = 0; i < n; i++) {
    tot += arr[i];
  }

  d.mean = tot/n;

  double sum = 0.0;
  for(int i= 0; i < n; ++i){
    sum += pow(arr[i] - d.mean, 2);
  }
  d.std = sqrt(sum / n);
  return d;
}


/*
 * Computes the background distribution assuming a random, uniform
 * read density
 *
 * size: total genome size
 * ncuts: total read weight
 * w: window size
 * b: feature bandwidth
 */
norm compute_background(double size, int ncuts, double* distr, int b, int w, char assay, int read_len) {
  
  int totalWindow = 1 + (int)(w * 2);     
  printf("Background window: %i\n", totalWindow);
  int cutDensity = (int)((ncuts / size) * totalWindow);
  printf("Cut density: %i\n", cutDensity);
  int thresholdIterations = 10000;
  double densities[thresholdIterations];
  read cuts[cutDensity];
  
  for(int i = 0; i < thresholdIterations; ++i){
    for(int j = 0; j < cutDensity; ++j) {
      cuts[j].pos = (int)(rand() % totalWindow);
      cuts[j].weight = 1;
    }
    ks_mergesort(by_pos, cutDensity, cuts, 0); // sort by pos
    densities[i] = density_by_read_list(totalWindow/2, distr, b, w, cuts, 0, cutDensity, assay, read_len);
  }
  
  return norm_distr(densities, thresholdIterations);
} 


int main(int argc, char *argv[]) {

  // seed random
  srand(time(NULL));

  if(argc < 6) {
    printf("Usage: feature_density <BAM> <feature length> <std threshold> <assay {d,a,c,f}> <output BED> [output WIG]\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  int feature_length = atoi(argv[2]); // feature length (default: 600)
  int b = feature_length / 6.0; // bandwidth - this is the way it's done in F-seq
  int std_threshold = atoi(argv[3]);
  char assay = argv[4][0];
  char *output_bed = argv[5];
  char *output_wig = NULL;
  if(argc >= 7) {
    output_wig = argv[6];
  }

  if(assay != 'd' && assay != 'c' && assay != 'a' && assay != 'f') {
    fprintf(stderr, "Unknown assay type '%c', should be one of 'd' (DNase), 'c' (ChIP), 'f' (FAIRE), or 'a' (ATAC)\n", assay);
    return -1;
  }

  int w = compute_window(b); // computed window size

  printf("Bandwidth (b): %i\n", b);
  printf("Window size (w): %i\n", w);


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

  // NOTE: we're going to precompute coverage across the genomes
  //   - on the other hand, we could just keep a list of reads and compute
  //     the coverage on the fly. IDK which way is faster
  readList *reads = malloc(sizeof(readList) * header->n_targets); // array of read vectors, one for each chromosome
  //float **coverage = malloc(sizeof(float*) * header->n_targets); // ~12Gb

  for (i = 0; i < header->n_targets; i++) {
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);

    genome_size += header->target_len[i];

    //coverage[i] = (float*)calloc(header->target_len[i], sizeof(float));
    printf("chrom %i\n", i);
    kv_init(reads[i]);

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

    read r;
    r.weight = weight;

    if(assay == 'd' || assay == 'a') {
      if (bam_is_rev(aln)) {
        //coverage[tid][endpos - offset] += weight;
        r.pos = endpos - offset;
      } else {
        //coverage[tid][pos + offset] += weight;
        r.pos = pos + offset;
      }
    } else if(assay == 'c' || assay == 'f') {
      r.pos = pos;
      /*
      for(i = pos; i <= endpos; i++) {
        coverage[tid][i] += weight;
      }
      */
    }

    kv_push(read, reads[tid], r); // these are no longer ordered by position! (b/c the reverse strand)
  }

  // clean up BAM file
  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing BAM file.\n");
  }

  double *distr = precompute(b, w, total_mapped);

  // compute background
  norm bg = compute_background(genome_size, total_mapped, distr, b, w, assay, read_len);
  printf("Background mean: %f\n", bg.mean);
  printf("Background std: %f\n", bg.std);

  FILE *bed_fout, *wig_fout;
  // output bed format may be nonstandard: "chrom startPos endPos avgDensityScore"
  bed_fout = fopen(output_bed, "w");
  if(output_wig != NULL) {
    wig_fout = fopen(output_wig, "w");
  }
  
  int32_t tid;
  double threshold = bg.mean + bg.std * std_threshold;
  double d, tot_density = 0.0;
  uint8_t high = 0; // flag indicating if we're in a region with density over the threshold
  int start;
  int r; // read index on a chromosome
  for(tid = 0; tid < header->n_targets; tid++) {
    r = 0; // start at the first read on this chromosome kv_A(reads[tid], r)
    if(output_wig != NULL) {
      fprintf(wig_fout, "fixedStep chrom=%s start=0 step=1\n", header->target_name[tid]);
    }
    printf("Chrom %s, %i bp...\n", header->target_name[tid], header->target_len[tid]);

    for(i = 0; i < header->target_len[tid]; i++) {
      if(i % 1000000 == 0) {
        printf("At pos %i\n", i);
      }

      // increment r until the read is within striking distance (read_len) of the target position
      // we have to leave a margin of read_len because some reads indicate their start and some their end (on the reverse strand)
      // so they aren't strictly ordered by position
      while(r < kv_size(reads[tid]) && kv_A(reads[tid], r).pos + read_len < i) {
        r++;
      }

      // compute density - this is the magic
      //d = density_by_coverage(i, distr, b, w, coverage[tid], header->target_len[tid]);
      d = density_by_read_list(i, distr, b, w, reads[tid].a, r, kv_size(reads[tid]), assay, read_len);

      if(output_wig != NULL) {
        fprintf(wig_fout, "%f\n", d);
      }
      if(d > threshold) {
        if(!high) {
          start = i;
          high = 1;
        }
        tot_density += d;
      } else if(high) {
        fprintf(bed_fout, "%s\t%i\t%i\t%f\n", header->target_name[tid], start, i, tot_density/(i - start)); // end position is EXclusive
        high = 0;
        tot_density = 0.0;
      }
    }
    break;
  }

  ret_val = fclose(bed_fout);
  if(output_wig != NULL) {
    ret_val = fclose(wig_fout);
  }

  bam_hdr_destroy(header);

  return 0;
}
