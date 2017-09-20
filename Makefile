CC=gcc
CFLAGS=-lhts -std=c99

OBJECTS = dc hc af td ex bf fd

all: $(OBJECTS)

dc: dynamic_correct.c
	$(CC) $(CFLAGS) dynamic_correct.c -o dc

hc: hexamer_correct.c
	$(CC) $(CFLAGS) hexamer_correct.c -o hc

af: allele_frequency.c
	$(CC) $(CFLAGS) allele_frequency.c -o af

td: tag_density.c
	$(CC) $(CFLAGS) tag_density.c -o td

ex: expand.c
	$(CC) $(CFLAGS) expand.c -o ex

bf: bed_filter.c
	$(CC) $(CFLAGS) bed_filter.c -o bf

fd: feature_density.c
	$(CC) $(CFLAGS) feature_density.c -o fd -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
