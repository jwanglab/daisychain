
     )
    (O
     )  Daisychain
    (
    O)
    (

A DNAse-seq Toolkit
===================


Tools
-----

### Utilities
-  BED filter
-  Allele frequency

### DNAse hypersensitivity
-  Feature density

### Footprinting
-  Tag density

    Ooooo
    (   )
     ) /
    (_/

### Bias correction
-  Hexamer corrector
-  Dynamic k-mer corrector
-  Weight "expander"


Dependencies
------------
Both should be cloned and built in *this* directory

http://github.com/lh3/klib

http://github.com/samtools/htslib


Installation
------------

    git clone http://github.com/txje/daisychain
    make

Right now, each tool compiles to a separate executable. Run it to get a *very* brief listing of its parameters.

I'll include more documentation and probably a central executable as I go along.
