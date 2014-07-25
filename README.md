![PANDASEQ-SAM](https://rawgithub.com/neufeld/pandaseq/master/pandaseq.svg)
============

PANDASEQ-SAM is a program to align Illumina reads, optionally with PCR primers embedded in the sequence, and reconstruct an overlapping sequence. This version works on SAM/BAM formatted files.

INSTALLATION
------------

[![Build Status](https://travis-ci.org/neufeld/pandaseq-sam.png?branch=master)](https://travis-ci.org/neufeld/pandaseq-sam)

A working version of PANDAseq and [libhts](https://github.com/samtools/htslib)  must first be installed, including the development headers.

On Ubuntu, this can be installed via

    sudo apt-get install libhts-dev pandaseq-dev

After the support packages are installed, one should be able to do:

    ./autogen.sh && ./configure && make && sudo make install

If you receive an error that `libpandaseq-sam.so.0` is not found, try running:

    sudo ldconfig

USAGE
-----

Please consult the manual page by invoking

    man pandaseq-sam

though this is only summarises the differences. Consult

    man pandaseq

for the details.

The short version is

    pandaseq-sam -f seq.sam

or

    pandaseq-sam -b -f seq.bam

