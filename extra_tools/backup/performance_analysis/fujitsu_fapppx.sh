#!/bin/sh

rm -rf proftxt
mkdir proftxt

fapppx -A -d Fprofd_Stati -Impi,hwm -p0,limit=4 -o proftxt/prof_stati.txt
fapppx -A -d Fprofd_SIMD  -Impi,hwm -p0,limit=4 -o proftxt/prof_simd.txt
fapppx -A -d Fprofd_MEM   -Impi,hwm -p0,limit=4 -o proftxt/prof_mem.txt
