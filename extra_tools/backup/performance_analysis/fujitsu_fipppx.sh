#!/bin/sh

rm -rf proftxt
mkdir proftxt

fipppx -A -d Fprofd_Stati -Icpu -p0,limit=4 -o proftxt/prof_cpu.txt
fipppx -A -d Fprofd_Stati -Ibalance -p0,limit=4 -o proftxt/prof_balance.txt
#fipppx -A -d Fprofd_Stati -Icall -p0,limit=4 -o proftxt/prof_call.txt
fipppx -A -d Fprofd_Stati -Ihwm -p0,limit=4 -o proftxt/prof_hwm.txt
fipppx -A -d Fprofd_Stati -Isrc:./src -p0,limit=4 -o proftxt/prof_src.txt
