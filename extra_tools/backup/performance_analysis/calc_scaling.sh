#!/bin/sh

### subroutine ###
arrange_name(){
tail -n 58 $1/log/*.0000.log.001 > wrk.txt

awk '{name[NR]=$1}
     END{
       printf "# name";
         for(i= 3;i<=12;i++){printf " ["i-1"]"name[i]};
         for(i=16;i<=33;i++){printf " ["i-4"]"name[i]};
         for(i=37;i<=57;i++){printf " ["i-7"]"name[i]};
       printf "\n";
     }' wrk.txt >> scaling.dat

rm wrk.txt
}

arrange_time(){
tail -n 58 $1/log/*.000000.log.001 > wrk.txt

awk '{time[NR]=$3}
     END{
       printf " ";
         for(i= 3;i<=12;i++){printf "  "time[i]};
         for(i=16;i<=33;i++){printf "  "time[i]};
         for(i=37;i<=57;i++){printf "  "time[i]};
       printf "\n";
     }' wrk.txt >> scaling.dat

rm wrk.txt
}
### main ###
rm -f scaling.dat
touch scaling.dat

arrange_name "test01"

arrange_time "test01"
arrange_time "test02"
arrange_time "test03"


