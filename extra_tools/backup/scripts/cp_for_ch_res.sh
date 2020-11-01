#!/bin/sh

indir=test13_bkup
innum=001
outdir=test13
outnum=000

ls $indir/cnt/*.cnt.$innum > wk
cat wk | sed -e "s/$indir/cp $indir/" > wk1
cat wk | sed -e "s/$indir/ $outdir/" | sed -e "s/.cnt.$innum/.cnt.$outnum/" > wk2
echo "#!/bin/sh" > cp_for_ch_res.log
paste wk1 wk2 >> cp_for_ch_res.log

mkdir -p $outdir/cnt
chmod +x cp_for_ch_res.log
./cp_for_ch_res.log

rm wk wk1 wk2
