#!/bin/sh
for x in $(cat lists/runlist.txt)
  do rpt=""
  if echo ${x} | grep "repeat" > /dev/null; then rpt=".1"; fi
  y1=$(echo ${x} | perl -pe 's/^.*?(repeat)?([rR]un_?)([0-9]+).*$/$1run_$3/')
  for f in $(ls /mnt/ihbi_ngs/iontorrent-kgq4/results/analysis/output/Home/${x}*/*.bam \
                /data/disk4/CurrentlySittingonIontorrent/analysis/output/Home/${x}*/*.bam | grep 'IonXpress')
    do y2=$(basename $f | perl -pe 's/^.*(IonXpress_[0-9]+).*$/$1/')
    if [ -e ${y1}${rpt}_${y2}.fastq.gz ];
       then echo -n ".";
    else
       echo "$f -> ${y1}${rpt}_${y2}.fastq.gz";
       pv ${f} | samtools view - | awk -F '\t' '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > ${y1}${rpt}_${y2}.fastq.gz;
    fi
  done
done
