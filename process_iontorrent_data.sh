#!/bin/sh

mkdir -p results
for x in *IonXpress_0*.fastq.gz; do
    y=$(./scripts/xpress2UUID.r ${x});
    if [ -e results/${y}.vcf ]; then
        echo "VCF file results/${y}.vcf already exists. Yay!"
    else
        echo "'${x}' -> '${y}'";
        echo "%%%% Processing ${y} %%%%";
        echo "** mapping ${y} **";
        bowtie2 -x db/hs_mtDNA_ref --local --rg-id ${y} -U ${x} | \
            samtools view -S -b - | samtools sort - results/${y};
        samtools index results/${y}.bam;
        maxdepth=$(samtools depth results/${y}.bam | cut -f 3 | sort -rn | uniq | head -n 1);
        echo "** getting VCF for ${y} (max depth = ${maxdepth}) **";
        pv results/${y}.bam | samtools mpileup -L ${maxdepth} -uf db/hs_mtDNA_ref.fasta - | \
            bcftools view -v -g - > results/${y}.vcf
    fi
    echo "** extracting heteroploidy data for ${y} **";
    awk -F '\t' 'BEGIN{print "ID Pos Ref Alt DP HQRf HQRr HQAf HQAr";}
        {if(!(/^#/)){print FILENAME,$2,$4,$5,$8}}' results/${y}.vcf | \
      perl -pe 's/results.//;s/.vcf / /;s/ [^ ]*DP=/ /;s/;.*DP4=/ /;
        s/;.*$//;s/([ACGTN]+),/$1;/g;tr/ /,/' \
      > results/heteroploidy_${y}.csv;
    echo "** getting consensus FASTQ for ${y} **";
    perl scripts/vcf2fq.pl -Q 20 -L 20 -f db/hs_mtDNA_ref.fasta results/${y}.vcf > results/con_${y}.fastq
    echo "** converting ${y} consensus FASTQ to plain FASTA **";
    perl scripts/fastq2fasta.pl results/con_${y}.fastq | perl -pe 'tr/acgt/ACGT/' > results/con_${y}.fasta
    perl -i -pe "s/^>.*$/>${y}/" results/con_${y}.fasta;
done

echo "** hashing consensus sequences **";
(for y in $(./scripts/xpress2UUID.r *IonXpress_0*.fastq.gz); do
   if [ -e results/con_${y}.fasta ]; then
      grep -v '^>' results/con_${y}.fasta | sha1sum | perl -pe "s/-/${y}/";
   fi
done)> results/variant_hashes.txt

cat results/heteroploidy_NI*.csv | \
    sort -t ',' -k 1,1 -k 2,2n | uniq > results/heteroploidy_all.csv
ls results/heteroploidy_* | grep 'heteroploidy_run' >/dev/null && \
    cat results/heteroploidy_run*.csv | \
    sort -t ',' -k 1,1 -k 2,2n | uniq | grep -v '^ID' >> results/heteroploidy_all.csv

