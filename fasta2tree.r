#!/usr/bin/Rscript

library("ape");
library("seqinr");
library("cluster");

data.in <- read.dna("aligned_nonNImutsOnly_all_plus_NI.fa", format = "fasta");

data.mac <- apply(data.in,2,function(x){min(table(as.character(x)))});

select.rows <- sort(unique(c(sample(1:dim(data.in)[1],1000),grep("^(gi|NI)",rownames(data.in)))));
select.rows <- grep("^(gi|NI)",rownames(data.in));

select.cols <- ((data.mac > max(data.mac)/1000) & (data.mac < max(data.mac)/2));

data.sub <- data.in[select.rows,select.cols];

data.dist <- dist.dna(data.sub, model = "indelblock");

sub.tree <- nj(data.dist);

bs.tree <- boot.phylo(sub.tree, data.sub, function(x){ nj(dist.dna(x, model = "indelblock"))}, B = 100);

tipStarts <- substring(sub.tree$tip.label,1,2);
lookup <- sapply(list(gi = "red", NI = "blue"),paste);
tipColours <- lookup[tipStarts];

sub.tree$tip.label <- rep("o", length(sub.tree$tip.label));

pdf("out_dendrogram_NI_2014-Sep-11.pdf");
plot(sub.tree, type="fan", show.tip.label = TRUE, tip.color = tipColours);
graphics.off();
