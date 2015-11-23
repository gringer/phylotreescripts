#!/usr/bin/Rscript
library(Biostrings);

ref.fa <- readDNAStringSet("/data/all/david/mitochondria/db/hs_mtDNA_ref.fasta")[[1]];
data.df <- read.csv("/data/all/david/mitochondria/db/mtDNA_haplogroups.csv",
                    stringsAsFactors=FALSE);
data.df <- subset(data.df, mutations != "reserved");

## Tweaks to tree to make it more correct
data.df$parent[(data.df$level == 1) & (data.df$name != "mt-MRCA")] <- "mt-MRCA";
data.df$level[data.df$name == "mt-MRCA"] <- 0;

## get number of [new] variants for each haplogroup
data.df$ncounts <- sapply(sapply(data.df$mutations,strsplit,split=" "),length);

data.melted.df <- data.frame(name=rep(data.df$name,data.df$ncounts),
                             parent=rep(data.df$parent,data.df$ncounts),
                             variant=unlist(strsplit(data.df$mutations," ")));
data.melted.df <- subset(data.melted.df, variant != "reserved");
data.melted.df$weight <- ifelse(grepl("\\(",data.melted.df$variant),0.5,1);
data.melted.df$oldVariant <- data.melted.df$variant;
data.melted.df$variant <- chartr("acgt","ACGT",gsub("(\\(|\\))","",data.melted.df$variant));
data.melted.df$pos <- as.numeric(sub("^[^0-9]*([0-9]+).*$","\\1",data.melted.df$variant));

## look at INDELs
unique(data.melted.df[grep("d",data.melted.df$variant),"variant"]); # deletions
unique(data.melted.df[grep("\\.",data.melted.df$variant),"variant"]); # insertions

## Correct deletion variants to be of the form <ref><pos><var>
simpleDel.pos <- grepl("^[ACGT][0-9]+d$",data.melted.df$variant);
data.melted.df$variant[simpleDel.pos] <-
    paste0(sapply(data.melted.df$pos[simpleDel.pos],function(x){as.character(ref.fa[(x-1):(x)])}),
           as.numeric(sub("^[ACGT]([0-9]+)d$","\\1",data.melted.df$variant[simpleDel.pos]))-1,
           sapply(data.melted.df$pos[simpleDel.pos],function(x){as.character(ref.fa[(x-1):(x-1)])}));
data.melted.df$pos[simpleDel.pos] <- data.melted.df$pos[simpleDel.pos] - 1;
## special case variant correction for inconsistent back insertion notation
data.melted.df$variant <- sub("^([ACGT]?)([0-9]+)\\..*?([ACGT]?)d!","\\1\\3\\2\\1\\3",data.melted.df$variant);
## remaining deletion variants should be of the form <start>-<end>d
delPoss <- grep("d",data.melted.df$variant);
delVars <- sapply(strsplit(sub("d$","",data.melted.df$variant[delPoss]),"-"), as.numeric);
delVars[1,] <- delVars[1,]-1;
data.melted.df$variant[delPoss] <- paste0(apply(delVars,2,function(x){as.character(ref.fa[x[1]:x[2]])}),
                                          data.melted.df$pos[delPoss]-1,
                                          sapply(delVars[2,],function(x){as.character(ref.fa[x])}));
data.melted.df$pos[delPoss] <- data.melted.df$pos[delPoss] - 1;
## sanity check for deletions
head(data.melted.df[simpleDel.pos,]);
head(data.melted.df[delPoss,]);

## Correct insertion variants to be of the form <anc><pos><var>
## .X[ACGT] seems to be the same as .1[ACGT]
## .2[ACGT] seems to be a transcription error
data.melted.df$variant <- sub("\\.[X2]([ACGT])",".1\\1",data.melted.df$variant);
## Insertions appear to insert *after* the specified position (see 8289 variant)
insPoss <- grep("\\.",data.melted.df$variant);
data.melted.df$variant[insPoss] <-
    paste0(sapply(data.melted.df$pos[insPoss],function(x){as.character(ref.fa[x])}),
           data.melted.df$pos[insPoss],
           sapply(data.melted.df$pos[insPoss],function(x){as.character(ref.fa[x])}),
           sub("^.*?\\.1","",data.melted.df$variant[insPoss]));
## sanity check for insertions
head(data.melted.df[insPoss,]);

## correct back mutations to be <reference><pos><alternate>
data.melted.df$variant <- sub("!!","",data.melted.df$variant);
data.melted.df$variant <- sub("^([ACGT])([0-9]+)([ACGT])!$","\\3\\2\\3",data.melted.df$variant);

## left-normalise INDELs
## derived from pseudocode from http://genome.sph.umich.edu/wiki/Variant_Normalization
leftNorm <- function(ref.seq, ref, pos = NULL, alt = NULL){
    if(any(grepl("[0-9]",ref))){
        leftNorm(ref.seq=ref.seq,
                 ref=sub("^(.*?)([0-9]+)(.*)$","\\1",ref),
                 pos=as.numeric(sub("^(.*?)([0-9]+)(.*)$","\\2",ref)),
                 alt=sub("^(.*?)([0-9]+)(.*)$","\\3",ref));
    } else {
        equalRights <- ((substring(ref,nchar(ref)) == substring(alt,nchar(alt))) & ((nchar(ref) * nchar(alt)) != 1));
        while(sum(equalRights) > 0){
            substring(ref,nchar(ref))[equalRights] <- "";
            substring(alt,nchar(alt))[equalRights] <- "";
            emptyAlleles <- ((nchar(ref) == 0) | (nchar(alt) == 0));
            pos[emptyAlleles] <- pos[emptyAlleles] - 1;
            ref[emptyAlleles] <- paste0(substring(as.character(ref.seq),pos[emptyAlleles],pos[emptyAlleles]),
                                        ref[emptyAlleles]);
            alt[emptyAlleles] <- paste0(substring(as.character(ref.seq),pos[emptyAlleles],pos[emptyAlleles]),
                                        alt[emptyAlleles]);
            equalRights <- ((substring(ref,nchar(ref)) == substring(alt,nchar(alt))) & ((nchar(ref) * nchar(alt)) != 1));
        }
        leftChangeable <- (substring(ref,1,1) == substring(alt,1,1)) & (nchar(ref) >= 2) & (nchar(alt) >= 2);
        while(sum(leftChangeable) > 0){
            substring(ref[leftChangeable],1,1) <- "";
            substring(alt[leftChangeable],1,1) <- "";
            pos[leftChangeable] <- pos[leftChangeable] + 1;
            leftChangeable <- (substring(ref,1,1) == substring(alt,1,1)) & (nchar(ref) >= 2) & (nchar(alt) >= 2);
        }
        paste0(ref,pos,alt);
    }
}

## sanity checks for leftNorm function
all(leftNorm("GGGCACACACAGGG",c("CA8","CAC6C","GCACA3GCA","GGCA2GG","GCA3G")) == "GCA3G");
leftNorm(ref.fa,c("ACCCCCTCTA8280A")) == "CACCCCCTCT8270C";
leftNorm(ref.fa,"G185G") == "G185G";

## do the actual left normalisation
data.melted.df$variant <- leftNorm(ref.fa,data.melted.df$variant);

## extract out variant characteristics
data.melted.df$ref <- sub("^(.*?)([0-9]+)(.*)$","\\1",data.melted.df$variant);
data.melted.df$pos <- as.numeric(sub("^(.*?)([0-9]+)(.*)$","\\2",data.melted.df$variant));
data.melted.df$alt <- sub("^(.*?)([0-9]+)(.*)$","\\3",data.melted.df$variant);

data.melted.df$vcfOrder <- sprintf("%05d;%s;%s",data.melted.df$pos,data.melted.df$ref,data.melted.df$alt);

## look at deletions
unique(grep("[ACGT][ACGT][0-9]+[ACGT]",data.melted.df$variant, value=TRUE))
## look at insertions
unique(grep("[ACGT][0-9]+[ACGT][ACGT]",data.melted.df$variant, value=TRUE))


## generate unadjusted count matrix (ignoring hierarchy)
hap.unadj.mat <- unclass(xtabs(weight ~ name + vcfOrder, data = data.melted.df));
hap.unadj.mat <- rbind(hap.unadj.mat,"mt-MRCA"=rep(0,ncol(hap.unadj.mat)));

hapParents <- data.df$parent;
names(hapParents) <- data.df$name;

getChain <- function(name, lookup){
    if(is.na(lookup[name])){
        name;
    } else {
        c(name, getChain(lookup[name], lookup));
    }
}

hap.unadj.mat[1:10,1:10];

## generate adjusted count matrix (accounting for hierarchy)
## [takes about 2 mins]
hap.adj.mat <- hap.unadj.mat * 0;
hapGroupCounts <- hap.unadj.mat[,1] * 0;
names(hapGroupCounts) <- rownames(hap.unadj.mat);
mutCounts <- hap.unadj.mat[1,] * 0;
names(mutCounts) <- colnames(hap.unadj.mat);
for(hap in rownames(hap.unadj.mat)){
    hapChain <- getChain(hap,hapParents);
    hap.adj.mat[hapChain,] <- t(t(hap.adj.mat[hapChain,,drop=FALSE]) + hap.unadj.mat[hap,]);
    hapGroupCounts[hapChain] <- hapGroupCounts[hapChain] + 1;
    mutCounts <- mutCounts + hap.unadj.mat[hap,];
    cat(".");
}

## sanity check on leaf node
hap <- "A17";
a <- hap.adj.mat[getChain(hap,hapParents),];
rowSums(a);
## make sure all haplotypes have at least one distinguishing mutation
min(rowSums(hap.adj.mat));

## p(mutation | haplogroup)
hap.adj.rowprop <- hap.adj.mat / rowSums(hap.adj.mat);

## note: p(haplogroup | mutation) probably isn't reliable, as most
## novel haplogroups are not discovered by a hypothesis-free
## approach. That's okay, as Bayes theorem can be used to approximate
## this (and work out how reliable it is)
#### not working...
####hap.adj.colprop <- t(t(hap.adj.mat) / colSums(hap.adj.mat));

## p(haplogroup)
hap.prop <- hapGroupCounts / hapGroupCounts["mt-MRCA"];

## p(mutation)
mut.prop <- mutCounts / sum(mutCounts);

