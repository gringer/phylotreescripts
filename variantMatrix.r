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
        c(getChain(lookup[name], lookup), name);
    }
}

hap.unadj.mat[1:10,1:10];

## generate adjusted count matrix (accounting for hierarchy)
## [takes about 2 mins]
## note: back mutations are not properly modelled
hap.adj.mat <- hap.unadj.mat * 0;
hapGroupCounts <- hap.unadj.mat[,1] * 0;
names(hapGroupCounts) <- rownames(hap.unadj.mat);
mutCounts <- hap.unadj.mat[1,] * 0;
names(mutCounts) <- colnames(hap.unadj.mat);
for(hap in rownames(hap.unadj.mat)){
    hapChain <- getChain(hap,hapParents);
    hap.adj.mat[hap,] <- colSums(hap.unadj.mat[hapChain,,drop=FALSE]);
    hapGroupCounts[hapChain] <- hapGroupCounts[hapChain] + 1;
    mutCounts <- mutCounts + colSums(hap.unadj.mat[hapChain,,drop=FALSE]);
    cat(".");
}

## order by row counts
hap.adj.mat <- hap.adj.mat[order(-rowSums(hap.adj.mat)),];

hap.adj.mat[1:10,1:10];

## sanity check on leaf node
hap <- "A17";
a <- hap.adj.mat[getChain(hap,hapParents),];
rowSums(a);
a[,colSums(a) > 0];
## make sure all haplotypes have at least one distinguishing mutation
min(rowSums(hap.adj.mat[-grep("mt-MRCA",rownames(hap.adj.mat)),]));

## p(mutation | haplogroup)
## -- I'm not quite sure what this means
##hap.adj.rowprop <- hap.adj.mat / rowSums(hap.adj.mat);

## p(haplogroup | mutation) probably isn't reliable, but this is an approximation
hap.adj.colprop <- t(t(hap.adj.mat) / colSums(hap.adj.mat));

## p(haplogroup)
##hap.prop <- hapGroupCounts / hapGroupCounts["mt-MRCA"];
## start with the simple definition
hap.prop <- 1 / length(hapGroupCounts);
hap.nums <- length(hapGroupCounts);

## p(mutation)
mut.prop <- mutCounts / sum(mutCounts);
mut.nums <- length(mutCounts);


## it is assumed (for now) that undeclared variants are not indicative of a haplotype

## okay, now to try it on KCCG data

inFile.KCCG <- "/data/all/david/mitochondria/bam_fromKCCG/OVLNormalised_STARout_KCCG_called.vcf.gz";
KCCG.vcf.df <- read.delim(inFile.KCCG, stringsAsFactors=FALSE,
                          skip=grep("^#CHROM",readLines(inFile.KCCG))-1);
colnames(KCCG.vcf.df)[1] <- "CHROM";
colnames(KCCG.vcf.df)[grep("^FR",colnames(KCCG.vcf.df))] <- sub("_.*$","",colnames(KCCG.vcf.df)[grep("^FR",colnames(KCCG.vcf.df))]);
frNames <- grep("^FR",colnames(KCCG.vcf.df), value=TRUE);
KCCG.vcf.df[,paste0(frNames,".GT")] <- sapply(KCCG.vcf.df[,frNames],sub,pattern=":.*$",replacement="");
KCCG.vcf.df[,paste0(frNames,".DPRr")] <- as.numeric(sapply(KCCG.vcf.df[,frNames],sub,pattern="^.*:(.*),.*$",replacement="\\1"));
KCCG.vcf.df[,paste0(frNames,".DPRa")] <- as.numeric(sapply(KCCG.vcf.df[,frNames],sub,pattern="^.*,",replacement=""));
KCCG.vcf.df[,"AAF"] <- rowSums(KCCG.vcf.df[,paste0(frNames,".DPRa")]) /
    (rowSums(KCCG.vcf.df[,paste0(frNames,".DPRa")]) + rowSums(KCCG.vcf.df[,paste0(frNames,".DPRr")]));
table(unlist(KCCG.vcf.df[,paste0(frNames,".GT")]))
head(unlist(KCCG.vcf.df[,paste0(frNames,".DPRr")]))
head(KCCG.vcf.df);
str(KCCG.vcf.df);
KCCG.vcf.df$refSig <- paste(sprintf("%05d",KCCG.vcf.df$POS),KCCG.vcf.df$REF,KCCG.vcf.df$REF,sep=";");
KCCG.vcf.df$altSig <- paste(sprintf("%05d",KCCG.vcf.df$POS),KCCG.vcf.df$REF,KCCG.vcf.df$ALT,sep=";");
KCCG.vcf.df[,paste0(frNames,".Sig")] <- sapply(KCCG.vcf.df[,paste0(frNames,".GT")],
                                               function(x){ifelse(x == "1/1",KCCG.vcf.df$altSig, KCCG.vcf.df$refSig)});
KCCG.vcf.df$cM <- 0;
KCCG.vcf.df[,paste0(frNames,".SigGT")] <- sapply(KCCG.vcf.df[,paste0(frNames,".GT")],
                                                 function(x){ifelse(x == "1/1"," C C",
                                                                    ifelse(x == "0/0"," A A",
                                                                           ifelse(x %in% c("0/1","1/0")," A C", " 0 0")))});
KCCG.vcf.df[["#CHROM"]] <- "26";
sum(KCCG.vcf.df$altSig %in% names(mut.prop));
(novel.vars.KCCG <- KCCG.vcf.df$altSig[!KCCG.vcf.df$altSig %in% names(mut.prop)]);
##write.table(KCCG.vcf.df[!KCCG.vcf.df$altSig %in% names(mut.prop),c("altSig","AAF")], row.names=FALSE, sep="\t",
##            quote=FALSE);

write.table(KCCG.vcf.df[,c("#CHROM","altSig","cM","POS",paste0(frNames,".SigGT"))], file="/data/all/david/mitochondria/KCCG_gt_matrix.tped",
            row.names=FALSE, col.names=TRUE, quote=FALSE);
write.table(cbind(frNames,1,0,0,0,0), file="/data/all/david/mitochondria/KCCG_gt_matrix.tfam", quote=FALSE,
            row.names=FALSE, col.names=FALSE);

test.vars <- KCCG.vcf.df[,paste0(frNames,".Sig")];

KCCG.haplotypes <- apply(test.vars,2,function(x){
    priors <- rowSums(log(hap.adj.colprop[,x[x %in% colnames(hap.adj.colprop)]]+0.5/mut.nums));
    ## find the most likely haplogroup(s) given expected mutation frequencies
    best <- names(which(priors == max(priors)));
    ## when there's disagreement, pick the haplogroup with the most members
    names(which(hapGroupCounts[best] == max(hapGroupCounts[best])));
});
names(KCCG.haplotypes) <- sub(".Sig$","",names(KCCG.haplotypes));

KCCG.lnlike <- apply(test.vars,2,function(x){
    priors <- rowSums(log(hap.adj.colprop[,x[x %in% colnames(hap.adj.colprop)]]+0.5/mut.nums));
    ## find the most likely haplogroup(s) given expected mutation frequencies
    best <- names(which(priors == max(priors)));
    ## when there's disagreement, pick the haplogroup with the most members
    priors[names(which(hapGroupCounts[best] == max(hapGroupCounts[best])))];
});
names(KCCG.lnlike) <- sub(".Sig$","",names(KCCG.lnlike));


KCCG.haps.df <- read.csv("/mnt/ihbi_ngs/KCCG_deccles/KCCGid_UUID_mapping.csv");
KCCG.haps.df <- KCCG.haps.df[order(KCCG.haps.df$KCCGID),];
KCCG.haps.df$patID[is.na(KCCG.haps.df$patID)] <- 0;
KCCG.haps.df$matID[is.na(KCCG.haps.df$matID)] <- 0;
KCCG.haps.df$WGS.hap <- KCCG.haplotypes[KCCG.haps.df$KCCGID];
KCCG.haps.df$WGS.lnlike <- round(KCCG.lnlike[KCCG.haps.df$KCCGID],2);
write.csv(KCCG.haps.df, file="/data/all/david/mitochondria/KCCG_WGS_haplotypes.csv", row.names=FALSE);


inFile.IT <- "/data/all/david/mitochondria/STARout_bams/OVLNormalised_STARout_IonTorrentRuns_1-25_called.vcf.gz";
IT.vcf.df <- read.delim(inFile.IT, stringsAsFactors=FALSE,
                          skip=grep("^#CHROM",readLines(inFile.IT))-1);
colnames(IT.vcf.df)[1] <- "CHROM";
itNames <- grep("\\.bam$",colnames(IT.vcf.df), value=TRUE);
IT.vcf.df[,paste0(itNames,".GT")] <- sapply(IT.vcf.df[,itNames],sub,pattern=":.*$",replacement="");
IT.vcf.df$refSig <- paste(sprintf("%05d",IT.vcf.df$POS),IT.vcf.df$REF,IT.vcf.df$REF,sep=";");
IT.vcf.df$altSig <- paste(sprintf("%05d",IT.vcf.df$POS),IT.vcf.df$REF,IT.vcf.df$ALT,sep=";");
IT.vcf.df[,paste0(itNames,".Sig")] <- sapply(IT.vcf.df[,paste0(itNames,".GT")],
                                             function(x){ifelse(x == "1/1",IT.vcf.df$altSig, IT.vcf.df$refSig)});
IT.vcf.df[,paste0(itNames,".DPRr")] <- as.numeric(sapply(IT.vcf.df[,itNames],sub,pattern="^.*:(.*),.*$",replacement="\\1"));
IT.vcf.df[,paste0(itNames,".DPRa")] <- as.numeric(sapply(IT.vcf.df[,itNames],sub,pattern="^.*,",replacement=""));
IT.vcf.df[,"AAF"] <- rowSums(IT.vcf.df[,paste0(itNames,".DPRa")]) /
    (rowSums(IT.vcf.df[,paste0(itNames,".DPRa")]) + rowSums(IT.vcf.df[,paste0(itNames,".DPRr")]));
IT.vcf.df$cM <- 0;
IT.vcf.df[,paste0(itNames,".SigGT")] <- sapply(IT.vcf.df[,paste0(itNames,".GT")],
                                                 function(x){ifelse(x == "1/1"," C C",
                                                                    ifelse(x == "0/0"," A A",
                                                                           ifelse(x %in% c("0/1","1/0")," A C", " 0 0")))});
IT.vcf.df[["#CHROM"]] <- "26";

test.vars <- IT.vcf.df[,paste0(itNames,".Sig")];


IT.haplotypes <- apply(test.vars,2,function(x){
    priors <- rowSums(log(hap.adj.colprop[,x[x %in% colnames(hap.adj.colprop)]]+0.5/mut.nums));
    ## find the most likely haplogroup(s) given expected mutation frequencies
    best <- names(which(priors == max(priors)));
    ## when there's disagreement, pick the haplogroup with the most members
    names(which(hapGroupCounts[best] == max(hapGroupCounts[best])));
});
names(IT.haplotypes) <- sub(".Sig$","",names(IT.haplotypes));

IT.lnlike <- apply(test.vars,2,function(x){
    priors <- rowSums(log(hap.adj.colprop[,x[x %in% colnames(hap.adj.colprop)]]+0.5/mut.nums));
    ## find the most likely haplogroup(s) given expected mutation frequencies
    best <- names(which(priors == max(priors)));
    ## when there's disagreement, pick the haplogroup with the most members
    priors[names(which(hapGroupCounts[best] == max(hapGroupCounts[best])))];
});
names(IT.lnlike) <- sub(".Sig$","",names(IT.lnlike));

IT.runmap.df <- read.csv("/data/all/david/mitochondria/db/UUID_Run_Barcodes_Mitochondria.csv");
head(IT.runmap.df);
IT.runmap.df$sig <- paste0("run",IT.runmap.df$Sequencing.run,"_IonXpress_",IT.runmap.df$Barcode);

IT.haps.df <- data.frame(IT.runID = names(IT.lnlike), IT.hap = IT.haplotypes[names(IT.lnlike)], IT.lnlike = IT.lnlike);
rownames(IT.haps.df) <- NULL;
IT.haps.df$IT.runID <- sub(".bam$","",sub("IonXpress_0+","IonXpress_",sub("run0","run",IT.haps.df$IT.runID)));
IT.haps.df$IT.ID <- sub("\\.[12]_","_",IT.haps.df$IT.runID);
## Sanity check -- number of non-matching IDs should be low (~9)
IT.haps.df$IT.ID[which(!(IT.haps.df$IT.ID %in% IT.runmap.df$sig))];
## Add in UUID
IT.haps.df$UUID <- IT.runmap.df$UUID[match(IT.haps.df$IT.ID,IT.runmap.df$sig)];
uuid.df <- read.csv("/data/all/david/genabel/NI_UUID_Ped_2012-Oct-23.csv");
## Add in patID, matID
IT.haps.df$patID <- uuid.df$patID[match(IT.haps.df$UUID,uuid.df$UUID)];
IT.haps.df$matID <- uuid.df$matID[match(IT.haps.df$UUID,uuid.df$UUID)];
IT.haps.df$patID[is.na(IT.haps.df$patID)] <- 0;
IT.haps.df$matID[is.na(IT.haps.df$matID)] <- 0;
## reorder columns
IT.haps.df <- IT.haps.df[order(IT.haps.df$UUID),c("IT.runID","IT.ID","UUID","patID","matID","IT.hap","IT.lnlike")];

write.csv(IT.haps.df, file="/data/all/david/mitochondria/IT_run1-25_haplotypes.csv", row.names=FALSE);

str(IT.vcf.df);
sum(IT.vcf.df$altSig %in% names(mut.prop));
novel.vars.IT <- IT.vcf.df$altSig[!IT.vcf.df$altSig %in% names(mut.prop)];

cbind(unlist(IT.vcf.df[2,paste0(itNames,".Sig")]),
      unlist(IT.vcf.df[2,paste0(itNames,".GT")]),
      unlist(IT.vcf.df[3,paste0(itNames,".Sig")]),
      unlist(IT.vcf.df[3,paste0(itNames,".GT")]));

cat(novel.vars.KCCG,"\n");
cat(novel.vars.IT,"\n");

write.table(IT.vcf.df[,c("#CHROM","altSig","cM","POS",paste0(itNames,".SigGT"))], file="/data/all/david/mitochondria/IT_gt_matrix.tped",
            row.names=FALSE, col.names=TRUE, quote=FALSE);
write.table(cbind(itNames,1,0,0,0,0), file="/data/all/david/mitochondria/IT_gt_matrix.tfam", quote=FALSE,
            row.names=FALSE, col.names=FALSE);


write.table(IT.vcf.df[(!IT.vcf.df$altSig %in% names(mut.prop)) & (IT.vcf.df$AAF > 0.05),c("altSig","AAF")], row.names=FALSE, sep="\t",
            quote=FALSE);

##(IT.vcf.df[!IT.vcf.df$altSig %in% names(mut.prop),c("altSig","QUAL")]);
