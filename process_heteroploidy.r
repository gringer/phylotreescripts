#!/usr/bin/Rscript

library(digest);
library(cluster);
library(gdata);

dateStr <- format(Sys.Date(), "%Y-%b-%d");

het.df <- read.csv("results/heteroploidy_all.csv", stringsAsFactors = FALSE);
## remove unlabelled sequences
het.df <- het.df[-grep("IonXpress_0(61|74|89)", het.df$ID),]
het.df$Alt[het.df$Alt == "."] <- het.df$Ref[het.df$Alt == "."];

het.unique.df <- unique(het.df[,c("Pos","Ref","Alt")]);
het.unique.df <- het.unique.df[order(het.unique.df$Pos,
                                     het.unique.df$Ref,het.unique.df$Alt),];

## load barcode -> UUID mappings
mapping.labid.df <- read.xls("lists/sample_barcodes_2014-04-14.xls");
mapping.labid.df$UUID[mapping.labid.df$UUID == 0] <-
    paste("LAB",mapping.labid.df$LabID[mapping.labid.df$UUID == 0],sep="_");
rownames(mapping.labid.df) <- paste("run",mapping.labid.df$Sequencing.run,
                                    "IonXpress",
                                    sprintf("%03d",mapping.labid.df$Barcode),
                                    sep = "_");

mapping.labid.df$code <-
    paste(sep="","NI_",mapping.labid.df$UUID,"-run",
          mapping.labid.df$Sequencing.run,"-ix",
          mapping.labid.df$Barcode);

## remove barcodes that cannot be determined
het.df <- het.df[-grep("run_11_IonXpress_006",het.df$ID),];

## convert unconverted runs to new ID
het.df$newID <- het.df$ID;
het.df$newID[grep("IonXpress",het.df$ID)] <-
    mapping.labid.df[grep("IonXpress",het.df$ID,value = TRUE),"code"];


## Calculate total number of high-likelihood reads
het.df$HQT <- rowSums(het.df[,grep("HQ",colnames(het.df))]);
## exclude low-scoring variant positions
het.df <- subset(het.df, HQT > 25); # based on mutMins density / hist plot
## Calculate total number of high-likelihood alternate reads
het.df$HQAT <- rowSums(het.df[,grep("HQA",colnames(het.df))]);
het.df$mutSig <- paste(het.df$Pos,het.df$Ref,sep=":");
het.df$mutStart <- het.df$Pos;
het.df$mutEnd <- het.df$Pos + nchar(as.character(het.df$Ref)) - 1;
het.df$mutRange <- paste(het.df$mutStart,het.df$mutEnd, sep="-");
mutSums <- tapply(het.df$HQT,het.df$mutSig,sum);
mutMins <- tapply(het.df$HQT,het.df$mutSig,min); # best distribution
mutMaxes <- tapply(het.df$HQT,het.df$mutSig,max);
## flag reference variants
het.df$isVariant <- !(het.df$Alt == het.df$Ref);
## SNPs have a reference sequence length of 1, and 1 base per additional ';'
het.df$isSNP <- (nchar(het.df$Ref) == 1) &
    ((nchar(het.df$Alt)+1) / 2 == nchar(gsub(";","",het.df$Alt)));
het.df$Alt[het.df$isSNP & !het.df$isVariant] <- ".";
## calculate alternate allele frequency
het.df$altFreq <-
  rowSums(het.df[,grep("HQA",colnames(het.df))]) / het.df$HQT;

## get ranges for "good" variants
stats.df <- read.csv("results/Marker_stats_NImt_variants_all_2014-May-08.csv",
                     stringsAsFactors = FALSE);
stats.df$isSNP <- (nchar(stats.df$REF) == 1) & (nchar(stats.df$ALT) == 1);
stats.df$REFALT <- paste(stats.df$REF,stats.df$ALT,sep=";");
## set up signature function
getSig <- function(x){
    if(length(x) > 1){
        return(sapply(x,getSig));
    } else {
        vals <- unlist(strsplit(x,";"));
        refLen <- nchar(vals[1]);
        altLens <- nchar(vals[-1]);
        diffs <- refLen - altLens;
        inss <- (diffs < 0);
        dels <- (diffs > 0);
        retStr <- "";
        if(sum(diffs == 0) > 0){
            retStr <- paste((vals[-1])[diffs == 0], sep = ";");
        }
        if(sum(inss) >= 1){
            if(nchar(retStr) > 0){
                retStr = sprintf("%s;",retStr);
            }
            if(sum(inss) == 1){
                retStr <- sprintf("%sI%d", retStr, -diffs[inss]);
            }
            if(sum(inss) > 1){
                retStr <-
                    sprintf("%sI%s", retStr,
                            paste(range(-diffs[inss]),sep="-"));
            }
        }
        if(sum(dels) >= 1){
            if(nchar(retStr) > 0){
                retStr = sprintf("%s;",retStr);
            }
            if(sum(dels) == 1){
                retStr <- sprintf("%sD%d", retStr, diffs[dels]);
            }
            if(sum(dels) > 1){
                retStr <-
                    sprintf("%sD%s", retStr, paste(range(diffs[dels]),sep="-"));
            }
        }
        return(retStr);
    }
}
# get signature for variants
stats.df$signature <-
    paste(stats.df$POS,"{",getSig(stats.df$REFALT),"}", sep = "");
stats.df$signature[stats.df$INDEL] <-
    sprintf("%d-%d{%s}", stats.df$POS[stats.df$INDEL],
            stats.df$POS[stats.df$INDEL] +
            sapply(strsplit(stats.df$REFALT[stats.df$INDEL],";"),
                   function(x){max(nchar(x))-1}),
            getSig(stats.df$REFALT[stats.df$INDEL]));


## find contiguous regions for INDELs
hetSub.df <- subset(het.df, !isSNP & isVariant);
## combine ranges from good variants and all variants
mutPoss <-
    unique(c(het.df[!het.df$isSNP,"mutRange"],
             sub("\\{.*$","",
                 stats.df$signature[grepl("-",stats.df$signature)])));
mutPoss.df <- data.frame(sig = mutPoss);
mutPoss.df$s <- as.numeric(sub("-.*$","",mutPoss.df$sig));
mutPoss.df$e <- as.numeric(sub("^.*?-","",mutPoss.df$sig));
## order by position
mutPoss.df <- mutPoss.df[order(mutPoss.df$s,mutPoss.df$e),];

## merge adjacent ranges
cR <- 1; # currentRange
mutRanges <- mutPoss.df[1,-1];
for(i in 2:dim(mutPoss.df)[1]){
    ## single base for reference cannot modify pre-existing range
    if(mutPoss.df$s[i] != mutPoss.df$e[i]){
        if(mutPoss.df$s[i] <= mutRanges$e[cR]){
            mutRanges$e[cR] <- max(mutRanges$e[cR],mutPoss.df$e[i]);
        } else {
            cR <- cR + 1;
            mutRanges[cR,] <- mutPoss.df[i,-1];
        }
    }
}
mutRanges$sig <- paste(mutRanges$s,mutRanges$e,sep="-");

## label INDELs with the corrected extended range
het.df$mutRange <- as.character(cut(het.df$mutStart, right = FALSE,
                                     c(mutRanges$s,max(het.df$Pos+1)),
                                     labels = mutRanges$sig));
het.df$mutRange[het.df$isSNP] <- het.df$mutStart[het.df$isSNP];
stats.df$mutRange <- as.character(cut(stats.df$POS, right = FALSE,
                                      c(mutRanges$s, max(stats.df$POS+1)),
                                      labels = mutRanges$sig));
stats.df$mutRange[stats.df$isSNP] <- stats.df$POS[stats.df$isSNP];

## determine relative length change
het.df$refLen <- nchar(het.df$Ref);
het.df$mutDesc <- paste(het.df$Ref,het.df$Alt,sep="->");
het.df$altLen <- sapply(strsplit(het.df$Alt,";"),
                        function(x){paste(nchar(x), collapse = ",")});
het.df$diffStr <- apply(het.df[,c("refLen","altLen")], 1, function(x){
    rl <- as.numeric(x["refLen"]);
    al <- as.numeric(strsplit(x["altLen"],",")[[1]]);
    paste(sort(unique(al - rl)),collapse=";");
});
## for no length change, fill in alternate SNPs
het.df$diffStr[het.df$isSNP] <- het.df$Alt[het.df$isSNP];


## combine variants with the same reference position
muts.all <-
    tapply(het.df$diffStr,het.df$mutRange,function(x){
        variants <- sort(unique(unlist(strsplit(x,";"))));
        alts <- paste(variants[variants != "."],collapse=";");
    });
## add in SNPs from "good" variants (where no variants are observed)
acGoodSNPs <- as.character(stats.df[stats.df$isSNP,"POS"]);
muts.all[acGoodSNPs[!(acGoodSNPs %in% names(muts.all))]] <-
    "";
muts.all[acGoodSNPs] <-
    sub("^;","",
        paste(muts.all[acGoodSNPs],stats.df[stats.df$isSNP,"ALT"], sep = ";"));
muts.all[acGoodSNPs] <-
    sapply(strsplit(muts.all[acGoodSNPs],";"),
           function(x){paste(unique(x),collapse=";")});
## convert numeric length differences into I<min>-<max>;D<min>-<max> notation
muts.all <- sapply(muts.all, function(x){
    if(!grepl("[A-Z]",x)){
        comps <- as.numeric(strsplit(x,";")[[1]]);
        outVec <- NULL;
        if(any(comps>0)){
            outVec <- sprintf("I%s",paste(unique(range(comps[comps>0])),
                                       collapse="-"));
        }
        del <- "";
        if(any(comps<0)){
            outVec <- c(outVec,
                        sprintf("D%s",paste(unique(range(-comps[comps<0])),
                                            collapse="-")));
        }
        return(paste(outVec,collapse=";"));
    }
    return(x);
});
muts.all[1:length(muts.all)] <- sprintf("%s{%s}",names(muts.all),muts.all);
## remove non-variant changes from list of mutations
muts.all <- muts.all[!grepl("\\{\\}",muts.all)];
muts.all <- muts.all[order(as.numeric(sub("-.*$","",names(muts.all))))];
het.df$mutStr <- muts.all[as.character(het.df$mutRange)];
# update "good" statistics file with correct mutation signature
stats.df$mutSig <- muts.all[stats.df$mutRange];
stats.df$mutSig[stats.df$QUAL != 999] <- NA;
stats.df$mutSig[stats.df$IAF > 0.1] <- NA;

mut.range.df <- data.frame(row.names = muts.all);
# determine start/end points for tagged mutations
mut.range.df$start <-
    sapply(strsplit(rownames(mut.range.df), split = "(-|\\{)"),
           function(x){x[1]});
mut.range.df$end <-
    sub("^[0-9]+-","",sapply(strsplit(rownames(mut.range.df), split = "\\{"),
           function(x){x[1]}));
mut.range.df$indel <- grepl("(I|D)",rownames(mut.range.df));

## Determine Alternate allele frequencies
het.HQT.mat <- unclass(xtabs(HQT ~ newID + mutStr,data = het.df));
het.HQAT.mat <- unclass(xtabs(HQAT ~ newID + mutStr,data = het.df));
het.HQAF.mat <- het.HQAT.mat / het.HQT.mat;
het.HQAF.mat[is.nan(het.HQAF.mat)] <- 0;

# generate count matrix
het.HQAC.mat <- het.HQT.mat;
het.HQAC.mat[1:length(het.HQAC.mat)] <-
    paste(het.HQAT.mat, "/", het.HQT.mat,sep="");
het.HQAC.mat[het.HQAC.mat == "0/0"] <- "";

# filter out variants that aren't in the "good" set
keep <- intersect(stats.df$mutSig, colnames(het.HQAF.mat));
het.HQAF.mat <- het.HQAF.mat[,keep];
het.HQAC.mat <- het.HQAC.mat[,keep];

## order columns based on base position
var.order <- order(as.numeric(sub("(-|\\{).*$","",colnames(het.HQAF.mat))));
het.HQAF.mat <- het.HQAF.mat[,var.order];
het.HQAC.mat <- het.HQAC.mat[,var.order];
## write to file
write.csv(round(het.HQAF.mat,3),
          (sprintf("results/AltFrac_positions_good_seqMerged_%s.csv",
                         dateStr)));
write.csv(het.HQAC.mat,
          (sprintf("results/AltCount_positions_good_seqMerged_%s.csv",
                         dateStr)));

hq.indels <- het.HQAF.mat[,grep("(I|D)",colnames(het.HQAF.mat))];
hq.snps <- het.HQAF.mat[,-grep("(I|D)",colnames(het.HQAF.mat))];


pdf(sprintf("results/AltFrac_hist_%s.pdf",
                         dateStr));
hist((hq.indels[(hq.indels >= 0) & (hq.indels <= 1)]),
             kernel = "rectangular", main =
     "Variant proportion for INDELs in NI population",
     breaks = 100, col = "red");
hist((hq.snps[(hq.snps > 0.2) & (hq.snps < 0.8)]),
             kernel = "rectangular", main =
     "Variant proportion for heteroplasmic SNPs in NI population",
     breaks = 100, col = "red");
graphics.off();

het.HQAF.trin <- het.HQAF.mat;
het.HQAF.trin[het.HQAF.trin < 0.2] <- 0;
het.HQAF.trin[het.HQAF.trin > 0.8] <- 2;
het.HQAF.trin[(het.HQAF.trin > 0) & (het.HQAF.trin < 1)] <- 1;
write.csv(round(het.HQAF.trin,3),
          sprintf("results/AltTrin_positions_good_seqMerged_%s.csv",
                  dateStr));

data.mat <- t(het.HQAF.trin);
# remove unknown IDs
data.mat <- data.mat[,!grepl("^(run|NI_LAB)", colnames(data.mat))];
## TODO: LAB IDs are actually NIES IDs
repeated.indivs <- unique(sub("\\..*$","",grep("run[0-9]+\\.1-ix",colnames(data.mat), value = T)));


repeated.indivs <- unique(sub("-.*$","",grep("run[0-9]+\\.1-ix",colnames(data.mat), value = T)));
## repeated individuals also have a non-repeated run (e.g. no data for run 13)
single.indivs <- unique(sub("-.*$","",grep("run[0-9]+-ix",colnames(data.mat), value = T)));
repeated.indivs <- intersect(repeated.indivs,single.indivs);
repeated.cols <- colnames(data.mat)[sub("-.*$","",colnames(data.mat)) %in% repeated.indivs];
## to exclude a bit of bias, only choose individuals who have two runs exactly
repeated.counts <- sort(table(sub("-.*$","",repeated.cols))) == 2;
repeated.indivs <- sort(names(repeated.counts)[repeated.counts]);
repeated.cols <- colnames(data.mat)[sub("-.*$","",colnames(data.mat)) %in% repeated.indivs];

## now generate statistics for variants
repeated.mat <- data.mat[,repeated.cols];
consistent.mat <-
    repeated.mat[,seq(1,length(repeated.cols),2)] ==
    repeated.mat[,seq(2,length(repeated.cols),2)];
colnames(consistent.mat) <- sub("-.*$","", colnames(consistent.mat));
prop.byVar <- apply(consistent.mat,1,function(x){sum(x) / length(x)});
## work out individual consistency by first excluding variants with low quality
plot(density(prop.byVar[prop.byVar > 0.8 & prop.byVar < 1]));
## okay, that's nice, but the filtering has already been done...
## [don't exclude SNPs with low consistency]
prop.byInd <- t(apply(consistent.mat,2,function(x){sum(x) / length(x)}));

## show proportions for repeats
table(round(prop.byVar,3));
table(prop.byVar == 1);
table(prop.byVar >= 0.9);

table(round(prop.byInd,3));

## merge genotypes for individuals
ind.ids <- sub("-.*$","",colnames(data.mat));
merged.df <- data.frame(row.names = rownames(data.mat));
for(id in unique(ind.ids)){
    cols = which(ind.ids == id);
    merged.df[[id]] <-
        apply(matrix(data.mat[,cols],ncol = length(cols)),1,function(x){
        if(all((x == x[1]))){
            x[1];
        } else if(all(((x != 0) == (x[1] != 0)))){
            ## all non-reference, but dissimilar
            -1;
        } else {
            ## all different
            -2;
        }
    });
}
