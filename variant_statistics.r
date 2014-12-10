#!/usr/bin/Rscript

library(digest);
library(cluster);
library(gdata);

if(grepl("scripts$",getwd())){
    setwd("..");
}

dateStr <- format(Sys.Date(), "%Y-%b-%d");
data.df <- read.delim("results/NImt_variants_all.vcf.gz", skip = 35);
rownames(data.df) <- paste(data.df$POS,data.df$REF,sep=".");
# remove unknown IDs
data.df <- data.df[,!grepl("^(run|NI_LAB)", colnames(data.df))];
## TODO: LAB IDs are actually NIES IDs
# trim off information columns and non-genotype statistics
data.gtonly.mat <- sapply(data.df[,-(1:9)],sub,pattern=":.*$",replacement="");

data.mat <- data.gtonly.mat
# convert to categorical / numerical variables
data.mat[data.mat == "./."] <- NA;
gtLevels <- levels(factor(data.mat));
gtVals <- 1:length(gtLevels);
names(gtVals) <- gtLevels;
data.mat[1:length(data.mat)] <- gtVals[data.mat];
data.mat <- matrix(as.numeric(data.mat),nrow = nrow(data.df));
rownames(data.mat) <- rownames(data.df);
colnames(data.mat) <- colnames(data.df)[-(1:9)];
data.mat[is.na(data.mat)] <- -1;

table(data.df$QUAL == 999);

table(data.df$QUAL == 999, (nchar(as.character(data.df$REF)) == 1) & (nchar(as.character(data.df$ALT)) == 1));

## repeated individuals have a repeat run
repeated.indivs <- unique(sub("\\..*$","",grep("run[0-9]+\\.1\\.ix",colnames(data.mat), value = T)));
## repeated individuals also have a non-repeated run (e.g. no data for run 13)
single.indivs <- unique(sub("\\..*$","",grep("run[0-9]+\\.ix",colnames(data.mat), value = T)));
repeated.indivs <- intersect(repeated.indivs,single.indivs);
repeated.cols <- colnames(data.mat)[sub("\\..*$","",colnames(data.mat)) %in% repeated.indivs];
## to exclude a bit of bias, only choose individuals who have two runs exactly
repeated.counts <- sort(table(sub("\\..*$","",repeated.cols))) == 2;
repeated.indivs <- sort(names(repeated.counts)[repeated.counts]);
repeated.cols <- colnames(data.mat)[sub("\\..*$","",colnames(data.mat)) %in% repeated.indivs];

## now generate statistics for variants
repeated.mat <- data.mat[,repeated.cols];
consistent.mat <-
    repeated.mat[,seq(1,length(repeated.cols),2)] ==
    repeated.mat[,seq(2,length(repeated.cols),2)];
colnames(consistent.mat) <- sub("\\..*$","", colnames(consistent.mat));
prop.byVar <- apply(consistent.mat,1,function(x){sum(x) / length(x)});
## work out individual consistency by first excluding variants with low quality
plot(density(prop.byVar[prop.byVar > 0.8 & prop.byVar < 1]));
## density plot suggest a critical point near 0.90, so use that as threshold
consistent.mat <- consistent.mat[data.df$QUAL == 999 & prop.byVar > 0.90,];
prop.byInd <- t(apply(consistent.mat,2,function(x){sum(x) / length(x)}));

## show proportions for repeats
table(round(prop.byVar,3));
table(prop.byVar == 1);
table(prop.byVar >= 0.9);

table(round(prop.byInd,3));

### Consistency statistics end here ###

## merge genotypes for individuals
ind.ids <- sub("\\..*$","",colnames(data.mat));
merged.df <- data.frame(row.names = rownames(data.mat));
for(id in unique(ind.ids)){
    cols = which(ind.ids == id);
    merged.df[[id]] <-
        apply(matrix(data.mat[,cols],ncol = length(cols)),1,function(x){
        if(all((x == x[1]))){
            x[1];
        } else if(all(((x != 1) == (x[1] != 1)))){
            ## all non-reference, but dissimilar
            -2;
        } else {
            ## all different
            -3;
        }
    });
}
## calculate minor allele frequency (only for consistent genotypes,
## lumping all non-reference variants as a non-reference allele)
variants.maf <-
    apply(merged.df,1,function(x){min(table(x[x>=1] == 1) / sum(x[x>=1]))});

marker.stats.df <- data.df[,2:6];
marker.stats.df$ID <- NULL;
marker.stats.df[names(prop.byVar),"IAF"] <- round(1-prop.byVar,3);
marker.stats.df[names(variants.maf),"MAF"] <- round(variants.maf,3);
marker.stats.df$INDEL <- grepl("INDEL",data.df$INFO);
marker.stats.df$ALT <- chartr(",",";",marker.stats.df$ALT);
infoData <- chartr(",;",";,",sub("INDEL;","",data.df$INFO));
names(infoData) <- rownames(data.df);
indel.data <- data.frame(t(sapply(strsplit(infoData[marker.stats.df$INDEL],","),
                                  sub,pattern="^.*?=",replacement="")));
colnames(indel.data) <-
    sub("=.*$","",unlist(strsplit((infoData[marker.stats.df$INDEL])[1],",")));
indel.data$RPB <- NA;
snp.data <- data.frame(t(sapply(strsplit(infoData[!marker.stats.df$INDEL],","),
                                sub,pattern="^.*?=",replacement="")));
colnames(snp.data) <-
    sub("=.*$","",unlist(strsplit((infoData[!marker.stats.df$INDEL])[1],",")));
snp.data$IS <- NA;
both.data <- rbind(snp.data, indel.data);
for(x in colnames(both.data)){
    marker.stats.df[rownames(both.data),x] <- both.data[,x];
}
marker.stats.df$RPB <- sub("e\\+","e",marker.stats.df$RPB);
marker.stats.df$HWE <- sub("e\\+","e",marker.stats.df$HWE);
write.csv(marker.stats.df,sprintf("results/Marker_stats_NImt_variants_all_%s.csv",dateStr),
          row.names = FALSE, quote = FALSE);

pdf(sprintf("results/NImt_variants_heatmap_dendrogram_all_%s.pdf",dateStr),
    width = 40, height = 8);
tCols.all = c("blue",rep("green",2),"yellow",rep("cyan",5),"white");
## display heatmap
tmp.mat <- data.mat;
tmp.mat[tmp.mat < 0] <- NA;
cluster.row <- diana(t(tmp.mat), metric = "manhattan");
tmp.mat[tmp.mat < 0] <- length(gtVals) + 1;
layout(matrix(1:2), heights = c(0.90,0.10));
par(mar = c(5,5,3,1));
image(x = as.numeric(sub("\\..*$","",rownames(data.mat)))+
      (1:dim(data.mat)[1])/(dim(data.mat)[1]),
      y = 1:dim(data.mat)[2], xlab = "Chromosome Location (bp)",
      z = tmp.mat[,cluster.row$order],
      ylab = "Individual",
      main = sprintf("Variant Plot (%s samples)",dim(data.mat)[2]),
      col = tCols.all);
par(mar = c(0,0,0,0));
plot.new();
legend("center",legend = c(gtLevels,"NA"), horiz = TRUE,
       fill = tCols.all);
## colour labels based on hash of indvidual ID
aD <- dendrapply(as.dendrogram(cluster.row), function(n) {
    if(is.leaf(n)){
        a <- attributes(n)
        labelCol <- paste("#",substring(digest(sub("\\..*$","",
                       attr(n,"label"))),1,6), sep="");
        attr(n, "edgePar") <- c(a$nodePar, list(col = labelCol));
        attr(n, "nodePar") <- list(pch = NA, lab.cex = 0.5, lab.col = labelCol);
    }
    n;
    });
layout(matrix(1));
par(mar = c(5,5,3,1));
# show dendrogram
plot(aD, main = "Dendrogram plot of Individual Samples");
ind.ids <- sub("\\..*$","",colnames(data.mat));
graphics.off();

pdf(sprintf("results/NImt_variants_heatmap_dendrogram_good_%s.pdf",dateStr),
    width = 40, height = 8);
markers.good.df <- subset(marker.stats.df,
                          (IAF <= 0.1) & (QUAL == 999));
tCols.all = c("blue",rep("green",2),"yellow",rep("cyan",5),"white");
names(tCols.all) <- 1:10;
## display heatmap
tmp.mat <- data.mat[rownames(markers.good.df),];
tmp.mat[tmp.mat < 0] <- NA;
cluster.row <- diana(t(tmp.mat), metric = "manhattan");
tmp.mat[is.na(tmp.mat)] <- length(gtVals) + 1;
#tCols.all <- tCols.all[sort(names(table(tmp.mat)))];
layout(matrix(1:2), heights = c(0.90,0.10));
par(mar = c(5,5,3,1));
image(x = as.numeric(sub("\\..*$","",rownames(tmp.mat)))+
      (1:dim(tmp.mat)[1])/(dim(tmp.mat)[1]),
      y = 1:dim(tmp.mat)[2], xlab = "Chromosome Location (bp)",
      z = tmp.mat[,cluster.row$order],
      ylab = "Individual",
      main = sprintf("Variant Plot (%s samples, 'good' variants only)",dim(tmp.mat)[2]),
      col = tCols.all);
par(mar = c(0,0,0,0));
plot.new();
legend("center",legend = c(gtLevels,"NA"), horiz = TRUE,
       fill = tCols.all);
## colour labels based on hash of indvidual ID
aD <- dendrapply(as.dendrogram(cluster.row), function(n) {
    if(is.leaf(n)){
        a <- attributes(n)
        labelCol <- paste("#",substring(digest(sub("\\..*$","",
                       attr(n,"label"))),1,6), sep="");
        attr(n, "edgePar") <- c(a$nodePar, list(col = labelCol));
        attr(n, "nodePar") <- list(pch = NA, lab.cex = 0.5, lab.col = labelCol);
    }
    n;
    });
layout(matrix(1));
par(mar = c(5,5,3,1));
# show dendrogram
plot(aD, main = "Dendrogram plot of Individual Samples (using only 'good' variants for clustering)");
ind.ids <- sub("\\..*$","",colnames(data.mat));
graphics.off();


pdf(sprintf("results/NImt_variants_heatmap_dendrogram_merged_all_%s.pdf",dateStr),
    width = 40, height = 8);
tCols.merge = c("blue",rep("green",2),"yellow",rep("cyan",5),"white", "red");
merged.mat <- as.matrix(merged.df);
merged.na.mat <- merged.mat;
merged.na.mat[merged.na.mat < 0] <- NA;
m.cluster.row <- diana(t(merged.na.mat), metric = "manhattan");
tmp.mat <- merged.mat;
tmp.mat[tmp.mat <= -2] <- length(gtVals) + 2;
tmp.mat[tmp.mat <= -1] <- length(gtVals) + 1;
layout(matrix(1:2), heights = c(0.90,0.10));
par(mar = c(5,5,3,1));
image(as.numeric(sub("\\..*$","",rownames(merged.mat)))+
      (1:dim(merged.mat)[1])/(dim(merged.mat)[1] * 2),
      1:dim(merged.mat)[2], xlab = "Chromosome Location (bp)",
      ylab = "Individual",
      main = sprintf("Merged Variant Plot (%s individuals)",
          dim(merged.mat)[2]),
      tmp.mat[,m.cluster.row$order], col = tCols.merge);
par(mar = c(0,0,0,0));
plot.new();
legend("center",legend = c(gtLevels,"NA", "INC"), horiz = TRUE,
       fill = c(tCols.merge, "red"));
layout(matrix(1));
par(mar = c(5,5,3,1));
## colour labels based on hash of indvidual ID
aD <- dendrapply(as.dendrogram(m.cluster.row), function(n) {
    if(is.leaf(n)){
        a <- attributes(n)
        labelCol <- paste("#",substring(digest(sub("\\..*$","",
                       attr(n,"label"))),1,6), sep="");
        attr(n, "edgePar") <- c(a$nodePar, list(col = labelCol));
        attr(n, "nodePar") <- list(pch = NA, lab.cex = 0.5, lab.col = labelCol);
    }
    n;
    });
# show dendrogram
plot(aD, main = "Dendrogram plot of Merged Individuals");
graphics.off();

## require QUAL=999, and <=10% inconsistent alleles for a marker to be considered good
markers.good.df <- subset(marker.stats.df,
                          (IAF <= 0.1) & (QUAL == 999));
merged.good.mat <- merged.mat[rownames(markers.good.df),];
apply(merged.good.mat,1,function(x){sum(x<=-2)})
merged.na.good.mat <- merged.good.mat;
merged.na.good.mat[merged.na.good.mat < 0] <- NA;
mg.cluster.row <- diana(t(merged.na.good.mat), metric = "manhattan");
tmp.mat <- merged.good.mat;
tmp.mat[tmp.mat <= -2] <- length(gtVals) + 2;
tmp.mat[tmp.mat <= -1] <- length(gtVals) + 1;
pdf(sprintf("results/NImt_variants_heatmap_dendrogram_goodmerged_all_%s.pdf", dateStr),
            width = 40, height = 8);
layout(matrix(1:2), heights = c(0.90,0.10));
par(mar = c(5,5,3,1));
image(as.numeric(sub("\\..*$","",rownames(merged.good.mat)))+
      (1:dim(merged.good.mat)[1])/(dim(merged.good.mat)[1] * 2),
      1:dim(merged.good.mat)[2], xlab = "Chromosome Location (bp)",
      ylab = "Individual",
      main = sprintf("Merged Variant Plot (%s individuals, 'good' variants only)",
          dim(merged.good.mat)[2]),
      tmp.mat[,m.cluster.row$order], col = tCols.merge);
par(mar = c(0,0,0,0));
plot.new();
legend("center",legend = c(gtLevels,"NA", "INC"), horiz = TRUE,
       fill = c(tCols.merge, "red"));
layout(matrix(1));
par(mar = c(5,5,3,1));
## colour labels based on hash of indvidual ID
aD <- dendrapply(as.dendrogram(mg.cluster.row), function(n) {
    if(is.leaf(n)){
        a <- attributes(n)
        labelCol <- paste("#",substring(digest(sub("\\..*$","",
                       attr(n,"label"))),1,6), sep="");
        attr(n, "edgePar") <- c(a$nodePar, list(col = labelCol));
        attr(n, "nodePar") <- list(pch = NA, lab.cex = 0.5, lab.col = labelCol);
    }
    n;
    });
# show dendrogram
plot(aD, main = "Dendrogram plot of Merged Individuals\n(using only 'good' variants for clustering)");
graphics.off();

plinkVals <- gtVals;
names(plinkVals) <- sub("/"," ",chartr("0123","1222",names(plinkVals)));
merged.na.plink.mat <- merged.na.good.mat;
colnames(merged.na.plink.mat) <- sub("^NI_","",colnames(merged.na.plink.mat))
merged.na.plink.mat[!is.na(merged.na.plink.mat)] <-
    names(plinkVals)[merged.na.plink.mat[!is.na(merged.na.plink.mat)]];
merged.na.plink.mat[is.na(merged.na.plink.mat)] <- "0 0";
## deal with duplicate positions (e.g. SNP and INDEL at same location)
## by adding a small amount to each position
poss <- as.numeric(sub("\\..*$","",(rownames(merged.na.plink.mat))));
poss <- poss + 1:length(poss) / 1000;
## convert to plink-like data frame
plink.df <- data.frame(row.names = rownames(merged.na.plink.mat),
                       chr = 25, name = rownames(merged.na.plink.mat),
                       cmPos = 0, pos = poss);
plink.df[,colnames(merged.na.plink.mat)] <-
    merged.na.plink.mat;
ped.df <- read.csv("data/NI_UUID_Ped_2012-Oct-23.csv", row.names = 1, stringsAsFactors = FALSE);
ped.df <- ped.df[colnames(merged.na.plink.mat),c("patID","matID","Gender")];
ped.df$Gender[ped.df$Gender == "Male"] <- "1";
ped.df$Gender[ped.df$Gender == "Female"] <- "2";
ped.df[is.na(ped.df)] <- "0";
plink.tfam.df <- data.frame(famID = 1, indID = colnames(merged.na.plink.mat),
                            patID = ped.df$patID, matID = ped.df$matID,
                            sex = ped.df$Gender, AS = -9);
## write out to file
write.fwf(plink.df, colnames = FALSE, sep = "  ",
          sprintf("results/NImt_goodVar_seqMerged_%s.tped", dateStr));
write.fwf(plink.tfam.df, colnames = FALSE,
          sprintf("results/NImt_goodVar_seqMerged_%s.tfam", dateStr));

## create BED format file
markers.good.df$chrom <- "RSRS";
markers.good.df$start <- markers.good.df$POS - 1;
markers.good.df$altLen <- sapply(strsplit(markers.good.df$ALT,";"),
    function(x){max(nchar(x))});
markers.good.df$maxLen <- pmax(nchar(as.character(markers.good.df$REF)),
                               markers.good.df$altLen);
markers.good.df$end <- markers.good.df$start + markers.good.df$maxLen;
markers.good.df$name <- rownames(markers.good.df);

write.fwf(markers.good.df[,c("chrom","start","end","name")],
          file = "results/NImt_goodVar.bed");
