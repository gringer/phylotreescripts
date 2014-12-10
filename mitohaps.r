#!/usr/bin/Rscript

if(grepl("scripts",getwd())){
    setwd("..");
}

data.df <- read.csv("results/NIhaplotypes_good.csv.gz", comment.char = "[",
                    stringsAsFactors = FALSE);
data.nmrca.df <- subset(data.df,(Best.Haplotype != "mt-MRCA") & (Match.Proportion > 0));

maxScore <- tapply(data.nmrca.df$Match.Proportion,data.nmrca.df$Name,max);

data.best.df <- subset(data.df, Match.Proportion == maxScore[Name]);
maxLevel <- tapply(data.best.df$Level, data.best.df$Name, max);

data.bestLevel.df <- subset(data.df, Level == maxLevel[Name]);
write.csv(data.bestLevel.df,"results/NIhaplotypes_highestLevelScore.csv",
          row.names = FALSE);
