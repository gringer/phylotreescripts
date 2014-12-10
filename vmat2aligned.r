#!/usr/bin/Rscript
library(Biostrings);

for(excludeNIMuts in c(TRUE, FALSE)){
    data.df <- read.delim("variant_matrix_all_plus_NI.tsv.gz", row.names = 1,
                          stringsAsFactors = FALSE);
    ids <- rownames(data.df);
    # exclude outlier-specific mutations
    if(excludeNIMuts){
        data.df <- data.df[,apply(data.df[-c(grep("^gi",ids),grep("^NI",ids)),],2,function(x){any(x != x[1])})];
    } else {
        data.df <- data.df[,apply(data.df[-grep("^gi",ids),],2,function(x){any(x != x[1])})];
    }
    maxIns <- paste(rep("-",max(apply(data.df,2,function(x){max(nchar(x))}))),collapse="");
    data.df <- sapply(data.df, function(x){
        ml <- max(nchar(x));
        if(ml > 1){
            x <- substring(paste(x,maxIns,sep=""),1,ml);
        }
        return(x);
    });
    rownames(data.df) <- ids;
    if(excludeNIMuts){
        cat(rbind(paste(">",ids,sep=""),apply(data.df,1,paste,collapse="")),
            sep="\n",
            file = "aligned_nonNImutsOnly_all_plus_NI.fa");
    } else {
        cat(rbind(paste(">",ids,sep=""),apply(data.df,1,paste,collapse="")),
            sep="\n",
            file = "aligned_allmuts_all_plus_NI.fa");
    }
}
