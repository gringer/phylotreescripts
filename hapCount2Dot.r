#!/usr/bin/Rscript

setwd("/bioinf/QUT-2014-Mar-10-GBIS/mitochondria");

hapcounts.df <- read.csv("HapCounts_mtDNA_all_plus_NI.csv", stringsAsFactor=FALSE);
hapcounts.df <- subset(hapcounts.df, count >= 5);

## all haplogroups
{
  graphFile <- "phyloTree.dot";
  graphSVG <- "phyloTree.svg";
  cat(file = graphFile,
      paste(
        "digraph PhyloTree {\n",
        "  node [shape=circle, labelloc=c, fill=lightgray, ",
        "fontcolor=gray, fontname=Helvetica];\n",
        sep=""));
  for (i in 1:(dim(hapcounts.df)[1])){
    id <- i;
    name <- hapcounts.df$haplotype[i];
    parent <- hapcounts.df$parent[i];
    count <- hapcounts.df$count[i];
    pid <- which(hapcounts.df$haplotype == parent);
    if(length(pid) > 0){
        cat(file = graphFile, append = TRUE, sprintf("  %d -> %d;\n", pid, id));
    }
  }
  for (i in 1:(dim(hapcounts.df)[1])){
      hapName <- hapcounts.df$haplotype[i];
      hapCount <- hapcounts.df$count[i];
      hapSize <- sqrt(hapCount)/max(sqrt(hapcounts.df$count)) * 20;
      cat(file = graphFile, append = TRUE,
          sprintf("  %d [label=\"%s\\n%d\", fontsize=%0.2f, width=%0.2f];\n",
                  i, hapName, hapCount, hapSize*20, hapSize));
  }
  cat(file = graphFile, append = TRUE, "}\n");
  cat("Creating SVG image for",graphFile,"... ");
  commandLine <- sprintf("dot -Tsvg %s -o%s", graphFile, graphSVG);
  system(commandLine);
  cat("done", fill = TRUE);
}


getDescendants <- function(count.df, hapName){
    toFind <- hapName;
    resLines <- which(count.df$haplotype == hapName);
    while(length(toFind) > 0){
        nextHaps <- NULL;
        for(hapToFind in toFind){
            hapResLines <- which(count.df$parent == hapToFind);
            resLines <- c(resLines,hapResLines);
            nextHapNames <- count.df$haplotype[hapResLines];
            nextHaps <- c(nextHaps,nextHapNames);
            #cat(sprintf("adding as children of %s: %s\n",hapToFind,
            #            paste(nextHaps, collapse=";")));
        }
        toFind <- nextHaps;
    }
    res.df <- count.df[resLines,];
    res.df$parent[res.df$haplotype == hapName] <- NA;
    return(res.df);
}

phyloGraph <- function(tree.df, fileOutName, col = NULL){
  graphFile <- fileOutName;
  graphSVG <- sub("(\\.dot)?$",".svg",fileOutName);
  cat(file = graphFile,
      paste(
        "digraph PhyloGraph {\n",
        "  node [shape=circle, labelloc=c, fillcolor=\"#E0E0E0\", ",
        "fontcolor=black, fontname=Helvetica, style=filled];\n",
        "  edge [style=bold];\n",
        sep=""));
  for (i in 1:(dim(tree.df)[1])){
    id <- i;
    name <- tree.df$haplotype[i];
    parent <- tree.df$parent[i];
    count <- tree.df$count[i];
    pid <- which(tree.df$haplotype == parent);
    if(length(pid) > 0){
        cat(file = graphFile, append = TRUE, sprintf("  %d -> %d;\n", pid, id));
    }
  }
  for (i in 1:(dim(tree.df)[1])){
      hapName <- tree.df$haplotype[i];
      hapCount <- tree.df$count[i];
      hapSize <- sqrt(hapCount)/max(sqrt(tree.df$count)) * 20;
      if(length(col) == 0){
          cat(file = graphFile, append = TRUE,
              sprintf("  %d [label=\"%s\\n%d\", fontsize=%0.3f, width=%0.3f];\n",
                      i, hapName, hapCount, hapSize*5, hapSize/5));
      } else {
          cat(file = graphFile, append = TRUE,
              sprintf("  %d [label=\"%s\\n%d\", fontsize=%0.3f, width=%0.3f, fillcolor=\"%s\"];\n",
                      i, hapName, hapCount, hapSize*5, hapSize/5, col[i]));
      }
  }
  cat(file = graphFile, append = TRUE, "}\n");
  cat("Creating SVG image for",graphFile,"... ");
  commandLine <- sprintf("dot -Tsvg %s -o%s", graphFile, graphSVG);
  system(commandLine);
  cat("done", fill = TRUE);
}

## All
hapcounts.df <- read.csv("HapCounts_mtDNA_all_plus_NI.csv", stringsAsFactor=FALSE);
hapcounts.df <- subset(hapcounts.df, count >=10);
phyloGraph(hapcounts.df, "AllHaps_phylotree.dot",
           col = rainbow(26)[as.integer(factor(substring(hapcounts.df$haplotype,1,1)))]);
phyloGraph(hapcounts.df, "AllHaps_phylotree_bw.dot");


## B group only
hapcounts.df <- read.csv("HapCounts_mtDNA_all_plus_NI.csv", stringsAsFactor=FALSE);
hapcounts.df <- subset(hapcounts.df, count > 0);
hapB.df <- getDescendants(hapcounts.df,"B");
phyloGraph(hapB.df, "hapB_phylotree.dot");

## Norfolk Island only
hapcounts.df <- read.csv("HapCounts_mtDNA_NIonly.csv", stringsAsFactor=FALSE);
hapcounts.df <- subset(hapcounts.df, count > 0);
hapNI.df <- getDescendants(hapcounts.df,"L3");
phyloGraph(hapNI.df, "NI_phylotree.dot");
