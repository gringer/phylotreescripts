#!/usr/bin/Rscript
library(XML);
out <- htmlTreeParse("http://www.phylotree.org/tree/main.htm", isURL = TRUE, useInternalNodes = TRUE);

# Get location for "official" mtDNA tree HTML file
links <- data.frame(t(xpathSApply(out,"//a",function(x){c(text=xmlValue(x),xmlAttrs(x)["href"])})));
buildLoc <- as.character(links$href[grep("mtDNA tree Build", links$text)]);
if(!grepl("^http",buildLoc)){
  buildLoc <- paste("http://www.phylotree.org/tree/",buildLoc, sep="");
}

# Extract tree data and create R XML structure
tempZipDir <- tempfile();
dir.create(tempZipDir);
zipFileName <- sprintf("%s/%s",tempZipDir,sub("^.*/","",buildLoc));
download.file(buildLoc,zipFileName);
unzip(zipFileName, exdir = tempZipDir);
zipFiles <- list.files(tempZipDir);
treeFileName <- sprintf("%s/%s",tempZipDir,grep("^mtDNA.*htm$", zipFiles, value = TRUE));
mtDoc <- htmlTreeParse(treeFileName, useInternalNodes = TRUE,
                       ignoreBlanks = FALSE);
unlink(tempZipDir, recursive = TRUE);

mtTree <- xpathApply(mtDoc, "//tr", function(x){
  tBits <- sub("^\\s+$","",gsub("\\s+"," ",gsub(" ","",getChildrenStrings(x))));
  return(tBits[names(tBits) == "td"]);
});
if(grep("mt-MRCA",mtTree) > 1){
  mtTree <- mtTree[-(1:(grep("mt-MRCA",mtTree)-1))];
}

lastNodeName <- "";
tree.df <- NULL;
nameVec <- NULL;
maxLen <- max(sapply(mtTree,length));
unknownInc <- rep(0,maxLen);
parentNodes <- rep(NA,maxLen);
for(x in mtTree){
  if(length(x) < maxLen){
    x <- c(x,rep("",maxLen - length(x)));
  }
  if(!all(x == "")){
    namePos <- min(which(x != ""));
    parentName <- NA;
    if((x[namePos+1] != "") || (namePos == 1)){
      lastNodeName <- x[namePos];
      unknownInc[namePos] <- 0;
      nodeName <- x[namePos];
      mutations <- x[namePos+1];
      if(namePos > 1){
        parentName <- parentNodes[namePos-1];
      }
    } else {
      mutations <- x[namePos];
      namePos <- namePos - 1;
      unknownInc[namePos-1] <- unknownInc[namePos-1] + 1;
      parentName <- parentNodes[namePos-1];
      nodeName <- sprintf("%s_X%d", parentName, unknownInc[namePos-1]);
    }
    parentNodes[namePos] <- nodeName;
    nameVec <- c(nameVec, nodeName);
    tree.df <- rbind(tree.df, data.frame(
      name = nodeName, parent = parentName, level = namePos,
      mutations = mutations, row.names = nodeName, stringsAsFactors = FALSE));
  }
}

## add 'normalised' mutation definitions
##  - mutations are normalised to make them easier to parse downstream
##  - normalised format is similar to VCF file
##  - <location>;<ref>;<alt>
tree.df$normalised <-
  sapply(tree.df$mutations, function(x){
    ## split by " ", remove back-mutation indicators and parentheses
    working.split <- unlist(strsplit(gsub("(\\(|\\)|!)","",toupper(x))," "));
    working.split <- working.split[working.split != "RESERVED"];
    locs.pos <- regexpr("[0-9]+", working.split);
    locs <- as.numeric(substr(working.split,locs.pos,
                              locs.pos + attr(locs.pos,"match.length")-1));
    ref <- substr(working.split,1,locs.pos-1);
    alt <- substring(working.split,locs.pos + attr(locs.pos,"match.length"));
    ## reorder by base location (stable sort to keep mutation order)
    working.split <- working.split[order(locs)];
    ref <- ref[order(locs)];
    alt <- alt[order(locs)];
    locs <- locs[order(locs)];
    ## deal with single-base deletions
    ref[substr(alt,1,1) == "D"] <- "X";
    alt[substr(alt,1,1) == "D"] <- "";
    ## deal with single-base back deletions
    ref[grepl("^\\..*D$",alt)] <- "X";
    alt[grepl("^\\..*D$",alt)] <- "";
    ## deal with multi-base deletions
    mb.pos <- grep("D",alt);
    if(length(mb.pos) > 0){
      alt[mb.pos] <- gsub("(\\-|D)+","",alt[mb.pos]);
      mb.lengths <- as.numeric(alt[mb.pos]) - locs[mb.pos] + 1;
      alt[mb.pos] <- "";
      ref[mb.pos] <- sapply(mb.lengths,function(x){paste(rep("X",x), collapse = "")});
    }
    ## deal with insertions
    ins.strs <- alt[substr(alt,1,1) == "."];
    ins.strs <- substring(ins.strs,gregexpr("[^X0-9\\.]+",ins.strs));
    ref[substr(alt,1,1) == "."] <- "";
    ## convert to X because it's difficult to determine from VCF file
    alt[substr(alt,1,1) == "."] <- chartr("ACGT","XXXX",toupper(ins.strs));
    return(paste(locs,ref,alt,sep=";", collapse = " "));
  });

tree.df$name <- gsub("\\\"","\\'\\'",tree.df$name);
tree.df$parent <- gsub("\\\"","\\'\\'",tree.df$parent);

write.csv(tree.df, "mtDNA_haplogroups.csv", row.names = FALSE);

