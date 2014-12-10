#!/usr/bin/Rscript

suppressMessages(library(gdata));
idMap.df <- read.xls("lists/sample_barcodes_2013-02-27.xls");
## add run to start of Sequencing run name
idMap.df$runTag <- paste("run",idMap.df$Sequencing.run,sep="_");
## use LAB_<ID> for individuals without UUID
idMap.df$UUID[is.na(idMap.df$UUID)] <- sprintf("LAB_%s",idMap.df$LabID[is.na(idMap.df$UUID)]);
idMap.df$UUID[idMap.df$UUID == 0] <- sprintf("LAB_%s",idMap.df$LabID[idMap.df$UUID == 0]);
###idMap.df$tag <- sprintf("IonXpress_%03d_R_%s_%s",idMap.df$Barcode,idMap.df$Run.Date,idMap.df$Run.Time);
## Generate NI string for individuals with UUID/LABID
idMap.df$idTag <- sprintf("NI_%s-run%s-ix%s",idMap.df$UUID,idMap.df$Sequencing.run,idMap.df$Barcode);
idMap.df$tag <- sprintf("%s_IonXpress_%03d",idMap.df$runTag,idMap.df$Barcode);

for(fileName in commandArgs(TRUE)){
  outName <- sub("\\.fast.*$","",fileName);
  if(outName %in% idMap.df$tag){
    outName <- idMap.df$idTag[which(outName == idMap.df$tag)];
  }
  cat(outName,fill = TRUE);
}
