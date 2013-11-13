#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------
# Link Amino acids to SIFT predictions
# Script written by Sofiane mezmouk (Ross-Ibarra laboratory)
#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------

for (chr in 1:10){
    # AA is the output of SNPtoAA.r
  AA <- read.table(paste("/group/jrigrp2/DiallelSofiane/AA/NucToAA-12_Chr",chr,".txt",sep=""),sep="\t", header=T)
  tmp <- read.table(paste("/home/smezmouk/DeleteriousMutations/SIFT/ResultSummary/resSIFT-chr",chr,sep=""),sep="\t", header=T)
  
  bof1 <- paste(as.vector(AA[,2]),as.vector(AA[,3]),sep="x")
  bof2 <- paste(as.vector(tmp[,2]),as.vector(tmp[,5]),sep="x")
  inter <- intersect(bof1,bof2)
  #print(length(inter))
  
  AA <- AA[match(inter,bof1),]
  sift <- tmp[match(inter,bof2),]; rm(inter,tmp,bof1,bof2)
  dimnames(sift)[[2]][32] <- "*"
  dimnames(sift)[[2]][33] <- "_"
  AA2 <- as.matrix(AA)
  
  for(i in 1:nrow(sift))
  {
    AA2[i,-c(1:6)] <- as.vector(as.matrix(sift[i,match(as.vector(as.matrix(AA2[i,-c(1:6)])),dimnames(sift)[[2]])]))
  }
  AA2[AA=="*"] <- "*"
  write.table(AA2,paste("/group/jrigrp2/DiallelSofiane/AA/SIFTpred_Chr",chr,sep=""), row.names=F, quote=F, sep="\t")
}

