#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------
# Link Amino acids to MAPP predictions
# Script written by Sofiane mezmouk (Ross-Ibarra laboratory)
#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# 1/ Link aa to MAPP predictions (good or bad)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

tmp <- read.table("resMAPP-allpart", header=T, sep="\t")
  #tmp format 
Gene	Prot	Position	Good.Amino.Acids	Bad.Amino.Acids
GRMZM2G059865	GRMZM2G059865_P01	1	M	ACDEFGHIKLNPQRSTVWY
GRMZM2G059865	GRMZM2G059865_P01	2	A	CDEFGHIKLMNPQRSTVWY
GRMZM2G059865	GRMZM2G059865_P01	3	M	ACDEFGHIKLNPQRSTVWY
GRMZM2G059865	GRMZM2G059865_P01	4	A	CDEFGHIKLMNPQRSTVWY


tmp[tmp=="N/A"] <- NA
tmp <- tmp[!((is.na(tmp[,4]))&(is.na(tmp[,5]))),]
first <- read.table("ListeProtFirstTranscrit", header=F, sep="\t")
  #first format
#AC147602.5_FGP004
#AC148152.3_FGP001
#AC148152.3_FGP005
#AC148152.3_FGP006
#AC148152.3_FGP008


dim(tmp)
tmp <- tmp[as.vector(tmp[,2])%in%as.vector(first[,1]), ]
dim(tmp)
length(unique(as.vector(tmp[,1])))

for(chr in 1:10)
{
  print(chr)
    # AA is the output of SNPtoAA.r
  AA <- read.table(paste("NucToAA_Chr",chr,".txt",sep=""),sep="\t", header=T)
  bof1 <- paste(as.vector(AA[,2]),as.vector(AA[,3]),sep="x")
  bof2 <- paste(as.vector(tmp[,2]),as.vector(tmp[,3]),sep="x")
  
  inter <- intersect(bof1,bof2)
  print(length(inter))
  
  AA <- as.matrix(AA[match(inter,bof1),])
  mapp <- tmp[match(inter,bof2),]; rm(inter,bof1,bof2)
  
  AA2 <- as.matrix(AA)
  for(i in 1:nrow(mapp))
  {
    AA2[i,-c(1:6)][AA2[i,-c(1:6)]!="*"] <- unlist(lapply(as.vector(AA[i,-(1:6)][AA2[i,-c(1:6)]!="*"]),function(x){grep(x,as.vector(as.matrix(mapp[i,4:5])))}))
  }      
  #AA2[AA=="*"] <- "*"
  
  #  write.table(AA2,paste("MAPP-Allpred_Chr",chr,sep=""), row.names=F, quote=F, sep="\t")
  rm (AA,AA2,mapp)
}




#-----------------------------------------------------------------------------------
# 3/ Link aa to MAPP quantitative scores
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

for(chr in 1:10){
  mapp <- read.table(paste("resMAPP-Chr", chr, sep=""), header=T, sep="\t")
  mapp <- mapp[,c(1:3,35:54)]

  AA <- read.table(paste("NucToAA_Chr",chr,".txt",sep=""),sep="\t", header=T)

  bof1 <- paste(as.vector(AA[,2]),as.vector(AA[,3]),sep="x")
  bof2 <- paste(as.vector(mapp[,2]),as.vector(mapp[,3]),sep="x")
  inter <- intersect(bof1,bof2)

  AA <- as.matrix(AA[match(inter,bof1),])
  mapp <- as.matrix(mapp[match(inter,bof2),]); rm(inter,bof1,bof2)

  AA2 <- AA
  for(i in 1:nrow(mapp))
  {
    AA2[i,-c(1:6)][AA2[i,-c(1:6)]!="*"] <- unlist(lapply(as.vector(AA[i,-(1:6)][AA2[i,-c(1:6)]!="*"]),function(x){mapp[i,grep(paste(x,".1", sep=""),as.vector(dimnames(mapp)[[2]]))]}))
  }

  write.table(AA2,paste("MAPP-scores-_Chr",chr,sep=""), row.names=F, quote=F, sep="\t")
}



