#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------
# This script gives the amino acid position (first transcript) of SNPs in exons
# It Also creat input files for polydNdS and a summary file for processing the results of 
# This script was written to carry out a particular analysis; it may not be applied to another case without changes 
# A header of each file is put under every file used 
# Script written by Sofiane mezmouk (Ross-Ibarra laboratory)
#-----------------------------------------------------------------------------------
#XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOX
#-----------------------------------------------------------------------------------


  #Choose a chromosome to work with (1 to 10)
Chr <- 10
#--------------
library(stringr)
library(gdata)
is.odd <- function(x) x %% 2 != 0
#--------------
#-----------------------------
# A function that converts a nucleotide sequence to an amino acid sequence
transcript <- function(x){xtm <- as.vector(str_sub(paste(as.vector(x),collapse="",sep=""),debut,fin)); as.vector(codon[tapply(xtm,rep(1:length(xtm)),function(x){match(x,codon[,1])}),2])}
#-----------------------------
#--------------

gffall <- read.table("ZmB73_5b_FGS.gff")
  # gffall format
#9	ensembl	chromosome	1	156750706	.	.	.	ID=9;Name=chromosome:AGPv2:9:1:156750706:1
#9	ensembl	gene	66347	68582	.	-	.	ID=GRMZM2G354611;Name=GRMZM2G354611;biotype=protein_coding
#9	ensembl	mRNA	66347	68582	.	-	.	ID=GRMZM2G354611_T01;Parent=GRMZM2G354611;Name=GRMZM2G354611_T01;biotype=protein_coding
#9	ensembl	intron	68433	68561	.	-	.	Parent=GRMZM2G354611_T01;Name=intron.1
#9	ensembl	intron	67142	67886	.	-	.	Parent=GRMZM2G354611_T01;Name=intron.2


codon <- read.table("Codons.txt", header=T, sep="\t")
  #codon format 
#Codon	AA_1	AA_3	AA_Full	AntiCodon
#TCA	S	Ser	Serine	TGA
#TCG	S	Ser	Serine	CGA
#TCC	S	Ser	Serine	GGA
#TCT	S	Ser	Serine	AGA


genelist <- read.table("GeneProtNames", header=F, sep="\t")
  #genelist format 
#AC147602.5_FG004	AC147602.5_FGT004	AC147602.5_FGP004
#AC148152.3_FG001	AC148152.3_FGT001	AC148152.3_FGP001
#AC148152.3_FG005	AC148152.3_FGT005	AC148152.3_FGP005
#AC148152.3_FG006	AC148152.3_FGT006	AC148152.3_FGP006
#AC148152.3_FG008	AC148152.3_FGT008	AC148152.3_FGP008


transc <- as.vector(read.table("ListeProtFirstTranscrit", header=F, sep="\t")[,1])
  #transc format
#AC147602.5_FGP004
#AC148152.3_FGP001
#AC148152.3_FGP005
#AC148152.3_FGP006
#AC148152.3_FGP008


genelist <- genelist[as.vector(genelist[,3]%in%transc),]; rm(transc)
geneposi <- read.table("GenePositions", header=T, sep="\t")
  #geneposi format 
#Genes	Chr	Start	End
#GRMZM2G059865	1	4854	9652
#GRMZM5G888250	1	9882	10387
#GRMZM2G093344	1	109519	111769
#GRMZM2G093399	1	136307	138929


geneposi <- geneposi[geneposi[,2]==Chr,]
genelist <- genelist[as.vector(genelist[,1]) %in% as.vector(geneposi[,1]),]
#---------------
geno <- read.table(paste("282_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_imp95_1024_chr",Chr,".hmp.txt", sep=""), header=T, sep="\t")
  # geno format 
#rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	33-16	38-11	...
#S1_2111	C/T	1	2111	+	NA	NA	NA	NA	NA	NA	C	C	...
#S1_10097	C/G	1	10097	+	NA	NA	NA	NA	NA	NA	C	C	...
#S1_10390	G/A	1	10390	+	NA	NA	NA	NA	NA	NA	G	G	...

geno <- as.matrix(geno[,-c(1:2,4:10)])
geno[is.element(geno, c("M","K","S","W","Y","R","V","H","D","B","H"))] <- "N"


#---------------

# Result file to keep track between the real positions and the positions in the sequence used for polydNdS estimates
RespolydNdS <- matrix(NA,1,8, dimnames=list(NULL, c("gene","SNP","Chr","Position","SeqPosition","Sens","LengthCDS","NbSeq")))

# Result file with the amino acid polymorphisms corresponding to the 
resaa <- matrix(NA,1,(6+ncol(geno)),dimnames=list(NULL,c("gene","transcript","AAposition","SNP1","SNP2","SNP3","B73ref",dimnames(geno)[[2]][-1])))
problemes <- vector()

#---------------
#---------------
#Loop over gene
for(i in 1:nrow(geneposi)){
  if(nrow(geno[as.numeric(as.vector(geno[,1]))%in%c(geneposi[i,3]:geneposi[i,4]),,drop=F])>0){ # if I have SNPs in the gene{
    gff <- gffall[grep(geneposi[i,1],gffall[,9]),]
    posgene <- as.vector(c(geneposi[i,3]:geneposi[i,4]))
    posgene <- posgene[order(posgene)]
    SENStransc <- as.vector(gff[grep("gene",gff[,3]),7])
    posi <- gffall[grep(as.vector(genelist[match(geneposi[i,1],genelist[,1]),2]),gffall[,9]),]
    posi <- posi[grep("CDS",posi[,3]),,drop=F]
    
    CDS <- c(posi[1,4]:posi[1,5])
    if (nrow(posi)>1)
    {
      for (j in 2:nrow(posi))
      {
        CDS <- c(CDS,c(posi[j,4]:posi[j,5]))
      }
      rm(j)
    }
    CDS <- CDS[order(CDS)]
    rm(posi)
    #----------------
    if(nrow(geno[as.numeric(as.vector(geno[,1]))%in%CDS,,drop=F])>0){
      geneseq <- readLines(paste("gene",geneposi[i,1],".fasta",sep=""))
        # geneseq format for geneAC147602.5_FG004.fasta
      #>AC147602.5_FG004 seq=gene; coord=3:178846540..178848005:-1
      #ATGGAGATCGTCGCCACGCGCTCCCCGGCTTGCTGCGCCGCCGTGTCCTTCTCCCAGTCG
      #TACAGGCCCAAGGTACGTACGGCACCTTCATATCTCGTGACTACTGTACGTAAGCGGAAA
      #GTAGCAGCAGCTCGTCGCGCACACGTGCAGAAGCCTTAAGTTTGCTGATGATGTTGATGA
      
      geneseq <- paste(geneseq[-1],collapse="", sep="")
      geneseq <- strsplit(geneseq,split=character(1),fixed=T)[[1]]
      
      tprot <- readLines(paste("tprot_",genelist[as.vector(genelist[,1])==as.vector(geneposi[i,1]),3],".fasta",sep=""))
        #tprot format for tprot_AC147602.5_FGP004.fasta
      #>AC147602.5_FGP004 seq=translation; coord=3:178846540..178848005:-1; parent_transcript=AC147602.5_FGT004; parent_gene=AC147602.5_FG004
      #MEIVATRSPACCAAVSFSQSYRPKASRPPTTFYGESVRVNTARPLSARRQSKAASRAALS
      #ARCEIGDSLEEFLTKATPDKNLIRLLICMGEAMRTIAFKVRTASCGGTACVNSFGDEQLA
      #VDMLANKLLFEALEYSHVCKYACSEEVPELQDMGGPVEGS
      
      
      tprot <- paste(tprot[-1],collapse="",sep="")
      tprot <- strsplit(tprot, split = "", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
      
      # Creat the nucleotide sequenc of every genotype 
      if(SENStransc=="-"){
        sequ <- matrix(rep(geneseq,ncol(geno)), length(geneseq),ncol(geno), dimnames=list(rev(posgene),c("B73ref",dimnames(geno)[[2]][-1])))
      }else
      {
        sequ <- matrix(rep(geneseq,ncol(geno)), length(geneseq),ncol(geno), dimnames=list(posgene,c("B73ref",dimnames(geno)[[2]][-1])))
      }
      rm(geneseq)
      
      sequ <- sequ[as.numeric(dimnames(sequ)[[1]])%in%CDS,,drop=F]
      
      tmp <- geno[as.numeric(as.vector(geno[,1]))%in%CDS,, drop=F]
      dimnames(tmp)[[1]] <- as.numeric(as.vector(tmp[,1])); tmp <- tmp[,-1,drop=F]
      
      if(SENStransc=="-")
      {
        tmp2 <- tmp[,,drop=F]
        tmp[tmp2=="A"] <- "T";tmp[tmp2=="T"] <- "A";tmp[tmp2=="C"] <- "G";tmp[tmp2=="G"] <- "C"
        tmp[tmp2=="M"] <- "K";tmp[tmp2=="K"] <- "M";tmp[tmp2=="Y"] <- "R";tmp[tmp2=="R"] <- "Y"
        rm(tmp2)
      }
      
      for(j in 1:nrow(tmp))
      {
        bof <- tmp[j,tmp[j,]!="N",drop=F]
        sequ[match(dimnames(bof)[[1]],dimnames(sequ)[[1]]),match(dimnames(bof)[[2]],dimnames(sequ)[[2]])] <- bof
        rm(bof)
      }
      rm(j)
      
      #-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-
      # write an input file for polydNdS
      bofseq <- apply(sequ, 2, function(x){paste(as.vector(x),collapse="",sep="")})
      bofseq <-unique(bofseq)
      bof <- vector()
      bof[is.odd(1:(length(bofseq)*2))] <-  paste("> sequenceNumber",c(1:length(bofseq)),sep="")
      bof[!is.odd(1:(length(bofseq)*2))] <- bofseq
      writeLines(bof,paste("seq_",as.vector(geneposi[i,1]),".fasta", sep="")); rm(bof)
      #---------
      bof <- cbind(as.vector(geneposi[i,1]),paste("S",Chr,"_",as.numeric(as.vector(dimnames(tmp)[[1]])),sep=""),Chr,as.numeric(as.vector(dimnames(tmp)[[1]])),
                   match(dimnames(tmp)[[1]],dimnames(sequ)[[1]]),SENStransc,nrow(sequ),length(bofseq))
      dimnames(bof)[[2]] <- c("gene","SNP","Chr","Position","SeqPosition","Sens","LengthCDS","NbSeq")
      RespolydNdS <- rbind(RespolydNdS, bof); rm(bof,bofseq)
      
      #-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-
      # nucleotide to aa
      debut <- seq(1,nrow(sequ),3)
      fin <- pmin(debut+2,nrow(sequ))
      
      AA <- matrix(apply(sequ,2,transcript),ncol=ncol(geno), byrow=F)
      AA <- cbind(c(1:nrow(AA)),dimnames(sequ)[[1]][debut],dimnames(sequ)[[1]][(debut+1)],dimnames(sequ)[[1]][fin],AA)
      
      # Put a warning if the aa sequence I transcript is different from the aa sequence from the files I upload for the reference B73
      if(sum(as.numeric(as.vector(AA[,5])[1:length(tprot)]!=tprot),na.rm=T)!=0){
        problemes[length(problemes)+1] <- as.vector(geneposi[i,1])
        #print("!!!problem"); print(as.vector(as.matrix(genelistmp[ii,])))
      }
      
      AA <- AA[(as.numeric(AA[,2])%in%as.numeric(as.vector(geno[,1])))|(as.numeric(AA[,3])%in%as.numeric(as.vector(geno[,1])))|(as.numeric(AA[,4])%in%as.numeric(as.vector(geno[,1]))),,drop=F]
      if (nrow(AA)>0){
        AA <- cbind(as.vector(geneposi[i,1]),as.vector(genelist[as.vector(genelist[,1])==as.vector(geneposi[i,1]),3]),AA)
        dimnames(AA) <- list(NULL,c("gene","transcript","AAposition","SNP1","SNP2","SNP3",dimnames(sequ)[[2]]))
        resaa <- rbind(resaa,AA)
      }
      rm(AA,debut,fin,tprot,sequ)
    }
    rm(gff,SENStransc,CDS)
  }
}

resaa <- resaa[-1,]
RespolydNdS <- RespolydNdS[-1,]

if(length(problemes)>0){write.table(problemes,paste("Problemes_Chr",Chr,sep=""), sep="\t", row.names=F, quote=F, col.names=F)}
write.table(RespolydNdS, paste("SummaryPolydNdS.Chr",Chr,sep=""), sep="\t", quote=F, row.names=F)
write.table(resaa,paste("NucToAA_Chr",Chr,".txt",sep=""), sep="\t", row.names=F, quote=F)




