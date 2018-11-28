################################################################################################
#
#                                   POLYGEN REGENOTYPING CALCULATIONS
#
################################################################################################
library(stringr)
library(ggplot2)
library(pbapply)

#PolyGen is broken into two functions
#The first function (polyGen) combines the locus table and read counts into a single file and submits to the second function
#The second funtion (genoSetup) genotypes the data and returns the final genotypes to polyGen where it undergoes final formatting
polyGen<-function(locusTable,readCounts){
  #combine locus information and read counts into single table
  combinedData<-cbind(locusTable,readCounts)
  #pbapply has integrated progress bar, use apply if pbapply cannot be installed
  results<-pbapply(combinedData,1,genoSetup,epsilon=0.01)
  #results<-apply(combinedData,1,genoSetup,epsilon=0.01)
  results<-t(results)
  return(results)
}

genoSetup<-function(genoData,epsilon=0.01){
  readList<-genoData[4:length(genoData)]
  n_alleles<-strsplit(genoData[3],",")
  n_alleles<-length(n_alleles[[1]])
  #use numbers in place of actual alleles to allow haplotypes and other allele coding (indel, etc)
  #convert back to actual alleles later
  alleles=as.character(seq(1:n_alleles))
  ploidy=as.numeric(genoData[2])
  #set up possible genotypes
  alleleList<-substr(alleles,1,n_alleles)
  alleleList<-replicate(ploidy, alleleList) 
  #convert to vector, then sort numerically and convert back to text
  alleleList<-as.vector(alleleList)    
  alleleList<-sort(as.numeric(alleleList))
  alleleList<-as.character(alleleList)
  genoCombos<-combn(alleleList,ploidy)
  GenotypeList<-t(genoCombos)
  possibleGenotypes<-unique(GenotypeList)
  
  #make conversion table for genotypes
  realAlleles<-genoData[3]
  realAllelesOrder<-unlist(strsplit(realAlleles,",",perl=TRUE))
  
  #make table to convert numeric genotype codes to real genotypes
  numericGeno<-apply(possibleGenotypes,1,function(x) paste(x,collapse=""))
  realGenotypes<-apply(possibleGenotypes,2,function(x) realAllelesOrder[as.numeric(x)])
  realGeno<-apply(realGenotypes,1,function(x) paste(x,collapse=","))
  genoConvert<-as.data.frame(t(rbind(numericGeno,realGeno)))
  
  
  #make matrix of relative dosage for each possible genotype
  relative_dosage<-matrix(NA,nrow=dim(possibleGenotypes)[1],ncol=n_alleles)
  
  for(i in 1:n_alleles){
    proportion<-sapply(numericGeno,function(x) str_count(x, alleles[i]))
    proportion<-proportion/ploidy
    relative_dosage[,i]<-proportion
  }
  
  #make matrix of read chances
  read_chances<-matrix(NA,nrow=1,ncol=ploidy)
  #read_chances<-(relative_dosage*(1-epsilon) + (1-relative_dosage)*epsilon)
  #updated error model dividing epsilon by 3 (n alleles -1) where n alleles is the number of possible bases (4, ATCG)
  read_chances<-(relative_dosage*(1-epsilon) + (1-relative_dosage)*epsilon/3)
  
  #make likelihood function
  likelihoodFunc<-function(reads){
    #convert reads to alleles
    #parse allele counts
    if(reads=="."){
      reads<-paste(as.character(replicate(n_alleles, 0)),collapse=",")
    }
    alleleCounts<-as.numeric(unlist(strsplit(reads,",")))
    initLikelihood<-log(1)
    likelihoodMatrix<-matrix(0,nrow=dim(possibleGenotypes)[1],ncol=ploidy)
    #likelihood calculation
    likelihoodMatrix<-t(initLikelihood+log(t(read_chances))*(alleleCounts))
    likelihood<-apply(likelihoodMatrix,1,sum)
    genotypes<-apply(possibleGenotypes,1,function(x) paste(x,collapse=""))
    
    #create likelihood results matrix
    like_of_geno<-as.data.frame(matrix(NA,nrow=dim(possibleGenotypes)[1],ncol=2))
    colnames(like_of_geno)<-c("genotype","likelihood")
    like_of_geno$genotype<-genotypes
    like_of_geno$likelihood<-likelihood
    
    #use a likelihood ratio test to determine the support for the 'best' genotype
    #order likelihoods
    geno_likes<-like_of_geno[order(-like_of_geno$likelihood),]
    #rename rows so they can be extracted in correct order
    rownames(geno_likes)<-1:nrow(geno_likes)
    #compare likelihood ratio of two most likely models
    LR<-2*(geno_likes[1,2]-geno_likes[2,2])
    #get p-value of likelihood ratio
    p<-1-pchisq(LR, 1)
    if(p<0.05){
      genoResult<-geno_likes[1,1]
      genoResult<-genoConvert$realGeno[genoConvert$numericGeno==genoResult]
    }else{
      genoResult<-"0"
    }
    result<-as.character(genoResult)
    return(result)
  }
  
  likelihoods<-sapply(readList,likelihoodFunc)
  return(likelihoods)
}


################################################################################################
#
#                                         FILE CONVERSION
#
################################################################################################
exportGenepop<-function(polygenResults,locusTable,exportParalogs=FALSE,filename="polygenResults.genepop"){
  #use ploidy and allele information from locusTable for genotype conversion and paralog filtering
  #combine locusTable and polyGenResults
  combinedData<-cbind(locusTable,polygenResults)
  #function to convert alleles to numeric code
  alleleConvert<-function(combinedData){
    alleles<-unlist(strsplit(as.character(combinedData[3]),","))
    genotypes<-combinedData[4:length(combinedData)]
    if(nchar(alleles[1])==1){
      numericAlleles<-str_replace_all(genotypes,"0",paste(rep("0",as.numeric(combinedData[2])*2),collapse=""))
      numericAlleles<-str_replace_all(numericAlleles,c("A"="01","C"="02","G"="03","T"="04","-"="05"))
      numericAlleles<-str_replace_all(numericAlleles,",","")
    }else{
      hapAlleles<-unlist(strsplit(as.character(combinedData[3]),","))
      hapAlleleNumbers<-seq(1,length(hapAlleles),by=1)
      #check if number allele for haplotype is already same length as ploidy, if not pad with a 0
      hapAlleleNumbers<-ifelse(hapAlleleNumbers<10,paste("0",hapAlleleNumbers,sep=""),hapAlleleNumbers)
      hapAlleleCode<-hapAlleleNumbers
      names(hapAlleleCode)<-hapAlleles
      numericAlleles<-str_replace_all(genotypes,"0",paste(rep("0",as.numeric(combinedData[2])*2),collapse=""))
      numericAlleles<-str_replace_all(numericAlleles,hapAlleleCode)
      numericAlleles<-str_replace_all(numericAlleles,",","")
    }
    return(numericAlleles)
  }
  #convert bases to genepop code
  genepopBases<-apply(combinedData,1,alleleConvert)
  rownames(genepopBases)<-colnames(polygenResults)
  rownames(genepopBases)<-paste(rownames(genepopBases),",",sep="")
  #exclude paralogs if exportParalogs flag is FALSE
  if(exportParalogs==FALSE){
    genepopBases<-genepopBases[,colnames(genepopBases) %in% locusTable$Locus_ID[locusTable$ploidy==2]]
  }
  #write output file
  write("Title Line:",filename)
  write(colnames(genepopBases),filename,ncolumns=1,append=TRUE)
  write("Pop",filename,append=TRUE)
  write.table(genepopBases,filename,quote=FALSE,sep="\t",col.names=FALSE,append=TRUE)
}

exportRubias<-function(polygenResults,locusTable,sampleMetaData=NULL,filename="polygenResults_rubias.txt"){
  #use ploidy and allele information from locusTable for genotype conversion and paralog filtering
  #combine locusTable and polyGenResults
  combinedData<-cbind(locusTable,polygenResults)
  #replace missing genotypes with NA
  combinedData[combinedData==0]<-NA
  #exclude loci that are paralogs, these are not currently accepted by rubias
  combinedData<-combinedData[combinedData$ploidy==2,]
  #remove Locus_ID, ploidy, and alleles columns and transpose
  genotypes<-t(combinedData[,4:dim(combinedData)[2]])
  #split genotypes into alleles
  splitGenotypes<-as.data.frame(apply(genotypes,2,function(x) data.frame(str_split_fixed(x,",",2))),check.names=FALSE)
  #format column names
  colnames(splitGenotypes)<-str_replace_all(colnames(splitGenotypes),c(".X1",".X2"),c("",".1"))
  #str_replace_all replaces NA with "", fill in missing NAs
  splitGenotypes[splitGenotypes==""]<-NA
  rownames(splitGenotypes)<-rownames(genotypes)
  
  #make new columns for sample type, reporting group, collection, and indivdual
  if(is.null(sampleMetaData)){
    rubiasFile<-data.frame(matrix(NA, nrow=dim(splitGenotypes)[1],ncol=4))
    colnames(rubiasFile)<-c("sample_type","repunit","collection","indiv")
    rubiasFile$indiv<-rownames(genotypes)
    rubiasFile<-cbind(rubiasFile,splitGenotypes)
  }else{
    splitGenotypes$indiv<-rownames(genotypes)
    #combine sample data with genotypes, retain all genotyped samples but include sample data only from genotyped samples
	rubiasFile<-merge(sampleMetaData,splitGenotypes,all.x=FALSE,all.y=TRUE)
    #reorder columns to fit rubias format
    rubiasFile<-rubiasFile[c(2:4,1,5:dim(rubiasFile)[2])]
  }
  write.table(rubiasFile,filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}



################################################################################################
#
#                                         DATA SUMMARIES
#
################################################################################################

#CALCULATE AVERAGE READ DEPTH, GENOTYPE RATE, MINOR ALLELE FREQUENCY, AND MAJOR ALLELE FREQUENCY
summarizeGTscore<-function(alleleReads, locusTable, genotypes){
	#calculate average read depth
	reads<-apply(alleleReads, 1:2, function(x) sum(as.numeric(unlist(str_split(x,",")))))
	totalReads<-apply(reads, 1, sum)
	#Average read depth is calculated by dividing my number of samples with non-zero reads
	#Is this appropriate?
	nonZeroSamples<-apply(reads, 1, function(x) {length(x)-sum(x=="0")})
	avgReads<-totalReads/nonZeroSamples
	calcReadDepth_results<-data.frame(keyName=names(avgReads), value=avgReads, row.names=NULL)
	colnames(calcReadDepth_results)<-c("Locus_ID","AvgReadDepth")
	
	#calculate genotype rate
	genotypeRate_results<-apply(genotypes,1,function(x){(length(x)-sum(x=="0"))/length(x)})
	genotypeRate_results<-data.frame(keyName=names(genotypeRate_results), value=genotypeRate_results,row.names=NULL)
	colnames(genotypeRate_results)<-c("Locus_ID","GenotypeRate")
	
	#calculate minor allele frequency and major allele frequency
	getAlleleFreqs<-function(alleleFreqsCombinedData){
		Locus_ID<-as.character(alleleFreqsCombinedData[1])
		#get unique alleles
		alleles<-unlist(str_split(unlist(alleleFreqsCombinedData[3]),","))
		#get counts and frequency for each unique allele
		alleleCounts<-unlist(lapply(alleles, function(x) sum(str_count(as.character(unlist(alleleFreqsCombinedData[4:length(alleleFreqsCombinedData)])),x))))
		alleleFreqs<-alleleCounts/sum(alleleCounts)
		#get minimum frequency and maximum frequency
		minFreq<-min(alleleFreqs)
		maxFreq<-max(alleleFreqs)
		#combine all frequencies into comma delimited string
		allFreqs<-paste(round(alleleFreqs,digits=2),collapse=",")
		#combine and return minFreq,maxFreq,allFreqs
		freqResults<-paste(Locus_ID,minFreq,maxFreq,alleleFreqsCombinedData[3],allFreqs,sep=" ")
		return(freqResults)
	}
	alleleFreqsCombinedData<-cbind(locusTable,genotypes)
	alleleFreqs_results<-as.data.frame(apply(alleleFreqsCombinedData,1,getAlleleFreqs),row.names=NULL)
	alleleFreqs_results<-data.frame(str_split_fixed(alleleFreqs_results[,1]," ",5))
	colnames(alleleFreqs_results)<-c("Locus_ID","minFreq","maxFreq","alleles","allFreqs")
	#alleleFreqs_results
	
	#combine results and return
	depth_geno<-merge(calcReadDepth_results, genotypeRate_results)
	all_combined<-merge(depth_geno,alleleFreqs_results)
	#ensure appropriate columns are numeric
	all_combined$AvgReadDepth<-as.numeric(as.character(all_combined$AvgReadDepth))
	all_combined$GenotypeRate<-as.numeric(as.character(all_combined$GenotypeRate))
	all_combined$minFreq<-as.numeric(as.character(all_combined$minFreq))
	all_combined$maxFreq<-as.numeric(as.character(all_combined$maxFreq))
	#set final column names
	colnames(all_combined)<-c("Locus_ID","AvgReadDepth","GenotypeRate","minAF","majAF","alleles","allFreqs")
	return(all_combined)
}

#calculate sample genotype rate
sampleGenoRate<-function(genotypes){
	#missing genotypes are coded as 0
	#calculate genotype rate as proportion of non-zero genotypes for each locus
	genotypeRate_results<-apply(genotypes,2,function(x){(length(x)-sum(x=="0"))/length(x)})
	genotypeRate_results<-data.frame(keyName=names(genotypeRate_results), value=genotypeRate_results,row.names=NULL)
	colnames(genotypeRate_results)<-c("sample","GenotypeRate")
	return(genotypeRate_results)
}

#PLOT GENOTYPE RESULTS
#This code works for bi-allelic data only
#Plotting ratios for more than two alleles requires three-dimensional space,
#Plotting in three-dimensional space is possible but requires significant revision

plotGenotypes<-function(locusTable,alleleReads,genotypes,type="ratio",savePlot="N",saveDir=""){
	#combine data to allow running through apply
	#add a column specifying the number of columns for each section
	#genotypes, total allele reads, and allele ratios should all have the same number of columns
	sectionColumns<-dim(genotypes)[2]
	combinedData<-cbind(locusTable,sectionColumns,genotypes,alleleReads)
	#use invisible to suppress stdout from dev.off() while allowing errors to print
	invisible(pbapply(combinedData,1,plotGenotypes_setup,type,savePlot,saveDir))
}

plotGenotypes_setup<-function(genoData,type,savePlot,saveDir){
	#separate data into categories
	locusTable<-genoData[1:3]
	sectionColumns<-as.numeric(genoData[4])
	genotypesStart<-5
	genotypesEnd<-5+sectionColumns-1
	alleleReadsStart<-genotypesEnd+1
	alleleReadsEnd<-alleleReadsStart+sectionColumns-1
	genotypes<-genoData[genotypesStart:genotypesEnd]
	alleleReads<-genoData[alleleReadsStart:alleleReadsEnd]
	#split reads by allele
	alleleReads_1<-str_split_fixed(alleleReads,",",2)[,1]
	alleleReads_2<-str_split_fixed(alleleReads,",",2)[,2]
	#combine results into data frame
	combinedData<-data.frame(genotypes,alleleReads_1,alleleReads_2,stringsAsFactors=FALSE)
	combinedData$alleleReads_1<-as.numeric(combinedData$alleleReads_1)
	combinedData$alleleReads_2<-as.numeric(combinedData$alleleReads_2)
	#calculate total reads
	combinedData$totalReads<-combinedData$alleleReads_1+combinedData$alleleReads_2
	#calculate read ratio
	combinedData$ratios<-combinedData$alleleReads_1/combinedData$totalReads
	
	#calculate summary statistics
	genotypeRate<-(length(combinedData$genotypes)-sum(combinedData$genotypes=="0"))/length(combinedData$genotypes)
	genotypeRate<-round(genotypeRate, digits=2)
	averageReadDepth<-mean(combinedData$totalReads)
	averageReadDepth<-round(averageReadDepth, digits=2)

	#plot results, exclude ratios with NA from plot to prevent warnings  
	locusID<-locusTable[1]
	summaryText=paste("average depth: ",averageReadDepth,"					","genotype_rate: ", genotypeRate,sep=" ")
	
	if(type=="ratio"){
		genoPlot<-ggplot()+geom_histogram(data=combinedData,aes(ratios,color=genotypes, fill=genotypes),binwidth=0.01,alpha=0.5,na.rm=TRUE)+
		xlim(-0.01,1.01)+labs(title=locusID,x="Allele Ratio",y="Frequency",subtitle=summaryText)+
		theme_bw()+theme(plot.title=element_text(hjust=0.5),plot.subtitle=element_text(hjust=0.5))
	}else if(type=="scatter"){
		genoPlot<-ggplot()+geom_point(data=combinedData,aes(x=alleleReads_1,y=alleleReads_2,color=genotypes))+
		xlim(range(combinedData$alleleReads_1,combinedData$alleleReads_2))+
		ylim(range(combinedData$alleleReads_1,combinedData$alleleReads_2))+
		labs(title=locusID,x="Allele 1 Reads",y="Allele 2 Reads",subtitle=summaryText)+coord_fixed(ratio=1)+
		theme_bw()+theme(plot.title=element_text(hjust=0.5),plot.subtitle=element_text(hjust=0.5))
	}
	#save plot if specified
	if(savePlot=="Y"){
		JPGoutfile<-paste(saveDir,"/",locusID,"_",type,".jpg",sep="")
		jpeg(file=JPGoutfile)
		print(genoPlot)
		#capture.output(dev.off(),file=NULL)
		dev.off()
		
	}else{
		print(genoPlot)
	}
}



summarizeMismatches<-function(mismatchData,saveDir){
  for(i in 1:dim(mismatchData)[1]){
    locusID<-mismatchData$Locus[i]
    locusSeq<-mismatchData$RefSeq[i]
    locusMismatches<-mismatchData$Mismatches[i]
    mismatchByPositions<-as.numeric(unlist(str_split(locusMismatches,",")))
    #mismatchByPositions<-as.numeric(mismatchByPositions)
    mismatchPositions<-as.numeric(seq(1,length(mismatchByPositions),by=1))
    #mismatchPositions<-as.numeric(mismatchPositions)
    seqBases<-unlist(str_split(locusSeq,""))
    mismatchPlotData<-data.frame(cbind(mismatchPositions,mismatchByPositions))
    locusLabel<-paste("Locus",locusID,sep=" ")
    annotateYpos<-(-1*max(mismatchByPositions)/100)
    plot<-ggplot()+geom_point(data=mismatchPlotData,aes(x=mismatchPositions,y=mismatchByPositions))+
      annotate("text",x=mismatchPositions,y=annotateYpos,label=seqBases,hjust=0.5,vjust=1,family="mono")+
      ggtitle(locusLabel)+xlab("Position")+ylab("Mismatches")+theme_bw()+theme(plot.title=element_text(hjust=0.5))+
      geom_line(data=mismatchPlotData,aes(x=mismatchPositions,y=mismatchByPositions))
    
    JPGoutfile<-paste(saveDir,"/",locusID,".jpg",sep="")
    jpeg(file=JPGoutfile, width=960, height=960)
    print(plot)
    dev.off()
  }
}



#Generate alignments of reference sequences, amplified sequences, primers, and probes
library(msa)
library(tools)
library(plyr)
library(utils)

alignMatchedSeqs<-function(referenceSeqs=NULL,primerProbes,matchedReads,minReads=1,maxAlignedSeqs=100,type="primerProbe",saveDir=""){
  outDir=paste(saveDir,"/",sep="")
  #loci<-unique(primerProbes$Locus)
  loci<-primerProbes
  for (x in 1:dim(loci)[1]){
    locus<-loci$Locus[x]
    if(type=="primerProbe"){
		locusSNP<-paste(loci$Locus[x],loci$SNP[x],sep="_")
		locusReads<-matchedReads[matchedReads$Locus==locusSNP,]
		#make new column of combined locus name and SNP position
		primerProbes$LocusSNP<-paste(primerProbes$Locus,primerProbes$SNPpos,sep="_")
		#filter to retain probes from one locus
		probes<-primerProbes[primerProbes$LocusSNP==locusSNP,]
		#convert problem characters in locus name
		safeLocusName<-gsub("-","_",locusSNP)
		safeLocusName<-gsub("\\.","_",safeLocusName)
		safeLocusName<-gsub("'","prime",safeLocusName)
	}else if(type=="primer"){
		locusReads<-matchedReads[matchedReads$Locus==locus,]
		#filter to retain probes from one locus
		probes<-primerProbes[primerProbes$Locus==locus,]
		#convert problem characters in locus name
		safeLocusName<-gsub("-","_",locus)
		safeLocusName<-gsub("\\.","_",safeLocusName)
		safeLocusName<-gsub("'","prime",safeLocusName)
	}
	
	
	if(dim(locusReads)[1]>0){
		#filter to retain unique reads and sum counts for each unique read
		uniqueLocusReads<-aggregate(Count ~ Sequence, data=locusReads, sum)
	}else{
		if(type=="primerProbe"){
			noReads<-paste(locusSNP,"No Reads", sep=" ")
		}else if(type=="primer"){
			noReads<-paste(locus,"No Reads", sep=" ")
		}
		print(noReads)
		next
	}
	if(max(uniqueLocusReads$Count>=minReads)){
		uniqueLocusReads<-uniqueLocusReads[uniqueLocusReads$Count>=minReads,]
		#restrict MSA alignment to N unique sequences ranked by number of reads (high to low)
		uniqueLocusReads<-head(uniqueLocusReads[order(-uniqueLocusReads$Count),], n=maxAlignedSeqs)
		
		#make fasta file of retained sequences
		makeFasta<-function(refSeq=NULL,sequences,probes,safeLocusName){
		  #initialize fasta output file
		  fastaOut<-paste(outDir,safeLocusName,".fasta",sep="")
		  cat(NULL,file=fastaOut,append=FALSE)
		  if(!is.null(refSeq)){
			  #add reference sequence to fasta file
			  cat(">Reference",as.character(refSeq$refSeq),sep="\n",file=fastaOut,append=TRUE)
		  }
		  #add matched sequence reads to fasta file
		  for (i in 1:dim(sequences)[1]){
			#locusName<-paste(">",sequences$Locus[i],"Reads: ", sequences$Count[i],sep="")
			locusName<-paste(">", sequences$Count[i], "Reads",sep=" ")
			locusSequence<-as.character(sequences$Sequence[i])
			cat(locusName,locusSequence,sep="\n",file=fastaOut,append=TRUE)
		  }
		  for (i in 1:dim(probes)[1]){
			probe1Name<-">probe1"
			probe2Name<-">probe2"
			#primer<-as.character(probes$Primer[i])
			probe1<-as.character(probes$Probe1[i])
			probe2<-as.character(probes$Probe2[i])
			#replace any brackets in probe with IUB code
			probe1<-str_replace_all(probe1,c("\\[AG\\]"="R","\\[GA\\]"="R","\\[CT\\]"="Y","\\[TC\\]"="Y","\\[GT\\]"="K","\\[TG\\]"="K","\\[AC\\]"="M","\\[CA\\]"="M","\\[GC\\]"="S","\\[CG\\]"="S","\\[AT\\]"="W","\\[TA\\]"="W"))
			probe2<-str_replace_all(probe2,c("\\[AG\\]"="R","\\[GA\\]"="R","\\[CT\\]"="Y","\\[TC\\]"="Y","\\[GT\\]"="K","\\[TG\\]"="K","\\[AC\\]"="M","\\[CA\\]"="M","\\[GC\\]"="S","\\[CG\\]"="S","\\[AT\\]"="W","\\[TA\\]"="W"))
			
			cat(probe1Name,probe1,sep="\n",file=fastaOut,append=TRUE)
			cat(probe2Name,probe2,sep="\n",file=fastaOut,append=TRUE)
		  }
		}
		
		if(is.null(referenceSeqs)){
			makeFasta(sequences=uniqueLocusReads,probes=probes,safeLocusName=safeLocusName)
		}else{
			#filter to retain reference sequence from one locus
			refSeq<-referenceSeqs[as.character(referenceSeqs$Locus)==as.character(locus),]
			#force uppercase to prevent errors in sequence alignment
			refSeq$refSeq<-toupper(refSeq$refSeq)
			
			makeFasta(refSeq,uniqueLocusReads,probes,safeLocusName)
		}
		
		#load fasta file of retained sequences
		locusSequences<-readAAStringSet(paste(outDir,safeLocusName,".fasta",sep=""))
		#align sequences
		#use capture.output to suppress stdout from msa alignment
		capture.output(locusAlignments<-msa(locusSequences),file=NULL)
		#convert alignments to dataframe
		seqCount<-names(locusAlignments@unmasked)
		seq<-paste(locusAlignments)
		seqDF<-data.frame(seqCount, seq)
		#head(seqDF)
		
		#write dataframe output
		dfOut<-paste(outDir,safeLocusName,".txt",sep="")
		write.table(seqDF,dfOut,sep="\t",quote=FALSE,row.names=FALSE)
		
		#call formatMSAalignments.pl to convert file to conditionally formatted excel file
		formatMSAcmd<-paste("perl formatMSAalignments.pl", dfOut, sep=" ")
		system(formatMSAcmd)
		
		#report that locus sequence alignment was completed
		if(type=="primerProbe"){
			completed<-paste(locusSNP,"completed", sep=" ")
		}else if(type=="primer"){
			completed<-paste(locus,"completed", sep=" ")
		}
		print(completed)
	  }else{
	  	if(type=="primerProbe"){
			noReads<-paste(locusSNP,"No Reads", sep=" ")
		}else if(type=="primer"){
			noReads<-paste(locus,"No Reads", sep=" ")
		}
		print(noReads)
	  }
	}
}

