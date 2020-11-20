################################################################################################
#
#                                   POLYGEN REGENOTYPING CALCULATIONS
#
################################################################################################
library(stringr)
library(ggplot2)
library(pbapply)
library(dplyr)

#PolyGen is broken into two functions
#The first function (polyGen) combines the locus table and read counts into a single file and submits to the second function
#The second funtion (genoSetup) genotypes the data and returns the final genotypes to polyGen where it undergoes final formatting
polyGen<-function(locusTable,readCounts,p_thresh=0.05,epsilon=0.01){
  #remove correction factor column from locusTable (if present)
  if("correctionFactors" %in% colnames(locusTable)){
    locusTable<-locusTable %>% select(-correctionFactors)
  }

  #combine locus information and read counts into single table
  combinedData<-cbind(locusTable,readCounts)
  #pbapply has integrated progress bar, use apply if pbapply cannot be installed
  results<-pbapply(combinedData,1,genoSetup,epsilon,p_thresh)
  #results<-pbapply(combinedData,1,genoSetup,epsilon=0.01)
  #results<-apply(combinedData,1,genoSetup,epsilon=0.01)
  results<-t(results)
  return(results)
}

genoSetup<-function(genoData,epsilon=0.01,p_thresh=0.05){
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
  #print(genoCombos)
  GenotypeList<-t(genoCombos)
  possibleGenotypes<-unique(GenotypeList)
  
  #make conversion table for genotypes
  realAlleles<-genoData[3]
  realAllelesOrder<-unlist(strsplit(realAlleles,",",perl=TRUE))
  
  #make table to convert numeric genotype codes to real genotypes
  numericGeno<-apply(possibleGenotypes,1,function(x) paste(x,collapse=","))
  realGenotypes<-apply(possibleGenotypes,2,function(x) realAllelesOrder[as.numeric(x)])
  realGeno<-apply(realGenotypes,1,function(x) paste(x,collapse=","))
  genoConvert<-as.data.frame(t(rbind(numericGeno,realGeno)))
  
  #function to generate allele dosage for each genotype
  generateDosage<-function(genos){
    #split genotype into each allele for dosage assignment
    #genoAlleles<-str_split_fixed(genos,",",2)
    genoAlleles<-str_split(genos,",",simplify=TRUE)
    
    #function to count number of times each possible allele is present in each possible genotype
    genoDosage<-function(singleGeno){
      #test vector of all possible alleles to see which allele matches the genotype
      #do this individually for each allele in the genotype
      matches<-t(sapply(singleGeno,function(x) as.numeric(alleles==x)))
      #sum vector of possible allele matches for each allele in genotype
      #this gives count of each possible allele in the genotype being tested
      matches<-apply(matches,2,function(x) sum(x))
      return(matches)
    }
    
    #count number of times each possible allele is present in each possible genotype
    genoMatches<-t(apply(genoAlleles,1,genoDosage))
    #divide count matrix by ploidy to get relative dosage of each allele
    dosage<-genoMatches/ploidy
    return(dosage)
  }
  
  relative_dosage<-generateDosage(numericGeno)
  
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
    #genotypes<-apply(possibleGenotypes,1,function(x) paste(x,collapse=""))
    genotypes<-apply(possibleGenotypes,1,function(x) paste(x,collapse=",")) #this is original for reference
    
    #create likelihood results matrix
    like_of_geno<-as.data.frame(matrix(NA,nrow=dim(possibleGenotypes)[1],ncol=2))
    colnames(like_of_geno)<-c("genotype","likelihood")
    like_of_geno$genotype<-genotypes
    like_of_geno$likelihood<-likelihood
    
    #print(like_of_geno)
    
    #use a likelihood ratio test to determine the support for the 'best' genotype
    #order likelihoods
    geno_likes<-like_of_geno[order(-like_of_geno$likelihood),]
    #rename rows so they can be extracted in correct order
    rownames(geno_likes)<-1:nrow(geno_likes)
    #compare likelihood ratio of two most likely models
    LR<-2*(geno_likes[1,2]-geno_likes[2,2])
    #get p-value of likelihood ratio
    p<-1-pchisq(LR, 1)
    #if(p<0.05){
    if(p<p_thresh){
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
exportGenepop<-function(polygenResults,locusTable,exportParalogs=FALSE,sampleWhitelist=NULL,locusWhitelist=NULL,sampleBlacklist=NULL,locusBlacklist=NULL,filename="polygenResults.genepop"){
  #use ploidy and allele information from locusTable for genotype conversion and paralog filtering
  #remove correction factor column from locusTable (if present)
  if("correctionFactors" %in% colnames(locusTable)){
    locusTable<-locusTable %>% select(-correctionFactors)
  }
  
  #filter locusTable and polyGenResults based on locus and sample whitelists and blacklists
  if(!is.null(locusWhitelist)){
    locusTable<-locusTable %>% filter(Locus_ID %in% locusWhitelist)
    polygenResults<-polygenResults[rownames(polygenResults) %in% locusWhitelist,]
  }
  if(!is.null(locusBlacklist)){
    locusTable<-locusTable %>% filter(!Locus_ID %in% locusBlacklist)
    polygenResults<-polygenResults[!rownames(polygenResults) %in% locusBlacklist,]
  }
  if(!is.null(sampleWhitelist)){
    polygenResults<-polygenResults[,colnames(polygenResults) %in% sampleWhitelist]
  }
  if(!is.null(sampleBlacklist)){
    polygenResults<-polygenResults[,!colnames(polygenResults) %in% sampleBlacklist]
  }
  
  #combine locusTable and polyGenResults
  combinedData<-cbind(locusTable,polygenResults)
  #print(combinedData[1:5,1:5])
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
  colnames(genepopBases)<-combinedData$Locus_ID
  #print (genepopBases[1:5,1:5])
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

exportRubias<-function(polygenResults,locusTable,sampleMetaData=NULL,sampleWhitelist=NULL,locusWhitelist=NULL,sampleBlacklist=NULL,locusBlacklist=NULL,filename="polygenResults_rubias.txt"){
  #use ploidy and allele information from locusTable for genotype conversion and paralog filtering
  #remove correction factor column from locusTable (if present)
  if("correctionFactors" %in% colnames(locusTable)){
    locusTable<-locusTable %>% select(-correctionFactors)
  }
  
  #filter locusTable and polyGenResults based on locus and sample whitelists and blacklists
  if(!is.null(locusWhitelist)){
    locusTable<-locusTable %>% filter(Locus_ID %in% locusWhitelist)
    polygenResults<-polygenResults[rownames(polygenResults) %in% locusWhitelist,]
  }
  if(!is.null(locusBlacklist)){
    locusTable<-locusTable %>% filter(!Locus_ID %in% locusBlacklist)
    polygenResults<-polygenResults[!rownames(polygenResults) %in% locusBlacklist,]
  }
  if(!is.null(sampleWhitelist)){
    polygenResults<-polygenResults[,colnames(polygenResults) %in% sampleWhitelist]
    if(!is.null(sampleMetaData)){
      sampleMetaData<-sampleMetaData %>% filter(indiv %in% sampleWhitelist)
    }
  }
  if(!is.null(sampleBlacklist)){
    polygenResults<-polygenResults[,!colnames(polygenResults) %in% sampleBlacklist]
    if(!is.null(sampleMetaData)){
      sampleMetaData<-sampleMetaData %>% filter(!indiv %in% sampleBlacklist)
    }
  }
  
  
  #combine locusTable and polyGenResults
  combinedData<-cbind(locusTable,polygenResults)
  rownames(combinedData)<-rownames(polygenResults)
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
	#remove correction factor column from locusTable (if present)
	if("correctionFactors" %in% colnames(locusTable)){
		locusTable<-locusTable %>% select(-correctionFactors)
	}
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
	
	####################################################################################################
	
	#calculate contamination score, this applies a binomial test to read counts from heterozygous genotypes
	#the contamination score is the proportion of heterozygous genotypes that failed the binomial test for each sample
	
	#function to create hetMatrix, 1 is heterozygous, 0 is homozygous
	idHets<-function(genos){
	  alleles<-str_split(genos,",")
	  hetCount<-unlist(lapply(alleles,function(x) length(unique(x))))
	  hetCount[hetCount!=2]<-0
	  hetCount[hetCount==2]<-1
	  return(hetCount)
	}

	#function to run binomial test
	binomTest<-function(reads){
	  ########binomial test doesn't work for haplotypes as written
	  #should be able to work for haplotypes as long as locus is diploid
	  #work on fixing, for now check if there are more than two alleles and set all reads to 0 if that is the case_when
	  #can't directly check for more than two alleles because splitting alleles previously into allele1 and allele2 forced
	  #counts into two columns.
	  #if more than two alleles, conversion to number in the second column results in NA.
	  #check for NA and if present reset all A1 and A2 in this code to 0
	  A1<-as.numeric(str_split_fixed(reads,",",2)[,1])
	  A2<-as.numeric(str_split_fixed(reads,",",2)[,2])
	  A1NA<-sum(is.na(A1))
	  A2NA<-sum(is.na(A2))

	  if(A1NA>0|A2NA>0){
	    A1<-0
		A2<-0
	  }

	  tot<-A1+A2
	  
	  if(tot!=0){
		result<-binom.test(A1, tot, 0.5, alternative=c("two.sided"), conf.level = 0.95)
		return(result$p.value)
	  }else{
		return(NA)
	  }
	}

	#function to manage binomial test
	testHet<-function(reads){
	  binomTest_results<-unlist(lapply(reads,binomTest))
	  numOutlier<-sum(binomTest_results<=0.05,na.rm=TRUE)
	  numHet<-sum(!is.na(binomTest_results))
	  #result is number of genotypes that failed binomial test over total number of heterozygous genotypes for each sample
	  result<-numOutlier/numHet
	  return(result)
	}

	#make matrix where heterozygous genotypes are 1 and homozygous genotypes are 0
	hetMatrix<-apply(genotypes,1,function(x){idHets(x)})

	#######################################################
	#this step causes problems if loci have more than two alleles, work on fixing...

	#split reads into allele 1 and allele 2
	allele1<-apply(alleleReads,1,function(x) as.numeric(str_split_fixed(x,",",2)[,1]))
	allele2<-apply(alleleReads,1,function(x) as.numeric(str_split_fixed(x,",",2)[,2]))

	#multiply by het matrix to keep only counts from heterozygous genotypes
	allele1HetCounts<-allele1*hetMatrix
	allele2HetCounts<-allele2*hetMatrix
	
	
	#recombine counts so only heterozygous genotypes from original count data are retained
	combinedHetCounts<-paste(allele1HetCounts,allele2HetCounts,sep=",")
	#reset dimensions and reassign column names
	dim(combinedHetCounts)<-c(dim(allele1HetCounts)[1],dim(allele1HetCounts)[2])
	colnames(combinedHetCounts)<-colnames(allele1HetCounts)


	#run function to get heterozygous contamination score for all samples
	hetContaminationResults<-apply(combinedHetCounts,2,testHet)
	hetContaminationResults<-data.frame(Locus_ID=names(hetContaminationResults),conScore=hetContaminationResults,stringsAsFactors=FALSE,row.names=NULL)
	
	
	####################################################################################################
	
	#combine results and return
	depth_geno<-merge(calcReadDepth_results, genotypeRate_results)
	all_combined<-merge(depth_geno,alleleFreqs_results)
	
	all_combined<-merge(all_combined,hetContaminationResults)
	
	all_combined$Locus_ID<-as.character(all_combined$Locus_ID)
	#ensure appropriate columns are numeric
	all_combined$AvgReadDepth<-as.numeric(as.character(all_combined$AvgReadDepth))
	all_combined$GenotypeRate<-as.numeric(as.character(all_combined$GenotypeRate))
	all_combined$minFreq<-as.numeric(as.character(all_combined$minFreq))
	all_combined$maxFreq<-as.numeric(as.character(all_combined$maxFreq))
	#set final column names
	colnames(all_combined)<-c("Locus_ID","AvgReadDepth","GenotypeRate","minAF","majAF","alleles","allFreqs","conScore")
	
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
	#remove correction factor column from locusTable (if present)
	if("correctionFactors" %in% colnames(locusTable)){
		locusTable<-locusTable %>% select(-correctionFactors)
	}
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
require(msa)
library(tools)
#library(plyr)
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



########################################################################
#
#
#                               DEVELOPMENT
#
#
########################################################################

#summarize sample genotype rate and heterozygosity
summarizeSamples<-function(genotypes){
	#missing genotypes are coded as 0
	#calculate genotype rate as proportion of non-zero genotypes for each locus
	genotypeRate_results<-apply(genotypes,2,function(x){(length(x)-sum(x=="0"))/length(x)})
	genotypeRate_results<-data.frame(keyName=names(genotypeRate_results), value=genotypeRate_results,row.names=NULL)
	colnames(genotypeRate_results)<-c("sample","GenotypeRate")
	
	#calculate heterozygosity
	#function to count number of heterozygous genotypes
	#need to split genotypes by , then compare alleles to estimate heterozygosity
	countHets<-function(genos){
		alleles<-str_split(genos,",")
		hetCount<-sum(unlist(lapply(alleles,function(x) length(unique(x))))>1)
		return(hetCount)
	}
	heterozygosity_results<-apply(genotypes,2,function(x){countHets(x)/length(x)})
	heterozygosity_results<-data.frame(keyName=names(heterozygosity_results), value=heterozygosity_results,row.names=NULL)
	colnames(heterozygosity_results)<-c("sample","Heterozygosity")
	
	#calculate contamination score, this applies a binomial test to read counts from heterozygous genotypes
	#the contamination score is the proportion of heterozygous genotypes that failed the binomial test for each sample
	
	#function to create hetMatrix, 1 is heterozygous, 0 is homozygous
	idHets<-function(genos){
	  alleles<-str_split(genos,",")
	  hetCount<-unlist(lapply(alleles,function(x) length(unique(x))))
	  hetCount[hetCount!=2]<-0
	  hetCount[hetCount==2]<-1
	  return(hetCount)
	}

	#function to run binomial test
	binomTest<-function(reads){
	  A1<-as.numeric(str_split_fixed(reads,",",2)[,1])
	  A2<-as.numeric(str_split_fixed(reads,",",2)[,2])
	  tot<-A1+A2

	  if(tot!=0){
		result<-binom.test(A1, tot, 0.5, alternative=c("two.sided"), conf.level = 0.95)
		return(result$p.value)
	  }else{
		return(NA)
	  }
	}

	#function to manage binomial test
	testHet<-function(reads){
	  binomTest_results<-unlist(lapply(reads,binomTest))
	  numOutlier<-sum(binomTest_results<=0.05,na.rm=TRUE)
	  numHet<-sum(!is.na(binomTest_results))
	  #result is number of genotypes that failed binomial test over total number of heterozygous genotypes for each sample
	  result<-numOutlier/numHet
	  return(result)
	}

	#make matrix where heterozygous genotypes are 1 and homozygous genotypes are 0
	hetMatrix<-apply(polyGenResults_singleSNP,2,function(x){idHets(x)})

	#split reads into allele 1 and allele 2
	allele1<-apply(singleSNP_alleleReads,2,function(x) as.numeric(str_split_fixed(x,",",2)[,1]))
	allele2<-apply(singleSNP_alleleReads,2,function(x) as.numeric(str_split_fixed(x,",",2)[,2]))

	#multiply by het matrix to keep only counts from heterozygous genotypes
	allele1HetCounts<-allele1*hetMatrix
	allele2HetCounts<-allele2*hetMatrix

	#recombine counts so only heterozygous genotypes from original count data are retained
	combinedHetCounts<-paste(allele1HetCounts,allele2HetCounts,sep=",")
	#reset dimensions and reassign column names
	dim(combinedHetCounts)<-c(dim(allele1HetCounts)[1],dim(allele1HetCounts)[2])
	colnames(combinedHetCounts)<-colnames(allele1HetCounts)

	#run function to get heterozygous contamination score for all samples
	hetContaminationResults<-apply(combinedHetCounts,2,testHet)
	hetContaminationResults<-data.frame(sample=names(hetContaminationResults),conScore=hetContaminationResults,stringsAsFactors=FALSE,row.names=NULL)
	
	sampleSummary<-merge(genotypeRate_results,heterozygosity_results)
	sampleSummary<-merge(sampleSummary,hetContaminationResults)
	return(sampleSummary)
}

#generate ratio and scatter plots of genotypes for samples
plotGenotypes_sample<-function(locusTable,alleleReads,genotypes,type="ratio",savePlot="N",saveDir=""){
  #remove correction factor column from locusTable (if present)
  if("correctionFactors" %in% colnames(locusTable)){
    locusTable<-locusTable %>% select(-correctionFactors)
  }
  #combine data to allow running through apply
  #add a column specifying the number of columns for each section
  #genotypes, total allele reads, and allele ratios should all have the same number of columns
  sectionColumns<-dim(genotypes)[2]
  combinedData<-cbind(locusTable,sectionColumns,genotypes,alleleReads)
  #use invisible to suppress stdout from dev.off() while allowing errors to print
  
  #invisible(pbapply(combinedData,1,plotGenotypes_sample_setup,type,savePlot,saveDir))
  
  #send allele reads only for plotting
  #remove paralogs since these are expected to have heterozygous ratios that differ from 1:1
  alleleReads<-alleleReads[rownames(alleleReads) %in% locusTable$Locus_ID[locusTable$ploidy==2],]
  #formattedReads<-rbind(colnames(alleleReads),alleleReads)
  
  #filter genotypes and add below allele reads
  genotypes<-genotypes[row.names(genotypes) %in% locusTable$Locus_ID[locusTable$ploidy==2],]

  formattedReads<-rbind(colnames(alleleReads),alleleReads,genotypes)
  
  invisible(pbapply(formattedReads,2,plotGenotypes_sample_setup,type,savePlot,saveDir))
}

plotGenotypes_sample_setup<-function(genoData,type,savePlot,saveDir){
  ##separate data into categories
  #locusTable<-genoData[1:3]
  #sectionColumns<-as.numeric(genoData[4])
  #genotypesStart<-5
  #genotypesEnd<-5+sectionColumns-1
  #alleleReadsStart<-genotypesEnd+1
  #alleleReadsEnd<-alleleReadsStart+sectionColumns-1
  #genotypes<-genoData[genotypesStart:genotypesEnd]
  #alleleReads<-genoData[alleleReadsStart:alleleReadsEnd]
  
  #need to separate alleleReads from genotypes
  #allele reads are first half of data, genotypes are second
  alleleEnd<-(length(genoData)/2)+0.5
  
  alleleReads<-genoData[2:alleleEnd]

  
  genotypes<-genoData[(alleleEnd+1):length(genoData)]
  

  #convert genotypes to 0 (homozygous) or 1 (heterozygous)
  #split genotypes by allele
  genoAllele_1<-str_split_fixed(genotypes,",",2)[,1]
  genoAllele_2<-str_split_fixed(genotypes,",",2)[,2]
  hetGenos<-genoAllele_1==genoAllele_2
  hetGenos<-hetGenos*1
  
  #split reads by allele
  alleleReads_1<-str_split_fixed(alleleReads,",",2)[,1]
  alleleReads_2<-str_split_fixed(alleleReads,",",2)[,2]
  ##combine results into data frame
  #combinedData<-data.frame(genotypes,alleleReads_1,alleleReads_2,stringsAsFactors=FALSE)
  combinedData<-data.frame(alleleReads_1,alleleReads_2,stringsAsFactors=FALSE)
  combinedData$alleleReads_1<-as.numeric(combinedData$alleleReads_1)
  combinedData$alleleReads_2<-as.numeric(combinedData$alleleReads_2)
  #calculate total reads
  combinedData$totalReads<-combinedData$alleleReads_1+combinedData$alleleReads_2
  #calculate read ratio
  combinedData$ratios<-combinedData$alleleReads_1/combinedData$totalReads
  
  #add het genotypes to combined data
  combinedData$hetGenos<-hetGenos
  
  ##calculate summary statistics
  #genotypeRate<-(length(combinedData$genotypes)-sum(combinedData$genotypes=="0"))/length(combinedData$genotypes)
  #genotypeRate<-round(genotypeRate, digits=2)
  averageReadDepth<-mean(as.numeric(combinedData$totalReads))
  averageReadDepth<-round(averageReadDepth, digits=2)
  
  #plot results, exclude ratios with NA from plot to prevent warnings  
  sampleID<-genoData[1]
  #summaryText=paste("average depth: ",averageReadDepth,"					","genotype_rate: ", genotypeRate,sep=" ")
  #summaryText=paste("average depth: ",averageReadDepth,"					","genotype_rate: TODO", sep=" ")
  
  if(type=="ratio"){
    #genoPlot<-ggplot()+geom_histogram(data=combinedData,aes(ratios,color=genotypes, fill=genotypes),binwidth=0.01,alpha=0.5,na.rm=TRUE)+
    genoPlot<-ggplot()+geom_histogram(data=combinedData,aes(ratios),binwidth=0.01,alpha=0.5,na.rm=TRUE)+
      #xlim(-0.01,1.01)+labs(title=sampleID,x="Allele Ratio",y="Frequency",subtitle=summaryText)+
      xlim(-0.01,1.01)+labs(title=sampleID,x="Allele Ratio",y="Frequency")+
      #theme_bw()+theme(plot.title=element_text(hjust=0.5),plot.subtitle=element_text(hjust=0.5))
      theme_bw()+theme(plot.title=element_text(hjust=0.5))
  }else if(type=="scatter"){
    #genoPlot<-ggplot()+geom_point(data=combinedData,aes(x=alleleReads_1,y=alleleReads_2,color=genotypes))+
    genoPlot<-ggplot()+geom_point(data=combinedData,aes(x=alleleReads_1,y=alleleReads_2,color=as.factor(hetGenos)))+
      xlim(range(combinedData$alleleReads_1,combinedData$alleleReads_2))+
      ylim(range(combinedData$alleleReads_1,combinedData$alleleReads_2))+
      #labs(title=sampleID,x="Allele 1 Reads",y="Allele 2 Reads",subtitle=summaryText)+coord_fixed(ratio=1)+
      labs(title=sampleID,x="Allele 1 Reads",y="Allele 2 Reads",color="Het Genotype")+coord_fixed(ratio=1)+
      #theme_bw()+theme(plot.title=element_text(hjust=0.5),plot.subtitle=element_text(hjust=0.5))
      theme_bw()+theme(plot.title=element_text(hjust=0.5))
  }
  #save plot if specified
  if(savePlot=="Y"){
    JPGoutfile<-paste(saveDir,"/",sampleID,"_",type,".jpg",sep="")
    jpeg(file=JPGoutfile)
    print(genoPlot)
    #capture.output(dev.off(),file=NULL)
    dev.off()
    
  }else{
    print(genoPlot)
  }
}


################################################################################################
#
#                SEX GENOTYPING USING SEX-SPECIFIC LOCUS AND AUTOSOMAL CONTROL LOCUS
#
################################################################################################

genotypeSex<-function(locusReads,sexLocus,controlLocus,p_thresh=0.05,epsilon=0.01,plotGenotypes="N",saveDir=NULL){
  #save locus name as concatenated sexLocus and controlLocus
  sexLocus_genoID<-paste(sexLocus,controlLocus,sep="_")
  #make locus table
  sexLocusTable<-data.frame(Locus_ID=sexLocus_genoID,ploidy=2,alleles="X,Y")
  #for each locus, combine reads per allele into total reads per locus
  #format data
  targetLocusReads<-locusReads[rownames(locusReads) %in% c(sexLocus,controlLocus),]
  #get total reads for each locus
  targetLocusReads_combined<-apply(targetLocusReads,1:2,function(x) sum(as.numeric(str_split_fixed(x,",",2))))
  #combine reads for each locus into single locus
  sexLocus_alleleReads<-paste(targetLocusReads_combined[1,],targetLocusReads_combined[2,],sep=",")
  #save in alleleReads format
  sexLocus_alleleReads<-t(data.frame(sexLocus_alleleReads,stringsAsFactors=FALSE))
  colnames(sexLocus_alleleReads)<-colnames(locusReads)
  rownames(sexLocus_alleleReads)<-sexLocus_genoID
  #genotype with polygen
  sexLocus_genotypes<-polyGen(sexLocusTable,sexLocus_alleleReads,p_thresh,epsilon)
  
  #combine results into data framealleleReads<-genoData[2:length(genoData)]
  #split reads by allele
  alleleReads_1<-str_split_fixed(sexLocus_alleleReads,",",2)[,1]
  alleleReads_2<-str_split_fixed(sexLocus_alleleReads,",",2)[,2]
  #combine results into data frame
  genotypes<-t(sexLocus_genotypes)
  genotypes<-genotypes[,1]
  combinedData<-data.frame(samples=colnames(locusReads),alleleReads_1,alleleReads_2,genotypes,stringsAsFactors=FALSE)
  combinedData$alleleReads_1<-as.numeric(combinedData$alleleReads_1)
  combinedData$alleleReads_2<-as.numeric(combinedData$alleleReads_2)
  
  #plot genotypes
  if(plotGenotypes=="Y"){
   plotGenotypes(sexLocusTable, sexLocus_alleleReads, sexLocus_genotypes, type='scatter', savePlot="Y", saveDir=saveDir)
  }
  #return results
  #return(t(sexLocus_genotypes))
  return(combinedData)
}


#adjust read counts using correction factors
correctReads<-function(locusTable,readCounts){
  combinedData<-cbind(locusTable,readCounts)
  #pbapply has integrated progress bar, use apply if pbapply cannot be installed
  results<-pbapply(combinedData,1,readCorrection)
  
  results<-t(results)
  #add sample names back
  colnames(results)<-colnames(readCounts)
  return(results)
}



readCorrection<-function(readData){
  #get alleles for locus
  allele1<-str_split_fixed(readData[3],",",2)[,1]
  allele2<-str_split_fixed(readData[3],",",2)[,2]
  #get genotypes
  A1hom<-paste(allele1,allele1,sep=",")
  Het<-paste(allele1,allele2,sep=",")
  A2hom<-paste(allele2,allele2,sep=",")
  
  #get correction factors
  A1corr<-as.numeric(str_split_fixed(readData[4],",",2)[,1])
  A2corr<-as.numeric(str_split_fixed(readData[4],",",2)[,2])
  
  #subset combined data to allele reads
  reads<-readData[5:length(readData)]
  #get reads for allele 1, allele 2, and combined reads
  allele1Reads<-as.numeric(str_split_fixed(reads,",",2)[,1])
  allele2Reads<-as.numeric(str_split_fixed(reads,",",2)[,2])
  sum_xy=allele1Reads+allele2Reads
  
  #make dataframe to store data
  genoData<-data.frame(allele1Reads,allele2Reads,sum_xy,stringsAsFactors=FALSE)
  head(genoData)
  #adjust allele counts using correction factor
  #counts are coerced to integers to repliate behavior of campbell pipeline.
  #not that this is not the same as rounding.  
  #values of 3.2 and 3.8 will both give integers of 3
  genoData<-genoData %>% mutate(allele1Count = as.integer(allele1Reads - (sum_xy/4*A1corr))) %>%
    mutate(allele2Count = as.integer(allele2Reads - (sum_xy/4*A2corr)))
  
  head(genoData)
  
  #set negative reads to 0
  genoData$allele1Count[genoData$allele1Count<0]<-0
  genoData$allele2Count[genoData$allele2Count<0]<-0
  
  #combine reads and return
  newReads<-paste(genoData$allele1Count,genoData$allele2Count,sep=",")
  return(newReads)
}

#genotype using same method as Campbell GTseq genotyper
campbellStyle<-function(locusTable,readCounts){
  
  combinedData<-cbind(locusTable,readCounts)
  #pbapply has integrated progress bar, use apply if pbapply cannot be installed
  results<-pbapply(combinedData,1,campbellType)
  
  results<-t(results)
  #add sample names back
  colnames(results)<-colnames(readCounts)
  return(results)
  
}



campbellType<-function(readData){
  #get alleles for locus
  allele1<-str_split_fixed(readData[3],",",2)[,1]
  allele2<-str_split_fixed(readData[3],",",2)[,2]
  #get genotypes
  A1hom<-paste(allele1,allele1,sep=",")
  Het<-paste(allele1,allele2,sep=",")
  A2hom<-paste(allele2,allele2,sep=",")
  
  #get correction factors
  A1corr<-as.numeric(str_split_fixed(readData[4],",",2)[,1])
  A2corr<-as.numeric(str_split_fixed(readData[4],",",2)[,2])

  #subset combined data to allele reads
  reads<-readData[5:length(readData)]
  #get reads for allele 1, allele 2, and combined reads
  allele1Reads<-as.numeric(str_split_fixed(reads,",",2)[,1])
  allele2Reads<-as.numeric(str_split_fixed(reads,",",2)[,2])
  sum_xy=allele1Reads+allele2Reads
  
  #make dataframe to store data
  genoData<-data.frame(allele1Reads,allele2Reads,sum_xy,stringsAsFactors=FALSE)
  head(genoData)
  #adjust allele counts using correction factor
  #counts are coerced to integers to repliate behavior of campbell pipeline.
  #not that this is not the same as rounding.  
  #values of 3.2 and 3.8 will both give integers of 3
  genoData<-genoData %>% mutate(allele1Count = as.integer(allele1Reads - (sum_xy/4*A1corr))) %>%
                         mutate(allele2Count = as.integer(allele2Reads - (sum_xy/4*A2corr)))
  
  head(genoData)
  
  #set negative reads to 0
  genoData$allele1Count[genoData$allele1Count<0]<-0
  genoData$allele2Count[genoData$allele2Count<0]<-0
  
  head(genoData)
  

  #Fix allele counts to non-zero number for division ratio calculation...
  genoData<-genoData %>% mutate(A1fix=case_when(allele1Count==0 ~ 0.1,
                                                allele1Count!=0 ~ allele1Count),
                         A2fix=case_when(allele2Count==0 ~ 0.1,
                                                allele2Count!=0 ~ allele2Count),
                         ratio=A1fix/A2fix)
  
  
  head(genoData)
  
  #initialize genotypes as "00"
  genoData$geno="00"
  #use TRUE as catch-all statement at end of case_when to replace unmatched conditions with the original value "00
  genoData<-genoData %>% mutate(geno=case_when((allele1Count+allele2Count)<10 ~ "00",
                                               ratio >= 10 ~ A1hom,
                                               ratio <= 0.1 ~ A2hom,
                                               ratio <= 0.2 ~ "00",
                                               ratio <= 5 ~ Het,
                                               TRUE ~ .$geno))
  
  #initialize genotype classification as NA
  genoData$genoClass="NA"
  genoData<-genoData %>% mutate(genoClass=case_when((allele1Count+allele2Count)<10 ~ "NA",
                                                    ratio >= 10 ~ "A1HOM",
                                                    ratio <= 0.1 ~ "A2HOM",
                                                    ratio <= 0.2 ~ "in-betweeners",
                                                    ratio <= 5 ~ "HET"))
  #return(genoData$ratio)
  return(genoData$geno)
  #return(genoData$allele1Count)
  #return(genoData$allele2Count)
}



#function to compare genotypes between all pairs of samples
#input format: rows are loci, columns are samples, alleles are 0 and 1, missing genotypes are NA
#heterozygous calls must be consistent within a locus, ie all 0/1, no 1/0
IDduplicateSamples<-function(genotypes,MAF=NULL){
  #function to calculate MAF
  calcMAF<-function(locusGenos){
    allele1Counts<-sum(str_count(locusGenos,"0"),na.rm=TRUE)
    allele2Counts<-sum(str_count(locusGenos,"1"),na.rm=TRUE)
    allele1Freq<-allele1Counts/sum(allele1Counts,allele2Counts)
    if(allele1Freq>0.5){
      MAF<-1-allele1Freq
    }else{
      MAF<-allele1Freq
    }
    return(MAF)
  }
  
  #filter loci using MAF if threshold is specified
  if(!is.null(MAF)){
    #calculate MAF
    message(paste("MAF threshold applied:",MAF,"MAF",sep=" "))
    message("calculating MAF")
    locusMAF<-pbapply(genotypes,1,calcMAF)
    #convert to dataframe
    locusMAF<-data.frame(locus_ID=names(locusMAF),MAF=locusMAF,row.names=NULL)
    locusMAF$locus_ID<-as.character(locusMAF$locus_ID)
    #filter loc based on MAF threshold
    genotypes<-genotypes[rownames(genotypes)%in%locusMAF$locus_ID[locusMAF$MAF>=MAF],]
  }else{
    message("No MAF threshold applied, using all loci")
  }
  
  #make matrix of called vs NA genotypes for faster counting of missing data
  genotypes_NAmatrix<-genotypes
  genotypes_NAmatrix[!is.na(genotypes_NAmatrix)]<-0
  genotypes_NAmatrix[is.na(genotypes_NAmatrix)]<-1
  class(genotypes_NAmatrix)<-"numeric"
  
  #identify all unique pairs of samples
  allPairs<-combn(dim(genotypes)[2], 2)
  ncombo<-dim(allPairs)[2]
  nloci<-dim(genotypes)[1]
  message(paste("number of samples:",dim(genotypes)[2],sep=" "))
  message(paste("number of loci:",nloci,sep=" "))
  message(paste("number of sample pairs:",ncombo,sep=" "))
  #reshape into 2xn matrix
  allPairs<-matrix(allPairs,nrow=2)
  
  #function to compare genotypes
  compareGenos<-function(samplePair){
    #count number of loci that are genotyped in both samples
    NAcounts<-genotypes_NAmatrix[,samplePair[1]]+genotypes_NAmatrix[,samplePair[2]]
    sharedCounts<-nloci-(sum(NAcounts)-sum(NAcounts[NAcounts==2])/2)
    #count number of loci with genotypes that match in both samples
    genotypeMatches<-sum(genotypes[,samplePair[1]]==genotypes[,samplePair[2]],na.rm=TRUE)
    #return counts of genotype matches and shared loci
    return(c(genotypeMatches,sharedCounts))
  }
  
  #do all pairwise sample comparisons
  message("comparing genotypes")
  matches<-pbapply(allPairs,2,compareGenos)
  
  #make dataframe of results
  comparisonResults<-data.frame(matrix(NA,nrow=dim(allPairs)[2],ncol=7))
  colnames(comparisonResults)<-c("Sample1","Sample2","matchedGenotypes","commonGenotypes","proportionMatch","proportionCommon","totalLoci")
  comparisonResults$Sample1<-colnames(genotypes)[allPairs[1,]]
  comparisonResults$Sample2<-colnames(genotypes)[allPairs[2,]]
  comparisonResults$matchedGenotypes<-matches[1,]
  comparisonResults$commonGenotypes<-matches[2,]
  comparisonResults$proportionMatch<-comparisonResults$matchedGenotypes/comparisonResults$commonGenotypes
  comparisonResults$proportionCommon<-comparisonResults$commonGenotypes/nloci
  comparisonResults$totalLoci<-nloci
  return(comparisonResults)
}





