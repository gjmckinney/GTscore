#/usr/bin/perl -w
use Algorithm::Combinatorics qw(combinations);
use strict;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#--------------------------------------------------------------------------------
#AmpliconReadCounter.pl
#
#Input:
#	Input consists of two files. The first file contains the primer/probe information
#	for each locus.  The second file contains a list of sequence files to count reads 
#	from.
#
#Usage:
#	perl AmpliconReadCounter.pl -p primerProbeFile.txt --files sampleList.txt
#
#Optional Arguments:
#	--printDiscarded	prints reads for each sample that did not have a primer and probe alignment for any locus
#	--useFullPrimer	uses the full primer for counting reads rather than the trimmed primer.
#	--prefix	optional prefix for output file names
#
#Input File Formats:
#	The primer/probe file is a tab delimited file containing the following information
#	for each locus: Locus ID, Ploidy, SNP position, Allele1, Allele2, Probe1, Probe2, and Primer.
#	If a locus has multiple SNPs an entry must be made for each SNP.
#
#	ex.
#	Locus	Ploidy	SNPpos	Allele1	Allele2	Probe1	Probe2	Primer
#	Locus1	4	53	C	G	CTGAGCGAAGG	CTGAGGGAAGG	AGGAAGTCTGAAGAGAGAACACTGA
#	Locus2	2	47	T	G	GAGGCTCTGAT	GAGGCGCTGAT	AGAGGAGCTGACCTGTGCTCT
#	Locus2	2	74	C	G	CAGAACAGCAT	CAGAAGAGCAT	AGAGGAGCTGACCTGTGCTCT
#
#Garrett McKinney
#University of Washington
#gjmckinn@uw.edu
#	5/9/2018
#--------------------------------------------------------------------------------

#record starting time so that total time to count reads can be reported
my$startTime=time;

#get options from command line
my$primerProbeFile='';
my$sequenceFiles='';
my$prefix='';
my$printDiscarded;
my$printMatched;
my$useFullPrimer;
my$alleleOrder="alphabetical";
my$inDir='';
my$outDir='';
my$inputType="fq";

GetOptions('p=s' => \$primerProbeFile, 'files=s' => \$sequenceFiles, 'prefix=s' =>\$prefix, 'printDiscarded' => \$printDiscarded, 
			'printMatched' => \$printMatched,'useFullPrimer' => \$useFullPrimer, 'alleleOrder=s' => \$alleleOrder, 'inDir=s' => \$inDir, 'outDir=s' => \$outDir,
			'i=s' => \$inputType);


#open primer probe file to get minimum primer length
#store this value for use in primer trimming below
open(PRIMERPROBE,"<$primerProbeFile")||die "cannot open $primerProbeFile:$!";
my$primerProbeHeader=<PRIMERPROBE>;


#initialize minimum primer length with ridiculously large value to guarantee real values are smaller
my$minPrimerLength=100000000;
while (my$line=<PRIMERPROBE>){
	chomp $line;
	my($locusID,$ploidy,$SNPpos,$allele1,$allele2,$probe1,$probe2,$primer,$A1corr,$A2corr)=split "\t", $line,10;
	my$primerLength=length($primer);
	if($primerLength<$minPrimerLength){
		$minPrimerLength=$primerLength
	}
}

#if($type eq "campbellStyle"){
#	$minPrimerLength=14;
#}

#$minPrimerLength=14;

open(PRIMERPROBE,"<$primerProbeFile")||die "cannot open $primerProbeFile:$!";

my$primerProbeHeader=<PRIMERPROBE>;
chomp $primerProbeHeader;

#get number of columns in primer probe file, use this to check if correction factors were used
my@primerProbeHeaderCols=split "\t", $primerProbeHeader;
my$corrFactorsEngaged=0;
if(scalar@primerProbeHeaderCols==10){
	$corrFactorsEngaged=1;
}

my%primerProbeAllele;
my%primerProbeLocus;
my%primerProbeHapAllele;
my%loci;
my%hapLoci;
my%ploidy;
my%hapPloidy;
#EXPERIMENTAL
my%primerLengths;
my%corrFactors;
my%primerProbeNum;
my%primerProbeAlleleNum;
my%lociNumOrder;
my%primerProbeHapLocus;
###experimental
my%SNPalleleOrder;

#store lists of loci to order later output of locus table and allele reads files
my@loci;
my@hapLoci;

while (my$line=<PRIMERPROBE>){
	chomp $line;
	my($locusID,$ploidy,$SNPpos,$allele1,$allele2,$probe1,$probe2,$primer,$A1corr,$A2corr)=split "\t", $line, 10;
	#store primer probe information in hash of hashes
	my$locus_SNP=$locusID."_".$SNPpos;
	#store ploidy for single SNP loci
	$ploidy{$locus_SNP}=$ploidy;
	#store ploidy for haplotype loci
	$hapPloidy{$locusID}=$ploidy;
	
	#store locus names to later order output of locus table and allele read files
	if(!exists $loci{$locus_SNP}){
		push@loci, $locus_SNP;
	}
	if(!exists $hapLoci{$locusID}){
		push@hapLoci, $locusID;
	}
	
	#if trimPrimers is enabled trim all primers to same length for substring sequence matching,
	#this greatly increase speed
	if(!$useFullPrimer){
		$primer=substr($primer, 0, $minPrimerLength);
	}
	#EXPERIMENTAL
	#my$trimmedPrimer=$primer;
	my$primerLength=length($primer);
	$primerLengths{$primerLength}++;
	
	#store allele order information for each SNP
	#if allele order is original then allele 1 is first, allele 2 is second
	#if allele order is alphabetical then sort alleles before storing
	if($alleleOrder eq "original"){
		$SNPalleleOrder{$locus_SNP}[1]=$allele1;
		$SNPalleleOrder{$locus_SNP}[2]=$allele2;
	}elsif($alleleOrder eq "alphabetical"){
		if($allele1 lt $allele2){
			$SNPalleleOrder{$locus_SNP}[1]=$allele1;
			$SNPalleleOrder{$locus_SNP}[2]=$allele2;
		}else{
			$SNPalleleOrder{$locus_SNP}[1]=$allele2;
			$SNPalleleOrder{$locus_SNP}[2]=$allele1;
		}
	}
	
	
	#store allele information
	$primerProbeAllele{$primer}{$probe1}=$allele1;
	$primerProbeAllele{$primer}{$probe2}=$allele2;
	#store SNP information
	$primerProbeLocus{$primer}{$probe1}=$locus_SNP;
	$primerProbeLocus{$primer}{$probe2}=$locus_SNP;
	
	#store information for haplotypes in format: primer,locus,snp position, probe, allele
	$primerProbeHapAllele{$primer}{$locusID}{$SNPpos}{$probe1}=$allele1;
	$primerProbeHapAllele{$primer}{$locusID}{$SNPpos}{$probe2}=$allele2;
	
	
	#store alleles for each locus
	#store single SNP allele information
	$loci{$locus_SNP}{$allele1}++;
	$loci{$locus_SNP}{$allele2}++;
	#store haplotype allele information
	$hapLoci{$locusID}{$SNPpos}{$allele1}++;
	$hapLoci{$locusID}{$SNPpos}{$allele2}++;
	
	
	######EXPERIMENTAL
	#store correction factor information
	$corrFactors{$locus_SNP}{$allele1}=$A1corr;
	$corrFactors{$locus_SNP}{$allele2}=$A2corr;
	#store probe information sorted by allele order
	$primerProbeNum{$primer}[1]=$probe1;
	$primerProbeNum{$primer}[2]=$probe2;
	#store allele information sorted by allele order
	$primerProbeAlleleNum{$primer}[1]=$allele1;
	$primerProbeAlleleNum{$primer}[2]=$allele2;
	#store allele information sorted by allele order
	$lociNumOrder{$locus_SNP}[1]=$allele1;
	$lociNumOrder{$locus_SNP}[2]=$allele2;
	#store haplotype locus name to access with primer in campbellStyle
	$primerProbeHapLocus{$primer}=$locusID;


}
close PRIMERPROBE;


#Store haplotype alleles for each locus.  Construct and save all possible combinations of SNPs ordered by position
my%lociHapAlleles;
my%haplotypeLengths;
foreach my$locus (keys %hapLoci){
	my%locusSNPposAlleles;
	my$SNPrank=0;
	#store SNP alleles based on ranked SNP position
	foreach my$SNPpos (sort {$a <=> $b} keys %{$hapLoci{$locus}}){
		$SNPrank++;
		my$alleleRank=0;
		foreach my$allele (sort keys %{$hapLoci{$locus}{$SNPpos}}){
			$alleleRank++;
			#print "$locus\t$SNPpos\t$allele\n";
			$locusSNPposAlleles{$SNPrank}{$alleleRank}=$allele;
		}
	}
	#store length of haplotypes for comparison against detected haplotypes in sequence data
	$haplotypeLengths{$locus}=$SNPrank;
	#Generate all possible haplotypes based on number of SNPs and alleles at each SNP.
	#The list of possible haplotypes will usually be larger than the number of observed haplotypes.
	my@alleleCodes;
	foreach my$i (1..$SNPrank){
		push@alleleCodes, (1,2);
	}
	my@all_combinations = combinations(\@alleleCodes, $SNPrank);
	#store unique combinations of alleles
	my%alleleCodeCombinations;
	foreach my$i (@all_combinations){
		my$alleleCombo;
		foreach my$j (@$i){
			$alleleCombo.=$j;
		}
		$alleleCodeCombinations{$alleleCombo}++;
	}
	#convert combinations of allele codes back to alleles and save as haplotypes
	foreach my$key (sort keys %alleleCodeCombinations){
		#print "$key\t";
		my$haplotype;
		my@SNPalleles=split "", $key;
		foreach my$i (0..$#SNPalleles){
			my$j=$i+1;
			$haplotype.=$locusSNPposAlleles{$j}{$SNPalleles[$i]};
		}
		$lociHapAlleles{$locus}{$haplotype}++;
		#print "$haplotype\n";
	}
}

#open list of sequence files and store in array
open(SEQFILES,"<$sequenceFiles")||die "cannot open $sequenceFiles:$!";
my@sequenceFiles;
while(my$line=<SEQFILES>){
	chomp $line;
	#skip empty lines
	if($line eq ''){
		next;
	}
	push@sequenceFiles,$line;
}
close SEQFILES;

my%alleleCounts;
my%haplotypeCounts;
my@samples;
my%sampleTotalReads;
my%sampleOffTarget;
my%samplePrimerOnlyMatched;
my%samplePrimerProbeMatched;
my%locusPrimerMatched;
my%locusPrimerProbeMatched;

my$sampleCount=0;
foreach my$seqFile (@sequenceFiles){
	my$sampleStartTime=time;
	$sampleCount++;
	my($sampleID,$ext)=split '\.', $seqFile, 2;
	push@samples,$sampleID;
	
	#if input directory is specified, add to seqFile
	if ($inDir){
		$seqFile=$inDir."/".$seqFile;
	}
	#Hash sequences for each sample, retain count of occurances for each unique sequence
	my%seqHash;
	my$sampleFH;
	if($inputType eq "fq"){
		open($sampleFH,"<$seqFile")||die "cannot open $seqFile:$!";
	}
	if($inputType eq "fastqgz"){
		$sampleFH = new IO::Uncompress::Gunzip $seqFile or die "gunzip failed: $GunzipError\n";
	}
	my$lineCount=0;
	while(my$line=<$sampleFH>){
		chomp $line;
		$lineCount++;
		if($lineCount%4==2){
			$seqHash{$line}++;
		}
	}
	close $sampleFH;
	
	#write unmatched reads to file if printDiscarded is enabled
	if ($printDiscarded){
		my$unmatchedReadFile=$sampleID."_unmatchedReads.txt";
		if ($outDir){
			$unmatchedReadFile=$outDir."/".$sampleID."_unmatchedReads.txt";
		}
		open(UNMATCHED,">$unmatchedReadFile")||die "cannot open $unmatchedReadFile:$!";
		print UNMATCHED "Sequence\tCount\n";
	}
	#write matched reads to file if printMatched is enabled
	if ($printMatched){
		my$matchedReadFile=$sampleID."_matchedReads.txt";
		if ($outDir){
			$matchedReadFile=$outDir."/".$sampleID."_matchedReads.txt";
		}
		open(MATCHED,">$matchedReadFile")||die "cannot open $matchedReadFile:$!";
		print MATCHED "Sequence\tCount\n";
	}
	
	my$matchedLocusID="";
	#iterate through hashed sequence and assign reads to loci based on primer/probe matches
	foreach my$seq (keys %seqHash){
		my$sequence=$seq;
		my$uniqueSeqCount=$seqHash{$seq};
		#add counts to total sequence counts for individual
		$sampleTotalReads{$sampleID}+=$uniqueSeqCount;
		#set flags for counting reads as matched or unmatched, use to count off-target, primer matched, and primer probe matched reads
		my$primerMatched=0;
		my$primerProbeMatched=0;
		#get first 14 bases of sequence for matching to primer sequences (also trimmed to 14 bases)
		#my$trimmedSeq=substr($sequence,0,14);
		foreach my$primerLength(sort {$b <=> $a} keys %primerLengths){
			
			if($primerProbeMatched==1){
				next;
			}
			my$trimmedSeq=substr($sequence,0,$primerLength);
			
			if(exists $primerProbeAllele{$trimmedSeq}){
				#set flag for primer matched reads
				$primerMatched=1;
				#Count alleles for single SNP genotypes
				foreach my$probeKey (keys %{$primerProbeAllele{$trimmedSeq}}){
					#fix indels first by changing to numeric code, then back with the ? in front
					my$revCompProbeKey=$probeKey;
					$revCompProbeKey=~s/A\?/1/g;
					$revCompProbeKey=~s/C\?/2/g;
					$revCompProbeKey=~s/G\?/3/g;
					$revCompProbeKey=~s/T\?/4/g;
					$revCompProbeKey=~s/1/\?A/g;
					$revCompProbeKey=~s/2/\?C/g;
					$revCompProbeKey=~s/3/\?G/g;
					$revCompProbeKey=~s/4/\?T/g;
					#test for reverse complement matches to probe
					#my$revCompProbeKey=reverse($probeKey);
					$revCompProbeKey=reverse($revCompProbeKey);
					$revCompProbeKey=~tr/ATCG/TAGC/;
					$revCompProbeKey=~tr/\]\[/\[\]/;
					
					if(($sequence=~/$probeKey/)||($sequence=~/$revCompProbeKey/)){
						#store matched allele
						my$alleleMatch=$primerProbeAllele{$trimmedSeq}{$probeKey};
						#store matched locus
						my$locusMatch=$primerProbeLocus{$trimmedSeq}{$probeKey};
						#add sequence count to allele for matched locus
						$alleleCounts{$sampleID}{$locusMatch}{$alleleMatch}+=$uniqueSeqCount;
						#set flag for primer probe matched reads
						$primerProbeMatched=1;
					}
				}
				#Count alleles for multi-SNP haplotypes
				#Iterate through loci associated with primer, then iterate through SNPs/probes for each locus
				foreach my$IDkey (keys %{$primerProbeHapAllele{$trimmedSeq}}){
					#EXPERIMENTAL
					$matchedLocusID=$IDkey;
					#EXPERIMENTAL
					#count number of reads that matched locus primer
					#$locusPrimerMatched{$IDkey}+=$uniqueSeqCount;
					my$haplotypeAlleleMatch;
					foreach my$SNPkey (sort {$a <=> $b} keys %{$primerProbeHapAllele{$trimmedSeq}{$IDkey}}){
						foreach my$probeKey (keys %{$primerProbeHapAllele{$trimmedSeq}{$IDkey}{$SNPkey}}){
							#fix indels first by changing to numeric code, then back with the ? in front
							my$revCompProbeKey=$probeKey;
							$revCompProbeKey=~s/A\?/1/g;
							$revCompProbeKey=~s/C\?/2/g;
							$revCompProbeKey=~s/G\?/3/g;
							$revCompProbeKey=~s/T\?/4/g;
							$revCompProbeKey=~s/1/\?A/g;
							$revCompProbeKey=~s/2/\?C/g;
							$revCompProbeKey=~s/3/\?G/g;
							$revCompProbeKey=~s/4/\?T/g;
							#test for reverse complement matches to probe
							#my$revCompProbeKey=reverse($probeKey);
							$revCompProbeKey=reverse($revCompProbeKey);
							$revCompProbeKey=~tr/ATCG/TAGC/;
							$revCompProbeKey=~tr/\]\[/\[\]/;
							if(($sequence=~/$probeKey/)||($sequence=~/$revCompProbeKey/)){
								#store matched single SNP allele
								my$alleleMatch=$primerProbeHapAllele{$trimmedSeq}{$IDkey}{$SNPkey}{$probeKey};
								#append matched single SNP allele to form haplotype 
								$haplotypeAlleleMatch.=$alleleMatch;
								#set flag for primer probe matched reads
								$primerProbeMatched=1;
							}
						}
					}
					
					#Test if the length of the haplotype allele matches the expected haplotype length
					#If matched, add sequence count to haplotype allele for the matched locus
					#look into storing haplotypes that do not have correct length.
					#haplotypes with incorrect lengths are expected due to sequencing errors
					#high proportions of incorrect lengths could be due to systemic problems with probes
					if(length($haplotypeAlleleMatch) == $haplotypeLengths{$IDkey}){
						$haplotypeCounts{$sampleID}{$IDkey}{$haplotypeAlleleMatch}+=$uniqueSeqCount;
					}
				}
			}
		}
		#for each sample assign read counts to primer probe matched, probe only matched, or off target reads
		if($primerProbeMatched==1){
			$samplePrimerProbeMatched{$sampleID}+=$uniqueSeqCount;
			if ($printMatched){
				print MATCHED "$sequence\t$seqHash{$seq}\n";
			}
		}elsif($primerMatched==1){
			$samplePrimerOnlyMatched{$sampleID}+=$uniqueSeqCount;
			if ($printDiscarded){
				print UNMATCHED "$sequence\t$seqHash{$seq}\n";
			}
		}else{
			$sampleOffTarget{$sampleID}+=$uniqueSeqCount;
			if ($printDiscarded){
				print UNMATCHED "$sequence\t$seqHash{$seq}\n";
			}
		}
		if($primerMatched==1){
			$locusPrimerMatched{$matchedLocusID}+=$uniqueSeqCount;
		}
		if($primerProbeMatched==1){
			$locusPrimerProbeMatched{$matchedLocusID}+=$uniqueSeqCount;
		}
	}
	if ($printDiscarded){
		close UNMATCHED;
	}
	if ($printMatched){
		close MATCHED;
	}
	my$sampleEndTime=time;
	my$sampleTotalTime=$sampleEndTime-$sampleStartTime;
	print "sample $sampleCount finished: $sampleID\tProcessing time: $sampleTotalTime seconds\n";
}

#################################################################################
#
#                      CREATE AND PRINT INDIVIDUAL SUMMARIES
#
#################################################################################
my$indSummaryFile=$prefix."GTscore_individualSummary.txt";
if ($outDir){
	$indSummaryFile=$outDir."/".$prefix."GTscore_individualSummary.txt";
}
open(INDSUMMARY,">$indSummaryFile")||die "cannot open $indSummaryFile:$!";
#open(INDSUMMARY,">GTscore_individualSummary.txt")||die "cannot open GTscore_individualSummary.txt:$!";
print INDSUMMARY "Sample\tTotal Reads\tOff-target Reads\tPrimer Only Reads\tPrimer Probe Reads\tOff-target Proportion\tPrimer Only Proportion\tPrimer Probe Proportion\n";
foreach my$sample (@samples){
	my$offTargetPercent;
	my$primerOnlyPercent;
	my$primerProbePercent;
	#if no reads exist for a sample set total reads to 0 and all other entries to NA
	#If reads exist for a sample, check each category and set reads to 0 no reads were recorded for that category
	if(not exists $sampleTotalReads{$sample}){
		$sampleTotalReads{$sample}=0;
		$offTargetPercent="NA";
		$primerOnlyPercent="NA";
		$primerProbePercent="NA";
		$sampleOffTarget{$sample}="NA";
		$samplePrimerOnlyMatched{$sample}="NA";
		$samplePrimerProbeMatched{$sample}="NA";
	}else{
		if(not exists $sampleOffTarget{$sample}){
			$sampleOffTarget{$sample}=0;
		}
		if(not exists $samplePrimerOnlyMatched{$sample}){
			$samplePrimerOnlyMatched{$sample}=0;
		}
		if(not exists $samplePrimerProbeMatched{$sample}){
			$samplePrimerProbeMatched{$sample}=0;
		}
		$offTargetPercent=$sampleOffTarget{$sample}/$sampleTotalReads{$sample};
		$offTargetPercent=sprintf("%.2f",$offTargetPercent);
		$primerOnlyPercent=$samplePrimerOnlyMatched{$sample}/$sampleTotalReads{$sample};
		$primerOnlyPercent=sprintf("%.2f",$primerOnlyPercent);
		$primerProbePercent=$samplePrimerProbeMatched{$sample}/$sampleTotalReads{$sample};
		$primerProbePercent=sprintf("%.2f",$primerProbePercent);
	}
	print INDSUMMARY "$sample\t$sampleTotalReads{$sample}\t$sampleOffTarget{$sample}\t$samplePrimerOnlyMatched{$sample}\t$samplePrimerProbeMatched{$sample}\t$offTargetPercent\t$primerOnlyPercent\t$primerProbePercent\n";
}
close INDSUMMARY;

#################################################################################
#
#                        CREATE AND PRINT LOCUS SUMMARIES
#
#################################################################################

my$locusSummaryFile=$prefix."GTscore_locusSummary.txt";
if ($outDir){
	$locusSummaryFile=$outDir."/".$prefix."GTscore_locusSummary.txt";
}
open(LOCUSSUMMARY,">$locusSummaryFile")||die "cannot open $locusSummaryFile:$!";
#open(LOCUSSUMMARY,">GTscore_locusSummary.txt")||die "cannot open GTscore_locusSummary.txt:$!";
print LOCUSSUMMARY "Locus\tPrimer Reads\tPrimer Probe Reads\tPrimer Probe Proportion\n";
foreach my$locus (@hapLoci){
	#set reads to 0 if read don't exist for a category
	if(not exists $locusPrimerMatched{$locus}){
		$locusPrimerMatched{$locus}=0;
	}
	if(not exists $locusPrimerProbeMatched{$locus}){
		$locusPrimerProbeMatched{$locus}=0;
	}
	#print LOCUSSUMMARY "$locus\t$locusPrimerMatched{$locus}\t$locusPrimerProbeMatched{$locus}\t";
	my$percPrimerProbe;
	if($locusPrimerMatched{$locus}==0){
		$percPrimerProbe="NA";
	}else{
		$percPrimerProbe=$locusPrimerProbeMatched{$locus}/$locusPrimerMatched{$locus};
		$percPrimerProbe=sprintf("%.2f",$percPrimerProbe);
	}
	#print LOCUSSUMMARY "$percPrimerProbe\n";
	print LOCUSSUMMARY "$locus\t$locusPrimerMatched{$locus}\t$locusPrimerProbeMatched{$locus}\t$percPrimerProbe\n";
}
close LOCUSSUMMARY;


#################################################################################
#
#                         CREATE AND PRINT LOCUS TABLES
#
#################################################################################

#print locus table for single-SNP genotypes
my$singleSNPlocusTable=$prefix."LocusTable_singleSNPs.txt";
if ($outDir){
	$singleSNPlocusTable=$outDir."/".$prefix."LocusTable_singleSNPs.txt";
}
open(LOCUSTABLE,">$singleSNPlocusTable")||die "cannot open $singleSNPlocusTable:$!";
#open(LOCUSTABLE,">LocusTable_singleSNPs.txt")||die "cannot open LocusTable_singleSNPs.txt:$!";
#print LOCUSTABLE "Locus_ID\trefAllele\taltAlleles\talleles\n";
if($corrFactorsEngaged==1){
	print LOCUSTABLE "Locus_ID\tploidy\talleles\tcorrectionFactors\n";
}else{
	print LOCUSTABLE "Locus_ID\tploidy\talleles\n";
}

foreach my$locus (@loci){
	#print LOCUSTABLE "$locus";
	print LOCUSTABLE "$locus\t$ploidy{$locus}";
	my$alleles;
	my$alleleCount=0;
	my$corrFactors;
	
	my$allele1=$SNPalleleOrder{$locus}[1];
	my$allele2=$SNPalleleOrder{$locus}[2];
	
	$alleles=$allele1.",".$allele2;
	$corrFactors=$corrFactors{$locus}{$allele1}.",".$corrFactors{$locus}{$allele2};
	if($corrFactorsEngaged==1){
		print LOCUSTABLE "\t$alleles\t$corrFactors\n";
	}else{
		print LOCUSTABLE "\t$alleles\n";
	}
}


#if($type eq "GTscore"){
	#print locus table for multi-SNP haplotypes
	my$haplotypeLocusTable=$prefix."LocusTable_haplotypes.txt";
	if ($outDir){
		$haplotypeLocusTable=$outDir."/".$prefix."LocusTable_haplotypes.txt";
	}
	open(HAPLOCUSTABLE,">$haplotypeLocusTable")||die "cannot open $haplotypeLocusTable:$!";
	#open(HAPLOCUSTABLE,">LocusTable_haplotypes.txt")||die "cannot open LocusTable_haplotypes.txt:$!";
	print HAPLOCUSTABLE "Locus_ID\tploidy\talleles\n";
	foreach my$locus (@hapLoci){
		#print HAPLOCUSTABLE "$locus";
		print HAPLOCUSTABLE "$locus\t$hapPloidy{$locus}";
		my$alleles;
		my$alleleCount=0;
		foreach my$alleleKey (sort keys %{$lociHapAlleles{$locus}}){
			$alleleCount++;
			#print HAPLOCUSTABLE "\t$alleleKey";
			if($alleleCount==1){
				$alleles=$alleleKey;
			}else{
				$alleles.=",$alleleKey";
			}
		}
		print HAPLOCUSTABLE "\t$alleles\n";
	}
#}

#################################################################################
#
#                     CREATE AND PRINT READ COUNT TABLES
#
#################################################################################

#print alleleReads table for single-SNP genotypes
my$alleleReadsFile=$prefix."AlleleReads_singleSNPs.txt";
if ($outDir){
	$alleleReadsFile=$outDir."/".$prefix."AlleleReads_singleSNPs.txt";
}
open(ALLELEREADS,">$alleleReadsFile")||die "cannot open $alleleReadsFile:$!";
#open(ALLELEREADS,">AlleleReads_singleSNPs.txt")||die "cannot open AlleleReads_singleSNPs.txt:$!";
foreach my$sample (@samples){
	print ALLELEREADS "\t$sample"
}
print ALLELEREADS"\n";

foreach my$locus (@loci){
	print ALLELEREADS "$locus";
	foreach my$sample (@samples){
		#print "$sample\n";
		print ALLELEREADS "\t";
		if(exists $alleleCounts{$sample}{$locus}){
			#my$allele1=$lociNumOrder{$locus}[1];
			#my$allele2=$lociNumOrder{$locus}[2];
			my$allele1=$SNPalleleOrder{$locus}[1];
			my$allele2=$SNPalleleOrder{$locus}[2];
			
			if(exists $alleleCounts{$sample}{$locus}{$allele1}){
				print ALLELEREADS "$alleleCounts{$sample}{$locus}{$allele1}";
			}else{
				print ALLELEREADS "0";
			}
			if(exists $alleleCounts{$sample}{$locus}{$allele2}){
				print ALLELEREADS ",$alleleCounts{$sample}{$locus}{$allele2}";
			}else{
				print ALLELEREADS ",0";
			}
		}else{
			print ALLELEREADS "0,0";
		}
	}
	print ALLELEREADS "\n";
}


#print alleleReads table for multi-SNP haplotypes
my$haplotypeAlleleReadsFile=$prefix."AlleleReads_haplotypes.txt";
if ($outDir){
	$haplotypeAlleleReadsFile=$outDir."/".$prefix."AlleleReads_haplotypes.txt";
}
open(HAPALLELEREADS,">$haplotypeAlleleReadsFile")||die "cannot open $haplotypeAlleleReadsFile:$!";
#open(HAPALLELEREADS,">AlleleReads_haplotypes.txt")||die "cannot open AlleleReads_haplotypes.txt:$!";
foreach my$sample (@samples){
	print HAPALLELEREADS "\t$sample"
}
print HAPALLELEREADS"\n";

foreach my$locus (@hapLoci){
	print HAPALLELEREADS "$locus";
	foreach my$sample (@samples){
		#print "$sample\n";
		print HAPALLELEREADS "\t";
		my$alleleCount=0;
		#$haplotypeCounts{$sampleID}{$IDkey}{$haplotypeAlleleMatch}++;
		if(exists $haplotypeCounts{$sample}{$locus}){
			#$lociHapAlleles{$locus}{$haplotype}++;
			foreach my$alleleKey (sort keys %{$lociHapAlleles{$locus}}){
				$alleleCount++;
				if(exists $haplotypeCounts{$sample}{$locus}{$alleleKey}){
					#print the number of reads found for an allele
					#print "$sample\t$locus\t$alleleKey\t$alleleCounts{$sample}{$locus}{$alleleKey}\n";
					if($alleleCount==1){
						print HAPALLELEREADS "$haplotypeCounts{$sample}{$locus}{$alleleKey}";
					}else{
						print HAPALLELEREADS ",$haplotypeCounts{$sample}{$locus}{$alleleKey}";
					}
				}else{
					#if there are no reads for an allele print 0
					if($alleleCount==1){
						print HAPALLELEREADS "0";
					}else{
						print HAPALLELEREADS ",0";
					}
				}
			}
		}else{
			#if a locus has no reads for a sample print 0 for each allele
			foreach my$alleleKey (sort keys %{$lociHapAlleles{$locus}}){
				$alleleCount++;
				#print "$sample\t$locus\t$alleleKey\t0\n";
				if($alleleCount==1){
					print HAPALLELEREADS "0";
				}else{
					print HAPALLELEREADS ",0";
				}
			}
		}
	}
	print HAPALLELEREADS "\n";
}


my$endTime=time;
my$totalTime=$endTime-$startTime;
print "Total Time: $totalTime\n";