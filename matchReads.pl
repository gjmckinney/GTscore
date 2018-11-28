#/usr/bin/perl -w
use Algorithm::Combinatorics qw(combinations);
use strict;
use Getopt::Long;

#--------------------------------------------------------------------------------
#matchReads.pl
#
#Input:
#	Input consists of two files. The first file contains the primer/probe information
#	for each locus.  The second file contains a list of sequence files to count reads 
#	from.
#
#Usage:
#	perl matchReads.pl -p primerProbeFile.txt --files sampleList.txt
#
#Optional Arguments:
#	--matchType	type of sequence matching to perform, primer matching or primer/probe matching
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
#	7/9/2018
#--------------------------------------------------------------------------------

#get options from command line
my$primerProbeFile='';
my$sequenceFiles='';
my$prefix='';
my$matchType='primerProbe';

GetOptions('p=s' => \$primerProbeFile, 'files=s' => \$sequenceFiles, 'prefix=s' =>\$prefix, 'matchType=s' =>\$matchType);


#get minimum primer length
open(PRIMERPROBE,"<$primerProbeFile")||die "cannot open $primerProbeFile:$!";
my$primerProbeHeader=<PRIMERPROBE>;
my$minPrimerLength=100000000;
while (my$line=<PRIMERPROBE>){
	chomp $line;
	my($locusID,$ploidy,$SNPpos,$allele1,$allele2,$probe1,$probe2,$primer)=split "\t", $line;
	my$primerLength=length($primer);
	if($primerLength<$minPrimerLength){
		$minPrimerLength=$primerLength
	}
}


#set minimum depth to report unique sequence
my$minSeqDepth=1;

open(PRIMERPROBE,"<$primerProbeFile")||die "cannot open $primerProbeFile:$!";

my$primerProbeHeader=<PRIMERPROBE>;
chomp $primerProbeHeader;
my%primerProbeAllele;
my%primerProbeLocus;
my%primerProbeHapAllele;
my%loci;
my%hapLoci;
while (my$line=<PRIMERPROBE>){
	chomp $line;
	my($locusID,$ploidy,$SNPpos,$allele1,$allele2,$probe1,$probe2,$primer)=split "\t", $line;
	#store information in hash of hashes
	my$locus_SNP=$locusID."_".$SNPpos;
	my$trimmedPrimer=substr($primer, 0, $minPrimerLength);
	#store allele information
	$primerProbeAllele{$trimmedPrimer}{$probe1}=$allele1;
	$primerProbeAllele{$trimmedPrimer}{$probe2}=$allele2;
	#store SNP information
	$primerProbeLocus{$trimmedPrimer}{$probe1}=$locus_SNP;
	$primerProbeLocus{$trimmedPrimer}{$probe2}=$locus_SNP;
	
	#store information for haplotypes in format: primer,locus,snp position, probe, allele
	$primerProbeHapAllele{$trimmedPrimer}{$locusID}{$SNPpos}{$probe1}=$allele1;
	$primerProbeHapAllele{$trimmedPrimer}{$locusID}{$SNPpos}{$probe2}=$allele2;
	
	#store alleles for each locus
	#store single SNP allele information
	$loci{$locus_SNP}{$allele1}++;
	$loci{$locus_SNP}{$allele2}++;
	#store haplotype allele information
	$hapLoci{$locusID}{$SNPpos}{$allele1}++;
	$hapLoci{$locusID}{$SNPpos}{$allele2}++;
}
close PRIMERPROBE;

#store list of loci to order later output of locus table
my@loci;
foreach my$key (sort keys %loci){
	push@loci,$key;
}

#store list of loci to order later output of haplotype locus table
my@hapLoci;
foreach my$key (sort keys %hapLoci){
	push@hapLoci, $key;
}

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
	push@sequenceFiles,$line;
}
close SEQFILES;

my%alleleCounts;
my%haplotypeCounts;
my%matchedSeqs;
my@samples;
my$sampleCount=0;
foreach my$seqFile (@sequenceFiles){
	$sampleCount++;
	my($sampleID,$ext)=split '\.', $seqFile, 2;
	push@samples,$sampleID;
	open(SAMPLE,"<$seqFile")||die "cannot open $seqFile:$!";
	while(my$seqID=<SAMPLE>){
		chomp $seqID;
		my$sequence=<SAMPLE>;
		chomp $sequence;
		my$plus=<SAMPLE>;
		chomp $plus;
		my$qual=<SAMPLE>;
		chomp $qual;
		#get first 14 bases of sequence for matching to primer sequences (also trimmed to 14 bases)
		my$trimmedSeq=substr($sequence,0,$minPrimerLength);
		if(exists $primerProbeAllele{$trimmedSeq}){
			#get primr aligned matches
			if($matchType eq "primer"){
				foreach my$locusID (keys %{$primerProbeHapAllele{$trimmedSeq}}){
					$matchedSeqs{$sampleID}{$locusID}{$sequence}++;
				}
			#get primer-probe aligned matches
			}elsif($matchType eq "primerProbe"){
				#Count alleles for single SNP genotypes
				foreach my$probeKey (keys %{$primerProbeAllele{$trimmedSeq}}){
					#test for reverse complement matches to probe
					my$revCompProbeKey=reverse($probeKey);
					$revCompProbeKey=~tr/ATCG/TAGC/;
					$revCompProbeKey=~tr/\]\[/\[\]/;
					if(($sequence=~/$probeKey/)||($sequence=~/$revCompProbeKey/)){
						#store matched allele
						my$alleleMatch=$primerProbeAllele{$trimmedSeq}{$probeKey};
						#store matched locus
						my$locusMatch=$primerProbeLocus{$trimmedSeq}{$probeKey};
						#add sequence count to allele for matched locus
						$alleleCounts{$sampleID}{$locusMatch}{$alleleMatch}++;
						#add sequence to matched sequences
						$matchedSeqs{$sampleID}{$locusMatch}{$sequence}++;
					}
				}
				#Count alleles for multi-SNP haplotypes
				#Iterate through loci associated with primer, then iterate through SNPs/probes for each locus
				foreach my$IDkey (keys %{$primerProbeHapAllele{$trimmedSeq}}){
					my$haplotypeAlleleMatch;
					foreach my$SNPkey (sort {$a <=> $b} keys %{$primerProbeHapAllele{$trimmedSeq}{$IDkey}}){
						foreach my$probeKey (keys %{$primerProbeHapAllele{$trimmedSeq}{$IDkey}{$SNPkey}}){
							#test for reverse complement matches to probe
							my$revCompProbeKey=reverse($probeKey);
							$revCompProbeKey=~tr/ATCG/TAGC/;
							$revCompProbeKey=~tr/\]\[/\[\]/;
							if(($sequence=~/$probeKey/)||($sequence=~/$revCompProbeKey/)){
								#store matched single SNP allele
								my$alleleMatch=$primerProbeHapAllele{$trimmedSeq}{$IDkey}{$SNPkey}{$probeKey};
								#append matched single SNP allele to form haplotype 
								$haplotypeAlleleMatch.=$alleleMatch;
							}
						}
					}
					#Test if the length of the haplotype allele matches the expected haplotype length
					#If matched, add sequence count to haplotype allele for the matched locus
					#look into storing haplotypes that do not have correct length.
					#haplotypes with incorrect lengths are expected due to sequencing errors
					#high proportions of incorrect lengths could be due to systemic problems with probes
					if(length($haplotypeAlleleMatch) == $haplotypeLengths{$IDkey}){
						$haplotypeCounts{$sampleID}{$IDkey}{$haplotypeAlleleMatch}++;
					}
				}
			}
		}
		
	}
	close SAMPLE;
	print "sample $sampleCount finished: $sampleID\n";
}

#print output
my$outfile=$prefix."matchedReads_".$matchType."Aligned.txt";
#open(MATCHREADS,">ampliconGenotyper_matchedReads_primerProbeAligned.txt")||die "cannot open ampliconGenotyper_matchedReads_primerProbeAligned.txt:$!";
open(MATCHREADS,">$outfile")||die "cannot open $outfile:$!";
print MATCHREADS "Sample\tLocus\tSequence\tCount\n";
foreach my$sample(sort keys %matchedSeqs){
	foreach my$locus (sort keys %{$matchedSeqs{$sample}}){
		foreach my$uniqueSeq (keys %{$matchedSeqs{$sample}{$locus}}){
			if($matchedSeqs{$sample}{$locus}{$uniqueSeq}>=$minSeqDepth){
				print MATCHREADS"$sample\t$locus\t$uniqueSeq\t$matchedSeqs{$sample}{$locus}{$uniqueSeq}\n";
			}
		}
	}
}
close MATCHREADS;