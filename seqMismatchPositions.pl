#/usr/bin/perl -w
use strict;
use Getopt::Long;

#--------------------------------------------------------------------------------
#seqMismatchPositions.pl
#
#Input:
#	Input consists of two files.  The first file contains reference sequence for each
#	amplicon in a tab delimited format.  The second file contains the matched sequences 
#	output by matchReaqds.pl.
#	
#
#Usage:
#	perl seqMismatchPositions.pl --amplicons ampliconRef.txt --matchedSeqs matchedSeqs.txt
#
#Optional Arguments:
#	--matchType	type of sequence matching to perform, primer matching or primer/probe matching
#	--prefix	optional prefix for output file names
#
#Input File Formats:
#	The amplicon reference sequence is a tab delimited file with two columns.  The first column
#	contains the locus name and the second column contains the reference sequence for the locus
#	amplicon.
#	
#	ex.
#	Locus	refSeq
#	Ots_RAD5584	TGCTACCCTCTGTTTACTGTTGTGCTATCATCAGATAATAGCKTCTTATGCTTTCTTTCSCCGAAAAGCCTTTTTAAAATCTGACATGTTGGCTGGATTC
#	Ots_RAD16850	CCTGCAGGTCAGGGTTTGCGGCAGCCTTCTGCAGYTCAGCGTTGATAATCTTTTCGGTCTTCTCCCG
#
#	The matchedSeqs files is a tab delimited file output by matchSeqs.pl; this file has four columns.
#	The first column contains the sample that the read was from, the second column contains the locus
#	that the sequence aligned to.  The third column contains the sequence.  The fourth column contains
#	the number of times that unique sequence was observed in the sequence file.
#
#	ex.
#	Sample1	Locus1	ATTGCACTCCCCACACCTGGCCGGGGGGCAGGGCGCTGTGGAGGAGGAGGGGTCAGAGTAGGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAA	51
#	Sample1	Locus2	GCAGGTAGATGGATTGGTTTGACCTATTTTCTATTAATAAAAAACTAAGTGTAAAAACTACACAGTGAAACACTGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA	310
#	Sample1	Locus2	GCAGGTAGATGGATTGGTTTGACTTATTTTCTATTAATAAAAAACTAAGTGTAAAAACTACACAGTGAAACACTGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA	300
#
#
#Garrett McKinney
#University of Washington
#gjmckinn@uw.edu
#	7/9/2018
#--------------------------------------------------------------------------------

my$ampliconFile='';
my$matchedSeqFile='';
my$prefix='';
my$matchType='primerProbe';

GetOptions('amplicons=s' => \$ampliconFile, 'matchedSeqs=s' => \$matchedSeqFile, 'prefix=s' =>\$prefix, 'matchType=s' =>\$matchType);

open(AMPLICON,"<$ampliconFile")||die "cannot open $ampliconFile:$!";
my$seqID;
my%amplicons;
while(my $line=<AMPLICON>){
	chomp $line;
	if($line=~/^>/){
		($seqID,my$discard)=split " ", $line, 2;
		$seqID=~s/^>//;
	}else{
		$amplicons{$seqID}=$line;
	}
}
close AMPLICON;

open(MATCHEDSEQS,"<$matchedSeqFile")||die "cannot open $matchedSeqFile:$!";
my%locusMismatches;
while(my$line=<MATCHEDSEQS>){
	chomp $line;
	my($sample,$locus,$seq,$seqCount)=split "\t", $line;
	my($tag,$SNP);
	my@ampliconBases;
	
	if($matchType eq "primer"){
		@ampliconBases=split "", $amplicons{$locus};
	}elsif($matchType eq "primerProbe"){
		($tag,$SNP)=split(/_([^_]+)$/, $locus);
		@ampliconBases=split "", $amplicons{$tag};
	}
	my@seqBases=split "", $seq;
	
	if(exists $locusMismatches{$locus}){
		foreach my$i (0..$#ampliconBases){
			if($ampliconBases[$i] ne $seqBases[$i]){
				$locusMismatches{$locus}[$i]+=$seqCount;
			}
		}
	}else{
		foreach my$i (0..$#ampliconBases){
			$locusMismatches{$locus}[$i]=0;
		}
	}
}
close MATCHEDSEQS;

#print output
my$outfile=$prefix."mismatchPositions_".$matchType.".txt";
open(MISMATCHPOSITIONS,">$outfile")||die "cannot open $outfile:$!";

print MISMATCHPOSITIONS "Locus\tRefSeq\tMismatches\n";
foreach my$locus (sort keys %locusMismatches){
	my($tag,$SNP);
	
	if($matchType eq "primer"){
		print MISMATCHPOSITIONS "$locus\t$amplicons{$locus}\t";
	}elsif($matchType eq "primerProbe"){
		($tag,$SNP)=split(/_([^_]+)$/, $locus);
		print MISMATCHPOSITIONS "$locus\t$amplicons{$tag}\t";
	}
	
	my$bases;
	print MISMATCHPOSITIONS "$locusMismatches{$locus}[0]";
	foreach my$position (1..$#{$locusMismatches{$locus}}){
		#$bases.=$locusMismatches{$locus}[$position];
		print MISMATCHPOSITIONS ",$locusMismatches{$locus}[$position]";
	}
	print MISMATCHPOSITIONS "\n";
}