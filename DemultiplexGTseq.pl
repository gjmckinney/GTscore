#-------------
# Begin script
#-------------

#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#--------------------------------------------------------------------------------
#Demultiplex.pl
#
#Input:
#	Input consists of two files. The first file contains barcodes and sample IDs, 
#	the second file contains the raw sequence data.
#
#Usage:
#	perl Demultiplex.pl -b barcodeFile.txt -s sequencefile.fastq
#
#Input File Formats:	
#	Barcode file should include each dual index barcode and sample ID in a separate line
#	The sample name is in the first column, the i7 barcode in the second column, and the i5 barcode in the third column
#	samples can be named whatever you like as long as there are no spaces in the sample name
#
#		ex. sample1	ATCACGAT	AAACGGTC
#			sample2	ATCACGAT	ACTCTTTC
#			sample3	ATCACGAT	ATTCCGTC
#
#	Sequence file should be a standard fastq format which consist of lines in 
#	repeating blocks of four.  
#		The first line includes ID information
#		The second line includes sequence information
#			(this is where barcodes are matched)
#		The third line contains more ID information?
#		The fourth line contains the quality score
#
#Output file Format:
#	A file containing demultiplexed sequences in fastq format is automatically
#	output for each sample using the sample names specified in the barcode file
#
#Garrett McKinney
#University of Washington
#gjmckinn@uw.edu
#	5/23/2018
#--------------------------------------------------------------------------------

my$barcodeFile='';
my$sequenceFile='';
GetOptions('b=s' => \$barcodeFile, 's=s' => \$sequenceFile);


#create empty file for discarded reads
my($seqFileName,$ext)=split '\.', $sequenceFile, 2;
my$discardFile=$seqFileName."_discarded.fastq";
open(DISCARDED, ">$discardFile")||die "cannot open $discardFile:$!";

#open barcode file and store information in hash
open(BARCODES,"<$barcodeFile")||die "cannot open $barcodeFile:$!";
my%barcodes;
my$barcodeLength;
#strip header from file
my$barcodeHeader=<BARCODES>;
while(my$line=<BARCODES>){
	chomp $line;
	#my($BC_i7,$BC_i5,$ID)=split "\t", $line, 3;
	my($ID,$BC_i7,$BC_i5)=split "\t", $line, 3;
	my$barcode=$BC_i7."+".$BC_i5;
	$barcodes{$barcode}=$ID;
	$barcodeLength=length($barcode);
}
close BARCODES;

#create sequence file for demultiplexing
#these will be written to when parsing the sequence data
foreach my$barcode (sort keys %barcodes){
	my$outFile=$barcodes{$barcode}.".fastq";
	open(OUTFILE,">$outFile")||die "cannot open $outFile:$!";
	close OUTFILE;
}

open(SEQUENCE,"<$sequenceFile")||die "cannot open $sequenceFile:$!";
my%barcode_seqs;
my$interval=0;
my$seqCount=0;
my$retainedReads=0;
my$discardedReads=0;
my%retainedBarcodeHash;
my%discardedBarcodeHash;
my$startTime=time;
my$lastTime=time;
while(my$seqID=<SEQUENCE>){
	my$sequence=<SEQUENCE>;
	my$plus=<SEQUENCE>;
	my$qual=<SEQUENCE>;
	chomp $seqID;
	chomp $sequence;
	chomp $plus;
	chomp $qual;
	
	#combinatorial barcode should be stored as last element of ID line, separated by colons 
	my@IDinfo=split ":", $seqID;
	my$seqBC=$IDinfo[9];
	
	#store sequence in barcode_seqs hash
	$barcode_seqs{$seqBC}.=$seqID."\n".$sequence."\n".$plus."\n".$qual."\n";
	
	#store 1,000,000 sequences, then print results and clear hash
	#last set of sequences won't print here, print after loop
	$seqCount++;
	$interval++;
	if($interval==1000000){
		#print sequences that matched sample barcode
		foreach my$barcode (sort keys %barcodes){
			my$outFile=$barcodes{$barcode}.".fastq";
			open(OUTFILE,">>$outFile")||die "cannot open $outFile:$!";
			print OUTFILE $barcode_seqs{$barcode};
			close OUTFILE;
			#count reads for retained sequences
			#reads have been concatenated so can't be counted directly
			#need to count number of newlines and divide by four
			my$counts=$barcode_seqs{$barcode}=~tr/\n//;
			$retainedReads+=$counts/4;
			#add counts of retained barcode to hash
			$retainedBarcodeHash{$barcode}+=$counts/4;
			#delete sequence from hash
			delete $barcode_seqs{$barcode};
		}
		#print sequences that did not match sample barcode
		foreach my$unmatched (sort keys %barcode_seqs){
			print DISCARDED $barcode_seqs{$unmatched};
			#count reads for retained sequences
			#reads have been concatenated so can't be counted directly
			#need to count number of newlines and divide by four
			my$counts=$barcode_seqs{$unmatched}=~tr/\n//;
			$discardedReads+=$counts/4;
			#add counts of discarded barcode to hash
			$discardedBarcodeHash{$unmatched}+=$counts/4;
		}
		
		#reset sequence counter and storage hash
		$interval=0;
		%barcode_seqs=();
		my$duration=time-$lastTime;
		$lastTime=time;
		print "Total Reads Processed: $seqCount\nInterval Time: $duration\n";
	}
	##printing progress for each read doubles run time
	#print "Reads Processed: $seqCount\r";
}
close SEQUENCE;


#print last sequences that matched sample barcode
foreach my$barcode (sort keys %barcodes){
	my$outFile=$barcodes{$barcode}.".fastq";
	open(OUTFILE,">>$outFile")||die "cannot open $outFile:$!";
	print OUTFILE $barcode_seqs{$barcode};
	close OUTFILE;
	#count retained reads
	my$counts=$barcode_seqs{$barcode}=~tr/\n//;
	$retainedReads+=$counts/4;
	#add counts of retained barcode to hash
	$retainedBarcodeHash{$barcode}+=$counts;
	#delete sequence from hash
	delete $barcode_seqs{$barcode};
}

#print last sequences that did not match sample barcode
foreach my$unmatched (sort keys %barcode_seqs){
	print DISCARDED $barcode_seqs{$unmatched};
	#count discarded reads
	my$counts=$barcode_seqs{$unmatched}=~tr/\n//;
	$discardedReads+=$counts/4;
	#add counts of discarded barcode to hash
	$discardedBarcodeHash{$unmatched}+=$counts;
}
close DISCARDED;

#print summary info
print "\n";
print "Total Reads: $seqCount\n";

#print number of reads retained and number discarded
print "Retained Reads: $retainedReads\n";
print "Discarded Reads: $discardedReads\n";

#print files with summary data for discarded reads
#include hash of barcodes?  Useful for identifying if samples were switched...
#include hash of sequences?  May be memory intensive, can always do later on discarded read file

#create file for retained reads
my$retainedSummaryFile=$seqFileName."_DemultiplexGTseq_log.txt";
open(RETAINEDSUMMARY,">$retainedSummaryFile")||die "cannot open $retainedSummaryFile:$!";
print RETAINEDSUMMARY "Barcode File: $barcodeFile\n";
print RETAINEDSUMMARY "Sequence File: $sequenceFile\n";
print RETAINEDSUMMARY "\n";
print RETAINEDSUMMARY "Total Reads: $seqCount\n";
print RETAINEDSUMMARY "Retained Reads: $retainedReads\n";
print RETAINEDSUMMARY "Discarded Reads: $discardedReads\n";
print RETAINEDSUMMARY "\n";
#print RETAINEDSUMMARY "Target Barcodes:\n";

print RETAINEDSUMMARY "Sample\tBarcode\tReads\n";
foreach my$key (sort {$retainedBarcodeHash{$b} <=> $retainedBarcodeHash{$a}} keys %retainedBarcodeHash){
	print RETAINEDSUMMARY "$barcodes{$key}\t$key\t$retainedBarcodeHash{$key}\n";
}

print RETAINEDSUMMARY "\n";
print RETAINEDSUMMARY "Top 100 Discarded Barcodes:\n";
print RETAINEDSUMMARY "Barcode\tSequenceCount\n";

#create file for discarded reads
my$discardedSummaryFile=$seqFileName."_discardedSummary.txt";
open(DISCARDEDSUMMARY,">$discardedSummaryFile")||die "cannot open $discardedSummaryFile:$!";
print "Barcode\tSequenceCount\n";
my$barcodeCount=0;
foreach my$key (sort {$discardedBarcodeHash{$b} <=> $discardedBarcodeHash{$a}} keys %discardedBarcodeHash){
	$barcodeCount+=1;
	if($barcodeCount<=100){
		print RETAINEDSUMMARY "$key\t$discardedBarcodeHash{$key}\n";
	}
	print DISCARDEDSUMMARY "$key\t$discardedBarcodeHash{$key}\n";
}
close RETAINEDSUMMARY;
close DISCARDEDSUMMARY;

#print run time
my$endTime=time;
my$totalTime=$endTime-$startTime;
print "Total Time: $totalTime\n";