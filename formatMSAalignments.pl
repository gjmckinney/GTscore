#/usr/bin/perl -w
use strict;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;

#--------------------------------------------------------------------------------
#formatMSAalignments.pl
#
#Input:
#	Input consists of the MSA alignments produced by the alignMatchedSeqs function
#	in GTscore.R.  
#
#Usage:
#	This script is automatically called by alignMatchedSeqs in GTscore
#
#Garrett McKinney
#University of Washington
#gjmckinn@uw.edu
#	6/30/2018
#--------------------------------------------------------------------------------

my$alignmentFile=$ARGV[0];
my($fileName,$ext)=split '\.', $alignmentFile, 2;
my$excelFile=$fileName.".xlsx";

#initialize workbook
my$workbook=Excel::Writer::XLSX->new($excelFile);
my$worksheet=$workbook->add_worksheet();
my$format = $workbook->add_format();
$format->set_align( 'center' );


#open file and write to excel
open(ALIGNED,"<$alignmentFile")||die "cannot open $alignmentFile:$!";
my$header=<ALIGNED>;
my@seqs;
my@seqs2;
my$maxColumns=0;
my$Nrows=0;
while(my$line=<ALIGNED>){
	$Nrows++;
	chomp $line;
	my($ID,$seq)=split "\t", $line;
	my@bases=split "", $seq;
	if((scalar@bases+1)>$maxColumns){
		$maxColumns=(scalar@bases+1);
	}
	if($line=~/^\w/){
		push@seqs, $line;
	}else{
		push@seqs2, $line;
	}
}
push@seqs, @seqs2;

#add formula to generate consensus sequence with IUB codes underneath aligned sequences
my$rowNum=0;
#my$maxRow=($Nrows+1);
$worksheet->write($rowNum, 0, "Consensus");
foreach my$i (1..$maxColumns-1){
	#get rows/columns in excel notation for formula range
	my$minRow=xl_rowcol_to_cell(1,$i);
	my$maxRow=xl_rowcol_to_cell($Nrows,$i);
	my$massiveFormula='IF(COUNTIF('.$minRow.':'.$maxRow.',"-")>0,"-",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"N",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0),"V",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"H",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"D",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"B",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"W",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0),"S",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0),"M",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"K",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"Y",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0),"R",IF(COUNTIF('.$minRow.':'.$maxRow.',"T")>0,"T",IF(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,"G",IF(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,"C",IF(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,"A",0))))))))))))))))';
	$worksheet->write_formula($rowNum, $i, $massiveFormula)
}
$rowNum++;


#my$rowNum=0;
#my$rowNum=1;
#my$maxColumns=0;
#while(my$line=<ALIGNED>){
foreach my$line (@seqs){
	chomp $line;
	my($counts,$bases)=split "\t", $line;
	$worksheet->write($rowNum, 0, $counts);
	#replace - in probe sequence alignments with empty spaces so it isn't confused with indel
	if($counts=~/^probe/){
		$bases=~s/-/ /g;
	}
	my@bases=split "", $bases;
	##get number of columns
	#my$baseColNum=scalar@bases+1;
	#if($rowNum==0){
	#	$maxColumns=$baseColNum;
	#}elsif($baseColNum>$maxColumns){
	#	$maxColumns=$baseColNum;
	#}
	#write bases to file
	foreach my$i (0..$#bases){
		my$colNum=$i+1;
		$worksheet->write($rowNum, $colNum, $bases[$i],$format); 
	}
	$rowNum++;
}


#set narrower width for columns with bases
$worksheet->set_column( 1, $maxColumns, 2 );
#center bases in columns

##add formula to generate consensus sequence with IUB codes underneath aligned sequences
#$worksheet->write($rowNum, 0, "Consensus");
#foreach my$i (1..$maxColumns-1){
#	#get rows/columns in excel notation for formula range
#	my$minRow=xl_rowcol_to_cell(0,$i);
#	my$maxRow=xl_rowcol_to_cell($rowNum-1,$i);
#	my$massiveFormula='IF(COUNTIF('.$minRow.':'.$maxRow.',"-")>0,"-",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"N",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0),"V",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"H",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"D",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"B",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"W",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0),"S",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"C")>0),"M",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"K",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,COUNTIF('.$minRow.':'.$maxRow.',"T")>0),"Y",IF(AND(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,COUNTIF('.$minRow.':'.$maxRow.',"G")>0),"R",IF(COUNTIF('.$minRow.':'.$maxRow.',"T")>0,"T",IF(COUNTIF('.$minRow.':'.$maxRow.',"G")>0,"G",IF(COUNTIF('.$minRow.':'.$maxRow.',"C")>0,"C",IF(COUNTIF('.$minRow.':'.$maxRow.',"A")>0,"A",0))))))))))))))))';
#	$worksheet->write_formula($rowNum, $i, $massiveFormula)
#}

#create format for each DNA bases
my$Aformat=$workbook->add_format(
	bg_color => '#C6EFCE',
	color    => '#006100');
my$Cformat=$workbook->add_format(
	bg_color => '#00B0F0',
	color    => '#0F243E');
my$Gformat=$workbook->add_format(
	bg_color => '#FFC7CE',
	color    => '#9C0006');
my$Tformat=$workbook->add_format(
	bg_color => '#FFEB9C',
	color    => '#9C6500');
my$indelFormat=$workbook->add_format(
	bg_color => '#F2F2F2');

#write conditional formats to cover the dataset for each base, get dimensions from intput file
$worksheet->conditional_formatting( 0,1,$rowNum,$maxColumns,
	{
		type		=> 'text',
		criteria	=> 'containing',
		value		=> 'A',
		format		=> $Aformat,
	}
);

$worksheet->conditional_formatting( 0,1,$rowNum,$maxColumns,
	{
		type		=> 'text',
		criteria	=> 'containing',
		value		=> "C",
		format		=> $Cformat,
	}
);

$worksheet->conditional_formatting( 0,1,$rowNum,$maxColumns,
	{
		type		=> 'text',
		criteria	=> 'containing',
		value		=> "G",
		format		=> $Gformat,
	}
);

$worksheet->conditional_formatting( 0,1,$rowNum,$maxColumns,
	{
		type		=> 'text',
		criteria	=> 'containing',
		value		=> "T",
		format		=> $Tformat,
	}
);

$worksheet->conditional_formatting( 0,1,$rowNum,$maxColumns,
	{
		type		=> 'text',
		criteria	=> 'containing',
		value		=> "-",
		format		=> $indelFormat,
	}
);

$workbook->close();