#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to filter the potential artefact polyA reads from a BAM file. Samtools must be installed; Most of the criteria were hard coded at the beginning
#	of the script
#
#	Input
#		--BAMPath=		file path; [obligatory]; Path of a BAM file;
#		--gffPath=		file path; [obligatory]; Path of the GFF that specify the introns;
#		--fastaPath		file path; [obligatory]; Path of the fasta of the reference genome;
#		--outDir=		dir path; output dir;
#
#	Usage
#		
#		perl BAMPolyAFilterer_v0.1.pl --BAMPath=/Volumes/C_Analysis/NGS/results/EHI_polyA_mapping_basic/EHI_polyA_mapping/tophat/allNoSec.bam --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v2_allFearues.forPileupCounter.gff  --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa
#
#	Assumption
#
#	History:
#		
#		v0.1
#		-debut, mainly based on BAM filterer;
#
#		v0.2
#			[Tue 27 Aug 2013 15:33:23 CEST] cleaned using perlScriptCleaner
#
#		v0.3
#			[09/10/2013 10:57] parameter $numTopALenRdForAvg added: the average length of A tails will be calculated using the top N reads with longest A tail instead of all read;
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-10-09 12:34]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/SAMFileHanding/BAMPolyAFilterer/v0.3/BAMPolyAFilterer_v0.3.pl --BAMPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/tophatPipeMapping/EHI_polyA_mapping/tophat/allNoSec.bam --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --outDir=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/tophatPipeMapping/EHI_polyA_mapping/BAMPolyAFilterer
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/SAMFileHanding/BAMPolyAFilterer/v0.3/BAMPolyAFilterer_v0.3.pl
#	--BAMPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/tophatPipeMapping/EHI_polyA_mapping/tophat/allNoSec.bam
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--outDir=/Volumes/A_MPro2TB/NGS/analyses/EHI_polyA/tophatPipeMapping/EHI_polyA_mapping/BAMPolyAFilterer
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: checkSamtoolsVersion|301, printCMDLogOrFinishMessage|452, readParameters|843
#	secondaryDependOnSub: currentTime|321, reportStatus|877
#
#<section ID="startingTasks" num="0">
#----------Read parameters ----------#
&printCMDLogOrFinishMessage("CMDLog");#->452
my ($BAMPath, $gffPath, $fastaPath, $outDir) = &readParameters();#->843
&checkSamtoolsVersion();#->301
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
#----hard coded the parameters
my $downANumRegion = 10; #---the size of the downstream region to check for the number of As 
my $downANumMax = 8; #---the maximum continous track of As in downANumRegion
my $downATrackRegion = 6; #---the size of the downstream region to check for continous track of As 
my $downATrackMax = 6; #---the maximum continous track of As in downATrackRegion
my $upANumRegion = 10; #---the size of the upstream region to check for the number of As 
my $upANumMax = 8; #---the maximum continous track of As in upANumRegion 
my $minMatch3EndBlock = 2; #---the minimum number of match in the 3'end block
my $minTerminalBlockLength = 5; #---the minimum size of the 3'end block
my $dynamicDownMaxAPct = 80; #---rm polyA read dynamically based on length of polyA tail;
my $readAPctMax = 70; #---maximum number of As within read
my $minimumLength = 20; #---read shorter than this length will be discarded
my $minRdATailLen = 5; #---read has polyA tail short than this length be discarded
my $nonRedundant = 'no'; #---will remove non-redudant read, keep only the read with longest A tail and smallest NM;
my $numTopALenRdForAvg = 3;#---only this number of reads with the longest A-tail length will be used to calculate the average A tail length

my $paramTag = join ".", (
	"polyA",
	"NR.$nonRedundant",
	"DR$downANumRegion.$downANumMax",
	"UR$upANumRegion.$upANumMax",
	"DT$downATrackRegion.$downATrackMax",
	"EM$minMatch3EndBlock",
	"TB$minTerminalBlockLength",
	"DD$dynamicDownMaxAPct",
	"RP$readAPctMax",
	"ML$minimumLength",
	"AL$minRdATailLen",
	"TA$numTopALenRdForAvg",
);
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag/"; push @mkDirAry, $resultDir;
my $resultBamDir = "$resultDir/bam/"; push @mkDirAry, $resultBamDir;
my $resultStorableDir = "$resultDir/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_readBasicGenomeInfo
#	primaryDependOnSub: generateTwoWaysJunctionIndex|388, readGff|671, readMultiFasta|789
#	secondaryDependOnSub: reportStatus|877
#
#<section ID="readBasicGenomeInfo" num="3">
my ($nameByGeneHsh_ref, $strandByGeneHsh_ref, $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $geneCtgryHsh_ref, $intronRngByGeneHsh_ref) = &readGff($gffPath);#->671
my $fastaHsh_ref = &readMultiFasta($fastaPath);#->789
my $junctIndxHsh_ref = &generateTwoWaysJunctionIndex($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref);#->388
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_filterBAM
#	primaryDependOnSub: readAndFilterBAMOnTheFly|485
#	secondaryDependOnSub: calculateSitePolyASeqInfo|206, dynmaicCheckPolyADnStrmRegionBasedOnReadLength|339, getRead3EndAndMatchBlockLenAndLongestGap|417, reportStatus|877
#
#<section ID="filterBAM" num="4">
&readAndFilterBAMOnTheFly($BAMPath, $junctIndxHsh_ref, $fastaHsh_ref, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax, $minMatch3EndBlock, $minTerminalBlockLength, $dynamicDownMaxAPct, $readAPctMax, $minimumLength, $minRdATailLen, $nonRedundant, $numTopALenRdForAvg, $resultBamDir, $resultStorableDir);#->485
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|452
#	secondaryDependOnSub: currentTime|321
#
#<section ID="finishingTasks" num="5">
&printCMDLogOrFinishMessage("finishMessage");#->452
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	checkTools [n=1]:
#		checkSamtoolsVersion
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=6]:
#		checkSamtoolsVersion, currentTime, printCMDLogOrFinishMessage
#		readMultiFasta, readParameters, reportStatus
#
#	reporting [n=1]:
#		currentTime
#
#	unassigned [n=6]:
#		calculateSitePolyASeqInfo, dynmaicCheckPolyADnStrmRegionBasedOnReadLength, generateTwoWaysJunctionIndex
#		getRead3EndAndMatchBlockLenAndLongestGap, readAndFilterBAMOnTheFly, readGff
#
#====================================================================================================================================================#

sub calculateSitePolyASeqInfo {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readAndFilterBAMOnTheFly|485
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_filterBAM|159
#	input: $cntg, $curntReadStart, $downANumMax, $downANumRegion, $downATrackMax, $downATrackRegion, $fastaHsh_ref, $junctIndxHsh_ref, $polyAPreCalSeqHsh_ref, $read3End, $readStrand, $upANumMax, $upANumRegion
#	output: none
#	toCall: &calculateSitePolyASeqInfo($fastaHsh_ref, $polyAPreCalSeqHsh_ref, $junctIndxHsh_ref, $curntReadStart, $readStrand, $cntg, $read3End, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax);
#	calledInLine: 218, 598
#....................................................................................................................................................#
	
	#&calculateSitePolyASeqInfo($fastaHsh_ref, $polyAPreCalSeqHsh_ref, $junctIndxHsh_ref, $curntReadStart, $readStrand, $cntg, $read3End, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax);#->206

	my ($fastaHsh_ref, $polyAPreCalSeqHsh_ref, $junctIndxHsh_ref, $curntReadStart, $readStrand, $cntg, $read3End, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax) = @_;

	my ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean);
	
	#------define the start and end region then get the sequences
	my $cntgSeq = $fastaHsh_ref->{$cntg};
	my $down20ntStart = my $up20ntStart = 0;
	my $down20ntGenSeq = my $up20ntGenSeq = "";
	if ($readStrand eq "+") {
		$down20ntStart = $read3End;
		$up20ntStart = $read3End - 20 - 1;
		$down20ntGenSeq = substr $cntgSeq, $down20ntStart, 20;
		$up20ntGenSeq = substr $cntgSeq, $up20ntStart, 20;

	} elsif ($readStrand eq "-") {
		$down20ntStart = $curntReadStart - 20 - 1;
		$up20ntStart = $curntReadStart;
		$down20ntGenSeq = substr $cntgSeq, $down20ntStart, 20;
		$up20ntGenSeq = substr $cntgSeq, $up20ntStart, 20;
		$down20ntGenSeq = reverse($down20ntGenSeq); $down20ntGenSeq =~ tr/ACGTacgt/TGCAtgca/;
		$up20ntGenSeq = reverse($up20ntGenSeq); $up20ntGenSeq =~ tr/ACGTacgt/TGCAtgca/;

	} else {
		die "undefined read strand\n";
	}

	my $downRegNumCheckSeq = substr $down20ntGenSeq, 0, $downANumRegion;
	my $downRegTrackCheckSeq = substr $down20ntGenSeq, 0, $downATrackRegion;
	my $upRegNumCheckSeq = substr $up20ntGenSeq, 0, $upANumRegion;
	
	#----check whether the read ends close to junction
	my $nextToJunct = 'no';
	my $down20ntcDNA = '';
	
	my @searchRngAry = ($read3End..($read3End+20));
	@searchRngAry = (($read3End-20)..$read3End) if ($readStrand eq "-");
	my $junctStart = my $junctEnd = 0;
	foreach my $srchPos (@searchRngAry) {#---check if the junction is around
		if (exists $junctIndxHsh_ref->{$cntg}{$readStrand}{$srchPos}) {
			$junctStart = $srchPos;
			$junctEnd = $junctIndxHsh_ref->{$cntg}{$readStrand}{$srchPos};
			$nextToJunct = "yes";
			last;
		}
	}

	#---get the "spliced" sequences, head half and the tail half
	if ($nextToJunct eq "yes") {
		if ($readStrand eq "+") {
			my $leftHalfLength = $junctStart - $read3End - 1;
			my $leftHalfSeq = substr $cntgSeq, $read3End, $leftHalfLength;
			my $rightHalfLength = 20 - $leftHalfLength;
			my $rightHalfSeq = substr $cntgSeq, $junctEnd, $rightHalfLength;
			$down20ntcDNA = $leftHalfSeq.$rightHalfSeq;
	
		} elsif ($readStrand eq "-") {
			my $rightHalfLength = $read3End - $junctStart - 1;
			my $rightHalfSeq = substr $cntgSeq, $junctStart, $rightHalfLength;
			my $leftHalfLength = 20 - $rightHalfLength;
			my $leftHalfSeq = substr $cntgSeq, ($junctEnd - $leftHalfLength-1) , $leftHalfLength;
			$down20ntcDNA = $leftHalfSeq.$rightHalfSeq;
			$down20ntcDNA = reverse($down20ntcDNA); $down20ntcDNA =~ tr/ACGTacgt/TGCAtgca/;
			
		} else {
			die "undefined read strand\n";
		}
	}

	$withinJunct20ntBoolean = $nextToJunct;
	$genomicCheckSeq20nt = $down20ntGenSeq;
	$cDNACheckSeq20nt = $down20ntcDNA;
	$downATrackBoolean = 'no';
	$downANumBoolean = 'no';
	$upANumBoolean = 'no';
	$downANumBoolean = 'yes' if ((($downRegNumCheckSeq =~ tr/A//) > $downANumMax) and ($downANumMax > 0));
	$upANumBoolean = 'yes' if ((($upRegNumCheckSeq =~ tr/A//) > $upANumMax) and ($upANumMax > 0));
	$downATrackBoolean = 'yes' if (($downRegTrackCheckSeq =~ m/([A]{$downATrackMax,})/) and ($downATrackMax > 0));
	
	@{$polyAPreCalSeqHsh_ref->{$cntg}{$read3End}{$readStrand}} = ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean);
}
sub checkSamtoolsVersion {
#....................................................................................................................................................#
#	subroutineCategory: checkTools, general
#	dependOnSub: reportStatus|877
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80
#	secondaryAppearInSection: >none
#	input: none
#	output: none
#	toCall: &checkSamtoolsVersion();
#	calledInLine: 88
#....................................................................................................................................................#

	my $samtoolsStdout = `samtools 2>&1`;
	if ($samtoolsStdout =~ m/\s+(Version: \S+)\s+/) {
		&reportStatus("Checking: samtools: $1", 0, "\n");#->877
	} else {
		die 'samtools not installed properly. Quitting.';
	}
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|452, reportStatus|877
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|80, 5_finishingTasks|169
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 472, 475, 480, 893
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub dynmaicCheckPolyADnStrmRegionBasedOnReadLength {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readAndFilterBAMOnTheFly|485
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_filterBAM|159
#	input: $cDNACheckSeq20nt, $downANumMax, $downANumRegion, $downATrackMax, $downATrackRegion, $dynamicDownMaxAPct, $genomicCheckSeq20nt, $rdATailLen, $rdName, $withinJunct20ntBoolean
#	output: $artefactPolyADynReadLength
#	toCall: my ($artefactPolyADynReadLength) = &dynmaicCheckPolyADnStrmRegionBasedOnReadLength($rdName, $dynamicDownMaxAPct, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $withinJunct20ntBoolean, $rdATailLen, $downANumRegion, $downATrackRegion, $downATrackMax, $downANumMax);
#	calledInLine: 607
#....................................................................................................................................................#

	my ($rdName, $dynamicDownMaxAPct, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $withinJunct20ntBoolean, $rdATailLen, $downANumRegion, $downATrackRegion, $downATrackMax, $downANumMax) = @_;

	my $downADynRegion = my $downADynMax = undef;
	
	my $artefactPolyADynReadLength = 'no';
	
	$rdATailLen = 20 if ($rdATailLen > 20);

	$downADynRegion = $rdATailLen;

	if ($rdATailLen > 10) {
		$downADynMax = sprintf "%.0f", $downADynRegion*($dynamicDownMaxAPct/100);
	} elsif (($rdATailLen > 5) and ($rdATailLen <=10)) {
		$downADynMax = sprintf "%.0f", $downADynRegion*(40/100);
	} else {#---<=5
		$downADynMax = sprintf "%.0f", $downADynRegion*(20/100);
	}
	
	#----check the dyn seq region
	my $downRegDynCheckSeq = substr $genomicCheckSeq20nt, 0, $downADynRegion;
	$artefactPolyADynReadLength = "yes" if ((($downRegDynCheckSeq =~ tr/A//) > $downADynMax) and ($downADynMax > 0));
	
	#----check the cDNA is close to junction
	if ($withinJunct20ntBoolean eq "yes") {
		my $downRegNumCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downANumRegion;
		my $downRegTrackCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downATrackRegion;
		my $downRegDynCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downADynRegion;
		if (((($downRegDynCheckcDNASeq =~ tr/A//) > $downADynMax) and ($downADynMax > 0))
		or ((($downRegNumCheckcDNASeq =~ tr/A//) > $downANumMax) and ($downANumMax > 0))
		or (($downRegTrackCheckcDNASeq =~ m/([A]{$downATrackMax,})/) and ($downATrackMax > 0))) {
			$artefactPolyADynReadLength = "yes";
		}
	}

	return $artefactPolyADynReadLength;
}
sub generateTwoWaysJunctionIndex {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_readBasicGenomeInfo|147
#	secondaryAppearInSection: >none
#	input: $cntgByGeneHsh_ref, $intronRngByGeneHsh_ref
#	output: $junctIndxHsh_ref
#	toCall: my ($junctIndxHsh_ref) = &generateTwoWaysJunctionIndex($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref);
#	calledInLine: 154, 400
#....................................................................................................................................................#
	
	#my $junctIndxHsh_ref &generateTwoWaysJunctionIndex($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref);#->388
	my ($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref) = @_;
	
	my $junctIndxHsh_ref = {};
	
	foreach my $geneID (keys %{$intronRngByGeneHsh_ref}) {
		my $cntg = $cntgByGeneHsh_ref->{$geneID};
		foreach my $intronID (keys %{$intronRngByGeneHsh_ref->{$geneID}}) {
			my $intronStart = $intronRngByGeneHsh_ref->{$geneID}{$intronID}{"start"};
			my $intronEnd = $intronRngByGeneHsh_ref->{$geneID}{$intronID}{"end"};
			$junctIndxHsh_ref->{$cntg}{"+"}{$intronStart} = $intronEnd;
			$junctIndxHsh_ref->{$cntg}{"-"}{$intronEnd} = $intronStart;
		}
	}
	return $junctIndxHsh_ref;
}
sub getRead3EndAndMatchBlockLenAndLongestGap {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: readAndFilterBAMOnTheFly|485
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_filterBAM|159
#	input: $cigarStr, $curntReadStart, $readStrand
#	output: $longestGap, $matchLengthAry_ref, $read3End
#	toCall: my ($read3End, $matchLengthAry_ref, $longestGap) = &getRead3EndAndMatchBlockLenAndLongestGap($cigarStr, $curntReadStart, $readStrand);
#	calledInLine: 575
#....................................................................................................................................................#

	my ($cigarStr, $curntReadStart, $readStrand) = @_;
	
	#---Get $read3End
	my $genomicLength = 0;
	my $matchLengthAry_ref = ();
	my $longestGap = 0;

	while ($cigarStr =~ /(\d+)M/g) {
		$genomicLength += $1;
		push @{$matchLengthAry_ref}, $1;
	}

	while ($cigarStr =~ /(\d+)([N|D])/g) {
		$genomicLength += $1;
		$longestGap = $1 if $1 > $longestGap;
	}

	my $read3End = $curntReadStart + $genomicLength - 1;
	$read3End = $curntReadStart if ($readStrand eq "-");
	
	return ($read3End, $matchLengthAry_ref, $longestGap);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|321
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80, 5_finishingTasks|169
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 86, 174
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->321
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->321
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->321
		print "=========================================================================\n\n";
	}
}
sub readAndFilterBAMOnTheFly {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: calculateSitePolyASeqInfo|206, dynmaicCheckPolyADnStrmRegionBasedOnReadLength|339, getRead3EndAndMatchBlockLenAndLongestGap|417, reportStatus|877
#	appearInSub: >none
#	primaryAppearInSection: 4_filterBAM|159
#	secondaryAppearInSection: >none
#	input: $BAMPath, $downANumMax, $downANumRegion, $downATrackMax, $downATrackRegion, $dynamicDownMaxAPct, $fastaHsh_ref, $junctIndxHsh_ref, $minMatch3EndBlock, $minRdATailLen, $minTerminalBlockLength, $minimumLength, $nonRedundant, $numTopALenRdForAvg, $readAPctMax, $resultBamDir, $resultStorableDir, $upANumMax, $upANumRegion
#	output: none
#	toCall: &readAndFilterBAMOnTheFly($BAMPath, $junctIndxHsh_ref, $fastaHsh_ref, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax, $minMatch3EndBlock, $minTerminalBlockLength, $dynamicDownMaxAPct, $readAPctMax, $minimumLength, $minRdATailLen, $nonRedundant, $numTopALenRdForAvg, $resultBamDir, $resultStorableDir);
#	calledInLine: 164
#....................................................................................................................................................#

	my ($BAMPath, $junctIndxHsh_ref, $fastaHsh_ref, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax, $minMatch3EndBlock, $minTerminalBlockLength, $dynamicDownMaxAPct, $readAPctMax, $minimumLength, $minRdATailLen, $nonRedundant, $numTopALenRdForAvg, $resultBamDir, $resultStorableDir) = @_;

	
	my $outBAMPath = "$resultBamDir/polyA.filtered.bam";
	open BAMOUT, "| samtools view -b -S - >$outBAMPath 2>/dev/null";

	&reportStatus("Checking BAM file size", 0, "\n");#->877
	my $flagStatOut = `samtools flagstat $BAMPath`;
	my ($totalReadNum) = split / /, $flagStatOut;

	#----print the headers
	open BAMHEADER, "samtools view -H $BAMPath |";#---Header only
	print BAMOUT $_	while <BAMHEADER>;
	close BAMHEADER;
	open BAMIN, "samtools view $BAMPath |"; #---no header

	&reportStatus("Start filtering BAM", 0, "\n");#->877

	my $procRead = my $passedRead = my $read3EndPosNum = 1;
	my $polyAPreCalSeqHsh_ref = {}; ###empty, will be filled along the way
	
	my $lastReadStart = 'inititation';
	my $nonRedundantLineHsh_ref = {};
	my $redundancyCountHsh_ref = {};
	my $tmpATailLenByPosAryHsh_ref = {};
	
	while (my $theLine = <BAMIN>) { 

		#last if ($passedRead >= 100000); #----adhoc debug

		chomp $theLine;
		my ($rdName, $flag, $cntg, $curntReadStart, $mapQ, $cigarStr, $cntgNext, $posNext, $len, $readSeq, $qual) = split /\t/, $theLine;
		$procRead++;
		
		if (($procRead % 100000 == 0) or ($procRead == $totalReadNum)) {
			my $passedReadPct = sprintf "%.5f", 100*$passedRead/$procRead;
			my $procReadPct = sprintf "%.5f", 100*$procRead/$totalReadNum;
			my $read3EndPosNumPct = sprintf "%.5f", 100*$read3EndPosNum/$procRead;
			&reportStatus("$passedReadPct\% passed || $procReadPct\% processed || $read3EndPosNumPct\% new site", 20, "\r");#->877
		}
		
		my $length = length $readSeq;
		my $readStrand = '+';
		$readStrand = '-' if ($flag & 16); ###---reference: http://seqanswers.com/forums/showthread.php?t=2301
		my $rdATailLen = undef;

		if ($rdName =~ m/\[Ax(\d+)\]/) {#---defined polyA
			$rdATailLen = $1;
		} else {
			die "Wrong polyA read naming format, no polyA length indicated. Qutting\n";
		}
		
		#############################################
		###############check for length###############
		next if $length < $minimumLength;

		##################################################
		###############check for rdPolyANum###############
		next if $rdATailLen < $minRdATailLen;

		###############################################
		###############check A rich read###############
		my $ARichRead = 'no';
		my $readANum = 0;
		
		if ($readStrand eq '+') {
			$readANum = $readSeq =~ tr/A//;
		} elsif ($readStrand eq '-') {
			$readANum = $readSeq =~ tr/T//;
		}
		
		my $readAPct = sprintf "%.2f", 100*($readANum/$length);
		$ARichRead = 'yes' if ($readAPct > $readAPctMax);
		next if $ARichRead eq 'yes';

		################################################
		###############shortTerminalBlock###############
		my $shortTerminalBlock = 'no';
		my ($read3End, $matchLengthAry_ref, $longestGap) = &getRead3EndAndMatchBlockLenAndLongestGap($cigarStr, $curntReadStart, $readStrand);#->417
		$shortTerminalBlock = 'yes' if (($matchLengthAry_ref->[-1] < $minTerminalBlockLength) or ($matchLengthAry_ref->[0] < $minTerminalBlockLength)); #----end match block < 5
		next if $shortTerminalBlock eq 'yes';

		###############################################
		###############end3BlockMismatch###############
		my $end3BlockMismatch = 'no';
		my $NM = 0;
		$NM = $1 if ($theLine =~ m/\tNM:i:(\d+)\t/);
		if ($NM > 0) {
			my $MD = $1 if ($theLine =~ m/\tMD:Z:(\S+)\t/);
			my @MDSplt = split /\D+/, $MD;
			if (($readStrand eq '+' and $MDSplt[-1] < $minMatch3EndBlock) or ($readStrand eq '-' and $MDSplt[0] < $minMatch3EndBlock)) {
				$end3BlockMismatch = 'yes';
				#print $readStrand."\t".$MD."\n";
			}
		}
		next if $end3BlockMismatch eq 'yes';

		###########################################
		###############artefactPolyA###############
		my $artefactPolyA = 'no';
		if (not $polyAPreCalSeqHsh_ref->{$cntg}{$read3End}{$readStrand}) {
			&calculateSitePolyASeqInfo($fastaHsh_ref, $polyAPreCalSeqHsh_ref, $junctIndxHsh_ref, $curntReadStart, $readStrand, $cntg, $read3End, $downANumRegion, $downATrackRegion, $upANumRegion, $downANumMax, $upANumMax, $downATrackMax);#->206
			$read3EndPosNum++;
		}
		my ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean) = @{$polyAPreCalSeqHsh_ref->{$cntg}{$read3End}{$readStrand}};
		$artefactPolyA = 'yes' if (($downATrackBoolean eq "yes") or ($downANumBoolean eq "yes") or ($upANumBoolean eq "yes"));
		next if $artefactPolyA eq 'yes';

		########################################################
		###############artefactPolyADynReadLength###############
		my $artefactPolyADynReadLength = &dynmaicCheckPolyADnStrmRegionBasedOnReadLength($rdName, $dynamicDownMaxAPct, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $withinJunct20ntBoolean, $rdATailLen, $downANumRegion, $downATrackRegion, $downATrackMax, $downANumMax);#->339
		next if $artefactPolyADynReadLength eq 'yes';

		$passedRead++;
		
		if ($nonRedundant eq 'yes') {

			$nonRedundantLineHsh_ref->{$curntReadStart}{$length}{$readStrand}{$rdATailLen}{$NM} = $theLine;

			#----change read start position, sort out the redudancy
			if (($curntReadStart ne $lastReadStart and $lastReadStart ne 'inititation') or ($procRead == $totalReadNum)){
				
				#---consider the last line, curnt pos have to print
				my $redundancy = 0;
				my @NRReadStartAry = $lastReadStart;
				push @NRReadStartAry, $curntReadStart if ($curntReadStart ne $lastReadStart and $procRead == $totalReadNum);
				foreach my $NRReadStart (@NRReadStartAry) {
					foreach my $NRLength (keys %{$nonRedundantLineHsh_ref->{$NRReadStart}}) {
						foreach my $NRStrand (keys %{$nonRedundantLineHsh_ref->{$NRReadStart}{$NRLength}}) {
							foreach my $NRRdPolyANum (sort {$b <=> $a} keys %{$nonRedundantLineHsh_ref->{$NRReadStart}{$NRLength}{$NRStrand}}) {
								foreach my $NRNM (sort {$a <=> $b} keys %{$nonRedundantLineHsh_ref->{$NRReadStart}{$NRLength}{$NRStrand}{$NRRdPolyANum}}) {
									push @{$tmpATailLenByPosAryHsh_ref->{$cntg}{$read3End}{$readStrand}}, $rdATailLen;
									print BAMOUT "$nonRedundantLineHsh_ref->{$NRReadStart}{$NRLength}{$NRStrand}{$NRRdPolyANum}{$NRNM}.\n";
									last;#---will only print the longest polyA
								}
								last;#---will only print the smallest NM
							}
						}
					}
					delete $nonRedundantLineHsh_ref->{$NRReadStart};
				}
			}

			$lastReadStart = $curntReadStart;

		} else {
			push @{$tmpATailLenByPosAryHsh_ref->{$cntg}{$read3End}{$readStrand}}, $rdATailLen;
			print BAMOUT "$theLine\n";
		}
	}
	close BAMOUT;
	
	my $polyACntgPosAvgTailLenHsh_ref = {};
	&reportStatus("Calculating average polyA tail length", 20, "\n");#->877
	foreach my $cntg (keys %{$tmpATailLenByPosAryHsh_ref}) {
		foreach my $read3End (keys %{$tmpATailLenByPosAryHsh_ref->{$cntg}}) {
			foreach my $readStrand (keys %{$tmpATailLenByPosAryHsh_ref->{$cntg}{$read3End}}) {
				my @sortedALenAry = sort {$b <=> $a} @{$tmpATailLenByPosAryHsh_ref->{$cntg}{$read3End}{$readStrand}};
				my $num = @sortedALenAry;
				my $maxIdx = $numTopALenRdForAvg-1;
				$maxIdx = $num -1 if $maxIdx > $num-1;
				my $ALen = sprintf "%.2f", sum(@sortedALenAry[0..$maxIdx])/($maxIdx+1);
				$polyACntgPosAvgTailLenHsh_ref->{$cntg}{$read3End}{$readStrand}{'num'} = $num;
				$polyACntgPosAvgTailLenHsh_ref->{$cntg}{$read3End}{$readStrand}{'ALen'} = $ALen;
				delete $tmpATailLenByPosAryHsh_ref->{$cntg}{$read3End}{$readStrand};
			}
		}
	}
	
	store($polyACntgPosAvgTailLenHsh_ref, "$resultStorableDir/polyACntgPosAvgTailLenHsh.pls");
	
	system "samtools index $outBAMPath";
}
sub readGff {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|877
#	appearInSub: >none
#	primaryAppearInSection: 3_readBasicGenomeInfo|147
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $geneCtgryHsh_ref, $intronRngByGeneHsh_ref, $nameByGeneHsh_ref, $strandByGeneHsh_ref
#	toCall: my ($nameByGeneHsh_ref, $strandByGeneHsh_ref, $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $geneCtgryHsh_ref, $intronRngByGeneHsh_ref) = &readGff($gffPath);
#	calledInLine: 152
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $strandByGeneHsh_ref = {};
	my $cntgByGeneHsh_ref = {};
	my $nameByGeneHsh_ref = {};
	my $exonRngByGeneHsh_ref = {};
	my $geneCtgryHsh_ref = {};
	my $intronRngByGeneHsh_ref = {};

	#---read the gff
	my $geneByRNAHsh_ref = {};
	open (GFF, $gffPath);
	&reportStatus("Reading: $gffPath", 0, "\n");#->877
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				$strandByGeneHsh_ref->{$geneID} = $geneStrd;
				$cntgByGeneHsh_ref->{$geneID} = $seq;
				$nameByGeneHsh_ref->{$geneID} = $geneName;

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				# The CDS is ignored at the moment, until it reaches the point that we are looking at UTRs
				#
				#my $mRNAID = $parent;
				#my $geneID = $geneByRNAHsh{$mRNAID};
				#$CDSCountHsh{$geneID}++;
				#my $CDSCount = $CDSCountHsh{$geneID};

				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"} = $featureStart;
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} = $featureEnd;
			 	#$geneCDSLenHsh{$geneID} = 0 if $CDSCount == 1; #---define the length hashfor the 1st time
			 	#$geneCDSLenHsh{$geneID} += ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} - ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"};
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $exonID = $unqID;
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				next if not $geneID;
				$exonRngByGeneHsh_ref->{$geneID}{$exonID}{"start"} = $featureStart;
				$exonRngByGeneHsh_ref->{$geneID}{$exonID}{"end"} = $featureEnd;

			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				next if not $geneID;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneCtgryHsh_ref->{$geneID} = $geneCategory;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;

	#---generate intron
	my $boundsForIntronHsh_ref = {};
	foreach my $geneID (keys %{$exonRngByGeneHsh_ref}) {
		my $exonCount = keys %{$exonRngByGeneHsh_ref->{$geneID}};
		my $seq = $cntgByGeneHsh_ref->{$geneID};
		foreach my $exonID (keys %{$exonRngByGeneHsh_ref->{$geneID}}) {
			if ($exonCount > 1) {
				push @{$boundsForIntronHsh_ref->{$geneID}}, $exonRngByGeneHsh_ref->{$geneID}{$exonID}{"start"};
				push @{$boundsForIntronHsh_ref->{$geneID}}, $exonRngByGeneHsh_ref->{$geneID}{$exonID}{"end"};
			}
		}
	}

	#---define the introns ranges
	my $geneIntronNum = keys %{$boundsForIntronHsh_ref};
	&reportStatus("$geneIntronNum gene were found to contain intron. Storing intron boundaries", 0, "\n");#->877

	foreach my $geneID (keys %{$boundsForIntronHsh_ref}) {
		my @sortedBounds = sort {$a <=> $b} @{$boundsForIntronHsh_ref->{$geneID}};
		my $boundNum = @sortedBounds;
		
		my $intronNum = 0;
		for (my $i = 1; $i < ($boundNum - 1); $i += 2) {
			$intronNum++;
			my $intronID = $geneID.$intronNum;
			$intronRngByGeneHsh_ref->{$geneID}{$intronID}{"start"} = $sortedBounds[$i] + 1;
			$intronRngByGeneHsh_ref->{$geneID}{$intronID}{"end"} = $sortedBounds[$i+1] - 1;
		}
	}
	
	return ($nameByGeneHsh_ref, $strandByGeneHsh_ref, $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $geneCtgryHsh_ref, $intronRngByGeneHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|877
#	appearInSub: >none
#	primaryAppearInSection: 3_readBasicGenomeInfo|147
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 153
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->877
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/ +/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|80
#	secondaryAppearInSection: >none
#	input: none
#	output: $BAMPath, $fastaPath, $gffPath, $outDir
#	toCall: my ($BAMPath, $gffPath, $fastaPath, $outDir) = &readParameters();
#	calledInLine: 87
#....................................................................................................................................................#
	
	my ($BAMPath, $gffPath, $fastaPath, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/BAMPolyAFilterer/";

	GetOptions 	("BAMPath=s" => \$BAMPath,
				 "gffPath=s"  => \$gffPath,
				 "fastaPath=s"  => \$fastaPath,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($BAMPath, $gffPath, $fastaPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($BAMPath, $gffPath, $fastaPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|321
#	appearInSub: checkSamtoolsVersion|301, readAndFilterBAMOnTheFly|485, readGff|671, readMultiFasta|789
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|80, 3_readBasicGenomeInfo|147, 4_filterBAM|159
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 315, 503, 513, 535, 650, 695, 771, 807
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->321

	return ();
}

exit;
