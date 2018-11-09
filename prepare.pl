use warnings;
use strict;

my $fasta_file="/well/htseq/Genomes/BWA/GRCh37";
my $mut_file="/data/ted/COAD/esophagael-sorted.txt";


open FASTA, $fasta_file or die "Can't open fasta file";

my @genome;
my $temp=[];
my $chr=0;

while(<FASTA>){
	chomp;
	if(/^>/){
		$chr++;
		push @genome, $temp ;
		$temp=[];
		last if $chr>24;
		next;
	}
	
	push @{$temp}, split "",$_;
	
}

open MUT,$mut_file or die "Can't open mutation file $mut_file";

my $dir="/data/ted/multi/eso/";
open SEQ, ">${dir}seq.txt" or die "Can't open output file ${dir}seq.txt";
open LABEL, ">${dir}label.txt" or die "Can't open output file ${dir}label.txt";

my ($old_pos,$pos);
my $old_chr=0;
my $Base;

sub Base_to_number($){
	$_=shift;
	return 0 if $_ eq "A";
	return 1 if $_ eq "G";
	return 2 if $_ eq "C";
	return 3 if $_ eq "T";
	return 4;
}


my ($sample_name,$ori,$base);


while(<MUT>){
	chomp;
	($sample_name,$chr,$pos,$ori,$base)=(split "\t", $_);
	$base=(split "",$base)[0];
	$Base=Base_to_number($base);
	next, if $Base>3;
	next, if Base_to_number($ori)>3;
	next, if ($chr eq 'X' or $chr eq 'Y');
	if($chr!=$old_chr){
		$old_chr=$chr;
		$old_pos=$pos;
		next;
	}
	print SEQ (join "",@{$genome[$chr]}[($pos-500)..($pos+500)]), "\n", if $pos!=$old_pos;
	
	print "ERROR: in $chr,$pos, $ori not ${$genome[$chr]}[$pos-1] \n", if ${$genome[$chr]}[$pos-1] ne $ori;
	
	print LABEL 4+$Base, "\n", if $pos-$old_pos<=1000 and $pos!=$old_pos;
	print LABEL $Base, "\n", if $pos-$old_pos>1000;
	$old_pos=$pos;
}











