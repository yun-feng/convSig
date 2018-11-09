use warnings;
use strict;

my $seq_file="/data/ted/multi/eso/seq_open.txt";
my $seq_bg_file="/data/ted/multi/eso/bg_open.txt";
my $out_file=">/data/ted/multi/eso/bg_cos400.txt";

sub Base_to_number($){
	$_=shift;
	return 0 if $_ eq "A";
	return 1 if $_ eq "G";
	return 2 if $_ eq "C";
	return 3 if $_ eq "T";
	return 4;
}


my $line;
my @bg_count;

open BG, $seq_bg_file or die "Can't open seq_bg file";
while(<BG>){
	chomp;
	$line=[];
	@{$line}=split /\t/, $_;
	push @bg_count, $line;
}

my $channel;
my $cos;
my $pos;
open SEQ, $seq_file or die "Can't open seq file";
open OUT, $out_file or die "Can't open out file";
while(<SEQ>){
	chomp;
	$line=[];
	@{$line}=split //, $_;
	$cos=0;
	foreach $pos (300..700){
		$channel=Base_to_number(${$line}[$pos]);
		$cos=$cos+${$bg_count[$channel]}[$pos];
	}
	$cos=$cos/401.0;
	print OUT $cos,"\n";
}
