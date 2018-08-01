#EMu parser
use warnings;
use strict;

my $dir="/home/ted/EMu/";
my $out="XRXS_mut.txt";
my $input="XRXS.txt";
open OUT,">$dir$out" or die "Can't open output $dir$out";
my @line;
my $mut;

open FH, "$dir$input" or die "Can't open input $dir$input";

while(<FH>){
	chomp;
	@line=(split "\t", $_);
	$mut=(join ">", @line[(-2,-1)]);
	print OUT "@line[(0,1,2)]\t$mut\n";
}