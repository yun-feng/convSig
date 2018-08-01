use warnings;
use strict;

my $seq="/home/ted/conv/deep/multi/seq.txt";

open SEQ, $seq or die "Can't open output file $seq";

my @channel=("A","C","G","T");

sub BaseToValue($$){
	my $array=shift;
	my $base=shift;
	my $temp=[];
	for(@{$array}){
		if ($_ eq $base){push @{$temp}, 1;}
		else {push @{$temp}, 0;}
	}
	return $temp;
}



my $out_file=">/home/ted/conv/deep/multi/channel.txt";

open OUT, $out_file or die "Can't open output file $out_file";

my $line=[];

my @seq;

while(<SEQ>){
	chomp;
	@seq=split "", $_;
	$line = \@seq;
	for(@channel){
		print OUT (join "\t", @{BaseToValue($line,$_)});
		print OUT "\n";
	}
	
	

}