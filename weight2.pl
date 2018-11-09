use warnings;
use strict;

my $seq_file="/data/ted/multi/eso/seq_open.txt";
my $seq_file2="/data/ted/multi/eso/seq_open2.txt";
my $seq_cos_file=">/data/ted/multi/eso/seq_open_cos_approx400.txt";

my $current=[];
my $comp=[];
my $count=0;
my $temp=0;
open SEQ_out,$seq_file2 or die "Can't open seq file";

#open OUT, $seq_count_file or die "Can't open out file";
open COS, $seq_cos_file or die "Can't open cos file";

my @seq;
my $frag=[];
open SEQ,$seq_file or die "Can't open seq file";
while(<SEQ>){
	chomp;
	$frag=[];
	@{$frag}=split //, $_;
	push @seq, $frag;
}




my $line=0;
my $randline;
my $c;
while(<SEQ_out>){
	chomp;
	@{$current}=split //, $_;
	$count=0.0;
	$c=10000;
	while($c>0){
		$temp=0;
		$randline=int(rand(scalar @seq));
		next, if $randline==$line;
		for(300..700){
			$temp++, if ${$seq[$randline]}[$_] eq ${$current}[$_];
		}
		$temp=$temp/401.0;
		$count=$count+$temp*$temp;
		$c--;
	}
	$count=$count/10000*((scalar @seq)-1)+1;
	print $line;
	print COS $count,"\n";
	$line++;
}
