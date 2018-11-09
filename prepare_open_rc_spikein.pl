use warnings;
use strict;

my $fasta_file="/well/htseq/Genomes/BWA/GRCh37";
my $mut_file="/data/ted/COAD/esophagael-de.txt";
my $chrom_file ="/data/ted/multi/Ubiquitous-sorted.bed";

my @rc_table=("T","C","G","A","N");

sub RC($){
	my $old_array=shift;
	my $new_array=[];
	my $temp;
	while(@{$old_array}){
	$temp=pop @{$old_array};
	$temp=Base_to_number($temp);
	push @{$new_array},$rc_table[$temp];
	}
	return $new_array;
}

sub Base_to_number($){
	$_=shift;
	return 0 if $_ eq "A";
	return 1 if $_ eq "G";
	return 2 if $_ eq "C";
	return 3 if $_ eq "T";
	return 4;
}

open FASTA, $fasta_file or die "Can't open fasta file";

my @genome;
my $temp=[];
my $chr=0;


print "fasta","\n";
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

print "mut","\n";
my $sample=-1;
my %sample_hash;
my @sample_table;
my $sample_name;
my @chrom_table;
my @pos_table;
my @mut_type;
my @ori_type;
my ($chrom,$pos,$ori,$mut);
open MUT,$mut_file or die "Can't open mutation file $mut_file";
while(<MUT>){
	chomp;
	($sample_name,$chrom,$pos,$ori,$mut)=(split "\t", $_);
	 $mut=Base_to_number($mut);
	next, if Base_to_number($ori)>3;
	next, if $mut>3;
	next, if ($chrom eq 'X' or  $chrom eq 'Y');
	push @chrom_table, $chrom;
	push @pos_table, $pos;
	push @ori_type, $ori;
	push @mut_type, $mut;
	$sample_hash{$sample_name}=++$sample, if(not exists($sample_hash{$sample_name}));
	push @sample_table, $sample_name;
}

print "open","\n";
my $dir="/data/ted/multi/eso/";
open SEQ, ">${dir}seq_open_rc_si6.txt" or die "Can't open output file ${dir}seq.txt";
open SEQ_bg, ">${dir}seq_bg_open_rc_si6.txt" or die "Can't open output file ${dir}seq.txt";
open LABEL, ">${dir}label_open_rc_si6.txt" or die "Can't open output file ${dir}label.txt";
open SAM, ">${dir}sample_open_rc_si6.txt" or die "Can't open output file ${dir}sample.txt";
open MID, ">${dir}mid_open_rc_si6.txt" or die "Can't open output file ${dir}mid.txt";

my $ochr;
my $os;
my $oe;
my $genome_pos;
my $flag=0;
my $seq=[];
my $spikein_pos;
$chrom=shift @chrom_table;
$pos=shift @pos_table;
$ori=shift @ori_type;
$mut=shift @mut_type;
$sample_name=shift @sample_table;
open CHROM,$chrom_file or die "Can't open mutation file $chrom_file";
while(<CHROM>){
	chomp;
	($ochr,$os,$oe)=(split "\t", $_);
	$genome_pos=$os;
	next, if $genome_pos+499> scalar @{$genome[$ochr]};
	while($genome_pos<$oe){
		while($chrom and ($chrom<$ochr or ($chrom==$ochr and $pos<$genome_pos))){
			$chrom=shift @chrom_table;
			$pos=shift @pos_table;
			$ori=shift @ori_type;
			$mut=shift @mut_type;
			$sample_name=shift @sample_table;
		}
		if($chrom and ($chrom==$ochr and $pos==$genome_pos)){
			print "ERROR: in $chrom,$pos, $ori not ${$genome[$ochr]}[$genome_pos-1] \n", if ${$genome[$ochr]}[$genome_pos-1] ne $ori;
			print SAM $sample_hash{$sample_name}+1, "\n";
			@{$seq}=@{$genome[$ochr]}[($genome_pos-501)..($genome_pos+499)];
			if (Base_to_number($ori)<2){
				$seq=RC($seq);
				$mut=3-$mut;
			}
			$mut=2, if $mut>2;
			$spikein_pos=3;
			while($spikein_pos>0){
				${$seq}[248+$spikein_pos]="A", if(${$seq}[500] eq 'C' and $mut==0 and rand()>0.6);
				$spikein_pos--;
			}
			if ($ori eq 'A' or $ori eq 'T'){
				print MID "0\t1\n";
			}
			else{
				print MID "1\t0\n";
			}
			print ${$seq}[500];
			print LABEL $mut+1, "\n";
			print SEQ (join "",@{$seq}), "\n";
			$chrom=shift @chrom_table;
                        $pos=shift @pos_table;
                        $ori=shift @ori_type;
                        $mut=shift @mut_type;
                        $sample_name=shift @sample_table;
			$flag=1;
			next;
		}
		
		else{
			if($flag==0){
			print SEQ_bg (join "",@{$genome[$ochr]}[($genome_pos-500)..($genome_pos+500)]), "\n";}
			
		}
	$genome_pos++;
	$flag=0;
	last , if $genome_pos+499> scalar @{$genome[$ochr]};
	}
	print $ochr," ",$os,"\n";
}

open OUT, ">${dir}sample_hash_open_rc_si6.txt" or die "Can't open output file ${dir}XRXS.sample.txt";
for (keys %sample_hash){
	print OUT $_,"\t",$sample_hash{$_},"\n";
}

