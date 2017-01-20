#!/usr/bin/perl -w

#piemon^2016.6.1.
use strict;
die "Usage:perl $0 <gff><genome.fasta><format.bs><12.out.txt>\n\n" unless @ARGV == 4;

$/ = "\tmRNA\t";
if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

open OUT,'>',"test.txt" or die $!;
my %position = ();
my %chromosome = ();
<IN>;
while(<IN>){
	my ($upstream,$exon1,$intron1,$exonL,$downstream,$chr);
	chomp;
	s/##sequence(.*\n)+//;
	my @line = split "\n",$_;
	next if ((scalar @line) < 5);
	$line[0] =~ m/(\d+)\t(\d+)\t/;
	$upstream = $1;
	$downstream = $2;
	my %hash = ();
	for my $line(1..$#line-1){
		my ($chromosome,$type,$start,$end)= (split "\t",$line[$line])[0,2,3,4];
		$chr = $chromosome;
		$chromosome{$chr} = 1;
		if (!defined $type){
			print "$chr\t$upstream\n";
		}
		next if ($type ne 'exon');
		$hash{$start} = $end;
	}
	my @array = sort {$a<=>$b} keys %hash;
	next if ((scalar @array) < 4);
	$exon1 = shift @array;
	$exon1 = $hash{$exon1};
	$intron1 = shift @array;
	$exonL = pop @array;
	print OUT "$chr\t$upstream\t$exon1\t$intron1\t$exonL\t$downstream\n";
	push @{$position{$chr}},(join "\t",($upstream,$exon1,$intron1,$exonL,$downstream));
}
$/ = "\n";

print "gff info have been read !!! \n";

my @chromosome = sort keys %chromosome;

$/ = ">";
if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
my %fasta_length = ();
<IN>;
while(<IN>){
	chomp;
	/^(\S+)/;
	my $head = $1;
	s/$head//;
	s/\s//g;
	$fasta_length{$head} = length($_);
	print "$head\t$fasta_length{$head}\n";
}
$/ = "\n";


my (@upstream,@exon1,@intron1,@inter_exon,@inter_intron,@exonL,@downstream,$number);

for my $chrom(@chromosome){

	if ($ARGV[2] =~ /\.gz$/){
		open IN,"gzip -dc $ARGV[1]|" or die $!;
	}
	else{
		open IN,$ARGV[2] or die $!;
	}
	my %hash = ();
	while(<IN>){
		my @line = split;#NC_003070.9     3       +       0       0       CHH     CTA     intergenic 
		next if (($line[3]+$line[4]) == 0);
		$hash{$line[0]}{$line[1]} = $line[3]/($line[3]+$line[4]);
	}
	print "format.bs $chrom part info have been read !!! \n";

	for (@{$position{$chrom}}){
		my (@upstream_tmp,@exon1_tmp,@intron1_tmp,@inter_exon_tmp,@inter_intron_tmp,@exonL_tmp,@downstream_tmp);
		my @line = split;
		if ($line[0]<100000) {# the distance of upstream or downstream;default is 10k.
			@upstream_tmp = &level_count($chrom,1,$line[0],40);#the piece of element;default is 40.
		}
		else{
			@upstream_tmp = &level_count($chrom,$line[0]-100000,$line[0],40);
		}
		@exon1_tmp = &level_count($chrom,$line[0],$line[1],20);#the piece of element;default is 20.
		@intron1_tmp = &level_count($chrom,$line[1],$line[2],20);
		@exonL_tmp = &level_count($chrom,$line[3],$line[4],20);
		if(($line[4]+100000) > $fasta_length{$chrom}){
			@downstream_tmp = &level_count($chrom,$line[4],$fasta_length{$chrom},40);
		}
		else{
			@downstream_tmp = &level_count($chrom,$line[4],$line[4]+100000,40);##
		}
		for(0..39){
			if (!defined $upstream[$_]){
				$upstream[$_] = $upstream_tmp[$_];
			}
			else{
				$upstream[$_] += $upstream_tmp[$_];
			}
			if(!defined $downstream[$_]){
				$downstream[$_] = $downstream_tmp[$_];
			}
			else{
				$downstream[$_] += $downstream_tmp[$_];
			}
		}
		for(0..19){
			if(!defined $exon1[$_]){
				$exon1[$_] = $exon1_tmp[$_];
			}
			else{
				$exon1[$_] += $exon1_tmp[$_];
			}
			if(!defined $intron1[$_]){
				$intron1[$_] = $intron1_tmp[$_];
			}
			else{
				$intron1[$_] += $intron1_tmp[$_];
			}
			if(!defined $exonL[$_]){
				$exonL[$_] = $exonL_tmp[$_];
			}
			else{
				$exonL[$_] += $exonL_tmp[$_];
			}
		}

	}

	$number += scalar @{$position{$chrom}};

sub level_count {

	warn "number of \@\_ error!" if (@_ < 4);
	my ($chr,$start,$end,$piece) = @_;

	my $length = $end-$start+1;
	my $bin = int($length/$piece);

	print "$chr\t$start\tbin:$bin\n" if ($bin < 10);

	my @array = ();
	for my $count(0..$piece-2){
		for my $pos(($start+$count*$bin)..($start+($count+1)*$bin)){
				my $sum = 0;
				if (exists $hash{$chr}{$pos}){
					$array[$count] += $hash{$chr}{$pos};
					$sum ++;
				}
				if ($sum == 0){
					$array[$count] = 0;
				}
				else{
					$array[$count] = $array[$count]/$sum;
				}
		}
	}
	for my $pos(($start+($piece-1)*$bin)..$end){
		my $sum = 0;
		if (exists $hash{$chr}{$pos}){
			$array[$piece-1] += $hash{$chr}{$pos};
			$sum ++;
		}
		if($sum == 0){
			$array[$piece-1] = 0;
		}
		else{
			$array[$piece-1] = $array[$piece-1]/$sum;
		}
	}

	warn "error!\n" if ((scalar @array) != $piece);
	return @array;
}

}

open OUT,'>',$ARGV[3] or die $!;

for(1..140){
	print OUT "$_\t";
}
print OUT "\n";

for(@upstream,@exon1,@intron1,@exonL,@downstream){
	$_ = $_/$number;
	print OUT "$_\t";
}
print OUT "\n";