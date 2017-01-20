#!/usr/bin/perl -w
#piemon^2016.6.27.
use strict;
die "Usage:perl $0 <gff> <genome.fasta> <bismark.bs/05.out.1> <12.out> <outdir>\n\n" unless @ARGV == 5;

$/ = "\tmRNA\t";
if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

my %position = ();
my %chromosome = ();
<IN>;
while(<IN>){
	my ($upstream,$downstream,$chr);
	chomp;
	s/##sequence(.*\n)+//;
	my @line = split "\n",$_;
	next if ((scalar @line) < 5);
	$line[0] =~ m/(\d+)\t(\d+)\t/;
	$upstream = $1;
	$downstream = $2;
	$line[1] =~ /^(\S+)\t/;
	$chr = $1;
	$chromosome{$chr} = 1;
	push @{$position{$chr}},(join "\t",($upstream,$downstream));
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


my (@upstream,@downstream,$number);
for(keys %position){
	$number += scalar @{$position{$_}};
}

for my $chrom(@chromosome){

	if ($ARGV[2] =~ /\.gz$/){
		open IN,"gzip -dc $ARGV[2]|" or die $!;
	}
	else{
		open IN,$ARGV[2] or die $!;
	}
	my %hash = ();
	while(<IN>){
		my @line = split;#NC_003070.9     3       +       0       0       CHH     CTA     intergenic 
		next if (($line[3]+$line[4]) == 0);
		$hash{$line[0]}{$line[1]} = $line[3]/($line[3]+$line[4])*100;
	}
	print "format.bs $chrom part info have been read !!! \n";

	for (@{$position{$chrom}}){
		my (@upstream_tmp,@downstream_tmp);
		my @line = split;
		@upstream_tmp = &level_count($chrom,$line[0]);
		@downstream_tmp = &level_count($chrom,$line[1]);

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

	#$number += scalar @{$position{$chrom}};

	sub level_count {

		my ($chr,$pos) = @_;
		my $bin = 100;
		my $start = ();
		if ($pos < 2000){
			$bin = int($pos/20);
			$start = 1;
		}
		elsif(($pos+2000) > $fasta_length{$chr}){
			$bin = int(($fasta_length{$chr}-$pos)/20);
			$start = $pos;
		}
		else{
			$start = $pos - 2000;
		}

		my @array = ();
		for my $count(0..39){
			my $sum = 0;
			for my $pos(($start+$count*$bin)..($start+($count+1)*$bin)){
				if (exists $hash{$chr}{$pos}){
					$array[$count] += $hash{$chr}{$pos};
					$sum ++;
				}
			}

			if ($sum == 0){
				$array[$count] = 0;
			}
			else{
				$array[$count] = $array[$count]/$sum;
			}
		}
		return @array;
	}

	}
}

open OUT,'>',"$ARGV[4]/$ARGV[3]" or die $!;

for(-40..-1){
	print OUT "$_\t";
}
print OUT "\n";
for(@upstream){
	$_ = $_/$number;
	print OUT "$_\t";
}
print OUT "\n";

for(1..40){
	print OUT "$_\t";
}
print OUT "\n";
for(@downstream){
	$_ = $_/$number;
	print OUT "$_\t";
}
print OUT "\n";

my $out = ();
my $pos = $ARGV[4];

$out .= "setwd(\"$pos\")\n";
$out .= "data<-read.table(\"12.out\")\n";
$out .= "x1<-data[1,]\n";
$out .= "y1<-data[2,]\n";
$out .= "x2<-data[3,]\n";
$out .= "y2<-data[4,]\n";
## png 
$out .= "png('12.png')\n";
$out .= "plot(t(x1),t(y1),type=\"l\",xaxt=\"n\",lwd=2,xlim=c(-50,50),xlab=\"Position\",ylab=\"Methylation Level(%)\",font=2,col=\"red\",main=\"sample\")\n";
$out .= "lines(t(x2),t(y2),type=\"l\",lwd=2,font=2,col=\"blue\")\n";
$out .= "axis(1,at=c(-40,-20,0,20,40),labels=c(\"-2k\",\"TSS\",\"2k......-2k\",\"TTS\",\"2k\"))\n";
$out .= "x1<-c(-20,-20,-20,-20,-20)\n";
$out .= "x2<-c(20,20,20,20,20)\n";
$out .= "z<-c(1,2,3,4,5)\n";
$out .= "lines(x1,z,type=\"l\",col=\"blue\")\n";
$out .= "lines(x2,z,type=\"l\",col=\"blue\")\n";
$out .= "dev.off()\n";
## pdf
$out .= "pdf('12.pdf')\n";
$out .= "plot(t(x1),t(y1),type=\"l\",xaxt=\"n\",lwd=2,xlim=c(-50,50),xlab=\"Position\",ylab=\"Methylation Level(%)\",font=2,col=\"red\",main=\"sample\")\n";
$out .= "lines(t(x2),t(y2),type=\"l\",lwd=2,font=2,col=\"blue\")\n";
$out .= "axis(1,at=c(-40,-20,0,20,40),labels=c(\"-2k\",\"TSS\",\"2k......-2k\",\"TTS\",\"2k\"))\n";
$out .= "x1<-c(-20,-20,-20,-20,-20)\n";
$out .= "x2<-c(20,20,20,20,20)\n";
$out .= "z<-c(1,2,3,4,5)\n";
$out .= "lines(x1,z,type=\"l\",col=\"blue\")\n";
$out .= "lines(x2,z,type=\"l\",col=\"blue\")\n";
$out .= "dev.off()\n";

open OUT,'>',"$pos/R_12.txt" or die $!;
print OUT $out;

system("Rscript $pos/R_12.txt");
