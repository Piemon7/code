#!/usr/bin/perl -w 

use strict;

die "Usage:perl $0 <CHG><CHH><CpG><depth.txt><genome><outfile>\n\n" unless @ARGV == 6;

my ($CG,$CHG,$CHH,$CN);
open IN,$ARGV[4] or die $!;
$/ = ">";<IN>;
my ($A,$C,$G,$T);
while(<IN>){
	chomp;
	$_ =~ s/^(.*)\n//;
	print "$1\t";
	$_ =~ s/\s+//g;
	my $length = length ($_);
	print "$length\n";
	$_ =~ tr/agctn/AGCTN/;
	$A     += $_ =~ tr/A//;
	$G     += $_ =~ tr/G//;
	$C     += $_ =~ tr/C//;
	$T     += $_ =~ tr/T//;
	my $line = $_;
	## C count
	for my $n(0..$length-3){
		my $str = substr ($line,$n,3);
		$CHG   += $str =~ s/(C[ATC]G)/$1/g;
		$CHH   += $str =~ s/(C[ACT][ACT])/$1/g;
		$CG    += $str =~ s/(CG[AGCT])/$1/g;
		$CN    += $str =~ s/(CN[AGCTN])/$1/g;
		$CN    += $str =~ s/(C[AGCT]N)/$1/g;
	}
	## G count (reverse C count)
	$line = reverse $line;
	$line =~ tr/ATCG/TAGC/;
	for my $m(0..$length-3){
		my $str = substr ($line,$m,3);
		$CHG   += $str =~ s/(C[ATC]G)/$1/g;
		$CHH  += $str =~ s/(C[ACT][ACT])/$1/g;
		$CG    += $str =~ s/(CG[AGCT])/$1/g;
	}
}
$/ = "\n";
print "whole genome C bases count is done!\n";

my @array = (0,1,2);
my @names = ("CHG","CHH","CpG");
my @values = ($CHG,$CHH,$CG,$CHG+$CHH+$CG);
my %values = ("CHG" => $CHG,"CHH" =>$CHH,"CpG"=>$CG);
my %Hash = ();

for my $i(@array){
	if($ARGV[$i] =~ /.gz$/){
		open IN,"gzip -dc $ARGV[$i]|" or die $!;
	}
	else {
		open IN,$ARGV[$i] or die $!;
	}

	my %hashCHG = ();
	while(<IN>){
		my $k = join "\t",(split "\t",$_)[2,3];
		$hashCHG{$k} = 1;
	}

	open IN,$ARGV[3] or die $!;
	my %count = ();
	while(<IN>){
		chomp;
		my @line = split;
		if (exists $hashCHG{join "\t",@line[0..1]}){
			$count{$line[2]} += 1;
		}
	}


	my $total = ();
	for (sort {$a<=>$b} keys %count){
		$total += $count{$_};
	}

	#print OUT $names[$i],"\t";
	for my $c(1..50){
		my $bases = $total;
	for (1..$c){
		$bases -= $count{$_};
	}
	#print OUT "$bases ";
	push @{$Hash{$names[$i]}},$bases;
	}
	#print OUT "\n";
}


open OUT,'>',$ARGV[5] or die $!;
print OUT "x\t";
for(0..50){
	print OUT "$_ ";
}
print OUT "\n";

for my $k(sort keys %Hash){
	print OUT "$k\t100 ";
	for my $v(@{$Hash{$k}}){
		print OUT $v/$values{$k}*100," ";
	}
	print OUT "\n";
}

print OUT "C\t100 ";
my @k = keys %Hash;
for my $x(0..49){
	my $sum = 0;
	for my $k(@k){
		$sum += ${$Hash{$k}}[$x];
	}
	print OUT $sum/$values[3]*100," ";
}
print OUT "\n";