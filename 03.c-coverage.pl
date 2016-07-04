#!/usr/bin/perl -w 

use strict;

die "Usage:perl $0 <reference><CHG><CpG><CHH>\n\n" unless @ARGV==4;

my ($total,$CG,$CHG,$CHH,$CN);
open IN,$ARGV[0] or die $!;
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
	for my $n(0..$length-3){
		my $str = substr ($line,$n,3);
		$CHG   += $str =~ s/(C[ATC]G)/$1/g;
		$CHH   += $str =~ s/(C[ACT][ACT])/$1/g;
		$CG    += $str =~ s/(CG[AGCT])/$1/g;
		$CN    += $str =~ s/(CN[AGCTN])/$1/g;
		$CN    += $str =~ s/(C[AGCT]N)/$1/g;
	}

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

if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1] |" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
<IN>;
my %CHG = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$x) = (split "\t",$_)[-3,-2,-1];
	if ($x ne "x" or $x ne "X"){
		$CHG{join "\t",($chr,$pos)} = 1;
	}
}


if ($ARGV[2] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[2] |" or die $!;
}
else{
	open IN,$ARGV[2] or die $!;
}
<IN>;
my %CpG = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$z) = (split "\t",$_)[-3,-2,-1];
	if ($z ne "z" or $z ne "Z"){
		$CpG{join "\t",($chr,$pos)} = 1;
		#print join "\t",($chr,$pos),"\n";
	}
}

if ($ARGV[3] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[3] |" or die $!;
}
else{
	open IN,$ARGV[3] or die $!;
}
<IN>;
my %CHH = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$h) = (split "\t",$_)[-3,-2,-1];
	if ($h eq "h" or $h eq "H"){
		$CHH{join "\t",($chr,$pos)} = 1;
	}
}

print "A\t$A\nG\t$G\nC\t$C\nT\t$T\n";
print "CG\t",(scalar keys %CpG),"\t",$CG,"\t",(scalar keys %CpG)*100/($CG),"\n";
print "CHG\t",(scalar keys %CHG),"\t",$CHG,"\t",(scalar keys %CHG)*100/($CHG),"\n";
print "CHH\t",(scalar keys %CHH),"\t",$CHH,"\t",(scalar keys %CHH)*100/($CHH),"\n";
#print "CN\t",$CN,"\n";
print "C\t",(scalar keys %CpG)+(scalar keys %CHG)+(scalar keys %CHH),"\t",($CG+$CHG+$CHH),"\t",((scalar keys %CpG)+(scalar keys %CHG)+(scalar keys %CHH))/($CG+$CHG+$CHH)*100,"\n";
