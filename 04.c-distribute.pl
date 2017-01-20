#!/usr/bin/perl -w

use strict;

die "Usage:perl $0 <CHG><CHH><CpG><genome>\n\n" unless @ARGV == 4;

open IN,"gzip -dc $ARGV[0] |" or die $!;
my %hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	${$hash{join "\t",($chr,$pos)}}{$stat} += 1;
}

my $CHG = ();
for my $k(keys %hash){
	next if (not exists ${$hash{$k}}{"X"});
	${$hash{$k}}{"x"} = 0 if (not exists ${$hash{$k}}{"x"});
	next if ((${$hash{$k}}{"X"} + ${$hash{$k}}{"x"}) < 0 );
	$CHG += 100*${$hash{$k}}{"X"}/(${$hash{$k}}{"X"} + ${$hash{$k}}{"x"});
}

open IN,"gzip -dc $ARGV[1] |" or die $!;
%hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	${$hash{join "\t",($chr,$pos)}}{$stat} += 1;
}

my $CHH = ();
for my $k(keys %hash){
	next if (not exists ${$hash{$k}}{"H"});
	${$hash{$k}}{"h"} = 0 if (not exists ${$hash{$k}}{"h"});
	next if ((${$hash{$k}}{"H"} + ${$hash{$k}}{"h"}) < 0 );
	$CHH += 100*${$hash{$k}}{"H"}/(${$hash{$k}}{"H"} + ${$hash{$k}}{"h"});
}

open IN,"gzip -dc $ARGV[2] |" or die $!;
%hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	${$hash{join "\t",($chr,$pos)}}{$stat} += 1;
}

my $CpG = ();
for my $k(keys %hash){
	next if (not exists ${$hash{$k}}{"Z"});
	${$hash{$k}}{"z"} = 0 if (not exists ${$hash{$k}}{"z"});
	next if ((${$hash{$k}}{"Z"} + ${$hash{$k}}{"z"}) < 0 );
	$CpG += 100*${$hash{$k}}{"Z"}/(${$hash{$k}}{"Z"} + ${$hash{$k}}{"z"});
}

open IN,$ARGV[3] or die $!;
$/ = ">";<IN>;
my ($A,$C,$G,$T);
my ($CHGt,$CHHt,$CpGt);
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
		$CHGt   += $str =~ s/(C[ATC]G)/$1/g;
		$CHHt  += $str =~ s/(C[ACT][ACT])/$1/g;
		$CpGt    += $str =~ s/(CG[AGCT])/$1/g;
		#$CN    += $str =~ s/(CN[AGCTN])/$1/g;
		#$CN    += $str =~ s/(C[AGCT]N)/$1/g;
	}

	$line = reverse $line;
	$line =~ tr/ATCG/TAGC/;
	for my $m(0..$length-3){
		my $str = substr ($line,$m,3);
		$CHGt   += $str =~ s/(C[ATC]G)/$1/g;
		$CHHt  += $str =~ s/(C[ACT][ACT])/$1/g;
		$CpGt    += $str =~ s/(CG[AGCT])/$1/g;
	}

}
$/ = "\n";

print "$C\t$G\n";
print "CHG:\t",$CHG/($CHGt),"\n";
print "CHH:\t",$CHH/($CHHt),"\n";
print "CpG:\t",$CpG/($CpGt),"\n";
print "C:\t",($CHG+$CHH+$CpG)/($CHGt+$CHHt+$CpGt),"\n";