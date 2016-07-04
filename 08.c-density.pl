#!/usr/bin/perl -w

use strict;

die "Usage:perl $0 <CHG><CHH><CpG><bedGraph.gz>\n\n" unless @ARGV == 4;

my $Mc = ();
open IN,"gzip -dc $ARGV[0] |" or die $!;
my %hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	next if ($stat =~ /[xhz]/);
	${$hash{join "\t",($chr,$pos)}}{$stat} = 1;
	$Mc++;
}

open IN,"gzip -dc $ARGV[1] |" or die $!;
%hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	next if ($stat =~ /[xhz]/);
	${$hash{join "\t",($chr,$pos)}}{$stat} = 1;
	$Mc++;
}

open IN,"gzip -dc $ARGV[2] |" or die $!;
%hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$stat) = (split "\t",$_)[-3,-2,-1];
	next if ($stat =~ /[xhz]/);
	${$hash{join "\t",($chr,$pos)}}{$stat} = 1;
	$Mc++;
}

my (%CHG,%CpG,%CHH);
open IN,"gzip -dc $ARGV[3]|" or die $!;
<IN>;
while(<IN>){
	chomp;
	my ($chr,$pos,$num) = (split "\t",$_)[0,2,3];
	if (exists ${$hash{join "\t",($chr,$pos)}}{"X"}){
		$CHG{$num} += 1;
	}
	if (exists ${$hash{join "\t",($chr,$pos)}}{"H"}){
		$CHH{$num} += 1;
	}
	if (exists ${$hash{join "\t",($chr,$pos)}}{"Z"}){
		$CpG{$num} += 1;
	}
}

for my $n(0..100){
	print "$n\t";
	if (exists $CHG{$n}){
		print "$CHG{$n}\t";
	}
	else {
		print "0\t";
	}
	if(exists $CHH{$n}){
		print "$CHH{$n}\t";
	}
	else{
		print "0\t";
	}
	if(exists $CpG{$n}){
		print "$CpG{$n}\n";
	}
	else{
		print "0\n";
	}
}