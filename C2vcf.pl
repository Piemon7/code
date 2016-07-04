#!/usr/bin/perl -w 

use strict;

die "Uasge:<genome><CHG><CHH><CpG><outdir>\n\n" unless @ARGV == 5;


my %hash = ();
if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}
$/ = ">";
<IN>;
while(<IN>){
	chomp;
	/^(.*)\n/;
	my $head = $1;
	s/^$head//;
	s/\s//g;
	$hash{$head} = $_;
}
$/ = "\n";

if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
open OUT,'>',"$ARGV[4]/CHG.vcf" or die $!;
while(<IN>){
	chomp;
	next if (/[hzx]$/);
	my ($chr,$pos)=(split "\t",$_)[2,3];
	my $base = substr($hash{$chr},$pos-1,1);
	print OUT "$chr\t$pos\t\.\t$base\n";
}

if ($ARGV[2] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[2]|" or die $!;
}
else{
	open IN,$ARGV[2] or die $!;
}
open OUT,'>',"$ARGV[4]/CHH.vcf" or die $!;
while(<IN>){
	chomp;
	next if (/[hzx]$/);
	my ($chr,$pos)=(split "\t",$_)[2,3];
	my $base = substr($hash{$chr},$pos-1,1);
	print OUT "$chr\t$pos\t\.\t$base\n";
}

if ($ARGV[3] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[3]|" or die $!;
}
else{
	open IN,$ARGV[3] or die $!;
}
open OUT,'>',"$ARGV[4]/CpG.vcf" or die $!;
while(<IN>){
	chomp;
	next if (/[hzx]$/);
	my ($chr,$pos)=(split "\t",$_)[2,3];
	my $base = substr($hash{$chr},$pos-1,1);
	print OUT "$chr\t$pos\t\.\t$base\n";
}