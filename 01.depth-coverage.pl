#!/usr/bin/perl -w
#piemon^2016.6.7.
use strict;

die "Usage:perl $0 <depth.txt> <reference> <depth-coverage.out> <outdir>\n\n" unless @ARGV == 4;

open IN,$ARGV[0] or die $!;
my %hash = ();
while (<IN>) {
	chomp;
	my @line = split;
	$hash{$line[2]} += 1;
}

if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
my $length = ();
while(<IN>){
	next if (/^[>#]/);
	$length += $_ =~ tr/[AGCTagct]//;
}

print "whole genome effect length:\t",$length,"\n";

open OUT,'>',$ARGV[2] or die $!;
for (sort {$a<=>$b} keys %hash){
	print OUT "$_\t",$hash{$_}/$length*100,"\n";
}

## R draw picture
my $pwd = $ARGV[3];
chomp $pwd;
my $outfile = "$ARGV[3]/R_01.txt";
open OUT,'>',$outfile or die $!;
print OUT "setwd(\"$pwd\")\n";
print OUT "data<-read.table(\"$ARGV[2]\")\n";

# pdf format
print OUT "pdf(\"01.pdf\")\n";
print OUT "plot(data\.frame(data[1],data[2]),type=\"l\",lwd=2,xlab=\"Depth\",ylab=\"Percentage(\%)\",font=2,col=\"red\",main=\"sample\",xlim=c(0,80),ylim=c(0,15))\n";
print OUT "dev\.off()\n";

# png format
#print OUT "png(\"01.png\")\n";
#print OUT "plot(data\.frame(data[1],data[2]),type=\"l\",lwd=2,xlab=\"Depth\",ylab=\"Percentage(\%)\",font=2,col=\"red\",main=\"sample\",xlim=c(0,80),ylim=c(0,15))\n";
#print OUT "dev\.off()\n";

system("R <$outfile --vanilla");