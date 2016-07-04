#!/usr/bin/perl -w 
#piemon^2016.6.7.
use strict;

die "Usage:perl $0 <CX_report.txt> <03.out> <04.out> <outdir>\n\n" unless @ARGV == 4;

my %hash = ();
my %ML = ();
my %count = ();
if ($ARGV[0] =~ /\.gz/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

while(<IN>){
	chomp;
	my ($yes,$no,$type) = (split "\t",$_)[3,4,5];
	my $depth = $yes + $no;
	if ($depth == 0){
		$hash{$type}{'all'}++;
		$hash{'C'}{'all'}++;
	}
	else{
		$hash{$type}{'yes'}++;
		$hash{'C'}{'yes'}++;
		$hash{$type}{'all'}++;
		$hash{'C'}{'all'}++;
	}

	next if ($depth < 2); ## only consider the site depth >= 2;
	my $percent = $yes/$depth*100;
	$ML{$type} += $percent;
	$ML{'C'}   += $percent;
	$count{$type}++;
	$count{'C'}++;
}

open OUT,'>',"$ARGV[3]/$ARGV[1]" or die $!;
for (sort keys %hash){
	print OUT $_,"\t",$hash{$_}{'yes'}/$hash{$_}{'all'}*100,"\n";
}

open OUT,'>',"$ARGV[3]/$ARGV[2]" or die $!;
for(sort keys %ML){
	print OUT "$_\t",$ML{$_}/$count{$_},"\n";
}


# R draw picture.
my $pwd = $ARGV[3];
chomp $pwd;
my $outfile = "$ARGV[3]/R_034.txt";
open OUT,'>',$outfile or die $!;
## 03.out
print OUT "setwd(\"$pwd\")\n";
print OUT "data<-read.table(\"03.out\")\n";
print OUT "y<-as.vector(t(data[2]))\n";
print OUT "x<-as.vector(t(data[1]))\n";
## pdf
print OUT "pdf(\"03.pdf\")\n";
print OUT "barplot(y,names.arg=x,border=\"black\",space = c(1,1,1,1),ylim=c(0,85),xlab=\"sample\",ylab=\"coverage(%)\",beside=TRUE,col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";
## png
print OUT "png(\"03.png\")\n";
print OUT "barplot(y,names.arg=x,border=\"black\",space = c(1,1,1,1),ylim=c(0,85),xlab=\"sample\",ylab=\"coverage(%)\",beside=TRUE,col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";

## 04.out
print OUT "data<-read.table(\"04.out\")\n";
print OUT "y<-as.vector(t(data[2]))\n";
print OUT "x<-as.vector(t(data[1]))\n";
## pdf 
print OUT "pdf(\"04.pdf\")\n";
print OUT "barplot(y,names.arg=x,border=\"black\",space = c(1,1,1,1),ylim=c(0,35),xlab=\"sample\",ylab=\"percentage(%)\",beside=TRUE,col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";
## png
print OUT "png(\"04.png\")\n";
print OUT "barplot(y,names.arg=x,border=\"black\",space = c(1,1,1,1),ylim=c(0,35),xlab=\"sample\",ylab=\"percentage(%)\",beside=TRUE,col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";
system ("Rscript < $outfile");