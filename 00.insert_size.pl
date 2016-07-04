#!/usr/bin/perl -w 
#piemon^2016.6.27
use strict;

die "Usage:perl $0 <outdir> <insert_size 1> <insert_size 2> .. \n\n" unless @ARGV > 1;

my $outdir = $ARGV[0];

my $xmax = 0;
my $xmin = 100;
my $ymax = 0;
for (1..$#ARGV){
	open IN,$ARGV[$_] or die $!;
	open OUT,'>',"$outdir/$_.insert_size.txt" or die $!;
	while(<IN>){
		my $line = $_;
		if ($line =~ /^(\d+)\s(\d+)\n$/){
			$xmin = $1 if ($1 < $xmin);
			$xmax = $1 if ($1 > $xmax);
			$ymax = $2 if ($2 > $ymax);
			print OUT $line;
		}
	}
}

$xmax = int(1.2*$xmax);
$xmin = int(0.8*$xmin);
$ymax = int(1.2*$ymax);
my %color = (
	'1' => 'red',
	'2' => 'green',
	'3' => 'blue',
	'4' => 'yellow',
	'5' => 'orange',
	'6' => 'purple',
	);

open OUT,'>',"$outdir/00.insert_size.R" or die $!;
print OUT "setwd(\'$outdir\')\n";
my $color = ();
my $sample = ();

## pdf 
for (1..$#ARGV){
	$sample .= "\"sample"."$_"."\"".",";
	$color .= "\"$color{$_}\"".",";
	if ($_ == 1){
		print OUT "data<-read.table(\'$_.insert_size.txt\')\n";
		print OUT "pdf(\'insert_size.pdf\')\n";
		print OUT "plot(data.frame(data[1],data[2]),type=\"l\",lwd=2,xlim=c($xmin,$xmax),ylim=c(0,$ymax),xlab=\"Length of Insert\",ylab=\"Count\",font=2,col=\"$color{$_}\")\n";
	}
	else{
		print OUT "data<-read.table(\'$_.insert_size.txt\')\n";
		print OUT "lines(data.frame(data[1],data[2]),type=\"l\",lwd=2,font=2,col=\"$color{$_}\")\n";
	}
}

$sample =~ s/,$//;
$color =~ s/,$//;

print OUT "legend(\"topright\",c($sample),lty=c(1),col=c($color))\n";
print OUT "dev.off()\n";

## png 
for (1..$#ARGV){
	if ($_ == 1){
		print OUT "data<-read.table(\'$_.insert_size.txt\')\n";
		print OUT "png(\'insert_size.png\')\n";
		print OUT "plot(data.frame(data[1],data[2]),type=\"l\",lwd=2,xlim=c($xmin,$xmax),ylim=c(0,$ymax),xlab=\"Length of Insert\",ylab=\"Count\",font=2,col=\"$color{$_}\")\n";
	}
	else{
		print OUT "data<-read.table(\'$_.insert_size.txt\')\n";
		print OUT "lines(data.frame(data[1],data[2]),type=\"l\",lwd=2,font=2,col=\"$color{$_}\")\n";
	}
}

print OUT "legend(\"topright\",c($sample),lty=c(1),col=c($color))\n";
print OUT "dev.off()\n";

## draw the picture
system('Rscript 00.insert_size.R');
