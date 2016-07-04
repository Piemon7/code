#!/usr/bin/perl -w 
#piemon^2016.6.27
use strict;

die "Usage:perl $0 <DMR.txt> <outdir> \n\n" unless @ARGV == 2;

my $outdir = $ARGV[1];

open OUT,'>',"$outdir/DMR_boxplot.txt" or die $!;
print OUT "type\tlen\n";
open IN,$ARGV[0] or die $!;
<IN>;
while(<IN>){
	chomp;
	my ($width,$diff) = (split "\t",$_)[3,5];
	#next if ($width eq 'width');
	if ($diff > 0){
		print OUT "1\t$width\n";
	}
	else{
		print OUT "2\t$width\n";
	}
}

my $outfile = "$outdir/DMR_boxplot.R";
open OUT,'>',$outfile or die $!;
my $out = ();
$out .= "setwd(\"$outdir\")\n";
$out .= "mydata=read.table('DMR_boxplot.txt',header=T)\n";
$out .= "png('DMR_boxplot.png',width=1000,height=800)\n";
$out .= "boxplot(len~type,data=mydata,varwidth=TRUE,col=c(\"red\",\"green\"),names=c(\"hyper\",\"hypo\"),notch=FALSE,outline=FALSE,main=\"sample1_vs_sample2\",ylab=\"DMR length\",ylim=c(0,120),cex=2,cex.lab=1.5,cex.axis=2,cex.mian=2,font=2)\n";
$out .= "dev.off()\n";
print OUT $out;

system("Rscript $outfile");
