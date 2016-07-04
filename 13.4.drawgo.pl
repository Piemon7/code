#!/usr/bin/perl -w 
#piemon^2016.6.28
use strict;
use File::Basename;

die "usage:perl $0 <gene1.out.go_enrichment> <outdir>\n\n" unless @ARGV == 2;
my $outdir = $ARGV[1];
open IN,$ARGV[0];
my $file = basename $ARGV[0];
my $png = $file.".png";
my $R_file = $file.".R";
$file .= '.temp';
open OUT,'>',"$outdir/$file" or die $!;
my $xlim = 50;
while(<IN>){
	next if (/^#/);
	s/biological\_process/red/g;
	s/cellular\_component/green/g;
	s/molecular\_function/orange/g;
	print OUT $_;
	$_ =~ /(\d+)$/;
	$xlim = $1 if ($xlim < $1);
}
close IN;
close OUT;
$xlim = (int($xlim/10)+2)*10;

my $R_commond = ();
$R_commond .= "setwd('$outdir')\n";
$R_commond .= "data<-read.table('$file',header=FALSE,sep=\"\t\")\n";
$R_commond .= "counts<-data[3]\n";
$R_commond .= "col<-data[1]\n";
$R_commond .= "png('$png',width=1200,height=1200)\n";
$R_commond .= "par(las=2)\npar(mar=c(5,20,10,10))\n";
$R_commond .= "barplot(as.matrix(counts),horiz=TRUE,cex.names=2,beside=TRUE,names.arg=as.matrix(data[2]),col=as.matrix(col),xlim=c(0,100))\n";
$R_commond .= "legend(\"topright\",bty='n',cex=2,c(\"biological_process\",\"cellular_component\",\"molecular_function\"),pch=c(15),col=c(\"red\",\"green\",\"orange\"))\n";
$R_commond .= "dev.off()\n";

open OUT,'>',"$outdir/$R_file" or die $!;
print OUT $R_commond;

system("Rscript $outdir/$R_file");

system("rm $outdir/$file");
system("rm $outdir/$R_file");
