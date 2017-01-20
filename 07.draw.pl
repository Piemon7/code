#!/usr/bin/perl -w 
#piemon^2016.6.28
use strict;
use File::Basename;

die "usage : perl $0 <07.out> <outdir> \n\n" unless @ARGV == 2;
my $outdir = $ARGV[1];
my $file = basename $ARGV[0];
$file .= ".R";

my $commond = ();
$commond .= "setwd('$outdir')\n";
$commond .= "data<-read.table('$ARGV[0]')\n";
$commond .= "col<-as.matrix(data[2])\n";
$commond .= "name<-as.matrix(data[1])\n";
$commond .= "png('07.png')\n";
$commond .= "counts<-as.matrix(data[3:5])\n";
$commond .= "par(las=1)\n";
#par(mar=c(5,5,10,10))
$commond .= "barplot(counts,width=0.05,ylim=c(0,100),col=col,ylab=\"Percentage (%)\",names.arg=c('sample1','sample2',NA))\n";
$commond .= "legend('right',name,bty='o',pch=c(15),col=col,cex=2)\n";
$commond .= "dev.off()\n";

open OUT,'>',"$outdir/$file" or die $!;
print OUT $commond;

system("Rscript $outdir/$file");
