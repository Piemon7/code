#!/usr/bin/perl -w 
#piemon^2016.7.1
use strict;
use File::Basename;

die "Usage:perl $0 <cy..1> <cy..2> <DMR> <outdir>\n\n" unless @ARGV == 4;

my %dmr = ();
my $outdir = $ARGV[3];
open OUT,'>',"$outdir/cluster.txt" or die $!;
open IN,$ARGV[2] or die $!;
<IN>;
while(<IN>){
	chomp;
	my ($chr,$start,$end) = (split "\t",$_)[0,1,2];
	$chr =~ s/\"//g;

	if ($ARGV[0] =~ /\.gz$/){
		open FL,"gzip -dc $ARGV[0]|" or die $!;
	}
	else{
		open FL,$ARGV[0] or die $!;
	}
	my $dmr_total_ml = 0;
	while(<FL>){
		chomp;
		my ($chr_1,$pos,$yes,$no) = (split "\t",$_)[0,1,3,4];
		next if ($chr_1 ne $chr);
		my $sum = $yes + $no;
		next if ($sum == 0);
		if ($start < $pos and $pos < $end){
			$dmr_total_ml += $yes/$sum;
		}
	}
	#print OUT "$chr\t$start\t$end\t$dmr_total_ml";
	print OUT "$dmr_total_ml";

	if ($ARGV[1] =~ /\.gz$/){
		open FH,"gzip -dc $ARGV[1]|" or die $!;
	}
	else{
		open FH,$ARGV[1] or die $!;
	}

	my $dmr_total_ml_2 = 0;
	while(<FH>){
		chomp;
		my ($chr_2,$pos,$yes,$no) = (split "\t",$_)[0,1,3,4];
		next if ($chr_2 ne $chr);
		my $sum = $yes + $no;
		next if ($sum == 0);
		if ($start < $pos and $pos < $end){
			$dmr_total_ml_2 += $yes/$sum;
		}
	}
	print OUT "\t$dmr_total_ml_2\n";
}

my $commond_R = ();
$commond_R .= "setwd('$outdir')\n";
$commond_R .= "library(reshape2)\n";
$commond_R .= "library(ggplot2)\n";
$commond_R .= "df<-read.table('cluster.txt')\n"
$commond_R .= "k <- kmeans(df,3)\n";
$commond_R .= "dfc <- cbind(df, id=seq(nrow(df)), cluster=k\$cluster)\n";
$commond_R .= "dfc\$idsort <- dfc\$id[order(dfc\$cluster)]\n";
$commond_R .= "dfc\$idsort <- order(dfc\$idsort)\n";
$commond_R .= "dfm <- melt(dfc, id.vars=c(\"id\",\"idsort\"))\n";
$commond_R .= "library(RColorBrewer)\n";
#heatmap(as.matrix(df)[order(k$cluster),],Rowv=NA,scale="none",labRow=rownames(df),labCol = rownames(df))
$commond_R .= "heatmap(as.matrix(df)[order(k\$cluster),],Rowv=NA,scale=\"none\",labCol=c('sample1','sample2'),labRow=NA,cexCol=1,)\n";

open OUT,'>',"$outdir/cluster.R" or die $!;
print OUT $commond_R;

system("/usr/bin/Rscript $outdir/cluster.R");
