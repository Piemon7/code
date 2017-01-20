#!/usr/bin/perl -w 

#piemon^2016.5.31.
use strict;

die "Usage:perl <05.out.1> <05.final.xls/05.out.2> <09.out> <outdir>\n\n" unless @ARGV == 4;

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}
my %hash = ();
my %info = ();

my %cylevel = ();

while(<IN>){
	chomp;
	my ($cyes,$cno,$type,$info) = (split "\t",$_)[3,4,5,7];
	my @info = split " ",$info;	
	my $cypercent = ();
	if (($cyes+$cno) == 0){
		for my $inf(@info){
			$info{$inf} = 1;
			$hash{$type}{$inf}->{'sum'} += 1;
		}
	}
	else{
		$cypercent = 100 * $cyes/($cyes+$cno);
		for my $inf(@info){
			$info{$inf} = 1;
			$hash{$type}{$inf}->{'sum'} += 1;
			$hash{$type}{$inf}->{'C'} += 1;
			$cylevel{$type}{$inf}->{'sum'} += 1;
			$cylevel{$type}{$inf}->{'C'} += $cypercent;
		}

	}
}

open OUT,'>',"$ARGV[3]/$ARGV[1]" or die $!;
open OUTPUT,'>',"$ARGV[3]/$ARGV[2]" or die $!;
print OUT " \tCG\tCHG\tCHH\tC\n";
print OUTPUT " \tCG\tCHG\tCHH\tC\n";

my @Ctype = ('CG','CHG','CHH');

for my $inf(sort keys %info){
	my ($sum,$C,$sump,$Cp);
	print OUT $inf,"\t";
	print OUTPUT $inf,"\t";
	for my $type(@Ctype){
	print OUT 100*$hash{$type}{$inf}->{'C'}/$hash{$type}{$inf}->{'sum'},"\t";
	print OUTPUT $cylevel{$type}{$inf}->{'C'}/$cylevel{$type}{$inf}->{'sum'},"\t";
	$sum += $hash{$type}{$inf}->{'sum'};
	$sump +=  $cylevel{$type}{$inf}->{'sum'};
	$C += $hash{$type}{$inf}->{'C'};
	$Cp +=  $cylevel{$type}{$inf}->{'C'};
	}
	print OUT $C/$sum*100,"\n";
	print OUTPUT $Cp/$sump,"\n";
}

## 09
open OUT,'>',"$ARGV[3]/R_59.txt" or die $!;
my $pos = $ARGV[3];
my $out = ();

$out .= "setwd(\"$pos\")\n";
$out .= "data<-read.table(\"09.out\",header=TRUE)\ncol=colnames(data)\n";
$out .= "row=rownames(data)\n";
$out .= "pdf(\"09.pdf\",width=800,height=500)\n";
$out .= "barplot(t(as.matrix(data)),beside=TRUE,xlab=\"sample\",ylab=\"Methylation Level(%)\",col=c(\"red\",\"yellow\",\"blue\",\"green\"),ylim=c(0,60))\n";
$out .= "legend(\"topright\",c(\"CG\",\"CHG\",\"CHH\",\"C\"),lty=c(1),col=c(\"red\",\"yellow\",\"blue\",\"green\"))\n";
$out .= "dev.off()\n";
$out .= "png(\"09.png\",width=800,height=500)\n";
$out .= "barplot(t(as.matrix(data)),beside=TRUE,xlab=\"sample\",ylab=\"Methylation Level(%)\",col=c(\"red\",\"yellow\",\"blue\",\"green\"),ylim=c(0,60))\n";
$out .= "legend(\"topright\",c(\"CG\",\"CHG\",\"CHH\",\"C\"),lty=c(1),col=c(\"red\",\"yellow\",\"blue\",\"green\"))\n";
$out .= "dev.off()\n";

## 05 
$out .= "setwd(\"$pos\")\n";
$out .= "data<-read.table(\"05.out.2\",header=TRUE)\ncol=colnames(data)\n";
$out .= "row=rownames(data)\n";
$out .= "pdf(\"05.pdf\",width=800,height=500)\n";
$out .= "barplot(t(as.matrix(data)),beside=TRUE,xlab=\"sample\",ylab=\"Methylation Level(%)\",col=c(\"red\",\"yellow\",\"blue\",\"green\"),ylim=c(0,80))\n";
$out .= "legend(\"topright\",c(\"CG\",\"CHG\",\"CHH\",\"C\"),lty=c(1),col=c(\"red\",\"yellow\",\"blue\",\"green\"))\n";
$out .= "dev.off()\n";
$out .= "png(\"05.png\",width=800,height=500)\n";
$out .= "barplot(t(as.matrix(data)),beside=TRUE,xlab=\"sample\",ylab=\"Methylation Level(%)\",col=c(\"red\",\"yellow\",\"blue\",\"green\"),ylim=c(0,80))\n";
$out .= "legend(\"topright\",c(\"CG\",\"CHG\",\"CHH\",\"C\"),lty=c(1),col=c(\"red\",\"yellow\",\"blue\",\"green\"))\n";
$out .= "dev.off()\n";

print OUT $out;
system("Rscript $ARGV[3]/R_59.txt");