#!/usr/bin/perl -w 
#piemon^2016.6.7
use strict;

die "Usage:perl $0 <CX_report.txt> <02.out> <outdir>\n\n" unless @ARGV == 3;
my %hash = ();
my %all = ();

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}
while(<IN>){
	chomp;
	my ($yes,$no,$type) = (split "\t",$_)[3,4,5];
	my $depth = $yes + $no;
	$hash{$type}{$depth}++;
	$all{$type}++;
	$hash{'C'}{$depth}++;
	$all{'C'}++;
}

open OUT,'>',$ARGV[1] or die $!;
for (0..50){
	print OUT "$_\t";
}
print OUT "\n";

for my $ktype(sort keys %all){ ## the order is C/CG/CHG/CHH.
	#print OUT $ktype,"\t";
	my $total = ();
	for (keys %{$hash{$ktype}}){
		$total += $hash{$ktype}{$_};
	}
	for my $num(0..50){
		print OUT $total/$all{$ktype}*100,"\t";
		$total -= $hash{$ktype}{$num};
	}
	print OUT "\n";
}

## R draw picture.
my $pwd = $ARGV[2];
chomp $pwd;
my $outfile = "$ARGV[2]/R_02.txt";
open OUT,'>',$outfile or die $!;

print OUT "setwd(\"$pwd\")\n";
print OUT "data<-read\.table(\"$ARGV[1]\")\n";
print OUT "name<-data[1]\n";
print OUT "x<-t(data[1,])\n";
print OUT "y1<-t(data[2,])\n";
print OUT "y2<-t(data[3,])\n";
print OUT "y3<-t(data[4,])\n";
print OUT "y4<-t(data[5,])\n";
## pdf 
print OUT "pdf(\"02\.pdf\")\n";
print OUT "plot(x,y1,type=\"l\",lwd=2,xlab=\"Depth\",ylab=\"Percentage(%)\",ylim=c(0,100),font=2,col=\"green\",main=\"sample\")\n";
print OUT "lines(x,y2,font=2,col=\"red\")\n";
print OUT "lines(x,y3,font=2,col=\"orange\")\n";
print OUT "lines(x,y4,font=2,col=\"blue\")\n";
#print OUT "legend(\"topright\",c(\"C\",\"CpG\",\"CHG\",\"CHH\"),pch=c(15),col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "legend(\"topright\",c(\"C\",\"CpG\",\"CHG\",\"CHH\"),lty=c(1),col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";
## png
=pub
print OUT "png(\"02\.png\")\n";
print OUT "plot(x,y1,type=\"l\",lwd=2,xlab=\"Depth\",ylab=\"Percentage(%)\",ylim=c(0,100),font=2,col=\"green\",main=\"sample\")\n";
print OUT "lines(x,y2,font=2,col=\"red\")\n";
print OUT "lines(x,y3,font=2,col=\"orange\")\n";
print OUT "lines(x,y4,font=2,col=\"blue\")\n";
print OUT "legend(\"topright\",c(\"C\",\"CpG\",\"CHG\",\"CHH\"),pch=c(15),col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "legend(\"topright\",c(\"C\",\"CpG\",\"CHG\",\"CHH\"),lty=c(1),col=c(\"green\",\"red\",\"orange\",\"blue\"))\n";
print OUT "dev.off()\n";
=cut

system("R < $outfile --vanilla");