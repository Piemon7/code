#!/usr/bin/perl -w
#piemon^2016.6.27
use strict;

die "Usage:perl $0 <cytosine.gz> <the min coverage> <outdir>\n\n" unless @ARGV == 3;

my $min_cover = $ARGV[1] or die $!;
if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else {
	open IN,$ARGV[0] or die $!;
}

my %hash = ();
my %density = ();
my %count = ();
while(<IN>){
	chomp;
	my ($A,$U,$type) = (split "\t",$_)[3,4,5];
	next if (($A+$U) < $min_cover);
	my $per = int($A/($A+$U)*100);
	$density{$per} = 1;
	${$hash{$type}}{$per} += 1;
	$count{$type} += 1;
}


open OUT,'>',"$ARGV[2]/08.out" or die $!;
for my $n(sort {$a<=>$b} keys %density){
#for my $n(0..100){
	print OUT "$n\t";
	if (exists ${$hash{"CHG"}}{$n}){
		print OUT ${$hash{'CHG'}}{$n}/$count{'CHG'},"\t";
	}
	else {
		print  OUT "0\t";
	}
	if(exists ${$hash{"CHH"}}{$n}){
		print OUT  ${$hash{'CHH'}}{$n}/$count{'CHH'},"\t";
	}
	else{
		print  OUT "0\t";
	}
	if(exists ${$hash{"CG"}}{$n}){
		print  OUT ${$hash{'CG'}}{$n}/$count{'CG'},"\n";
	}
	else{
		print  OUT "0\n";
	}
}

my $out = ();
my $pos = ();
$pos = $ARGV[2];

$out .= "setwd(\"$pos\")\n";
$out .= "data<-read.table(\"08.out\")\n";
$out .= "x<-data[,1]\n";
$out .= "y1<-data[,2]\n";
$out .= "y2<-data[,3]\n";
$out .= "y3<-data[,4]\n";
$out .= "png('08.png')\n";
$out .= "plot(data.frame(x,y1),type=\"l\",lwd=2,xlab=\"Methylation Level(%)\",ylab=\"Density\",font=2,col=\"orange\",main=\"sample\",xlim=c(0,100),ylim=c(0,1))\n";
$out .= "lines(data.frame(x,y2),type=\"l\",lwd=2,ylab=\"Density\",font=2,col=\"blue\",main=\"sample\")\n";
$out .= "lines(data.frame(x,y3),type=\"l\",lwd=2,ylab=\"Density\",font=2,col=\"red\",main=\"sample\")\n";
$out .= "legend(\"topright\",c(\"CHG\",\"CHH\",\"CpG\"),lty=c(1),col=c(\"orange\",\"blue\",\"red\"))\n";
$out .= "dev.off()\n";
$out .= "pdf('08.pdf')\n";
$out .= "plot(data.frame(x,y1),type=\"l\",lwd=2,xlab=\"Methylation Level(%)\",ylab=\"Density\",font=2,col=\"orange\",main=\"sample\",xlim=c(0,100),ylim=c(0,1))\n";
$out .= "lines(data.frame(x,y2),type=\"l\",lwd=2,ylab=\"Density\",font=2,col=\"blue\",main=\"sample\")\n";
$out .= "lines(data.frame(x,y3),type=\"l\",lwd=2,ylab=\"Density\",font=2,col=\"red\",main=\"sample\")\n";
$out .= "legend(\"topright\",c(\"CHG\",\"CHH\",\"CpG\"),lty=c(1),col=c(\"orange\",\"blue\",\"red\"))\n";
$out .= "dev.off()\n";

open OUT,'>',"$ARGV[2]/R_08.txt" or die $!;
print OUT $out;

system("Rscript $ARGV[2]/R_08.txt");
