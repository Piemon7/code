#!/usr/bin/perl -w 
#piemon^2016.7.4
use strict;
use File::Basename;
die "usage:perl $0 <individual.txt> <vcf> <pca count> <outdir>
the max number of group is 17. it depends on the \@col.\n\n" unless @ARGV == 4;

my $vcftools = "/lustre/project/og04/shichunwei/biosoft/vcftools_0.1.13/bin/vcftools";
my $plink = "/lustre/project/og04/shichunwei/biosoft/plink/plink";
my $gcta = "/lustre/project/og04/shichunwei/biosoft/gcta/gcta64";
my %hash = ();
my $outdir = $ARGV[3];
my $name = basename $ARGV[1];
my $vcf_file = $ARGV[1];
my $run_pca = ();

#opendir DIR,$outdir or die $!;
chdir $outdir;
#close DIR;

$run_pca .= "$vcftools --vcf $vcf_file --plink --out pear\n";
$run_pca .= "$plink --noweb --file pear --make-bed --out pear_bfile\n";
$run_pca .= "$gcta --bfile pear_bfile --make-grm --out pear\n";
my $tmp_pca = "tep_pca"."$ARGV[2]";
$run_pca .= "$gcta --grm pear --pca $ARGV[2] --out $tmp_pca\n";

print $run_pca;
system("$run_pca");

my @col = ('red','green','blue','orange','yellow','black','pink','gold','brown','wheat','azure','gray','aquamarine','coral','navy','purple','powderblue');

my $order = -1;
my %individual = ();
my @color3D = ();
open IN,$ARGV[0] or die $!;
open OUT,'>',"$outdir/$name.eigenvec.legend" or die $!;
while(<IN>){
	chomp;
	my ($individual,$group) = (split "\t",$_)[2,5];
	if (not exists $hash{$group}){
		$order++;
		#warn "color is not enough!" unless (-e $col[$order]);
		$hash{$group} = join "\t",($col[$order],$order+1);
		print OUT "$group\t$col[$order]\t",$order+1,"\n";
	}
	push @color3D,$col[$order];
	$individual{$individual} = $hash{$group};
}

open IN,"$outdir/$tmp_pca.eigenvec" or die $!;
open OUT,'>',"$outdir/$name.eigenvec.final" or die $!;
while(<IN>){
	chomp;
	/^(\S+)/;
	my $k = $1;
	print OUT $_,"\t",$individual{$k},"\n";
}

my $first = 0;
my $color3D_rep = ();
my $count = ();
my $sum = ();
for(0..$#color3D){
	if ($_ == $#color3D){
		$count = $_ - $first + 1;
		$color3D_rep .= "rep\(\"$color3D[$_]\"\,$count)".",";
		$sum += $count;
		next;
	}

	if ($color3D[$_] ne $color3D[$_+1]){
		$count = $_ - $first + 1;
		$color3D_rep .= "rep\(\"$color3D[$_]\"\,$count)".",";
		$first = $_ + 1;
		$sum += $count;
	}

}
$color3D_rep =~ s/,$//;

open R,'>',"$outdir/$name.R" or die $!;

print R "setwd('$outdir')\n";
print R "data<-read.table('$name.eigenvec.final')\n";

print R "col<-data[",$ARGV[2]+3,"]\n";
print R "pch<-data[",$ARGV[2]+4,"]\n";
print R "legend<-read.table('$name.eigenvec.legend')\n";

for my $pca1(1..$ARGV[2]-1){
	for my $pca2(($pca1+1)..$ARGV[2]){
		my $name = "pc".$pca1.$pca2;
		print R "## $pca1 $pca2\n";
		print R "$name<-data[,c(",$pca1+2,",",$pca2+2,")]\n";
		my $two = "2D_"."$name".".svg";
		print R "svg(\"$two\")\n";
		print R "plot($name,ylab=\"PC$pca2\",xlab=\"PC$pca1\",pch=as.matrix(pch),col=as.matrix(col),main=\"PCA analysis 2D plot\")\n";
		print R "legend(\"topright\",as.matrix(legend[1]),pch=as.matrix(legend[3]),col=as.matrix(legend[2]),cex=0.5)\n";
		print R "dev.off()\n";
	}
}

for my $pca1(1..$ARGV[2]-2){
	for my $pca2(($pca1+1)..($ARGV[2]-1)){
		for my $pca3(($pca2+1)..$ARGV[2]){
			my $name = "pc".$pca1.$pca2.$pca3;
			print R "$name<-data[,c(",$pca1+2,",",$pca2+2,",",$pca3+2,")]\n";
			my $three = "3D_"."$name".".svg";
			print R "svg(\"$three\")\n";
			print R "library('scatterplot3d')\n";
			print R "scatterplot3d($name,type='h',highlight.3d=FALSE,color=c($color3D_rep),xlab=\"PC$pca1\",ylab=\"PC$pca2\",zlab=\"PC$pca3\",main=\"PCA analysis 3D plot\")\n";
			print R "legend(\"topright\",as.matrix(legend[1]),lty=c(1),col=as.matrix(legend[2]),cex=0.5)\n";
			print R "dev.off()\n";
		}
	}
}

=pub

$commond .= "svg('2d.svg')\n";
$commond .= "plot(pc12,ylab=\"PC2\",xlab=\"PC1\",pch=as.matrix(pch),col=as.matrix(col),main=\"PCA analysis 2D plot\")\n";
$commond .= "legend<-read.table('tmp_pca6.eigenvec.legend')\n";
$commond .= "legend(\"topright\",as.matrix(legend[1]),pch=as.matrix(legend[3]),col=as.matrix(legend[2]),cex=0.5)\n";
$commond .= "dev.off()\n";

$commond .= "svg('3d.svg')\n";
$commond .= "pc123<-data[3:5]\n";
$commond .= "library('scatterplot3d')\n";

$commond .= "scatterplot3d(pc123,type='h',highlight.3d=FALSE,color=c($color3D_rep),xlab=\"PC1\",ylab=\"PC2\",zlab=\"PC3\",main=\"PCA analysis 3D plot\")\n";
$commond .= "legend(\"topright\",as.matrix(legend[1]),lty=c(1),col=as.matrix(legend[2]),cex=0.5)\n";
$commond .= "dev.off()\n";
=cut

system("/usr/bin/Rscript $outdir/$name.R");