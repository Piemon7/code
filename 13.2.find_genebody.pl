#!/usr/bin/perl -w 
#piemon^2016.6.27
use strict;
use File::Basename;
die "Usage:perl $0 <gff> <DMR.txt> <outdir>\n\n" unless @ARGV == 3;

my %hash = ();
open IN,$ARGV[1] or die $!;
while(<IN>){
	chomp;
	my ($chr,$start,$end) = (split "\t",$_)[0..2];
	$chr =~ s/\"//g;
	next if ($chr eq 'seqnames');
	push @{$hash{$chr}},(join "\t",($start,$end));
}
close IN;

$/ = "\tmRNA\t";
if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

<IN>;


my %result = ();

while(<IN>){
	chomp;
	$_ =~ s/(\n#.*)+//g;
	my @line = split "\n",$_;

	my @filter = ();
	for (0..$#line-1){
		push @filter,$line[$_];
	}

	$filter[0] =~ /^(\d+)\t(\d+)\t\S+\t(\S+)/;
	my ($mrna_1,$mrna_2,$strand) = ($1,$2,$3);

	$filter[1] =~ /^(\S+)/;
	my $chr = $1;
	next if (not exists $hash{$chr});
	
	my @array = @{$hash{$chr}};
	my %type = ();
	
	next if ((scalar @filter) < 2);
	
	for (1..$#filter){
		my ($type,$start,$end) = (split "\t",$filter[$_])[2..4];
		push @{$type{$type}},(join "\t",($start,$end));
	}

	my @anno = ();
	for my $dmr(@array){
		my ($pos1,$pos2) = split "\t",$dmr;
		next if ($pos2<$mrna_1 or $mrna_2<$pos1);
		my @cds = ();
		## cds 
		if (exists $type{'CDS'}){
			for my $cds(@{$type{'CDS'}}){
				my ($cds1,$cds2) = split "\t",$cds;
				push @cds,($cds1,$cds2);
				if (($cds1 <$pos1 and $pos1 < $cds2) or ($cds1<$pos2 and $pos2<$cds2)){
					push @anno,'CDS';
				}
			}
		}
		## exon
		my @intron = ();
		if (exists $type{'exon'}){
			for my $exon(@{$type{'exon'}}){
				my ($exon1,$exon2) = split "\t",$exon;
				push @intron,($exon1,$exon2);
				if (($exon1 <$pos1 and $pos1 < $exon2) or ($exon1<$pos2 and $pos2<$exon2)){
				push @anno,'exon';
				}
			}
		}

		## intron
		if (exists $type{'exon'}){
			shift @intron;
			pop @intron;
			for(my $i=0;$i<$#intron;$i+=2){
				if (($intron[$i]<$pos1 and $pos1<$intron[$i+1]) or ($intron[$i]<$pos1 and $pos2<$intron[$i+1])){
					push @anno,'intron';
				}
			}
		}
		## utr
		if (exists $type{'CDS'}){
			my @sort = sort {$a<=>$b} @cds;
			#$strand
			my $utr_left  = shift @sort;
			my $utr_right = pop @sort;
			if (($mrna_1<$pos1 and $pos1<$utr_left) or ($mrna_1 < $pos2 and $pos2<$utr_left)){
				if ($strand eq '+'){
					push @anno,'utr3';
				}
				else{
					push @anno,'utr5';
				}
			}
			if (($utr_right<$pos1 and $pos1<$mrna_2) or ($utr_right<$pos2 and $pos2<$mrna_2)){
				if ($strand eq '+'){
					push @anno,'utr5';
				}
				else{
					push @anno,'utr3';
				}
			}
		}
		my $out = ();
		if ((scalar @anno) >0){
			my %uniq = ();
			for (@anno){
				$uniq{$_} = 1;
			}
			$out = join " ",keys %uniq;
			$result{join "\t",($chr,$pos1)} = $out;
		}
	}
}

$/ = "\n";

my $outfile_1 = basename $ARGV[1];
$outfile_1 =~ s/\.txt$//g;
$outfile_1 .= ".anno";
my $outfile_2 = $outfile_1.'.txt';
my $outdir = $ARGV[2];

my %genebody = (
	'CDS' =>'0',
	'utr5' => '0',
	'utr3' => '0',
	'intron' =>'0',
	'exon' => '0',
	#'promoter' =>'0',
	);

open OUT,'>',"$outdir/$outfile_1" or die $!;
open IN,$ARGV[1] or die $!;
<IN>;
while(<IN>){
	chomp;
	my ($chr,$start) = (split "\t",$_)[0,1];
	$chr =~ s/\"//g;
	my $line = $_;
	if (exists $result{join "\t",($chr,$start)}){
		my $anno = $result{join "\t",($chr,$start)};
		print OUT $line,"\t",$anno,"\n";
		my @anno = split " ",$anno;
		for(@anno){
			$genebody{$_} += 1;
		}
		
	}
	else{
		print OUT $line,"\tNA\n";
	}
}


open OUT,'>',"$outdir/$outfile_2" or die $!;
for(sort keys %genebody){
	print OUT "$_\t$genebody{$_}\n";
}
