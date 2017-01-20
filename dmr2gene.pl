#!/usr/bin/perl -w 
#piemon^2016.6.22
use strict;

die "Usage:perl $0 <GFF> <DMR> <dmr2gene.out>\n\n" unless @ARGV == 3;

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

my %hash = ();
while(<IN>){
	next if (/^#/);
	chomp;
	my ($chr,$type,$start,$end,$string) = (split "\t",$_)[0,2..4,8];
	next if ($type ne 'gene');
	my $k = join "\t",($chr,$start,$end);
	$string =~ /;Name=(.*)/;
	my $gene = $1;
	$gene =~ s/;.*//g;
	$hash{$k} = $gene;
	#print $gene,"\n";
}

open IN,$ARGV[1] or die $!;
<IN>;
open OUT,'>',$ARGV[2] or die $!;
print OUT "#CHR\tSTART\tEND\tGENE\n";
while(<IN>){
	chomp;
	my ($chr,$start,$end) = (split "\t",$_)[0..2];
	$chr =~ s/\"//g;
	my $pos = int ($start + $end);
	my $anno = ();
	for my $k(keys %hash){
		my @line = split "\t",$k;
		next if ($line[0] ne $chr);
		if ($line[1] < $pos and $pos < $line[2]){
			$anno = $hash{$k};
			next;
		}
		if ($line[1] < $start and $start < $line[2]){
			$anno = $hash{$k};
			next;
		}
		if ($line[1] < $end and  $end < $line[2]){
			$anno = $hash{$k};
			next;
		}
	}

	$anno = 'NULL' if (! defined $anno);
	print OUT "$chr\t$start\t$end\t$anno\n";
}