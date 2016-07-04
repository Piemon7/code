#!/usr/bin/perl -w 
#piemon^2016.6.28
use strict;
use File::Basename;
die "usage:perl $0 <bam> <fasta> <05.out.1/bs> <outdir>\n\n" unless @ARGV == 4;

my $outdir = $ARGV[3];

my @chr = ();
my @array = `grep ">" $ARGV[1]`;
for(@array){
	chomp;
	/^>(\S+)/;
	push @chr,$1;
}

#my %chr = ();
open IN,"samtools view $ARGV[0]|" or die $!;
while(<IN>){
	chomp;
	my ($chr,$pos,$string) = (split "\t",$_)[2,3,13];
	#$chr{$chr} = 1;
	open OUT,'>>',"$outdir/$chr.chr.temp" or die $!;
	$string =~ s/^XM:Z://;
	print OUT "$pos\t$string\n";
}
#my @chr = keys %chr;
print "chr num:",scalar @chr,"\n";
print `date`,"\n";

for my $chr(@chr){
	my $judge = 0;
	if (-e "$outdir/$chr.chr.temp"){
		$judge = 1;
	}
	next if ($judge == 1);

	my %position = ();
	open IN,"$outdir/$chr.chr.temp" or die $!;
	while(<IN>){
		chomp;
		my @line = split;## pos,string
		my $len = length($line[1]);
		for my $count(0..$len){
			my $dot = substr($line[1],$count,1);
			next if ($dot eq '.');
			if ($dot eq 'h' or $dot eq 'x' or $dot eq 'z'){
				$position{$line[0]+$count+1}{'1'} += 1;
			}
			elsif($dot eq 'H' or $dot eq 'X' or $dot eq 'Z'){
				$position{$line[0]+$count+1}{'2'} += 1;
			}
			else{
				warn;
			}
		}
	}

	open OUT,'>',"$outdir/$chr.out.temp" or die $!;
	for(sort {$a<=>$b} keys %position){
		print OUT "$_\t",$position{$_}{'2'},"\t",$position{$_}{'1'},"\n";
	}
	`rm $outdir/$chr.chr.temp`;
	print "$chr out.temp done\n";
}

=pod 
if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
$/ = ">";
<IN>;
while(<IN>){
	my %hash = ();
	chomp;
	my $line = $_;
	my $head = (split "\n",$_)[0];
	$line =~ s/^$head//;
	$head = (split " ",$head)[0];
	open OUT,'>',"$outdir/$head.bs.temp" or die $!;
	$line =~ s/\s//g;
	$line =~ tr/agctn/AGCTN/;
	my $len = length($line);
	for my $pos(0..$len){
		my $string = substr($line,$pos,1);
		if ($string eq 'C' ){
			my $three = substr($line,$pos,3);
			print OUT "$head\t",$pos+1,"\t+\t$three\n";
		}
		elsif($string eq 'G'){
			my $three = substr($line,$pos,3);
			print OUT "$head\t",$pos+1,"\t-\t$three\n";
		}
	}
}
$/ = "\n";
=cut

open OUT,'>',$ARGV[2] or die $!;
for my $chr(@chr){
	my %uniq = ();
	if (-e "$outdir/$chr.out.temp"){
		open IN,"$outdir/$chr.out.temp" or die $!;
		while(<IN>){
			my @line = split;# pos,yes,no
			$uniq{$line[0]} = join "\t",@line[1,2];
		}		
	`rm $outdir/$chr.out.temp`;
	}

	open IN,"$outdir/$chr.bs.temp" or die $!;
	while(<IN>){
		my @line = split;# chr,pos,strand,base
		my $out1 = join "\t",@line[0..2];
		my $out2 = ();
		if (exists $uniq{$line[1]}){
			$out2 = $uniq{$line[1]};
		}
		else{
			$out2 = join "\t",('0','0');
		}
		print OUT $out1,"\t",$out2,"\t$line[3]\n";
	}
	`rm $outdir/$chr.bs.temp`;
}
