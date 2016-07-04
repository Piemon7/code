#!/usr/bin/perl -w 
#piemon^2016.6.15.
use strict;
die "Usage:perl $0 <gff><bed>\n\n" unless @ARGV == 2;

if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

open OUT,'>',$ARGV[1] or die $!;

## find utr 
$/ = "\tmRNA\t";
<IN>;
while(<IN>){
	my (%up,%hash,$direction,@utr);
	chomp;
	my @line = split "\n",$_;
	my @filter = ();
	for (@line){
		next if ($_ =~ /^#/);
		next if ($_ =~ /^\S+\t\S+\tregion/);
		push @filter,$_;
	}
	$filter[0] =~ /^(\d+)\t(\d+)\t\S+\t(\S+)/;
	push @utr,($1,$2);
	if($3 eq '+'){
		$direction = 'f';
	}
	elsif($3 eq '-'){
		$direction = 'r';
	}
	else{
		warn;
		print $filter[0],"\n\n";
	}

	$filter[0] =~ /gene=(\S+);/;
	my $gene = $1;
	$gene =~ s/;.*//g;
	$filter[1] =~ /^(\S+)\t/;
	my $chr = $1;
	#$chr =~ s/\_//g;
	#$chr =~ s/\.//g;

	my $cds_count = 0;
	for(1..$#filter-1){
		$filter[$_] =~ /^\S+\t\S+\t(\S+)\t/; 
		my $type = $1;
		$cds_count ++ if ($type eq 'CDS');
		push @{$hash{$type}},$filter[$_];
	}
	next if $cds_count == 0;

	for(@{$hash{'CDS'}}){
		my @line = split "\t",$_;
		push @utr,@line[3,4];
	}

	my @sorted_utr = sort {$a<=>$b} @utr;
	my $left = shift @sorted_utr;
	$left = $left."\t".(shift @sorted_utr);

	my $right = pop @sorted_utr;
	$right =(pop @sorted_utr)."\t".$right;
	
	if($direction eq 'f'){
		# NM_017829_utr3_*_*_chr22_18216906_r
		my @left = split "\t",$left;
		print OUT "$chr\t$left\t",('NM_123456_utr5_*_*_'.$chr.'_'.($left[0]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
		my @right = split "\t",$right;		
		print OUT "$chr\t$right\t",('NM_123456_utr3_*_*_'.$chr.'_'.($right[0]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
	}
	elsif($direction eq 'r'){
		my @right = split "\t",$right;
		print OUT "$chr\t$right\t",('NM_123456_utr5_*_*_'.$chr.'_'.($right[0]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
		my @left = split "\t",$left;
		print OUT "$chr\t$left\t",('NM_123456_utr3_*_*_'.$chr.'_'.($left[0]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
	}
	else{
		warn;
	}
	
}

$/ = "\n";

## up 5000

if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

while(<IN>){
	chomp;
	next if (/^#/);
	my @line = split "\t",$_;

	my ($region,$start,$end,$type) = @line[2..4,6];
	next if ($region ne 'gene');

	$line[-1] =~ /;Name=(.*);/;
	my $gene = $1;
	$gene =~ s/;.*//g;
	#$gene = shift(split ";",$gene);

	if (! defined $gene){
		print "$_\n\n";
	}

	my $chr = $line[0];
	#$chr =~ s/\_//g;
	#$chr =~ s/\.//g;

	my ($out1,$out2);
	if ($type eq '+'){
		my $direction = 'f';
		$out1 = $start-5000;
		$out1 = 0 if ($out1 < 0);
		$out2 = $start + 1;
		print OUT "$chr\t$out1\t$out2\t",('NM_123456_up_5000_'.$chr.'_'.$out1.'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
	}
	elsif($type eq '-'){
		my $direction = 'r';
		$out2 = $end + 5000;
		$out1 = $end + 1;
		print OUT "$chr\t$out1\t$out2\t",('NM_123456_up_5000_'.$chr.'_'.$out1.'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";		 
	}
	else{
		warn;
	}
}


## intron  and   CDS

if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

$/ = "\tmRNA\t";
<IN>;
while(<IN>){
	chomp;
	my @line = split "\n",$_;
	my @filter_exon = ();
	my @filter_cds  = ();
	my ($chr,$type,$start,$end,$direction);
	for (@line){
		next if ($_ =~ /^#/);
		#next if ($_ =~ /^\S+\t\S+\tregion/);
		if (/^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)/){
			($chr,$type,$start,$end,$direction) = ($1,$2,$3,$4,$5);
			if ($type eq 'exon'){
				push @filter_exon,($start,$end);
			}
			elsif($type eq 'CDS'){
				push @filter_cds,($start,$end);
			}
		}
	}
	next if (!defined $direction);
	$direction = 'f' if ($direction eq '+');
	$direction = 'r' if ($direction eq '-');

	$line[0] =~ /gene=(\S+)/;
	my $gene = $1;
	$gene =~ s/;.*//g;
	next if (!defined $gene);

	my @sorted_exon = sort {$a<=>$b} @filter_exon;
	my @sorted_cds  = sort {$a<=>$b} @filter_cds;
	for(my $i=1;$i<$#sorted_exon-1;$i+=2){
		my $out = join "\t",($sorted_exon[$i],$sorted_exon[$i+1]);
		print OUT "$chr\t$out\t",('NM_123456_intron_'.int($i/2).'_*_'.$chr.'_'.($sorted_exon[$i]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
	}

	for(my $i=0;$i<$#sorted_cds-1;$i+=2){
		my $out = join "\t",($sorted_cds[$i],$sorted_cds[$i+2]);
		print OUT "$chr\t$out\t",('NM_123456_cds_'.int($i/2).'_*_'.$chr.'_'.($sorted_cds[$i]+1).'_'.$direction),"\t*\t*\tNM_123456\t$gene\t*\t*\n";
	}
}

$/ = "\n";
