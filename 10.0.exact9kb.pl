#!/usr/bin/perl -w 

#piemon^2016.6.3.
use strict;

die "Usage:perl $0 <9kb.summary.report><CG|CHH|CHG><weblogo-prefile.out>\n\n" unless @ARGV == 3;

my $C_type = $ARGV[1] or die $!;

open IN,$ARGV[0] or die $!;

my $base_num = 10000;
my %hash = ();
while(<IN>){
	chomp;
	my ($type,$count,$A,$C,$G,$T) = split "\t",$_;
	next if ($type ne $C_type);
	$A = int($A*$base_num);
	$C = int($C*$base_num);
	$G = int($G*$base_num);
	$T = int($T*$base_num);
	
	for(0..$A-1){
		push @{$hash{$count}},"A";
	}
	for(0..$C-1){
		push @{$hash{$count}},"C";
	}
	for(0..$G-1){
		push @{$hash{$count}},"G";
	}
	for(0..$T-1){
		push @{$hash{$count}},"T";
	}

}

my @out = ();
for my $line(0..$base_num-1){
	my $out = ();
	for my $base(1..12){
		$out .= ${$hash{$base}}[$line];
	}
	push @out,$out;
}

open OUT,'>',$ARGV[2] or die $!;
for (@out){
	my $len = length $_;
	next if ($len < 12);
	print OUT "$_\n";
}