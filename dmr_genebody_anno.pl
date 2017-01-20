#!/usr/bin/perl -w 
#piemon^2016.6.22
use strict;

die "Usage:perl $0 <GFF> <DMR.txt> <dmr_genebody.anno> \n\n" unless @ARGV == 3;

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

my %hash = ();

