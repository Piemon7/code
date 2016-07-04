#!/usr/bin/perl -w 
use strict;

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|";
}
else {
	open IN,$ARGV[0] or die $!;
}

my @unsorted = <IN>;

my @new = map { $_->[0] } 
sort { $a->[1]->[0] <=> $b->[1]->[0] or
	$a->[1]->[1] <=> $b->[1]->[1] }
	map {[$_,[ $_ =~ /^NC\_(\S+)\t(\d+)\t/ ]]} @unsorted;


#print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t15a15 16a11 16a12 16a13 16a3 16a4 16a6 16a9 17a1 17a3 17a5 17a7 17a8 18a3 18a7 19a14 19a1 19a2 19a3 19a9 1a10 1a11 1a13 1a15 1a4 1a5 1a8 20a11 20a13 20a14 20a3 20a8 20a9 21a10 21a12 21a13 21a14 21a8 22a1 22a2 22a3 22a6 2a11 2a12 2a2 2a5 2a7 3a11 3a15 3a6 3a8 3a9 4a13 4a14 4a2 4a7 5a10 5a4 5a5 5a7 5a9 6a11 6a14 6a3 6a4 6a7 7a11 7a12 7a7 7a8 XY12a2 XY7\n";
for (@new){
	print $_;
}