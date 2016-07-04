#!/usr/bin/perl -w 
#piemon^2016.6.6

use strict;

die "Usage:perl <cytosine>
coverage > 4 and ml > 80%
\n\n" unless @ARGV == 1;

if($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}
my %hash = ();
while(<IN>){
	chomp;
	my ($chr,$pos,$yes,$no,$type) = (split "\t",$_)[0,1,3,4,5];
	next if (($yes+$no) == 0);
	next if (($yes+$no) < 2);
	my $percent = $yes/($yes+$no);
	$hash{$type}++ if ($percent > 0.8 and ($yes+$no)> 4);
}

my $mC = ();
for (keys %hash){
	print $_,"\t",$hash{$_},"\n";
	$mC += $hash{$_};
}
print "mC\t$mC\n";
for (keys %hash){
	print "$_\t",$hash{$_}/$mC,"\n";
}
