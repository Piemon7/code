#!/usr/bin/perl -w 
#piemon^2016.6.29.
use strict;

die "Usage:perl $0 <genome.fasta> <cytosine> <outdir>\n
weblogo site: http://weblogo.berkeley.edu/logo.cgi.
upload the file CHG.weblogo. 
select 'Frequency Plot' if necessary.
unselect 'Show fine print'.
unselect 'Label Sequence Ends'.
\n" unless @ARGV == 3;

my $outdir = "$ARGV[2]/9kb.weblogo";
`mkdir -p $outdir`;

my %chromosome = ();

$/ = ">";

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}

#my ($C,$G,$CG);

<IN>;
while(<IN>){
	chomp;
	my $line = $_;
	my $head = (split "\n",$_)[0];
	$line =~ s/$head//;
	$line =~ s/\s//g;
	$head =~ /^(\S+)/;
	$head = $1;
	print "$head\t",length($line),"\n";
	$chromosome{$head} = $line;
}
$/ = "\n";
print "genome info have been read!!!\n";

open IN,$ARGV[1] or die $!;
if ($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}

my %hash = ();
my %C_count = ();

while(<IN>){
	my ($chr,$pos,$string,$yes,$no,$type) = (split "\t",$_)[0..5];
	next if ($no == 0);
	my $percent = $yes/($yes + $no);
	next if ($percent < 0.8);# M level > 80%.
	
	my $ninebases = ();
	if ($string eq "+"){
		$ninebases = substr($chromosome{$chr},$pos-6,12);
		push @{$hash{$type}},$ninebases;
	}
	elsif($string eq "-"){
		$ninebases = substr($chromosome{$chr},$pos-7,12);
		$ninebases = reverse $ninebases;
		$ninebases =~ tr/AGCTagct/TCGAtcga/;
		push @{$hash{$type}},$ninebases;
	}
	else {
		warn "$string is wrong !!!\n";
	}

	if ($type eq 'CHH'){
		$C_count{$type} += 1;
	}
	elsif($type eq 'CHG'){
		$C_count{$type} += 1;
	}
	elsif($type eq 'CG'){
		$C_count{$type} += 1;
	}
	else{
		warn "wrong!!!\n";
	}
}

my $temp_file = "$outdir/9kb.temp";
open OUT,'>',$temp_file or die $!;

for my $k(sort keys %hash){
	#print OUT "$k\n";
	#print OUT "\tA\tC\tG\tT\n";
	for my $count(0..11){
		print OUT "$k\t";
		print OUT $count+1,"\t";
		my (%sum,$num);
		for my $line(@{$hash{$k}}){
			my $base = substr ($line,$count,1);
			$base =~ tr/agct/AGCT/;
			$sum{$base} += 1;
			$num++;
		}
		my @array = ('A','C','G','T');
		for (@array){
			
			if(exists $sum{$_}){
				print OUT $sum{$_}/$num,"\t";
			}
			else{
				print OUT "0\t";
			}
		}
		print OUT "\n";
	}
}

my $perl = "/lustre/project/og04/shichunwei/project/bisulfite_sequencing/script/10.0.exact9kb.pl";
for(keys %C_count){
	if ($C_count{$_} > 10000){
		system("perl $perl $temp_file $_ $outdir/$_\.weblog");
	}
	else{
		my @array = @{$hash{$_}};
		open OUT,'>',"$outdir/$_\.weblog" or die $!;
		for(@array){
			print OUT $_,"\n";
		}
	}
}

system("rm $temp_file");
