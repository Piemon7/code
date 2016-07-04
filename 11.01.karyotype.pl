#!/usr/bin/perl -w 
#piemon^2016.6.29.
use strict;
use File::Basename;

die "Usage:perl $0 <genome.fa> <cytosine.gz> <outdir> <bin_size>
bin_size is selectable.the default bin_size is 100000\n\n" unless @ARGV > 3;

my $outdir = $ARGV[2];

if ($ARGV[0] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else{
	open IN,$ARGV[0] or die $!;
}
my %hash = ();
$/ = ">";
<IN>;
while(<IN>){
	chomp;
	/^(\S+).*\n/;
	my $head = $1;
	$head =~ s/\_//g;
	s/^.*\n//;
	s/\s//g;
	my $len = length($_);
	$hash{$head} = $len;
	print $head,"\t",$len,"\n";
}
$/ = "\n";

my $name = basename $ARGV[0];
$name =~ s/\.gz$//g;
$name =~ s/\.fasta$//g;
my $karyotype = $name.".karyotype.txt";

# karyotype 
open OUT,'>',"$ARGV[2]/$karyotype"or die $!;
my $i = 1;
for (sort keys %hash){
	my $last = "chr".$i;
	$_ =~ s/\_//g;
	print OUT "chr - $_ ",$_," 0 $hash{$_} ",$last," \n";
	$i++;
}

my %histogram = ();
my $bin = ();
if((scalar @ARGV) > 3){
	$bin = $ARGV[3];
}
else{
	$bin = 100000;
}


if($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1] |" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
while(<IN>){
	chomp;
	my ($chr,$pos,$yes,$no,$type) = (split "\t",$_)[0,1,3,4,5];
	$chr =~ s/\_//g;
	next if (($yes + $no) == 0);
	my $percent = $yes/($yes+$no);
	$pos = int($pos/$bin);
	$histogram{$type}{$chr}->{$pos} += $percent;
}

# CG/CHG/CHH histogram.
for my $type("CG","CHH","CHG"){
	my $out = $type.".histogram.txt";
	open OUT,'>',"$ARGV[2]/$out" or die $!;
	for my $chr(sort keys %{$histogram{$type}}){
		for my $k(sort {$a<=>$b} keys %{$histogram{$type}{$chr}}){
			print OUT "$chr\t",$k*$bin,"\t",($k+1)*$bin,"\t",$histogram{$type}{$chr}->{$k},"\n";
		}
	}
}

# karyotype layout
my $layout = "karyotype.and.layout.".$name.".conf";
open OUT,'>',"$ARGV[2]/$layout" or die $!;
print OUT "karyotype = $ARGV[2]/$karyotype
chromosomes_order_by_karyotype = yes
chromosomes_units              = $bin
chromosomes_display_default    = yes";
`cp /lustre/project/og04/shichunwei/biosoft/circos-0.69-2/example/test/*.conf $outdir/`;
`cp /lustre/project/og04/shichunwei/biosoft/circos-0.69-2/example/test/circos.model.conf $outdir/circos.conf`;
`perl -p -i -e 's/model\.layout/$layout/g' $outdir/circos.conf`;
`perl -p -i -e 's/model\.CGhistogram/CG\.histogram\.txt/g' $outdir/circos.conf`;
`perl -p -i -e 's/model\.CHGhistogram/CHG\.histogram\.txt/g' $outdir/circos.conf`;
`perl -p -i -e 's/model\.CHHhistogram/CHH\.histogram\.txt/g' $outdir/circos.conf`;
`cp /lustre/project/og04/shichunwei/biosoft/circos-0.69-2/example/test/ticks.conf $outdir/ticks.conf`;

my $temp = `head -n 1 $outdir/genome.karyotype.txt |awk '{print \$4\"\t\"\$5\"\t\"\$6}'`;
chomp $temp;
`echo \"$temp	CG\" > $outdir/typeCG.txt`;
`echo \"$temp	CHG\" > $outdir/typeCHG.txt`;
`echo \"$temp	CHH\" > $outdir/typeCHH.txt`;

`perl -p -i -e 's/model\.typeCG/typeCG\.txt/g' $outdir/circos.conf`;
`perl -p -i -e 's/model\.typeCHG/typeCHG\.txt/g' $outdir/circos.conf`;
`perl -p -i -e 's/model\.typeCHH/typeCHH\.txt/g' $outdir/circos.conf`;
#model.typeCG

my $circos = "/lustre/project/og04/shichunwei/biosoft/circos-0.69-2/bin/circos";

system ("cd $outdir;$circos -conf $outdir/circos.conf -debug_group summary,timer");
