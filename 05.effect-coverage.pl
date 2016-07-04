#!/usr/bin/perl -w
#piemon^2016.5.31.
use strict;

die "Usage:perl <genome> <gff> <CX_report.txt> <05.out.1(bismark.bs)> << outdir.default is ./ >>

## the 05.out.1 add an annotation line at the last of CX_report file.
## the gff shoule contains 'gene','mRNA','CDS','exon','tRNA','repeat_region'.
## rewrite the script according to the gff format.\n\n" unless @ARGV >= 4;
## the 05.out.1 add an annotation line at the last of CX_report file.
my $outdir = ();
if (! defined $ARGV[4]){
	$outdir = "./";
}
else{
	$outdir = $ARGV[4];
}

my @chromosome = ();
my %chromosome = ();
$/ = ">";
if($ARGV[0] =~ /\.gz/){
	open IN,"gzip -dc $ARGV[0]|" or die $!;
}
else {
	open IN,$ARGV[0] or die $!;
}

<IN>;
while(<IN>){
	my $line = $_;
	chomp $line;
	$line =~ /^(\S+)/;
	my $chr = $1;
	push @chromosome,$chr;
	$line =~ s/^$chr//;
	$line =~ s/\s//g;
	$chromosome{$chr} = $line;
	print $chr,"\t",length($line),"\n";
}
$/ = "\n";
print "genome fasta file have been read!\n";


$/="\tmRNA\t";
if($ARGV[1] =~ /\.gz$/){
	open IN,"gzip -dc $ARGV[1]|" or die $!;
}
else{
	open IN,$ARGV[1] or die $!;
}
<IN>;

if (-e "$outdir/anno.txt"){
	`rm $outdir/anno.txt`;
}

while(<IN>){
	chomp;
	my @line = split "\n",$_;
	$line[0] =~ /^(\d+)\t(\d+)\t/;
	my $mRNA_start = $1;
	my $mRNA_end = $2;
	next if (!defined $mRNA_start);
	next if (!defined $mRNA_end);
	my (%exon,%intron,%cds,%five_UTR,%three_UTR);
	my $chromosome = ();
	my $strand = ();
	for (1..($#line-1)){
		next if ($line[$_] =~ /^#/);
		if ($line[$_] =~ /^(\S+)\t\S+\t\S+\t\d+\t\d+\t\S+\t(\S+)/){
			$chromosome = $1;
			$strand = $2;
			last;
		}
	}

	for my $line(1..($#line-1)){
		next if ($line[$line] =~ /^#/);
		$line[$line] =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t/;
		
		if ($2 eq 'exon'){
			for my $pos($3..$4){
				$exon{$pos} = 1;
			}
		}
		if ($2 eq 'CDS'){
			for my $pos($3..$4){
				$cds{$pos} = 1;
			}
		}
	}

	if (%cds){
		my @cds = keys %cds;
		@cds = sort {$a<=>$b} @cds;
		my $left = shift @cds;
		my $right = pop @cds;
		
		for ($mRNA_start..$left){
			if ($strand eq '+'){
				$five_UTR{$_} = 1;
			}
			else{
				$three_UTR{$_} = 1;
			}
		}
		for($right..$mRNA_end){
			if ($strand eq '+'){
				$three_UTR{$_} = 1;
			}
			else{
				$five_UTR{$_} = 1;
			}
		}
	}
	

	open OUT,'>>',"$outdir/anno.txt" or die $!;
	for($mRNA_start..$mRNA_end){
		my $base = substr ($chromosome{$chromosome},($_-1),1);
		if (! defined $base) {print $chromosome,"\t",$mRNA_start,"\t",$mRNA_end,"\n";}
		$base =~ tr/cg/CG/;
		next if ($base ne 'C' and $base ne 'G');
		print OUT "$chromosome\t$_\t";
		if (exists $five_UTR{$_}){
			print OUT " 5_utr";
		}
		if(exists $three_UTR{$_}){
			print OUT " 3_utr";
		}
		if(exists $exon{$_}){
			print OUT " exon";
		}
		if(not exists $exon{$_}){
			print OUT " intron";
		}
		if(exists $cds{$_}){
			print OUT " cds";
		}
			print OUT " mRNA\n";
	}

}

$/="\n";
print "UTR region have been found !\n";

if (-e $ARGV[3]){
	`rm $ARGV[3]`;
}

for my $chromosome(@chromosome){
	if($ARGV[1] =~ /\.gz$/){
		open IN,"gzip -dc $ARGV[1]|" or die $!;
	}
	else{
		open IN,$ARGV[1] or die $!;
	}
	my (%gene,%tRNA,%repeat);
	while(<IN>){
		next if (/^#/);
		/^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t/;
		my ($chr,$type,$start,$end)=($1,$2,$3,$4);		
		#next if ($type eq 'region');
		next if ($chr ne $chromosome);
		if ($type eq 'gene'){
			for ($start..$end){
				$gene{$_} = 1;
			}
		}
		if ($type eq 'tRNA'){
			for ($start..$end){
				$tRNA{$_} = 1;
			}
		}
		if ($type eq 'repeat_region'){
			for ($start..$end){
				$repeat{$_} = 1;
			}
		}
	}

	print "gff $chromosome part have been all read !\n";

	die "anno.txt is not exists !\n" unless (-e "$outdir/anno.txt");
	
	open IN,"$outdir/anno.txt" or die $!;
	my %anno = ();
	while(<IN>){
		chomp;
		my ($chr,$pos,$info)=(split "\t",$_)[0,1,2];
		next if ($chr ne $chromosome);
		$anno{$pos} = $info;
	}


	if($ARGV[2] =~ /\.gz$/){
		open IN,"gzip -dc $ARGV[2]|" or die $!;
	}
	else{
		open IN,$ARGV[2] or die $!;
	}

	open OUT,'>>',$ARGV[3] or die $!;
	while(<IN>){
		chomp;
		my $row = $_;
		my @line = split "\t",$_;
		my $pos = $line[1];
		next if ($line[0] ne $chromosome);
		print OUT "$row\t";
		if (exists $gene{$pos}){
			print OUT "gene ";
		}
		else{
			print OUT "intergenic ";
		}

		if(exists $anno{$pos}){
			print OUT "$anno{$pos} ";
		}
		if(exists $tRNA{$pos}){
			print OUT "tRNA ";
		}
		if(exists $repeat{$pos}){
			print OUT "repeat ";
		}
		print OUT "\n";
	}
}
print "ALL done!!!\n";
