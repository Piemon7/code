#/usr/bin/perl -w 
#piemon^2016.6.22.
use strict;
use File::Basename;
use Getopt::Long;

my ($outdir,$ref,$list,$gff);
GetOptions(
	'list=s' => \$list,	
	'outdir|o=s' => \$outdir,
	'ref=s' => \$ref,
	'gff=s' => \$gff,
);

my $usage = <<USAGE;
Usage :  perl $0 -l -outdir -ref.\n
options:
	-list     a list(phred33-quals). pair-end format data.
	-outdir   outdir. default is ./
	-ref	  reference (fasta or fa).
	-gff      gff3
	all need absolute path. 

	list format 
	sample1 /path/sample1_1.fastq /path/sample1_2.fastq
	sample2 /path/sample2_1.fastq /path/sample2_2.fastq
	the name shoule be same.

	the insert size range: 50 ~ 500.

	-help     show this help.\n
USAGE

print $usage and exit unless ( $list and $ref and $gff);

my $bismark_dir = "/lustre/project/og04/shichunwei/biosoft/bismark_v0.16.1/";
my $bowtie2_dir = "/lustre/project/og04/shichunwei/biosoft/bowtie2-2.2.9/";
my $perl        = "/lustre/project/og04/shichunwei/bin/perl";
my $samtools    = "/lustre/project/og04/shichunwei/biosoft/samtools-1.3/samtools";
my $picard_dir  = "/lustre/project/og04/shichunwei/biosoft/picard-tools-1.119/";
my $java        = "/lustre/project/og04/shichunwei/biosoft/jre1.8.0_91/bin/java";
my $script      = "/lustre/project/og04/shichunwei/project/bisulfite_sequencing/script/";
my $qsub        = "/lustre/project/og04/shichunwei/project/temp/GATK/qsub-sge.pl";

$outdir = $outdir || `pwd`;

my %list = ();
open IN,$list or die $!;
while(<IN>){
	chomp;
	my @line = split;
	$list{$line[0]}{'1'} = $line[1];
	$list{$line[0]}{'2'} = $line[2];
}

## bismark call.
`mkdir -p $outdir/ref`;
`ln -s $ref $outdir/ref/genome.fasta`;
my $ref_dir = "$outdir/ref";
my $commond_bismark_preparation = ();
$commond_bismark_preparation .= "$perl $bismark_dir/bismark_genome_preparation --path_to_bowtie $bowtie2_dir --verbose --bowtie2  $ref_dir\n";
my $commond_bismark_call = ();
my $commond_bismark_extract = ();
for (keys %list){
	$commond_bismark_call .= "mkdir -p $outdir/$_\n";
	$commond_bismark_call .= "$perl $bismark_dir/bismark -q --phred33-quals -p 8 --multicore 8 -N 1 -I 50 -X 500 -1 $list{$_}{'1'} -2 $list{$_}{'2'} $ref_dir --samtools_path $samtools -o $outdir/$_\n";
	my $bam_file = "$_"."_1"."_bismark_bt2_pe.bam";
	$commond_bismark_extract .=	"$perl $bismark_dir/bismark_methylation_extractor -p --no_overlap --comprehensive --no_header --report --gzip --bedGraph --counts --CX --buffer_size 10G --cytosine_report $outdir/$_/$bam_file --genome_folder $ref_dir -o $outdir/$_ --multicore 8\n";
}

print $commond_bismark_preparation;
print $commond_bismark_call;
print $commond_bismark_extract;

## bam 2 sort sam 
my $commond_sort = ();
## 00 插入片段长度检测
my $commond_insert_size = ();
## 01 不同reads测序深度下的覆盖度分布
my $commond_allbases_depth = ();
## 02 C碱基有效测序深度的累积分布
my $commond_Csite_distribute = ();
## 03 C位点的有效覆盖度分布图
## 04 样品的甲基化水平分布图
my $commond_Csite_effect = ();

## 06 mC count & percent
my $commond_mC_count = ();

## 05 不同基因组区域的有效覆盖度
## 09 不同基因组区域内 CG、CHG 和 CHH 中 C 的甲基化水平
my $commond_region_coverage = ();

## 08 CG、CHG 和 CHH 中的所有 C 的甲基化水平
my $commond_methy_level = ();

## 10 9kb sequence 
my $commond_9kb = ();
## 11 circle 图
my $commond_circle = ();

## 12 基因组不同转录元件中的 DNA 甲基化水平
my $commond_element_ml = ();

## R methylKit 
my $commond_methylKit = ();

for (keys %list){
	my $bam_file = "$_"."_1"."_bismark_bt2_pe.bam";
	my $bam_sort = $bam_file.".sort";
	my $sam_sort = $bam_sort.".sam";
	$commond_sort .= "$samtools sort -o $outdir/$_/$bam_sort $outdir/$_/$bam_file -@ 8\n";
	$commond_sort .= "$samtools view -O sam -o $outdir/$_/$sam_sort $outdir/$_/$bam_sort -@ 8\n";
	$commond_insert_size .= "$java -jar $picard_dir/CollectInsertSizeMetrics.jar I=$outdir/$_/$bam_sort O=$outdir/$_/pe.insertSize.txt H=$outdir/$_/pe.insert.pdf\n";
	$commond_insert_size .= "$script/00.insert_size.pl $outdir/$_/ $outdir/$_/pe.insertSize.txt\n";
	$commond_allbases_depth .= "$samtools depth -a $outdir/$_/$bam_sort > $outdir/$_/depth.txt\n";
	$commond_allbases_depth .= "$perl $script/01.depth-coverage.pl $outdir/$_/depth.txt $ref $outdir/$_/01.out $outdir/$_\n";
	my $CX_report_file = $_."_1"."_bismark_bt2_pe.CX_report.txt.gz";
	$commond_Csite_distribute .= "$perl $script/02.Csite_distribute.pl $outdir/$_/$CX_report_file $outdir/$_/02.out $outdir/$_\n";
	$commond_Csite_effect .= "$perl $script/03.C_effect_coverage.pl  $outdir/$_/$CX_report_file 03.out 04.out $outdir/$_\n";
	
	$commond_mC_count .= "$perl $script/06.cytosine_count.pl $CX_report_file > $outdir/$_/06.out\n";

	open METH,'>',"$outdir/$_\.R" or die $!;
	my $commond_methylKit = (); 
	$commond_methylKit .= "setwd(\'$outdir/$_\')\nlibrary(methylKit)\n";
	$commond_methylKit .= "myobj=read.bismark(\"$sam_sort\",\"$_\",assembly=\"hg19\",save.context=c(\"CpG\",\"CHG\",\"CHH\"),read.context=\"none\",mincov=5,minqual=20,save.folder = getwd())\n";
	print METH $commond_methylKit;

	$commond_region_coverage .= "$perl $script/05.effect-coverage.pl $ref $gff $outdir/$_/$CX_report_file $outdir/$_/05.out.1 $outdir/$_\n";
	$commond_region_coverage .= "$perl $script/05.effect-coverage-2.pl $outdir/$_/05.out.1 05.out.2 09.out $outdir\n";
	$commond_methy_level .= "$perl $script/08.cytosine-density.pl $outdir/$_/$CX_report_file 4 $outdir/$_\n";

	$commond_9kb .= "$perl $script/10.exact9kb.pl $ref $CX_report_file $outdir/$_\n";
	$commond_circle .= "$perl $script/11.01.karyotype.pl $ref $CX_report_file $outdir/$_ 100000\n";	
	$commond_element_ml .= "$perl $script/12.TSS_C_level.pl $gff $ref $outdir/$_/05.out.1 12.out $outdir/$_\n";

}
print $commond_sort;
print $commond_insert_size;
print $commond_allbases_depth;
print $commond_Csite_distribute;
print $commond_Csite_effect;
print $commond_mC_count;
print $commond_region_coverage;
print $commond_methy_level;
print $commond_9kb;
print $commond_circle;
print $commond_element_ml;

## R methylKit
open METH,'>',"$outdir/save.contest.sh" or die $!;
for (keys %list){
	print METH "Rscript $outdir/$_\.R\n";
}
close METH;
