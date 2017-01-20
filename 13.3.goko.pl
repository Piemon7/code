#!/usr/bin/perl -w 
#piemon^2016.6.28
use strict;
use Getopt::Long;
my ($protein,$species,$outdir,$nr,$kegg,$gff,$dmrfile);
GetOptions(
	'p=s' => \$protein,
	's=s' => \$species,
	'o=s' => \$outdir,
	'nr=s' => \$nr,
	'kegg=s' => \$kegg,
	'gff=s' => \$gff,
	'i=s' => \$dmrfile,
);

my $usage = <<USAGE;
usage:
	-p 	protein,fa/fasta
	-s 	species.default is 'species'.
	-o 	outdir.default is ./
	-nr 	plant,animal,bacteria or mammals
	-kegg 	eukaryote,prokaryote,all
	-gff 	gff3
	-i 	dmrfile
eg.
	perl $0 -p protein.fa -s human -o outdir -nr mammals -kegg eukaryote -gff human.gff -i DMR.txt
	perl $0 -p protein.fa -nr mammals -kegg eukaryote -gff human.gff -i DMR.txt
USAGE

die $usage unless ($protein and $nr and $kegg and $gff and $dmrfile);

$outdir  ||= `pwd`;
$species ||= "species";

my $perl = "/lustre/project/og04/shichunwei/bin/perl";
my $ko_go = "/nfs/pipe/RNA/RNA-seq/version1/functional/kogo_annotation/ko_go.pl";
my $script = "/lustre/project/og04/shichunwei/project/bisulfite_sequencing/script/kogo";

my %nr = (
	'plant'    => '/nfs/database/db/Pub/nr/DNA/20150603/Plants.fa',
	'animal'   => '/nfs/database/db/Pub/nr/DNA/20150603/animal.fa',
	'bacteria' => '/nfs/database/db/Pub/nr/DNA/20150603/Bacteria.fa',
	'mammals'  => '/nfs/database/db/Pub/nr/DNA/20150603/Mammals.fa',
	);

my %kegg = (
	'eukaryote'  => "/nfs/database/db/Pub/kegg/DNA/current/kegg_eukaryote_clean.fa",
	'prokaryote' => "/nfs/database/db/Pub/kegg/DNA/current/kegg_prokaryote_clean.fa",
	'all'        => "/nfs/database/db/Pub/kegg/DNA/current/kegg_all_clean.fa",
	);

my $nr_db   = $nr{$nr};
my $kegg_db = $kegg{$kegg};

my $commond_kogo = ();
$commond_kogo .="$perl $ko_go -ko -go -input $protein -species $species -blast blastp -evalue 1e-5 -rank 5 -kegg $kegg_db -nr $nr_db -obo /nfs/onegene/user/1gene/liangq/database/go-basic.obo -outdir $outdir\n";
#system("$commond_kogo");
print $commond_kogo;

my $commond = ();
$commond .= "$perl $script/03.get_gene_protein_info.pl $gff $outdir/$species.gene2pro\n";
$commond .= "sh /nfs/pipe/RNA/RNA-ref/version1/functional/kogo_annotation/DealannotTr2Gene.sh  $outdir/$species.gene2pro $outdir/$species.annot $outdir/$species.ko $outdir $species $outdir/$species.nr.desc\n";
$commond .= "path=/nfs/pipe/RNA/RNA-ref/version1/functional\n";
$commond .= "export LD_LIBRARY_PATH=/nfs/config/boost/lib:\$LD_LIBRARY_PATH\n";
$commond .= "$perl \$path/dealGOObo.pl -go /nfs/onegene/user/1gene/liangq/database/go-basic.obo  -prefix $outdir/go\n";
$commond .= "cat $outdir/$species.C $outdir/$species.F $outdir/$species.P |awk -F '\t' '{a[\$3]=a[\$3]\",\"\$5} END {for (i in a) {print i\"\t\"a[i]}}' |sed -e 's/,//1' > $outdir/$species.gene2goMapping\n";
#system("$commond");
print $commond;

my $anno = ();
$anno .= "$perl $script/00.extract_genes_from_gff.pl $gff $outdir/gene_from_gff.info\n";
$anno .= "awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' >  $outdir/final.info\n";
## dmr2gene,select genes
$anno .= "$perl $script/dmr2gene.pl $gff $dmrfile $outdir/dmr2gene.txt\n";
$anno .= "$perl $script/select_gene.pl $outdir/final.info $outdir/dmr2gene.txt $outdir/gene_select.txt\n";

## anno
$anno .= "$perl $script/01.go_function_enrichment.pl $outdir/go.anno $outdir/gene_select.txt $outdir/out.go_function $outdir/out.go_enrichment\n";
$anno .= "$perl $script/02.keggpathway.pl $outdir/$species.ko $outdir/$species.kegg.path $outdir/gene_select.txt $outdir/out.keggpathway\n";

#system("$anno");
print $anno;
