#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values

unless (@ARGV) {
	die "perl $0 <SNP file> <gff> <ref.fa> <out dir><1 or 0 [1:sample is arabidopsis and gff ID is different with ref ID(Chrchloroplast<->ChrC);0:others]>\n\nNote:if the format of gff file is not gff-version 3(not containing 'ID=' or 'Parent='),you should change the program according to your gff file.\n";
}

my ($Verbose,$Help,$Outdir);
GetOptions(
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);

$Verbose = 1;

my $Outdir=$ARGV[3];
my $species=$ARGV[4];
$Outdir ||= "./";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $snp_file = shift;
my $refGene_file = shift;
my $hs_file = shift;
    

$Outdir =~ s/\/$//;
my $snp_file_basename = basename($snp_file);
$snp_file_basename=~ s/\.xls.gz//g;
my $output_info_file = "$Outdir/$snp_file_basename.CDS.xls";
#my $output_ref_cds_file = "$Outdir/$snp_file_basename.ref.cds";
#my $output_yh_cds_file = "$Outdir/$snp_file_basename.cds";
#my $output_axt_file = "$Outdir/$snp_file_basename.axt";

my %SNP; ##store SNP data
#my %CDS; ##store gene data


my %CODE = (
	'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
	'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
	'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
	'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
	'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
	'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
	'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
	'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
	'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
	'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
	'ATG' => 'M',                                                                         # Methionine
	'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
	'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
	'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
	'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
	'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
	'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
	'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
	'TGG' => 'W',                                                                         # Tryptophan
	'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
	'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
);

my %Abbrev = (
		'A' => [ 'A' ],
		'C' => [ 'C' ],
		'G' => [ 'G' ],
		'T' => [ 'T' ],
		'M' => [ 'A', 'C' ],
		'R' => [ 'A', 'G' ],
		'W' => [ 'A', 'T' ],
		'S' => [ 'C', 'G' ],
		'Y' => [ 'C', 'T' ],
		'K' => [ 'G', 'T' ],
		'V' => [ 'A', 'C', 'G' ],
		'H' => [ 'A', 'C', 'T' ],
		'D' => [ 'A', 'G', 'T' ],
		'B' => [ 'C', 'G', 'T' ],
		'X' => [ 'A', 'C', 'G', 'T' ],  
		'N' => [ 'A', 'C', 'G', 'T' ]
);

## chromosome01    51      T       4543    3.453115        T       A       49      1       3       16.68   1587    3285 
my %capui = ("AC","M", "CA","M", "GT","K", "TG","K", "CT","Y", "TC","Y", "AG","R", "GA","R", "AT","W", "TA","W", "CG","S", "GC","S");
if ($snp_file=~/.gz$/) {
	open(IN,"gzip -dc $snp_file |");
} else {
	open(IN,$snp_file) || die "can't open the $snp_file\n";
}
while(<IN>){
	chomp;
	my @temp=split;
	my $chr = $temp[0];
	if ($species==1) {
		if ($chr=~/mitochondria/) {
        	        $chr=~s/mitochondria/M/g;
	        } elsif ($chr=~/chloroplast/) {
                        $chr=~s/chloroplast/C/g;
                }
	}
	my $pos = $temp[1];
	#my $ref_type = $temp[2];
	my $ref_type = $temp[3]; #shichunwei.

#	my $yh_type = $capui{"$temp[5]$temp[9]"};
#	my $orign_type = "$temp[3]\t$temp[5]\t$temp[9]";
	#my $yh_type = $temp[3];#$capui{"$temp[3]"};
	my $yh_type = $temp[4];#shichunwei
#print"$temp[5] $temp[6] $yh_type\n";
	$SNP{$chr}{$pos} = [$ref_type,$yh_type];
}
close(IN);

print "read snp_file done" if($Verbose);
#print Dumper \%SNP;

my %UTR_3;
my %UTR_5;
my %mRNA;
my %INTRON;
my %CDS;
my %miRNA;
my %miRNA_star;
my %rRNA;
my %TE;
my %tRNA;
my %chromosome;
my %exon;
my %gene;
my %ncRNA;
my %protein;
my %pseudogene;
my %pseudogenic_exon;
my %pseudogenic_transcript;
my %snoRNA;
my %snRNA;
my %transposable_element_gene;
my %transposable_element;
my %transposon_fragment;
my %start_codon;
my %stop_codon;
my %hash_type;
my %transcript;
# chr1    mRNA    rep_CDS 177642  179462  .       +       0       seq_id "AK070557";locus_id "Os01g0103100"
if ($refGene_file=~/.gz$/) {
	open(REF,"gzip -dc $refGene_file |");
} else {
	open(REF,$refGene_file) || die"can't open the $refGene_file\n";
}
while(<REF>){
	chomp;
	my @temp=split(/\t/);
        my $chr=$temp[0];
	my $type=$temp[2];
	my $start=$temp[3];
	my $end=$temp[4];
	my $strand=$temp[6];
	my $info=$temp[8];

	$hash_type{$type}=1;

        my $groups=$temp[8];
        my @groups = split(/\s*;\s*/, $groups);
        my (%groups,$name);
        for my $group (@groups) {
                my ($tag,$value) = split /=/,$group;
                $groups{$tag}=$value;   # patch for those alter-splices
        }
        my @name_order=qw/Parent ID/;
        @name_order=qw/ID Parent/ if $type =~ /mRNA/;
        for (@name_order) {
	        if ($groups{$_}) {$name=$groups{$_};last;}
        }

        my $ID = $name;
	
	if ($type =~ "CDS") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$CDS{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "pseudogenic_transcript") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$pseudogenic_transcript{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "transcript") {
		push @{$transcript{$chr}{$ID}}, [$start,$end,$strand,$info];
	}
	if ($type =~ "snoRNA") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$snoRNA{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "snRNA") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$snRNA{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "transposon_fragment") {
#                 my ($ID) = $temp[8]=~/(\S+);/;
                 push @{$transposon_fragment{$chr}{$ID}}, [$start,$end,$strand,$info];
	}

       if ($type =~ "start_codon") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$start_codon{$chr}{$ID}}, [$start,$end,$strand,$info];

       }

       if ($type =~ "stop_codon") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$stop_codon{$chr}{$ID}}, [$start,$end,$strand,$info];

       }


	if ($type =~ "transposable_element_gene") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$transposable_element_gene{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "transposable_element") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$transposable_element{$chr}{$ID}}, [$start,$end,$strand,$info];

        }



	if ($type =~ "chromosome") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$chromosome{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "exon") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$exon{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "gene") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$gene{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "ncRNA") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$ncRNA{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "protein") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$protein{$chr}{$ID}}, [$start,$end,$strand,$info];

        }
	if ($type =~ "pseudogene") {
#                my ($ID) = $temp[8]=~/(\S+);/;
                push @{$pseudogene{$chr}{$ID}}, [$start,$end,$strand,$info];

        }

        if ($type =~ "pseudogenic_exon") {
#		my ($ID) = $temp[8]=~/(\S+);/;
		push @{$pseudogenic_exon{$chr}{$ID}}, [$start,$end,$strand,$info];

	}
	if($type =~ "three_prime_UTR" || $type =~ "3'-UTR") {
#		my ($ID) = $temp[8]=~/(\S+);/;
		push @{$UTR_3{$chr}{$ID}},[$start,$end,$strand,$info];	
		
	}
	if($type =~ "five_prime_UTR" || $type =~ "5'-UTR") {
#		my ($ID) = $temp[8]=~/(\S+);/;
                push @{$UTR_5{$chr}{$ID}},[$start,$end,$strand,$info];

        }
	if($type =~ "mRNA") {
#		my ($ID) = $temp[8]=~/(\S+);/;
               push @{$mRNA{$chr}{$ID}},[$start,$end,$strand,$info];

        }
	if($type =~ "miRNA") {
#		my ($ID) = $temp[8]=~/(\S+);/;
                push @{$miRNA{$chr}{$ID}},[$start,$end,$strand,$info];

        }
	if ($type =~ "intron") {
#		my ($ID) = $temp[8]=~/(\S+);/;
		push @{$INTRON{$chr}{$ID}},[$start,$end,$strand,$info];
	}
	if ($type =~ "miRNA-star") {
#		my ($ID)=$temp[8]=~/(\S+);/;
		push @{$miRNA_star{$chr}{$ID}},[$start,$end,$strand,$info];
	}
	if ($type =~ "rRNA") {
#		my ($ID)=$temp[8]=~/(\S+);/;
		push @{$rRNA{$chr}{$ID}},[$start,$end,$strand,$info];
	}
	if ($type =~ "TE") {
#		my ($ID)=$temp[8]=~/(\S+);/;
		push @{$TE{$chr}{$ID}},[$start,$end,$strand,$info];
	}
	if ($type =~ "tRNA") {
#                my ($ID)=$temp[8]=~/(\S+);/;
                push @{$tRNA{$chr}{$ID}},[$start,$end,$strand,$info];
        }

}
close(REF);

print "read refGene done" if($Verbose);
##print Dumper \%CDS;


open INFO, "| gzip >$output_info_file.gz" || die "fail $output_info_file";
print INFO "Chr\tpos\tRef_base<->sample_base\tsnp_status\tstrand\tannotation_type\tFeature_type\tcodon_phase(the_phase_of_the_codon_mutation)\tcodon_mutate(the_codon_before_mutation<->the_codon_after_mutation)\tAA_mutate(Amino_acids_before_mutation<->Amino_acids_after_mutation\tSynonymoud_mutations\tnon_synonymous_mutations\tstat_pos\tend_pos\tGene_id_and_its_functions\n";
#open REFCDS, ">".$output_ref_cds_file || die "fail $output_ref_cds_file";
#open YHCDS, ">".$output_yh_cds_file || die "fail $output_yh_cds_file";
#open AXT, ">".$output_axt_file || die "fail $output_axt_file";
#print INFO  "Chromosome-name\tposition\tref_base<->sample_base\tsnp_status\tstrand\tAnnotation_type\tFeature_type\tcodon_phase\tcodon_mutate\tAA_mutate\tsynonymous\tnonsynonymous\tpos_start\tpos_end\tgene_id_and_functions\n";
##read human genome file, we must first get correct codon
##the phase method can't work near the splicing site
##the input must be correct gene models
if ($hs_file=~/.gz$/) {
	open(IN, " gzip -dc $hs_file  |");
} else { 
	open(IN, $hs_file) || die ("can not open $hs_file\n");
}

$/=">"; <IN>; $/="\n";
if (exists $hash_type{"mRNA"}){
                open O,"| gzip >$Outdir/$snp_file_basename.mRNA.xls.gz" or die $!;
                print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
                close O ;
}
if (exists $hash_type{"three_prime_UTR"}){
                open O,"| gzip >$Outdir/$snp_file_basename.three_prime_UTR.xls.gz" or die $!;
                print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
                close O ;
        }
if (exists $hash_type{"five_prime_UTR"}){
                open O,"| gzip >$Outdir/$snp_file_basename.five_prime_UTR.xls.gz" or die $!;
                print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
                close O ;
        }
if (exists $hash_type{"exon"}){
                open O,"| gzip >$Outdir/$snp_file_basename.exon.xls.gz" or die $!;
                print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
                close O ;
        }
if (exists $hash_type{"gene"}){
                open O,"| gzip >$Outdir/$snp_file_basename.gene.xls.gz" or die $!;
                print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
                close O ;
        }
while (<IN>) {
	my $chr = $1 if(/^(\S+)/);
	if ($species==1) {
	        if ($chr=~/mitochondria/) {
                        $chr=~s/mitochondria/M/g;
                } elsif ($chr=~/chloroplast/) {
                        $chr=~s/chloroplast/C/g;
                        }
	}
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$seq = uc($seq); ## all upper case
	$/="\n";
	
	warn "read hs $chr done" if($Verbose);
	
	next if(!exists $SNP{$chr} || !exists $CDS{$chr});

	my $chr_snp_p = $SNP{$chr};
	my $chr_cds_p = $CDS{$chr};

	my $chr_utr_3_p = $UTR_3{$chr};
	my $chr_utr_5_p = $UTR_5{$chr};
	my $chr_mrna_p  = $mRNA{$chr};
	my $chr_intron_p = $INTRON{$chr};
	my $chr_miRNA_p = $miRNA{$chr};
	my $chr_miRNA_star_p = $miRNA_star{$chr};
	my $chr_rRNA_p = $rRNA{$chr};
	my $chr_TE_p = $TE{$chr};
	my $chr_tRNA_p = $tRNA{$chr};
	my $chr_chromosome_p = $chromosome{$chr};
        my $chr_start_codon_p = $start_codon{$chr};
        my $chr_stop_codon_p = $stop_codon{$chr};

	my $chr_gene_p = $gene{$chr};
	my $chr_ncRNA_p = $ncRNA{$chr};
	my $chr_protein_p = $protein{$chr};
	my $chr_pseudogene_p = $pseudogene{$chr};
	my $chr_pseudogenic_exon_p = $pseudogenic_exon{$chr};
	my $chr_pseudogenic_transcript_p = $pseudogenic_transcript{$chr};
	my $chr_transcript_p =$transcript{$chr};
	my $chr_snoRNA_p = $snoRNA{$chr};	
	my $chr_snRNA_p = $snRNA{$chr};
	my $chr_transposable_element_gene_p = $transposable_element_gene{$chr};
	my $chr_exon_p = $exon{$chr};
	my $chr_transposable_element_p = $transposable_element{$chr};
	my $chr_transposon_fragment_p = $transposon_fragment{$chr};
	loop_region($chr_mrna_p,$chr_snp_p,"mRNA",$chr) if (exists $hash_type{"mRNA"});
	loop_region($chr_utr_3_p,$chr_snp_p,"three_prime_UTR",$chr) if (exists $hash_type{"three_prime_UTR"});
	loop_region($chr_utr_5_p,$chr_snp_p,"five_prime_UTR",$chr) if (exists $hash_type{"five_prime_UTR"});
	loop_region($chr_utr_3_p,$chr_snp_p,"3'-UTR",$chr) if (exists $hash_type{"3'-UTR"});
	loop_region($chr_utr_3_p,$chr_snp_p,"5'-UTR",$chr) if (exists $hash_type{"5'-UTR"});
	loop_region($chr_intron_p,$chr_snp_p,"intron",$chr) if (exists $hash_type{"intron"});
	loop_region($chr_miRNA_p,$chr_snp_p,"miRNA",$chr) if (exists $hash_type{"miRNA"});
	loop_region($chr_miRNA_star_p,$chr_snp_p,"miRNA_star",$chr) if (exists $hash_type{"miRNA_star"});
	loop_region($chr_rRNA_p,$chr_snp_p,"rRNA",$chr) if (exists $hash_type{"rRNA"});
	loop_region($chr_TE_p,$chr_snp_p,"TE",$chr) if (exists $hash_type{"TE"});
	loop_region($chr_tRNA_p,$chr_snp_p,"tRNA",$chr) if (exists $hash_type{"tRNA"});	
	loop_region($chr_chromosome_p,$chr_snp_p,"chromosome",$chr) if (exists $hash_type{"chromosome"});
        loop_region($chr_start_codon_p,$chr_snp_p,"start_codon",$chr) if (exists $hash_type{"start_codon"});
        loop_region($chr_stop_codon_p,$chr_snp_p,"stop_codon",$chr) if (exists $hash_type{"stop_codon"});
	loop_region($chr_exon_p,$chr_snp_p,"exon",$chr) if (exists $hash_type{"exon"});
	loop_region($chr_gene_p,$chr_snp_p,"gene",$chr) if (exists $hash_type{"gene"});
	loop_region($chr_ncRNA_p,$chr_snp_p,"ncRNA",$chr) if (exists $hash_type{"ncRNA"});
	loop_region($chr_protein_p,$chr_snp_p,"protein",$chr) if (exists $hash_type{"protein"});
	loop_region($chr_pseudogene_p,$chr_snp_p,"pseudogene",$chr) if (exists $hash_type{"pseudogene"});
	loop_region($chr_pseudogenic_exon_p,$chr_snp_p,"pseudogenic_exon",$chr) if (exists $hash_type{"pseudogenic_exon"});
	loop_region($chr_pseudogenic_transcript_p,$chr_snp_p,"pseudogenic_transcript",$chr) if (exists $hash_type{"pseudogenic_transcript"});
	loop_region($chr_transcript_p,$chr_snp_p,"transcript",$chr) if (exists $hash_type{"transcript"});
	loop_region($chr_snoRNA_p,$chr_snp_p,"snoRNA",$chr) if (exists $hash_type{"snoRNA"});
	loop_region($chr_snRNA_p,$chr_snp_p,"snRNA",$chr) if (exists $hash_type{"snRNA"});
	loop_region($chr_transposable_element_gene_p,$chr_snp_p,"transposable_element_gene",$chr) if (exists $hash_type{"transposable_element_gene"});
	loop_region($chr_transposable_element_gene_p,$chr_snp_p,"transposable_element",$chr) if (exists $hash_type{"transposable_element"});	
	loop_region($chr_transposon_fragment_p,$chr_snp_p,"transposon_fragment",$chr) if (exists $hash_type{"transposon_fragment"});
#	my ($output_info,$output_refcds,$output_yhcds,$output_axt);
		my ($output_refcds,$output_yhcds,$output_axt);
	foreach my $gene_id (keys %$chr_cds_p){
		my $gene_p = $chr_cds_p->{$gene_id};
		my $strand = $gene_p->[0][2];
		my $gene_cds_str;
		my $gene_cds_len = 0;
		my $mRNA_pos = 0;
		my $ref_cds_str;
		my $axt_cds_str;
		@$gene_p = (sort {$a->[0] <=> $b->[0]} @$gene_p);
		for (my $i = 0; $i < @$gene_p; $i++){
			my $cds_p = $gene_p->[$i];
			$gene_cds_str .= substr($seq,$cds_p->[0]-1,$cds_p->[1]-$cds_p->[0]+1);
		}
		$gene_cds_len = length($gene_cds_str);
		Complement_Reverse(\$gene_cds_str) if ($strand eq "-");

		$ref_cds_str = $gene_cds_str;
		$axt_cds_str = $gene_cds_str;
		my (%output_info, $ref_base, $yh_base, $snp_status, $codon_num, $phase_num, $codon_phase_str, $ref_codon, $yh_codon);
		if($strand eq "+"){
			for(my $i = 0; $i < @$gene_p; $i++){
				my $cds_p = $gene_p->[$i];
				for(my $j = $cds_p->[0]; $j <= $cds_p->[1]; $j++){
					$mRNA_pos++;
					next unless (exists $chr_snp_p->{$j});
					$ref_base = substr($seq, $j - 1, 1);
					$chr_snp_p->{$j}[0] =~ tr/agct/AGCT/;
					warn "ref_base error on $chr:$j, please check: $ref_base $chr_snp_p->{$j}[0]" if ($ref_base ne $chr_snp_p->{$j}[0]);
					$yh_base = $chr_snp_p->{$j}[1];
					my $site = $mRNA_pos;
					substr($gene_cds_str, $site - 1, 1) = $yh_base;
					substr($axt_cds_str, $site - 1, 1) = different_base($ref_base,$yh_base);
					$snp_status = hom_het($ref_base, $yh_base);
					($codon_num, $phase_num) = codon_phase($site);
					$ref_codon = substr($ref_cds_str, $site - $phase_num - 1, 3);
					$yh_codon = substr($gene_cds_str, $site - $phase_num - 1, 3);
					warn "$yh_codon has N, please check: $gene_id\t$j\n" && next if $yh_codon =~ /N/;
					my ($codon_mutate_str, $aa_mutate_str, $synonymous, $nonsynonymous) = site_info($yh_codon, $ref_codon, $phase_num, $snp_status);
					push @{$output_info{$j}}, ($ref_base, $yh_base, $snp_status, $site, $phase_num, $codon_mutate_str, $aa_mutate_str, $synonymous, $nonsynonymous,$cds_p->[0],$cds_p->[1],$cds_p->[3]);
				}
			}
			my @k = (sort {$a <=> $b} keys %output_info);
			for(my $key = $k[0]; $key <= $k[-1]; $key++){
				next unless exists $output_info{$key};
				my $key1 = $key - 1;
				my $key2 = $key - 2;
				if($output_info{$key}->[4] == 0){
					next if exists $output_info{$key + 1} || exists $output_info{$key + 2};
					print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
				}elsif($output_info{$key}->[4] == 1){
					next if exists $output_info{$key + 1};
					if(!exists $output_info{$key1}){
						print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}else{
						print INFO "$chr\t$key1/$key\t$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}
				}else{
					if(!exists $output_info{$key1} && !exists $output_info{$key2}){
						print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}elsif(exists $output_info{$key1} && !exists $output_info{$key2}){
						print INFO "$chr\t$key1/$key\t$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}elsif(!exists $output_info{$key1} && exists $output_info{$key2}){
						print INFO "$chr\t$key2/$key\t$output_info{$key2}->[0]<->$output_info{$key2}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key2}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key2}->[3]:$output_info{$key2}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}else{
						print INFO "$chr\t$key2/$key1/$key\t$output_info{$key2}->[0]<->$output_info{$key2}->[1]/$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key2}->[2]/$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key2}->[3]:$output_info{$key2}->[4]/$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}
				}
			}
		}else{
			for(my $i = @$gene_p - 1; $i >= 0; $i--){
				my $cds_p = $gene_p->[$i];
				for(my $j = $cds_p->[1]; $j >= $cds_p->[0]; $j--){
					$mRNA_pos++;
					next unless (exists $chr_snp_p->{$j});
					$ref_base = substr($seq, $j - 1, 1);
					$chr_snp_p->{$j}[0] =~ tr/agct/AGCT/;#shichunwei.
					warn "ref_base error on $chr:$j, please check: $ref_base $chr_snp_p->{$j}[0]" if ($ref_base ne $chr_snp_p->{$j}[0]);
					$yh_base = $chr_snp_p->{$j}[1];
					$ref_base =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					$yh_base =~ tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
					#my $site = $gene_cds_len - $mRNA_pos + 1;
					my $site = $mRNA_pos;
					substr($gene_cds_str, $site - 1, 1) = $yh_base;
					substr($axt_cds_str, $site - 1, 1) = different_base($ref_base,$yh_base,$chr,$gene_id,$j);
					$snp_status = hom_het($ref_base, $yh_base);
					($codon_num, $phase_num) = codon_phase($site);
					$codon_phase_str = "$site:$phase_num";
					$ref_codon = substr($ref_cds_str, $site - $phase_num - 1, 3);
					$yh_codon = substr($gene_cds_str, $site - $phase_num - 1, 3);
					warn "$yh_codon has N, please check: $gene_id\t$j\n" && next if $yh_codon =~ /N/;
					my ($codon_mutate_str, $aa_mutate_str, $synonymous, $nonsynonymous) = site_info($yh_codon, $ref_codon, $phase_num, $snp_status);
					push @{$output_info{$j}}, ($ref_base, $yh_base, $snp_status, $site, $phase_num, $codon_mutate_str, $aa_mutate_str, $synonymous, $nonsynonymous,$cds_p->[0],$cds_p->[1],$cds_p->[3]);
				}
			}
			my @k = (sort {$b <=> $a} keys %output_info);
			for(my $key = $k[0]; $key >= $k[-1]; $key--){
				next unless exists $output_info{$key};
				my $key1 = $key + 1;
				my $key2 = $key + 2;
				if($output_info{$key}->[4] == 0){
					next if exists $output_info{$key - 1} || exists $output_info{$key - 2};
					print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
				}elsif($output_info{$key}->[4] == 1){
					next if exists $output_info{$key - 1};
					if(!exists $output_info{$key1}){
						print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}else{
						print INFO "$chr\t$key1/$key\t$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}
				}else{
					if(!exists $output_info{$key1} && !exists $output_info{$key2}){	
						print INFO "$chr\t$key\t$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}elsif(exists $output_info{$key1} && !exists $output_info{$key2}){
						print INFO "$chr\t$key1/$key\t$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}elsif(!exists $output_info{$key1} && exists $output_info{$key2}){
						print INFO "$chr\t$key2/$key\t$output_info{$key2}->[0]<->$output_info{$key2}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key2}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key2}->[3]:$output_info{$key2}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}else{
						print INFO "$chr\t$key2/$key1/$key\t$output_info{$key2}->[0]<->$output_info{$key2}->[1]/$output_info{$key1}->[0]<->$output_info{$key1}->[1]/$output_info{$key}->[0]<->$output_info{$key}->[1]\t$output_info{$key2}->[2]/$output_info{$key1}->[2]/$output_info{$key}->[2]\t$strand\tGene\tCDS\t$output_info{$key2}->[3]:$output_info{$key2}->[4]/$output_info{$key1}->[3]:$output_info{$key1}->[4]/$output_info{$key}->[3]:$output_info{$key}->[4]\t$output_info{$key}->[5]\t$output_info{$key}->[6]\t$output_info{$key}->[7]\t$output_info{$key}->[8]\t$output_info{$key}->[-3]\t$output_info{$key}->[-2]\t$output_info{$key}->[-1]\n";
					}
				}
			}
		}
#		$output_axt .="$gene_id\_ref&$gene_id\_sample\n".$ref_cds_str."\n".$axt_cds_str."\n\n";
		Display_seq(\$ref_cds_str);
		Display_seq(\$gene_cds_str);
#		$output_yhcds .= ">$gene_id\n".$gene_cds_str;
	}
#	print YHCDS $output_yhcds;
#	print AXT $output_axt;
	
}
close(IN);


close INFO;
#close REFCDS;
#close YHCDS;
#close AXT;

sub site_info{
        my($yh_codon, $ref_codon, $phase_num, $snp_status) = @_;
        my ($tem_ref, $tem_yh);
        my @tem1 = split //, $yh_codon;
        my @tem2 = split //, $ref_codon;
        if ($phase_num == 0){
                $tem_yh = $tem1[1] . $tem1[2];
                $tem_ref = $tem2[1] . $tem2[2];
        }elsif($phase_num == 2){
                $tem_yh = $tem1[0] . $tem1[1];
                $tem_ref = $tem2[0] . $tem2[1];
        }else{
                $tem_yh = $tem1[0] . $tem1[2];
                $tem_ref = $tem2[0] . $tem2[2];
        }
        my ($yh_codon1, $yh_codon2);
        my ($synonymous, $nonsynonymous) = (0,0);
        my ($codon_mutate_str, $aa_mutate_str);
        if($tem_yh ne $tem_ref){
		if ($snp_status =~ /Hom/) {
	                $codon_mutate_str = "$ref_codon<->$yh_codon;";
        	        $aa_mutate_str = "$CODE{$ref_codon}<->COMPLEX;";
               		$synonymous = $nonsynonymous = "NA";
		} elsif ($snp_status =~ /Het-one/) {
			my $base = substr($yh_codon,$phase_num,1);
			my $baseref=substr($ref_codon,$phase_num,1);
			my @all;
		        foreach my $replace (@{$Abbrev{$base}}){
				if ($replace ne $baseref) {				
	                		my $new_codon = $yh_codon;
	                		substr($new_codon,$phase_num,1) = $replace;
                			$codon_mutate_str = "$ref_codon<->$new_codon;";
					$aa_mutate_str = "$CODE{$ref_codon}<->COMPLEX;";
					$synonymous = $nonsynonymous = "NA";
				}
			}
		} else {
			my $base = substr($yh_codon,$phase_num,1);
			my @all;
			foreach my $replace (@{$Abbrev{$base}}){
				my $new_codon = $yh_codon;
				substr($new_codon,$phase_num,1) = $replace;
				push @all,$new_codon;
			}
			$codon_mutate_str = "$ref_codon<->$all[0];$ref_codon<->$all[1];";
			$aa_mutate_str = "$CODE{$ref_codon}<->COMPLEX;";
			$synonymous = $nonsynonymous = "NA";
		}
				
        
	
        }else{
                ($yh_codon1,$yh_codon2) = convert_codon($yh_codon, $ref_codon, $phase_num) if ($snp_status =~ /Het/);
                $yh_codon1 = $yh_codon2 = $yh_codon if ($snp_status =~ /Hom/);
                if ($ref_codon ne $yh_codon1){
                        if ($CODE{$ref_codon} eq $CODE{$yh_codon1}) {
                                $synonymous += 0.5;
                        }else{
                                $nonsynonymous += 0.5;
                        }
			$codon_mutate_str .= "$ref_codon<->$yh_codon1;";
                        $aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$yh_codon1};";
                }
                if ($ref_codon ne $yh_codon2) {
                        if ($CODE{$ref_codon} eq $CODE{$yh_codon2}) {
                                $synonymous += 0.5;
                        }else{
                                $nonsynonymous += 0.5;
                        }
			$codon_mutate_str .= "$ref_codon<->$yh_codon2;" if($yh_codon2 ne $yh_codon1);
                        $aa_mutate_str .= "$CODE{$ref_codon}<->$CODE{$yh_codon2};" if($yh_codon2 ne $yh_codon1);
                }
#                $codon_mutate_str = "$ref_codon<->$yh_codon;";
        }
        return ($codon_mutate_str, $aa_mutate_str, $synonymous, $nonsynonymous);
}

## loop each region
###########################
sub loop_region{
	my $p = shift;
	my $chr_snp_p = shift;
	my $type = shift;
	my $chr = shift;

 	open O,"| gzip >>$Outdir/$snp_file_basename.$type.xls.gz" or die $!;
#	print O "Chr\tPosition\tRef\tSample_base\tFeature_type\tStart_position\tEnd_position\tStrand\tGene_id_and_its_functions\n" ;
#	print O "Chromosome-name\tposition\tref\tsnp\tFeature_type\tStart_pos\tEnd_positon\tStrand\tGene_id_and_functions\n";
	foreach my $gene_id (%$p){

		my $gene_p = $p->{$gene_id};
                my $strand = $gene_p->[0][2];

		for (my $i=0; $i<@$gene_p; $i++) { ##loop for each region
                        my $region_p = $gene_p->[$i];
                        for (my $j=$region_p->[0]; $j<=$region_p->[1]; $j++) {

			if (exists $chr_snp_p->{$j}) {
		
				#my $ref_base = substr($seq,$j-1,1);	
				print O "$chr\t$j\t$chr_snp_p->{$j}[0]\t$chr_snp_p->{$j}[-1]\t$type\t$region_p->[0]\t$region_p->[1]\t$region_p->[2]\t$region_p->[3]\n";



			}

			}


		}


	}

	close O;

}

##convert complex codon into explict codons
#############################################
sub convert_codon{
	my $codon = shift;
	my $ref_codon = shift;
	my $phase = shift;
	my @all;
	
	my $base = substr($codon,$phase,1);
	foreach my $replace (@{$Abbrev{$base}}) {
		my $new_codon = $ref_codon;
		substr($new_codon,$phase,1) = $replace;
		push @all,$new_codon;
	}

	return @all;
}



##caculate codon and phase
#############################################
sub codon_phase {
	my $pos = shift;

	my $phase = ($pos-1)%3;
	my $codon = ($pos-1-$phase)/3 + 1;
	
	##print "$codon,$phase,,\n";
	return ($codon,$phase);

}



sub different_base {
	my ($ref,$yh) = @_;
	my $base;
#print"$yh";exit;
	if ($yh =~ /[ACGT]/) {##add "acgt" by shichunwei.2016.5.30.
		$base = $yh;
	}elsif($yh =~ /[MRWSYK]/){
		$base = $Abbrev{$yh}[0] if($Abbrev{$yh}[0] ne $ref);
		$base = $Abbrev{$yh}[1] if($Abbrev{$yh}[1] ne $ref);
		
	}else{
#print"$yh";exit;
	die "different_base unknown complex character, please check the snp data";
	}
	##print "$ref,$yh,$base;\n";
	return $base;
}

##get the status of a snp
#############################################
sub hom_het {
	my ($ref,$yh) = @_;
	
	return "same" if($ref eq $yh);
	
	my $status;
	

	if ($yh =~ /[ACGT]/) { # add by shichunwei.2016.5.30.
		$status = "Hom";
	}elsif($yh =~ /[MRWSYK]/){
		if ($Abbrev{$yh}[0] ne $ref && $Abbrev{$yh}[1] ne $ref) {
			$status = "Het-two";
		}else{
			$status = "Het-one";
		}
	}else{
		die "hom_het unknown complex character, please check the snp data";
	}

	return $status;
}

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


##complement and reverse the given sequence
#usage: Complement_Reverse(\$seq);
#############################################
sub Complement_Reverse{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
		$$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
		$$seq_p=reverse($$seq_p);  
	}
}
#############################################

