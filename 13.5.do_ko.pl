#!/usr/bin/perl -w 
#piemon^2016.6.30
use strict;
use File::Basename;

die "usage:perl $0 <ko> <dmr2gene.out> <komap> <outdir>
	komap: one of animal,plant,prokaryote,fungi,microorganism,other.
\n\n" unless @ARGV == 4;

my $name = basename $ARGV[1];
my $bg = $ARGV[0];
my $outdir = $ARGV[3];
my $glist = "$outdir/$name.glist";
my $fg    = "$outdir/$name.ko";

my $komap = "/nfs/database/db/Pub/kegg/RNA/59.3/komap/$ARGV[2]\_ko_map.tab";

open IN,$ARGV[1] or die $!;
open OUT,'>',$glist or die $!;
print OUT "GeneID\n";
while(<IN>){
	next if (/^#/);
	chomp;
	my $last = (split "\t",$_)[-1];
	print OUT $last,"\n";
}

my $commond = ();
$commond .= "path=/nfs/pipe/RNA/RNA-ref/version1/functional\n";
$commond .= "export LD_LIBRARY_PATH=/nfs/config/boost/lib:\$LD_LIBRARY_PATH\n";
$commond .= "bg=$bg\n";
$commond .= "glist=$glist\n";
$commond .= "perl /nfs/pipe/RNA/RNA-ref/version1/functional/getKO.pl -glist $glist -bg $bg -outdir $outdir\n";
$commond .= "perl /nfs/pipe/RNA/RNA-ref/version1/functional/pathfind.pl -fg $fg -komap $komap -bg $bg -output $outdir/$name.path\n";
$commond .= "perl /nfs/pipe/RNA/RNA-ref/version1/functional/pathway_enrichFigure.pl $outdir\n";

#print $commond;
system ("$commond");
