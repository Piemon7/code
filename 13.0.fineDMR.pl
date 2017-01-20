#!/usr/bin/perl -w 
#piemon^2016.6.22
use strict;

die "Usage: perl $0 <CpG1> <CpG2> <outdir> 
INFO:	get dmr region. 
	maybe should merge CG,CHG,CHH first.
	it takes lots of time for large file.
\n\n" unless @ARGV == 3;

my $commond = ();
$commond .= "setwd(\"$ARGV[2]\")\nlibrary(methylKit)\n";
$commond .= "file.list=list(\"$ARGV[0]\",\"$ARGV[1]\")\n";
$commond .= "myobj=read(file.list,sample.id=list(\"test1\",\"test2\"),assembly=\"hg19\",treatment=c(1,0))\n";
$commond .= "meth=unite(myobj,destrand=FALSE)\nmyDiff=calculateDiffMeth(meth,mc.cores=6)\n";
$commond .= "write.table(as.data.frame(myDiff),file=\"myDiff.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.01)\n";
$commond .= "write.table(as.data.frame(myDiff25p),file=\"myDiff25p.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "myDiff25pHypo=get.methylDiff(myDiff,difference=25,qvalue=0.01,type=\"hypo\")\n";
$commond .= "write.table(as.data.frame(myDiff25pHypo),file=\"myDiff25pHypo.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "myDiff25pHyper=get.methylDiff(myDiff,difference=25,qvalue=0.01,type=\"hyper\")\n";
$commond .= "write.table(as.data.frame(myDiff25pHyper),file=\"myDiff25pHyper.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "library(edmr)\nlibrary(GenomicRanges)\nlibrary(IRanges)\nlibrary(mixtools)\nlibrary(data.table)\n";
$commond .= "mydmr=edmr(myDiff, mode=1, ACF=TRUE)\nwrite.table(as.data.frame(mydmr),file=\"DMR.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "mysigdmr=filter.dmr(mydmr)\nwrite.table(as.data.frame(mysigdmr),file=\"sigDMR.txt\",sep = \"\\t\", row.names = FALSE)\n";

$commond .= "mydmr2=edmr(myDiff, mode=2, ACF=TRUE)\nwrite.table(as.data.frame(mydmr2),file=\"DMR2.txt\",sep = \"\\t\", row.names = FALSE)\n";
$commond .= "mysigdmr2=filter.dmr(mydmr2)\nwrite.table(as.data.frame(mysigdmr2),file=\"sigDMR2.txt\",sep = \"\\t\", row.names = FALSE)\n";

my $outfile = "$ARGV[2]/edmr.R";
open OUT,'>',$outfile or die $!;
print OUT $commond;
#system("Rscript $outfile");
