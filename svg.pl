#!/usr/bin/perl -w 

#piemon^2016.6.2.
use strict;
use SVG;

die "Usage:perl $0 <9kb.summary.report>\n\n" unless @ARGV == 1;

open IN,'<',$ARGV[0] or die $!;
my %hash = ();
while(<IN>){
	chomp;
	my @line = split;#the order is ACGT;
	$hash{$line[0]}{$line[1]} = join "\t",@line[2..5];
}



my $svg = SVG->new(width=>1400, height=>1400);
#$svg->circle(id => "This is a circle", cx => 250, cy => 250, r => 150, fill => "red");
#$svg->rect(x => 150, y => 150, width => 200, height => 200);
$svg->line(x1 => 200, y1 => 200, x2 => 200 , y2 => 1200, stroke=>'black',"stroke-width"=>2);
$svg->line(x1 => 200, y1 => 1200,x2 => 1200, y2=>1200, stroke=>'black',"stroke-width"=>2);
my $bin = 80;

my %col = ('0'=>"red",'1'=>"blue",'2'=>"green",'3'=>"orange");
my %base = ('0'=>"A",'1'=>"C",'2'=>"G",'3'=>"T");

for my $count(1..12){
	my @ACGT = split "\t",$hash{'CG'}{$count};
	my $y_count = 1;
	for (0..3){
		next if ($ACGT[$_] == 0);
		$svg->text(x => 200+$bin*$count, y => 1000*$y_count,width => 500, height => 500,"text-anchor"=>"middle","font-family"=>"Arial","font-size"=>300, "-cdata" => $base{$_},fill=>$col{$_});
		$y_count -= $ACGT[$_];
	}
}



#$svg->text(x => 800, y => 500,"font-family"=>"Arial", "text-anchor"=>"end","font-size"=>200, "-cdata" => "T",fill=>"red");
#$svg->text(x => 600, y => 300, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>200, "-cdata" => "C",fill=>"blue");
#$svg->text(x => 400, y => 600, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>200, "-cdata" => "T",fill=>"orange");
#$svg->title(id=>'document-title')->cdata('This is the title');#

my $out = $svg->xmlify;
open SVGFILE, ">temp.svg";
print SVGFILE $out;
