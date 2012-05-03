#! usr/bin/env perl -w
use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {

&usage if(@ARGV < 1);
my $command = shift(@ARGV); 
my %func = (covbedhist2sites=>\&covbedhist2sites,
            covbedhist2cov=>\&covbedhist2cov,
            getbedintervals=>\&getbedintervals,
            gff2bed=>\&gff2bed);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}};

}

sub gff2bed {

#checked 1.30.12, all correct

my $usage = (qq/
'
Usage:       bedtk.pl gff2bed <gff> [options]

Description: extracts features of specified type from gff
             and returns BED3 file with chromosome, 
             and coordinates for all features of that type 
       
             Types are typically gene, cds, mrna etc. Start 
             coordinate is always less than end coordinate 
             in output as in gff.  

	     Output (tab-delimited):
             col1=chromosome
             col2=5p coordinate (0-based)
             col3=3p coordinate (1-based)

             Writes to stdout

Options:     -t feature type 
             -i <gff>

Note: Output file coordinates are 0-based for 5' coordinates
      and 1-based for 3' coordinates.
'
\n/);

my %args = ();
getopts('t:i:',\%args);

die($usage) unless(defined($args{t}) and defined($args{i}));
die("error: file $args{i} not found\n") unless(-e $args{i});
 
open F,$args{i}, or die $!;
while(<F>){
  chomp;
  next if(/^$/);
  next if(/^#/);

  my @s = split /\t/;
  my $chr = $s[0];
  my $feat = $s[2];
  my $fivep = $s[3];
  my $threep = $s[4];
 
  print join("\t",$chr,$fivep - 1, $threep),"\n" if($feat eq $args{t});
}
close F;

}


sub covbedhist2cov {

#covbedhist2cov output checked carefully, all correct (1.29.12)

my $usage = (qq/
'
Usage:   bedtk.pl covbedhist2cov <coverageBed -hist output file>
	 or
	 coverageBed -hist -abam <bam file> -b <bedfile> | bedtk.pl covbedhist2cov -
	 
Description: Calculates coverage depth in intervals. Writes to stdout 
             tab-delimited feature, 5' boundary, 3' boundary, coverage depth.

Options: none

Note: using in pipe requires bedtools installation
'
\n/); 

die($usage) if (@ARGV == 0 && -t STDIN);

#extract data blocks, one block for each interval in original bedfile
my %int2hist = ();
my @intorder = ();
while(<>){
	chomp;
	my @s = split;
 	
	next if(/^$/);
	next if(/^all\s+/); #coverageBed -hist by default will 
                                    #have coverage for all intervals, ignore this
	die("error: unexpected number of columns in coverageBed -hist output\n") unless(@s==7);
		
	if (/^(\S+\s+\d+\s+\d+)/){
		push @{$int2hist{$1}}, $_;
		push @intorder,$1 if(@intorder==0 or $intorder[-1] ne $1); #store input interval order	
	}else{
		die("error: unexpected coverageBed -hist output format\n");	
	}
}

#iterate over each interval
foreach my $interval (@intorder){

	#iterate over each bin in histogram
	my %intsize = (); #="interval size"
	my $sum = 0;
	foreach my $bin (@{$int2hist{$interval}}){
		my @s = split /\s+/,$bin;
		$sum += ($s[3]*$s[4]);
		$intsize{$s[5]}++;
	}

	#get coverage depth for interval
	die("error: unexpected window size in coverageBed -hist output\n") unless(scalar keys %intsize == 1);
	my @intsize = keys %intsize;
	my $cov = $sum/$intsize[0];
	print $interval,"\t",$cov,"\n"; 
}
}

sub getbedintervals {

#checked 1.13.2012, all correct
#changed flag options from -b (start) to -s, 
#                          -s (step) to -t
#                          -i (input bed) to -b, checked 2.1.12 allcorrect
                     
my $usage = (qq/
'
Usage: 	     perl bedtk.pl getbedintervals [<arguments>]

Description: getbedintervals generates a file in BED3* format 
             with coordinates for sliding window analysis of a region.  

             The input is the feature name (e.g., chr10), the 
             0-based start, and 1-based end coordinates of the feature 
             [ie, the region or chromosome within which windows are to be 
             defined], the window size and step size.
              
             Alternatively, a BED3 file with feature name, 
             0-based start coordinate and 1-based end coordinate can be 
             specified, if intervals desired for many features (e.g.,
             all chromosomes).

             Output is BED3 (0-based start coords and 1-based end).

             Writes to stdout     

Options:     -f feature name*
             -w window size
             -t step size
             -s start coordinate (0-based)  
             -e end coordinate (1-based) 
             -b <input bedfile>

Note: *the output file strictly meets UCSC BED3 definition only if the feature 
      specified by -f is a chromosome. Any feature name is acceptable for 
      this script provided that no white space in the name.

Example: perl bedtk.pl getbedintervals -f chr1 -s 0 -e 50 -t 5 -w 10
chr1	0	10
chr1	5	15
chr1	10	20
chr1	15	25
chr1	20	30
chr1	25	35
chr1	30	40
chr1	35	45
chr1	40	50
'
\n/);

	die($usage) if(@ARGV < 1);
 
	my %args = ();
	getopts('f:w:t:s:e:b:',\%args);
	die($usage) if(!defined($args{w}) or !defined($args{t}));
	die("error: stepsize > windowsize which is unconventional\n") if($args{t} > $args{w});

	if(defined($args{b})){
		die("filename in -b does not exist\n") unless(-e $args{b});
		open F,$args{b}, or die $!;
		while(<F>){
			chomp;
			next if(/^$/);

			if(/^(\S+)\t(\d+)\t(\d+)$/){
				&_getinterval($1,$2,$3,$args{w},$args{t});
			}else{
				die("error: check that input file is in BED3 format\n");
			}
		}
		close F;
	}else{
		die($usage) unless( defined($args{f}) and defined($args{s}) and defined($args{e}) );                 
		&_getinterval($args{f},$args{s},$args{e},$args{w},$args{t});
	}
}

sub _getinterval {

	my ($feat,$s,$e,$windowsize,$stepsize) = @_;

	die("error: windowsize [$windowsize] larger than feature [$s,$e]\n") if($s + $windowsize > $e);
        die("error: start coord ($s) must be less than end ($e)\n") if($s > $e);

	while( ($s + $windowsize) <= $e){
		print join ("\t",$feat,$s,$s + $windowsize),"\n"; 
		$s += $stepsize;
	}
}


sub covbedhist2sites {

#checked carefully 3.6.12 -> all correct
my $usage = (qq/
'
Usage:       bedtk.pl covbedhist2sites [options] <coverageBed -hist output file>
              or
             coverageBed -hist -abam <bam file> -b <bedfile> | bedtk.pl covbedhist2sites [options] -
	 
Description: Calculates the number of sites with less than, equal to, or greater than a specified 
             coverage depth in coverageBed -hist stdout or output file.
             
             Writes to stdout

Options:     -l [depth]  number of sites less than specified depth 
             -e [depth]  number of sites equal to specified depth
             -g [depth]  number of sites greater than specified depth

Note: using in pipe requires bedtools installation

Example:
cat tmp.bed
chromosome_1	0	10000
chromosome_1	10000	20000

coverageBed -hist -abam tmp.bam -b tmp.bed > tmp.covhist

cat tmp.covhist
chromosome_1	10000	20000	3	3	10000	0.0003000
chromosome_1	10000	20000	4	18	10000	0.0018000
chromosome_1	10000	20000	5	13	10000	0.0013000
chromosome_1	0	10000	4	2	10000	0.0002000
chromosome_1	0	10000	5	1	10000	0.0001000

perl bedtk.pl covbedhist2sites -l 5 tmp.covhist
chromosome_1	10000	20000	21
chromosome_1	0	10000	2

coverageBed -hist -abam tmp.bam -b tmp.bed | perl bedtk.pl 
covbedhist2sites -l 5 -  
chromosome_1	10000	20000	21
chromosome_1	0	10000	2
'
\n/); 

die($usage) if (@ARGV == 0 && -t STDIN);

my %args = ();
getopts('l:e:g:',\%args);

die($usage) unless(exists($args{l}) or exists($args{e}) or exists($args{g}));
die($usage) unless(scalar (keys %args) == 1);

#extract data blocks, one block for each interval in original bedfile
my %int2hist = ();
my @intorder = ();
while(<>){
	chomp;
	my @s = split;
 	
	next if(/^$/);
	next if(/^all\s+/); #coverageBed -hist by default will 
                            #have coverage for all intervals, ignore this
	die("error: unexpected number of columns in coverageBed -hist output\n") unless(@s==7);
		
	#expected columns in output for coverageBed -hist:
        #feature 5' 3' coveragedepth numberofsites@covdepthlevel totalsitesinterval
	if (/^(\S+\s+\d+\s+\d+)/){
		push @{$int2hist{$1}}, $_;
		push @intorder,$1 if(@intorder==0 or $intorder[-1] ne $1); #store input interval order	
	}else{
		die("error: unexpected coverageBed -hist output format\n");	
	}
}

die("error: unknown error reading input\n") unless( scalar(keys %int2hist) == scalar(@intorder) );

#iterate over each interval
foreach my $interval (@intorder){

	#iterate over each bin in histogram
	my %intsize = (); #="interval size"
	my %depth2sites = ();
	my $totsites = 0;

	#expected columns in output for coverageBed -hist:
        #feature 5' 3' coveragedepth numberofsites@covdepthlevel totalsitesinterval
	foreach my $bin (@{$int2hist{$interval}}){
		my @s = split /\s+/,$bin;
		my $depth = $s[3];
		my $sitesatdepth = $s[4]; #number of sites at depth == $depth
		$depth2sites{$depth} += $sitesatdepth;

		$totsites += $sitesatdepth;
		$intsize{$s[5]}++;
	}

	my @intsize = keys %intsize;
	die("error: unexpected window size in coverageBed -hist output\n") unless(@intsize == 1);
	die("error: number of sites doesnt sum to interval size\n") unless($intsize[0] == $totsites);

	#count sites meeting criteria defined by -l, -e, or -g for current interval
	my $sitecount = 0; 	
	map { $sitecount += $depth2sites{$_} if($_ < $args{l}) } keys %depth2sites if(defined $args{l});
	map { $sitecount += $depth2sites{$_} if($_ == $args{e}) } keys %depth2sites if(defined $args{e});
	map { $sitecount += $depth2sites{$_} if($_ > $args{g}) } keys %depth2sites if(defined $args{g});

	print join("\t",$interval,$sitecount),"\n";
}
}
sub usage {

die(qq/
Version: 1.2 (3.6.12)
Usage:   bedtk.pl <command> [<arguments>]\n
Command: covbedhist2cov   get the coverage depth from bedtools
                                coverageBed -hist output
         covbedhist2sites get the number of sites that meet a specific
                                coverage depth criterion (e.g., >,=, or < 
                                a specified depth) from coverageBed -hist 
                                output
         getbedintervals  get bed3-formatted file of intervals (e.g.,
                                for sliding window analysis)
         gff2bed          get bed3-formatted coordinates for all 
                                features of a specified type from a gff
\n/);
}







