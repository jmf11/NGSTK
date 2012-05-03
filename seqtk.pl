#! usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure('bundling');


&main;
exit;

sub main {

&usage if(@ARGV < 1);
my $command = shift(@ARGV); 
my %func = (seq2gc=>\&seq2gc,subseq=>\&subseq,fq2fa=>\&fq2fa,seq2len=>\&seq2len);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}};

}

sub seq2len {
my $usage = '
Usage:       seqtk.pl seq2len [file1.fas] [file2.fas] ... [fileN.fas]

Description: Return lengths of sequences in fasta format.
             Default output is a tab-delimited, 2 column table
             with seq id and length

Options:     --bed return BED3-formatted file
';



}
sub fq2fa {

#checked 2.3.12, all correct	
my $usage = '

Usage:       seqtk.pl fq2fa -i <fastqfile>

Description: Converts fastq to fasta. Fastq can be either
             standard 4-line (id,seq,+,qual) or wrapped 
             (ie, id,multi-line seq, +, multi-line qual) version.
             
	         Requires that line with + must match /^\+$/, that is,
             must not contain additional characters as was typical
             of some early fastq files.

             Writes to stdout

Options:     -i <fastqfile>

Note:        fq2fa uses the parsing strategy suggested by
             Cock et al. 2010 to account for the fact that
             @ and + are not restricted to the identifier
             and separator lines. 

Example:
cat tmp.fq
@chr1
AAAAAAA
+
#$%@#$$
@chr2
TTTTTTT
TTTTTTT
+
#%^$%^$
@@@@@@@
@chr3
CCCCCCCCCC
+
helloworld

perl seqtk.pl fq2fa -i tmp.fq
>chr1
AAAAAAA
>chr2
TTTTTTT
TTTTTTT
>chr3
CCCCCCCCCC
';

my %args = ();
GetOptions(\%args,'in=s');

die($usage) if(!defined($args{in}));
die("error: $args{in} not found") unless(-e $args{in});

my $t = 3; #$t == toggle
my %id2qlen = ();
my %id2slen = ();
my $id = '';

open F,$args{in}, or die $!;
while(<F>){
	chomp;
	next if(/^$/);

	if($t == 1){
		if(/^\+$/){
			$t = 2;
		}else{
			print $_,"\n"; 
			$id2slen{$id} += length($_);
		}
	}elsif($t == 2){
		$id2qlen{$id} += length ($_);
		$t = 3 if($id2qlen{$id} == $id2slen{$id});
	}elsif($t == 3){
		if(/^@(.*)$/){
			$id = $1;
			if( defined($id2slen{$id}) or defined($id2qlen{$id}) ){
				die("error: sequence $id appears multiple times\n");
			}
			print '>'.$id,"\n"; 
			$t = 1;
		}else{
			die("error: first character in line 1 should be @ symbol\n") if($. == 1);
			die("error: first character in line line after qual block isnt @\n");
		}
	}else{
		die("error: unexpected error\n");
	}
}
close F;

if(scalar(keys %id2slen) != scalar(keys %id2qlen)){
	die("error: unexpected error\n");
}

}

sub seq2gc {

#checked carefully, all correct 2.1.12
my $usage = (qq/
'
Usage:       seqtk.pl seq2gc <fasta file>
              or 
             seqtk.pl subseq [options] | perl seqtk.pl seq2gc -

Description: Determines the GC content of a set of sequences
             in fasta format (either from stream or input file)

             Writes to stdout

Options:     none

Note: *Amibiguous and heterozygous sites ignored.
      *Sequences can be lowercase, uppercase or mixed

Example:
>seq1
ATGC
GGGG
>seq2
ATNN
NNGC
>seq3
atnn
nngc
>seq4
ATYN
Gtra

perl seqtk.pl seq2gc tmp.fas
seq1	0.75
seq2	0.5
seq3	0.5
seq4	0.2
'
\n/);

die($usage) if (@ARGV == 0 && -t STDIN);

my $id = '';
my @order = ();
my %id2gc = ();
my %id2tot = ();
while(<>){
	chomp;
	next if(/^$/);
  if(/^>(.*)$/){
		$id = $1;
		die("error: duplicate fasta sequence headers\n") if(exists($id2gc{$id}));
		$id2gc{$id} = 0;
    $id2tot{$id} = 0;
		push @order,$id;
	}else{
    my $seq = uc $_; #upper case all sequence first
		$id2gc{$id} += $seq =~ tr/GC//;
	  $id2tot{$id} += $seq =~ tr/[ATGC]//;
	}
}

foreach my $id (@order){
	if($id2tot{$id}==0){
		print join("\t",$id,"NA"),"\n";
	}else{
		print join("\t",$id,$id2gc{$id}/$id2tot{$id}),"\n"; 
	}
}

}

sub subseq {

#checked to here, all correct 2.3.12
my $usage = (qq/
'
Usage: 	     perl seqtk.pl subseq [<arguments>]

Description: subseq extracts DNA sub-sequences from a fasta file.
             Sequences can be extracted by specifying a single id
             (using options --id,--start,--end, --in) or in batch using a 
             BED3-formatted input file (--bed, --in).

             Writes to stdout.	
          
Options:     --id       [STR] sequence id
             --start 	[INT] start coordinate (0-based, unless -1 is set)
             --1 	      1-based start coordinate (default= 0-based)
             --end	[INT] end coordinate (1-based)
             --in 	[STR] fasta input file
             --bed 	[STR] tab-delmited BED3 input file without a header. 
                              Three columns are:
                              (1) sequence id, 
                              (2) 5 prime coordinate (0-based),
                              (3) 3 prime coordinate (1-based)
                              **--bed incompatible with --id,--start, --end, --1 

Notes: *coordinates following ":" in output header all 1-based
       *use quotes if fasta header has whitespace (e.g., --id "my id")
       *--start should be < --end (or --start <= --end if --1). 
       *Revcomp not implemented for --start > --end 
       *number of characters per line from input is preserved. 
       
Example:
cat tmp.fas
>chr2
AATT
GGCC
>chr3
TTGG
CCAA

perl seqtk.pl subseq --id chr2 --start 4 --end 8 --in tmp.fas
>chr2:5..8
GGCC

perl seqtk.pl subseq --id chr2 --start 5 --end 8 --in tmp.fas --1
>chr2:5..8
GGCC

cat tmp.bed
chr2	0	6
chr2	3	7
chr3	2	6

perl seqtk.pl subseq --bed tmp.bed --in tmp.fas
>chr2:1..6
AATT
GG
>chr2:4..7
TGGC
>chr3:3..6
GGCC

'
\n/);
 
die($usage) if(@ARGV < 1);

my %args = (1=>0);
GetOptions(\%args,"id=s","start=i","end=i","in=s","bed=s","1");

die("error: file $args{in} not found\n") unless(-e $args{in});

#get hash of input sequences (key=id,value=seq)
#*retrieve linewidth from original input file so it can be preserved in output
my ($id2seq,$linewidth) = _getid2seq($args{in}); #%id2seq contains all sequences in $args{i} (could optimize)
my %id2seq = %{$id2seq};

if(defined($args{id})){

	#command line input feature [could be 0-based or 1-based start coordinate]

	die($usage) unless(defined($args{start}) and defined($args{end}));
	die($usage) if(defined $args{bed});
	die("error: $args{id} not found in $args{in}\n") unless(exists($id2seq{$args{id}}));

	$args{start} = ($args{start} - 1) if($args{1}==1); #if start coord 1-based, change to 0-based

	my $header = '>'.$args{id}.':'.($args{start} + 1).'..'.$args{end};
	my ($subseq) = _extractsubseq($args{start},$args{end},$id2seq{$args{id}}); 
	&_printseq($header,$subseq,$linewidth);

}elsif(defined $args{bed}){ 

	#BED input [start coord already 0-based, end coord 1-based]

	die($usage) if(defined($args{id}) or defined($args{start}) or defined($args{end}));
	die("error: --1 and --bed cant both be set (BED-file start coords must be 0-based)\n") if($args{1}==1);
	die("error: file $args{bed} not found\n") unless(-e $args{bed});
    
	open F,$args{bed}, or die $!;
	while(<F>){
		chomp;
		next if(/^$/);
		die("error: BED not tab-delimited\n") unless(/\t/);

		my @s = split /\t/;
		die("error: file specified in --bed is not BED3\n") if(@s != 3); #howl if not 3 columns
		
		my $f = $s[0];
		my $s = $s[1]; #0-based (per BED specification)
		my $e = $s[2]; #1-based (per BED specification)

		die("error: found white space in BED column 1 id\n") if($f =~ /\s+/);
		die("error: found non-digit in BED column 2 or 3\n") if($s =~ /[^0-9]/ or $e =~ /[^0-9]/);
		die("error: $f not found in $args{in}\n") unless(exists($id2seq{$f}));
       
		my $header = '>'.$f.':'.($s + 1).'..'.$e;
		my ($subseq) = _extractsubseq($s,$e,$id2seq{$f});

		&_printseq($header,$subseq,$linewidth);
	}
	close F;
}else{
print 'uhoh',"\n"; exit; 
	die($usage);
}
}

sub _extractsubseq {

my ($s,$e,$seq) = @_;	#$s = 0-based start
											#$e = 1-based end
my $len = length ($seq); 
my $sublen = $e - $s; 	#$s is zero-based, $e is 1-based
die("error: cant extract sequence of length < 1\n")if($sublen < 1);
die("error: -e $e is beyond end of sequence\n")if( ($s + $sublen) > $len); #correct
my $subseq = substr $seq, $s, $sublen; #zero-indexed offset
return($subseq);	

}

sub _printseq {
	my ($header,$seq,$linewidth) = @_;
	$seq =~ s/(\w{1,$linewidth})/$1\n/g; #add eol character (\n) for writing output 
	print join("\n",$header,$seq);
}

sub _getid2seq {
	my ($file) = @_;

	my %id2seq = ();
	my $id = '';
	my %linewidth = ();
	open F,$file, or die $!;
	while(<F>){
		chomp;
		next if(/^$/);
        if($. == 1){
			die("error: $file is not in fasta format\n")unless(/^>/); 
		}

		if(/^>(.*)$/){
			$id = $1;
			die("error: multiple entries with same id in $file\n") if(exists $id2seq{$id});
			$id2seq{$id} = '';
		}else{
			$id2seq{$id} .= $_;
			$linewidth{length ($_) }++;
		}
	}
	close F;

	my @linewidth = sort { $b <=> $a } keys %linewidth;
	return(\%id2seq,$linewidth[0]);
}

sub usage {
die(qq/
Version: 1.1 (2.3.12)
Usage:   seqtk.pl <command> [<arguments>]\n
Command: seq2gc           get GC content
         subseq           get sub-sequence(s)
         fq2fa            convert fastq to fasta
\n/);
}

