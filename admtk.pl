#written by Jonathan M. Flowers
#! usr/bin/env perl -w
use strict;
use warnings;
use Getopt::Long;
use Digest::MD5;

&main;
exit;

sub main {

&usage if(@ARGV < 1);
my $command = shift(@ARGV); 
my %func = (checksums=>\&checksums);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}};

}

sub checksums {

my $usage = '

Usage:		admtk.pl checksums [file 1] [file 2] ... [file n]

Description:	Checks the md5 sum of specified file(s) (e.g. filename) against
		the md5 sum in a file with the name filename.md5

Options:	none

Note:		contenst of filename.md5 must have one line matching the 
		regex /^\S+\s+\S+$/.
';

die("$usage\n") unless(@ARGV > 0);

my @files = @ARGV;
map { print "error: $_ not found" unless(-e $_) } @files;
map { print "error: $_.md5 not found" unless(-e $_.'.md5') } @files;
    
foreach my $file (@files){
	my ($digest) = _md5_checksum($file);
	my ($source_digest) = _read_md5_file($file);
	my $match = ($digest eq $source_digest) ? 1 : 0;
	if($match){
		print "$file md5sum validated\n";
	}else{
		print "*** error: $file md5sum did not validate ***\n"; exit;
	}
}

}

sub _md5_checksum {

my $file = shift @_;
open(FILE, $file), or die $!;
binmode(FILE);
my $digest = Digest::MD5->new->addfile(*FILE)->hexdigest;
close(FILE);

return ($digest);

}

sub _read_md5_file {
    
my $file = shift @_;
open(MD5FILE,$file.'.md5'), or die $!;
$_ = <MD5FILE>;
close(MD5FILE);
my @source_digest = split;

die("error: unrecognized format in file $file.md5\n") if(@source_digest != 2);
die("error: unrecognized format in file $file.md5\n") if($source_digest[1] ne $file);

return $source_digest[0];

}

sub usage {

die(qq/
Version: 1.1 (4.27.12)
Usage:   admtk.pl <command> [options]
Command: checksums	check md5sums

\n/);

}
