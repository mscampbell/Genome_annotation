#!/usr/bin/perl -w 
use strict;
#use lib ('/home/mcampbell/lib');
#use PostData;
use FileHandle;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\nnucleotide_replacer.pl: 
converts non A,T,G,C neucleotides to \'Ns\ '
input is a fasta file\n\n
USAGE: nucleotide_replacer.pl <file.fasta>\n\n";

die($usage) unless $ARGV[0];

my $FILE = $ARGV[0];

parse($FILE);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while(defined(my $line = <$fh>)){
	chomp $line;
	if ($line =~ /^>/){
	print "$line\n";
    }
	if ($line =~ /^\w/){
	    $line =~ s/D/N/g;
	    $line =~ s/Y/N/g;
	    $line =~ s/R/N/g;
	    $line =~ s/W/N/g;
	    $line =~ s/M/N/g;
	    $line =~ s/K/N/g;
	    $line =~ s/S/N/g;
	    $line =~ s/H/N/g;
	    $line =~ s/D/N/g;
	    $line =~ s/V/N/g;
	    $line =~ s/B/N/g;
	    $line =~ s/a/A/g;
	    $line =~ s/t/T/g;
	    $line =~ s/c/C/g;
	    $line =~ s/g/G/g;
	    print "$line\n"; 
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

