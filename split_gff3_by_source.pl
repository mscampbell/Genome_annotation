#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{
    
    my $file = shift;       
    my %out_files;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	chomp($line);
	my @col = split(/\t/, $line);
	
	my $filename = $col[1];
	$filename =~ s/:/_/;
	$filename .= '.gff';
	
	if (!defined($out_files{$filename})){
	    open (my $fh2, ">>", $filename);
	    print $fh2 $line ."\n";
	    $out_files{$filename}=$fh2;
	}
	else {
	    my $fh2 = $out_files{$filename};
	    print $fh2 $line ."\n";
	}
    }
    $fh->close();
    foreach my $file (keys %out_files){
	my $fh2 = $out_files{$file};
	close($fh2);
    }
}
#-----------------------------------------------------------------------------

