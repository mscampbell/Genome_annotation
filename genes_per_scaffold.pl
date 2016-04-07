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

my $data = parse($FILE);
report($data);
#PostData($data);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    foreach my $key (keys %$data){
	print $key."\t". $data->{$key}."\n";;
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my %data;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'gene'){
	    $data{$cols[0]}++;
	}
	else{next;}

    }
    return (\%data);
    $fh->close();
}
#-----------------------------------------------------------------------------

