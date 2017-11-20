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
    my %data;
    my $tot =0;
    my $tot_fun =0;
    my $tot_psu =0;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	
	if ($cols[2] eq 'tRNA'){
	    my ($id) = $cols[8] =~ /ID=(\S+?);/;
	    my @stuff = split(/\-/, $id);
	    my ($aa, $cod) = split(/\_/, $stuff[3]);
	    $data{$aa}{'cod'}{$cod}++;
	    $data{$aa}{'count'}++;
	    
	}
    }
    $fh->close();
    
    foreach my $aa (sort keys %data){
	print "\n", $aa ,"\t $data{$aa}{'count'}\n";
	if ($aa eq 'Pseudo'){
	    $tot_psu += $data{$aa}{'count'};
	}
	else {
	    $tot_fun += $data{$aa}{'count'};
	}

	foreach my $cod (keys %{$data{$aa}{'cod'}}){
	    print "$cod\t$data{$aa}{'cod'}{$cod}\n";
	   
	}
    }
    print "\ntotal_functional_tRNAs\t$tot_fun\n";
    print "total_pseudo_tRNAs\t$tot_psu\n";
}
#-----------------------------------------------------------------------------

