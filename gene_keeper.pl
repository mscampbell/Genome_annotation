#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_s);
getopts('iegpcms');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
takes a list of gene ids and prints the genes from a gff3_file
gene_keeper.pl <id_file.txt> <gff3_file.gff>
\n\n";

my $FILE1 = $ARGV[0];
my $FILE2 = $ARGV[1];
die($usage) unless $ARGV[1];

my %LU_G;
my %LU_T;

build_lu_gid($FILE1);
build_lu_tid($FILE2);
filter($FILE2);
#PostData(\%LU_T);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub filter{
    my $file = shift;

    my $fh = new FileHandle;
    $fh->open($file);
    print "##gff-version 3\n";
    while (defined(my $line = <$fh>)){
        chomp($line);

        last if $line =~ /^\#\#FASTA/;
        next if $line =~ /^\#/;
        my @array = split(/\t/, $line);
	
        if ($array[2] eq 'gene'){
	    my ($id) = $array[8] =~ /ID=(\S+?);/;
	    print $line."\n" if defined($LU_G{$id}); 
	}
	elsif ($array[2] eq 'mRNA' || $array[2] eq 'transcript'){
	    
	    my ($id) = $array[8] =~ /ID=(\S+?);/;
	    $line =~ s/\ttranscript\t/\tmRNA\t/ if $opt_s;
	    print $line."\n" if defined($LU_T{$id}); 
	}
	elsif ($array[2] eq 'exon'||
	    $array[2] eq 'CDS'||
	    $array[2] eq 'three_prime_UTR'||
	    $array[2] eq 'five_prime_UTR'){
	    my $bool = 0;
	    if ($array[8] =~ /Parent=(\S+?);/ || $array[8] =~ /Parent=(\S+)$/){
		my $ids = $1;
		my @ids_array = split(/,/, $ids);
		foreach my $x (@ids_array){
		    if (defined($LU_T{$x})){
			$bool++;
		    }	
		}
	    }
	    print $line."\n" if $bool;
	}
	
    }
    $fh->close();
    
}
#-----------------------------------------------------------------------------
sub build_lu_gid{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	$LU_G{$line}=1;
    }
    $fh->close();
}
#-----------------------------------------------------------------------------
sub build_lu_tid{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	
	last if $line =~ /^\#\#FASTA/;
        next if $line =~ /^\#/;
        my @array = split(/\t/, $line);

        if ($array[2] =~ 'mRNA' || $array[2] eq 'transcript'){
	    my ($tid) = $line =~ /ID=(.+?);/;
            my ($gid) = $line =~ /Parent=(.+?);/;
	    if (defined($LU_G{$gid})){
		$LU_T{$tid}=1;
	    }
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

