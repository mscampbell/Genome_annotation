#!/usr/bin/perl -w 
use strict;
#use lib ('/home/mcampbell/lib');
#use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
takes a list of transcript ids and prints the complete gene models from a gff3_file
transcript_keeper.pl <id_file.txt> <gff3_file.gff>
\n\n";

my $FILE1 = $ARGV[0]; #id file
my $FILE2 = $ARGV[1]; #gff3 file
die($usage) unless $ARGV[1];

my %LU_G;
my %LU_T;


build_lu_tid($FILE1);
build_lu_gid($FILE2); #the gff3 file gets read twice once to build a gene id
filter($FILE2);       #to transcript id lookup and then once to filter based 
#PostData(\%LU_T);    #on the list of transcript ids in the id file.
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
	elsif ($array[2] eq 'mRNA'){
	    my ($id) = $array[8] =~ /ID=(\S+?);/;
	    print $line."\n" if defined($LU_T{$id}); 
	}
	elsif ($array[2] eq 'exon' ||
	    $array[2] eq 'CDS' ||
	    $array[2] eq 'three_prime_UTR' ||
	    $array[2] eq 'five_prime_UTR'){
	    my @gone;
	    my $bool = 0;
	    my ($ids) = $array[8] =~ /Parent=(\S+);?/; #i'm using an array
		my @ids_array = split(/,/, $ids);      #because some of these 
		foreach my $x (@ids_array){            #features have multiple 
		    if (defined($LU_T{$x})){           #parents 
			$bool++;
		    }
		    else {
			push(@gone, $x);
		    }
		}
	    foreach my $toss (@gone){ #This gets rid of the parent ids that 
		                      #are not in the transcripts to keep list
		if ($line =~ /,$toss,/){ 
		    $line =~ s/,$toss//;
		}
		elsif ($line =~ /,$toss$/){
		    $line =~ s/,$toss$//;
		}
		elsif ($line =~ /=$toss,/){
		    $line =~ s/$toss,//;
		}
	    }
	    print $line."\n" if $bool;
	}
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
	$LU_T{$line}=1;
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
	
	last if $line =~ /^\#\#FASTA/;
        next if $line =~ /^\#/;
        my @array = split(/\t/, $line);

        if ($array[2] =~ 'mRNA'){
	    my ($tid) = $line =~ /ID=(.+?);/;
            my ($gid) = $line =~ /Parent=(.+?);/;
	    if (defined($LU_T{$tid})){
		$LU_G{$gid}=1;
	    }
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

