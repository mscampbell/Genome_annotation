#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;
use PostData;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
my %DATA;
my %PRINTED;
parse($FILE);
report();
#PostData(\%DATA);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    foreach my $gid (keys %DATA){
	foreach my $tid (sort {$DATA{$gid}{$b} <=> $DATA{$gid}{$a}  } keys %{$DATA{$gid}}){
	    print "$tid\n"unless defined($PRINTED{$gid});
	    $PRINTED{$gid} = 1;
	}
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my %lu;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;

	my @col = split(/\t/, $line);
	if ($col[2] eq 'mRNA'){
	    my ($t) = $col[8] =~ /ID=(\S+?);/;
	    my ($g) = $col[8] =~ /Parent=(\S+?);/;
	    $lu{$t} = $g; #this only workd if the mRNA lines are bfore the cds lines
	}
	if ($col[2] eq 'CDS' && ($col[8] =~ /Parent=(\S+?);/ || $col[8] =~ /Parent=(\S+)/)){
	    my $id = $1;
	    my $gid = $lu{$id};
	    #if($col[8] =~ /Parent=(\S+?)-mRNA/){
	    #	($gid)  = $col[8] =~ /Parent=(\S+?)-mRNA/;
	    #}
	    #elsif ($col[8] =~/Parent=(transcript:A\S+?);/){
	    #	
	    #	($gid) = $col[8]=~/Parent=(transcript:A\S+?);/;
	    #	$gid =~ s/transcript/gene/;
	    #	$gid =~ s/_FGT/_FG/;
	    #}
	    #elsif ($col[8] =~ /Parent=(\S+?)-/){
	    #	($gid)  = $col[8] =~ /Parent=(\S+?)-/;
	    #}
	    #elsif ($col[8] =~ /Parent=(\S+?)\_/){
	    #	($gid)  = $col[8] =~ /Parent=(\S+?)\_/;
	    #}
	    #else {die "this happened $col[8]\n";}
	    my $len = $col[4] - $col[3];
	    $DATA{$gid}{$id} += $len;
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

