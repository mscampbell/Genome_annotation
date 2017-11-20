#!/usr/bin/perl -w 
use strict;
use Getopt::Std;
use vars qw($opt_d);
getopts('d');

use warnings;
use lib ('/Users/mcampbell/mcampbell/lib');
use PostData;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "\n\nscript to print out the top 100 pfam domnains from a GFF3 file
if iprscan info has been pushed on to Dbxref in column 9\n\n
USAGE: pfam_domain_report.pl <GFF3 file>\n\n";

die($usage) unless $ARGV[0];

my $FILE1 = $ARGV[0];
my %GFF3;
my %DATA;
my %PFDS;
my @GENES;


parse_gff3($FILE1);
load_data();
pfd();
print "pfam domains\n";
pfds();
#print "\n\nThank you for choosing Mike\'s lamprey_script.pl\n\n";
print "Frac genes with \"True\" iprscan hits:". $DATA{t_w_iprs_h}/$DATA{t_trans}."\n";
#PostData(\%DATA);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
#prints the top 100 pfam ids for file
sub pfds{
    my $number = 0;
    foreach my $dom (sort {$PFDS{$b} <=> $PFDS{$a}} keys %PFDS){
        print $dom."\t".$PFDS{$dom}."\n" unless $number >= 100;
        $number++;
    }
 }
#-----------------------------------------------------------------------------
sub parse_gff3{
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
        chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	next if $line =~ /^\s*$/;    
	next if $line =~ /^\>/;    
	my @GFF_col = split(/\t/, $line);
	
	if ($GFF_col[2] eq 'mRNA'){
	    my @mrna9 = split(/\;/, $GFF_col[8]); 
	    $mrna9[0] =~ s/ID=//;
	    $mrna9[0] =~ s/\-RA//;
	    foreach my $pair (@mrna9){
		my ($k, $v) = split(/\=/,$pair);
		$GFF3{$mrna9[0]}{$k}=$v;
	    }           
	}
	
    }
    close($FH);
}
######################### example data structure ########################
#{PMZ_0024127-RA} = HASH
#{Alias} = maker-scaffold_23256.1-1073-augustus-gene-0.0-mRNA-1
#{Dbxref} = InterPro:IPR001593,InterPro:IPR018281,PANTHER:PTHR11830,Pfam:PF01015,Prosite:PS01191
#{Name} = PMZ_0024127-RA
#{Note} = Similar to rps3a: 40S ribosomal protein S3a (Danio rerio)
#{Ontology_term} = GO:0003735,GO:0005622,GO:0005840,GO:0006412
#{PMZ_0024127-RA} = UNDEFINED VALUE
#{Parent} = PMZ_0024127
#{_AED} = 0.00
#{_QI} = 0|-1|0|1|-1|1|1|0|274

#-----------------------------------------------------------------------------
sub pfd{
    foreach my $x (@GENES){
        if (defined($GFF3{$x}{Dbxref})){
            my @stuff = split(/\,/,$GFF3{$x}{Dbxref});
            foreach my $y (@stuff){
                if ($y =~ /^Pfam:/){
                    my ($k, $v) = split(/\:/,$y);
                    $PFDS{$v}++;
                }
            }
        }
    }
}
#-----------------------------------------------------------------------------
########################
sub load_data{
    my @keys= keys %GFF3;
    foreach my $key (@keys){
	$DATA{t_trans}++;
	if ($GFF3{$key}{Dbxref}){
	    push(@GENES, $key);
	    $DATA{t_w_iprs_h}++;
	}
    }	
}
#-----------------------------------------------------------------------------

