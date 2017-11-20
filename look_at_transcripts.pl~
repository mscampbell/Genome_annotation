#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;
use gff3_annotation_stuff_b;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
my %DATA;
my ($gff_hr, $lu_hr) = gff3_annotation_stuff_b::parse_gff($FILE);
my %GID_UTR;
build_data($gff_hr);
count_transcripts_and_cdss();
UTR_presence_check();
#get the median CDS length
#get the median cds lenth
#get the median 5' and 3' utr lengths
#get the number of gene s with 3' or 5' utr 
#PostData($gff_hr);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub UTR_presence_check{
    #PostData(\%GID_UTR);
    my $genes_with_tputr=0;
    my $genes_with_fputr=0;
    my $trans_with_tputr=0;
    my $trans_with_fputr=0;
    
    
    foreach my $gid (keys %{$GID_UTR{'g'}}){
	$genes_with_tputr++ if defined($GID_UTR{'g'}{$gid}{'tp'}) && $GID_UTR{'g'}{$gid}{'tp'} > 0 ;
	$genes_with_fputr++ if defined($GID_UTR{'g'}{$gid}{'fp'}) && $GID_UTR{'g'}{$gid}{'fp'} > 0;
    }
    foreach my $tid (keys %{$GID_UTR{'t'}}){
	$trans_with_tputr++ if defined($GID_UTR{'t'}{$tid}{'tp'}) && $GID_UTR{'t'}{$tid}{'tp'} > 0;
	$trans_with_fputr++ if defined($GID_UTR{'t'}{$tid}{'fp'}) && $GID_UTR{'t'}{$tid}{'fp'} > 0;
    }
    print "number of genes with five prime UTR\t$genes_with_fputr\n";
    print "number of genes with three prime UTR\t$genes_with_tputr\n";
    print "number of trans with five prime UTR\t$trans_with_fputr\n";
    print "number of trans with three prime UTR\t$trans_with_tputr\n";

}
#-----------------------------------------------------------------------------
sub count_transcripts_and_cdss{
    #this takes all of the flattened CDSs aned exons and drops them into a 
    #hash so the redundent ones get collapsed
    my $t_t  = 0;
    my $t_ut = 0;
    my $t_uc = 0;

    my %ct;
    foreach my $rs (keys %DATA){
	foreach my $tid (keys %{$DATA{$rs}}){
	    $t_t++;
	    foreach my $ut (keys %{$DATA{$rs}{$tid}{'exons'}}){
		$ct{$rs}{'exons'}{$ut}++;
	    }
	    foreach my $uc (keys %{$DATA{$rs}{$tid}{'cds'}}){
		$ct{$rs}{'cds'}{$uc}++;
	    }
	}
    }

    #this loop counts the numbner of unique CDS and unique transcripts
    foreach my $rs (keys %ct){
	foreach my $ex (keys  %{$ct{$rs}{'exons'}}){
	    $t_ut++; #it is importantt that these counters don't
	}            #add the values of the keys
	foreach my $cs (keys  %{$ct{$rs}{'cds'}}){
	    $t_uc++;
	}
	
    }
    print "total transcripts\t $t_t\n";
    print "total unique transcripts\t $t_ut\n";
    print "total unique cds\t $t_uc\n\n";
}
#-----------------------------------------------------------------------------
sub flatten{
    # This takes an array of arrays and reduces it to a scalar that
    # can be used in string comparisons
    my $mda = shift;
    my @f_array;
    foreach my $ar (@$mda){
	push(@f_array , $ar->[0]);
	push(@f_array , $ar->[1]);
    }
    my $fas = join("_", @f_array);
    return $fas;
}
#-----------------------------------------------------------------------------
sub median{
    my $array = shift;
    @{$array} = sort {$a <=> $b} @{$array};
    my $count = scalar (@{$array});
    if ($count % 2) {
        return $array->[int($count/2)];
    } 
    else {
        return ($array->[$count/2] + $array->[$count/2 - 1]) / 2;
    }
}

#-----------------------------------------------------------------------------
sub build_data{
    # this takes the hash reference from the gff3_annotation_stuff::parse_gff
    # call in main and puts it into another structure that is a littel easier
    # to work with
    my $hr = shift;
    #PostData($hr);
    my @cds_lengths;
    my @trans_lengths;
    my @tpu_lengths;
    my @fpu_lengths;
    my $t_g = 0;
   foreach my $gid (keys %{$hr}){
       $t_g++;
	my $rs = $hr->{$gid}->{'gl'}->{'refseq'};
	print $gid, "mmm\n" if !defined($rs);
	foreach my $tid (keys %{$hr->{$gid}->{'tl'}}){
	    
	    #build the structure to help with the string comparisons
	    my $f_cds = flatten($hr->{$gid}->{'tl'}->{$tid}->{'CDS'});
	    my $f_exons = flatten($hr->{$gid}->{'tl'}->{$tid}->{'exons'});
	    $DATA{$rs}{$tid}{'cds'}{$f_cds}++;
	    $DATA{$rs}{$tid}{'exons'}{$f_exons}++;
	    
	    #while you are hear loop through and build an array of cds lengths
	    my $cds_len = 0;
	    foreach my $cds_chunks (@{$hr->{$gid}->{'tl'}->{$tid}->{'CDS'}}){
		my ($beg, $end) = @{$cds_chunks};
		$cds_len += $end - $beg +1;
	    }
	    #and do the same thing for the spliced transcript
	    my $trans_len = 0;
	    foreach my $exon_chunks (@{$hr->{$gid}->{'tl'}->{$tid}->{'exons'}}){
		my ($beg, $end) = @{$exon_chunks};
		$trans_len += $end - $beg +1;
	    }
	    #and for the UTRS
	    my $tpu_len = 0;
	    foreach my $exon_chunks (@{$hr->{$gid}->{'tl'}->{$tid}->{'three_prime_UTR'}->{'coordinates'}}){
		my ($beg, $end) = @{$exon_chunks};
		#print "$tid\tbeg:$beg\tend:$end\n";
		$tpu_len += $end - $beg +1;
		$GID_UTR{'g'}{$gid}{'tp'} += $tpu_len;
		$GID_UTR{'t'}{$tid}{'tp'} += $tpu_len;
	    }
	    my $fpu_len = 0;
	    foreach my $exon_chunks (@{$hr->{$gid}->{'tl'}->{$tid}->{'five_prime_UTR'}->{'coordinates'}}){
		my ($beg, $end) = @{$exon_chunks};
		#print "$tid\tbeg:$beg\tend:$end\n";
		$fpu_len += $end - $beg +1;
		$GID_UTR{'g'}{$gid}{'fp'} += $fpu_len;
		$GID_UTR{'t'}{$tid}{'fp'} += $fpu_len;
	    }
	    push(@fpu_lengths, $fpu_len)     unless $fpu_len   == 0;
	    push(@tpu_lengths, $tpu_len)     unless $tpu_len   == 0;
	    push(@cds_lengths, $cds_len)     unless $cds_len   == 0;
	    push(@trans_lengths, $trans_len) unless $trans_len == 0;
	    #PostData(\@fpu_lengths);
	}
   }
    my $median_fpu_len = median(\@fpu_lengths);
    my $median_tpu_len = median(\@tpu_lengths);
    my $median_cds_len = median(\@cds_lengths);
    my $median_trans_len = median(\@trans_lengths);
    print "median cds length\t$median_cds_len\n";
    print "median transcript length\t$median_trans_len\n";
    print "median five prime UTR length\t$median_fpu_len\n";
    print "median three prime UTR length\t$median_tpu_len\n";
    print "total gene count\t$t_g\n";
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
#	last if $line =~ /\#\#FASTA/;
#	next if $line =~ /\#/;
#	
#	my @cols = split(/\t/, $line);
#	if $cols[
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

