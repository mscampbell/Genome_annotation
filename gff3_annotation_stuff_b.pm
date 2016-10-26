#!/usr/bin/perl -w 
#------------------------------------------------------------------------
#----                  gff3_annotation_stuff.pm                      ---- 
#------------------------------------------------------------------------
package gff3_annotation_stuff_b;
use strict;
#use PostData;
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#----------------------------------FUNCTIONS----------------------------------
#-----------------------------------------------------------------------------
sub parse_gff{
    #this function returns a hash ref with the genes in it and a hash ref with
    #gid tid look up information 
    #$gid_tid_lu{'g2t'}{$gid} = $tid; 
    #$gid_tid_lu{'t2g'}{$tid} = $gid;
    #$gid_tid_lu{'tn2g'}{$tname} = $gid; helpfull if you got the ids from a 
    #FASTA file

    #print STDERR "Now Parsing GFF\n";
    my $file = shift;
    my %gff;
    my %gid_tid_lu;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
        chomp($line);
        last if $line =~ /^\#\#FASTA/;
        next if $line =~ /^\#/;
        my @cols = split(/\t/ , $line);

#parse the gene lines
        if ($cols[2] eq 'gene'){
            my $gid;
            if ($cols[8] =~ /ID=(\S+?);/ || $cols[8] =~ /ID=(\S+)/){
                $gid = $1;
            }
            else {die "Not valid gff3";}
            #load everything before column 9
            $gff{$gid}{'gl'}{'begin'}   = $cols[3];
            $gff{$gid}{'gl'}{'end'}     = $cols[4];
            $gff{$gid}{'gl'}{'strand'}  = $cols[6];
            $gff{$gid}{'gl'}{'refseq'}  = $cols[0];
            $gff{$gid}{'gl'}{'source'}  = $cols[1];
            $gff{$gid}{'gl'}{'score'}   = $cols[5];
            $gff{$gid}{'gl'}{'phase'}   = $cols[7];
            
            #load column 9
            my @col9 = split(/;/, $cols[8]);
            foreach my $kvp (@col9){
                my ($k, $v) = split(/=/, $kvp);
                $gff{$gid}{'gl'}{$k} = $v;
            }
        }
#parse the mRNA lines
        elsif ($cols[2] eq 'mRNA'){
            my $tid;
            my $gid;
	    my $tname;
            if ($cols[8] =~ /ID=(\S+?);/){
                $tid = $1;
	    }
            else {die "Not valid gff3";}
            
            if ($cols[8] =~ /Parent=(\S+?);/ || $cols[8] =~ /Parent=(\S+)/){
                $gid = $1;
            }
            else {die "Not valid gff3";}
	    if ($cols[8] =~ /Name=(\S+?);/ || $cols[8] =~ /Name=(\S+)/){
		$tname = $1
	    }
            #build a look up for later
            $gid_tid_lu{'g2t'}{$gid} = $tid;
            $gid_tid_lu{'t2g'}{$tid} = $gid;
	    $gid_tid_lu{'tn2g'}{$tname} = $gid;
	    
            #load everything before column 9
            $gff{$gid}{'tl'}{$tid}{'begin'}   = $cols[3];
            $gff{$gid}{'tl'}{$tid}{'end'}     = $cols[4];
            $gff{$gid}{'tl'}{$tid}{'strand'}  = $cols[6];
            $gff{$gid}{'tl'}{$tid}{'refseq'}  = $cols[0];
            $gff{$gid}{'tl'}{$tid}{'source'}  = $cols[1];
            $gff{$gid}{'tl'}{$tid}{'score'}   = $cols[5];
            #$gff{$gid}{'tl'}{$tid}{'phase'}   = $cols[7];
            
            #load column 9
            my @col9 = split(/;/, $cols[8]);
            foreach my $kvp (@col9){
                my ($k, $v) = split(/=/, $kvp);
                
                if ($k eq '_QI'){
                    my @qi = split(/\|/, $v);

                    $gff{$gid}{'tl'}{$tid}{'five_prime_utr_length'} = $qi[0]; 
                    $gff{$gid}{'tl'}{$tid}{'frac_ss_confirmed_by_est'} = $qi[1]; 
                    $gff{$gid}{'tl'}{$tid}{'frac_exons_matching_est'} = $qi[2]; 
                    $gff{$gid}{'tl'}{$tid}{'frac_exons_with_est_or_pro_overlap'} = $qi[3]; 
                    $gff{$gid}{'tl'}{$tid}{'frac_ss_confirmed_by_gene_pred'} = $qi[4]; 
                    $gff{$gid}{'tl'}{$tid}{'frac_exons_with_gene_pred_overlap'} = $qi[5]; 
                    $gff{$gid}{'tl'}{$tid}{'number_of_exons'} = $qi[6]; 
                    $gff{$gid}{'tl'}{$tid}{'three_prime_utr_length'} = $qi[7]; 
		    $gff{$gid}{'tl'}{$tid}{'lenght_of_protein_product'} = $qi[8];
                }
                else{
                    $gff{$gid}{'tl'}{$tid}{$k} = $v;
                }
            }
            
        }
#parse the exon, cds, and utr lines
        elsif ($cols[2] eq 'exon' ||
               $cols[2] eq 'CDS' ||
               $cols[2] eq 'three_prime_UTR' ||
               $cols[2] eq 'five_prime_UTR'){
               
            my @tids; #becasue of multiple parenting this my have more
                     #than one id in it (comma separated). If not this
                     #will jsut turn into an array with one value
            if ($cols[8] =~ /Parent=(\S+?);/ || $cols[8] =~ /Parent=(\S+)/){
                @tids = split(/,/, $1);

            }
            else {die "Not valid gff3";}
            #get the featrue ids of things that are not exons
            my $fid;
            #if (($cols[2] eq 'CDS' || #this works for MAKER output but the ensemble dump doesn't put IDs on CDSs
	    if (($cols[2] eq 'three_prime_UTR' ||
                 $cols[2] eq 'five_prime_UTR') && 
		($cols[8] =~ /ID=(\S+?);/      || 
		 $cols[8] =~ /ID=(\S+)/)){
                $fid = $1;
            }
	    elsif (($cols[2] eq 'three_prime_UTR') && #another edit for enesembl
                ($cols[8] =~ /Parent=(\S+?);/      || $cols[8] =~ /Parent=(\S+)/)){
		$fid = $1;
		$fid .= "tputr";
	   }
	    elsif (($cols[2] eq 'five_prime_UTR') && #another edit for ensembl
                ($cols[8] =~ /Parent=(\S+?);/      || $cols[8] =~ /Parent=(\S+)/)){
		$fid = $1;
		$fid .= "fputr";
	   }
            
            else {die " I don't think this is valid gff3. I need ids for cds and utr features\n$cols[2] $cols[8]\n" unless $cols[2] eq 'exon' || $cols[2] eq 'CDS';}
            
            foreach my $tid (@tids){
                my $gid = $gid_tid_lu{'t2g'}{$tid};
                if ($cols[2] eq 'exon'){
                    push (@{$gff{$gid}{'tl'}{$tid}{'exons'}}, [$cols[3], $cols[4]]);
                }
                elsif ($cols[2] eq 'CDS'){
                    push (@{$gff{$gid}{'tl'}{$tid}{'CDS'}}, [$cols[3], $cols[4]]);
                }
                else {
                    $gff{$gid}{'tl'}{$tid}{$cols[2]}{'id'} = $fid;
                    if ($cols[2] eq 'CDS'){
                        push (@{$gff{$gid}{'tl'}{$tid}{$cols[2]}{'coordinates_and_pahse'}}, [$cols[3], $cols[4], $cols[7]]);
                    }
                    else{
                        push (@{$gff{$gid}{'tl'}{$tid}{$cols[2]}{'coordinates'}}, [$cols[3], $cols[4]]);
                    }
                }
            }
        }
        
        else {
            #print STDERR "WARNING unexpected typ: $cols[2]\n";
        }
    }
    $fh->close();
    return (\%gff, \%gid_tid_lu);
}


#-----------------------------------------------------------------------------
1;
