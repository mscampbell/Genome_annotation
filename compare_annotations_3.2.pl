#!/usr/bin/perl -w 
$| = 1;
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use Getopt::Std;
use vars qw($opt_t $opt_c $opt_u $opt_s $opt_o $opt_v $opt_i);
getopts('tc:u:sovi:');
#use PostData;
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\tMIKE_compare_annotations_3.x.pl: Compares 2 GFF3 files and returns 
\tcomparrison statistics OR takes one GFF3 file and returns basic statistics 
\tinput is one or two GFF3 files.\n\n

VERSION 3 NOTE: I added sensitivity, specificity, and accuracy calculations 
                along with median exons per mRNA.  I sped it up a little and 
                removed the necessity to modify the parse subroutine to match 
                your specific contig, scaffold, or chromosome. I also added 
                the option to run the script on a specific scaffold.
             
USAGE:          compare_annotations_3.pl <GFF3_file_a> <GFF3 file_b>
                compare_annotations_3.pl <GFF3_file_a>

OPTIONS:        -t (trim) standardize scaffold names. 
                Note: you may need to alter the sub trim.
                -c <scaffold name> limits calculations to a scaffold of interest.
                -u <a or b> prints out the genes unique to the specified file
                -i <a or b> prints out the genes not unique to the specified file
                \n\n";

die($usage) unless $ARGV[0];

my $A = $ARGV[0];
my $B = $ARGV[1];

my %LU;
my %EX;

my %G;
my %I;
my %R;

my %S;
my %MU;
my %MI;

my %EXLEN;
my %LEN_U;
my %LEN_I;

my %COL;
my @FOO;

parse($A,'a');
parse($B,'b');

build_G('a');
build_G('b');

intersection();

report();
stacked_column(@FOO) if $opt_s; #opt_o modifies the print out here
one_to_one() if $opt_v;
get_overlap($opt_i, 'mRNA') if $ARGV[1] && $opt_i;
#PostData(\%LU);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub one_to_one{ #this was a hack to get at the one to one overlapping genes. This gives ids for all of the overlapping genes I wrote and accessory script to get the one to ones called get_one2ones.pl. this has only been tested on arabidopsis data
    my $out_a = 'intersect_a';
    my $out_b = 'intersect_b';
    my $fh_a = new FileHandle;
    $fh_a->open(">$out_a");
    
    foreach my $scaff (keys %I) {
	foreach my $a_gid (keys %{$I{$scaff}{INT}{'a'}}){
	    print $fh_a "$a_gid\t" unless $I{$scaff}{INT}{'a'}{$a_gid}[1];
	    print $fh_a join("\t", @{$I{$scaff}{INT}{'a'}{$a_gid}}) unless $I{$scaff}{INT}{'a'}{$a_gid}[1];
	    print $fh_a "\n" unless $I{$scaff}{INT}{'a'}{$a_gid}[1];
	}
    }
    $fh_a->close();
}
#-----------------------------------------------------------------------------
sub stacked_column{
    my $a = shift;
    my $b = shift;
    my $sh = shift;
    my $sn = shift;
    my $sp = shift;
    my $ac = shift;
    
    print "\nstacked_column\n\n";
    if ($opt_o){
	print "$a\n";
	print "$b\n";
	print "$sh\n";
	print "$ARGV[0]\n";
	print "$ARGV[1]\n";
	print "$sn\n";
	print "$sp\n";
	print "$ac\n";
    }

    else {
	print "a_only\t$a\n";
	print "b_only\t$b\n";
	print "shared\t$sh\n";
	print "file_a:\t$ARGV[0]\n";
	print "file_b:\t$ARGV[1]\n";
	print "sn:\t$sn\n";
	print "sp:\t$sp\n";
	print "ac:\t$ac\n";
    }
}
#-----------------------------------------------------------------------------
sub report {


        my $both_scaff   = $R{scaff_report}{both};
        my $a_only_scaff = $R{scaff_report}{a_only} || 0;
        my $b_only_scaff = $R{scaff_report}{b_only} || 0 if $ARGV[1];
        my $tot_scaff_2  = $both_scaff + $a_only_scaff +  $b_only_scaff if $ARGV[1];
        my $tot_scaff_1  = $a_only_scaff  unless $ARGV[1];
	my $tot_bases_a  = $MU{'a'}{t_len};
	my $tot_bases_b  = $MU{'b'}{t_len};
	print "*********************************************\n";
	print "********** SCAFFOLD SUMMARY *****************\n";
	print "*********************************************\n";
	print "TOTAL BASES ANALYZED in A:$tot_bases_a\n" if $MU{a}{t_len};
	print "TOTAL BASES ANALYZED in A:NA\n" unless $MU{a}{t_len};
	print "TOTAL BASES ANALYZED in B:$tot_bases_b\n" if $ARGV[1] && $MU{b}{t_len};
	if ($ARGV[1]) {
	    print "TOTAL BASES ANALYZED in B:NA\n" unless $MU{b}{t_len};
	}
	print "TOTAL SCAFFOLDS ANALYZED:$tot_scaff_2\n" if $ARGV[1];
	print "TOTAL SCAFFOLDS ANALYZED:$tot_scaff_1\n" unless $ARGV[1];
	print "NUM. UNIQUE TO A:$a_only_scaff \n" if $ARGV[1];
	print "NUM. UNIQUE TO B:$b_only_scaff \n" if $ARGV[1];
	
	my $tot_genes_a_pc = get_tot_genes('a', 'mRNA');
	my $tot_genes_b_pc = get_tot_genes('b', 'mRNA') if $ARGV[1];
        my $tot_genes_a_nc = get_tot_genes('a', 'ncRNA');
        my $tot_genes_b_nc = get_tot_genes('b', 'ncRNA') if $ARGV[1];

        print "*********************************************\n";
        print "************** GENE SUMMARY *****************\n";
        print "*********************************************\n";
	print "NUM. PC genes in A:$tot_genes_a_pc\n";
	print "NUM. PC genes in B:$tot_genes_b_pc\n" if $ARGV[1]; 	
        print "NUM. NC genes in A:$tot_genes_a_nc\n";    
        print "NUM. NC genes in B:$tot_genes_b_nc\n" if $ARGV[1];
	print "Gene density in  A:" . 
	    ($tot_genes_a_pc + $tot_genes_a_nc)/($tot_bases_a/1000) 
	    ." per kb\n" if $tot_bases_a;
	print "Gene density in  B:" . 
	    ($tot_genes_b_pc + $tot_genes_b_nc)/($tot_bases_b/1000) 
	    ." per kb\n" if $tot_bases_b;


        my $num_genes_in_a_w_o_in_b = get_num_w_overlap('a') if $ARGV[1];
        my $num_genes_in_b_w_o_in_a = get_num_w_overlap('b') if $ARGV[1];
	my $or = $num_genes_in_a_w_o_in_b/$num_genes_in_b_w_o_in_a if $ARGV[1]; 
	my $a_uniq_pc = get_unique_count('a', 'mRNA') if $ARGV[1];
	my $a_uniq_nc = get_unique_count('a', 'ncRNA') if $ARGV[1];
        my $b_uniq_pc = get_unique_count('b', 'mRNA') if $ARGV[1];
        my $b_uniq_nc = get_unique_count('b', 'ncRNA') if $ARGV[1];

	my $tp_a2b = $num_genes_in_b_w_o_in_a if $ARGV[1];
	my $tp_plus_fn_a2b = $tp_a2b + $b_uniq_pc if $ARGV[1];
	my $tp_plus_fp_a2b = $tp_a2b + $a_uniq_pc if $ARGV[1];
	my $sn_a2b = $tp_a2b/$tp_plus_fn_a2b if $ARGV[1];
	my $sp_a2b = $tp_a2b/$tp_plus_fp_a2b if $ARGV[1];
	my $ac_a2b = ($sn_a2b + $sp_a2b)/2 if $ARGV[1];

	my $tp_b2a = $num_genes_in_a_w_o_in_b if $ARGV[1];
	my $tp_plus_fn_b2a = $tp_b2a + $a_uniq_pc if $ARGV[1];
	my $tp_plus_fp_b2a = $tp_b2a + $b_uniq_pc if $ARGV[1];
	my $sn_b2a = $tp_b2a/$tp_plus_fn_b2a if $ARGV[1];
	my $sp_b2a = $tp_b2a/$tp_plus_fp_b2a if $ARGV[1];
	my $ac_b2a = ($sn_b2a + $sp_b2a)/2 if $ARGV[1];

        print "*********************************************\n" if $ARGV[1];
        print "********** GENE OVERLAP SUMMARY *************\n" if $ARGV[1];
        print "*********************************************\n" if $ARGV[1];
	print "NUM. genes in A W/one or more overlaping gene(s) in B:$num_genes_in_a_w_o_in_b\n" if $ARGV[1];
	print "NUM. genes in B W/one or more overlaping gene(s) in A:$num_genes_in_b_w_o_in_a\n" if $ARGV[1];
	print "Within this intersection, the A:B overlap ratio is:$or \n" if $ARGV[1];	
	print "NUM. PC genes unique to A:$a_uniq_pc\n" if $ARGV[1];
	print "NUM. NC genes unique to A:$a_uniq_nc\n" if $ARGV[1];
        print "NUM. PC genes unique to B:$b_uniq_pc\n" if $ARGV[1];
        print "NUM. NC genes unique to B:$b_uniq_nc\n" if $ARGV[1];
	print "SN. a2b:$sn_a2b\n" if $ARGV[1];
	print "SP. a2b:$sp_a2b\n" if $ARGV[1];
	print "AC. a2b:$ac_a2b\n" if $ARGV[1];
	print "SN. b2a:$sn_b2a\n" if $ARGV[1];
	print "SP. b2a:$sp_b2a\n" if $ARGV[1];
	print "AC. b2a:$ac_b2a\n" if $ARGV[1];
	@FOO = ($a_uniq_pc, $b_uniq_pc, $num_genes_in_b_w_o_in_a, $sn_a2b, $sp_a2b, $ac_a2b);
ex_path_u();	

	my $a_ave_g_len   = $MU{a}{t_len_gene}/($tot_genes_a_pc + $tot_genes_a_nc);
	my $a_ave_ex_len  = $MU{a}{t_len_ex}/$MU{a}{num_ex};
	my $a_ave_in_len  = $MU{a}{t_len_in}/$MU{a}{num_in};
	my $a_ave_ex_mRNA = $MU{a}{num_ex}/$MU{a}{num_mRNA};
	my $a_med_ex_mRNA = median_exons('a', 'u');
	my $b_ave_g_len   = $MU{b}{t_len_gene}/($tot_genes_b_pc + $tot_genes_b_nc) if $ARGV[1];
	my $b_ave_ex_len  = $MU{b}{t_len_ex}/$MU{a}{num_ex} if $ARGV[1];
	my $b_ave_in_len  = $MU{b}{t_len_in}/$MU{b}{num_in} if $ARGV[1];
	my $b_ave_ex_mRNA = $MU{b}{num_ex}/$MU{b}{num_mRNA} if $ARGV[1];
	my $b_med_ex_mRNA = median_exons('b', 'u') if $ARGV[1];
	my $a_med_exlen_u = median($LEN_U{'exons'}{'a'});
	my $b_med_exlen_u = median($LEN_U{'exons'}{'b'}) if $ARGV[1];
	my $a_med_in_len  = median($LEN_U{'ints'}{'a'});
	my $b_med_in_len  = median($LEN_U{'ints'}{'b'}) if $ARGV[1];
	my $a_med_g_len   = median($LEN_U{'trans'}{'a'});
	my $b_med_g_len   = median($LEN_U{'trans'}{'b'}) if $ARGV[1];
	my $a_max_in_len  = max($LEN_U{'ints'}{'a'});
	my $b_max_in_len  = max($LEN_U{'ints'}{'b'}) if $ARGV[1];

	print "*********************************************\n";
        print "*********** EXON/INTRON SUMMARY *************\n";
        print "*********************************************\n";
	print "UNION\n" if $ARGV[1];
	print "AVE. gene length    A:$a_ave_g_len\n";
	print "MED. gene length    A:$a_med_g_len\n";
	print "AVE. gene length    B:$b_ave_g_len\n" if $ARGV[1];	
	print "MED. gene length    B:$b_med_g_len\n" if $ARGV[1];
	print "AVE. exons per mRNA A:$a_ave_ex_mRNA\n";
	print "MED. exons per mRNA A:$a_med_ex_mRNA\n";
	print "AVE. exons per mRNA B:$b_ave_ex_mRNA\n" if $ARGV[1];	
	print "MED. exons per mRNA B:$b_med_ex_mRNA\n" if $ARGV[1];	
	print "AVE. exon length    A:$a_ave_ex_len\n";
	print "MED. exon length    A:$a_med_exlen_u\n";
	print "AVE. exon length    B:$b_ave_ex_len\n" if $ARGV[1];
        print "MED. exon length    B:$b_med_exlen_u\n" if $ARGV[1];
	print "AVE. intron length  A:$a_ave_in_len\n";
	print "MED. intron length  A:$a_med_in_len\n";
	#print "MAX. intron length  A:$a_max_in_len\n";
        print "AVE. intron length  B:$b_ave_in_len\n" if $ARGV[1];        
        print "MED. intron length  B:$b_med_in_len\n" if $ARGV[1];        
        #print "MAX. intron length  B:$b_max_in_len\n" if $ARGV[1];        

#	PostData(\%EXLEN);

ex_path_int() if $ARGV[1];
i_gene_path() if $ARGV[1];
	
	if ($ARGV[1]){
	my $ai_ave_g_len   = $MI{a}{t_len_gene}/($num_genes_in_a_w_o_in_b);
        my $ai_ave_ex_len  = $MI{a}{t_len_ex}/$MI{a}{num_ex};
        my $ai_ave_in_len  = $MI{a}{t_len_in}/$MI{a}{num_in};
        my $ai_ave_ex_mRNA = $MI{a}{num_ex}/$MI{a}{num_mRNA};
        my $bi_ave_g_len   = $MI{b}{t_len_gene}/($num_genes_in_b_w_o_in_a);
        my $bi_ave_ex_len  = $MI{b}{t_len_ex}/$MI{a}{num_ex};
        my $bi_ave_in_len  = $MI{b}{t_len_in}/$MI{b}{num_in};
        my $bi_ave_ex_mRNA = $MI{b}{num_ex}/$MI{b}{num_mRNA};
	my $ai_med_exlen   = median($LEN_I{'exons'}{'a'});
	my $bi_med_exlen   = median($LEN_I{'exons'}{'b'});
	my $ai_med_in_len  = median($LEN_I{'ints'}{'a'});
	my $bi_med_in_len  = median($LEN_I{'ints'}{'b'});
	my $ai_med_g_len   = median($LEN_I{'trans'}{'a'});
        my $bi_med_g_len   = median($LEN_I{'trans'}{'b'}) if $ARGV[1];
	my $ai_med_ex_mRNA = median_exons('a', 'i');
	my $bi_med_ex_mRNA = median_exons('b', 'i');

	print "INTERSECTION\n";
	print "AVE. gene length    A:$ai_ave_g_len\n";
        print "MED. gene length    A:$ai_med_g_len\n";
	print "AVE. gene length    B:$bi_ave_g_len\n";
	print "MED. gene length    B:$bi_med_g_len\n";
        print "AVE. exons per mRNA A:$ai_ave_ex_mRNA\n";
	print "MED. exons per mRNA A:$ai_med_ex_mRNA\n";
        print "AVE. exons per mRNA B:$bi_ave_ex_mRNA\n";
	print "MED. exons per mRNA B:$bi_med_ex_mRNA\n";
        print "AVE. exon length    A:$ai_ave_ex_len\n";
        print "MED. exon length    A:$ai_med_exlen\n";
	print "AVE. exon length    B:$bi_ave_ex_len\n";
        print "MED. exon length    B:$bi_med_exlen\n";
	print "AVE. intron length  A:$ai_ave_in_len\n";
	print "MED. intron length  A:$ai_med_in_len\n";
        print "AVE. intron length  B:$bi_ave_in_len\n";
	print "MED. intron length  B:$bi_med_in_len\n";
    }
}
#-----------------------------------------------------------------------------
sub get_overlap {
    my $flag = shift;
    my $type = shift;
    
    foreach my $scaff (keys %I){
	foreach my $gid (keys %{$I{$scaff}{INT}{$flag}}){
	    print "Overlapping FILE:$flag\tGID:$gid\tSCAF:$scaff\n";
	}
    }
}
#-----------------------------------------------------------------------------
sub get_unique_count {
	my $flag = shift;
	my $type = shift;
	my $tot = 0;
	
	foreach my $scaff (keys %I){
		foreach my $gid (@{$I{$scaff}{UNQ}{$flag}}){
			$tot++ if $G{$flag}{$scaff}{$gid}{t} eq $type;
                        
                        #print out the uniqe gene information if -u was used
			if ($opt_u && ($flag eq $opt_u) 
                                   && 
                            ($G{$opt_u}{$scaff}{$gid}{t} eq $type)){
			    print "Unique FILE:$flag\tGID:$gid\tSCAF:$scaff";
			    print "\tbigin:$G{$flag}{$scaff}{$gid}{'b'}";
			    print "\tend:$G{$flag}{$scaff}{$gid}{'e'}\n";
			}
		}
	}
	return $tot;
}
#-----------------------------------------------------------------------------
sub get_num_w_overlap {
	my $flag = shift;

	my $tot = 0;
	foreach my $scaff (keys %I){
		$tot += keys %{$I{$scaff}{INT}{$flag}};	
	}

	return $tot;
}
#-----------------------------------------------------------------------------
sub get_tot_genes {
	my $flag = shift;
	my $type = shift;


	my $tot = 0;
	foreach my $scaff (keys %{$G{$flag}}){
		foreach my $g_id (keys %{$G{$flag}{$scaff}}){
		    my $crap = $G{$flag}{$scaff}{$g_id}{t};
			print "CRAP:$crap\n";
		    print "$g_id\n";
		    $tot++ if $G{$flag}{$scaff}{$g_id}{t} eq $type;
		}
	}

	return $tot;
}
#-----------------------------------------------------------------------------
sub intersection {


	while (defined (my $scaff = each %S)){
		if    (exists($S{$scaff}{'a'}) && exists($S{$scaff}{'b'})){
			$R{scaff_report}{both}++;
		}
		elsif (exists($S{$scaff}{'a'})){
			$R{scaff_report}{a_only}++;
		}
		elsif (exists($S{$scaff}{'b'})){
			$R{scaff_report}{b_only}++;
		}

	
		my @a_gids = keys %{$G{'a'}{$scaff}};
		my @b_gids = keys %{$G{'b'}{$scaff}};
		

		list_logic(\@a_gids, \@b_gids, $scaff);
		t_len(\@a_gids, \@b_gids, $scaff, 'u');
	}
}
#-----------------------------------------------------------------------------
#MIKE#                                                                     
sub i_gene_path{    
    foreach my $scaff (keys %I){
	get_len($scaff);
    }
    #t_len(\@a_gids, \@b_gids, $scaff, 'i');
}
#-----------------------------------------------------------------------------
#MIKE#
sub get_len{
    my $scaff = shift;
    my @a_gids = keys %{$I{$scaff}{INT}{'a'}};
    my @b_gids = keys %{$I{$scaff}{INT}{'b'}};
    t_len(\@a_gids, \@b_gids, $scaff, 'i');
}
#-----------------------------------------------------------------------------
sub t_len{
    my $a_gids = shift;
    my $b_gids = shift;
    my $scaff  = shift;
    my $x      = shift;
    foreach my $a_gid (@{$a_gids}){
	my $a_gene = $G{'a'}{$scaff}{$a_gid};
	my $a_b = $a_gene->{b};
	my $a_e = $a_gene->{e};
	my $l_t = $a_gene->{e} - $a_gene->{b} +1;
	if ($x eq 'u'){	
	    $MU{a}{t_len_gene} += $l_t;
	    push(@{$LEN_U{'trans'}{'a'}}, $l_t);
	}  
	if ($x eq 'i'){
            $MI{a}{t_len_gene} += $l_t;
	    push(@{$LEN_I{'trans'}{'a'}}, $l_t);
	    $MI{a}{t_genes}++;
        }    
    }
    foreach my $b_gid (@{$b_gids}){
        my $b_gene = $G{'b'}{$scaff}{$b_gid};
	my $b_b = $b_gene->{b};
        my $b_e = $b_gene->{e};
        my $l_t = $b_gene->{e} - $b_gene->{b} +1;
	if ($x eq 'u'){     
	    $MU{b}{t_len_gene} += $l_t;
	    push(@{$LEN_U{'trans'}{'b'}}, $l_t);
	}
	if ($x eq 'i'){
            $MI{b}{t_len_gene} += $l_t;
	    push(@{$LEN_I{'trans'}{'b'}}, $l_t);
	    $MI{b}{t_genes}++;
        }
    }
}
#-----------------------------------------------------------------------------
sub list_logic {
	my $a_gids = shift;
	my $b_gids = shift;
	my $scaff  = shift;

	
	foreach my $a_gid (@{$a_gids}){
		my $a_gene = $G{'a'}{$scaff}{$a_gid};
		my $a_b = $a_gene->{b};
		my $a_e = $a_gene->{e}; 
		my $a_s = $a_gene->{s};
		my $toggle = 0;
		foreach my $b_gid (@{$b_gids}){
			my $b_gene = $G{'b'}{$scaff}{$b_gid};
			my $b_b = $b_gene->{b};
			my $b_e = $b_gene->{e};
			my $b_s = $b_gene->{s};	

			next unless span_int($a_s, $a_b, $a_e, $b_s, $b_b, $b_e);
			#print "AAAAAAAAAAAAA\n";
			next unless e_int($a_gene, $b_gene);			
			#print "AAAAAAAAAAAAA\n";

			push(@{$I{$scaff}{INT}{'a'}{$a_gid}}, $b_gid);

			$toggle = 1;
		}

		push(@{$I{$scaff}{UNQ}{'a'}}, $a_gid) unless $toggle;
	}

        foreach my $b_gid (@{$b_gids}){
                my $b_gene = $G{'b'}{$scaff}{$b_gid};
                my $b_b = $b_gene->{b};
                my $b_e = $b_gene->{e};
                my $b_s = $b_gene->{s};
		my $toggle = 0;
                foreach my $a_gid (@{$a_gids}){
                        my $a_gene = $G{'a'}{$scaff}{$a_gid};
                        my $a_b = $a_gene->{b};
                        my $a_e = $a_gene->{e};
                        my $a_s = $a_gene->{s};
			#print "stuff:$b_s, $b_b, $b_e, $a_s, $a_b, $a_e\n";
                        #print "CCCCCCCCCCCC\n";
			next unless span_int($b_s, $b_b, $b_e, $a_s, $a_b, $a_e);
			#print "AAAAAAAAAAAAA\n";
                        next unless e_int($b_gene, $a_gene);
			#print "BBBBBBBBBBBBB\n";

                        push(@{$I{$scaff}{INT}{'b'}{$b_gid}}, $a_gid);

			$toggle = 1;
                }

		 push(@{$I{$scaff}{UNQ}{'b'}}, $b_gid) unless $toggle;
        }

}
#-----------------------------------------------------------------------------
sub e_int {

	my $a_gene = shift;
	my $b_gene = shift;

	my $a_s = $a_gene->{s};
	my $b_s = $b_gene->{s};

	return 0 unless $a_s eq $b_s;

	foreach my $a_exon (@{$a_gene->{x}}){
		my $a_b = $a_exon->[0];
		my $a_e = $a_exon->[1];

		foreach my $b_exon (@{$b_gene->{x}}){
			my $b_b = $b_exon->[0];
			my $b_e = $b_exon->[1];

			return 1 if span_int($a_s, $a_b, $a_e, $b_s, $b_b, $b_e);
			
		}
	}

	return 0;
}
#-----------------------------------------------------------------------------
sub span_int {
	my $a_s = shift;
	my $a_b = shift;
	my $a_e = shift;
	my $b_s = shift;
	my $b_b = shift;
	my $b_e = shift;
	
	#print "entered_span_int\n";
	return 0  unless $a_s eq $b_s;
	#print "passed_step_one $a_s:$b_s\n";
	#print "passed\n" if (($b_b <= $a_e && $b_b >= $a_b) 
	 #        ||  ($b_b <= $a_b && $b_e >= $a_e));

	return 1 if (($b_b <= $a_e && $b_b >= $a_b) 
	         ||  ($b_b <= $a_b && $b_e >= $a_e)
	         ||  ($b_b <= $a_b && $b_e >= $a_b));
	#print "got through step two FIX_THIS\n";
	#print "bb:$b_b,be:$b_e,ab:$a_b,ae:$a_e\n";
	return 0;
}
#-----------------------------------------------------------------------------
sub build_G {
	my $flag = shift;

	foreach my $scaff (keys %{$EX{$flag}}){
		foreach my $tid (keys %{$EX{$flag}{$scaff}}){

			my $exons = get_exons($flag, $scaff, $tid);

			my $type;
			my $gid;
			if     (defined($LU{$flag}{mRNA_to_gene}{$tid})){
				$type = 'mRNA';
				$gid = $LU{$flag}{mRNA_to_gene}{$tid}[0];	
			}
			elsif (defined($LU{$flag}{ncRNA_to_gene}{$tid})){
				$type = 'ncRNA';
				$gid = $LU{$flag}{ncRNA_to_gene}{$tid}[0];
			}

			$G{$flag}{$scaff}{$gid}{x} = $exons if defined $gid;

			$G{$flag}{$scaff}{$gid}{t} = $type  if defined $gid;
			print "GID:$gid\n";
			print "TID:$tid\n";
		    }
	}
}
#-----------------------------------------------------------------------------
sub get_exons {
	my $flag  = shift;
	my $scaff = shift;
	my $tid   = shift;


	my $t_exons = $EX{$flag}{$scaff}{$tid}; 
	my @exons;
	foreach my $e (@{$t_exons}){
		push(@exons, $e);

	}
	return \@exons;

}
#-----------------------------------------------------------------------------
sub ncRNA {
	my $term = shift;

	if    ($term eq 'ncRNA'){
		return 1;
	}
	elsif ($term eq 'transcript'){
		return 1;
	}
        elsif ($term eq 'rRNA'){
                return 1;
        }
        elsif ($term eq 'miRNA'){
                return 1;
        }
        elsif ($term eq 'tRNA'){
                return 1;
        }
        elsif ($term eq 'lincRNA'){
                return 1;
        }
        elsif ($term eq 'snoRNA'){
                return 1;
        }
        elsif ($term eq 'piRNA'){
                return 1;
        }
        elsif ($term eq 'snRNA'){
                return 1;
        }
	else {
		return 0;
	}

	
	
}
#-----------------------------------------------------------------------------
sub parse {

	my $file = shift;	
	my $flag = shift;
	
	return unless $file && $flag;
	print STDERR "NOW PARSING:$file $flag\n";

	my $fh = new FileHandle;
	   $fh->open($file);

	my $i = 0;
	while (defined(my $line = <$fh>)){
		chomp($line);
		#skip appended fasta sequence
		last if $line =~ /^\#\#FASTA/;
		#skip comment lines
		next if $line =~ /^\#/;
		
		my @stuff = split(/\t/, $line);

                #if -c is used only lines matching the specified
                #scaffold will be used 
		if ($opt_c){
		    next unless $stuff[0] eq $opt_c;
		}

		my $scaffold = trim(shift(@stuff));

		$S{$scaffold}{$flag}++;

		if    ($stuff[1] eq 'gene'){
			add_gene($scaffold, $flag, \@stuff);
		}
		elsif ($stuff[1] eq 'mRNA'){
			add_mrna($scaffold, $flag, \@stuff);
		
		}
                elsif (ncRNA($stuff[1])){
                        add_ncRNA($scaffold, $flag, \@stuff);
                }
		elsif ($stuff[1] eq 'exon'){
			add_exon($scaffold, $flag, \@stuff);
		}
		elsif ($stuff[1] eq 'contig'){
			add_contig($scaffold, $flag, \@stuff);
		}
		my $tot_scaffold = keys %S;

		print STDERR "TOT SCAFFOLD PROCESSED:$tot_scaffold\n" if $i/100000 == int $i/100000;
		$i++;

	}

	$fh->close();

}
#-----------------------------------------------------------------------------
sub tv {
	my $key   = shift;
	my $field = shift;

	my @tv = split(/\;/, $field);


	my ($tag, $value);
	foreach my $tv (@tv){
		next unless $tv =~ /$key/i;
		($tag, $value) = split(/\=/, $tv);
		last;
	}


	
	return $value;
}
#-----------------------------------------------------------------------------
sub mike_tv{
    my $key   = shift;
    my $field = shift;

    my @tv = split(/\;/, $field);
    my @values;

    my ($tag, $value);
    foreach my $tv (@tv){
	next unless $tv =~ /$key/i;
	($tag, $value) = split(/\=/, $tv);
	@values = split(/,/, $value);
	last;
    }



    return \@values;

}
#-----------------------------------------------------------------------------
sub add_contig {
        my $scaff = shift;
        my $flag  = shift;
        my $stuff = shift;


	$MU{$flag}{t_len} += $stuff->[3];
       
	   
}
#-----------------------------------------------------------------------------
# modified by Mike to allow alternative splicing
sub add_exon {
    my $scaff = shift;
    my $flag  = shift;
    my $stuff = shift;
    
    if ($stuff->[7] =~ /\S+,\S+/){
	my $parent_array = mike_tv('Parent',$stuff->[7]);
	my $eID          = tv('ID', $stuff->[7]);
	foreach my $parent (@{$parent_array}){
	    push(@{$LU{$flag}{exon_to_trans}{$eID}}, $parent);
	    push(@{$LU{$flag}{trans_to_exon}{$parent}}, $eID);
	    
	    push(@{$EX{$flag}{$scaff}{$parent}}, [$stuff->[2], $stuff->[3]]);
	    $MU{$flag}{num_ex}++;
	    load_EXLEN($scaff, $flag, $parent, $stuff);#MIKE#             
	}

    }

    else{
	my $parent = tv('Parent',$stuff->[7]);
	my $eID    = tv('ID', $stuff->[7]);
	
	push(@{$LU{$flag}{exon_to_trans}{$eID}}, $parent);
	push(@{$LU{$flag}{trans_to_exon}{$parent}}, $eID);
	
	push(@{$EX{$flag}{$scaff}{$parent}}, [$stuff->[2], $stuff->[3]]);
	$MU{$flag}{num_ex}++;
	load_EXLEN($scaff, $flag, $parent, $stuff);#MIKE#   
	}
}
#-----------------------------------------------------------------------------
#MIKE#
sub load_EXLEN{
    my $scaff = shift;
    my $flag  = shift;
    my $tid = shift;
    my $stuff = shift;
   
    push(@{$EXLEN{$flag}{$scaff}{$tid}}, $stuff->[2]);
    push(@{$EXLEN{$flag}{$scaff}{$tid}}, $stuff->[3]);
}
#-----------------------------------------------------------------------------
#MIKE#
sub ex_path_u{
    foreach my $flag (keys %EXLEN){
	foreach my $scaff (keys %{$EXLEN{$flag}}){
	    foreach my $tid (keys %{$EXLEN{$flag}{$scaff}}){
				
		ex_len($scaff, $flag, $tid, 'u');    
		in_len($scaff, $flag, $tid, 'u');    
	    }
	}
    }
}
#-----------------------------------------------------------------------------
#MIKE#
sub ex_path_int{
    foreach my $scaff (keys %I){
        foreach my $flag (keys %{$I{$scaff}{INT}}){
            foreach my $gid (keys %{$I{$scaff}{INT}{$flag}}){
		my $tids = get_tid($flag, $gid);
                foreach my $tid (@{$tids}){
		    $MI{$flag}{num_mRNA}++;    
		    ex_len($scaff, $flag, $tid, 'i');
		    in_len($scaff, $flag, $tid, 'i');
		}        
	    }
        }
    }
}
#-----------------------------------------------------------------------------
sub get_tid{
    my $flag = shift;
    my $gid  = shift;

    return unless defined(@{$LU{$flag}{gene_to_mRNA}{$gid}});

    my @tids = @{$LU{$flag}{gene_to_mRNA}{$gid}};
    return \@tids;
}
#-----------------------------------------------------------------------------
#MIKE#
sub ex_len{
    my $scaff = shift;
    my $flag  = shift;
    my $tid   = shift;
    my $x     = shift;
    
    my @exons = sort(@{$EXLEN{$flag}{$scaff}{$tid}});
    while (@exons >=2){
	my $exlen = $exons[1] - $exons[0] + 1;
	if ($x eq 'u'){	
	    $MU{$flag}{t_len_ex} += $exlen;
	    push(@{$LEN_U{'exons'}{$flag}}, $exlen);
	}
	if ($x eq 'i'){
	    $MI{$flag}{t_len_ex} += $exlen;
	    $MI{$flag}{num_ex}++;
	    push(@{$LEN_I{'exons'}{$flag}}, $exlen);    
	}	
	my $removed1 = shift(@exons);
	my $removed2 = shift(@exons);
    }
}
#-----------------------------------------------------------------------------
#MIKE#
sub in_len{
    my $scaff = shift;
    my $flag  = shift;
    my $tid   = shift;
    my $x     = shift;

    my @introns = sort(@{$EXLEN{$flag}{$scaff}{$tid}});
    my $remove_beg = shift(@introns);
    while (@introns >=2){
	my $inlen = $introns[1] - $introns[0] - 1;
	if ($x eq 'u'){
	    $MU{$flag}{t_len_in} += $inlen;
	    push(@{$LEN_U{'ints'}{$flag}}, $inlen);
	    $MU{$flag}{num_in}++;
	    #print "tid:$tid intlen:$inlen\n"; #devel line
	}
	if ($x eq 'i'){
	    $MI{$flag}{t_len_in} += $inlen;
	    push(@{$LEN_I{'ints'}{$flag}}, $inlen);
	    $MI{$flag}{num_in}++;
	}
	my $removed1 = shift(@introns);
	my $removed2 = shift(@introns);
    }
}
#-----------------------------------------------------------------------------
sub add_mrna {
        my $scaff = shift;
        my $flag  = shift;
        my $stuff = shift;

        my $parent = tv('Parent',$stuff->[7]);
	my $mid     = tv('ID',$stuff->[7]);
	
	$MU{$flag}{num_mRNA}++;
        
	push(@{$LU{$flag}{mRNA_to_gene}{$mid}}, $parent);
        push(@{$LU{$flag}{gene_to_mRNA}{$parent}}, $mid);

}
#-----------------------------------------------------------------------------
sub add_ncRNA {
        my $scaff = shift;
        my $flag  = shift;
        my $stuff = shift;

        my $parent = tv('Parent',$stuff->[7]);
        my $tid     = tv('ID',$stuff->[7]);

        push(@{$LU{$flag}{ncRNA_to_gene}{$tid}}, $parent);
        push(@{$LU{$flag}{gene_to_ncRNA}{$parent}}, $tid);


}
#-----------------------------------------------------------------------------
sub add_gene {
	my $scaff = shift;
	my $flag  = shift;
	my $stuff = shift;

	my $name = tv('ID',$stuff->[7]);

	$G{$flag}{$scaff}{$name}{s} = $stuff->[5];
	$G{$flag}{$scaff}{$name}{b} = $stuff->[2];
	$G{$flag}{$scaff}{$name}{e} = $stuff->[3];
	
}
#-----------------------------------------------------------------------------
sub trim {
	my $s = shift;

	$s =~ s/\..+$// if $opt_t;

	return $s;
}
#-----------------------------------------------------------------------------
#MIKE#
sub median{
    my $array = shift;
    @{$array} = sort {$a <=> $b} @{$array};
    my $count = scalar (@{$array});
    if ($count % 2) {
	return $array->[int($count/2)];
    } 
    else {
	#print " first thing: $array->[$count/2]\n";
	#print " second thing: $array->[$count/2 - 1]\n";
	return ($array->[$count/2] + $array->[$count/2 - 1]) / 2;
    }
    #PostData($array);
}
#-----------------------------------------------------------------------------
#MIKE#
sub max{
    my $array = shift;
    my @arrayb = sort {$a <=> $b} @{$array};
    return (pop @arrayb);

}
#-----------------------------------------------------------------------------
#MIKE#
sub median_exons{
    my $file = shift;
    my $flag = shift;
    my @exons;
    
    if ($flag eq 'u'){
	foreach my $chr (keys %{$EX{$file}}){
	    foreach my $mRNA (keys %{$EX{$file}{$chr}}){
		my $x = @{$EX{$file}{$chr}{$mRNA}};
		push (@exons, $x);
	    }
	}
    }
    if ($flag eq 'i'){
	foreach my $scaff (keys %I){
	    foreach my $gid (keys %{$I{$scaff}{INT}{$file}}){
		my $tids = get_tid($file, $gid);
		foreach my $mRNA (@{$tids}){    
		    my $x = @{$EX{$file}{$scaff}{$mRNA}};
		    push (@exons, $x);
		}    
	    }
	}
    }
    #PostData(\@exons);
    return median(\@exons);
}
#-----------------------------------------------------------------------------
