#!/usr/bin/perl -w 
use strict;
#use lib ('/home/mcampbell/lib');
#use PostData;
#use Getopt::Std;
#use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
#getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
This script takes a populated MAKER opts file and outputs an editable 
maker_evm.ctl file. 
 
make_evm_weights.pl maker_opts.ctl > maker_evm.ctl
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
my %DATA;
my %LU;

parse($FILE);
make_weights();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_possible_types{
    my $type = shift;
    if ($type eq 'est_pass'){
	return ('est2genome', 'blastn');
    } 
    if ($type eq 'altest_pass'){
	return ('cdna2genome', 'tblastx');
    } 
    if ($type eq 'protein_pass'){
	return ('protein2genome', 'blastx');
    }
    if ($type eq 'model_pass'){
	return ('maker');
    } 
    if ($type eq 'pred_pass'){
	return ('snap', 'augustus', 'fgenesh', 
		'genemark', 'pred_gff', 'model_gff');
    } 
	

}
#-----------------------------------------------------------------------------
sub get_gff_labels{
    my $type = shift;
    my $maker_gff_status = shift;
    my %sources;
    my @labels;
    my @gff3_files;
    my @posibilities; #needs a better name it is possible sources from maker_gff
    if (!$maker_gff_status){
	@gff3_files = split(/\,/, $DATA{$type});
    }
    if ($maker_gff_status){
	@posibilities = get_possible_types($type);
	@gff3_files = split(/\,/, $DATA{'maker_gff'});
    }
    foreach my $gff3_file (@gff3_files){
	$gff3_file =~ s/\s//g;

        my $fh = new FileHandle;
	$fh->open($gff3_file);
	
	while (defined(my $line = <$fh>)){
	    chomp($line);
	    last if $line =~ /^\#\#FASTA/;
	    next if $line =~ /^\#/;
	    my @stuff = split(/\t/, $line);
	    $sources{$stuff[1]}++;
	}
	$fh->close();
    }
    
    foreach my $x (keys %sources){
	if ($maker_gff_status){

	    foreach my $pos (@posibilities){
            #this loop is to deal with source gotchas when 
	    #using maker_gff
		if ($x =~ /^$pos/){
		    my $y = lc $x;
		    push(@labels, "$y");
		}
		else {next;}
	    }
	}
	else{
	    my $y = lc $x;
	    push(@labels, "$type:$y");
	}
    }
    return(@labels);
}
#-----------------------------------------------------------------------------
sub get_source_label{
    my $type = shift;
    my @labels;
    return('est2genome','blastn') if $type eq 'est';
    return('cdna2genome','tblastx') if $type eq 'altest';
    return('protein2genome','blastx') if $type eq 'protein';
    return('snap') if $type eq 'snaphmm';
    return('augustus') if $type eq 'augustus_species';
    return('fgenesh') if $type eq 'fgenesh_par_file';
    return('genemark') if $type eq 'gmhmm';
    if ($type eq 'est_gff' || 
	$type eq 'altest_gff' ||
	$type eq 'protein_gff'||
	$type eq 'model_gff'||
	$type eq 'pred_gff'){ 
    @labels = get_gff_labels($type, '0');
    }
    if ($type eq 'est_pass' ||
	$type eq 'altest_pass' ||
	$type eq 'protein_pass' ||
	$type eq 'model_pass' ||
	$type eq 'pred_pass'){
    @labels = get_gff_labels($type, '1');
    }
    return(@labels) if $type eq 'est_gff';
    return(@labels) if $type eq 'altest_gff';
    return(@labels) if $type eq 'protein_gff';
    return(@labels) if $type eq 'pred_gff';
    return(@labels) if $type eq 'model_gff';

    return(@labels) if $type eq 'est_pass';
    return(@labels) if $type eq 'altest_pass';
    return(@labels) if $type eq 'protein_pass';
    return(@labels) if $type eq 'model_pass';
    return(@labels) if $type eq 'pred_pass';
}
#-----------------------------------------------------------------------------
sub gene_finder_input{
    my $type = shift;
    my $category = shift;
    my $masked = shift;
    my $unmask = shift;
    $masked = 0 if $type eq 'gmhmm';
    my $source_label = get_source_label($type);

    my @abinit_files = split(/\,/, $DATA{$type});
    my $has_had_tag = 0;
    print "#-----Abinitio Prediction weights\n$category=1\n" unless defined($LU{$category});
    $LU{$category}=1;
    foreach my $abinit_file (@abinit_files){
	my ($file, $tag) = split(/\:/, $abinit_file);
	if (defined($tag)){
	    $tag =~ s/\s//g;
	    $tag = ':'.$tag;
	    if ($masked){ 
		print "$category".":"."$source_label"."_masked$tag"."="."1\n";
		$has_had_tag = 1;
	    }
	    if(!$masked || $DATA{'unmask'}){
		print "$category".":"."$source_label$tag"."="."1\n";
		$has_had_tag = 1;
	    }
	}
	if(!defined($tag) && $has_had_tag == 1){
	    die("FAILURE:every $type file needs to have a tag\n");
	}
	else {next;} 
    }
    if($has_had_tag == 0 && $masked){
	print "$category".":"."$source_label"."_masked"."="."1\n";
    }
    if($has_had_tag == 0 && (!$masked || $DATA{'unmask'})){
	print "$category".":"."$source_label"."="."1\n";
    }
}
#-----------------------------------------------------------------------------
sub tagged_input{
    my $type = shift;
    my $category = shift;
    my @source_label = get_source_label($type);
    my @files = split(/\,/, $DATA{$type});
    my $has_had_tag = 0;
    my $counter=0;
    my $tag_count=0;
    if ($category eq 'evmtrans'){
	print "#-----Transcript weights\n$category=1\n" unless defined($LU{$category});
	$LU{$category}=1;
    }
    elsif ($category eq 'evmprot'){
	print "#-----Protein weights\n$category=1\n" unless defined($LU{$category});
	$LU{$category}=1;
    }


#    print $source_label[0]."\n";
    foreach my $file (@files){
	$counter++;
	my ($path, $tag) = split(/\:/, $file);
	if (defined($tag)){
	    $tag =~ s/\s//g;
	    $tag = ":".$tag;
	    $tag_count++;
	    if (@source_label > 1){
		foreach my $x (@source_label){
		    print "$category".":"."$x$tag"."="."1\n";
		    $has_had_tag = 1;
		}
	    }
	    else{print "$category".":"."$source_label[0]$tag"."="."1\n";
		 $has_had_tag = 1;
		 #$tag_count++;
	    }
	}
	if(!defined($tag) && $has_had_tag == 1){
	 #   print "$type:";
	    die("FAILURE:every $type file needs to have a tag\n");
	}
	else {next;} 
    }
    if($has_had_tag == 0 && ($counter != $tag_count)){
	if (@source_label >= 1){
	    foreach my $x (@source_label){
		print "$category".":"."$x"."="."1\n";
	    }
	}
	else{die("$tag_count:FAILURE:every $type file needs to have a tag\n");}
    }
}
#-----------------------------------------------------------------------------
sub make_weights{
    my $masked = 0;
    my $unmask = $DATA{'unmask'};

#if there is repeat masking the gene finders labels will change
    if (defined($DATA{model_org}) ||
	defined($DATA{rmlib}) ||
	defined($DATA{repeat_protein}) ||
	defined($DATA{rm_gff}) ||
	$DATA{rm_pass}==1){
	$masked = 1;
    }
#these pass through options will take some investigating
    if (defined($DATA{'maker_gff'})){
	if($DATA{'est_pass'}){
	    tagged_input('est_pass', 'evmtrans');
	}
	if($DATA{'altest_pass'}){
	    tagged_input('altest_pass', 'evmtrans');
	}
	if($DATA{'protein_pass'}){
	    tagged_input('protein_pass', 'evmprot');
	}
	if($DATA{'model_pass'}){
	    tagged_input('model_pass', 'evmab');
	}
	if($DATA{'pred_pass'}){
	    tagged_input('pred_pass', 'evmab');
	}
    }
    if (defined($DATA{'est'})){
	tagged_input('est', 'evmtrans');
    }
    if (defined($DATA{'altest'})){
	tagged_input('altest', 'evmtrans');
    }
    if (defined($DATA{'est_gff'})){
	tagged_input('est_gff','evmtrans');
    }
    if (defined($DATA{'altest_gff'})){
	tagged_input('altest_gff','evmtrans');
    }
    if (defined($DATA{'protein'})){
	tagged_input('protein', 'evmprot');
    }
    if (defined($DATA{'protein_gff'})){
	tagged_input('protein_gff','evmprot');
    }
    if (defined($DATA{'snaphmm'})){
	gene_finder_input('snaphmm','evmab',$masked, $unmask);
    }
    if (defined($DATA{'gmhmm'})){
	gene_finder_input('gmhmm','evmab',$masked, $unmask);
    }
    if (defined($DATA{'augustus_species'})){
	gene_finder_input('augustus_species','evmab',$masked, $unmask);
    }
    if (defined($DATA{'fgenesh_par_file'})){
	gene_finder_input('fgenesh_par_file','evmab',$masked, $unmask);
    }
    if (defined($DATA{'pred_gff'})){
	tagged_input('pred_gff','evmab');
    }
    if (defined($DATA{'model_gff'})){
	tagged_input('model_gff','evmab');
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	return if $line =~ /^\#-----MAKER Behavior Options/;
	next if $line =~ /^\#/;
	$line =~ s/\s*\#.*//;
	my @kvp = split(/\=/, $line);
	$DATA{$kvp[0]} = $kvp[1] if defined($kvp[1]) && $kvp[1] !~ /^\s$/;

    }
    $fh->close();
}
#-----------------------------------------------------------------------------

