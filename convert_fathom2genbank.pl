#!/usr/bin/perl -w
use strict;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "Usage: 
convert_fathom2genbank.pl uni.ann uni.dna number_of_genes_to_use(or use 'all')

Note:
The two files genome.ann and genome.dna is recommended to derive from SNAP's
fathom package. Run:
fathom genome.ann genome.dna -categorize 1000
And take uni.dna and uni.ann in the output.



";

die $usage unless $ARGV[2];

my $ann = $ARGV[0];
my $dna = $ARGV[1];
my $total_num = $ARGV[2];

$total_num = 1000000 if $total_num eq 'all';

my $annotation = get_ann($ann);
my $seq = get_seq($dna);

my $count = 0;
while (my ($id, $gene) = each %$annotation) {
	my $length = length($seq->{$id});

	print "LOCUS       $id   $length bp  DNA\n";

	print "FEATURES             Location/Qualifiers\n";

	print "     source          1..$length\n";
	my @cors;

	foreach my $cor (@{$gene->{exons}}) {
		push @cors, join '..', @$cor;
	}
	my $cors = join ',', @cors;
	if (scalar @cors >1) {
		$cors = 'join('.$cors.')';
	}
	
	if ($gene->{str} eq '-') {
		$cors = 'complement('.$cors.')';
	}

	my $join = $cors;
	print "     CDS             ";

	my @join = split /,/, $join;
	my $line_length = 0;
	for (my $i = 0; $i<=$#join; $i++) {
		if ($line_length + length($join[$i]) >59) {
			print "\n                     ";
			$line_length = 0;
		}
		print $join[$i];
		if ($i!=$#join) {
			print ",";
		}
		else {
			print "\n";
		}
		$line_length += length($join[$i]);
		
	}
	print_seq($seq->{$id});

	print "//\n";
	if (++$count >=$total_num) {
		last;
	}
}
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub print_seq {
	my $seq = shift;

	print "ORIGIN";
	my $header = 1;

	my @seq = $seq =~ /(.{1,10})/g;

	my $count = 1;
	foreach my $column (@seq) {
		if ( ($count-1)/6 == int (($count-1)/6)) {
			my $space_l = 9-length($header);
			print "\n";
			for (my $i=1; $i<=$space_l;$i++) {
				print " ";
			}

			print "$header ";
			$header += 60;
		}
		
		print lc($column)." ";
		$count ++;
	}
	print "\n";
}
#-----------------------------------------------------------------------------
sub get_ann {
	my $file = shift;

	my %ann;
	open FH, $file;

	my $id;
	while (<FH>) {
		chomp;

		if (/^>(\S+)/) {
			$id = $1;
		}
		else {
			my @i = split /\t/, $_;

			if ($i[1] < $i[2]) {
				$ann{$id}->{str} = '+';
			}
			else {
				$ann{$id}->{str} = '-';
			}
                        my ($l , $r) = ($i[1] < $i[2])? ($i[1], $i[2]):
                                                ($i[2], $i[1]);
                        push @{$ann{$id}->{exons}}, [$l, $r];
		}
	}
	close FH;

	return \%ann;
}	
#-----------------------------------------------------------------------------
sub get_seq {
	my $file = shift;
	my %seq;

	open FH, $file ;
	my $id;

	while (<FH>) {
                chomp;

                if (/^>(\S+)/) {
                        $id = $1;
                }
		else {
			$seq{$id} .= $_;
		}
	}
	close FH;
	return \%seq;
}
#-----------------------------------------------------------------------------
