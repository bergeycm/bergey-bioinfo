#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

use lib '/home/cmb433/local_perl/';

use Bio::Restriction::Analysis;
use Bio::Restriction::EnzymeCollection;
use Bio::SeqIO;

# ----------------------------------------------------------------------------------------
# --- Get user input
# ----------------------------------------------------------------------------------------

my $usage = "Usage: perl $0 --input_fasta input.fa [--min_length 100000]\n";

my $input_fasta;
my $min_length = 0;
my $verbose;

GetOptions ("input_fasta=s" => \$input_fasta,
			"min_length:i"  => \$min_length,
			"verbose"       => \$verbose)
or die("ERROR: Problem with passed arguments.\n$usage\n");

if (!$input_fasta) {
	die "ERROR: Missing input FASTA file.\n$usage\n";
}

# Open and load DNA sequence in FASTA file
chomp $input_fasta;
my $seqIO_obj = Bio::SeqIO->new(-file => $input_fasta, -format => "fasta");

# Abort if results folder already exists
my $out_dir = $input_fasta . '_enzymes/';
if (-d $out_dir) {
	die "ERROR: Results directory [$out_dir] already exists. Aborting.\n";
}

# Create folder to hold results
mkdir $out_dir;

# ----------------------------------------------------------------------------------------
# --- Load data on enzymes
# ----------------------------------------------------------------------------------------

# Create collection of just Type II enzymes
my $complete_collection = Bio::Restriction::EnzymeCollection->new();
my $type_ii_collection  = Bio::Restriction::EnzymeCollection->new( -empty => 1 );
$type_ii_collection->enzymes( grep { $_->type() eq 'II' } $complete_collection->each_enzyme() );

# Could also be run on nonambiguous enzymes only:
# my $type_ii_nonambig = Bio::Restriction::EnzymeCollection->new( -empty => 1 );
# $type_ii_nonambig->enzymes( grep { $_->is_ambiguous == 0 } $type_ii_collection->each_enzyme() );

# Put these guys in arrays as Bio::Restriction::Enzyme objects
# One array is all Type II, second array (unused) is just the unambiguous Type II
my @enzymes = $type_ii_collection->each_enzyme();

# Could also be run on nonambiguous enzymes only:
# my @na_enzymes = $type_ii_nonambig->each_enzyme();

my $num_enz = scalar @enzymes;

# ----------------------------------------------------------------------------------------
# --- Do digestion for each sequence
# ----------------------------------------------------------------------------------------

# Keep track of how much sequence is analyzed
my $total_length_digested = 0;

while (my $seq_obj = $seqIO_obj->next_seq) {

	my $seq_string;
	if ($seq_obj->desc() ne "") {
		$seq_string = $seq_obj->desc();
	} else {
		$seq_string = $seq_obj->display_id();
	}
	print STDERR "Processing $seq_string...\n";
	
	# Skip if sequence length is less than minimum length
	my $seq_len = $seq_obj->length();
	if ($seq_len < $min_length) {
		print STDERR "\tSkipping sequence of length $seq_len.\n";
		next;
	}
	
	$total_length_digested += $seq_len;
	
	my $e_count = 1;

	# Start restriction enzyme analysis using our subset of enzymes
	my $ra = Bio::Restriction::Analysis->new(	-seq=>$seq_obj, 
												-enzymes=>$type_ii_collection);

	# For each enzyme, cut DNA and save fragments
	
	for my $e (@enzymes) {
		
		my @cut_positions = $ra->positions($e->name());
		
		# Abort if there's only one fragment for this enzyme and for this contig
		# Avoid an unneccessary opening of the file
		if (scalar @cut_positions == 1) {
			next;
		}
		
		push @cut_positions, 0;
		push @cut_positions, length $seq_obj->seq();
		
		@cut_positions = sort { $a <=> $b } @cut_positions;
		
		open (ENZ, '>>' . $input_fasta . '_enzymes/Enz_' . $e->name() . '.txt')
			or die "Could not open enzyme file for " . $e->name() . "\n";

		# Lock file, since many instances of this program will be writing to it
		##flock(ENZ, 2);
		
		# Write fragment info to file, but skip the first and last fragment
		for (my $cp = 1; $cp < scalar (@cut_positions) - 1; $cp++) {
		
			my $start = 1 + $cut_positions[$cp];
			my $end   = $cut_positions[$cp + 1];
			my $frag_len = 1 + $end - $start;
		
			print ENZ $e->name() . "\t";					# Print enzyme
			if ($seq_obj->desc() ne "") {					# Print contig description
				print ENZ $seq_obj->desc() . "\t";
			} else {
				print ENZ $seq_obj->display_id() . "\t";	# ...or contig ID
			}
			print ENZ $start . "\t";						# Print start bp
			print ENZ $end   . "\t";						# Print end bp
			print ENZ $frag_len;							# Print fragment length
			print ENZ "\n";
		}
		
		close (ENZ);
		print "\tFinished enzyme $e_count of $num_enz.\n" if $verbose;
		$e_count++;
	}

}

print STDERR "Finished simulated digestion.\n";
print STDERR "\tTotal sequence digested:\t$total_length_digested bp\n";
print STDERR "\tTotal number of enzymes:\t$num_enz\n";

exit(0);
