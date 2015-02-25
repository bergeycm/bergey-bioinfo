#!/usr/bin/env perl

use strict;
use warnings;
use LWP::Simple;
use XML::XPath;
use XML::XPath::XMLParser;

# Use DAS of UCSC to fetch specific sequence by its given chromosome position
# From here: https://www.biostars.org/p/6156/

my $chr  = shift;
my $pos  = shift;
my $size = shift;

my $usage = "Example: perl extract_seq_from_ucsc.pl 14 482780 1000\n";

if (! $size) {
	die "ERROR: You must pass three arguments: chr. num., position, and size.\n$usage";
	
}

chomp $size;

my $start = $pos - ($size/2);
my $end   = $pos + ($size/2);

# Figure out URL for the DAS server. Example:
# http://genome.ucsc.edu/cgi-bin/das/calJac3/dna?segment=chr14:482280,483280

my $URL_gene ="http://genome.ucsc.edu/cgi-bin/das/papAnu2/dna?segment=chr";
$URL_gene .= $chr . ":" . $start . "," . $end;

my $xml = get($URL_gene);

my $xp = XML::XPath->new(xml=>$xml);

my $nodeset = $xp->find('/DASDNA/SEQUENCE/DNA/text()'); # find all sequences
# there should be only one node, anyway:    
foreach my $node ($nodeset->get_nodelist) {

	my $seq = $node->getValue;
	$seq =~ s/\s//g; # remove white spaces
	print ">papAnu2_chr" . $chr . ":" . $start . "-" . $end . "\n";
	print $seq, "\n";
	
}
