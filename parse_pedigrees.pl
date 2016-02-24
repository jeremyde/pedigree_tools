#!/usr/bin/perl

=pod

=head1 NAME

extract_pedigree_names.pl

=head1 SYNOPSIS

extract_pedigree_names.pl -i pedigrees.csv > names.txt

=head1 DESCRIPTION

This script extracts unique names from pedigrees.

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut

use strict;
use Getopt::Std;
use lib "$ENV{PWD}/lib";
use Bio::GeneticRelationships::ParsePedigree;

our ($opt_i, $opt_h);
getopts('i:h');
if ($opt_h){
    help();
    exit;
}
if (!$opt_i) {
    print STDERR "\nInput filename required\n\n\n";
    help();
}

my $line=0;
open FILE, "<", $opt_i or die "No such file $opt_i";
while (<FILE>) {
    $line++;
    chomp $_;
#    $_ =~ y/\n//d;
#    $_ =~ s/\m//d;
    $_ =~ s/\r//g;
    my @row =  split('\t', $_);
    my $accession = $row[0];
    my $pedigree = $row[1];
    #chomp $accession;
    #chomp $pedigree;
    my $pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $pedigree, accession => $accession);
    my $pedigree_objs = $pedigree_parse->parse_pedigrees();

    my $cross_data_ref = $pedigree_parse->get_cross_data(); 

    #if ($pedigree_parse->is_root()) {
#	print "$line\t$accession"."\t".$pedigree."\t".$pedigree."\t".$pedigree."\t".$pedigree."\t0\tA\n";
#	next;
#    }

    if ($pedigree_parse->has_parse_error()) {
	print "$line\t$accession\tParse error\t$pedigree\n";
    } elsif (!$cross_data_ref){
	print "$line\t$accession\tParse error\t$pedigree\n";
    } else {
	my %cross_data = %$cross_data_ref;

	foreach my $cross_key (keys %cross_data) {
	    my $cross_info = $cross_data{$cross_key};
	    print "$line\t$accession"."\t".$cross_info->{'progeny'}."\t".$cross_info->{'progeny_pedigree'}."\t".$cross_info->{'parent_a'}."\t".$cross_info->{'parent_b'}."\t".$cross_info->{'level'}."\t".$cross_info->{'a_or_b'}."\n";
#		print "prog(".$cross_info->{'progeny_pedigree'}.")\n";
#		print "prog(".$cross_info->{'progeny'}.")\n";
	}
    }
}
close FILE;


sub help {
  print STDERR <<EOF;
  $0:

    Description:
     This script extracts unique names from pedigrees.

    Usage:

      extract_pedigree_names.pl -i pedigrees.csv > names.txt

    Flags:

      -i <input_pedigree_file>      input pedigree file (mandatory)

      -h <help>                     help

EOF
exit (1);
}
