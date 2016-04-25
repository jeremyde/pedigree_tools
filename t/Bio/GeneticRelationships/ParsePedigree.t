# A test for parsing pedigrees
## Jeremy D. Edwards (jde22@cornell.edu) 2016

use strict;
use warnings;
use Test::More "no_plan";

use Data::Dumper;

BEGIN {require_ok('Moose');}
BEGIN {use_ok('Bio::GeneticRelationships::ParsePedigree');}

my $pedigree_parse;
my $pedigree_objs;

my $test_accession='testaccession';
my %test_pedigrees;

#pedigrees to test
$test_pedigrees{'simple'}='accession_1/accession_2';
$test_pedigrees{'2_deep_right'}='accession_1//accession_2/accession_3';
$test_pedigrees{'2_deep_left'}='accession_1/accession_2//accession_3';
$test_pedigrees{'3_deep_slashes'}='accession_1//accession_2/accession_3///accession_3';
$test_pedigrees{'3_deep_number'}='accession_1/2/accession_2/accession_3/3/accession_3';
$test_pedigrees{'3_deep_parentheses'}='(accession_1/(accession_2/accession_3))/accession_3';
$test_pedigrees{'named_cross_left_name_first'}='namevar(accession_1/2/accession_2/accession_3)/3/accession_3';
$test_pedigrees{'named_cross_left_name_last'}='(accession_1/2/accession_2/accession_3)namevar/3/accession_3';
$test_pedigrees{'backcross'}='accession_1/accession_2//accession_3*4///accession_4';
$test_pedigrees{'mixed'}='accession_1/3/accession_2/accession_3//(accession_4//accession_5/accession_6/3/accession_7)';


my $cross_data_ref;
my $parent_a_name;
my $parent_b_name;
my $parent_b_name2;

# Test simple pedigree
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'simple'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
is($cross_data_ref->{'testaccession'}->{'parent_a'},'accession_1');
is($cross_data_ref->{'testaccession'}->{'parent_b'},'accession_2');

# Test pedigree two deep on right side
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'2_deep_right'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
is($cross_data_ref->{'testaccession'}->{'parent_a'},'accession_1');
ok($parent_b_name = $cross_data_ref->{'testaccession'}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree two deep on left side
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'2_deep_left'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
is($cross_data_ref->{'testaccession'}->{'parent_b'},'accession_3');
ok($parent_a_name = $cross_data_ref->{'testaccession'}->{'parent_a'});
is($cross_data_ref->{$parent_a_name}->{'parent_a'},'accession_1');
is($cross_data_ref->{$parent_a_name}->{'parent_b'},'accession_2');

# Test pedigree three deep with three slashes
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'3_deep_slashes'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_a_name = $cross_data_ref->{'testaccession'}->{'parent_a'});
ok($parent_b_name = $cross_data_ref->{$parent_a_name}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree three deep with number between slashes
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'3_deep_number'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_a_name = $cross_data_ref->{'testaccession'}->{'parent_a'});
ok($parent_b_name = $cross_data_ref->{$parent_a_name}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree three deep with parentheses
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'3_deep_parentheses'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_a_name = $cross_data_ref->{'testaccession'}->{'parent_a'});
ok($parent_b_name = $cross_data_ref->{$parent_a_name}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree with a named variety and corresponding pedigree in parentheses with name first
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'named_cross_left_name_first'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_b_name = $cross_data_ref->{'namevar'}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree with a named variety and corresponding pedigree in parentheses with name last
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'named_cross_left_name_last'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_b_name = $cross_data_ref->{'namevar'}->{'parent_b'});
is($cross_data_ref->{$parent_b_name}->{'parent_a'},'accession_2');
is($cross_data_ref->{$parent_b_name}->{'parent_b'},'accession_3');

# Test pedigree containing backcross
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'backcross'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
ok($parent_a_name = $cross_data_ref->{'testaccession'}->{'parent_a'});
ok(my $parent_BC2_name = $cross_data_ref->{$parent_a_name}->{'parent_a'});
is($cross_data_ref->{$parent_a_name}->{'parent_b'},'accession_3'); #error
ok(my $parent_BC1_name = $cross_data_ref->{$parent_BC2_name}->{'parent_a'});
is($cross_data_ref->{$parent_BC2_name}->{'parent_b'},'accession_3'); #error

# Test pedigree with a mixed use of parentheses and an equally deep pedigree inside the parentheses as outside
ok($pedigree_parse = Bio::GeneticRelationships::ParsePedigree->new( pedigree => $test_pedigrees{'mixed'}, accession=> $test_accession));
ok($pedigree_parse->parse_pedigrees());
ok($cross_data_ref = $pedigree_parse->get_cross_data()); 
is($cross_data_ref->{'testaccession'}->{'parent_a'},'accession_1');
ok($parent_b_name = $cross_data_ref->{'testaccession'}->{'parent_b'});
ok($parent_b_name2 = $cross_data_ref->{$parent_b_name}->{'parent_b'});
is($cross_data_ref->{$parent_b_name2}->{'parent_b'},'accession_7');
