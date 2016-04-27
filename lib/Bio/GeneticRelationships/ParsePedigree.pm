
=head1 NAME

Bio::GeneticRelationships::ParsePedigree - a module to parse pedigrees from strings.

=head1 USAGE

 my $pedigree_parse = Bio::GeneticRelationships->new({ pedigree => $pedigree} );
 my $validated = $pedigree_parse->validate(); #is true when the pedigree is valid
 $pedigree_parse->get_pedigree_objects();

=head1 DESCRIPTION

Parses a pedigree string. This module is intended to be used in independent loading scripts and interactive dialogs.

=head1 AUTHORS

 Jeremy D. Edwards (jde22@cornell.edu)

=cut

package Bio::GeneticRelationships::ParsePedigree;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Try::Tiny;
use Data::GUID qw(guid);

has 'pedigree' => (isa => 'Str', is => 'ro', predicate => 'has_pedigree', required => 1);
has 'accession' => (isa => 'Str', is => 'ro', predicate => 'has_accession', required => 1);
has 'cross_data' => (isa => 'HashRef[HashRef[Str]]', is => 'rw', predicate => 'has_cross_data');
has 'leaf_data' => (isa => 'HashRef[HashRef[Str]]', is => 'rw', predicate => 'has_leaf_data');
has 'parse_error' => (isa => 'Str', is => 'rw', predicate => 'has_parse_error');

sub parse_pedigrees {
    my $self = shift;
    my $pedigree_str = $self->get_pedigree();
    my $accession = $self->get_accession();
    print STDERR "\n\nProcessing $pedigree_str\n\n\n";
    my $parse_success = _split_pedigree($self, $pedigree_str, "A", 0, $accession);
    return 1;
}

sub is_root {
    my $self = shift;
    my $pedigree_str = $self->get_pedigree();
    my $accession = $self->get_accession();
    my $parse_success = _split_pedigree($self, $pedigree_str, "A", 0, $accession);
    if (!$parse_success) {
	return;
    }
    if ($parse_success = $pedigree_str) {
	return 1;
    } else {
	return;
    }
}

sub _mask_inside {
    my $pedigree_str = shift;
    my $parenth_open = 0;
    my $mask_pos = 0;
    my $was_masked;
    my $parenth_first;
    my $parenth_last;
    my $parenth_sets = 0;
    my @mask_chars = split('', $pedigree_str);
    for my $char (@mask_chars) {
	if ($char eq "(") {
	    $was_masked = 1;
	    #mask parenthesis inside parenthesis
	    if ($parenth_open > 0) {
		$mask_chars[$mask_pos] = 0;
	    }
	    $parenth_open++;
	    if (!$parenth_first) {
		$parenth_first = $mask_pos;
	    }
	} elsif ($char eq ")") {
	    # mask nested parentheses
	    if ($parenth_open > 1) {
		$mask_chars[$mask_pos] = 0;
	    }
	    # count non-nested sets of parentheses
	    if ($parenth_open == 1) {
		$parenth_sets++;
	    }
	    $parenth_open--;
	    $parenth_last = $mask_pos;
	} elsif ($parenth_open > 0) {
	    $mask_chars[$mask_pos] = 0; #mask using zeros
	}
	$mask_pos++;
    }
    $pedigree_str = join('',@mask_chars);
    return ($pedigree_str, $was_masked, $parenth_first, $parenth_last, $parenth_sets);
}

sub _split_pedigree {
    my $self = shift;
    my $pedigree_str = shift;
    my $a_or_b = shift;
    my $level = shift;
    my $unmasked_pedigree_str = $pedigree_str;
    my $accession = shift;
    my $parent_A;
    my $parent_B;
    my $longest = 0;
    my $next_longest = 0;
    my $longest_start;
    my $longest_end;
    chomp $pedigree_str;

    # remove spaces around slash
    $pedigree_str =~ s/\/ /\//g;
    $pedigree_str =~ s/ \//\//g;

    # detect mismatched parenthesis
    my $offset = 0;
    my $forward_count = 0;
    my $reverse_count = 0;
    my $forward_result = index($pedigree_str,'(',$offset);
    while ($forward_result != -1) {
	$forward_count++;
	$offset = $forward_result + 1;
	$forward_result = index($pedigree_str,'(',$offset);
    }
    $offset = 0;
    my $reverse_result = index($pedigree_str,')',$offset);
    while ($reverse_result != -1) {
	$reverse_count++;
	$offset = $reverse_result + 1;
	$reverse_result = index($pedigree_str,')',$offset);
    }
    if ($forward_count != $reverse_count) {
	print STDERR "Error:mismatched parenthesis in $pedigree_str\n";
	$self->set_parse_error('-1');
	return -1;
    }

    # convert string to array
    my @chars = split('', $pedigree_str);

    # remove parenthesis that begin and end string
    my $trim_done = -1;
    while ($trim_done == -1) {
	if ((substr($pedigree_str,0,1) eq '(') && (substr($pedigree_str,length($pedigree_str)-1,length($pedigree_str)) eq ')')) {
	    my $parenth_logic = 0;
	    my $logic_bad;
	    my $pedigree_str_trim = substr($pedigree_str, 1, length($pedigree_str) - 2);
	    my @chars_trim = split('', $pedigree_str_trim);
	    for my $char (@chars_trim) {
		if ($char eq '(') {
		    $parenth_logic++;
		}
		elsif ($char eq ')') {
		    $parenth_logic--;
		}
		if ($parenth_logic < 0) {
		    $logic_bad = 1;
		}
	    }
	    if (!$logic_bad){ 
		$pedigree_str = $pedigree_str_trim;
		$unmasked_pedigree_str = $pedigree_str_trim;
		@chars = @chars_trim;
	    } else {
		$trim_done = 1;
	    }
	} else {
	    $trim_done = 1;
	}
    }

    # mask inside parenthesis
    my $was_masked;
    my $parenth_first;
    my $parenth_last;
    my $parenth_sets;
    $unmasked_pedigree_str = $pedigree_str;
    ($pedigree_str, $was_masked, $parenth_first, $parenth_last, $parenth_sets) = _mask_inside($pedigree_str); 

    #check for cross name outside of parenthesis 
    if ($was_masked && $parenth_sets == 1) {
	my $name_side;
	my $name_pos;
	my $possible_cross_name;
	my $possible_named_cross;
	my $possible_named_cross_unmasked;
	if (substr($unmasked_pedigree_str,0,1) eq "(") {
	    $name_pos = $parenth_last;
	    $possible_cross_name = substr($unmasked_pedigree_str,$name_pos+1, (length($unmasked_pedigree_str)-$name_pos));
	    $possible_named_cross = substr($pedigree_str, 1, $name_pos-1);
	    $possible_named_cross_unmasked = substr($unmasked_pedigree_str, 1, $name_pos-1);
	}
	elsif (substr($unmasked_pedigree_str,length($unmasked_pedigree_str)-1,1) eq ")") {
	    $name_pos = $parenth_first;
	    $possible_cross_name = substr($unmasked_pedigree_str,0, $name_pos);
	    $possible_named_cross = substr($pedigree_str, $name_pos+1, (length($pedigree_str)-$name_pos)-2);
	    $possible_named_cross_unmasked = substr($unmasked_pedigree_str, $name_pos+1, (length($unmasked_pedigree_str)-$name_pos)-2);
	}
	
	#check possible cross name for characters indicating a cross and not a name
	if ($possible_cross_name) {
	    unless ($possible_cross_name =~ m/\//) {
		$accession = $possible_cross_name;
		$pedigree_str = $possible_named_cross;
		$unmasked_pedigree_str = $possible_named_cross_unmasked;
		#need to feed unmasked back through masking if using cross name
		($pedigree_str, $was_masked, $parenth_first, $parenth_last) = _mask_inside($unmasked_pedigree_str); 
	    }
	}
    }

    # set accession if none provided
    if (!$accession) {
	# a unique code to use as an accession name if one isn't provided
	my $guid_code = Data::GUID->new()->as_hex();
	$accession = $pedigree_str."_".$guid_code;
    }

    #find longest stretch of slashes
    #and check that there isn't another stretch of slashes of equal length
    my $current = 0;
    my $current_pos = 0;
    @chars = split('', $pedigree_str);
    for my $char (@chars) {
	$current_pos++;
	if ($char eq '/') {
	    $current++;
	    if ($current > $longest) {
		$next_longest = $longest;
		$longest = $current;
		$longest_start = $current_pos - ($current-1);
		$longest_end = $current_pos - 1;
	    }
	} else {
	    $current = 0;
	}
    }
    if ($longest > 0 && ($longest <= $next_longest)) {
	print STDERR "Pedigree structure error, split at equal levels in $pedigree_str\n";
	$self->set_parse_error('-1');
	return -1;
    }

    # search for slashes using the /3/ format
    my @slashes;
    my $masked_pedigree_str  = join( '' , @chars ) ;
    push @slashes, pos($masked_pedigree_str) while ($masked_pedigree_str =~ m/\/\d{1,2}\//g);
    foreach my $slash (@slashes) {
	my $slash_offset = 3;
	my $cross_str = substr($masked_pedigree_str, $slash-$slash_offset,3); 
	if (substr($cross_str,0,1) ne '/') {
	    $slash_offset = 4;
	    $cross_str = substr($masked_pedigree_str, ($slash-$slash_offset),4);
	}
	if (substr($cross_str,0,1) ne '/') {
	    print STDERR "Slash detection error in $pedigree_str\n";
	    $self->set_parse_error('-1');
	    return -1;
	}
	my $crossnum = substr($cross_str,1,length($cross_str)-2);
	if ($crossnum > $longest) {
	    $next_longest = $longest;
	    $longest = $crossnum;
	    $longest_start = $slash-$slash_offset+1;
	    $longest_end = ($slash-$slash_offset)+$slash_offset-1;
	}
    }

    # return string if nothing left to split
    if ($longest == 0) {
	my $leaf_return;
	$leaf_return=_process_leaf($self, $unmasked_pedigree_str, $accession);
	if (!$leaf_return) {
	    print STDERR "Pedigree structure error (leaf): $pedigree_str\n";
	    $self->set_parse_error('-1');
	    return -1;
	}
	return $leaf_return;
    }

    # Deal with backcrosses
    my $right_side = substr($unmasked_pedigree_str, $longest_end+1, length($unmasked_pedigree_str));
    my $left_side = substr($unmasked_pedigree_str, 0, $longest_start-1);
    # check if there are any slashes and if it ends with a backcross
    my $backcross_number;
    my $backcross_accession;
    if (!($right_side =~ m/\//) && ($right_side =~ m/\*\d{1,2}$/)) {
	my @split_backcross = split('\*', $right_side);
	$backcross_accession = $split_backcross[0];
	$backcross_number = $split_backcross[1];
	$backcross_number--;
	my $new_backcross_name;
	if ($backcross_number > 1) {
	    $new_backcross_name = "(($left_side)/$backcross_accession*".$backcross_number.")/$backcross_accession";
	    $longest_start = length($left_side)+length($backcross_accession)+length($backcross_number)+7;
	} else {
	    $new_backcross_name = "(($left_side)/$backcross_accession)/$backcross_accession";
	    $longest_start = length($left_side)+length($backcross_accession)+6;
	}
	$longest_end = $longest_start-1;
	$unmasked_pedigree_str=$new_backcross_name;
	($pedigree_str, $was_masked, $parenth_first, $parenth_last) = _mask_inside($unmasked_pedigree_str);
    } 
    elsif (!($left_side =~ m/\//) && ($left_side =~ m/\*\d{1,2}$/)) {
	my @split_backcross = split('\*', $left_side);
	$backcross_accession = $split_backcross[0];
	$backcross_number = $split_backcross[1];
	$backcross_number--;
	my $new_backcross_name;
	if ($backcross_number > 1) {
	    $new_backcross_name = $backcross_accession."/(".$backcross_accession."*".$backcross_number."/(".$right_side."))";
	    $longest_start = length($backcross_accession)+1;#or 2???
	} else {
	    $new_backcross_name = "$backcross_accession/($backcross_accession/($right_side))";
	    $longest_start = length($backcross_accession)+1;
	}
	$longest_end = $longest_start-1;
	$unmasked_pedigree_str=$new_backcross_name;
	($pedigree_str, $was_masked, $parenth_first, $parenth_last) = _mask_inside($unmasked_pedigree_str);
    }

    # split pedigree by parents 
    my $parent_ids;
    $parent_ids = _do_next_split_pedigree($self, $unmasked_pedigree_str, $longest_start, $longest_end, $level+1);
    if (!$parent_ids || $parent_ids == -1) {
	print STDERR "\n\nPedigree structure error in backcross $pedigree_str\n";
	$self->set_parse_error('-1');
	return -1;
    }
    $parent_A = @$parent_ids[0];
    $parent_B = @$parent_ids[1];

    #write to hash
    my %cross_info;
    #store values in hash
    $cross_info{'progeny'} = $accession;
    $cross_info{'parent_a'} = $parent_A;
    $cross_info{'parent_b'} = $parent_B;
    $cross_info{'progeny_pedigree'} = $unmasked_pedigree_str;
    $cross_info{'a_or_b'} = $a_or_b;
    $cross_info{'level'} = $level;
    print STDERR "Parent_A: $parent_A\n";
    print STDERR "Parent_B: $parent_B\n";
    my $cross_data_hash = $self->get_cross_data();
    $cross_data_hash->{$accession} = \%cross_info;
    $self->set_cross_data($cross_data_hash);
    return $accession;
}

sub _do_next_split_pedigree {
    my $self = shift;
    my $pedigree_str = shift;
    my $longest_start = shift;
    my $longest_end = shift;
    my $level = shift;
    my $parent_A = substr($pedigree_str, 0, $longest_start - 1);
    my $parent_B = substr($pedigree_str, $longest_end + 1, length($pedigree_str));
    my $parent_A_id;
    my $parent_B_id;
    $parent_A_id = _split_pedigree($self, $parent_A, "A", $level);
    $parent_B_id = _split_pedigree($self, $parent_B, "B", $level);
    if (!$parent_A_id || !$parent_B_id ) {
	print STDERR "Pedigree structure error: _do_split_pedigree\n";
	$self->set_parse_error('-1');
	return -1;
    }
    my @return_data = ($parent_A_id, $parent_B_id);
    return \@return_data;
}

sub _process_leaf {
    #currently doing nothing with this other than storing leaf names
    my $self = shift;
    my $pedigree_str = shift;
    my $accession = shift;
    my %leaf_info;
    #store values in hash
    $leaf_info{'accession'} = $accession;
    $leaf_info{'leaf'} = $pedigree_str;
    print STDERR "Leaf: $pedigree_str\n";
    my $leaf_data_hash = $self->get_leaf_data();
    $leaf_data_hash->{$accession} = \%leaf_info;
    $self->set_leaf_data($leaf_data_hash);
    return $pedigree_str;
}


#######
1;
#######
