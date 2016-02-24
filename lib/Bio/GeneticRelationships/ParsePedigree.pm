
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
use Data::Dumper;
use Bio::GeneticRelationships::Pedigree;
use Bio::GeneticRelationships::Individual;
use Data::GUID qw(guid);

has 'pedigree' => (isa => 'Str', is => 'ro', predicate => 'has_pedigree', required => 1);
has 'accession' => (isa => 'Str', is => 'ro', predicate => 'has_accession', required => 1);
has 'cross_data' => (isa => 'HashRef[HashRef[Str]]', is => 'rw', predicate => 'has_cross_data');
has 'parse_error' => (isa => 'Str', is => 'rw', predicate => 'has_parse_error');

sub parse_pedigrees {
    my $self = shift;
    my $pedigree_str = $self->get_pedigree();
    my $accession = $self->get_accession();

    my $parse_success = _split_pedigree($self, $pedigree_str, "A", 0, $pedigree_str, $accession);

    return 1;
}

sub is_root {
    my $self = shift;
    my $pedigree_str = $self->get_pedigree();
    my $accession = $self->get_accession();
    my $parse_success = _split_pedigree($self, $pedigree_str, "A", 0, $pedigree_str, $accession);
    if (!$parse_success) {
	return;
    }
    if ($parse_success = $pedigree_str) {
	return 1;
    } else {
	return;
    }
}


sub _split_pedigree {
    my $self = shift;
    my $pedigree_str = shift;
    my $a_or_b = shift;
    my $level = shift;
    my $unmasked_pedigree_str = shift;
    my $accession = shift;
    my $success = -1;
    my $parent_A;
    my $parent_B;

    my $longest = 0;
    my $next_longest = 0;
    my $longest_start;
    my $longest_end;

    print STDERR "\n_split_pedigree begin\n";
    print STDERR "Pedigree unprocessed: $pedigree_str\n";

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
	print STDERR "Error:mismatched parenthesis\n";
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
	    print STDERR "Pedigree begins and ends with ()$pedigree_str_trim\n";
	    
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
		print STDERR "Leading and trailing parenthesis removed\n";
		$pedigree_str = $pedigree_str_trim;
		$unmasked_pedigree_str = $pedigree_str_trim;
		@chars = @chars_trim;
		print STDERR "Trimmed pedigree is: $pedigree_str\n";
	    } else {
		$trim_done = 1;
	    }
	} else {
	    $trim_done = 1;
	}
    }

    if (substr($pedigree_str,0,1) eq '(') {

	# find closing parenthesis
	my $split_pos;
	my $parenth_balance = 1;
	my $str_pos = 0;
	while ($parenth_balance != 0) {
	    $str_pos++;
	    my $char = $chars[$str_pos];
	    if ($char eq '(') {
		$parenth_balance++;
	    }
	    elsif ($char eq ')') {
		$parenth_balance--;
	    }
	    if ($parenth_balance == 0) {
		$split_pos = $str_pos;
	    }
	}
	
	if (!$split_pos) {
	    print STDERR "Error:mismatched parenthesis\n";
	    $self->set_parse_error('-1');
	    return -1;
	}

	# split at closing parenthesis if at a /
	if (scalar(@chars) <= $split_pos+1) {
	    print STDERR "Error: parenthesis process error\n";
	    $self->set_parse_error('-1');
	    return -1;
	}

	
	if ($chars[$split_pos+1] eq '/') {
	    print STDERR "() followed by a / so splitting there\n";
	    print STDERR "sending to _do_next_split \nPedigree:$pedigree_str\nUnmasked:$unmasked_pedigree_str\n";
	    #_do_next_split_pedigree($self, $pedigree_str, $split_pos+2, $split_pos+1, $unmasked_pedigree_str, $level+1);
	    my $rightofslash = substr($pedigree_str, $split_pos+1, length($pedigree_str));
	    print STDERR "rightofslash: $rightofslash\n";
	    
	    my $crosschars =  $rightofslash =~ m/\/+/;

	    print STDERR "Crosschars: $crosschars\n";
	    
	    $longest = length($crosschars);
	    $longest_start = $split_pos + 2; 
	    $longest_end = $split_pos + length($crosschars);
	    
	}   
    }
    

    # discover / outside of any ()
    my $free_slash;
    my $parenth_balanth = 0;
    for my $char (@chars) {
    	if ($char eq '(') {
    	    $parenth_balanth++;
    	}
    	if ($char eq ')') {
    	    $parenth_balanth--;
    	}
    	if ($char eq '/' && $parenth_balanth == 0) {
    	    $free_slash = 1;
    	}
    }

    if ($free_slash) {
	print STDERR "Contains free slash\n";
    }

    if (!$longest_start) {
	if (!$free_slash) {
	    # find parenthesis used to name a cross of the form (parentA/parentB)crossname
	    my $mask_pos = 0;
	    my $parenth_open = 0;
	    my $parenth_close_pos;
	    my $parenth_open_pos;
	    for my $char (@chars) {
		if ($char eq "(") {
		    $parenth_open++;
		    if (!$parenth_open_pos) {
			$parenth_open_pos = $mask_pos;
		    }
		} elsif ($char eq ")") {
		    $parenth_open--;
		    #find closure of first open
		    if (($parenth_open == 0) && (!$parenth_close_pos)) {
			$parenth_close_pos = $mask_pos;
		    }
		}
		$mask_pos++;
	    }
	    # when cross name is at the beginning
	    if ($parenth_close_pos &&  (substr($pedigree_str,0,1) eq '(')) {
		my $after_close = substr($pedigree_str, $parenth_close_pos+1, length($pedigree_str));
		print STDERR "String following closing ) is:$after_close\n";
		#if no more crosses after close, this is a crossname
		if (!($after_close =~ m/\//)) { 
		    # if no accession was provided, use this cross name as the accession
		    if (!$accession) {
			$accession = $after_close;
			print STDERR "Set accession name to: $after_close\n";
		    }
		    print STDERR "The string follwing closing ) appears to be a cross name: $after_close\n";
		    $pedigree_str = substr($pedigree_str, 1, $parenth_close_pos - 1);
		    @chars = split('', $pedigree_str);
		    print STDERR "The new pedigree with the cross name removed is:$pedigree_str\n";
		}
	    }
	    # when cross name is at the end
	    # problem is here.  need don't find cross names when there is a '/' outside ()
	    if ($parenth_close_pos &&  (substr($pedigree_str,length($pedigree_str)-1,length($pedigree_str)) eq ')')) {
		my $before_close = substr($pedigree_str, 0, $parenth_open_pos);
		print STDERR "String before opening ( is:$before_close\n";
		if (!($before_close =~ m/\//)) { 
		    # if no accession was provided, use this cross name as the accession
		    if (!$accession) {
			$accession = $before_close;
			print STDERR "Set accession name to: $before_close\n";
		    }
		    print STDERR "The string before the opening ( appears to be a cross name: $before_close\n";
		    $pedigree_str = substr($pedigree_str, $parenth_open_pos+1, (length($pedigree_str)-length($before_close))-2);
		    @chars = split('', $pedigree_str);
		    print STDERR "The new pedigree with the cross name removed is:$pedigree_str\n";
		}
	    }
	    
	    # mask between parenthesis
	    $unmasked_pedigree_str = $pedigree_str;
	    $parenth_open = 0;
	    $mask_pos = 0;
	    for my $char (@chars) {
		if ($char eq "(") {
		    $chars[$mask_pos] = 0;
		    $parenth_open++;
		} elsif ($char eq ")") {
		    $chars[$mask_pos] = 0;
		    $parenth_open--;
		} elsif ($parenth_open > 0) {
		    $chars[$mask_pos] = 0;
		}
		$mask_pos++;
	    }
	    my $maskedseq = join('',@chars);
	    print STDERR "Masked between ():$maskedseq\n";

	}


	#find longest stretch of slashes
	#and check that there isn't another stretch of slashes of equal length

	my $current = 0;
	my $current_pos = 0;
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
	    print STDERR "Pedigree structure error:  split at equal levels\n";
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
		print STDERR "Slash detection error\n";
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
    }
    print STDERR "deepest cross is level: $longest\n";

    # return string if nothing left to split
    if ($longest == 0) {

	my $leaf_return;
	$leaf_return=_process_leaf($pedigree_str, $level);

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
    print STDERR "Right side is: $right_side\n";
    print STDERR "Left side is: $left_side\n";
    # check if there are any slashes and if it ends with a backcross
    my $backcross_number;
    my $backcross_accession;
    if (!($right_side =~ m/\//) && ($right_side =~ m/\*\d{1,2}$/)) {
	print STDERR "Right side does not contain slash and is backcross\n";
	my @split_backcross = split('\*', $right_side);
	$backcross_accession = $split_backcross[0];
	$backcross_number = $split_backcross[1];
	print STDERR "Backcross name: $backcross_accession\nBackcross number: $backcross_number\n";
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
	$pedigree_str = $unmasked_pedigree_str;

	    
	print STDERR "Backcross modified pedigree: $new_backcross_name\n";

    }


    my $parent_ids;

    print STDERR "/ based splitter sending to _do_next_split \nPedigree:$pedigree_str\nUnmasked:$unmasked_pedigree_str\n";
    $parent_ids = _do_next_split_pedigree($self, $pedigree_str, $longest_start, $longest_end, $unmasked_pedigree_str, $level+1);
    if (!$parent_ids || $parent_ids == -1) {
	print STDERR "\n\nPedigree structure error: parent ids\n";
	$self->set_parse_error('-1');
	return -1;
    }

    $parent_A = @$parent_ids[0];
    $parent_B = @$parent_ids[1];

 
    if (!$accession) {
	print STDERR "No accession defined, so creating one\n";

	# a unique code to use as an accession name if one isn't provided
	my $guid_code = Data::GUID->new()->as_hex();

	$accession = $pedigree_str."_".$guid_code;
    }


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
    $success = 1;


    return $accession;


}

sub _do_next_split_pedigree {
    my $self = shift;
    my $pedigree_str = shift;
    my $longest_start = shift;
    my $longest_end = shift;
    #my $cross_name = shift;
    #my $accession = shift;
    my $unmasked_pedigree_str = shift;
    my $level = shift;


    my $parent_A = substr($pedigree_str, 0, $longest_start - 1);
    my $parent_B = substr($pedigree_str, $longest_end + 1, length($pedigree_str));

    my $parent_A_id;
    my $parent_B_id;

    $parent_A_id = _split_pedigree($self, $parent_A, "A", $level, $parent_A);
    $parent_B_id = _split_pedigree($self, $parent_B, "B", $level, $parent_B);

 
    if (!$parent_A_id || !$parent_B_id ) {
	print STDERR "Pedigree structure error: _do_split_pedigree\n";
	$self->set_parse_error('-1');
	return -1;
    }
    my @return_data = ($parent_A_id, $parent_B_id);

    return \@return_data;
}

sub _process_leaf {
    my $pedigree_str = shift;
    my $level = shift;
    return $pedigree_str;
}


#######
1;
#######
