use strict;

my %onetothree = (A => "ALA", R => "ARG", N => "ASN", D => "ASP", C => "CYS", Q => "GLN",
		  E => "GLU", G => "GLY", H => "HIS", I => "ILE", L => "LEU", K => "LYS",
		  M => "MET", F => "PHE", P => "PRO", S => "SER", T => "THR", W => "TRP",
		  Y => "TYR", V => "VAL", B => "ASX", Z => "GLX", X => "XAA");
my @keys = keys %onetothree;
my @values = values %onetothree;
my %threetoone = ();
foreach my $i (0..$#keys)
{
    $threetoone{$values[$i]} = $keys[$i];
}

my %nametothree = (ALANINE => "ALA", ARGININE => "ARG", ASPARAGINE => "ASN", "ASPARTIC ACID" => "ASP",
		   CYSTEINE => "CYS", GLUTAMINE => "GLN", "GLUTAMIC ACID" => "GLU", GLYCINE => "GLY", 
		   HISTIDINE => "HIS", ISOLEUCINE => "ILE", LEUCINE => "LEU", LYSINE => "LYS",
		   METHIONINE => "MET", PHENYLALANINE => "PHE", PROLINE => "PRO", SERINE => "SER", 
		   THREONINE => "THR", TRYPTOPHAN => "TRP", TYROSINE => "TYR", VALINE => "VAL", STOP => "XAA");

my %codons = (TTT => "F", TTC => "F", TTA => "L", TTG => "L", TCT => "S",
	      TCC => "S", TCA => "S", TCG => "S", TAT => "Y", TAC => "Y",
	      TAA => "X", TAG => "X", TGT => "C", TGC => "C", TGA => "X",
	      TGG => "W", CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	      CCT => "P", CCC => "P", CCA => "P", CCG => "P", CAT => "H",
	      CAC => "H", CAA => "Q", CAG => "Q", CGT => "R", CGC => "R",
	      CGA => "R", CGG => "R", ATT => "I", ATC => "I", ATA => "I",
	      ATG => "M", ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	      AAT => "N", AAC => "N", AAA => "K", AAG => "K", AGT => "S",
	      AGC => "S", AGA => "R", AGG => "R", GTT => "V", GTC => "V",
	      GTA => "V", GTG => "V", GCT => "A", GCC => "A", GCA => "A",
	      GCG => "A", GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	      GGT => "G", GGC => "G", GGA => "G", GGG => "G");

sub iscodon
{
    my $aa = shift;
    $aa = uc $aa;
    return $codons{$aa};
}

sub threetoone
{
    my $aa = shift;
    $aa = uc $aa;
    my $amino = iscodon($aa);
    return $amino if (defined($amino));
    return $threetoone{$aa};
}

sub onetothree
{
    my $aa = shift;
    $aa = uc $aa;
    return $onetothree{$aa};
}

sub nametothree
{
    my $aa = shift;
    $aa = uc $aa;
    return $nametothree{$aa};
}


sub mutation_conversion
{
	my @muts = @_;
	my @results;
	my $protein = 0;
	foreach my $aa (@muts)
	{
		unless (($aa =~ /[ATGC]{1}/i) && ($aa =~ /\w{1}/))
		{
			$protein = 1;
			last;
		}
	}
	
	foreach my $aa (@muts)
	{
		my $mut = $aa;
		my $changed;
		if ($aa =~ /\w{4,}/)
		{
			$changed = nametothree($aa);
		}
		elsif ($aa =~ /\w{3}/)
		{
			my $amino = iscodon($aa);
			if (defined($amino))
			{
				$changed = onetothree($amino);
			}
			else
			{
				$changed = $aa; 
			}
		}
		else
		{
			if ($protein == 1)
			{
				$changed = onetothree($aa);
			}
			else
			{
				$changed = $aa;
			}
		}
		$changed = $mut unless(defined($changed));
		push(@results, uc $changed);
	}
	return @results;
}

=pod Subroutines

=head1 iscodon

This subroutine takes one parameter and checks if it is a 
3 nucleic acid codon. If it is the one letter amino acid is 
returned, else nothing is returned.

=head1 onetothree

This subroutine takes a single letter amino acid 
(case-insensitive) and returns the three letter amino acid.
If the input cannot be converted nothing is returned. 

=head1 threetoone

This subroutine takes a three letter amino acid abbreviation
and returns a one letter abbreviation. The subroutine also 
uses iscodon() in case the 3 letter string is a codon. If it
is a codon the appropriate amino acid for humans is returned.
The subroutine is case-insensitive. If the input can't be
recognized nothing is returned.

=head1 nametothree

This subroutine takes an amino acid name (case-insensitive)
and returns the three letter abbreviation. If it fails
nothing is returned. 

=head1 mutation_conversion

This is a very specific subroutine. It is intended to take
two parameters (wtaa and mtaa) and return and array of the 
three letter abbreviations unless both the wtaa and mtaa are
in the character class [ATGC]. However, it can actually take 
any size array of suspected amino acids and unless all of 
them are in the aforementioned character class the 3 letter 
abbreviations will be returned. Because it is crucial that
the number of returned elements are the same in number and 
order any bad data (i.e. stop, ins, del, etc.) will simply
be returned in the appropriate postion. This subroutine can 
take ANYTHING (full name, three letter abbreviation, one 
letter abbreviation) and give back a result.

=head2 EXAMPLE

my ($wtaa, $mtaa) = mutation_conversion($wt, $mt);

=cut


1;
