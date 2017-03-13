use strict;
use Cwd;
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin); 

my $root = $Bin; 
require $root."/AACONVERSION.pm";

my $print_without_positions = 0;
my $DEBUG                   = 0;

my $ver    = "1.19";
#my $infile = $ARGV[0] 
  ; #list of abstracts. Each line must contain tab delimited pubmed ID \t Title \t Abstract

use vars qw( $opt_h $opt_s $opt_f $header $run_sent $infile);

## Check command line
my $prog = basename($0);

getopt('fhs');

$opt_h = uc($opt_h);
$opt_s = uc($opt_s);

if($opt_f){
   $infile = $opt_f;
}

if ($opt_h eq "YES") {
  $header = "TRUE";
} 

if ($opt_s eq "YES") {
  $run_sent = "YES";
}

#my $run_sent = $ARGV[1];
	#takes the given text and runs EMU on each sentence.  flag: -y for sentence run or -n for using full input text
my $HugoGeneNames = $root."/HUGOGeneNames.txt"
  ; #file that contains the list of the HUGO gene names. Each row must contain one gene name.
my $CellLineNames = $root."/Cell_line_list_short.txt"
  ; #file that contains the list of the Cell Line names. Each row must contain one gene name.
#my $AbGeneNames;   #"AbGene_PCa.txt.genes"; #list of gene names instead of HUGO gene names
my $OrganismFile; #"pc_mut-ct_humans_animals.txt"; #organism to the abstract
my $argc = @ARGV;
=pod
if ($argc > 1){
  $AbGeneNames = $ARGV[1];    #file that contains the list of the AbGene names. Each row must contain one gene name.
}
=cut
if ($argc > 2){
  $OrganismFile = $ARGV[2];   # file that contains the organism information
}


# Variables for storing pattern-matching components
my ( $left_border, $right_border, $symbol_connective, $symbol_connective_2,
    $rs_connective, $word_connective, $na_pre_match_set_char, $indel_connective );
my ( $residue_set1, $residue_set2, $residue_set3, $residue_set4, $residue_set5,
     $residue_set6, $residue_set7 )
  ;    # Amino Acid subject & object declarations

# Patterns
my ( $wtaa, $mtaa, $pos1, $pos2, $pos3, $pos4, $pos, $position, $rs_mutation, $pattern_body,
    $position_patterns, $position_patterns2, $pattern, $multi, $conn, $mutation_pos, $pos_patt,
	$variant_type );

my ( $mutations, $genes, $patt, $pos, $tmp1, $tmp2, $type );
my ( $position23, $position1, $position4 );
my $cell_lines;

my %normalize_indel = ();
# Variables for storing false-positive patterns
#my ($fallible_patterns);

#initializes the mutation patterns
load_positive_patterns();

#Reading additional informations
my %organisms = (); #hash
my $organism;
open( FILE, "$OrganismFile" );
while (<FILE>) {
    if (/^\d+/)
    {    # only need to read lines that start with numbers not empty lines
        chomp;
        my ( $pmid, $organism ) = split /\|/;
        $organisms{$pmid} .= "$organism;";
    }
}

close FILE;

my @Hugo_gene_name = ();
open( FILE, "$HugoGeneNames" ) || die "can't open file Hugo: $!\n\n";
while (<FILE>) {
    chomp;
    my ( $gene, $other ) = split /\~/;
    push @Hugo_gene_name, $gene if ( length($gene) > 1 && !($gene =~ m/($residue_set4)/g)); #skips gene of length less than two and if gene name can be confused with codon 
}
#add the gene names p53 manually to the gene name list.
push @Hugo_gene_name, "p53";
push @Hugo_gene_name, "P53";

close FILE;


#skips cell lines in abstract
my @Cell_line_name = ();
open( FILE, "$CellLineNames" ) || die "can't open file Cell Lines: $!\n\n";
while (<FILE>) {
    chomp;
    my ( $gene, $other ) = split /\~/;
    chop ($gene) if ($gene =~ /\s$/);
    push @Cell_line_name, $gene;
}
close FILE;

=pod
#opens Ab gene names
my %Ab_gene_name = ();
open( FILE, "$AbGeneNames" );
#if (open( FILE, "$AbGeneNames" )){
  my $patt;
  my $pmid;
  while (<FILE>) {

      chomp;
      $patt = "-----\\D+(\\d+)\\D+-----";

      if (/$patt/g) {
          $pmid = $1;
      }
      else {
          $Ab_gene_name{$pmid} .= "$_;";
      }
  }
#}
close FILE;
=cut
print "FINISH READING ADDITIONAL DATA\n";

#load_fallible_patterns();

#Open output files
my $file1 = "EMU_$ver" . "_HUGO_$infile";
#my $file2 = "EMU_$ver" . "_ABG_$infile";
#my $file3 = "EMU_$ver" . "_fall_$infile";

open( OUTPUT_F,  ">$file1" ) || die "can't open file $file1: $!\n\n";
#open( OUTPUT_F2, ">$file2" ) || die "can't open file $file2: $!\n\n";
#open( OUTPUT_Fa, ">$file3" ) || die "can't open file $file3: $!\n\n";
printf OUTPUT_F "pmid\torganism\tmut_pat1\tpos_patt\tmutation type\twtaa\tmtaa\tpos\tgenes\ttype\n";
#printf OUTPUT_F2 "pmid\torganism\tmut_pat1\tpos_patt\twtaa\tmtaa\tpos\tgenes\ttype\n";

#Reading input abstracts
if($run_sent eq "YES"){
	my $input = mk_sent($infile);
	open( FILE, "$input" ) || die "can't open file $infile: $!\n\n";
}else{
	open( FILE, "$infile" ) || die "can't open file $infile: $!\n\n";
}
 
 my $tmp = <FILE> if ($header eq "TRUE");

while (<FILE>) {

    chomp;
    my ( $pmid, $title, $abstract ) = split /\t/;

    my $text            = " $title $abstract";
    my %mutations_table = ();
    my %mutation_type   = ();
	my %mutation_variant = ();
    $genes = "";
    $cell_lines = "";
    
    foreach my $gene (@Hugo_gene_name) {

        if ( $text =~ m/\W$gene\W/ ) {    # find all genes matchs per article
            $genes .= "$gene;";           #
        }
    }
    chop($genes);
	
	#for sent, if no genes, next sentence
	next if ($genes eq "" && $run_sent eq "-y");
	
	
    foreach my $gene (@Cell_line_name) {
        if ( $text =~ m/\W$gene\W/ ) {    # find all cell lines matchs per article
            $cell_lines .= "$gene;";           #
        }
    }
    my %mut_store;

    $patt  = "$left_border"."(p\\.$residue_set2(\\d+)$residue_set2)"."$right_border";
    while ( $text =~ s/($patt)/ / ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $pos      = $4;
        ( $wtaa, $mtaa ) = mutation_conversion( $3, $5 );

        print("First1 a Print; $1;$2;$3;$4;$5;$6;\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }
        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt";
        $mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN";
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
        print("Second1 Print: $wtaa;$pos;$mtaa\n") if ($DEBUG);
    }
    $patt  = "$left_border"."(p\\.($residue_set1)(\\d+)($residue_set1))"."$right_border";
    while ( $text =~ s/($patt)/ /) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $pos      = $4;
        ( $wtaa, $mtaa ) = mutation_conversion( $3, $5 );

        print("First1 b Print; $1;$2;$3;$4;$5;$6;\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }
        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt";
        $mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN" if (length($wtaa) > 1);
		$mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN;DNA" if (length($wtaa) == 1);
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
        print("Second1 Print: $wtaa;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt  = "$left_border"."(([cgmr])\\.([+-\\d]+)($residue_set1)($symbol_connective)($residue_set1))"."$right_border";
    while ( $text =~ s/($patt)/ / ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $type     = $3;
        $pos      = $4;
        ( $wtaa, $mtaa ) = mutation_conversion( $5, $7 );

        print("First2 c Print; $1;$2;$3;$4;$5;$6;\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;
        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }
        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt";
        $mutation_type{"$wtaa;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$wtaa;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$wtaa;$pos;$mtaa"} = "MITHO" if ($type eq "m");
        $mutation_type{"$wtaa;$pos;$mtaa"} = "RNA" if ($type eq "r");
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
        print("Second2 Print: $wtaa;$pos;$mtaa\n") if ($DEBUG);
    }


    $patt =    # analyzing pattern for one letter code
      "$left_border"
      . "(($position1)"
      . "($residue_set1)"
      . "($position23)"
      . "($symbol_connective)"
      . "($position23)"
      . "($residue_set1)"
      . "($position4))"
      . "$right_border";
    while ( $text =~ m/(?=($patt))/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $conn     = $9;
        $pos1     = $4;
        $tmp1     = $5;
        $pos2     = $8;
        $pos3     = $10;
        $tmp2     = $14;
        $pos4     = $15;
        ( $wtaa, $mtaa ) = mutation_conversion( $6, $12 );

        if ($tmp1 =~ m/[\.,;%]/  || $tmp2 =~ m/[\.,;%]/) {
          next;
        }

        print(
            "First Print d; $1;$2;$3;$4;$5;$6;$7;$8;$9;$10;$11;$12;$13;$14;$15;$16;$17\n")
          if ($DEBUG);
        ( $pos, $multi ) = one_position( $pos1, $pos2, $pos3, $pos4 );
        unless ( $pos2 || $pos3 || $conn ) { $pos = ""; }

        $mutation_pos =
          length($`);    #stores the position where the mut pattern was found

        if ( !$pos && $conn && !$multi )
        {    #find positions when there wasn't a position in the pattern

            ( $pos, $pos_patt ) =
              search_closest_position( $mutation_pos, $text );
        }
        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;
        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

		#adding types
		$mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN" if (length($wtaa) > 1);
        $mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN;DNA" if (length($wtaa) == 1);
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
		#end addition

        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt"
          if ( ($mutation_type{"$wtaa;$pos;$mtaa"}) && ( $print_without_positions || $pos ) && ( $pos || $conn ) );
        print("Second Print: $wtaa;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
      . "(($position1)"
      . "($residue_set2)"
      . "($position23)"
      . "($symbol_connective|$word_connective)"
      . "($position23)"
      . "($residue_set2)"
      . "($position4))"
      . "$right_border";
    while ( $text =~ m/(?=($patt))/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $pos1     = $4;
        $tmp1     = $5;
        $pos2     = $9;
        $pos3     = $11;
        $tmp2     = $16;
        $pos4     = $17;
        ( $wtaa, $mtaa ) = mutation_conversion( $6, $13 );

        if ($tmp1 =~ m/[\.,;%]/){
          $pos1 = "";
        } elsif ($tmp2 =~ m/[\.,;%]/) {
          $pos4 = "";
        }

        print(
            "Third print e: $1;$2;$3;$4;$5;$6;$7;$8;$9;$10;$11;$12;$13;$14;$15;$16;$17\n")
          if ($DEBUG);
        ( $pos, $multi ) = one_position( $pos1, $pos2, $pos3, $pos4 );

        $mutation_pos = length($`);

        if ( !$pos && !$multi ) {
            ( $pos, $pos_patt ) =
              search_closest_position( $mutation_pos, $text );
        }
        print("Printing into table hash: $wtaa;$pos;$mtaa\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;
        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/"($pattern_body)"/);
        }
		$mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN";
        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt"
          if ( ($mutation_type{"$wtaa;$pos;$mtaa"}) && $print_without_positions || $pos );
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
        print("Fourth print: $wtaa;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
      . "(($residue_set4)"
      . "($position23)"
      . "($symbol_connective)"
      . "($position23)"
      . "($residue_set4))"
      . "$right_border";
    while ( $text =~ m/(?=($patt))/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $pos1     = $6;
        $conn     = $7;
        $pos2     = $9;
        ( $wtaa, $mtaa ) = mutation_conversion( $3, $10 );

        print("f $1;$2;$3;$4;$5;$6;$7;$8;$9;$10;$11;$12;$13;$14;$15\n")
          if ($DEBUG);

        ( $pos, $multi ) = one_position( $pos1, $pos2 );

        $mutation_pos = length($`);

        if ( !$pos && $conn && !$multi ) {
            ( $pos, $pos_patt ) =
              search_closest_position( $mutation_pos, $text );
        }
        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;
        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }
		$mutation_type{"$wtaa;$pos;$mtaa"} = "PROTEIN";
        $mutations_table{"$wtaa;$pos;$mtaa"} = "$pattern\t$pos_patt"
          if ( ( $print_without_positions || $pos ) && ( $pos || $conn ) );
		$mutation_variant{"$wtaa;$pos;$mtaa"} = "MISSENSE";
        print("$wtaa;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)" . "$left_border" . "($rs_connective)" . "$right_border";
    while ( $text =~ m/(?=($patt))/g ) {    # find all pattern matchs per article

        $pattern     = $1;
        $rs_mutation = $2;
        print("$rs_mutation\n") if ($DEBUG);
		$mutation_type{"$rs_mutation;;"} = "PROTEIN;DNA;RNA";
		$mutation_variant{"$rs_mutation;;"} = "RSID";
        $mutations_table{"$rs_mutation;;"} = "$pattern\t";
    }

	#############################
	#INS and DEL patterns
	#############################
	$patt = "(?i)"
      . "$left_border"
	  . "((p|c|g|m|r)?\\.?"
      . "([\\+\\-_\\d]+)"
      . "\s?($indel_connective)"
      . "($residue_set5))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;
        $pos      = $4;
		$variant_type = $5;
        $mtaa = $6;

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

        $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($pos < 0);
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN" if ($type eq "p" || ( $mtaa =~ /[RNDQEHILKMFPSWYV]/ && $mtaa !~ /[ATUGC]/ ));
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "RNA" if ($type eq "r");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA" if ($pos > 0 && $type eq "" && $mtaa =~ /T/i && $mtaa !~ /U/i);
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($pos > 0 && $type eq "" && $mtaa !~ /[UT]/i);
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }

	$patt = "(?i)"
      . "$left_border"
	  . "((p)\\."
      . "([\\+\\-_\\d]+)"
      . "\s?($indel_connective)"
      . "($residue_set6))"
      . "$right_border";
    while ( $text =~ m/(?=($patt))/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;
        $pos      = $4;
		$variant_type = $5;
        $mtaa = $6;

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN";
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
	  . "(([gcrmp]?)\\.?"
      . "([+\\-_\\d]+)"
      . "\s?($indel_connective)"
      . "(\\d+))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;        
		$pos      = $4;
		$variant_type = $5;
        $mtaa = $6;

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

        $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN" if ($type eq "p");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "RNA" if ($type eq "r");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($type eq "");
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
	  . "(([gcrmp]?)\\.?"
      . "([\\+\\-_\\d]+)"
      . "($indel_connective)"
      . "(\\d+)"
      . "([\\-_])"
      . "($indel_connective)"
      . "(\\d+))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;
        $pos      = $4;
		$variant_type = "$5|$8";
        $mtaa = "$6|$9";

		$pos =~ s/_/-/g;
		#my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8;$9\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

		$mutation_type{"$variant_type;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$variant_type;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$variant_type;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$variant_type;$pos;$mtaa"} = "PROTEIN" if ($type eq "p");
        $mutation_type{"$variant_type;$pos;$mtaa"} = "RNA" if ($type eq "r");
        $mutation_type{"$variant_type;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($type eq "");
        $mutations_table{"$variant_type;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$variant_type;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $variant_type;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
      . "(([+\\-_\\d]+)"
      . "(\\s+)"
      . "(\\d+N)"
      . "(\\s+)"
      . "($indel_connective))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
        $pos      = $3;
		$variant_type = "$7";
        $mtaa = "$5";

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8;$9\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

		$mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA";
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }

    $patt = "(?i)"
      . "$left_border"
	  . "(([gcrmp]?)\\.?"
      . "\\s*([\\+\\-_\\d]+)"
      . "(\\s+)"
      . "($indel_connective)"
      . "(\\s+)"
      . "($residue_set6))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;
        $pos      = $4;
		$variant_type = "$6";
        $mtaa = "$8";

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8;$9\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

		$mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN" if ($type eq "p" || ( $mtaa =~ /[RNDQEHILKMFPSWYV]/ && $mtaa !~ /[ATUGC]/ ));
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "RNA" if ($type eq "r");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA" if ($type eq "" && $mtaa =~ /[T]/i);
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($type eq "" && $mtaa !~ /[RNDQEHILKMFPSTWYVUT]/i && $mtaa =~ /[ACG]/i);
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }

    #adding for pattern R541-E543delinsK
	$patt = "(?i)"
      . "$left_border"
	  . "((p|c|g|m|r)?\\.?"
      . "($residue_set6)"
      . "([\\+\\-_\\d]+)"
      . "($residue_set6)"
      . "([\\+\\-_\\d]+)"
      . "($indel_connective)"
      . "($residue_set6))"
      . "$right_border";
    while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
        $pos_patt = "";
        $pattern  = $1;
        $pattern_body = $2;
		$type = $3;
        $pos      = "$5-$7";
		$variant_type = $8;
        $mtaa = $9;

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

        print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8\n") if ($DEBUG);

        $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
        $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

        if ($cell_lines ne "") {
            print("$pattern_body \t $cell_lines\n") if ($DEBUG);
            next if ($cell_lines =~ m/($pattern_body)/);
        }

        $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($type eq "c");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "GENOM" if ($type eq "g");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN" if ($type eq "p" || ( $mtaa =~ /[RNDQEHILKMFPSWYV]/ && $mtaa !~ /[ATUGC]/ ));
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "RNA" if ($type eq "r");
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA" if ($pos > 0 && $type eq "" && $mtaa =~ /T/i && $mtaa !~ /U/i);
        $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($pos > 0 && $type eq "" && $mtaa !~ /[UT]/i && $mtaa !~ /[RNDQEHILKMFPSTWYVUT]/i );
        $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
        print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
    }


	$patt = "(?i)"
      . "$left_border"
	  . "((p|c|g|m|r)?\\.?"
      . "IVS([\\+\\-_\\d]+)"
      . "($indel_connective)"
      . "($residue_set5|\\d+))"
      . "$right_border";
	 
	while ( $text =~ m/($patt)/g ) {   # find all pattern matchs per article
	    $pos_patt = "";
	    $pattern  = $1;
	    $pattern_body = $2;
		$type = $3;
	    $pos      = $4;
		$variant_type = $5;
	    $mtaa = $6;

		$pos =~ s/_/-/g;
		my $norm_variant = $normalize_indel{$variant_type};

	    print("First1 Print: $1;$2;$3;$4;$5;$6;$7;$8\n") if ($DEBUG);

	    $pattern_body =~ s/\(/\\(/g; $pattern_body =~ s/\)/\\)/g;
	    $pattern_body =~ s/\[/\\[/g; $pattern_body =~ s/\]/\\]/g;

	    if ($cell_lines ne "") {
	        print("$pattern_body \t $cell_lines\n") if ($DEBUG);
	        next if ($cell_lines =~ m/($pattern_body)/);
	    }

	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($pos < 0);
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "DNA" if ($type eq "c");
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "GENOM" if ($type eq "g");
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "MITO" if ($type eq "m");
		$mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN" if ($type eq "p" || ( $mtaa =~ /[RNDQEHILKMFPSWYV]/ && $mtaa !~ /[ATUGC]/ ));
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "RNA" if ($type eq "r");
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA" if ($pos > 0 && $type eq "" && $mtaa =~ /T/i && $mtaa !~ /U/i);
	    $mutation_type{"$norm_variant;$pos;$mtaa"} = "PROTEIN;DNA;RNA" if ($pos > 0 && $type eq "" && $mtaa !~ /[UT]/i);
	    $mutations_table{"$norm_variant;$pos;$mtaa"} = "$pattern\t$pos_patt";
		$mutation_variant{"$norm_variant;$pos;$mtaa"} = "INDEL";
	    print("Second1 Print: $norm_variant;$pos;$mtaa\n") if ($DEBUG);
	}

    chop( $organisms{$pmid} );
    #chop( $Ab_gene_name{$pmid} );
    foreach my $mutations ( keys %mutations_table ) {
        my ( $wtaa, $pos, $mtaa ) = split( /;/, $mutations );
		print "found $mutations\n" if ($mutations eq "G;1141;T");
        next if ($wtaa eq $mtaa);
        print OUTPUT_F
"$pmid\t$organisms{$pmid}\t$mutations_table{$mutations}\t$mutation_variant{$mutations}\t$wtaa\t$mtaa\t$pos\t$genes\t$mutation_type{$mutations}\n";
        #print OUTPUT_F2
#"$pmid\t$organisms{$pmid}\t$mutations_table{$mutations}\t$wtaa\t$mtaa\t$pos\t$Ab_gene_name{$pmid}\t$mutation_type{$mutations}\n";
    }
#     print OUTPUT_F "$pmid\t$organisms{$pmid}\t\t\t\t\t$genes\n" unless (%mutations_table);
}
close OUTPUT_F;
#close OUTPUT_F2;
print "FINISHED\n";

### LOADS PATTERNS INTO ARRAYS ###
sub load_positive_patterns {

    # Variations of possible left border characters
    $left_border =
"(?:(?:\\s+)|(?:,)|(?:\\.)|(?:-)|(?:\\/)|(?:;)|(?:=)|(?::)|(?:\\*)|(?:\\[)|(?:\\])|(?:\\{)|(?:\\})|(?:\\()|(?:\\)))";

    # Variations of possible right border characters
    $right_border =
#"(?:(?:\\s+)|(?:,)|(?:\\.)|(?:;)|(?::)|(?:\\/)|(?:\\*)|(?:\\[)|(?:\\])|(?:\\{)|(?:\\})|(?:\\()|(?:\\)))";
 "(?:(?:\\s+)|(?:,)|(?:-)|(?:\\.)|(?:;)|(?::)|(?:\\/)|(?:\\*)|(?:\\])|(?:\\})|(?:\\)))";

# Variations of possible connectives linking wild type and mutation
#	$symbol_connective = "(?|(?:\\s+to\\s+)|(?:-+to-+)|(?:\\s+(?:to|into|for|of|by)\\s+(?:a\\s+))|(?:\\s*-*(?:>)\\s*)|(?:\\s*\\d+\\s*-*\\s*)|(?:---+))";
    $symbol_connective =
"(?:|(?:\\s*)|(?:\\s+to\\s+)|(?:-+to-+)|(?:\\s+(?:to|into|for|of|by)\\s+(?:a\\s+))|(?:\\s*-*(?:>)\\s*)|(?:---+))";

    # Specific RS id number for known mutations
    $rs_connective = "(?:[r|s]s\\d+(?:[a-zA-Z])*)";

    $position_patterns =
"(?:amino acid|position|base pair|residue|nucleotide|codon)\\s*(\\d{2,})|((?:intron|exon)\\s*\\d{1,})";
    $position_patterns2 =
"(?:amino acid|position|base pair|residue|nucleotide|codon|)\\s*((?:intron|exon|)\\s*\\d{2,})";
#"(?:codon|amino acid|position|base pair|residue|nucleotide|intron|exon)\\s*(\\d{1,})";

    # Templates for positions.
    $position23 = "[-\\(\\{\\[\\s]*(\\d{2,})[\\)\\}\\]\\-\\s]*|";

#    $position1 = "\\W*(\\d{2,})(?:\\W*)|";
    $position1 = "\\W*$position_patterns2(\\W*)|";
#    $position1 = "\\W*(\\d{2,})(?:[[\\W]&&[^\\.,;]]*|\\-*|\\s+)|";
#    $position4 = "(?:[[\\W]&&[^\\.,;]]*|\\-*|\\s+|[\\(\\{\\[\\-])(\\d{2,})\\W*|";
    $position4 = "(\\W*)(\\d{2,})\\W*|";


    # Amino Acid match set assignments
    $residue_set1 = "[ARNDCQEGHILKMFPSTWYV]";
    $residue_set2 =
"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Alanine|Arginine|Asparagine|Aspartic Acid|Cysteine|Glutamic Acid|Glutamine|Glycine|Histidine|Isoleucine|Leucine|Lysine|Methionine|Phenylalanine|Proline|Serine|Threonine|Tryptophan|Tyrosine|Valine)";

#"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)";
#    $residue_set3 =
#"(Alanine|Arginine|Asparagine|Aspartic Acid|Cysteine|Glutamic Acid|Glutamine|Glycine|Histidine|Isoleucine|Leucine|Lysine|Methionine|Phenylalanine|Proline|Serine|Threonine|Tryptophan|Tyrosine|Valine)";
    $residue_set4 = "([ACGT]{3})";

	$residue_set5 = "([ACGTU]{1,})";

	$residue_set6 = "[ARNDCQEGHILKMFPSTWYV]{1,}";

	$residue_set7 = "(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Alanine|Arginine|Asparagine|Aspartic Acid|Cysteine|Glutamic Acid|Glutamine|Glycine|Histidine|Isoleucine|Leucine|Lysine|Methionine|Phenylalanine|Proline|Serine|Threonine|Tryptophan|Tyrosine|Valine){1,}";

 # List of possible words and adverbs used for connecting wild type and mutation
    my @verb_List = (
        "substit\\w*", "lead\\w*",   "replac\\w*", "chang\\w*",
        "mutat\\w*",   "devia\\w*",  "modif\\w*",  "alter\\w*",
        "switch\\w*",  "exhang\\w*", "variat\\w*", "instead\\w*"
    );
    my @adverb_List = ( "to", "into", "for", "of", "by", "with" );
    my $verb_List   = join( "|", @verb_List );
    my $adverb_List = join( "|", @adverb_List );

# Variations of possible words and phrases linking wild type and mutation
#	$word_connective = "(?:\\s*\\d*\\s+)" . "(?:" . "(?:(?:\\w+\\s*){0,5}(?:$verb_List)\\s+)" . "|" . "(?:(?:\\w+\\s*){0,5}(?:$adverb_List)(?:\\s+a)?\\s+)" . "|" . "(?:(?:\\w+\\s*){0,5}(?:$verb_List)\\s+(?:$adverb_List)(?:\\s+a)?\\s+)" . "|" . "(?:\\s+(?:$verb_List)\\s+(?:$adverb_List)(?:\\s+a)?\\s+)" . ")";
    $word_connective =
        "(?:\\s+)" . "(?:"
      . "(?:(?:\\w+\\s*){0,5}(?:$verb_List)\\s+)" . "|"
      . "(?:(?:\\w+\\s*){0,5}(?:$adverb_List)(?:\\s+a)?\\s+)" . "|"
      . "(?:(?:\\w+\\s*){0,5}(?:$verb_List)\\s+(?:$adverb_List)(?:\\s+a)?\\s+)"
      . "|"
      . "(?:\\s+(?:$verb_List)\\s+(?:$adverb_List)(?:\\s+a)?\\s+)" . ")";

	$indel_connective = "ins\/del|del\/ins|del\/del|insdel|delins|ins|del|insertion|deletion";

	%normalize_indel = (insdel => "INDEL", delins => "INDEL", "ins\/del" => "INDEL", "del\/ins" => "INDEL", "ins-del" => "INDEL", "del-ins" => "INDEL", ins => "INS", del => "DEL", "del/del" => "DEL");
}

=pod
### Load fallible patterns Into Arrays

sub load_fallible_patterns {

    # Amino Acid false positives, Nucleic Acid false positives
    my (@fallible_patterns);

    # Case Insensitive Option
    push( @fallible_patterns, "(?i)" );

    #  Sequence of three or more Amino Acids (ie. Tye-Cys-Gly or T-C-G )
    push( @fallible_patterns,
"(?:$residue_set1\\s*(?:-|(?:>))\\s*$residue_set1\\s*(?:-|(?:>))\\s*$residue_set1)
|(?:$residue_set2\\s*(?:-|(?:>))\\s*$residue_set2\\s*(?:-|(?:>))\\s*$residue_set2)
|(?:$residue_set4\\s*(?:-|(?:>))\\s*$residue_set4\\s*(?:-|(?:>))\\s*$residue_set4)"
    );
#|(?:$residue_set3\\s*(?:-|(?:>))\\s*$residue_set3\\s*(?:-|(?:>))\\s*$residue_set3)

    #OR
    push( @fallible_patterns, "|" );

# Digit(s) follwed by Amino Acid  followed by digit(s) follwed by Amino Acid (ie. 45Cys56Arg or 45C56A)
    push( @fallible_patterns,
"(?:\\d+$residue_set1\\d+$residue_set1)
|(?:\\d+$residue_set2\\d+$residue_set2)
|(?:\\d+$residue_set4\\d+$residue_set4)"
    );
#|(?:\\d+$residue_set3\\d+$residue_set3)

    # OR
    push( @fallible_patterns, "|" );

# Amino Acid  followed by digit(s) follwed by Amino Acid followed by digit(s)(ie. Cys45Arg56 or C45A56)
    push( @fallible_patterns,
"(?:$residue_set1\\d+$residue_set1\\d+)
|(?:$residue_set2\\d+$residue_set2\\d+)"
#|(?:$residue_set3\\d+$residue_set3\\d+)"
    );

    # OR
    push( @fallible_patterns, "|" );

# Amino Acid followed by a zero followed by digit(s) followed by Amino Acid (ie. Thr035Trp or R0556C)
    push( @fallible_patterns,
"(?:(?:$residue_set1)(?:0\\d+)(?:$residue_set1))
|(?:(?:$residue_set2)(?:0\\d+)(?:$residue_set2))"
#|(?:(?:$residue_set3)(?:0\\d+)(?:$residue_set3))
    );

    # OR
    push( @fallible_patterns, "|" );

# Amino Acid (Three Letter Name) followed by digit(s) followed by hyphen followed by mino Acid (Three Letter Name) followed by digit
    push( @fallible_patterns,
        "(?:(?:$residue_set2)\\d+-(?:$residue_set2)\\d)" );

    # OR
    push( @fallible_patterns, "|" );

# Left border character followed by Amino Acid followed by digit followed by an optional hyphen followed by Amino Acid followed by right border character
    push( @fallible_patterns,
"(?:(?:$left_border)(?:$residue_set1)(?:\\d-?)(?:$residue_set1)(?:$right_border))"
    );

    # OR
    push( @fallible_patterns, "|" );

# Sequence of Amino Acid followed by digit(s) followed by symbol connective followed by Amino Acid followed by digit(s) (ie. Arg45-->Cys98 or A45-->C98)
    push( @fallible_patterns,
"(?:(?:$residue_set1)(?:\\d+)(?:$symbol_connective)(?:$residue_set1)(?:\\d+))
|(?:(?:$residue_set2)(?:\\d+)(?:$symbol_connective)(?:$residue_set2)(?:\\d+))"
    );

    push( @fallible_patterns, "|" );
    push( @fallible_patterns,
"(?:(?:$residue_set4\\s*(?:-|(?:>))\\s*$residue_set4\\s*(?:-|(?:>))\\s*$residue_set4))"
    );

    # Join Amino Acid false positives
    $fallible_patterns = join( '', @fallible_patterns );
}
=cut

sub search_closest_position {
    my $closest      = 50000;
    my $mutation_pos = shift;
    my $text         = shift;
    my ( $pos, $pos_patt );

    while ( $text =~ m/($position_patterns)/g )
    {    # find all pattern matchs per article

        if ( abs( length($`) - $mutation_pos ) < $closest ) {
            $closest  = abs( length($`) - $mutation_pos );
            if ($2) {
              $pos      = $2;
            } else {
              $pos      = $3;
            }
            $pos_patt = $1;
        }
    }
    return ( $pos, $pos_patt );
}

sub one_position {
    my $pos        = "";
    my $multi      = 0;
    my (@position) = @_;
    my $len        = scalar(@position);

    for ( my $i = 0 ; $i < $len ; $i++ )
    {    #removing empty positions from the array
        print "$i\t$position[$i]\n" if ($DEBUG);
        if ( $position[$i] eq "" ) {
            splice( @position, $i, 1 );
            $len--;
            $i--;
        }
    }
    print "LEN:$len\n" if ($DEBUG);
    if ( $len == 1 || ( $len == 2 && $position[0] == $position[1] ) ) {
        $pos = $position[0];
    }
    elsif ( $len > 1 ) {
        $multi = 1;
    }
    print "POS:$pos\tMULTI:$multi\n" if ($DEBUG);
    return ( $pos, $multi );
}

sub mk_sent {
	my $file = shift;
	my $output = "$file.out.txt" ;
	open (FILE , "<$file")|| die "could not open the input $file \n";
	open (OUT, ">$output") || die "could not open the sentence output file\n";
	my $tmp = <FILE> if ($header eq "TRUE");
	print OUT "pmid\ttitle\n";
	#my $temp_sent = "";
	while(<FILE>){
		chomp;
		my $row = $_;
		my($pmid, $title, $text) = split(/\t/, $row);
		#print "$row\n";
		print OUT "$pmid\ttitle\t$title\n";
		my $temp_sent = "";
		my @sentences = split(/\.\s/, $text);
		
		
		foreach my $sent (@sentences){
			if($sent =~ /^[a-zA-Z]\.\s$/ || $sent =~ /\s[a-zA-Z]\.\s$/ || $sent =~ /\s[a-zA-Z][a-zA-Z]\.\s$/ || $sent =~ /^[a-zA-Z][a-zA-Z]\.\s$/){
				#print "found $row\n";
				$temp_sent .= "." + $sent;
				next;
			}
			
			if($temp_sent eq ""){
				#print "in if\n";
				print OUT "$pmid\ttitle\t$sent\n";
			
			}else{
				#print "in else\n";
				$temp_sent .= "." + $sent;
				print OUT "$pmid\ttitle\t$temp_sent\n";
				$temp_sent = "";
			}
		}
	}
	close FILE;
	close OUT;

	return $output;	
}	
