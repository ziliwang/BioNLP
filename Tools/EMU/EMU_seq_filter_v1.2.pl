#! /usr/bin/perl -w
use strict;
use Getopt::Std;
use LWP::Simple;
use POSIX qw(ceil floor);
use AACONVERSION qw/ threetoone onetothree nametothree iscodon mutation_conversion /;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
unless (@ARGV == 2)
{
	print "Argument error there should be only 2 arguments\n";
	usage();
	exit;
}

my (%seen_genes, %seen_gids, %seen_prots, %gene_list) = ();

open(IN, "$infile") || die "Can't open file: $infile";
open(OUT, ">$outfile") || die "Can't open file: $outfile";
my $line = 1;

printf OUT "pmid\torganism\tmut_pat1\tpos_patt\tmutation type\twtaa\tmtaa\tpos\tgenes\ttype\tfasta_check\tseq_filter type\tgi\tgene_name\tprot_id\n";

my $tmp = <IN>;

while(<IN>)
{

    chomp;

	my ($pmid, $org, $mutpat, $pospat, $variant_type, $wtaa, $mtaa, $position, $genepat, $mut_type, @others) = split/\t/;
	
    unless ($pmid =~ /\d+/) #Skips rows with non_digit values (intended to skip the header)
    {
		next;
    }
    ($wtaa, $mtaa) = mutation_conversion($wtaa, $mtaa) unless ($variant_type eq "INDEL");
	# col 11 - 15 coorespond to Type, Fasta_check, Gene IDs, Gene Names, Proteins
    my ($col11, $col12, $col13, $col14, $col15);
    $col11 = getcol11($wtaa,$mtaa);
    if($position && $genepat && $wtaa && $mtaa && ($position =~ /^\d+$/))
    {
		print "$line-->$genepat<--\n";
		my @genes = split(/;/, $genepat);
		my %good_genes = ();
		my %good_proteins = ();
		my %good_gid = ();
		my $overallcheck = 0;
		my $max_len = 0; #must remain in the outer loop because of multiple genes per record
		
		foreach my $gene (@genes)
		{
			my @array = collected_data(\%seen_genes, $gene, \&getgid);
			#print "results have: ".@array." elements: @array for $gene\n";
			foreach my $gid (@array)
			{
				my @prots = collected_data(\%seen_gids, $gid, \&getprots);
				   
				foreach my $prot (@prots)
				{
					my @fasta = collected_data(\%seen_prots, $prot, \&getfasta);
					
					my $fastcheck = '';
					foreach my $fast (@fasta)
					{
						$fast =~ s/(\w{70})\n/$1/g; #There are 70 amino acids per line, this makes it so the is only one line of amino acids
						my ($desc, $aa) = split(/\n/, $fast);
						my ($junk, $gi, @others) = split(/\|/, $desc);
						my $len = length $aa;
						
						#print "LENGTH $len\n";
						$max_len = $len if ($len > $max_len);
						
						my $wt = convert($wtaa); #  ensures the aa is the one letter abbreviation as this is what checkrecord requires
						my $mt = convert($mtaa); #  ensures the aa is the one letter abbreviation as this is what checkrecord requires
						$fastcheck = checkrecord($position, $wt, $mt, $aa);
						$overallcheck = 1 if ($fastcheck !~ /^NO$/);
						
						$good_proteins{$gi} = "${gi}|${fastcheck}" if ($fastcheck !~ /^NO$/);
						$seen_prots{$prot}{$fast} = "nada";
					}
					
					$good_genes{$gene} = "nada" if ($fastcheck !~ /^NO$/);
					$good_gid{$gid} = "nada" if ($fastcheck !~ /^NO$/);
					$seen_gids{$gid}{$prot} = "nada";
				}
				if ($position > $max_len)
				{
					$col11 = "DNA" unless ($col11 =~ /PROTEIN/i);
				}
				#print "MAX LENGTH $max_len\n";
				$seen_genes{$gene}{$gid} = "nada";
			}
			$gene_list{$gene} = "nada"
		}
		$col12 = 'YES' if ($overallcheck == 1);
		$col12 = 'NO' if ($overallcheck == 0);
		$col15 = join(";", values %good_proteins);
		$col14 = join(";", keys %good_genes);
		$col13 = join(";", keys %good_gid);
    }
    else
    {
		$col12 = "BAD DATA";
    }
    
	$col11 = $others[-1] if (defined($others[-1]) && ($others[-1] =~ /\w+/));
	my @tojoin = ($pmid, $org, $mutpat, $pospat, $variant_type, uc $wtaa, uc $mtaa, $position, $genepat, $mut_type, $col12, $col11, $col13, $col14, $col15);
    

    
    my $size = @tojoin;
    my $count = 1;
    foreach my $ele (@tojoin) # I didn't use the join function because it was printing out errors about some of the fields being uninitialized (which is acceptable) 
    {
		if (defined($ele))
		{
			print OUT "$ele"
		}
		print OUT "\n" if ($count == $size);
		print OUT "\t" unless ($count == $size);
		$count++;
    }
    
    $line++;
}


sub usage
{
	print << "STOP";

USAGE: perl $0 inputfile.txt outputfile.txt

For additional help see: perldoc $0

STOP
}

sub collected_data
{
    my $hashref = shift;
    my $key = shift;
    my $subref = shift;
    my @array;
    #print "the key in collected data is -->$key<--\n";
    #unless (defined($hashref->{$key}))
    unless (defined $hashref->{$key})
    {
		#print "Collecting Data\n";
		@array = $subref->($key);
		foreach my $ele (@array)
		{
			$hashref->{$key}{$ele} = "nada";
		}
    }	
    else
    {
		#print "Using stored data\n";
		@array = keys %{ $hashref->{$key} }
    }
    return @array;
}

sub getcol11
{
    my ($wtaa,$mtaa) = @_;
    if ($wtaa && $mtaa)
    {
		my $aminowt = iscodon($wtaa);
		my $aminomt = iscodon($mtaa);
		$wtaa = $aminowt if (defined($aminowt));
		$mtaa = $aminomt if (defined($aminomt));
		unless (($wtaa =~ /^(A|T|G|C)$/) && ($mtaa =~ /^(A|T|G|C)$/))
		{
			return "PROTEIN";
		}
		else
		{
			##
			##
			#
			# CHECK FASTA HERE (Decided not to check fasta here)
			#
			##
			##
			return "UNKNOWN";
		}
    }
    else
    {
		return "UNKNOWN";
    }
}

sub convert
{
    my $wtaa = shift;
    my $wt;
    if ($wtaa =~ /\w{4,}/)
    {
	$wt = nametothree($wtaa);
	$wt = threetoone($wt);
    }
    elsif ($wtaa =~ /\w{3}/)
    {
        $wt = threetoone($wtaa);
    }
    else
    {
	$wt = $wtaa;
    }
    return $wt;
    
}

sub getprots
{
	my $genetorefseq = "refseq.txt";
	my $gid = shift;
	open(REFER, "$genetorefseq") || die "Can't open file: $genetorefseq";
	my @results = ();
	while(<REFER>)
	{
		chomp;
		unless (/^#/)
		{
			my ($taxid, $geneid, $stat, $NucAcc, $NucGI, $ProtAcc, $ProtGI, @others) = split/\t/;
			if ($geneid == $gid)
			{
				push(@results, "$ProtGI");
			}
		}
	}
	close REFER;
	return @results;
}

sub getfasta
{
    my @gis = @_;
	
	my ($url,$base,$id2,$gi,$result, $email, @array);
	$email = "tsuznic1\@jhu.edu"; #used with efetch so that NCBI can contact me if there are problems                
	$base="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=";
	if (@gis > 1)
	{
		my $num_of_GIs = 200;
		my $int = 1; 
		foreach my $gi (@gis)
		{
		
			$base=$base.$gi.","; # Ultimately adds a comma separated list of GIs                                         

			#Requests of multiple genpept files (>1000) always fail, This breaks it up into calls of $num_of_GIs per URL
			if (($int % $num_of_GIs) == 0)
			{
				chop $base;
				my $url=$base."&retmode=text&rettype=fasta&email=$email";
				my $result = get($url);
				sleep 1; #NCBI requires a minimum of 3 seconds between large requests so that the system is not overloaded                                                                                                               
				my $fail_count = 0;
				until (defined($result))
				{
					if($fail_count > 9)
						{
								die "ERROR: $url failed at $fail_count times";
						}
						$result = get($url);
						sleep 1;
						$fail_count++;
				}
				my @fastas = split(/\n\n/, $result);
				push(@array, @fastas);
				$base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=";
			}
			$int++;
		}

		# This is required so that a call isn't made with no GIs                                                         
		if (($int % $num_of_GIs) != 1)
		{
			chop $base;
			my $url=$base."&retmode=text&rettype=fasta&email=$email";
			my $result = get($url);
			unless (defined($result))
			{
				print "$url failed in $0\n";
			}
			my @fastafiles = split(/\n\n/, $result);
			push(@array, @fastafiles);
		}
		return @array;
	}
	else
	{
		my $url=$base.$gis[0]."&retmode=text&rettype=fasta&email=$email";
		my $result = get($url);
		sleep 3; #NCBI requires a minimum of 3 seconds between large requests so that the system is not overloaded
		my $fail_count = 0;
		until (defined($result))
		{
			if($fail_count > 9)
				{
						die "ERROR: $url failed at $fail_count times";
				}
				$result = get($url);
				sleep 3;
				$fail_count++;
		}
		chomp($result); chomp($result);
		return $result;
	}
}

sub checkrecord
{
    my ($loc, $wtaa, $mtaa, $seq) = @_;
    return "NO" unless ($loc && ($wtaa =~ /^\w$/) && ($mtaa  =~ /^\w$/) && $seq);
    my @aa   = split //, $seq;
	return "NO" if ($loc > @aa);
    my $loc2 = $loc - 1;
	my @results = ();
    if (($aa[$loc2] =~ /$wtaa/i) )
    {
        push @results, "YES";
    }
    if (($aa[0] =~ /M/i) && ($aa[$loc] =~ /$wtaa/i))
    {
        push @results, "PLUS";
    }
	if (($aa[$loc2] =~ /$mtaa/i))
    {
        push @results, "REV";
    }
	if (($aa[0] =~ /M/i) && ($aa[$loc] =~ /$mtaa/i))
    {
        push @results, "REVPLUS";
    }
    
	if (@results)
	{
		return join("|", @results);
	}
	else
    {
        return "NO";
    }
}


sub getgid
{
	my $gname = shift;
	my %seen = ();
	my @answer = ();
		unless (defined($seen{$gname}))
		{
			if($gname =~ /\b.*\b/)
			{
				#print "$gname<--\n";
				my @queries = ("${gname}[Preferred symbol]+AND+human[Organism]", "${gname}[Gene Name]+AND+human[Organism]",
						"${gname}[Gene Full Name]+AND+human[Organism]", "${gname}[Protein Name]+AND+human[organism]+AND+refseq[Filter]",
						"${gname}[Protein Name]+AND+human[organism]+AND+swissprot[Filter]");
				foreach my $query (@queries)
				{
					unless ($query =~ /Protein Name/)
					{
						@answer = getquery($query, $gname, "gene");
					}
					else
					{
						@answer = getquery($query, $gname, "protein");
					}
					last if (@answer != 0);
				}

				$seen{$gname} = "nothing";
			}
		}
	return @answer;
}

sub getquery
{
        my $query = shift;
        my $gname = shift;
		my $id;
        my $max = 5000;
        my $start = 0;
        my $database = shift;
        my $ebase="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
        my $url = $ebase."esearch.fcgi?db=$database&term=$query&retmax=$max&retstart=$start&usehistory=y";
        my $result = get($url);
        sleep 1;
        $result =~ m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
        my $Count    = $1;
        #print "$Count\n";
        my $QueryKey = $2;
        my $WebEnv   = $3;
        while($result=~/[^<]+<Id>(\d+)<\/Id>/sg){ $id.="$1,"; }
        #print "NOW id is $id\n";
        my $loop_counter = ceil($Count/$max);
        until ($loop_counter <= 1)
        {
                $start += $max;
                $url = $ebase."esearch.fcgi?db=$database&term=$query&retmax=$max&retstart=$start&usehistory=y";
                $result = get($url);
                sleep 1;
                while($result=~/[^<]+<Id>(\d+)<\/Id>/sg){ $id.="$1,"; }
                $loop_counter--;
        }
		my @return;
		return @return unless (defined($id)); 
        chop $id;
        #Fetch the ids we found and  input into an array
        my @arrayids=();
        if($id=~m/\,/){ (@arrayids)=split/\,/,$id; } else{ @arrayids=($id);}
		
        foreach my $protid(@arrayids)
        {
                unless ($protid eq "") {push(@return, "$protid")};
        }
		
		return @return;
}



=pod
    
=head1 Parse Table Script for EMU Output



=head1 Parameters

=head3 First: EMU output file (used for input)

The input file is expected to have the following tab delimited format:

PMID, Organism, mutation pattern, position pattern, wild-type amino acid, mutant-type amino acid, position, genes

Only the PMID, wild-type amino acid, mutant-type amino acid, position and genes columns are used by the program. 
The others can contain anything as long as the field exists. Also any column past the genes column will be handled 
appropriately and will not affect the program, however they will NOT be printed out in the results. This is so that
the output of the program will always have the same number of columns.  

=head3 Second: Name of output file

The format of the output file has the same format as the input file with 6 additional tabe delimited columns. Namely:

Sequence Type, Fasta Check, Gene IDs that pass fasta check, Gene names that pass fasta check, proteins that pass fasta check, OMIM check

Sequence Type can have the values DNA, prot, UNKNOWN or ERROR if something goes wrong.
 


=cut
