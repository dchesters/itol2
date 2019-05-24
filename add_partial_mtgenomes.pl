#!/usr/bin/perl

# 
# 
# to run:
# perl fetch_entries_from_genbank_flatfile.pl [ncbi_flatfile_db] [output_name] [search_string]
# 
# 
# for example:
# perl ~/usr_scripts/fetch_entries_from_genbank_flatfile.pl proteobacteria_orders.gb Desulfobulbus 320352146
# 
# 
# 
# script simpy prints entries which contain the string given by the user.
# required since relevent blast+ tools dont use ncbi flatfile format (curiously)
# 
# 
# 
# CHANGE LOG
# 2018-JAN-26: requires user specifiy taxonomic level, earlier i worked at genus level,
#   seems for some reason i last used this script a species level 
#   
# 
# 
# 
# 
# 
# 
# 




#my $search_query 	= $ARGV[2];



my $file1 		= $ARGV[0];
my $flatfile 		= $ARGV[1];
my $tax_level 	= $ARGV[2];

unless($tax_level =~ /genus/)	{die "\nmissing switch, whihc taxonoic level (only genus supported)\n"};


open(IN, $file1) || die "\nerror 37\n";

print "opened partials log file:$file1\n";

while (my $line = <IN>)
	{
	# print $line;


	if($tax_level =~ /species/) # UNtested
		{
		if($line =~ / tax:([A-Z][a-z]+[\_\s]{0,1}[a-z]+) ac:/)
			{
			my $tax = $1;$tax =~ s/\_/ /;
			$store_genera_with_partial_genomes{$tax}=1;$tax_parsed++;#	print "tax:$tax\n";
			};
		};

	if($tax_level =~ /genus/) # UNtested
		{
		if($line =~ / tax:([A-Z][a-z]+).+ ac:/)
			{
			my $tax = $1;$tax =~ s/\_/ /;
			$store_genera_with_partial_genomes{$tax}=1;$tax_parsed++;#	print "tax:$tax\n";
			};
		};


	}
close IN;

print "tax_parsed:$tax_parsed\n";
unless($tax_parsed >= 1){die "\nerror, nothing parsed.\n"};


@keys = keys %store_genera_with_partial_genomes;




$how_many_to_print = 1;# of the number in the db matching your query, how many do you want?



#########################################################################################################






$entries_matching_search_term=0;

print "reading .gb file of completes\n";

#######################################
parse_genbank_flatfileNEW($flatfile);#
#######################################




		open(OUT5, ">partial_mtgenomes_needed");
		close OUT5;

my @print_array;

foreach $genus(@keys)
	{

	my $count = $partial_genus_found_in_compeltes_file{$genus};
	print "genus with partial:$genus count in completes:$count\n";

	unless($count >= 1)
		{
		$partials_to_add++;
		if($genus =~ /\s/){$genus = "'$genus'"};

	#	open(OUT5, ">>partial_mtgenomes_needed");
	#	print OUT5 "$genus ";
	#	close OUT5;

		push @print_array , $genus;

		}
	};


@print_array = sort @print_array;
my $printstring = join ' ', @print_array;
		open(OUT5, ">>partial_mtgenomes_needed");
		print OUT5 "$printstring\n";
		close OUT5;



print "
partials_to_add:$partials_to_add
\n\nfin.\n";
exit;


##############################################################################################



###########################################################################

sub parse_genbank_flatfileNEW
{
my $inputfile = shift;

open(IN2, $inputfile) || die "cant open infile:$inputfile\n";
print "opened $inputfile\n";
my $fileasstring = "";
while(my $line = <IN2>)
	{$fileasstring .= $line}
close(IN2);


my @file_as_array = split(/\n\/\/\n/ , $fileasstring);


# open(MAIN_OUT, ">$outfile_name") or die "cantopen\n";


for my $index(0 .. ($#file_as_array-1))
	{
	my $entry = $file_as_array[$index];
	parse_this_entry($entry);	
	}

# close(MAIN_OUT);

}





###########################################################################



sub parse_this_entry
{
my $current_entry = shift;
my $current_entry_copy = $current_entry;
$counttotal++;

if($counttotal =~ /00$/)
	{
	print "entry number $counttotal\n";
	};

my $accession = "";



foreach $genus(@keys)
	{
#	print "\tgenus:$genus";

if($current_entry =~ /$genus/)
	{
#	print MAIN_OUT "$current_entry\n\/\/\n"; 
	$entries_matching_search_term++;

	$partial_genus_found_in_compeltes_file{$genus}++;

	if($entries_matching_search_term >= $how_many_to_print)
		{
	#	close MAIN_OUT;
#		print "\n in the script you only specified $how_many_to_print should be printed,\n" , 
#		"this number has been found, so quitting now before reading rest of flatfile DB\n";
#		die"";
		}

	}

	}


#die "";


}; # sub parse_this_entry


###########################################################################





