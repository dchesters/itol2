#!/usr/bin/perl

# 
# 
# to run:
# perl fetch_entries_from_genbank_flatfile.pl [ncbi_flatfile_db] [output_name] [search_string]
# 
# 
# for example:
#  perl ~/usr_scripts/fetch_entries_from_genbank_flatfile.pl proteobacteria_orders.gb Desulfobulbus 320352146
# 
# if your search includes spaces, run like:
#  perl ~/usr_scripts/fetch_entries_from_genbank_flatfile.pl sequence.gb.partials "tax.partial_mt" "tax"
# 
# script simpy prints entries which contain the string given by the user.
# required since relevent blast+ tools dont use ncbi flatfile format (curiously)
# 
# 
# 
# 
# change log 
# 2017-05-18: allows strings containing spaces to be used as search term
# 
# 
# 
# 
# 
# 
# 
# 
###############################################################################################################



my $flatfile 		= $ARGV[0];
my $outfile_name 	= $ARGV[1];
my $search_query 	= $ARGV[2];


$how_many_to_print 	= 1;# of the number in the db matching your query, how many do you want?
$case_sensitive 	= 1;


#########################################################################################################






unless($flatfile 	=~ /[\w\d]/ && 
	$outfile_name 	=~ /[\w\d]/ && 
	$search_query 	=~ /[\w\d]/)
	{
	die "\nerror, something missing in your command.\n"
	};

if($outfile_name =~ s/\s+/_/g){print "\nspaces in output file name, changed to underscores\n"};

print "you speciefied:
	flatfile:($flatfile) 
	outfile_name:($outfile_name) 
	search_query:($search_query)\n";



$entries_matching_search_term=0;

#######################################
parse_genbank_flatfileNEW($flatfile);#
#######################################



print "
entries_matching_search_term:$entries_matching_search_term

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


open(MAIN_OUT, ">$outfile_name") or die "cantopen\n";


for my $index(0 .. ($#file_as_array-1))
	{
	my $entry = $file_as_array[$index];
	parse_this_entry($entry);	
	}

close(MAIN_OUT);

}





###########################################################################



sub parse_this_entry
{
my $current_entry = shift;
my $current_entry_copy = $current_entry;
$counttotal++;

my $accession = "";

my $query_found =0;
if($case_sensitive == 0)
	{
	if($current_entry =~ /$search_query/i)
		{$query_found=1};
	};
if($case_sensitive == 1)
	{
	if($current_entry =~ /$search_query/)
		{$query_found=1};
	};

if($query_found == 1)
	{
	print MAIN_OUT "$current_entry\n\/\/\n"; 
	$entries_matching_search_term++;

	if($entries_matching_search_term >= $how_many_to_print)
		{
		close MAIN_OUT;
		print "\n in the script you only specified $how_many_to_print should be printed,\n" , 
		"this number has been found, so quitting now before reading rest of flatfile DB\n";
		die"";
		}

	}




}


###########################################################################





