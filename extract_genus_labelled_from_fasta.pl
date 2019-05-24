


$database_file = $ARGV[0];


$starting_node  = 6960;


###############################################################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object




###############
store_nodes();#
###############



###################
parse_namesfile();#
###################


$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################



###############################################################################################################################################



	open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";

open(OUT1 ,  ">genus_labelled.fas") || die "";
open(OUT2 ,  ">not_genus_labelled.fas") || die "";

	#open(OUT1, ">$database_file.m_rm") || die"";
	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
				$entry_number++;
				if($entry_number =~ /000$/)
					{print "$entry_number\n"};

				########################################
				process_entry($fasta_entry);#
				########################################
				}
			$fasta_entry = $fasta_line;
			}else{
			$fasta_entry .= $fasta_line;
			}
		};
	close(FASTA_IN);

	# do the last entry!:
	unless(length($fasta_entry)<=2)
		{
		$fasta_entry =~ s/^>//;
		########################################
		process_entry($fasta_entry);#
		########################################
		}

close OUT1; close OUT2;



	sub process_entry
		{
		my $line = shift;
		my $current_id = "";
		if ($line =~ /^(.+)\n/ )
			{
			$current_id = $1;
			$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
			}else{
			die "\nerror 5631.\n"
			};

		# need to know if current label is a genus , or not.

		if($current_id =~ /^([A-Z][a-z]+)_/)
			{
			my $label = $1;
			my $taxon_info = $complete_lineage_for_this_species{$label};

			if($taxon_info =~ / genus:$label /)
				{
				print OUT1 ">$current_id\n$line\n";
				# print "\nlabel:$label taxon_info:$taxon_info\n";
				}else{
				print OUT2 ">$current_id\n$line\n";
				
				};

			}else{
			# not labelled with genus ...
			print OUT2 ">$current_id\n$line\n";

			}

	#	print OUT1 ">$current_id\n$line";
		}










#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub store_nodes
{

my %rank_hash;

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "cant open nodes.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database .... nodes.dmp .... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

			$ncbi_nodes{$tax_id}{rank} = $rank;
		#	$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode;
			$ncbi_nodes{$tax_id}{parent} = $parent_tax_id;
			$ncbi_nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

#my @ranks = keys %rank_hash;@ranks = sort @ranks;
#print "ranks found in nodes.dmp:\n";
#print LOG "ranks found in nodes.dmp:\n";

#foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

my @all_nodes = keys %ncbi_nodes;@all_nodes = sort @all_nodes;

print scalar @all_nodes , " nodes have been read.\n";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "cant open names.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


print "\nnames.dmp, parsing 'scientific name', ignoring others ... ";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;
		$ncbi_nodes{$tax_id}{name} = $name;

		$names_line_counter++;#print "$names_line_counter\n";


		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);

print "$names_line_counter names parsed.\n";


}



#####################################################################################################
#
#
#
#####################################################################################################




sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	}



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/; ####################### sep2013
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	# print "child_complete_lineage:$child_complete_lineage\n";
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################







