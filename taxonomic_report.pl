#!/usr/bin/perl

###################################################################################################################################
#
#
#	taxonomic_report.pl
#
#    	Copyright (C) 2019 Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
# ##################################################################################################################################
# 
#
# 
# original output format for fasta file seems to be that format they use in prokaryote profilling, and SAP
# 
# -seqfile XX -output XX -node XX -fasta or -newick
# 	newick reading disc
# 
# 
# 
# CHANGE LOG
# 2015-09-29: does newick termnials 
# 2015-10-07: does it better
# 2015-11-03: improved format of fasta ID
# 2016-01-12: prints counts of each rank
# 2017-09-08: improved genus name parsing of fasta IDs
# 2019-04-24: prints comm-like table so can caluate PD for different higher ranks
# 
# 
# 
# 
# curiously, newick only works if adding more than one taxonomic string
# 
# 
# 
# 





my $arg_string  = join ' ', @ARGV;

# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose
$filter_level 		= 1;# 1= remove member if it is of a genus not found in tax db. 2= rm if member is of species ....
				# note, currently does not allow option 0
$process_backbone_tree_terminal_IDs = 0;
	# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis
	# if == 2, then trim off species names

$truncate_new_taxonomic_strings = 0;
$fasta_ID_format =2; # applies to fasta IDs, not newick. 1= SAP, 2 = regular

#@ranksarray 	= ("kingdom" , " phylum" , " order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
#my @ranksarray 	= ("kingdom" , " phylum" , " order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
#@ranksarray 	= ("kingdom" , " phylum" , " order", " family");
#@ranksarray 	= (" class" , " order");
#@ranksarray 	= (" order", " family");#, " family"
#@ranksarray 	= (" order", " family");#, " family"



##########################################################################################################################################


# globals:
%check_hash4;

$treefile;
$keyfile;
$fas_file;
$reference_file;
$output_filename;# 	= "$treefile.diet_clades";
%species_names;
$starting_node;
$support_cutoff = "NA";


#####################################
read_command_arguments($arg_string);#
#####################################




###############################################################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object

%assign_these_new_terminals_to_taxon;


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
print "traverse nodes 1 ... ";
######################################################
traverse_nodes($starting_node , $starting_name); # 
######################################################
print " is done\n";





#######################################################################################################################



if($assign_to_tree==1)
	{
	open(INTREE, $fas_file) || die "\nerror ... input file you have given ($fas_file) cannot be opened. \n";
	while (my $line = <INTREE>)	{if($line =~ /(.+\(.+)/){$newick_string = $1}};
	close(INTREE);
	$newick_string =~ s/ //g;$tree_length = length($newick_string);
	print "\nnewick string of length $tree_length has been read from file:$treefile\n\n";

	while($newick_string =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/){$count_scientific_notation_branchlengths_removed++}; 
	# remove regular branchlengths: 0.02048
	while($newick_string =~ s/\:\-*\d+\.\d+//){$count_regular_format_branchlengths_removed++}; 
	# remove 0 length branchlengths
	while($newick_string =~ s/\:\d+//){$count_zero_length_branchlengths_removed++};
	while($newick_string =~ s/(\))\d\.\d+/$1/){$count_proportional_node_support_removed++}; 

#############################################################################################################

print "parsing newick string ....\n";

	while ($newick_string =~ /\(([^\(\)]+)\)/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
		{
		my $node = $1; my $nodeID = "INTERNAL_NODE_$interal_node"; my $length_newick_string = length($newick_string);
		my @child_nodes = split /\,/ , $node;
		if(scalar @child_nodes == 1)# an error, internal node with no split
			{
			if($node =~ /^INTERNAL_NODE_\d+$|INTERNAL_NODE_\d+\:\d+\.\d+$/)
			{
			$newick_string =~ s/\(([^\(\)]+)\)/$1/;
			print "why only one child?\n";
			}else{
			print "some other .... length_newick_string:$length_newick_string ID:$nodeID\n";
			if($length_newick_string <= 500){print "\tnode:$node\n";};
			$interal_node++;
			$root_node = $nodeID;

			};

		}
	elsif(scalar @child_nodes == 2)
		{
	#	print "bifurcation $#child_nodes\n";

		$child_counts{$nodeID} = $#child_nodes;

		for $i(0 .. $#child_nodes)
			{
			$nodes{$nodeID}{$i} = $child_nodes[$i];$nodes{$child_nodes[$i]}{parent} = $nodeID;
			unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++}
			};

		$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;

		$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
		$interal_node++;
		}else{
		print "POLYTOMIE, child nodes:$#child_nodes, IN:$interal_node\n";

		while($node =~ s/([^\,]+)\,([^\,]+)/INTERNAL_NODE_$interal_node/)
			{
			my $member1 = $1; my $member2 = $2;			
		#	print "\tmember1:$member1 member2:$member2, repl:INTERNAL_NODE_$interal_node\n";
			$child_counts{$nodeID} = 1;
				$nodes{$nodeID}{0} = $member1;$nodes{$member1}{parent} = $nodeID;
				$nodes{$nodeID}{1} = $member2;$nodes{$member2}{parent} = $nodeID;
			if($node =~ /([^\,]+)\,([^\,]+)/)
				{
				print "\t\tnothing\n";
				}else{
				$newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/;
				print "\t\tnewick string replacment with INTERNAL_NODE_$interal_node\n";
				};
			$root_node = $nodeID;
			$interal_node++;$nodeID = "INTERNAL_NODE_$interal_node";
			};

			
		};
	}#while ($newick_string =~



my @terminal_labels = keys %terminals;@terminal_labels = sort @terminal_labels;

$fas_file = "$fas_file.fas";

open(TEMPFILE, ">$fas_file");
foreach my $tax(@terminal_labels)
	{print TEMPFILE ">$tax\n"};
close TEMPFILE;

	};

###################################################################################################################################

# my %terminals = ();
# my $root_identity = "";
# my $date = localtime time;

%query_IDs;	# the fasta file contains names of new entries to be placed onto backbone tree. 
		# each are stored in the hash %query_IDs.


open(LOG, ">$output_filename.LOG") || die "\nerror opening output\n";

print "reading fasta file ... ";
	############
	read_fas();# reads fasta IDs, parses a name and retrieves lineages, counts all taxon at all ranks
	############
print " done\n";




# 2019-04: need to print a comm-like table, higher taxon vs all members


my @all_orders = keys %all_orders_hash;@all_orders = sort @all_orders;
my @all_fasta_members = keys %query_IDs; @all_fasta_members = sort @all_fasta_members;


open(OUTCOMM, ">commlike_table");
my $printstrnig = "order\t";
foreach my $sp(@all_fasta_members)
	{$printstrnig .= "$sp\t"};
$printstrnig =~ s/\t$//;print OUTCOMM "$printstrnig\n";

foreach my $order(@all_orders)
	{
	my $printstrnig = "$order\t";	# print "order:$order\n";
	foreach my $sp(@all_fasta_members)
		{
		if($order_of_member{$sp} eq $order)
			{
			$printstrnig .= "1\t";#	print "$sp belongs to $order\n";	
			}else{
			$printstrnig .= "0\t";				
			};
		};	
	$printstrnig =~ s/\t$//;print OUTCOMM "$printstrnig\n";
	};
#	if($rank eq "family"){$family_of_member{$id} = $tax};
close OUTCOMM;

my @all_orders = keys %all_families_hash;@all_orders = sort @all_orders;

open(OUTCOMM, ">commlike_table.fam");
my $printstrnig = "family\t";
foreach my $sp(@all_fasta_members)
	{$printstrnig .= "$sp\t"};
$printstrnig =~ s/\t$//;print OUTCOMM "$printstrnig\n";

foreach my $order(@all_orders)
	{
	my $printstrnig = "$order\t";	# print "order:$order\n";
	foreach my $sp(@all_fasta_members)
		{
		if($family_of_member{$sp} eq $order)
			{
			$printstrnig .= "1\t";#	print "$sp belongs to $order\n";	
			}else{
			$printstrnig .= "0\t";				
			};
		};	
	$printstrnig =~ s/\t$//;print OUTCOMM "$printstrnig\n";
	};
#	if($rank eq "family"){$family_of_member{$id} = $tax};
close OUTCOMM;




$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "
going through tax heirarchy and printing counts...
name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;
print "traverse nodes 2 ... ";



open(OUT, ">$output_filename") || die "\nerror opening output\n";
######################################################
traverse_nodes2($starting_node , $starting_name);#  taxonomic tree, not phylogeny
######################################################
close OUT;




my @ranks = keys %rank_count; 
foreach my $rank(@ranks)
	{
	my $string = $rank_count{$rank};$string =~ s/^\t+(.+)\t+$/$1/;
	my @array = split /\t+/ , $string;
	my $counttax = scalar @array;

		print LOG "rank:$rank count tax:$counttax\n";	
	};

print " is done\n";


# for making standard file of all ranks as columns and listed for each sp

open(OUT7, ">tabulated_taxa");

my $printstring = "member\t";
foreach my $rank(@ranks)
	{
	unless($rank eq "no_rank")
		{
$printstring .= "$rank\t";
		}
	};		
$printstring =~ s/\t$//;
		print OUT7 "$printstring\n";

foreach my $id(@all_fasta_members)
	{
	print "\n$id\n";my $printstring = "$id\t";

	foreach my $rank(@ranks)
		{
		unless($rank eq "no_rank")
		{
		my $tax_at_rank = "NA";
		if($higher_tax_of_member{$id}{$rank} =~ /\w/)
			{
			$tax_at_rank = $higher_tax_of_member{$id}{$rank};
			print "\trank:$rank\t$tax_at_rank\n";

			};
		$printstring .= "$tax_at_rank\t";
		};
		};
	$printstring =~ s/\t$//;print OUT7 "$printstring\n";
	};
close OUT7;

close LOG;


	print "user speified base node for taxonomies:$starting_name\n";
	print "\n\n\nFASTA read and taxonomies assigned. end of script\n";
	exit;

#	};







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";

my $count_diet_IDs = 0;
#open(OUT1, ">$output_filename") || die"\nerror 140\n";
my $fasta_entry = "";
while(my $fasta_line= <IN>)
	{
	if($fasta_line =~ /^>.+/)
		{
		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;
			$process_entry_number++;
			if($process_entry_number =~ /000$/)
				{print "process_entry_number:$process_entry_number\n"};

			########################################
			process_entry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}
	};


unless(length($fasta_entry)<=2)
	{
	$fasta_entry =~ s/^>//;
	########################################
	process_entry($fasta_entry);#
	########################################
	}

close IN;
#close OUT1;




sub process_entry
{
my $line = shift;
my $id = "";
if ($line =~ /^(.+)\n/ )
	{
	$id = $1;
	$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
	}else{
	die "\nerror 5631.\n"
	};

	# if($id =~ s/^Q([A-Z][a-z]+_[a-z]+)/$1/){print "WARNING, leading character Q trimmed from query\n"};

	my $binomial = $id; $binomial =~ s/^([A-Z][a-z]+_[a-z]+)[\|\_].+/$1/;
	$store_binomials{$binomial} = 1;

	my $genus = $id;#print "id:$id\n";
		if($genus =~ s/^([A-Z][a-z]+)_[a-z]+$/$1/i)
			{
			$count_binomials++;
		#	}elsif($genus =~ s/^[a-z]+_([A-Z][a-z]+)_[a-z]+/$1/i)
		#		{
		#		$extended_ids++;
				}elsif($genus =~ s/^([A-Z][a-z]+)_[a-z]+_[0-9]+$/$1/){
				# binomial followed by GI: Lasiopogon_cinctus_158267666
				}elsif($genus =~ s/^([A-Z][a-z]+)_[a-z]+_[A-Z]{2}[0-9]{6}$/$1/){
				# binomil followed by accession: Abagrotis_dickeli_KJ383120

				}elsif($genus =~ s/^([A-Z][a-z]+)_sp_.+|^([A-Z][a-z]+)_cf.+|^([A-Z][a-z]+)_nr.+/$1/)
				{
				}elsif($genus =~ s/^([A-Z][a-z]+)_.+/$1/) # or_ID:(Abraxas_QMA414513_BD00944935)
				{


			#	print  "\nWARNING 153. not parsed iD:$genus\n"	
				};
		$count_diet_IDs++;
		$query_IDs{$id}= 1;



		my $look_for_ID; # 1= remove member if it is of a genus not found in tax db. 2= rm if member is of species ....
		if($filter_level == 1){$look_for_ID = $genus};
		if($filter_level == 2){$look_for_ID = $id};

		if(exists($complete_lineage_for_this_species{$look_for_ID}))
			{
			my $completelineage = $complete_lineage_for_this_species{$look_for_ID};
		#	print "completelineage:$completelineage\n";

			while($completelineage =~ s/(\S+):(\S+)//)
					{
					my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
					$tax_count{$tax}++;
					unless($rank_count{$rank} =~ /\t$tax\t/)
						{$rank_count{$rank} .= "\t$tax\t"};

					$higher_tax_of_member{$id}{$rank} = $tax;

				#	print "rank:$rank tax:$tax\n";
					if($rank eq "order"){
						$order_of_member{$id} = $tax;
						$all_orders_hash{$tax} = 1;
						 # print "\t$tax\n";
						};
					if($rank eq "family")
						{
						$family_of_member{$id} = $tax;
						$all_families_hash{$tax} = 1;

						};

					};


			$present++;#print "fasta ID is PRESENT from taxonomy DB:($id)\n";
			}else{
			$absent++; 
			print LOG "fasta ID is MISSING from taxonomy DB:($id), look_for_ID:($look_for_ID)\n";
			}

};



print "

NUMBER OF FASTA ENTRIS:$process_entry_number
	fasta IDs present in tax DB:$present
	fasta IDs absent  in tax DB:$absent
";


my @k8 = keys %store_binomials;

print "
count entries which are binomials:$count_binomials 
";
print "number of binoaisl: " , scalar @k8 , "
";




}#sub read_fas





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub get_species_list
{
my $list_of_terminals = shift;

my %hash_of_terminals = ();

if($verbose == 1){print "list_of_terminals:$list_of_terminals\n"};

my @array_of_terminals = split /\t/, $list_of_terminals;

foreach my $termnal(@array_of_terminals)
	{
	$termnal =~ s/^([^_]+)_.+/$1/;
	$hash_of_terminals{$termnal}= 1;
	}

my @array_of_terminals2 = keys %hash_of_terminals;

if($verbose == 1){
#print "list_of_species:@array_of_terminals2\n";
}

my $list_of_species = join ' ', @array_of_terminals2;

return($list_of_species);

}



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

# -seqfile -output -node -fasta

if($arguments =~ /-input\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	die "\ncommand error, no -input [input_file_name]\n"
	}

if($arguments =~ /-fasta/)
	{
	$inputformat = "fas";
	}else{
	if($arguments =~ /-newick/)
		{
		$assign_to_tree=1;
		}else{
		die "error";
		};
	}


if($arguments =~ /-output\s+(\S+)/)
	{
	$output_filename = $1;
	}else{
	}



if($arguments =~ /-node\s+(\d+)/)
	{
	$starting_node = $1;
	}else{
	print "user did not given NCBI taxonomic number. using default 33208 (Metazoa)
if you have references outside this default, you need to get the appropriate NCBI taxonomy number from http://www.ncbi.nlm.nih.gov/taxonomy/
and input this here with the option -node taxon_number
";
	$starting_node = 33208;
	}


#$output_filename	= "$treefile.query_clades";


print "\n
user options have been read.....
support_cutoff:			$support_cutoff
taxon node encompassing refs:	$starting_node
treefile:			$treefile
query fasta file:		$fas_file
reference_file:			$reference_file

output will be written to file:	$output_filename


";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub read_keyfile
{
open(IN, $keyfile) || die "\n\nerror 91. cant open $keyfile\n";
my $species_strings_read=0;

while (my $line = <IN>)
	{
	#print "\n$line";

	#MeeRBBegrac Berberis_gracilis 258166 species no_rank:eudicotyledons no_rank:eudicotyledons order:Ranunculales family:Berberidaceae genus:Berberis species:gracilis

	if($line =~ /^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(.+)/)
		{
		my $tobycode = $1;my $species_name = $2;my $ncbi_number = $3; my $rank = $4; my $taxonomic_path = $5;
		#print "\ttobycode:$tobycode species_name:$species_name ncbi_number:$ncbi_number rank:$rank taxonomic_path:$taxonomic_path\n";
		$species_names{$tobycode} = $species_name;$taxonomic_paths{$tobycode} = $taxonomic_path;
		$species_strings_read++;
		$tobycodes{$species_name} = $tobycode;

	#	if($tobycode eq "MLAA4"){print "$tobycode $species_name ... quit\n\n";die ""}

		}
	}

close(IN);

if($species_strings_read == 0){die "\n\nerror 112. seems to be problem reading keyfile, lines in that file dont match regex.\n"}

print "$species_strings_read species strings read from file, plus associated taxonomic information\n";

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
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;

	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

	if($name_string =~ /Zorochros/)
		{
		#print "name_string:($name_string) child:$child complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";
		};


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

















sub record_tree3
{
my $t5 = shift;

my $tree1= "";
open(IN, $t5) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at tree:$treefile\n";

if($tree1 =~ s/\:[0-9\.]+//g){print "\nwarning, branhc lengths removed\n"};
$tree1 =~ s/ //g;


my $newick_string 	= $tree1;
my $interal_node	= 0;



# new newick parser ....

while ($newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";

	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas


	my @child_nodes = split /\,/ , $node;
	$child_counts3{$nodeID} = $#child_nodes;



	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes3{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes3{$child_nodes[$i]}{parent} 	= $nodeID;

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1}

		}
	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";

	$root_node3 = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;

	print "\nnewick_string:$newick_string\n";
$newick_string =~ s/\((INTERNAL_NODE_\d+)\)/$1/;

	}


print "your tree has been read, it has $interal_node nodes.\n";
#print "newick_string:$newick_string\n";die;


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################

sub traverse_nodes2
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
		my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
		my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
		my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
		my $column = scalar split /\s/ , $child_complete_lineage;

		if($tax_count{$name_string}>=1)
			{
			for my $j(1 .. $column){print OUT "\t"};
			my $counttax = $tax_count{$name_string};
			print OUT "$rank:$name_string:$counttax\n"
			};

			###########################################
			traverse_nodes2($child , $child_taxstring);#
			###########################################
		}
	}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################




