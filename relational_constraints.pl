#!/usr/bin/perl

###################################################################################################################################
# 

# 
# 
# 
# 
# sub store_nodes reads file nodes.dmp, which contains ncbi taxonomic heirachy, 
# 	with nodes represrted by unique numbers
# sub parse_namesfile parses names.dmp for 'scientific name' of each node
# sub traverse_nodes, recurses through all nodes of the taxonomic heirachy, from user specified position
# 	constructs lineage inforamtion for each species
# sub read_fas, reads users fasta file, stores all taxa included,
# 	retrives lineage for each, and stores all higher taxa.
# sub record_tree2, read user tree in newick format,
# 	some basic processing if required, stores structure of tree,
# 	reads tip taxa, again rretreives their higher taxa and records these.
# sub backbone_constraints2
# 	for each terminal of teh backbone tree,
# 	retrives lineage, then proceeds through all ranks present, starting from genus, in a root-ward direction  
# 	higher taxon names are assigned to the terminal, until the point is readch that the higher taxon name 
# 	is possessed also by a terminal anywhere else on the tree. 
# 	note, in the unlikely event of complete monophylys on the backbone tree,
# 	it would be neccessary only to consult the taxa of the terminals sister to the termnial being assigned.
# the sub prints a newick string just of these assignments, for information.
# next , for the taxon at each terminal, fasta IDs belonging to that terminal are retrived ,
# and , where >1, are joined in comma-sepeatred parenthseis,
# and inserted into the newick string , to replace its encompasing taxon name
# finally, any terminals from taxa which are absent in the sequence data,
# are pruned from the constraint tree.
#
#
# 
# 
# 
# 










# fork of read_deep_level_constraints.pl
# doesnt set constrAINTS that require floating taxa (see fasttree constraints),
# so can be set in strict newick format

# make sure old version of this is kept


# ((Sialis,Stenosialis)_lutaria:1.0,






#
#	backbone_constraints_newick.pl
#	perl script which reads the topology of a phylogeny and the taxonomic groups of the contained clades,
#	then makes an equivelent constraint tree for the species of a set of DNA sequences 
#	suitable for RAXML!
#
#		  
#    	Copyright (C) 2015-2018  Douglas Chesters
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
#	Insecta = 50557, Endopterygota 	= 33392
#
#
#
#
##########################################################################################################################################
#
#
#
#
#	uses code from parse_clades.pl / bagpipe_phylo.pl 
# 	(Papadopoulou et al. 2014 / Chesters et al. 2015)
#	bagpipe_phylo.pl is written to be a general, standalone script for 
#	
#	you need to download the NCBI taxonomy database, go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/
#	download the file taxdump.tar.gz, and unzip in the working directory. 
#	the script will read 2 of the resulting files: names.dmp and nodes.dmp.
# 	it will get from these the higher taxonomies to be used for assignment to species given in your tree.
#	only other input is a fasta file of the queries (-seqfile), the script simply identifies query names from this file 
#
# 	to run:
#
# 	3 OUTPUTS: $treefile.rax_constraint; $fas_file.$treefile.tnt; $fas_file.$treefile.ft
#	
#	
#	
#	
#	
#	
#	Change Log:
#	21 Jan 2015: version 1.001
#	22 oct 2015: modified some bacterial labels when reading ncbi taxonomy db:/^Candidatus\s(\w+)/$1/
#		these are often not used, at least in hte tree  
#		if cant find taxonomy for binomials, tries genus name
#	01 dec 2015: runs without outgroup specified
#		finds taxonomies for genus name when lineage for species is not found
#	28 Jan 2018: implemented another algorithm,
#		an advance over the 'basic constraint method'
#		makes newick style constraint tree, will usually omit non-overlapping taxa,
#		maintains backbone topology, inferrs higher ranks for backbone terminals to the extent possible
#	22 Mar 2018: more helpful error message when user uses out of date ncbi taxonomy database
#		got it working on backbone trees with just genus names as labels
#	01 Nov 2018: user specifies ranks at which constraints are made, in command.
#		avoids 'bug', early version was hacked so order constraints not made,
#		for reasons that might make sense at the time,
#		but which meant inadequate constraining since.
#	21 Dec 2018: prints fasttree format constraints of same set printed to newick.
#		although exsiting scripts available for fasttree format relational constraints,
#		this script might work better in cases of low overlap and missing taxonomies.
#		renamed, backbone_constraints_newick.pl to relational_constraints.pl
#	
#	
#	
#	
#	
##########################################################################################################################################


# back to the drawing board.
# as opposed to selecting a shared taxon, there would often be a few ranks above that that node would represent,
# need to draw these on, then select a criteria for assigning, mybe there would be ambiguity.
# 
# 
# second type of constraint made, monophyly of the assigned name, 
# at first instance root->tip, if contains no other taxa, and >1 of the taxa of rank below.
# 

# 
# 
# 
# SUBS:
# read_command_arguments, store_nodes, parse_namesfile, traverse_nodes, read_fas, record_tree2
# traverse_backbone_tree_and_infer_taxonomic_node_labels {get_terminals_from_this_node2, get_shared_taxa, get_all_taxa}
# traverse_backbone_tree_and_collapse
# 
# 
# 
# 
# 


my $arg_string  = join ' ', @ARGV;



# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose

#$tree_ID_format = 2;# i think 1=Order_Genus, 2= Order_Genus_species


# $process_backbone_tree_terminal_IDs = 2;
# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis
# if == 2, then binomial to genus

$fasttree_format_constraints_file 	= "ftcons"; # output.
$fasttree_constraints_tabdelim 		= "relational_constraints_tabdelim";


$process_new_consraint = 0; # default 0, 
				# somtimes there will be taxa inferred from the input backbone tree
				# which are not in the subsequent dataset,
				# usually rare, but if == 1 , script will try prune them from the 
				# subsequent constraint tree


$ugly_hack = 1;

$unconstrain_broadly_labelled_data = 1;
# much data only ident for order level, usually ignored, although many slip through, this is an extra chaeck
# it tests if the species has a family designation, if not, then not constrained.



$newick_reformat = 1;# 1 = backbone ID are binomial, convert to genus name only






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
$outgroup = "NULL";
$outfile_prefix;



#####################################
read_command_arguments($arg_string);#
#####################################

$output_filename = "$fas_file.$treefile";


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




$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


# hack, some missing species
# in other scripts, the presence of these when they actually are in the ncbi tax DB (e.g. newer db is used), causes error.

#$complete_lineage_for_this_species{"Lemyra_melli"} 		= $complete_lineage_for_this_species{"Lemyra_jankowskii"};
#$complete_lineage_for_this_species{"Promalactis_suzukiella"} 	= $complete_lineage_for_this_species{"Promalactis_jezonica"};
#$complete_lineage_for_this_species{"Hylaeus_dilatatus"} 		= $complete_lineage_for_this_species{"Hylaeus_difficilis"};
#$complete_lineage_for_this_species{"Lethe_dura"} 			= $complete_lineage_for_this_species{"Lethe_diana"};
#$complete_lineage_for_this_species{"Eucryptorrhynchus_brandti"} 	= $complete_lineage_for_this_species{"Eucryptorhynchus_brandti"};




###############################################################################################################################################


open(BCN_LOG, ">backbone_constraints_newick_LOG") || die "";
close BCN_LOG;


my %terminals = ();
my $root_identity = "";
my $date = localtime time;

%count_tax_in_new_data;
%query_IDs;	# the fasta file contains names of new entries to be placed onto backbone tree. 
		# each are stored in the hash %query_IDs.


############
read_fas();#
############



####################################################################################################################################


%nodes;		# user phylogeny stored in this hash
$root_node;	# this global will keep the value of the last in the loop, which will be the root
%count_taxa_represented_in_terminals;	# will use this to determine whether taxa are monophyletic.
$print_internal_nodes=0;


# read backbone tree into hash. node identities for hash key, 
# then parent/child identities as hash entry
# gets complete lineage for each tip, and for tax at all ranks, 
# counts number of tips for each into object count_taxa_represented_in_terminals
# optional processing of tip labels. otherwise, this sub just does bog-standard reading of tree

#########################
record_tree2($treefile);#
#########################

print "
your backbone tree has been read, TIPS:$count_the_terminal_nodes ALL NODES:$print_internal_nodes.
\n";





####################################################################################################################################

# jan 2018 , another constraint method

#########################
backbone_constraints2();#
#########################

print "\nprinting fasttree format relational constraints\n";
$ft_relational_constraints_subcalls =0;

###############################################################
print_fasttree_format_relational_constraints($root_node , 0);#
###############################################################

print "left sub\n";

# relational constraint for each relevent backbone node, now in this object:
my @rel_constraints = keys %relational_contraints;@rel_constraints = sort @rel_constraints;

# something wierd happening, bunce of hash keys of fasta ID appear as 1.
# cant fingure out why, jsut read in again:
open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. \n";
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1; if($id == 1){die "\nWTF\n"};
		$fasta_IDs_yet_again{$id} = 1;
		}
	}
close IN;


my @the_fasta_IDs = %fasta_IDs_yet_again;@the_fasta_IDs = sort @the_fasta_IDs;

	foreach my $id5(@the_fasta_IDs)
		{
	#	print "id5:$id5 ";
		};
# die "";

my $printing_limit = 0;
for my $i(0 .. $#rel_constraints)
	{
	my $rel_constraint_ID = $rel_constraints[$i];
	if($printing_limit <= 6){print "FT relational constraint $i of $#rel_constraints ($rel_constraint_ID)\n"};
	if($printing_limit =~ /00$/){print "FT relational constraint $i of $#rel_constraints ($rel_constraint_ID)\n"};

	$printing_limit++;

	my $rel_constraint = $relational_contraints{$rel_constraint_ID};
#	print "\ni:$i rel_constraint_ID:$rel_constraint_ID contraint:$rel_constraint\n";
	my @split_csontr = split /RelationalConstraint\d/ , $rel_constraint;
	# RelationalConstraint0

	foreach my $id(@the_fasta_IDs)
		{
		unless($id == 1) # im stumped why bunch of 1's appeared and persisting as hash keys,
		{
	#	print "id:$id\n";
		my $state = "-";
		foreach my $j(0 .. $#split_csontr) 
			{
			my $contraint7 = $split_csontr[$j];
			if($contraint7 =~ /\t$id\t/){$state = $j};
			};
		if($state =~ /\d/){$count_states_each_member{$id}++};
		$relational_constraints_per_fasta_ID{$id} .= $state;
		$relational_constraints_per_fasta_ID2{$id} .= "$state\t";
		if($id == 1){die "\nHUH!\n"}
		};
		};

	};

##################

print "\nfasta_IDs:@the_fasta_IDs\n";

open(OUT_REL_CONS, ">relational_constraints.ft_format") || die "";
open(OUT_REL_CONS2, ">relational_constraints.ft_format.tdt") || die "";
$printystring = "member\t";
for my $i(0 .. $#rel_constraints)
	{
	my $rel_constraint_ID = $rel_constraints[$i];$printystring .= "$rel_constraint_ID\t";
	};
print OUT_REL_CONS2 "$printystring\n";

foreach my $id(@the_fasta_IDs)
	{
	unless($id == 1) # im stumped why bunch of 1's appeared and persisting as hash keys,
	{
	my $current_line = $relational_constraints_per_fasta_ID{$id};
	my $current_line2 = $relational_constraints_per_fasta_ID2{$id};
	my $fam = "NA";
	if($family_of_fasta_ID{$id} =~ /\w/){$fam =  $family_of_fasta_ID{$id}};

	print OUT_REL_CONS "$id\t$current_line\n";
	print OUT_REL_CONS2 "$fam $id\t$current_line2\n";
	}};

close OUT_REL_CONS;
close OUT_REL_CONS2;

foreach my $id(@the_fasta_IDs)
	{
	if($count_states_each_member{$id} =~ /\d/){$members_with_ft_rel_constraint++};
	};

my $count_ft_rel_constraints = scalar @rel_constraints;
my $count_fasta_ids = scalar @the_fasta_IDs;

print "
for fasttree relational constraints,
 number of constraints:$count_ft_rel_constraints
 number of members in fasta file:$count_fasta_ids
 of which, number with at least one constraint:$members_with_ft_rel_constraint
";
print "\n\nFIN.\n\n";
exit;

####################################################################################################################################


%all_taxa_on_backbone_tree;
%monophylys;

	# traverse nodes of the user genome tree, 
	# for use in constraining treesearch in the barcode data


$new_newick_string = "($root_node)";
$new_newick_string2 = "($root_node)";
$newick_print_string = "($root_node)";



print "\n
calling sub traverse_backbone_tree_and_infer_taxonomic_node_labels, 
which will traverse nodes of the backbone tree, 
infer taxonomic names for each node
";



#########################################################################
traverse_backbone_tree_and_infer_taxonomic_node_labels($root_node , 0);#
#########################################################################

my @alltax9 = keys %all_taxa_on_backbone_tree;@alltax9 = sort @alltax9;

open(NEWICK1 , ">$outfile_prefix.newick_string1") || die "\nerror 193.\n";
print NEWICK1 "$new_newick_string\n";
close NEWICK1;
open(NEWICK2 , ">$outfile_prefix.newick_string2") || die "\nerror 193.\n";#
print NEWICK2 "$new_newick_string2\n";
close NEWICK2;


my @mps = keys %monophylys;@mps = sort @mps;
foreach my $mp(@alltax9)
	{
	if(exists($monophylys{$mp}))
		{
		print "$mp IS monophyletic\n";
		}else{
		#print "$mp NOT monophyletic\n";
		
		}
}
print "\n$#mps name monophylys found on backbone tree, these will be set for next tree search if relevent 
, in addition to the other things\n";


# taxonomic node labels assigned to backbone tree.
# now traverse tree again and decide at which point to stop (collpase all decendents)
# this must be done where nodes are reached at which not all taxon represented in the fasta file 
# are included as decendents.
# also might need to collapse at where polyphyletic taxa are reached.



print "\ncalling sub traverse_backbone_tree_and_collapse
which will do something.
";


$yet_another_newick_string = "($root_node)";


#########################################################################
traverse_backbone_tree_and_collapse($root_node);#
#########################################################################


#while($yet_another_newick_string =~ s/\,*\(\)//){};
#while($yet_another_newick_string =~ s/\(\,\(/"(("/e){};
while($yet_another_newick_string =~ s/(\,)\,+/$1/){};
while($yet_another_newick_string =~ s/\(\,\)//){};
while($yet_another_newick_string =~ s/\(\,\(/"(("/e){};
while($yet_another_newick_string =~ s/\(\,([A-Z])/($1/){};
while($yet_another_newick_string =~ s/\,\)/)/){};


	# e.g. $yet_another_newick_string = ((Lepismatidae-Nicoletiidae)Zygentoma);
open(NEWICK2000 , ">$outfile_prefix.newick_string2000") || die "\nerror 193.\n";#
print NEWICK2000 "$yet_another_newick_string;\n";
close NEWICK2000;



# this object contains all names of taxa written onto the constraint tree:
# constraints_written


my @constraints_writtenn = keys  %constraints_written;@constraints_writtenn = sort @constraints_writtenn;
foreach my $constraint(@constraints_writtenn)
	{
#	print "constraint:$constraint\n";
	my @split = split /\-/ , $constraint;
	foreach my $item(@split)
		{
	#	print "\t$item\n";
		$all_constraint_taxa{$item} =1;
		};
	};





#my @ft_constraints_again = keys  %ft_constraints;
#@ft_constraints_again = sort @ft_constraints_again;
#open(CONS_TABLE , ">$fasttree_constraints_tabdelim") || die "\nerror 258\n";
#foreach my $constraint(@ft_constraints_again)
#	{
#	#print "constraint:($constraint)\n";
#	$constraint =~ s/\n/\t/g;#
#	$constraintno++;
#	print CONS_TABLE "constraint:$constraintno\t$constraint\n";
#	}
#close CONS_TABLE;




my @new_terminals = keys %query_IDs;@new_terminals = sort @new_terminals;

# my @all_taxa = keys %all_taxa_on_backbone_tree;@all_taxa = sort @all_taxa;


# open(OUT5 , ">$fasttree_format_constraints_file") || die "\nerrro\n";
# print OUT5 " " , scalar @new_terminals , " " , scalar @ft_constraints_again + scalar @mps + 1, "\n";

# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1

# comorpha)Brachycera)Diptera)Diptera-Siphonaptera-Mecoptera)Diptera-Siphonaptera-Mecoptera-Amphiesmenoptera,((S


my $coutn = 0;

open(OUT_LIST, ">$outfile_prefix.list_of_constrained_members") || die "\nerror 410\n";

foreach my $new(@new_terminals)
	{

	$coutn++;
	if($coutn=~ /0000$/)
		{
		print "$coutn of $#new_terminals, ID:$new\n";
		};
	my $genus = $new;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
	
	my $lineage = "NA";

	if(exists($complete_lineage_for_this_species{$new}))
		{
		# complete_lineage_for_this_species, not just species, 
		# but has lineage for every taxon in ncbi, only mod is space->underscore
		$lineage = $complete_lineage_for_this_species{$new};

		}elsif(exists($complete_lineage_for_this_species{$genus}))
			{
			$lineage = $complete_lineage_for_this_species{$genus};
			};

		if($unconstrain_broadly_labelled_data == 1)
			{
			unless($lineage =~ /family\:[A-Z][a-z]+/){$lineage = "NA"; $new_members_without_family_names_ignored++};
			};

		unless($lineage eq "NA")
			{

			# SEP2016, unconstrain these
		#	}else{
		#	die "WARNING ... cant find lineage.\n";
			
		#	print "\n$coutn of $#new_terminals, ID:$new\nlineage:$lineage\n";

			my $current_member_to_be_assigned = "NA";
			while($lineage =~ s/\:(\w+)//)
				{
				my $test_tax = $1;# print "\t\ttax in current lineage:$test_tax\n";
				if(exists($all_constraint_taxa{$test_tax}))
					{
					$current_member_to_be_assigned = $test_tax; # print "\t\t\texist as constraint taxon [re]assigning\n";
					};
				};

			$new_member_assigned_to_constraint_taxon { $new } = $current_member_to_be_assigned;
			print OUT_LIST "$new\n";
			};

	};#$new(@new_terminals)

my @testkeys = keys %new_member_assigned_to_constraint_taxon;

close OUT_LIST;

print "
constraint tax determined for each of $#testkeys new members 
new_members_without_family_names_ignored:$new_members_without_family_names_ignored

";




# finally need to read constraint tree and replace tax with new member lists
# $yet_another_newick_string


place_specieslevel_data_into_taxon_constraints($yet_another_newick_string);


# all labels on input constrqint whic have been string replaced.:
# add_these_to_this

# all_constraint_taxa contains all seperate tax names assigned .. NOT relational constraints types (hyphenated)


$xy = 0;
if($xy == 1)
{
my @tax_on_in_constrt = keys %all_constraint_taxa;@tax_on_in_constrt = sort @tax_on_in_constrt;
print "\nlook from $#tax_on_in_constrt tax on the input constraint tree, check thses are all replaced\n";
foreach my $tax(@tax_on_in_constrt)
	{
	if(exists($add_these_to_this{$tax}))
		{
	#	print "\t$tax replaced\n";
		# taxa as node labeles will be removed at next step so dont need worry about these
		# ony worry abotu input taxa which remain which are terminals.
		}elsif($yet_another_newick_string =~ /\)[A-Za-z\-]*$tax/)
		{}elsif($yet_another_newick_string =~ /$tax[^a-zA-Z\_]/)
		{
		my $context = "NA";
		if($yet_another_newick_string =~ /(........$tax........)/)
			{$context = $1};
		print "warning. terminal taxon ($tax) in your INPUT constraint tree might be retained in output\n";
		print "\tin some cases beacuse no seqs of this tax in your fasta file.\n";
		print "\tafter removing node labels, grep both fasta file and newick, rm if neccessary\n";
		print "\tcontext in new constraint tree:$context\n";
		push @prune_tax , $tax;
		};
	};
};

# try different way::::


my @tax_on_in_constrt = keys %all_internal_node_strings_assigned_to_input_constraint;
@tax_on_in_constrt = sort @tax_on_in_constrt;
foreach my $tax(@constraints_writtenn)
	{
	if($yet_another_newick_string =~ /$tax/)
		{
		print "constraint:$tax\n";$constraints_klinging_on++;

#		if($yet_another_newick_string =~ /\)$tax/)
#			{
		#	print "\tnever mind, just node label\n";
#			}els
			if(
			# $yet_another_newick_string =~ /[^a-zA-Z\_\-\)]$tax[^a-zA-Z\_\-]/

			$yet_another_newick_string =~ /[^a-zA-Z\_\-\)]$tax[^a-zA-Z\_\-]XXXXXXXXXX/


			)
			{
			print "\tmust be terminal.... needs pruning\n";
			$tax =~ s/\-/"\\-"/e;
			push @prune_tax , $tax;
			}elsif($yet_another_newick_string =~ /$tax[\_a-zA-Z0-9]{6}/){
		#	print "\tno problem, its something else\n";
			}else{
			my $context = "NA";
			if($yet_another_newick_string =~ /(.{20}$tax.{20})/){$context = $1};
		#	print "\thuh!\n\tcontext:$context\n";
			};
		};
	};

if($constraints_klinging_on =~ /[\w\d]/)
	{
	print "\nconstraints_klinging_on:$constraints_klinging_on\n";
	};


# if( $ugly_hack == 1){ @prune_tax = ( "Dytiscoidea\\-Haliploidea" , "Cucujoidea\\-Chrysomeloidea") };


if($#prune_tax >= 0 && $process_new_consraint == 1)
	{
print "\nyou have opted to prune ... \n";

	foreach my $tax(@prune_tax)
		{
		print "terminal to prune:$tax\n";
		# context in new constraint tree:lisonia,Archaeopsyllinae,(Polype


	my $testsrting = $yet_another_newick_string;
	my $count_instances_of_tax = 0;
	while( $testsrting =~ s/[\,\(\)]$tax\,//){$count_instances_of_tax++};
	if($count_instances_of_tax >= 2){print "\nwarmning, multipple instances ($count_instances_of_tax) of tax ($tax) in newick string \$yet_another_newick_string\n"};

		if($yet_another_newick_string =~ s/\,$tax(\,)/$1/)		
			{
			print "\teasy to prune, was in a polytomie\n";
			}elsif($yet_another_newick_string =~ /\,\(($tax\,\(.+)/)
			{
			print "\tfound tax in bifurcation\n";
			# nt tree:alinae,(Chrysidini,(Poli	

			my $string = $1;
			my $count_nestedness = 1;my $string_to_replace = "(";my $regex_to_replace = "\\(";
			my $parentheses_encountered =0;

			my $x = 1;
			while($x == 1)
				{
				my $preceeding_string; my $next_parenthesis;
				if($string =~ s/^([^\(\)]*)([\(\)])//)
					{
					$preceeding_string = $1; $next_parenthesis = $2;
					}else{die "\nerror 522\n"};

				#unless($parentheses_encountered == 0){
				$regex_to_replace .= "$preceeding_string" . "\\" . "$next_parenthesis";
				$string_to_replace .= "$preceeding_string$next_parenthesis";
					#};
				$parentheses_encountered++;

				if($next_parenthesis eq "("){$count_nestedness += 1};
				if($next_parenthesis eq ")"){$count_nestedness -= 1};

			#	print "parentheses_encountered:$parentheses_encountered nestedness:$count_nestedness string:$preceeding_string next:$next_parenthesis\n";
				if($count_nestedness == 0)
					{
					$x = 0; # break loop
				#	print "break loop\n";
					};
				};


	# 	(Chrysidini,
	#		(Polistes,Abispa,Vespa,Wallacidia,
	#			(Linepithema,Leptomyrmex,Solenopsis,Camponotus,Polyrhachis,Formica,Myrmica,Pristomyrmex,Vollenhovia,
	#				(Nomada,Colletes,Osmia,Hylaeus,Melipona,Megachile,Lasioglossum,Halictus,Andrena,Seladonia,Sphecodes,Bombus
	#				)Philanthus
	#			)Vespoidea-Apoidea-Sphecoidea
	#		)


			my $replace_with = $string_to_replace;
			unless($replace_with =~ s/(\()$tax\,\((.+)\)/$1$2/){die "\nerror 555.\n"};

			my $following_string = "NA"; 
			if($yet_another_newick_string =~ /$regex_to_replace(.{40})/)
				{$following_string = $1}

		#	print "\nstring_to_replace:$string_to_replace\n";
		#	print "replace_with:$replace_with\n";
		#	print "REGEX_to_replace:$regex_to_replace\n";
		#	print "following 40 chars:$following_string\n";

			# here need to make sure non a-z0-9 characters are appropriate for use in regex,
			# for example parentheses need two backslashes before (done already)
			# maybe some hyphens . others?

			# sometimes have new tax as node label for some reason, just rm
			# us,Haliplus_fulvus,Haliplus_fasciatus)Gyrinus_sp_BOLD_AAG0707_943417329,Dineut
			$yet_another_newick_string =~ s/($regex_to_replace)[A-Z][a-z]+[A-Za-z0-9\_]+/$1/;

			if($yet_another_newick_string =~ s/$regex_to_replace/$replace_with/)#
				{
				print "\nsuccess\n";
				if($yet_another_newick_string =~ /(.{20}$tax.{20})/){my $error = $1;print "\nerror, still there!\n($error)\n"};
				}else{
				print "\nstring_to_replace:$string_to_replace\n";
				print "replace_with:$replace_with\n";
				print "within:$yet_another_newick_string\n";
				die "\nerror 560 trying to prune in newick string\n"
				};

		#	die;

			}else{
			if($yet_another_newick_string =~ /(.{20}$tax[^a-zA-Z\_].{20})/)
				{$contenxt = $1; print "\ncontext:$contenxt\n"}else{print "\ncant even find it\n"};	
			die "\nerror .. this was only a quick hack .. not written for current configuration\n"
			}

		};

	};


open (OUT67 , ">$outfile_prefix.newick_constraint_tree2") || die "\n\n";
print OUT67 "$yet_another_newick_string;\n";
close OUT67;




# do basic constraints, at the genus level .. 

my @backbone_genera = keys %terminals;@backbone_genera = sort @backbone_genera;

# used to assist pruning
# open(TEST_FILE , ">$outfile_prefix.needs_pruning") || die "\nerror 759\n";

my $j = 0;
foreach my $genus(@backbone_genera)
	{
	$j++;
	my $count_new_members_assigned_to_this_genus = 0;
	my $replace_backbone_genus_with_string = "";
	foreach my $new(@new_terminals)
		{
		if($new =~ /$genus[_]/)
			{
			$count_new_members_assigned_to_this_genus++;
			$replace_backbone_genus_with_string .= "$new,";
		#	$print "\tnew member $new\t"
			};
		};
	$replace_backbone_genus_with_string =~ s/\,$//;
	print "Genus $j of $#backbone_genera backbone genera. Setting basic constraints named $genus, " , 
		"found $count_new_members_assigned_to_this_genus\n";
	if($count_new_members_assigned_to_this_genus >= 1)
		{
		print TEST_FILE ">$genus\nACTG\n"
		}else{print "\n\tfound nothing for backbone genus $genus .. will need pruning. quitting.\n"};

	$backbone_copy_newick =~ s/$genus([^a-z])/$replace_backbone_genus_with_string$1/;

	};

close TEST_FILE;


open (OUT699 , ">$outfile_prefix.basic_constraint_tree") || die "\n\n";
print OUT699 "$backbone_copy_newick\n";
close OUT699;





print "

line 801,
FIN!
";


exit;









open(OUT57 , ">$outfile_prefix.taxon_constraints_output_LOG") || die "\nerrro\n";
foreach my $tax(@mps)
	{
	if($record_taxon_constraints_ouput_MPa{$tax}==1 && $record_taxon_constraints_ouput_MPb{$tax}==1)
		{
		print OUT57 "backbone monophyly $tax\tapplied to new matrix \n";
		}else{
		print OUT57 "backbone monophyly $tax\tnot applied to new matrix \n";
		}
	}
close OUT57;



#unless($outgroup_found == 1){die "\nerror assigning outgroup constraint.\n"}


close OUT5;

print "out of $#new_terminals new seqs, $assigned22 can be unambiguously placed on the backbone tree, 
$not_assigned22 cannot

\n";

if($no_lineage_found >= 1)
	{
	print "no lineage found for $no_lineage_found, therefor also cant assign\n";
	}else{
	print "lineage was found for all ... something\n ";
	};



%store_node_assigned_newIDs;

$newick5 = "($root_node)";

#############################################
get_terminals_from_this_node30($root_node);#
#############################################

$newick5 =~ s/^\((.+)\)$/$1/;

$raxml_constraint_tree = "$treefile.rax_constraint";

open(NEWICK5 , ">$outfile_prefix.raxml_format_constraint_tree") || die "\nerror 193.\n";
print NEWICK5 "$newick5;\n";
#print "$newick5\n";
close NEWICK5;




# tnt: given up trying to get it running.

# open(TNT_FILE_HANDLE , ">$fas_file.$treefile.tnt") || die "\nerror 366\n";

%assigned_already;

# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string

$tnt_command_string = "tnt64 'mxram 2000, rseed 1,p tnt_input.nex timeout 0:20:00, echo=, taxname=, ";


# for fasttree constraint file (uses format similar to phylip)
$constraint_number;
%phylip_constraints;

# will be stored in hash under keys 1..$constraint_number
# ids you can get from my @all_new_IDs = keys %query_IDs;


# this works out tnt commands. 

#############################################
get_terminals_from_this_node31($root_node);#
#############################################



# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1



my @all_new_IDs = keys %query_IDs;@all_new_IDs = sort @all_new_IDs;

print "\nmaking constraint file for fasttree.
dimensions: taxa:" , scalar @all_new_IDs,
"
constraints:$constraint_number
\n";

# open(OUT_FT , ">$output_filename.ft") || die "\nerror 421\n";
print OUT_FT " " , scalar @all_new_IDs , " $constraint_number\n";

foreach my $id(@all_new_IDs)
	{
	print OUT_FT "$id        ";
		
#	for my $constaint_no(1..$constraint_number)
#		{
		my $state = $phylip_constraints{$id};print OUT_FT "$state";
#		}

	print OUT_FT "\n";
	};




#constrain =, outgroup 0, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;

$tnt_command_string .= "constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;'";



#print "\nrunnning tnt\n";
#system($tnt_command_string);

print TNT_FILE_HANDLE "$tnt_command_string";
close TNT_FILE_HANDLE;

print "
printed tnt
";



#$new_newick_string =~ s/\)(\w+)/ $assign_these_new_terminals_to_taxon{$1} . ")" . $1 /ge;

		if($new_newick_string =~ /(........Pierini.......)/)
			{print "dolar 1:$1\n"};


my @taxa_to_which_new_species_are_assigned = keys %assign_these_new_terminals_to_taxon;
@taxa_to_which_new_species_are_assigned = sort @taxa_to_which_new_species_are_assigned;
print "\n" , scalar @taxa_to_which_new_species_are_assigned , 
" taxon names on the backbone tree, to which new terminals will be added\n\n";

print  "new newick string:$new_newick_string;\n";


foreach my $taxon(@taxa_to_which_new_species_are_assigned)
	{
	if($taxon =~ /\w/)
	{
	my $new_terminals_for_this = $assign_these_new_terminals_to_taxon{$taxon};
	$new_terminals_for_this =~ s/\,$//;

	# hack to cut down tree for vewiing:
	$cut_down_for_viewing = 1;
	if($cut_down_for_viewing ==1)
		{
		$new_terminals_for_this =~ s/^([\_A-Za-z]+)\,.+\,([\_A-Za-z]+)$/$1 . "," . $2/e;#################
		};
	print "\ntax:($taxon)\tnew terminals:($new_terminals_for_this)\n";

	if($new_newick_string =~ /(.......[\(\)\.]$taxon[\(\)\.]......)/)
		{print "looks like:$1\n"};

	my $taxname = $taxon;
	$inculde_tax_branch_labels = 0;
unless($inculde_tax_branch_labels ==1){$taxname = ""};

	if($new_newick_string =~ s/\(\)$taxon([\,\(\)])/ "(" . $new_terminals_for_this . ")" . $taxname . $1/e)
		{
			#if($new_newick_string =~ /(........Pierini.......)/)
			#	{print "dolar 1:$1\n"};
		}elsif($new_newick_string =~ s/(\w)\)$taxon([\,\(\)])/ $1 . "," . $new_terminals_for_this . ")" . $taxname . $2/e)
			{
			}elsif($new_newick_string =~ s/\)\)$taxon([\,\(\)])/ ")," .  $new_terminals_for_this . ")" . $taxname . $1/e){#....tata))Cucujiformia)Polyp....
			}else{
#print  "\nERROR. newick string:$new_newick_string;\n";
		
		die "\nerror 266:cant find $taxon on tree\n"
			
			};

unless($inculde_tax_branch_labels ==1)
	{
#ta)))Nymphalidae,((Ant

while($new_newick_string =~ s/([\,\(\)])$taxon([\,\(\)])/$1$2/g){};

	};



	}};





open(NNS, ">$outfile_prefix.new_newick_string4");
print NNS "$new_newick_string;\n";
close NNS;

# this doesnt make sense:
# record_tree3("new_newick_string2");#
#########################





$new_newick_string3 = "($root_node3)";

###############################################
get_terminals_from_this_node3($root_node3 , 0);#
###############################################


open(NNS3, ">$outfile_prefix.new_newick_string3");
print NNS3 "$new_newick_string3;\n";
close NNS3;






die;


my @test_query_IDs = keys %query_IDs;
foreach my $id(@test_query_IDs)
	{
	#print "query:$id\n";
	unless (exists($terminals{$id}))
		{print "\n\nWARNING:member found in fasta file ($id) is absent from tree. are you sure you are using the correct pair of files?\n\n"}
	}


$terminals_belonging_to_current_node;




print "\n\n\nend of script\n";
exit;





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub record_tree2
{


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at tree:$treefile\n";

$tree1 =~ s/ //g;

$second_backbone_copy_newick = $tree1;
$third_backbone_copy_newick = $tree1;

# this one will have internal taxonomic node labesl:
$forth_backbone_copy_newick = $tree1;

# thus , rm supprto node labesl:
$forth_backbone_copy_newick =~ s/(\))[0-9\.]+/$1/g;



$tree_parse = 0; 	# assuming sensible tree.
if($tree_parse ==1)
	{
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	$tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 
	# remove regular branchlengths: 0.02048
	$tree1 =~ s/\:\-*\d+\.\d+//g; 
	# remove 0 length branchlengths
	$tree1 =~ s/\:\d+//g;
	# whole number node support
	# estris)100)100,UCE_v2_Bomb    &   ris))100,UC
	while($tree1 =~ s/(\))100([\)\,])/$1$2/){};
	# 1 or 2 digits
	while($tree1 =~ s/(\))[0-9]{1,2}([\)\,])/$1$2/){};
	};


my $newick_string 	= $tree1;
my $interal_node	= 0;
$backbone_copy_newick = $tree1;






#####################################################################################

# new newick parser ....

#while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
while ($newick_string =~ s/\(([^\(\)]+)\)([0-9\.]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";


		# seems raxml constraint needs to be unrooted
	if($newick_string =~ s/^\((INTERNAL_NODE_\d+)\:[\d\.]+(\,[A-Za-z_]+\:[\d\.]+)\)(\;)/$1$3/)
		{
		$node .= $2;
		print "\nwarning, since raxml constraint needs to be unrooted, " ,
			"tried to unroot your rooted tree. may or may not work.\n\troot node now:$node\n";
		};


	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas

	if($boot =~ /\d/){$nodes{$nodeID}{support} = $boot}else{$nodes{$nodeID}{support} = 100};#print "\tnodeID:$nodeID boot:$boot\n";

	my @child_nodes = split /\,/ , $node;#print "\@child_nodes:@child_nodes\n";
	$child_counts{$nodeID} = $#child_nodes;

	if($interal_node >= 6)
		{}elsif($interal_node == 5)
		{print "........\n.......\n"}else{
		print "nodeID:$nodeID\n\t$#child_nodes child_nodes:@child_nodes\n";
open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";
print BCN_LOG "NEWICK nodeID:$nodeID\tcount_child_nodes$#child_nodes child_nodes:@child_nodes\n";
close BCN_LOG;
		
		};


	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";
		my $bl = "NA"; 	
		if($child_nodes[$i] =~ /\:(.+)/)	# branchlength found
			{$bl = $1};
		$child_nodes[$i] =~ s/\:(.+)//;
		#print "node:$child_nodes[$i]\tbl:$bl\n ";

		# optional processing of tip labels
		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($remove_accessions_from_reference_species == 1)
				{$child_nodes[$i] =~ s/^([A-Z][a-z]+_[a-z]+)_.+$/$1/}
			if($process_backbone_tree_terminal_IDs == 1)
				{$child_nodes[$i] =~ s/.+_([A-Z][a-z]+_[a-z]+)$/$1/}
			if($process_backbone_tree_terminal_IDs == 2)
				{
				unless($child_nodes[$i] =~ s/^([A-Z][a-z]+)_.+$/$1/){die "\nERROR, process_backbone_tree_terminal_IDs == 2, however cant remove species from terminal:$child_nodes[$i]\n"};
				
				}
			};


		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }
		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;

		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($terminals{$child_nodes[$i]} == 1)
				{
				if($process_backbone_tree_terminal_IDs == 2){print "\nwarning, you have set script to remove species names, this assumes no conspecifics, otherwise expect crash here.\n"};
				die "\nfatal error, found duplicate tip in backbone tree:$child_nodes[$i]\n"
				};
			$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++;
			my $genus_name = $child_nodes[$i];$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			my $complete_lineage = $complete_lineage_for_this_species{$child_nodes[$i]};
			unless($complete_lineage =~ /\w/){$complete_lineage= $complete_lineage_for_this_species{$genus_name}};
			unless($child_nodes[$i] =~ /^[\w\_]+$/)
				{
				print "\nerror 100. no terminal label:$child_nodes[$i]\n";
				# die "";
				}
			if( $complete_lineage =~ /\w\s/)
				{
				$complete_lineages_retreived++;

				#print "child_nodes[i]:$child_nodes[$i] complete_lineage:$complete_lineage\n";
				# need to count how many terminals there are for each taxa,
				# for a given taxa, if all of them are decended from one node
				# it means they are monophlyetic, and will be constrained.

			#	print "\n$child_nodes[$i]\n";
				while($complete_lineage =~ s/^([^:]+):(\w+)//)
					{
					my $current_rank = $1;my $current_taxname = $2;
				#	print "\tcurrent_rank:($current_rank) current_taxname:($current_taxname)\n";
				#	current_rank:( family) current_taxname:(Bothrideridae)
					$count_taxa_represented_in_terminals{$current_taxname}++;
				#	$tax_is_of_this_rank
					}

				}else{
			#	print "\$terminal, child $i:$child_nodes[$i], genus_name:($genus_name) complete_lineage:$complete_lineage\n";
				print "\nERROR (1195). sub record_tree2. no complete lineage retrieved for genus $genus_name\n" , 
				"\t$complete_lineages_retreived have sucessfully been retrvied previously\n" , 
				"\ttip: check you are not using outdated NCBI taxonomy database\n";
				print "check its in there with command grep \"$genus_name\" names.dmp\n";
				if($genus_name =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/)
					{
					}else{
				#	die ""
					};
				# die "\tquitting.\n\n";
				};



			}#unless($child_nodes[$i] =~ /INTERNAL_NODE_/)

		}

	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";
	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
	$interal_node++;$print_internal_nodes++;
#	if(length($newick_string)<= 500){print "newick_string:$newick_string\n"};

	}; # while ($newick_string =~ s

#####################################################################################

 print "remaining of newick_string :$newick_string\n";


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




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
print "sucessfully opened fasta file $fas_file\n";


open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";
print BCN_LOG "sucessfully opened fasta file $fas_file\n";
# close BCN_LOG;


my $count_diet_IDs = 0;
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1; $fas_ID_Print++; if($fas_ID_Print <= 5){print BCN_LOG "native fasta ID format:$id\n"};
		if($id =~ /^[a-z]+_[a-z]+$/i)
			{
			$count_binomials++;
			}elsif($id =~ /^[a-z]+_[a-z]+_[a-z]+/i)
				{
				$extended_ids++;
			#	die "\nfasta ID:$id. not writenn for non-binomials \n";
				}elsif($id =~ /^[a-z]+$/i){
				$genus_names_only++;
				};
		$index_of_fastaID{$id} = $count_diet_IDs;
		$count_diet_IDs++;
		$query_IDs{$id}= 1;
		$fasta_IDs_again{$id}= 1; # print "stored $id\n";

		# get all tax names of this species, then record +1 for each.
		my $complete_lineage = $complete_lineage_for_this_species{$id};
		unless($complete_lineage =~ /\:\w/)
			{
			my $genusname = $id;

			$genusname =~ s/AntiopalaX/Antiopala/;
			$genusname =~ s/BareaX/Barea/; # BareaX_ectadia_342779467
			$genusname =~ s/^([A-Z][a-z]+)BIOUG\d+.+/$1 . "_"/e; # ArchaeognathaBIOUG07747D02)
			$genusname =~ s/^([A-Z][a-z]+)BOLD\w+\d+.+/$1 . "_"/e; # ArchaeognathaBOLDAAG6186
			$genusname =~ s/^([A-Z][a-z]+)CCDB\d+.+/$1 . "_"/e; # BlattodeaCCDB23299D03
			$genusname =~ s/^([A-Z][a-z]+)DJM2.+/$1 . "_"/e; # BlattodeaDJM2individual3


			unless($genusname =~ s/^([A-Z][a-z]+)_.*/$1/)
				{
			#	print "\nerror trying to extrat genus name from ($id)\n";
				};
		#	print "\n no lineage found for id ($id)\ntrying genus name ($genusname)\n";

			$complete_lineage = $complete_lineage_for_this_species{$genusname};
			unless($complete_lineage =~ /\:\w/)
				{
			#	print "warning ... cant find lineage for sequence id:$id genusname:$genusname\n";
				# as of SEP2016, will try without constraints for some unlabeled BOLD data
				# i note mispellings slip though here and will be unconstrained
				# but these are mosrly the order labelled BOLD seqs
				};
			}# unless($complete_lineage =~ /\:\w/)
		if($print_ID_count <= 10)
			{
			print "fasta ID:$id\n";
			}elsif($print_ID_count == 11){print ".......\n"};
			$print_ID_count++;

		#print "complete_lineage:$complete_lineage\n";

		# sometimes identical string (name) used for different ranks
		my %tax_used = ();

		while ($complete_lineage =~ s/ (\w+)\:(\w+) / /)
			{
			my $rank = $1; my $tax = $2;$count_tax_in_new_data{$tax}++; #print "tax:$tax\n";
			if($rank eq "family"){
			$family_of_fasta_ID{$id} = $tax};

			unless($tax_used{$tax} == 1)
				{
				$fasta_ids_for_each_taxon{$tax} .= "\t$id\t";
				};	
			$tax_used{$tax} = 1;

			};		

		}else{
		if($line =~ /^>/){die "\nerror 646.\n"}
		}
	}
close(IN);

print BCN_LOG "
count_binomials:$count_binomials 
extended_ids (Order_Genus_species):$extended_ids
genus_names_only:$genus_names_only
";






$autodetect_fasta_IDs = 0;
if($count_diet_IDs == $count_binomials)
	{$autodetect_fasta_IDs = "binomial";print "\nautodetected strict binomial used in your fasta file. good!\n"};
if($count_diet_IDs == $extended_ids)
	{$autodetect_fasta_IDs = "extended";print "\nwarning. you are not using strict binomials in your fasta file, not guarunteed to work\n"};

my @check_species_filtered = keys %query_IDs;
# print "check_species_filtered:@check_species_filtered\n";die "";

if($count_diet_IDs == scalar @check_species_filtered)
	{
	print "\npass check, all strings in fasta file are unique (species filtered)\n"
	}else{
	print "\nwarning. duplicate IDs found in your fasta file.\n"
	}

print "\n*** $count_diet_IDs fasta IDs *** in your fasta file ($fas_file). will be looking for these in the tree ...\n\n";

unless($count_diet_IDs >= 1)
	{die "\nerror ... no fasta ID's found in the file ($fas_file). quitting.\n"}



close BCN_LOG;


}#sub read_fas





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_terminals_from_this_node
{
my $next = $_[0];my $sum_branchlength = $_[1];


#print "\nNEW NODE ($next)";

my $child_nodes = $child_counts{$next};
my @next1 = ();
my @non_terminal_child_nodes = ();
my @child_nodes_test = ();


my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();
my $tax_for_all_child_nodes = "";
my $all_child_nodes_are_internal =1;

my $tax_for_parent = "";
					# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{
	$all_terminals_from_this_node ="";
	get_terminals_from_this_node2($next, 0);
	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);
	#print "current node ($next) is internal, shared tax name ($tax_for_parent)\n";
	};


for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	$all_terminals_from_this_node ="";# for current child node only

	if($nodes{$next}{$i} =~ /INTERNAL_NODE_\d/)
		{
		# if current child is internal node, go and get the decending terminals
			# this sub does this for terminal nodes only:
			# $all_terminals_from_this_node .= "$test\t";
			# gets binomial names
		#print "line561. getting terminal binomials for child node $nodes{$next}{$i}\n";
		get_terminals_from_this_node2($nodes{$next}{$i}, 0);
		
		}else{
		# otherwise, just append the node name itself
		my $established_name = $nodes{$next}{$i};
		$all_terminals_from_this_node .= "$established_name\t";
		$all_child_nodes_are_internal =0
		};



		$all_terminals_from_both .= $all_terminals_from_this_node;
		$all_terminals_from_this_node =~ s/\t$//;

				##############################################################
		my $tax = 	get_shared_taxonomic_name($all_terminals_from_this_node);#
				##############################################################

		push @taxon_for_child_nodes, $tax;
		#print "child $i of node $next has ID $nodes{$next}{$i}. $count_terminal_nodestermianls got, shared tax:$tax.\n";

		# here assess if child node is same taxonomy as node.
		# if so, this node is collapsed.
		# this is done by stepping ahead to getting the two grandchild nodes, 
		# then write both these in instead of the child node.


		if($tax eq $tax_for_parent)
			{
			print "this child has same taxon ($tax,$tax_for_parent), collapsing.\n";
			my $grandchild_nodes = $child_counts{$nodes{$next}{$i}};
			for $j(0 .. $grandchild_nodes)
				{
				my $push_node = $nodes{$nodes{$next}{$i}}{$j};
				print "push grandchild to array j:$j, push_node:$push_node\n";
				push @next1 , $push_node;
				if($push_node =~ /INTERNAL_NODE_\d/)	# for makeing the newick string with termainals for new data,
					{
					push @non_terminal_child_nodes, $push_node;
					push @child_nodes_test, $push_node;
					}elsif# first need to remove existing terminals
					($push_node =~ /\w/)
						{
						push @child_nodes_test, $push_node;
						}else{die "\n\nerror 788\n"}	
				};
			}else{
			print "this child has diff taxon ($tax,$tax_for_parent)\n";
			my $push_node = $nodes{$next}{$i};
			push @next1 , $push_node;	
			if($push_node =~ /INTERNAL_NODE_\d/)
				{push @non_terminal_child_nodes, $push_node};
			push @child_nodes_test, $push_node;
			};

		if($i == $child_nodes)
			{
			$all_terminals_from_both=~ s/\t$//;
			print "line793, calling sub get_shared_taxonomic_name\n";
			$tax_for_all_child_nodes = get_shared_taxonomic_name($all_terminals_from_both);
			unless($tax =~ /\w/){die "\nno tax name retreived for this set of terminals\n"}
			};



	}



	# this is the default way to building up the newick string with new nodes.
	# it is overwritten later in situations where the node needs to be collapsed.


# this is collapsed where appropriate, but includes terminal nodes. 
# printed to an additional newick string so user can see better how tree has been collapsed 
my $join_child_nodes_test = join ',', @child_nodes_test;

# removed terminals:
my $join_child_nodes = join ',',@non_terminal_child_nodes;
my $swap_string ="";

#if($#non_terminal_child_nodes <= 0)
#	{
#	print "no decendents. dead end.\n";
##$swap_string = "$tax_for_parent";
#$swap_string = "($join_child_nodes)$tax_for_parent";
#	}else{
$swap_string = "($join_child_nodes)$tax_for_parent";
my $swap_string_test = "($join_child_nodes_test)$tax_for_parent";
#	}



#print "tax_for_parent:$tax_for_parent\n";
#die;
$all_internal_taxa_on_constraint_tree{$tax_for_parent}=1;


if($next =~ /INTERNAL_NODE_\d/ && $all_child_nodes_are_internal == 1)
	{
#	$all_terminals_from_this_node ="";
#	get_terminals_from_this_node2($next, 0);
#	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);


	my $all_branches_from_this_node_are_the_same_taxon = 1;
	foreach my $test_taxon(@taxon_for_child_nodes)
		{
		if($tax_for_parent eq $test_taxon)
			{}else{$all_branches_from_this_node_are_the_same_taxon=0}
		};


	# collapse node method 1
	# turns not really what i want to do here.
	# but i keep the code in case it is what i want to do in a future situation
	if($all_branches_from_this_node_are_the_same_taxon ==1)
		{
		print "node $next will be COLLAPSED. tax_for_parent:$tax_for_parent taxon_for_child_nodes:@taxon_for_child_nodes\n";
	#	$swap_string = "$join_child_nodes";
		# collapse node. get identity of parent of current node (which is itself parent)
		# then add the child nodes to this  
		# no. instead you need to get the grandchild node of the matching taxon branch, and pull it up to be a child node.

		#$child_counts{$nodeID} = $#child_nodes;
		#$nodes{$child_nodes[$i]}{parent} 	= $nodeID;

		}


	};

my $join_the_child_nodes;
$collapse_nodes = 0;
if($collapse_nodes == 1)
	{
	# build the newick string, replace current internal node with child nodes
	print "node:$next becomes ($swap_string)\n";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;
	$newick_print_string=~ s/$next(\W)/$swap_string_test$1/;
	}else{
	@next1 = ();
	for $i2(0 .. $child_nodes)
		{push @next1, $nodes{$next}{$i2}};
	# default:
	$join_the_child_nodes = join ',', @next1;
	$swap_string = "($join_the_child_nodes)$tax_for_parent";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;
	$newick_print_string=~ s/$next(\W)/$swap_string$1/;

	};

# collapsing is acheived by not looping to child nodes, but skipping strait to grandchild nodes.
# similarly, grandchild nodes are printed to the growing newick string instead of child nodes.


for my $index(0 .. $#next1)
	{
############	my $test = @next1[$index];#print "test:$test\n";####################################### surely this is wrong?
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################
		}
	}

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub get_shared_taxa
{
my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $shared_tax_substring = "";


# gets lineage (series of taxa) shared by all members. also option to get only lowest shared taxa



my $index=0;

my $test_species = $tax_array[0];

	# hash complete_lineage_for_this_species, has not just species, 
	# but has lineage for every taxon in ncbi, only modification is space->underscore

my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
#print "line928. test_species:$test_species test_taxonomy:$test_taxonomy\n";


unless($test_taxonomy =~ /\w+/)
	{
	my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
	$test_taxonomy = $complete_lineage_for_this_species{$genusname};
	print "cant find taxonomy for species ($test_species) trying genus ($genusname)\n";
	if($test_taxonomy =~ /\w+/)
		{
		print "\tphew, found using genus name ($genusname), continuing ...\n";
		}else{
		if($test_species =~ /optera_[A-Z]/)
			{print "\nlooks like user format was not processed to an established tax name\n"};
		die "\nerror 1701, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
		}
	}

my $most_inclusive_name = "NA";
my $most_inclusive_lineage = "";


while($test_taxonomy =~ s/^([^:]+):(\w+)//)
	{
	my $current_rank = $1;my $current_taxname = $2;
	# over counts, should not count taxa at parent node:
	#$number_of_nodes_to_which_this_taxon_has_been_assigned{$current_taxname}++;
	#print "\tcurrent_rank:$current_rank current_taxname:$current_taxname\n";

	if($#tax_array >= 1)
		{
		# this is the regular situation, there is more than one terminal you have got a lineage for,
		# and you are finding the shared taxa for these

		my $all_members_have_this_name =1;
		foreach my $tax(@tax_array)
			{
			my $test_taxonomy2 = $complete_lineage_for_this_species{$tax};
			my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			unless($test_taxonomy2 =~ /\w/){$test_taxonomy2 = $complete_lineage_for_this_species{$genus_name}};
			#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";
			if($test_taxonomy2 =~ /\:$current_taxname\s/)
				{print ""}else{
				$all_members_have_this_name =0;#print "0";
				}
			}

		#print "\nall_members_have_this_name:$all_members_have_this_name\n";

		if($all_members_have_this_name == 1)
			{
			$most_inclusive_name = $current_taxname;
			$most_inclusive_lineage .= "$current_taxname\t";
			};

		}else{

		# alternativly, you only have one taxa, so just parse the whole lineage for this		

		$most_inclusive_name = $current_taxname;
		$most_inclusive_lineage .= "$current_taxname\t";

		}

	#print "\tmost_inclusive_lineage:$most_inclusive_lineage\n";
	}


#if($tax_list eq "Apis_florea"){
#print "tax_list:($tax_list) tax_array:(@tax_array) test_taxonomy:($test_taxonomy)\n";
#die;
#}


# print name only:
#return($most_inclusive_name);
# or lineage:
#print "\treturning:($most_inclusive_lineage)\n";
return($most_inclusive_lineage);

}





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

print "
\n\n
     ***** backbone_constraints_newick.pl *****
		  
\n";




if($arguments =~ /-treefile\s+(\S+)/)
	{
	$treefile = $1;
	}else{
	print "\nerror reading command arguments  (-treefile)\n\n";die""
	}

if($arguments =~ /-outfile_prefix\s+(\S+)/)
	{
	$outfile_prefix = $1;
	}else{
	print "\nerror reading command arguments  (-outfile_prefix)\n\n";die""
	}


if($arguments =~ /-support\s+(\S+)/)
	{
	$support_cutoff = $1;
	}else{
	$support_cutoff = 50;
	print "user did not give cutoff for acceptable boot support (not relevent if you are not using bootstrapped trees)\nusing default ($support_cutoff).\n"
	}

if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	}

if($arguments =~ /-references\s+(\S+)/)
	{
	$reference_file = $1;
	}else{
	}

if($arguments =~ /-backbone_terminal_format\s+(\d)/)
	{
	$process_backbone_tree_terminal_IDs = $1;

# $process_backbone_tree_terminal_IDs = 2;


	}else{
	die "\n\ncommand error, please give format you have used for terminal labels of backbone tree\n" , 
	"  -backbone_terminal_format [0,1,2]\n" , 
	"  if == 0, then just genus names\n", 
	"  if == 1 (unlikely), then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis\n", 
	"  if == 2, then Genus_speciesID (which is converted to genus)\n\n";

	}


if($arguments =~ /-outgroup\s+(\S+)/)
	{
	$outgroup = $1;
	unless($outgroup =~ /[A-Z][a-z]+/){die "\noutgroup switch used but cannot parse name \n"}
	}else{
	print "\nyou did not specify the outgroup.\n";
	}


if($arguments =~ /-constrain_ranks\s+(\S[^\-]+)/)
	{
	$constrain_ranks1 = $1;
	@constrain_ranks = split /\s+/ , $constrain_ranks1;
	print "\nuser specified ranks at which to make constraints:\n@constrain_ranks\n";
	}else{
	@constrain_ranks = ("order","suborder","infraorder","superfamily","family","subfamily","genus");
	print "\nyou did not specify ranks at which to make constraints, applying default:\n@constrain_ranks\n";
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

		my $current_rankcode = 0;
		foreach my $rank_test("order","suborder","infraorder","series","superfamily",
			"family","subfamily","tribe","subtribe","genus",
			"subgenus","species group","species","subspecies")
				{
				$current_rankcode++;
				if($rank eq $rank_test){$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode}
				};


			$ncbi_nodes{$tax_id}{rank} = $rank;
			
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

		$name =~ s/^Candidatus\s(\w+)$/$1/;

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
	};



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;

	
	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name


	# NOV 2018: user decides which ranks are constrained	
	my $constrain_this_node = 0; # print "\n\n";
	foreach my $user_constrain_rank ( @constrain_ranks )
		{
		if($rank eq $user_constrain_rank)
			{$constrain_this_node = 1};
	#	print "rank:$rank user:$user_constrain_rank constrain:$constrain_this_node\n";
		};
	if($constrain_this_node == 1){$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage};

# forgot about this, explains some imperfect species level trees
#	unless($rank eq "order") # 31 aug 2016: way to stop seqs only given order ident from being constrained.
#		{
#		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};


	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	$ncbi_taxnumber_for_taxname{$originalname} = $child;
#	print "ncbi_number:$child tax name:$originalname\n";

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




sub get_terminals_from_this_node2
{
my $next = shift;

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	
	#print "line1423node:$next i:$i child:$nodes{$next}{$i}\n";
	}



for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test 1430:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node2($test );#	# recurse
		#####################################

		}else{
		$count_terminal_nodes++;

		if($test eq "")
			{
			print "\nnext:($next)\n";
			die "\nerror 1443. looks like sub get_terminals_from_this_node2 was called with a terminal node.\n"
			}


		### here set the format used in the previous (mtgenome) tree.
		# currently im using Order_Genus in the tree, and parsing the genus name.
#		unless($test =~ s/.+_//)
#			{die "strange format:($test)\n"}

		# 8mar, Coleoptera_Tetraphalerus_bruchi and parsing genus
	#	unless($test =~ s/.+_(.+_.+)/$1/)
	#		{die "error 1229. expecting Order_Genus_species, got strange format:($test)\n"}
		# now process ids in intial reading of the tree, so they just have a regular taxon name

		$all_terminals_from_this_node .= "$test\t";

		#print "line1459. test:$test\t";
		}
	}

return();

}























sub record_tree3  # sub not used
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


sub get_terminals_from_this_node3
{
my $next = $_[0];my $sum_branchlength = $_[1];

my $child_nodes = $child_counts3{$next};
my @next1 = ();

for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	my $push_node = $nodes3{$next}{$i};
	push @next1 , $push_node;	
	}

# default:
my $join_child_nodes = join ',', @next1;
$swap_string = "($join_child_nodes)$tax_for_parent";
$new_newick_string3 =~ s/$next(\W)/$swap_string$1/;

for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node3($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################
		}
	}

return();

}








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub traverse_backbone_tree_and_infer_taxonomic_node_labels
{

my $next = $_[0];my $sum_branchlength = $_[1];


print "\nNEW NODE ($next)\n";

my @non_terminal_child_nodes = ();
my @child_nodes_test = ();
my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();

my $lowest_taxname_for_node = "";
my $lineage_for_node = "";
$all_terminals_from_this_node ="";
$taxa_for_terminals_of_node = "";

# need a list of taxa derived from currnet node. ($all_terminals_from_this_node .= "$test\t";)
# get tax names for these.
# then go through all remaining members of the tree, (%terminals, not incl $all_terminals_from_this_node)
# if none of these have these tax names (get_shared_taxa), they are defined as a constraint.
# then all new members from these, are constrained.			


		# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{

	# traversing backbone tree, sub run for current node.
	# get members (terminals) of this node

	##########################################
	get_terminals_from_this_node2($next, 0);#
	##########################################


	# following gets lineage (series of taxa) shared by all members of clade. 
	# also option to get only lowest shared taxa

				####################################################
	$lineage_for_node = 	get_shared_taxa($all_terminals_from_this_node);#
				####################################################
	unless($lineage_for_node =~ /\w/){die "error 1944 could not parse lineage for node ($next)\n"};
	#print "get_all_taxa 1\n";

					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($all_terminals_from_this_node);#
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/){die "error 1950 could not parse taxa for terminals of node ($next)\n"};


	#########################################################
	my $monophyly_name = test_for_monophylies($all_terminals_from_this_node);#
	#########################################################

	# if monophyly found, record the taxonomic name on the backbone node.
	unless($monophyly_name eq "NA")
		{$node_monophylys{$next} = $monophyly_name};

	}else{#if($next =~ /INTERNAL_NODE_\d/)

	# for tips, you still need to know the lineage. but it wont be a shared lineage.
	# in this case you dont have a list of terminals, theres just one, the curent node

				#########################
	$lineage_for_node = 	get_shared_taxa($next);#
				#########################
	unless($lineage_for_node =~ /\w/){die "error 1963 could not parse lineage for node ($next)\n"};

	#print "get_all_taxa 2\n";
					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($next);# bug fix 22 oct 2015
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/){die "error 1969 .could not parse taxa for terminals of node ($next)\n"};


	};

# parse all names, put into this object
my $lineage_string = $lineage_for_node;
while($lineage_string =~ s/([^\t]+)\t$//)
	{
	my $taxon_for_node = $1;$all_taxa_on_backbone_tree{$taxon_for_node}=1;
	#print "taxname1:$taxon_for_node lineage_string:$lineage_string\n";
	};


# for current node, parse lowest shared taxa from the lineage:
$lowest_taxname_for_node = $lineage_for_node;
unless($lowest_taxname_for_node =~ s/.+\t([^\t]+)\t$/$1/)
	{
	unless($lowest_taxname_for_node =~ s/^$starting_name_ff	/$starting_name_ff/ )# basal TOL will have shared taxon being cellular orgs, same as root specified by user
		{
		print "next:$next\nlineage_for_node:$lineage_for_node\ntaxa_for_terminals_of_node:$taxa_for_terminals_of_node\nerror 1982\n";
		die "\neror 1550, lowest_taxname_for_node:($lowest_taxname_for_node)\n"
		};
	};

# check for error:
unless($lowest_taxname_for_node =~ /^[A-Za-z\_]+$/)
	{die "\nerror 1638. lowest shared taxa for current node ($next) has unexpected format:($lowest_taxname_for_node).\n"}
# and store the LOWEST name in this hash (no intermediates)
$node_label{$next} = $lowest_taxname_for_node;

print "lowest_taxname_for_node:$lowest_taxname_for_node\n";

$all_taxnames_for_terminals_decended_from_this_node{$next} = $taxa_for_terminals_of_node;

# 
# next compare the complete shared lineage of the current node, with that of the parent.
# in the example below, the node in which all decendents are Lucanidae,
# in addition to being the position at which new Lucanids are added
# can also accomodate the above ranks Scarabaeoidea, Polyphaga,
# but not Coleoptera, which is best assigned to the parent node.
# therefor find (if present) all intermediate ranks.
#        
#                    -------
#                    |
#           -----Lucanidae
#           |        |
#           |        ------
# ---Coleoptera      
#           |        -----
#           |        |
#           -----Adephaga
#                    |
#                    -----


my $parent_node = $nodes{$next}{parent};
my $parent_taxon = "";
my $all_potential_taxa ="";

if($node_label{$parent_node} =~ /\w/)
	{
	$parent_taxon = $node_label{$parent_node};
	#print "node($next), tax($lowest_taxname_for_node) parent($parent_node) Tax($parent_taxon)\n";


	if($lowest_taxname_for_node eq $parent_taxon)	# compare the lowest tax name of current node
		{						# to lowest tax name of parent node

		# Node and parent both have the same assigned taxon, nothing to do.
		# Although whilst here, need to count the number of nodes each taxon has been assigned to 
		$number_of_nodes_to_which_this_taxon_has_been_assigned{$lowest_taxname_for_node}++;

		}else{
		while($lineage_for_node =~ s/([^\t]+)\t$//)
			{
			# go through whole lineage of current node 
			# (including lowest, which is not trimmed from this particular variable), 
			# and stop when reaching the taxon assigned to parent
			# this direction:Corcyra_cephalonica->Corcyra->Galleriinae->Pyralidae

			my $intermediate_taxon_for_node = $1;
			unless($intermediate_taxon_for_node =~ /^[A-Za-z\_]+$/){die "\nerror 1691. unexpected taxon name:$intermediate_taxon_for_node\n"}
			#print "\ttaxname1:$intermediate_taxon_for_node. comparing to that assigned to parent ($parent_taxon)\n";

			# the first element is the lowest taxon of current node
			# if that is the same as that of the parent node, 
			# there will be nothing put in $all_potential_taxa
			# if there is only one item in $all_potential_taxa, it will be same as lowest taxon
			if($intermediate_taxon_for_node eq $parent_taxon)
				{
				#print "same, break.\n";
				$lineage_for_node = "";# delete lineage string to break loop.
				$all_potential_taxa =~ s/[\s\.]+$//;# if names have been assigned, there will be terminal char
				}else{

				# taxon assigned to parent node has not been reached yet, so append taxon name 
				$all_potential_taxa .= "$intermediate_taxon_for_node.";#print "appending\n";

				# And whilst here, keep counting the number of nodes each taxon has been assigned to 
				$number_of_nodes_to_which_this_taxon_has_been_assigned{$intermediate_taxon_for_node}++;
				};
			};
		}

	}else{#if($node_label{$parent_node} =~ /\w/)
	
	unless($next eq $root_node)# only the root node wont have a parent, so no taxon will be found for that one.
		{
		die "\nerror 1703. no taxon assigned to parent ($parent_node) of node ($next)\n";
		};
	}


my $label_string_to_assign = "";
if($all_potential_taxa =~ /\w[\.\s]\w/)# taxa found which are intermediate between those assigned to node and parent
	{
	$label_string_to_assign = $all_potential_taxa	# $all_potential_taxa has fullstop seperated taxon names
	}else{							# so they can be printed by tree viewers.

	$label_string_to_assign = $lowest_taxname_for_node # otherwise, just use lowest shared taxon name
	};
unless($label_string_to_assign =~ /\w\w+/){die "\nwhy no taxon assigned to current node ($next). \n"}


my $child_nodes = $child_counts{$next};
my @next1 = ();

# record the taxonomic name(s) given to this node
$nodes{$next}{node_taxonomic_label}=$label_string_to_assign;

#print "\n\nNEW BACKBONE NODE ($next) assigned name(s):$label_string_to_assign\n";







for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}

#print "node:($next) child nodes:(@next1) label_string_to_assign:($label_string_to_assign) lineage_for_node:($lineage_for_node)\n";

#if($next eq "Apis_florea"){die""};

# get to the tip,
#NEW NODE (Apis_florea)
#node:(Apis_florea) child nodes:() label_string_to_assign:()

# not used any more ?:
# build newick string. this string is an intermediate, 
# it contains empty parentheses ready to be filled with the species-dense data,
#thus not readable by tree viewers

my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;
$swap_string 		= "($join_the_child_nodes)$label_string_to_assign";
$new_newick_string 	=~ s/$next(\W)/$swap_string$1/;

# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


$all_internal_taxa_on_constraint_tree{$tax_for_parent}	= 1;

# the approach in this sub, is called all the way up to tip nodes, so tax lineages etc can be easier assessed.
# which means it needs to determine the end has been reached below (no child nodes for this tip)
# its a bit different from the usual approach, in which the sub is not called for tip nodes.

if(exists($child_counts{$next}))
	{

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_infer_taxonomic_node_labels($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################

		}

	};

return();

}






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


sub get_terminals_from_this_node30
{
my $next = $_[0];
my $child_nodes = $child_counts{$next};

#print "\nNEW NODE:($next)\n";

my @next1 = ();


my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

for $i(0 .. $child_nodes)
	{


	# get child nodes as normal
	my $push_node = $nodes{$next}{$i};

	unless(exists($query_IDs{$push_node}))
		{
		unless($push_node =~ /INTERNAL_NODE_/)
			{
			$backbone_child_nodes_absent_from_new_IDs++;
			print "warning 1875. node $push_node is not internal, yet is not found in list of new IDs\n";
			};
		}

	push @next1 , $push_node;	
	#print "child node:$i node ID:$push_node\n";

	my $node_taxonomic_label = $nodes{$push_node}{node_taxonomic_label};

	my @split_node_tax_labels = split /\./ , $node_taxonomic_label;
	foreach my $candidate_tax(@split_node_tax_labels)
		{
		#print "\tchild node is assigned candidate_tax:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

		unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

		if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
			{
			my $string_length_new = length($assign_these_new_terminals_to_taxon{$candidate_tax});

			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added, length:$string_length_new\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_members =~ s/\,$//;

			# each new member come across, count in hash.
			# later check duplicates are not met.
			my @split_new_members= split /\,/ , $new_members;
			foreach my $newID(@split_new_members)
				{
				$store_node_assigned_newIDs{$newID}++;$count_new_IDs_added_to_current_node++;
				unless($newID eq $push_node){$add_to_current_node .= $newID . ",";}
				};

			# append, this will be a new child of the current node
			#$add_to_current_node .= $new_members . ",";
			my @check_split= split /\,/ , $add_to_current_node;

			#if($add_to_current_node =~ /pisumAcyr/){die "\n1887.\n"}

			#$check1++;

			}#if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))

		};#foreach my $candidate_tax(@split_node_tax_labels)

	};

#print "
#count_new_IDs_added_to_current_node:$count_new_IDs_added_to_current_node
#backbone_child_nodes_absent_from_new_IDs:$backbone_child_nodes_absent_from_new_IDs
#";
# check any havnt been assigned elsewhere (shouldnt have)
my @all_new_IDs_assigned_at_this_point = keys %store_node_assigned_newIDs; 
foreach my $newID(@all_new_IDs_assigned_at_this_point)
	{
	if($store_node_assigned_newIDs{$newID}>1)
		{
		print "\nerror 1880. new ID:($newID), to be assigned to node:($next),\n";
		print "which has tax:(...), was assigned previously. quit.\n";die};
	#unless($newID eq ""){};# terminal comma needs keeping
				# so last object will be empty
	};




my $copy = $add_to_current_node;#$copy .= join ',', @next1;
#$copy =~ s/\,$//;
my @check_split= split /\,/ , $copy;
my %check_hash=();
foreach my $check4(@check_split)
	{
#unless($check4 =~ /^[A-Z][a-z]+_[a-z]+$/){die "\nerror 1828:$check4\n"}
	if(exists($query_IDs{$check4}))
		{
	$check_hash{$check4}++;#if($check_hash{$check4}>= 2){die "\nerror 1822, $check4 assigned already\n"}
	#$check_hash4{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}
		}
	}

my @check_split2= keys %check_hash , @next1;
my %check_hash2=();
foreach my $check4(@check_split2)
	{

	$check_hash2{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}

	}


my @new_nextnode = keys %check_hash2;

#if($check1>= 3){die "\nerror 1816\n"}


#print "node:$next child nodes:@next1 node_taxonomic_label:$node_taxonomic_label\n";
if($add_to_current_node =~ /./)
	{
	#print "new members will be added to this node, ($#check_split)\n"
	};
#if($next1[1] =~ /Trachypachus_holmbergi/i){die}

#	


$howmany+= $backbone_child_nodes_absent_from_new_IDs;

if(
$backbone_child_nodes_absent_from_new_IDs >= 1 && 
$count_new_IDs_added_to_current_node < 2
){
	# this adds to 2, so ok
unless($count_new_IDs_added_to_current_node==1 && $backbone_child_nodes_absent_from_new_IDs == 1)
	{
unless($count_new_IDs_added_to_current_node == 0 && $backbone_child_nodes_absent_from_new_IDs == 1)
{
#$howmany++;
#die "\nyou have a problem. current node ... \n";
}	}
}

# default:
my $join_child_nodes = $add_to_current_node . join ',', @next1;

$swap_string = "($join_child_nodes)";
unless($newick5 =~ s/$next(\W)/$swap_string$1/)
	{print "\nerror 1962, cant build string\n"};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node30($test  );#	# recurse
		#####################################
		}
	}

return();

}








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub get_terminals_from_this_node31
{
my $next 	= $_[0];
my $child_nodes = $child_counts{$next};


my @next1 = ();
my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

my $node_taxonomic_label = $nodes{$next}{node_taxonomic_label};

# keep this print: 
print "\n*** NEW NODE of backbone *** 
Making TNT / FT constraint no:$constraint_number; node:($next); tax:($node_taxonomic_label)\n";


 # sometimes they have multiple names a la: Diabrotica.Diabroticites.Diabroticina.Luperini.Galerucinae
my @split_node_tax_labels = split /\./ , $node_taxonomic_label;

my $new_IDs_going_to_this_node = "";
my $count_tax_assigned_to_node =0;
my $total_tax_assigned_to_node  =scalar @split_node_tax_labels;
foreach my $candidate_tax(@split_node_tax_labels)
	{
	$count_tax_assigned_to_node++;
	print "\ttax:$count_tax_assigned_to_node of total:$total_tax_assigned_to_node assigned. taxon name:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

	unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

	if(exists($assign_these_new_terminals_to_taxon_TNT{$candidate_tax}))
		{
		print "\tIS new members to assign to this name\n";

		if(exists($assigned_already{$candidate_tax}))
			{
			print "\twarning, taxon($candidate_tax) was come across on previous node. ignoring.\n";
			}else{
			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_IDs_going_to_this_node .= $new_members;
			$new_members =~ s/\,$//;
			my @new_members_to_node = split /\,/ , $new_members;
			print "\tcount new members to this node:" , scalar @new_members_to_node , "\n";
			};

		$assigned_already{$candidate_tax} = 1;
		}else{#if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
		print "\tNO new members to assign to this name\n";
		}

	};#foreach my $candidate_tax(@split_node_tax_labels)






# constrinat is made of the new taxa, and the list of terminasl in backbone tree, so get the terminals:
$all_terminals_from_this_node = "";
# does this:$all_terminals_from_this_node .= "$test\t";

#######################################
get_terminals_from_this_node2($next);#
#######################################





unless($next eq $root_node)
	{



# 3 components to each constraint: 1)the taxa monophyletic 2) the taxa excluded 3) the taxa that float.
# err on the side of floating!

my %floaters = ();
my %constrained = ();

# all backbone IDs:		%terminals
# backbone IDs to constrain:	$all_terminals_from_this_node


my @all_backbone_IDs = keys %terminals;@all_backbone_IDs = sort @all_backbone_IDs;
foreach my $backbone_ID(@all_backbone_IDs)
	{
	#print "backbone_ID:$backbone_ID\n";	
	if($all_terminals_from_this_node =~ /$backbone_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$backbone_ID}=1;# taxa in the backbone tree decended from this node will be constrained MP
		}else{
		$constrained{$backbone_ID}=2;# taxa in backbone tree not decended from this node will be excluded from MP
		#print "\texclude from MP\n";#die;
		};
	};

#$all_terminals_from_this_node

my @all_new_IDs = keys %query_IDs;
foreach my $new_ID(@all_new_IDs)
	{
	my $check_overlapping = $constrained{$new_ID};
	#print "new_ID:$new_ID check_overlapping:$check_overlapping\n";
	unless($check_overlapping =~ /\d/)
	{
	if($new_IDs_going_to_this_node =~ /$new_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$new_ID}=1;# new taxa assigned to this node will be constrained MP
		}else{
		#print "\tfloat\n";
		$constrained{$new_ID}=3;# new taxa not assigned to this node will float
		};
	}};



my @everyone = keys %constrained;@everyone = sort @everyone;

my $constrain_string="";
my $float_string="";
my %check_assigned = ();
$constraint_number++;


foreach my $member(@everyone)
	{
	#print "member:$member\n";
	my $contr = $constrained{$member};
	my $fasta_index;
	if($index_of_fastaID{$member} =~ /\d/)
		{
		$fasta_index = $index_of_fastaID{$member};
		if(exists($check_assigned{$fasta_index}))
			{die "\nerror 2218. member:$member fasta_index:$fasta_index\n"};
		$check_assigned{$fasta_index}=1;
		}else{die "\nerror 2190. constr:$contr\tmem:$member\t\n"}
	#print "$contr\t$member\tfile index:$fasta_index\n";

	if($contr == 1){$phylip_constraints{$member} .= "1";$constrain_string .= "$fasta_index "}
	if($contr == 2){$phylip_constraints{$member} .= "0"};
	if($contr == 3){$phylip_constraints{$member} .= "-";$float_string .= "$fasta_index "}
	}

$constrain_string	=~ s/\s$//;
$float_string		=~ s/\s$//;

my $tnt_command;
$tnt_command = "[ $constrain_string ($float_string) ]";

$tnt_command_string .= "force + $tnt_command, ";
#print TNT_FILE_HANDLE "$tnt_command_string";


#		$index_of_fastaID{$id} = $count_diet_IDs;
# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string
	# crashes unless you have single quotes around force command, due to parentheses confusing bash.
#	tnt64 mxram 2000, rseed 1,p tnt_input.nex timeout 0:00:30, echo=, taxname=, 'force + [ 14 46 67 (1 2) ],' constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / , export *tnt_output, log / ,quit;
#print "backbone constraint (list of all_terminals_from_this_node):$all_terminals_from_this_node\n";
#print "new_IDs_going_to_this_node:$new_IDs_going_to_this_node\n";


# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1




	};#unless($next eq $root_node)




for $i(0 .. $child_nodes)# get child nodes as normal
	{
	push @next1 , $nodes{$next}{$i};
	#print "child node:$i node ID:$nodes{$next}{$i}\n";
	};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node31($test  );#	# recurse
		#####################################
		}
	}

return();

}


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





sub traverse_backbone_tree_and_collapse
{

my $next = shift;

# this has species names, species groups, genera, families etc. and has the outgroup on first call.
my $all_taxnames_for_terminals = $all_taxnames_for_terminals_decended_from_this_node{$next};
#print "all_taxnames_for_terminals:($all_taxnames_for_terminals)\n";


# taxonomic name assigned to node of backbone tree:
my $label = $node_label{$next};# just one name in this
# or this one:
#my $label = $nodes{$next}{node_taxonomic_label};
#print "\nNEW NODE ($next), label:$label\n";


my @non_terminal_child_nodes = ();my @child_nodes_test = ();
my $all_terminals_from_both = "";my @taxon_for_child_nodes = ();
my $lowest_taxname_for_node = "";my $lineage_for_node = "";
$all_terminals_from_this_node ="";


#my @tax_labels_current_node_of_backbone_tree = split /\./ , $label;



#foreach my $taxlab(@tax_labels_current_node_of_backbone_tree)
#	{
my $taxlab = $label;
my $tax_number = $ncbi_tax_number_for_this_species{$taxlab};

# all taxa at rank below the name:
my $child_nodes = $ncbi_nodes{$tax_number}{child_nodes};$child_nodes =~ s/^\t//;



#print "
#label assigned to current node of backbone:$taxlab, NCBI tax_number:$tax_number
#tax_number:$tax_number has child_nodes:$child_nodes\n";

my @child_nodes_array = split(/\t/, $child_nodes);
my $current_constraint = "";my $current_constraintB = "";

my $mpA =0;my $mpB =0;my $data_missing_reached=0;
my $count_floaters = 0;my $count_constrained = 0;

# each child node (taxonomic name) derived from name that was given to current node of tree
foreach my $child(@child_nodes_array)
	{
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g;	
	# print "\tchild taxon ID:$child name:$name_string\n";

	# for each child taxon (1 rank below) derived from the 
	# taxonomic label assigned to the current node of the backbone tree.

	my $found_in_backbone_decendent=0;
	if($all_taxnames_for_terminals =~ /\t$name_string\t/)
		{
		# child taxon is found in the lineage of terminals decended from this backbone node			
		# so a monophyletic constraint may be made			
		$found_in_backbone_decendent =1;#print "\tIS in list of taxnames of decendents\n";

		#print "\t\t\ttax found in a terminal of this node of backbone tree\n"
		}else{
		#print "\tNOT in list of taxnames of decendents\n";
		#print "\t\t\ttax NOT found in a terminal of this node of backbone tree\n";
		};

	# how many times has this taxon been assigned, ANYWHERE on the backbone tree
	my $number_anywhere_in_backbone = $number_of_nodes_to_which_this_taxon_has_been_assigned{$name_string};

	# is this taxon found in the fasta file, if not, ignore.
	my $how_many_in_fasta = $count_tax_in_new_data{$name_string};

	# if child taxon if found in backbone tree and new fasta members, but is not in decendents of current node, 
	# will constrained outside of monophyly
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 0)
		{$current_constraint .= "$name_string:0\n";$mpA=1;
		$count_constrained += $how_many_in_fasta;
	#	print "\t\tmpA:$mpA, constrain $child:$name_string ($how_many_in_fasta) outside MP\n"
		}

	# constrained MP
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 1)
		{$current_constraint .= "$name_string:1\n";$mpB=1;$current_constraintB .= "$name_string-";
		$count_constrained += $how_many_in_fasta;
	#	print "\t\tmpB:$mpB, constrain $child:$name_string ($how_many_in_fasta) inside MP\n"
		}

	# no information on where these should be placed (not on backbone). so they float where they want.
	if($number_anywhere_in_backbone <1)
		{
		if( $how_many_in_fasta >= 1)
			{$current_constraint .= "$name_string:-\n";
			$data_missing_reached =1;
			$count_floaters += $how_many_in_fasta;
		#	print "\t\tnot in backbone, BUT IS in fasta ($how_many_in_fasta), unknown for constraining\n";
			}else{
		#	print "\t\tnot in backbone, AND NOT in fasta, so, inconsequential.\n";
			}
		};
#	print "\t$child:$name_string how_many_in_fasta:$how_many_in_fasta\n";
#	print "\tnbr_anywhere_in_backbone:$number_anywhere_in_backbone found_in_backbone_decendent:$found_in_backbone_decendent\n";


	};#foreach my $child(@child_nodes_array)


my $constraint_length = length($current_constraint);


# mpA and mpB check that there are taxa defined within and external , i.e. a monophyly to be made.

#####################################################################

my $test1 = 0; if($mpA ==1 && $mpB == 1 && $count_floaters == 0){$test1 = 1};

my $current_constraint_string = "";
if($node_monophylys{$next} =~ /\w/ || $test1 == 1)
	{
	print "\ndefine cosntraint for BACKBONE NODE NUMBER:$next\n";	
	};

my $monophyletic_taxon_defined_from_this_node = "NA";
if($node_monophylys{$next} =~ /\w/)
	{
	$monophyletic_taxon_defined_from_this_node = $node_monophylys{$next};
	print "\tmonophyletic_taxon:$monophyletic_taxon_defined_from_this_node\n";
	$current_constraint_string = $monophyletic_taxon_defined_from_this_node;
	};

 if($mpA ==1 && $mpB == 1 && $count_floaters == 0)
	{
	$current_constraintB =~ s/-$//;
	print "\trelational constraint: $current_constraintB\n";
	$current_constraint_string = $current_constraintB;
	$ft_constraints{$current_constraint}=1;
	};

# some duplicates:
if($constraints_written{$current_constraint_string} == 1)
	{$current_constraint_string = ""};
$constraints_written{$current_constraint_string} = 1;

#####################################################################



my $child_nodes = $child_counts{$next};
my @next1 = ();

for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}


my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;

if($current_constraint_string =~ /\w/)
	{	
	$swap_string 		= "($join_the_child_nodes)$current_constraint_string";
	
	}else{
	$swap_string 		= "$join_the_child_nodes";
	
	};

$yet_another_newick_string 	=~ s/$next(\W)/$swap_string$1/;


# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	#$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	#$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


if(exists($child_counts{$next}))
	{
	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_collapse($test );#	# recurse
		#####################################
		}
	};

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_all_taxa
{


my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $all_tax_names_for_these_terminals = "";

unless($tax_list =~ /./){die "\nerror 2830, why has sub get_all_taxa been called with nothing\n"};

foreach my $tax(@tax_array)
	{
	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
	unless($test_taxonomy =~ /\w/){$test_taxonomy = $complete_lineage_for_this_species{$genus_name}};

	#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";

	while($test_taxonomy =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		unless($all_tax_names_for_these_terminals =~ /\t$current_taxname\t/)
			{$all_tax_names_for_these_terminals .= "	$current_taxname	";
			}
		}
	}
#print "\ncount tax_array:" , scalar @tax_array , ", all_members_have_this_name:$all_members_have_this_name\n";


return($all_tax_names_for_these_terminals);


}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub test_for_monophylies
{
my $tax_list = shift; # input is list of all terminals from current node
my @tax_array = split /\t/ , $tax_list;
my %all_tax_names_for_these_terminals = ();
my $number_input = scalar @tax_array;
print "\nsub test_for_monophylies, number terminals input:$number_input\n";


#count_taxa_represented_in_terminals


foreach my $tax(@tax_array)
	{
	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	unless($test_taxonomy =~ /\w/)
		{
		my $genusname = $tax;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
		$test_taxonomy = $complete_lineage_for_this_species{$genusname};
		#print "cant find taxonomy for species ($tax) trying genus ($genusname)\n";
		if($test_taxonomy =~ /\w+/)
			{
			#print "\tphew, found using genus name ($genusname), continuing ...\n";
			}else{
			die "\nerror 2889, no taxonomic informartion found for ($tax) nor ($genusname)\n";		
			}
		};
	#print "\tspecies:$tax\n";# complete taxonomy:$test_taxonomy\n";
	while($test_taxonomy =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		$all_tax_names_for_these_terminals {$current_taxname}++;
	#	print "\t\tcurrent_rank:$current_taxname:" , $all_tax_names_for_these_terminals{$current_taxname}, "\n";
		}
	}



#count_taxa_represented_in_terminals

# for a given node of the backbone tree, there may be several tax names for which the terminlas are monophyletic,
# only the lowest (tip-ward) rank is valid,
# for example just 2 terminals are input into the monophyly test sub,
# these 2 belong to Chrysidini tribe, and the only representatives of the superfamily Chrysidoidea,
# thus each of Chrysidini, Chrysidinae, Chrysididae, Chrysidoidea
# would 'monophyletic' in loosest test, but for Chrysidinae and above this would be too permissive a test,
# as for the most part there would be only one (Chrysidini in this case) of 
# presumably several fundamental representatives (other tribes),
# so, here need to get rank of each, and select lowest where several are present.

my @all_taxnames_from_terminals_of_current_node = keys %all_tax_names_for_these_terminals;
@all_taxnames_from_terminals_of_current_node = sort @all_taxnames_from_terminals_of_current_node;
my $current_no_taxnames = scalar @all_taxnames_from_terminals_of_current_node;
my $monophyly_found=0;
my $putative_monophhyly_rank = 0;my $putative_monophhyly_name = "NA";

foreach my $tax(@all_taxnames_from_terminals_of_current_node)
	{
	my $count_total_in_backbone_terminals = $count_taxa_represented_in_terminals{$tax};
	my $count_decended_from_current_node = $all_tax_names_for_these_terminals{$tax};
	my $how_many_terminals_from_current_node = scalar @tax_array;
	my $rank_number = $rank_codes{$tax};
	unless($rank_number =~ /\d/)
		{$rank_number = 1;print "wanriing, couldtn find rank number for tax:$tax\n"};

#	print "tax:$tax rank_number:$rank_number total_in_backbone_terminals:$count_total_in_backbone_terminals\n";

	if($count_total_in_backbone_terminals == $count_decended_from_current_node && 
		$count_total_in_backbone_terminals == $how_many_terminals_from_current_node)
		{
		# this was too liberal a definition:
	#	$monophylys{$tax}=1; 

		$monophyly_found = 1; print "putative monophyly:$tax\n";

		if($rank_number >= $putative_monophhyly_rank)
			{
			$putative_monophhyly_rank = $rank_number;
			$putative_monophhyly_name = $tax;
			};
		
		}
	}

 if($monophyly_found == 1 )
	{
	$monophylys{$putative_monophhyly_name} = 1;
	print "DEFINED monophyly:$putative_monophhyly_name\n";
	};


# preiovusly was this, though seems the returned varable not used
# return($all_tax_names_for_these_terminals);



# this is now required retunred:
return($putative_monophhyly_name);







};





#####################################################################################################
#
#
#
#####################################################################################################


sub place_specieslevel_data_into_taxon_constraints
{

my $newick_constraint_string 	= shift;
my $newick_length 		= length($newick_constraint_string);

print "\n
sub place_specieslevel_data_into_taxon_constraints
newick_length:$newick_length
\n";

open(VERBOSE_LOG, ">backbone_constraints_newick_verbose_LOG");


my $interal_node	= 0;

while ($newick_constraint_string =~ s/\(([^\(\)]+)\)([a-zA-Z\-]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{

	my $node = $1;my $internal_label = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";
	my $node_length = length($node);
	# contains multi-furcations.
	my @child_nodes = split /\,/ , $node;my $number_child_nodes = scalar @child_nodes;

	

	push @child_nodes , $internal_label;
	$all_internal_node_strings_assigned_to_input_constraint{$internal_label} = 1;


#	print "\nfirst adding new members of child nodes\n";
	for $i ( 0 .. $#child_nodes )
		{
		print VERBOSE_LOG "nodeID:$nodeID node_length:$node_length count_child_nodes:$#child_nodes current_child_no:$i, child:$current_child, count added $count_added_to\n";

		#print "$child_nodes[$i]\n";
		my $current_child = "-$child_nodes[$i]-";
		if($current_child =~ /[\,]/){die "\nerror 3084:$current_child\n"};
		# $new_member_assigned_to_constraint_taxon { $new } = $current_member_to_be_assigned;
		# my @testkeys = keys %new_member_assigned_to_constraint_taxon;
		my @add = ();

		foreach my $new ( @testkeys )
			{
			my $is_assigned_to = $new_member_assigned_to_constraint_taxon { $new };
			print VERBOSE_LOG "\tnew:$new is_assigned_to:$is_assigned_to\n";
			if($current_child =~ /-$is_assigned_to-/)
				{
				unless($this_added_already{$new} == 1) 
					{
					print VERBOSE_LOG "\t\tADDING\n";
					 push @add , $new 
					};
				$this_added_already{$new} = 1;
				};
			};
		my $count_added_to = scalar @add;
		$current_child =~ s/^-(.+)-$/$1/;
		if($count_added_to >= 1){$add_these_to_this{$current_child} = join ',' , @add};

		};

#	print "node:$interal_node, number_child_nodes:$number_child_nodes int_label:$internal_label\n\tchild nodes:@child_nodes\n";

	# node:4, int_label:Thysanoptera
	# 	child nodes:Thripinae
	# node:5, int_label:Psocoptera-Thysanoptera-Phthiraptera
	# 	child nodes:Psocoptera-Phthiraptera INTERNAL_NODE_4





	$interal_node++;
	}


print "
remaining \$newick_constraint_string (should be nothing):$newick_constraint_string
";




# finally, string replace newick contraints with new data

my @keys = keys %add_these_to_this;

# open(LOG2 , ">$outfile_prefix.bcn.LOG") || die "\nerror 3441\n";

foreach my $key(@keys)
	{
	my $replace_with = $add_these_to_this{$key};
	my $count_members = scalar split /\,/ , $replace_with;

	my $screen_print = $replace_with;$screen_print =~ s/^(.{120}).+/$1........../;
	#if(length($replace_with)<= 1000){
		print "replcace assign tax $key with $count_members members:$screen_print\n";
	#	print LOG2 "replcace assign tax $key with $count_members members:$screen_print\n";
print VERBOSE_LOG "replcace assign tax $key with $count_members members:$screen_print\n";

	#};

		# put comma after bracket otherwise will end up as node label:
	if($yet_another_newick_string =~ s/([\)])$key([\(\)\,])/$1,$replace_with$2/) # bugfix 20160928
		{}elsif($yet_another_newick_string =~ s/([\(\)\,])$key([\(\)\,])/$1$replace_with$2/)
		{
		my $test1 = $1;my $test2 = $2;#print "\tpreceednig char:$test1 following char:$test2\n";
		}else{
		print "\nerror 3441???\n";
		};

	};

# close LOG2;
close VERBOSE_LOG;

# correct clade with onyl one tip. a few node lables will be removed by this, but most remain.

# (Challia)Orthopteroidea-Neoptera_incertae_sedis,
$yet_another_newick_string =~ s/\(([A-Z][a-z]+)\)[a-zA-Z\-\_]+/$1/g;
# ERROR: Expecting ',' in tree; found: character ')'
# chnura)),((((Challia),(Mesocapnia,Aptero
$yet_another_newick_string =~ s/\(([a-z]+)\)\,/$1,/ig;






} # sub

#####################################################################################################
#
#
#
#####################################################################################################



sub backbone_constraints
{


	# THIS SUB NOT RUN



# this didnt work, if backbone trees were perfect then maybe it would

print "
sub backbone_constraints
";

my @all_backbone_terminals = keys %terminals;@all_backbone_terminals = sort @all_backbone_terminals;


	# THIS SUB NOT RUN


foreach my $termnial(@all_backbone_terminals)
	{
	my $parent = $nodes{$termnial}{parent};
#	print "termnial:$termnial parent:$parent\n";
	my $sister_node;
	for my $i(0,1)
		{
		my $child_node_of_parent = $nodes{$parent}{$i};
		if($child_node_of_parent eq $termnial)
			{}else{
			$sister_node = $child_node_of_parent;
			};
	#	print "\ti:$i child_node_of_parent:$child_node_of_parent\n";
		};
	unless($termnial =~ /\w/ && $parent =~ /\w/ && $sister_node=~ /\w/){ die "\nerror 3552\n"};
	print "\tsister_node:$sister_node\n";

	my $most_inclusive_taxa_of_sister_node = "";
	if($sister_node =~ /INTERNAL_NODE_/)
		{
		$all_terminals_from_this_node = "";
		get_terminals_from_this_node2($sister_node, 0);#
		$most_inclusive_taxa_of_sister_node = get_shared_taxa($all_terminals_from_this_node);#
		unless($most_inclusive_taxa_of_sister_node =~ s/.+\t(\w+)\t$/$1/)
			{
			print "\nwarning, expecting lineage here, algorthim won work for this terminal\n";
			$most_inclusive_taxa_of_sister_node = "NA";
			};
		}else{
		$most_inclusive_taxa_of_sister_node = $sister_node;
		};
	print "\tmost_inclusive_taxa_of_sister_node:$most_inclusive_taxa_of_sister_node\n";

	# THIS SUB NOT RUN


	################################################################################

	unless($most_inclusive_taxa_of_sister_node eq "NA")
		{
		# for genus of current terminal, walk up the taxonomic ranks, 
		# until taxa is shared by sister of backbone tree.
		
		my $test_species = $termnial;
		my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		unless($test_taxonomy =~ /\w+/)
			{
			my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
			$test_taxonomy = $complete_lineage_for_this_species{$genusname};
			if($test_taxonomy =~ /\w+/)
				{}else{
				if($test_species =~ /optera_[A-Z]/)
					{print "\nlooks like user format was not processed to an established tax name\n"};
				die "\nerror 3781, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
				};
			};


	# THIS SUB NOT RUN


		print "test_taxonomy:($test_taxonomy)\n";
		my $upper_limit=0;my $upper_taxon = "NA";
		while($test_taxonomy =~ s/\s+(\w+):(\w+)\s+$/ /)
			{
			my $current_rank = $1;my $current_taxname = $2;

			unless($upper_limit == 1)
			{
			print "\t\tNEW rank:$current_rank taxname:$current_taxname\n";
			my $ncbi_taxnumber = $ncbi_taxnumber_for_taxname{$current_taxname};
			unless($ncbi_taxnumber =~ /[\d\w]/){die "\ndidnt find tax number for name\n"};

			# traversing heirachy from this node, 
			# checking no offspring are sister taxa
			$all_taxa_from_this_node="";
			traverse_heirachy($ncbi_taxnumber);
		#	print "all_taxa_from_this_node:$all_taxa_from_this_node\n";
			if($all_taxa_from_this_node =~ /\t$most_inclusive_taxa_of_sister_node\t/)
				{
				$upper_limit = 1;
				};
			if($upper_limit == 0){$upper_taxon = $current_taxname};

			};
			};
			print "highest taxon not shared by backbone sister:$upper_taxon\n";
		$backbone_terminal_rename{ $termnial } = $upper_taxon;
		};
};



	# THIS SUB NOT RUN


foreach my $termnial(@all_backbone_terminals)
	{
	my $replace_backbone_genus_with_string = "";

	if($backbone_terminal_rename{ $termnial } =~ /\w/)
		{
		$replace_backbone_genus_with_string = $backbone_terminal_rename{ $termnial };

		if($process_backbone_tree_terminal_IDs == 0)
			{
			$second_backbone_copy_newick =~ s/$termnial([^a-zA-Z])/$replace_backbone_genus_with_string-$termnial$1/;
			}else{
			$second_backbone_copy_newick =~ s/$termnial([\_])/$replace_backbone_genus_with_string-$termnial$1/;
			};		

		};


	};

	# THIS SUB NOT RUN

close TEST_FILE;


open (OUT6999 , ">$outfile_prefix.less_basic_constraint_tree") || die "\n\n";
print OUT6999 "$second_backbone_copy_newick\n";
close OUT6999;



	# THIS SUB NOT RUN




die "";




};



#####################################################################################################
#
#
#
#####################################################################################################






sub traverse_heirachy
{
my $current_node = $_[0]; #my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated
my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);

if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;
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
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	my $name_assignment_to_taxnumber = "";
	$all_taxa_from_this_node .= "\t$originalname\t";

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

		###########################################
		traverse_heirachy($child );#
		###########################################
	}}

	
}#sub traverse_heirachy





#####################################################################################################
#
#
#
#####################################################################################################




#####################################################################################################
#
#
#
#####################################################################################################



sub backbone_constraints2
{


# attempt 2 ...............







print "
sub backbone_constraints2
";

my @all_backbone_terminals = keys %terminals;@all_backbone_terminals = sort @all_backbone_terminals;

foreach my $termnial(@all_backbone_terminals)
	{
	$terminal_count++;
	$terminal_print2++;if($terminal_print2 <= 6){print "$terminal_count of $#all_backbone_terminals\n"}
	if($terminal_count =~ /00$/){$terminal_print++;if($terminal_print >= 7 && $terminal_print <= 12){print "$terminal_count of $#all_backbone_terminals\n"}};

		my $test_species = $termnial;my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		unless($test_taxonomy =~ /\w+/)
			{
			my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;$test_taxonomy = $complete_lineage_for_this_species{$genusname};
			if($test_taxonomy =~ /\w+/)
				{}else{
				if($test_species =~ /optera_[A-Z]/){print "\nlooks like user format was not processed to an established tax name\n"};

				print "\nerror 4000, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
				if($genusname =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/){}else{
					#die ""
					};

				};
			}; # print "\n\nTERMINAL:$termnial test_taxonomy:($test_taxonomy)\n";
		my $upper_limit=0;my $upper_taxon = "NA";my $break_loop = 0;
		while($test_taxonomy =~ s/\s+(\w+):(\w+)\s+$/ /)
			{
			my $current_rank = $1;my $current_taxname = $2;
		#	print "\trank:$current_rank tax:$current_taxname\n";

			unless($break_loop == 1)
			{
			foreach my $termnial2(@all_backbone_terminals)
				{
				unless($termnial2 eq $termnial)
					{
					my $test_species2 = $termnial2;
					my $test_taxonomy2 = $complete_lineage_for_this_species{$test_species2};
					unless($test_taxonomy2 =~ /\w+/)
						{
						my $genusname2 = $test_species2;$genusname2 =~ s/^([A-Z][a-z]+)_.+/$1/;
						$test_taxonomy2 = $complete_lineage_for_this_species{$genusname2};
						if($test_taxonomy2 =~ /\w+/)
							{}else{
							if($test_species2 =~ /optera_[A-Z]/){print "\nlooks like user format was not processed to an established tax name\n"};
							print "\nerror 923, no taxonomy found for species ($test_species2) nor genus ($genusname2). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
							if($genusname2 =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/)
								{
								}else{
								# die ""
								};
							};
						};#print "\ttermnial2:$termnial2 test_taxonomy2:($test_taxonomy2)\n";
					#######################
					while($test_taxonomy2 =~ s/\s+(\w+):(\w+)\s+$/ /)
						{
						my $current_rank2 = $1;my $current_taxname2 = $2;
						if($current_taxname eq $current_taxname2)
							{
							$break_loop = 1;
							};

						};

					#####################
					};
				};# foreach my $termnial2(@all_backbone_terminals)

			unless($break_loop == 1){$upper_taxon = $current_taxname};

			};# unless($break_loop == 1) # this stops at upper taxon
			}; # while($test_taxonomy =~


	

		#	print "terminal:$termnial new_taxon:$upper_taxon\n";
		$backbone_terminal_rename{ $termnial } = $upper_taxon;


			
	}; # foreach my $termnial(@all_backbone_terminals)


print "done ... looking through terminals again, replacing strings ... \n";

open(LOG_OUT , ">Backbone_Constraints_Newick.list_of_constraints") || die "";

$terminal_count = 0;$terminal_print2 = 0;$terminal_print = 0;
foreach my $termnial(@all_backbone_terminals)
	{
	$terminal_count++;
	$terminal_print2++;
	if($terminal_count =~ /00$/){$terminal_print++;if($terminal_print >= 7 && $terminal_print <= 12){print "$terminal_count of $#all_backbone_terminals\n"}};

	my $replace_backbone_genus_with_string = "";

	if($backbone_terminal_rename{ $termnial } =~ /\w/)
		{

		$replace_backbone_genus_with_string = $backbone_terminal_rename{ $termnial };
		$backbone_substrings_original_length += $termnial;
		$backbone_substrings_replace_length += $replace_backbone_genus_with_string;

		if($process_backbone_tree_terminal_IDs == 0)
			{
			unless($second_backbone_copy_newick =~ s/([\,\(])$termnial([^a-zA-Z])/$1$replace_backbone_genus_with_string-$termnial$2/)
				{die "\nwarning, could not replace string in newick:$termnial\n"};
			
			}else{
			unless($second_backbone_copy_newick =~ s/([\,\(])$termnial([\_])/$1$replace_backbone_genus_with_string-$termnial$2/)
				{die "\nwarning, could not replace string in newick:$termnial\n"};
			};


		unless($backbone_terminal_rename{ $termnial } eq "NA")
			{
		#	$store_all_constraints{$replace_backbone_genus_with_string} = "$replace_backbone_genus_with_string-$termnial"
		# try and make simple as possible , replace with new tax lists, on original newick string, not taxon assigned string:
			$store_all_constraints{$termnial} = $replace_backbone_genus_with_string;
			};

		if($terminal_print2 <= 6)
			{
			print "$terminal_count of $#all_backbone_terminals:$termnial\n";
			print "  replace terminal ($termnial) with ($replace_backbone_genus_with_string)\n";
			};



		}; # 	if($backbone_terminal_rename{ $termnial } =~ /\w/)

	 print LOG_OUT "$termnial\t$replace_backbone_genus_with_string\n";

	};

close TEST_FILE;
close LOG_OUT;

print "
backbone_substrings_original_length:$backbone_substrings_original_length
backbone_substrings_replace_length:$backbone_substrings_replace_length
";



open (OUT69999 , ">$outfile_prefix.less_basic_constraint_tree.backbone_taxa_assigned") || die "\n\n";
print OUT69999 "$second_backbone_copy_newick\n";
close OUT69999;

open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";


my @list_constraints = keys %store_all_constraints; @list_constraints = sort @list_constraints;

my $cosntr_number = 0;
foreach my $termnial(@list_constraints)# this list of terminals is not neccessary complete,
	{					# will omit those in which upper_taxon inference did not work

	my $constreaint = $store_all_constraints{$termnial};$cosntr_number++;
#	print "termnial:$termnial constreaint:$constreaint\n";

	if($fasta_ids_for_each_taxon{$constreaint} =~ /\S/)	
		{
		my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
		$new_members_included_in_constraint =~ s/^\t+//;$new_members_included_in_constraint =~ s/\t+$//;
		$new_members_included_in_constraint =~ s/\t\t+/	/g;
		my @split_memers = split /\t+/ , $new_members_included_in_constraint;
		my $count_new_members = scalar @split_memers;
		print BCN_LOG "CONSTR number:$cosntr_number, termnial:$termnial, taxon_assigned:$constreaint, count_new:$count_new_members new_members_incl:$new_members_included_in_constraint\n";

#		my $find_string = $store_all_constraints{$constreaint};
		my $find_string = $termnial;
		my $replace_with = join ',' , @split_memers;

		# $third_backbone_copy_newick

		if($count_new_members >= 2)
			{

			if($newick_reformat == 1)# 1 = backbone ID are binomial, convert to genus name only
				{

				if($third_backbone_copy_newick =~ /(.{10}$find_string.{20})/)
					{
					my $context = $1;$context_print++;
					if($context_print <= 6){print BCN_LOG "\tfound location in newick to replace:$context\n"};
					};

				# if($count_new_members <= 9){$constreaint = ""}; # !!!!

				$forth_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . "(" . $replace_with . ")" . $constreaint /e;

					# this is the important string replace:
				if($process_backbone_tree_terminal_IDs == 0)
					{
					if($third_backbone_copy_newick =~ s/([\(\,])$find_string(\:)/ $1 . "(" . $replace_with . ")" . $2/e)
						{
						$sucessful_data_insertions++;
						$list_all_inserted_things .= ",$replace_with,";
						}else{$unsucessful_data_insertions++};			
					
					}else{
					if($third_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . "(" . $replace_with . ")" /e)
						{$sucessful_data_insertions++}else{$unsucessful_data_insertions++};			
					};
				}else{
				die "\nnot implemented\n";
				};			


			}else{
			if($find_string eq $replace_with)
				{
				$list_all_inserted_things .= ",$replace_with,";

				}else{

				if($newick_reformat == 1)# 1 = backbone ID are binomial, convert to genus name only
					{
					# remember, higher taxa might have been assigned.
					# and if only one genus is in the new data,
					# it may not neccessary be the same genus as the one in the backbone.
					# hence these might be different here, still need to replcae
					if($process_backbone_tree_terminal_IDs == 0)
						{
						$third_backbone_copy_newick =~ s/([\(\,])$find_string(\:)/ $1 . $replace_with . $2/e;			
						$list_all_inserted_things .= ",$replace_with,";

						}else{
						$third_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . $replace_with /e;			
						};
					}else{
					die "not impleented\n";
					};
				};

			};

		foreach my $mem(@split_memers)
			{
			if($new_members_assigned{$mem} =~ /./){die "\nerror, $mem assigned twice\n"};
			$new_members_assigned{$mem} = $constreaint;
		#	print "$mem contrained to $constreaint\n";
			};


		}else{
		push @needs_pruning , $constreaint;
	#	print "CONSTR:$constreaint, NO new membres\n";
		};
	};

my @new_members_assigned_to_backbine_constraint = keys %new_members_assigned;
print "

sucessful_data_insertions:$sucessful_data_insertions++
unsucessful_data_insertions:$unsucessful_data_insertions++

$#needs_pruning need pruning:@needs_pruning
new_members_assigned_to_backbine_constraint:$#new_members_assigned_to_backbine_constraint

";
close BCN_LOG;



 open(TEST_FILE , ">$outfile_prefix.needs_pruning") || die "\nerror 759\n";
foreach my $taxas(@needs_pruning)
	{
	print TEST_FILE ">$taxas\nACTG\n"
	}
close TEST_FILE;



open (OUT69999999 , ">$outfile_prefix.less_basic_constraint_tree") || die "\n\n";
print OUT69999999 "$third_backbone_copy_newick\n";
close OUT69999999;

open (OUT699999998 , ">$outfile_prefix.less_basic_constraint_tree.internal_tax_labeled") || die "\n\n";
print OUT699999998 "$forth_backbone_copy_newick\n";
close OUT699999998;


$list_all_inserted_things =~ s/^\,+//;$list_all_inserted_things =~ s/\,+$//;
my @splitit = split /\,+/ , $list_all_inserted_things;
open(FLE, ">$outfile_prefix.list_constrained_members") || die "";
print FLE "@new_members_assigned_to_backbine_constraint\n";
foreach my $thing(@new_members_assigned_to_backbine_constraint) # @splitit)
	{
	print FLE "$thing\t$new_members_assigned{$thing}\n";	
	};
close FLE;


print "
although script has lots of code from various method implementations, currently only doing one thing.
	check taxa assigned to backbone tree in file:
	$outfile_prefix.less_basic_constraint_tree.backbone_taxa_assigned

	and constraint tree is
	$outfile_prefix.less_basic_constraint_tree
	which will probably require some pruning before doing anything else.
	
	backbone_constraints_newick_LOG list each constraint, and taxa assigned to it


";


};



#####################################################################################################
#
#
#
#####################################################################################################



sub print_fasttree_format_relational_constraints
{
my $next = $_[0];my $sum_branchlength = $_[1];
$ft_relational_constraints_subcalls++;
if($ft_relational_constraints_subcalls <= 6){print "ft_relational_constraints_subcalls:$ft_relational_constraints_subcalls\n"};
if($ft_relational_constraints_subcalls =~ /000$/){print "ft_relational_constraints_subcalls:$ft_relational_constraints_subcalls\n"};

# print "\nNEW NODE ($next)\n";

my @non_terminal_child_nodes = ();
my @child_nodes_test = ();
my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();

my $lowest_taxname_for_node = "";
my $lineage_for_node = "";
$taxa_for_terminals_of_node = "";


if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{



	}else{#if($next =~ /INTERNAL_NODE_\d/)

	# for tips, you still need to know the lineage. but it wont be a shared lineage.
	# in this case you dont have a list of terminals, theres just one, the curent node

#	$taxa_for_terminals_of_node = get_all_taxa($next);# bug fix 22 oct 2015


	};

if($node_label{$parent_node} =~ /\w/)
	{
	$parent_taxon = $node_label{$parent_node};
	#print "node($next), tax($lowest_taxname_for_node) parent($parent_node) Tax($parent_taxon)\n";
	}else{#if($node_label{$parent_node} =~ /\w/)
	
	}



my $child_nodes = $child_counts{$next};
my @next1 = ();
my @current_node_relational_constraint = ();

# record the taxonomic name(s) given to this node
$nodes{$next}{node_taxonomic_label}=$label_string_to_assign;


my $terminal_members_found = 0;

for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	$all_terminals_from_this_node ="";

	if($push_node =~ /INTERNAL_NODE_\d/) 
		{
		##########################################
		get_terminals_from_this_node6($push_node, 0);#
		##########################################
		}else{

		my $constreaint = $store_all_constraints{$push_node}; # print "the termnial:$push_node the constreaint:$constreaint\n";
		if($fasta_ids_for_each_taxon{$constreaint} =~ /\S/)	
			{
			my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
			$all_terminals_from_this_node .= "	$new_members_included_in_constraint	";
			}
		};

	if($all_terminals_from_this_node =~ /[\w\d]/)
		{
		$terminal_members_found++;
		push @current_node_relational_constraint , $all_terminals_from_this_node
		};

	};

if($terminal_members_found >= 2)
	{
#	print "\nFT relational contraint:\n";
	for my $j(0 .. $#current_node_relational_constraint)
		{
		my $contr98 = $current_node_relational_constraint[$j]; # print "$j $contr98\n";
		$relational_contraints{$next} .= "\t$contr98\tRelationalConstraint$j\t";
		};
	};





if(exists($child_counts{$next}))
	{

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		print_fasttree_format_relational_constraints($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################

		}

	};

return();

}


#####################################################################################################
#
#
#
#####################################################################################################



sub get_terminals_from_this_node6
{
my $next = shift;

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	
	#print "line1423node:$next i:$i child:$nodes{$next}{$i}\n";
	}



for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test 1430:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node6($test );#	# recurse
		#####################################

		}else{
		$count_terminal_nodes++;

		if($test eq "")
			{
			print "\nnext:($next)\n";
			die "\nerror 4435. looks like sub get_terminals_from_this_node6 was called with a terminal node.\n"
			}



		my $constreaint = $store_all_constraints{$test};$cosntr_number++;
	#	print "termnial:$termnial constreaint:$constreaint\n";

		if($fasta_ids_for_each_taxon{$constreaint} =~ /\S/)	
			{
			my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
		#	$new_members_included_in_constraint =~ s/^\t+//;
		#	$new_members_included_in_constraint =~ s/\t+$//;
		#	$new_members_included_in_constraint =~ s/\t\t+/	/g;
			$all_terminals_from_this_node .= "	$new_members_included_in_constraint	";

		#	my @split_memers = split /\t+/ , $new_members_included_in_constraint;
		#	my $count_new_members = scalar @split_memers;
			}

			}
	}

return();

} # sub get_terminals_from_this_node6



#####################################################################################################
#
#
#
#####################################################################################################



