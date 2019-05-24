#!/usr/bin/perl

###################################################################################################################################
#
#
#		  
#    	Copyright (C) 2015-2018 Douglas Chesters
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
#	
#	
##########################################################################################################################################
#
#	
#	
#	
#	
# 	script reads the IDs in the fasta file. 
#	any ID in the tree not in the fasta file, is pruned.
#	
#	
#	
#	
#	prune tree arguments:
#		-treefile XXX, -seqfile XXX, -output XXX -start_prune_list tax1 tax2 -end_prune_list
#	
#
#	
#	
#	change log
#	08 Apr 2015: uses branch lengths.
#	29 Jan 2018: some minor details printed to screen.
#	22 Mar 2018: if pruning makes extra pair of parentheses near root, these are removed
#	28 Apr 2018: taxa to be pruned can be provided as list in the command:
#		-start_prune_list TR_Tetragonula_carbonaria -end_prune_list
#	20 Nov 2018: still prints output even if nothing found to prune, in case used in pipeline (as i do use this)
#	08 Dec 2018: bufix in processing start and end of output newick string
#	04 May 2019: workaround for reading newick strings with common error usually introduced during rerooting: 
#		missing comma before outgroup.
#	07 May 2019: prints list of terminal in pruned tree, which may or may not be known when run, 
#		depends on type of list input.
#	09 May 2019: option -process_tree_labels, when == 1, trim off species names for binomials in input tree
#		also, prints usefull message if unrooted tree input by mistake.
#	
#	
#	
#	
#	
#
##########################################################################################################################################



my $arg_string  = join ' ', @ARGV;



# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose

# $tree_ID_format = 2;# i think 1=Order_Genus, 2= Order_Genus_species


# $process_backbone_tree_terminal_IDs = 1;
# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis
# if == 2

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


open(LOG, ">prune_tree_LOG") || die "";


###############################################################################################################################################



my $root_identity = "";
my $date = localtime time;




# user might give terminal to prune according to overlap with fasta file

if($fas_file =~ /\w/)
	{

	%fasta_IDs;	# ids from fasta file are stored in this hash
	############
	read_fas();#
	############

	};



#######################################################################

print LOG "
reading newick....
";


%terminals = ();
%nodes;		# user phylogeny stored in this hash
$root_node;	# this global will keep the value of the last in the loop, which will be the root
%child_counts;


# read backbone tree into hash. node identities for hash key, then parent/child identities as hash entry

#########################
record_tree2($treefile);#
#########################

my @backbone_tree_terminals = keys %terminals;@backbone_tree_terminals = sort @backbone_tree_terminals;
my $count_to_prune=0;my $count_not_to_prune=0;my $count_all=0;


foreach my $tax(@backbone_tree_terminals)
	{
	$count_all++;#print "backbone tax:$tax\n";

	if($fas_file =~ /\w/)
		{
		if(exists($fasta_IDs{$tax}))
			{
			$count_not_to_prune++;#print "\texists\n";
			}else{
			$prune_these{$tax}=1;$count_to_prune++;#print "\tPRUNE $tax\n";
			};
		};
	};


# or prune terminals given in command

if($prune_string =~ /\w/)
	{
	my @split_list = split /\s+/ , $prune_string;
	foreach my $tax5(@split_list)
		{
		$item_number++; $prune_these{$tax5}=1;$count_to_prune++; 
		print "item $item_number from prune list:$tax5\n";
		unless(exists($terminals{$tax5}))
			{die "\nerror, cant find prune taxon $tax5 in tree terminals.\n"};
		};
	};



if($count_all == $count_to_prune){die "\nsome error, pruning everything!\n"};


if($count_to_prune >=1)
	{
	print "will prune $count_to_prune from backbone tree, will retain:$count_not_to_prune\n";

	if($count_to_prune <=10)
		{
		my @prune_array = keys %prune_these;@prune_array = sort @prune_array;
		print "prune list:@prune_array\n";
		};

	###############
	prune_nodes();#
	###############


	}else{
	print "\n\nNOTHING to prune. HOWEVER, in case this is being run within a pipeline, will print unpruned tree to output.\n";	

	};






$new_newick_string3 = "($root_node)";

print "writing new newick string ... \n";
print LOG "\n\n  .... writing new newick string .... \n";

###############################################
get_terminals_from_this_node3($root_node);#
###############################################


open(OUTLIST, ">list_terminals_pruned_tree");
foreach my $terinal(@store_terminals_of_pruned_tree)
	{
	print OUTLIST "$terinal\n";
	};
close OUTLIST;




if($count_not_to_prune >= 500 && $newick_substring_replacements <= 4)
	{
	print "\n\nWARNING. large tree should have been printed but seems not done so. implies fail.\n\n";
	}else{
	print "\tnewick built, substring replacements:$newick_substring_replacements\n";
	
	};





 if($new_newick_string3 =~ s/^\((\(.+\))\)$/$1/)
	{print "extra parentheses removed\n"};


open(NEWICK5 , ">$output_filename") || die "\nerror 193.\n";

 print NEWICK5 "$new_newick_string3;\n";
# print NEWICK5 "($new_newick_string3:1.0);\n";


#print "$newick5\n";
close NEWICK5;


print "pruned tree writen to $output_filename

depending on whats pruned, may have an extra pair of parentheses at extremes of string, 
	you can manually remove these
	unforunately, the opposite can also occur ! you may need to add a pair to ends.
	thus, please check file 
";

close LOG;

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


$tree_parse = 1; 	
if($tree_parse ==1)
	{

	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	$tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 
	# remove regular branchlengths: 0.02048
	$tree1 =~ s/\:\-*\d+\.\d+//g; 
	# remove 0 length branchlengths
	$tree1 =~ s/\:\d+//g;

	# sumtrees node supports
	$tree1 =~ s/\)[01]\.[0-9]+/\)/g; 

	# whole number node support
	# estris)100)100,UCE_v2_Bomb    &   ris))100,UC
	while($tree1 =~ s/(\))100([\)\,])/$1$2/){};
	# 1 or 2 digits
	while($tree1 =~ s/(\))[0-9]{1,2}([\)\,])/$1$2/){};

	};


# print "tree1:$tree1\n";die;

my $newick_string 	= $tree1;
my $interal_node	= 0;



# new newick parser ....

#while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
#	{
	#my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";

#while ($newick_string =~ s/\(([^\(\)]+)\)([^\(\)\,]+)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
#	{
#	my $node = $1;my $nodelab = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; 

# try leaving node lable to follwoing iterations:
while ($newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
	{
	my $node = $1;my $nodeID = "INTERNAL_NODE_$interal_node";
#	print "\nNEW NODE($node)\n";

	my @child_nodes = split /\,/ , $node;
	$child_counts{$nodeID} = $#child_nodes;

	for $i(0 .. $#child_nodes)
		{
		my $bl = "1.0"; 	
		if($child_nodes[$i] =~ s/\:(.+)//)	# branchlength found
			{$bl = $1};

		#print "parent($nodeID) child $i of $#child_nodes:($child_nodes[$i]), bl:$bl\n";

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($process_tree_labels == 1){$child_nodes[$i] =~ s/^([A-Z][a-z]+)_[a-z]+$/$1/};
			};

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;

		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++;
			$tips_printed_to_screen++;
			if($tips_printed_to_screen < 6){print "tree tip:$child_nodes[$i]\n"};
			if($tips_printed_to_screen == 6){print " .....\n"};

			}

		}
#	print "node:$interal_node\n\tchild nodes:@child_nodes\n";
#	print "node:$interal_node\tchild nodes:$#child_nodes\n";

	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;

	if(length($newick_string) <= 1000)
		{
		print LOG "newick_string:$newick_string\n";

		if($newick_string =~ s/^(\(\(INTERNAL_NODE_\d+\))([A-Z][a-z_]+\)\;)$/$1,$2/)
			{
			print "\n\nWARNING: 
	problem with your newick string, missing comma before outgroup, probably done during rerooting. 
	hacked correction which may or may not work\n\n"
			};


		};
	}



print "
your backbone tree has been read, TIPS:$count_the_terminal_nodes INTERNAL NODES:$interal_node.
\n";


if(length($newick_string)>=40)
	{
	print "unexpected length of newick string remaining after parse:$newick_string\n";
	die "\nseems tree not read correctly\n"
	
	}else{
	print "seems tree is read correctly\n";
	};


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




sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";

my $count_diet_IDs = 0;
my $items_printed_to_screen = 0;
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1;#print "id:$id\n";
		$items_printed_to_screen++;if($items_printed_to_screen < 6){print "fasta ID:$id\n"};
		if($items_printed_to_screen == 6){print "....\n"};

		if($id =~ /^[a-z]+_[a-z]+$/i)
			{
			$count_binomials++;
			}elsif($id =~ /^[a-z]+_[a-z]+_[a-z]+/i)
				{
				$extended_ids++;
				};
	#	$index_of_fastaID{$id} = $count_diet_IDs;
		$count_diet_IDs++;
		$fasta_IDs{$id}= 1;
		}else{
		if($line =~ /^>/){die "\nerror 646.\n"}
		}
	}
close(IN);

print "
count_binomials:$count_binomials 
extended_ids:$extended_ids
";
$autodetect_fasta_IDs = 0;
if($count_diet_IDs == $count_binomials)
	{$autodetect_fasta_IDs = "binomial";print "\nautodetected strict binomial used in your fasta file. good!\n"};
if($count_diet_IDs == $extended_ids)
	{$autodetect_fasta_IDs = "extended";print "\nwarning. you are not using strict binomials in your fasta file, not guarunteed to work\n"};

my @check_species_filtered = keys %fasta_IDs;

if($count_diet_IDs == scalar @check_species_filtered)
	{
	print "\npass check, all strings in fasta file are unique (species filtered)\n"
	}else{
	print "\nwarning. duplicate IDs found in your fasta file.\n"
	}

print "\n*** $count_diet_IDs fasta IDs *** in your fasta file ($fas_file). will be looking for these in the tree ...\n\n";

unless($count_diet_IDs >= 1)
	{die "\nerror ... no fasta ID's found in the file ($fas_file). quitting.\n"}




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




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
\n\n
     ***** prune_tree.pl *****
		  
\n";


# perl parse_clades.pl -treefile RAxML_bipartitions.diet_group.0.reformatted -seqfile diet_group.0.fas -keyfile key_Mar2013_Magnoliophyta
#$treefile;
#$fas_file;
#$output_filename 	= "$treefile.diet_clades";



if($arguments =~ /-treefile\s+(\S+)/)
	{
	$treefile = $1;
	}else{
	print "\nerror reading command arguments  (-treefile)\n";die""
	}


if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	unless($arguments =~ /\-start_prune_list\s+(.+)/)
		{die "\ninsufficient input\n"};
	}

if($arguments =~ /\-start_prune_list\s+(.+)\s+\-end_prune_list/)
	{
	$prune_string = $1;
	}else{
	}


if($arguments =~ /\-output\s+(\S+)/)
	{
	$output_filename = $1;
	}else{
	$output_filename = "prune_tree_output.nwk";
	};

if($arguments =~ /\-process_tree_labels\s+(\d)/)
	{
	$process_tree_labels = $1;
	}else{
	};


#$output_filename	= "$treefile.query_clades";


print "\n
user options have been read.....
treefile:			$treefile
query fasta file:		$fas_file
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








sub prune_nodes
{
print "sub prune_nodes\n";

my @prune_tax = keys %prune_these;@prune_tax = sort @prune_tax;

foreach my $tax(@prune_tax)
	{
#	print "will prune tax:$tax\n";
	my $parent  = $nodes{$tax}{parent};
	unless($parent =~ /[\d\w]/){die "\nerror 485. no parent node for $tax\n"};
	my $grandparent  = $nodes{$parent}{parent};
	unless($grandparent =~ /[\d\w]/)
		{
		print "\nerror 487. no grandparent node for $tax (parent:$parent) .... tip, " , 
			"this script probably wont work on unrooted trees, perhaps your tree is rooted.\n";
		die "";
		};

	my $count_decendents_of_parent = $child_counts{$parent};# will be -1, $child_counts{$nodeID} = $#child_nodes;
	unless($count_decendents_of_parent =~ /\d/){die "\nerror 2235\n"}
	#print "parent_node:$parent has decendents:$count_decendents_of_parent and grandparent:$grandparent\n";

#	unless($count_decendents_of_parent == 1){die "\nnot written for non-totally bifurcating. (\$count_decendents_of_parent == $count_decendents_of_parent)\n"}

	my $children_of_grandparent = $child_counts{$grandparent};
	unless($children_of_grandparent =~ /\d/){die "\nerror 2242\n"}

	for my $i(0 .. $children_of_grandparent)
		{
		# go through all children of grantparent, need to find which $i it is
		my $child  = $nodes{$grandparent}{$i};
		#print "\tchild of grandparent:$grandparent i:$i is:$child\n";
		if($child eq $parent)	
			{
			# correct lineage found
			# parent node is no longer required, so link sister node to grandparent node
			#delete($nodes{$parent_id});
			for my $j(0 .. $count_decendents_of_parent)
				{
				# go through all children of parent, need to find all sisters of node to be deleted.
				my $child2  = $nodes{$parent}{$j};
				if($child2 eq $tax)
					{
					# arrived at node to be deleted. so ignore.
					delete($nodes{$tax});$nodes_deleted++;
					}else{
					# arrived at sister of node to be deleted.
					# this is now the child of the grandparent
					$nodes{$grandparent}{$i} = $child2;
					
					# the branchlength of the child of the grandparent is
					# now the length of the previous child branch plus the length of the grandchild 
					$nodes{$grandparent}{branchlength} += $nodes{$parent}{branchlength};

					# grandparent is now parent of sister node:
					$nodes{$child2}{parent} = $grandparent;
					
					};
				};
			};
			
			
		}


	
	};#foreach my $tax(@prune_tax)


print "\tnodes_deleted:$nodes_deleted\n"

};#sub prune_nodes


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub get_terminals_from_this_node3
{
my $next = shift;
my $child_nodes = $child_counts{$next};
my @next1 = ();my @new_node = ();

for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	my $push_node = $nodes{$next}{$i};
	my $bl = $nodes{$push_node}{branchlength};
	push @next1 , $push_node;
	push @new_node , "$push_node:$bl";	
	}

# default:
my $join_child_nodes = join ',', @new_node;
$swap_string = "($join_child_nodes)";
if($new_newick_string3 =~ s/$next(\W)/$swap_string$1/){$newick_substring_replacements++};
print LOG "current node:$next child_count:$child_nodes child nodes:(@next1)\n";

for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node3($test);#	# recurse
		#####################################
		}else{
		push @store_terminals_of_pruned_tree , $test;
		};
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







