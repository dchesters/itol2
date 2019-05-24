

#
# 	sample_db.pl
#
#    	Copyright (C) 2014  Douglas Chesters
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
#	############################################################################################## 
#
#
#
#
# 	to take a random sample of 60 from a fasta file: 
# 	perl ~/usr_scripts/sample_db.pl in.fas 60 2 out.fas NULL
#
#
#
# 	change log
#
# 	02OCT2013: taxonomic sampling using keyfile not tobycodes.
# 	26AUG2014: cleaned up a little
# 
# 
# 
# 
# 
# 






my $sample_this_db 			= $ARGV[0];
my $number_per_subfamily 		= $ARGV[1];
my $taxonomic_or_random_sampling	= $ARGV[2]; # 1 = per family, 2 = random 
my $outfile				= $ARGV[3];
my $keyfile 				= $ARGV[4];
#my $keyfile = "key_Aug2013_Euka.reduced_according_to.Euka.fas.ng";



my $taxonomic_level_for_sampling	= 4;		#4=family, 3=superfam
$sequence_length_limit = 3000;


if($sample_this_db eq $outfile){die "\nerror same filename for input and output.\n"}













%selected_ids;


if($taxonomic_or_random_sampling == 2)
	{
	print "\n\n  **  SAMPLE_DB.PL  **\n\nyou have selected random sampling....\n";

############
read_ids();#
############

#my @selcted_sequence_ids = keys %selected_ids;
#print "\n$#selcted_sequence_ids\n";

#############
print_ids();# #
#############

	}


if($taxonomic_or_random_sampling == 1)
	{

%which_family;

####################
read_taxkey_file();#
####################

my @array_of_fams = keys %count_family;

print scalar @array_of_fams , " groups at tax level $taxonomic_level_for_sampling\n";


#############
print_ids2();# #
#############

print "\ntotal_printed:$total_printed\n";


	}

print "\n\nFIN.\n";
exit;


####################################################
sample_for_locus_definition($number_per_subfamily);#
####################################################




print "\nscript sample_db.pl finished\n";
exit;



##############################################################################################################
#
#
#
#
#
##############################################################################################################




sub sample_for_locus_definition
{
print "\nsub sample_for_locus_definition\n";

my $resample_number_for_each_family = shift;

my %how_many_for_this_family = ();
my $how_many_currently_sampled = 0;
my $file_as_string = "";
open(IN_FILTER, $sample_this_db) || die "\n\nerror 1454. cant open db:$sample_this_db\n";
while (my $line = <IN_FILTER>)	{$file_as_string .= $line};close(IN_FILTER);
my @all_lines = split />/, $file_as_string;
print scalar @all_lines , " seqs in file $sample_this_db\n";

fisher_yates_shuffle( \@all_lines );

open(OUT, ">$outfile") || die "\n\nerror 3270:$sample_this_db\n";
my $count_orinted=0;
for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];
	if($line =~ /^(.+)/)
		{
		my $speciesid = $1;	#print "$speciesid\n";
		$line =~ s/^.+\n//;
		$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
		$line =~ s/[\s\t]//g;


		if ($taxonomic_or_random_sampling == 1)
			{

			if($speciesid =~ /^([^7_$taxonomic_level_for_sampling]+[a-zA-Z]$taxonomic_level_for_sampling).+[_]/)
				{
				my $current_family = $1;

				unless($speciesid =~ /^([^$taxonomic_level_for_sampling]+BOLD.*$taxonomic_level_for_sampling).+[_]/)
					{
					unless($how_many_for_this_family{$current_family}>=$resample_number_for_each_family)
						{print OUT ">$speciesid\n$line\n";$count_orinted++}
						$how_many_for_this_family{$current_family}++;
					}
				}else{
				#print "$speciesid\n";
				}



			}else{
			
			unless($how_many_currently_sampled>=$resample_number_for_each_family)
			{print OUT ">$speciesid\n$line\n";$count_orinted++}
			$how_many_currently_sampled++;
			

			}



		}
	}


close(OUT);

my @allfams = keys %how_many_for_this_family;@allfams = sort @allfams;

open(OUTLOG, ">listoffams") || die "";

foreach my $fam(@allfams){
print OUTLOG "$fam\t$how_many_for_this_family{$fam}\n";
}
close(OUTLOG);

print "\n\ntotal resampled:$count_orinted\n";
open(OUT , ">>annotation_results") || die "\n\nerror 3327\n";
print OUT "\ntotal resampled:$count_orinted\n";
close(OUT);



}



##############################################################################################################
#
#
#
#
#
##############################################################################################################





######################################################################################



	 sub fisher_yates_shuffle 	# http://perldoc.perl.org/perlfaq4.html
		{
	 	my $deck = shift;my $i = @$deck;
 		while (--$i) 
			{my $j = int rand ($i+1);
 			@$deck[$i,$j] = @$deck[$j,$i]}
 		}

######################################################################################




sub read_ids
{

my %all_species_ids = ();

open(IN_FILTER, $sample_this_db) || die "\n\nerror 1454. cant open db:$sample_this_db\n";

print "reading db ($sample_this_db) and storing IDs\n";

while (my $line = <IN_FILTER>)	
	{
	if($line =~ /^>(.+)/)
		{
		my $speciesid = $1;	#print "$speciesid\n";
		$speciesid =~ s/\n//;$speciesid =~ s/\r//;
		$all_species_ids{$speciesid}  = "";
		}
	}

my @all_ids = keys %all_species_ids;

print scalar @all_ids , " in input file. randomizing these IDs. then selecting the first $number_per_subfamily of them.\n";

fisher_yates_shuffle( \@all_ids );

my $how_many_currently_sampled = 0;

foreach my $id(@all_ids)
	{
	#print "$id\n";
	unless($how_many_currently_sampled>=$number_per_subfamily)
		{
		$selected_ids{$id}= 1;
	#	print "storing $id\n";
		}
	$how_many_currently_sampled++;
	}




}



##############################################################################################################
#
#
#
#
#
##############################################################################################################



sub print_ids
{

print "reading DB again, printing seqeunces for selected IDs to output.\n";

open(FASTA_IN, $sample_this_db) || die "error 275, Cant open database $database_file.\n";

open(OUT, ">$outfile") || die "\n\nerror 3270:$sample_this_db\n";


	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
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

close OUT;
}


##############################################################################################################
#
#
#
#
#
##############################################################################################################



sub process_entry
{
my $line = shift;
my $current_id = "";


if ($line =~ /^(.+)\n/ )
	{
	$current_id = $1;
	$line =~ s/^.+\n//;
	$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);


	if(exists($selected_ids{$current_id}))
		{
		print OUT ">$current_id\n$line\n";
		}

	}else{
	die "\nerror 5631.\n"
	}







}

##############################################################################################################
#
#
#
#
#
##############################################################################################################



sub read_taxkey_file
{

my $count_genus_annotation=0;my $count_order_annotation=0;my $count_family_annotation=0;
my $count_genus_nonannotation=0;my $count_order_nonannotation=0;my $count_family_nonannotation=0;
my $count_species_level_tobycodes=0;

open(IN1, $keyfile) || die "\nerr 336\n";



while (my $line = <IN1>)
	{
	# find species information:
	if($line =~ /^(\S+)\s(\S+)\s\S+\sspecies\s(.+)/)
		{
		my $tobycode = $1;my $binomial = $2;my $string  = $3;$count_species_level_tobycodes++;
	
#		if($id_format == 2){$tobycode = $binomial}
	
#		my $order = "NA";	if($string =~ /\sorder\:(\w+)\s/){$order = $1;$count_order_annotation++}else{$count_order_nonannotation++}
#my $family = "NA";
		if($taxonomic_level_for_sampling == 4)
			{
			if($string =~ /\sfamily\:(\w+)\s/)
				{
				my $family = $1;$which_family{$binomial} = $family;$count_family{$family}=1;
				}
			}
		if($taxonomic_level_for_sampling == 3)
			{
			if($string =~ /\ssuperfamily\:(\w+)\s/)
				{
				my $family = $1;$which_family{$binomial} = $family;$count_family{$family}=1;
				}
			}



#		my $genus = "NA";	if($string =~ /\sgenus\:(\w+)\s/){$genus = $1;$count_genus_annotation++}else{$count_genus_nonannotation++}

#		unless($genus eq "NA"){$which_genus{$tobycode} = $genus}
#		unless($family eq "NA"){$which_family{$tobycode} = $family}
#		unless($order eq "NA"){$which_order{$tobycode} = $order}
#		unless($genus eq "NA"){$which_parent{$tobycode} = $genus}#
#		unless($genus eq "NA"){unless($family eq "NA"){$which_parent{$genus} = $family}}
#		unless($family eq "NA"){unless($order eq "NA"){$which_parent{$family} = $order}}

		}
	
	
	}


print "\nfile information for $taxkey_file\n";
print "\tspecies level tobycodes:$count_species_level_tobycodes
\n";

close IN1;

}


###############################################################################################################################

###############################################################################################################################


sub print_ids2
{



	open(FASTA_IN, $sample_this_db) || die "error 426 Cant open db:$database_file.\n";

open(OUT, ">$outfile") || die "\n\nerror 3270:$sample_this_db\n";


	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
				########################################
				process_entry2($fasta_entry);#
				########################################
				}
			$fasta_entry = $fasta_line;
			}else{
			$fasta_entry .= $fasta_line;
			}
		};
	close(FASTA_IN);

close OUT;
}


##############################################################################################################
#
#
#
#
#
##############################################################################################################


sub process_entry2
{
my $line = shift;
my $current_id = "";


if ($line =~ /^(.+)\n/ )
	{
	$current_id = $1;
	$line =~ s/^.+\n//;
	$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);

	my $speciesID = $current_id;$speciesID =~ s/_[^_]+$//;
	if(exists($which_family{$speciesID}) && $current_seq_length <= $sequence_length_limit)
		{
		my $seq_is_from_this_family = $which_family{$speciesID};

		if($printed_for_this_family{$seq_is_from_this_family} >= $number_per_subfamily)	
			{}else{

			$printed_for_this_family{$seq_is_from_this_family}++;$total_printed++;
			print OUT ">$current_id\n$line\n";
			#print "species id from DB ($speciesID) found in hash of taxkey IDs\n";
			
			}


		}else{
	#	print "species id from DB ($speciesID) NOT found in hash of taxkey IDs\n";
		
		}


	}else{
	die "\nerror 5631.\n"
	}







}



