


# 
# 
# 
# 
# perl ~/usr_scripts/remove_fasta_entries.pl inv.fas.parsed.ng.rr.new misidentified_entries
# 
#	format in remove_file, just a list:
#
# 	if($line =~ /(.+)/)
#		{my $rmID = $1;
# 
# change log
# 14 june 2014: in addition to looking for remove entries as full fasta id, 
#	option to look using a regex, useful if you want to remove all members of a species without specifying accessions
# 14 jan 2015: also retain fasta entries
# 05 Oct 2015: retain fasta entries works if user uses list of regex's
# 
# 
# 
# 



$database_file 		= $ARGV[0];
$rm_file 		= $ARGV[1];
$out_file 		= $ARGV[2];

	
$regex 			= 0;	# default 0.
$remove_or_retain	= 2;# 1=remove. 2=retain
$blastclust_list 	= 0; # default = 0; 1= select one per OTU




##############################################################################################################



unless($database_file =~ /[\d\w]/ && $rm_file =~ /[\d\w]/ && $out_file =~ /[\d\w]/  )
	{
	die "\ncommand error, something missing.\n"
	};




open(FASTA_IN, $rm_file) || die "Cant open $rm_file.\n";
print "opened $rm_file\n";
while(my $line= <FASTA_IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
#	print "line:$line\n";

	if($blastclust_list == 1){$line =~ s/^(\S+)\s.*$/$1/};

	if($line =~ /(.+)/)
		{
		my $rmID = $1; 
		my $lines_parsed++;
		if($lines_parsed =~ /000$/){print "lines_parsed:$lines_parsed\trmID:$rmID\n"};
	#	 print "ID:($rmID)\n";

		#12S	Hylemya_variata_FJ025383
		#12S	Thambemyia_pagdeni_FJ808133

		$search_IDs{$rmID} = 1;
		}
	}
close(FASTA_IN);


my @array_of_search_IDs = keys %search_IDs;
@array_of_search_IDs = sort @array_of_search_IDs;


print scalar @array_of_search_IDs , "entries parsed, to be looked for in db for removal\n";

#die;

######################################################################################





	open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";
	open(OUT1, ">$out_file") || die"";
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

	unless(length($fasta_entry)<=2)
		{
		$fasta_entry =~ s/^>//;
		########################################
		process_entry($fasta_entry);#
		########################################
		}

	
close OUT1;

print "\nbeen removed from file:$rm
retained:$rt
\n\n";



###################################################################################################################

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

#print "current_id:$current_id\n";
	
if($regex == 1)
	{
	
	my $current_entry_matches_user_string=0;
	foreach my $search_regex(@array_of_search_IDs)
		{
		if($current_id =~ /$search_regex/)
			{
			#print "found entry to be removed ($current_id), matches remove regex ($search_regex)\n";
			$current_entry_matches_user_string =1;
			}
		}

	if($current_entry_matches_user_string == 1)
		{
		if($remove_or_retain == 1)# 1=remove. 2=retain
			{
			# current sequence matches string specified by user to be removed
			$rm++;
			}else{	
			# ... to be retained
			$rt++;print OUT1 ">$current_id\n$line"
			}
		
		}else{

		if($remove_or_retain == 1)
			{
			# current seq ID does not match string in users list of items to be removed
			$rt++;print OUT1 ">$current_id\n$line";
			}else{
			# ... to be retained
			$rm++;
			}
		}

	}else{	

#my $current_idCOPY = $current_id;
#$current_idCOPY =~ s/.+_(\d+)$/$1/;
#if($rm_these{$current_idCOPY} == 1)	{$rm_these{$current_id}=1};

	if($search_IDs{$current_id} == 1)
		{
		if($remove_or_retain == 1)
			{
			$rm++;
			$lines_parsed2++;
			if($lines_parsed2 =~ /000$/){print "will rm $current_id\n"};
		
			}else{
			#print "will retain $current_id\n";
			$rt++;print OUT1 ">$current_id\n$line";
			}
		}else{

		if($remove_or_retain == 1)
			{
			#print "will retain $current_id\n";
			$rt++;print OUT1 ">$current_id\n$line";
			}else{
			$rm++;#print "will rm $current_id\n";
			};

		}

	}


}

###################################################################################################################











