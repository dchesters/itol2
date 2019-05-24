#!/usr/local/bin/perl


# 
# 
# 
# 
# 
# Script adapted from rm_column_from_hmmalignment and dougseq.pl. 
# It removes columns which are mostly missing data.
# 
# 
# 
# 
# 
# 
#
# change log
# oct2013: option to remove only terminal columns
# jul2014: checks for entries with no remaining data in output
#
# 
# 
# 
# 


$infile 			= $ARGV[0];
$outfile 			= $ARGV[1];
unless($infile =~ /./ && $outfile =~ /./){die "\nerror 30.\n"};


$gene_name 			= "16S_hmmer_gb";
$trim_terminal_gaps_only	= 0;



$subscript_infile = $gene_name .  "_columns_rm";
$remove_columns_from_hmmalignment_logfile = $subscript_infile . ".column_removal_logfile";



# $column_removal_cutoff = 0.98; 

# 0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,


					# high cutoff number = few columns removed
for $column_removal_cutoff(0.9)	# low number=few sites remaining
	{

$entry_counter = 0;
my $current_seq_length;
my @seqs;
my @ids;
my %seq_hash;
my %ids_r;
@site_array =();
$sequence_length =0;
@remove_these =();


print "script name: column_removal_analysis\n";
print "infile:$infile outfile:$outfile\n";
print "settings, cutoff for column removed:$column_removal_cutoff percent gapped\n";



##############
read_fasta();#
##############


if($trim_terminal_gaps_only == 1)
	{

my $start_main_alignment =0;my $start_position=0;
for $i(0..$sequence_length)
	{
	if( ($site_array[$i]/$entry_counter) < $column_removal_cutoff && $start_main_alignment ==0)
		{$start_main_alignment =1;$start_position = $i}
	};

my $end_main_alignment =0;my $end_position=0;
for $k(reverse(0 .. ($sequence_length-1)))
	{
	if( ($site_array[$k]/$entry_counter) < $column_removal_cutoff && $end_main_alignment ==0)
		{$end_main_alignment =1;$end_position = $k;
	#print "$k 0\t"
	}else{
	#print "$k 1\t"	
	}
	};

print "\nstart_position:$start_position end_position:$end_position\n";
$how_many_columns=0;
$how_many_removed=0;

for $i(0..$sequence_length)
	{
	$how_many_columns++;
	if( ($site_array[$i]/$entry_counter) >= $column_removal_cutoff)
		{
		if(  $i <=  $start_position || $i >=  $end_position)
			{
		$remove_these[$i]=1;$how_many_removed++;
		print "0";
			}else{
		$remove_these[$i]=0;
		print "1";

			}
		}else{
		$remove_these[$i]=0;
		print "1";
		};
	
	};


	}else{ # if($trim_terminal_gaps_only == 1)

$how_many_columns=0;
$how_many_removed=0;

for $i(0..$sequence_length)
	{
	$how_many_columns++;
	if( ($site_array[$i]/$entry_counter) >= $column_removal_cutoff)
		{
		$remove_these[$i]=1;$how_many_removed++;
		#print "0";
		}else{
		$remove_these[$i]=0;
		#print "1";
		};
	
	};


	}




close(IN);



###############
write_fasta();#
###############



$remaining = $how_many_columns - $how_many_removed;
print "$how_many_columns columns, of which $how_many_removed removed. $remaining remaining\n";

#open(LOGFILE, ">>$remove_columns_from_hmmalignment_logfile") || die "cant open logfile:$remove_columns_from_hmmalignment_logfile\n";
#$date = localtime time;
#print LOGFILE "RUNNING_SCRIPT:column_removal_analysis.pl $date\n";
#print LOGFILE "infile:$infile outfile:$outfile\nnumber_sequences:$seq_counter total_columns:$how_many_columns ";
#print LOGFILE "columns_removed:$how_many_removed columns_remaining:$remaining\n";
#print LOGFILE "cutoff_for_column_removed:$column_removal_cutoff percent_gapped\n";
#close(LOGFILE);


# system("perl assess.pl $subscript_infile");

	};

print "end of script:column_removal_analysis.pl\n";
die;



##############################################################################


##############################################################################






sub read_fasta
	{

	open(FASTA_IN, $infile) || die "Cant open input:$infile.\n";

	# local variables:

	my $line;
	my $current_id;
	my $current_sequence;
	my $current_name;
	my $current_seq;
	my $i4;

	my $current_char;

	while($line = <FASTA_IN>)
		{ #cb1
		$line =~ s/\n//;
		$line =~ s/\r//;

		if ($line=~/^\s{0,2}>\s{0,2}(.+)$/ )
			{


			if($entry_counter==1)
				{$sequence_length=length($current_sequence);
				for $i4(0..$sequence_length){$site_array[$i4]=0}};

				for $i4(0..$sequence_length)
					{
					$current_char=substr($current_sequence,$i4,1);
					if($current_char=~/[nN\?\-]/){$site_array[$i4]=$site_array[$i4]+1};
					}

			

			if ($entry_counter >= 1)
				{
				$current_seq_length=length($current_sequence);
				$seqs[$entry_counter] = $current_sequence;
				};

			$entry_counter ++;# print "ec:$entry_counter line:$line cid:$current_id\n";
			$current_id = $1;
			$current_id =~ s/(\S+)\s\d{1,4}\s+bp\s*$/$1/; # remove 688 bp from end of id if present
			$ids[$entry_counter] = $current_id;
			$ids_r{$current_id} = 1;
			$current_sequence = "";
			}else{
			$current_sequence = $current_sequence . $line;
			};

		}; #cb1

	$ids[$entry_counter] = $current_id;
	$seqs[$entry_counter] = $current_sequence;
	$ids_r{$current_id} = 1;

	close (FASTA_IN);


	print "$entry_counter entries. length:$current_seq_length\n";


	};



####################################################################################
#
#
#
#
#
####################################################################################





sub write_fasta
	{

	open(FASTA_OUT ,">$outfile") || die "cant open output file:$outfile\n";
	

	for $i(1 .. $entry_counter)
		{
		$current_name = $ids[$i];
		$current_seq = $seqs[$i];
		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}
		

		my $printstring = "";
		for $i(0..$sequence_length)
			{
			my $current_char=substr($current_seq,$i,1);
			if($remove_these[$i]==0)
				{
				$printstring .= $current_char;
				#print FASTA_OUT "$current_char";
				};
			};

		if($printstring =~ /^[n-]+$/i)
			{
			print "current entry ($current_name) now has no characters, not printing \n";
			}else{
		print FASTA_OUT ">$current_name\n";
		print FASTA_OUT "$printstring";
		print FASTA_OUT "\n";
			}
		};
	close(FASTA_OUT);
	

	};







