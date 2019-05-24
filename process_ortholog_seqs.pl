









$process_misof_orthologs 	= 1;
$misof 				= "/home/douglas/scripted_analyses/insect_TOL_analysis/data/Misof_et_al_Supplementary_Archives/Supplementary_Archive_1/1kite_orthologset_hexapoda/1kite_100taxa_hexapoda_2_HMMer3/1kite_100taxa_hexapoda_2_HMMer3.fa";

if($process_misof_orthologs == 1)
	{

	process_misof_orthologs();
	die "";

	};





@files = (
# aug 2016:
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Daphnia.prot.b",

# 
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Halyomorpha.prot.b",
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Zootermopsis.prot.b",
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Apis.prot.b",
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Tribolium.prot.b",
"/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/Papilio.prot.b"
);


foreach my $file(@files)
{
print "file:$file\n";

	open(FASTA_IN, $file) || die "Cant open $database_file.\n";
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

	# do the last entry!:
	unless(length($fasta_entry)<=2)
		{
		$fasta_entry =~ s/^>//;
		########################################
		process_entry($fasta_entry);#
		########################################
		}



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
			}

		if(exists($store_seqs{$current_id})){print "warning $current_id foudn already\n"};
		$store_seqs{$current_id} = $line;
	#	print OUT1 ">$current_id\n$line";
		}






} #foreach my $file(@files)








$file2="/home/douglas/scripted_analyses/insect_TOL_analysis/orthologs/congruent_orths";

open(IN, $file2) || die "cantoen $file2\n";
open(OUT , ">InsCoreOrth") || die "";

my $orth_numebr = 100000;
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	my @split = split /\t+/ , $line;
		$orth_numebr++;
	foreach my $current_id(@split)
		{
		if($store_seqs{$current_id} =~ /\w/)
			{
			print OUT ">$orth_numebr\n$store_seqs{$current_id}";
			}else{
			print "warnig cant find $current_id\n";
			};
		};


#print $line;


	};


close IN;
close OUT;






########################################################




sub process_misof_orthologs
{

open(MISOF, $misof) || die "\ncant open $misof\n";
open(OUT3 , ">Misof_CO") || die "\nerror\n";
print "opened MISOF orthologs file\n";

while (my $line  = <MISOF>)
	{
	if($line  =~ /^>(.+)/)
		{
		my $id = $1;unless($id =~ s/^([^\|]+)\|.+/$1/){die "\nparse error\n"};

		$store_misof_ortholog_ids{$id} = 1;

		print OUT3 ">$id\n";
		}else{
		print OUT3 $line;
		};	
	};



close MISOF;
close OUT3;

my @keys = keys %store_misof_ortholog_ids;
print "misof has $#keys ortholgos \n";

};








########################################################









