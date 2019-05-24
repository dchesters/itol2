



# insectaNUCL.EOG5R228R.clo_pruned.fa.RMmulti




my @file_list = split /\n/ , `ls /home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/insectaNUCL.*.clo_pruned.fa.RMmulti`;

print "file_list:$#file_list\n";


foreach my $file(@file_list)
	{
#	print "file:$file\n";
	my $orth_ID = $file;

# file:/home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/insectaNUCL.EOG5ZGMV9.clo_pruned.fa.RMmulti

	$orth_ID =~ s/.+insectaNUCL\.([^\.]+)\.clo_pruned.+/$1/;
#	print "\torth_ID:$orth_ID\n";


# ls /home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/RAxML_result.insectaNUCL.*.clo_pruned.fa.sed.phy

# /home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/RAxML_result.insectaNUCL.EOG5ZPCC0.clo_pruned.fa.sed.phy

	my $command = "cat "  . 
	"/home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/RAxML_result.insectaNUCL." . 
		$orth_ID . ".clo_pruned.fa.sed.phy" . " >> filtered_genetrees";

	print "$command\n";
	system($command);



	};






