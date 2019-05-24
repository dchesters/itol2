

my $file = $ARGV[0];
my $job = $ARGV[1];

# 
# 
# 
# 
# 	sumtrees out has redundent format, support given both in normal way, plus in square box:)1.00[&support=1.0],
# 		default for this script is retain regular support
# 
# 
# 
# 
# 
# 
# 2015-11-18: do similar with sumtree output, remove these labels:
# 	s,Campodea_augens[&support=1.0],((Pogo
# 2016-03-30: some other string added to end of queries (probably by raxml) is removed
# 2016-06-24: more strings on sumtrees output when its run on bootstrapped trees
# 2018-03-05: processes output of raxml IC calculation
# 2019-05-02: processes raxml SH-like node supprots
# 2019-05-09: prints out node supports to table that can be read into R
# 
# 
# 

unless($file =~ /./ && $job =~ /\d/){die "\ncommand error\n"};

$basic_process = 1; # if ==1 . just rm QUERY___ and node labels, if file is too big to open and do manually.

$terminal_format = 2; # 1 = $raxID\t$termnial; 2 = 

# $print_support = 0;# ignore



if($job == 4)
	{
	print "\nJOB == 4, processing output of Raxml Internode Certainty.\n";

	# 	internode certainty out tree

	# 	[0.935,0.935]):0.90000000000000002220[0.977,0.977]):0.90000000000000002220[0.959,0.959],Diaphorina_citri

# length[x,y] 
# length   corresponds to the original branch length, x corresponds to the IC score and  y  to the ICA score.


	open(IN, $file) || die "\n\nerror, cant open input file ($file). quitting\n\n";
	open(OUT, ">$file.procd") || die "\nerror 6.\n";
	while (my $line = <IN>)
		{
		if($line =~ /..+/)
			{

			while($line =~ s/\)\:([0-9\.]+)\[(\-*[0-9\.]+)\,(\-*[0-9\.]+)\]/ ")" . $2 /e)
				{
				my $branchlength = $1; my $IC1 = $2;my $IC2 = $3;
				$sum_IC += $IC1; $sum_ICA += $IC2; $count_obervatins++; 
			#	print "branchlength:$branchlength IC:$IC1 ICA:$IC2\n";
				};


			# bracnhlengths of internal nodes were removed while processing node labels,
			# thus shouls also remove bracnh lengths of terminals
			while($line =~ s/\:[\d\.]+//)
				{$terminal_branch_lengths_removed++};

			if($line =~ /(.{30}\[.{30})/)
				{my $remaining = $1; print "\nwarning: $remaining\n"};
			print OUT $line;

			};
		};

my $av_IC = $sum_IC / $count_obervatins;
my $av_ICA = $sum_ICA / $count_obervatins;
print "
terminal_branch_lengths_removed:$terminal_branch_lengths_removed

node scores read:$count_obervatins
	average IC:$av_IC
	average ICA:$av_ICA



";


	close OUT;
	exit;
	};

################################
if($job == 5)
	{
	print "\nJOB == 5, processing output of Raxml SH-like node support.\n";

# 04:0.13536897429145425265):0.05258078540271623330[79]):0.15298495426669769803[100]):0.02337187937838788171[69],Mesopsocidae_GMHLN17715_BD00880966:0.11867037726425760935):0.01951065893646064878[81]):0.02874651955601923575[93]):0

	open(IN, $file) || die "\n\nerror, cant open input file ($file). quitting\n\n";
	open(OUT, ">$file.procd") || die "\nerror 6.\n";
	open(OUT7, ">$file.node_supports_tabulated") || die "\nerror 105.\n";

	while (my $line = <IN>)
		{
		if($line =~ /..+/)
			{

			while($line =~ s/\)\:([0-9\.]+)\[(\-*[0-9\.]+)\]/ ")" . $2 /e)
				{
				my $branchlength = $1; my $SH = $2;
				$sum_SH += $SH; $count_obervatins++; 
				print OUT7 "$count_obervatins\t$SH\n";
			#	print "branchlength:$branchlength IC:$IC1 ICA:$IC2\n";

				if($count_obervatins =~ /0000$/){print "intermin. supports processed so far:$count_obervatins\n"};
				};


			# bracnhlengths of internal nodes were removed while processing node labels,
			# thus shouls also remove bracnh lengths of terminals
			while($line =~ s/\:[\d\.]+//)
				{$terminal_branch_lengths_removed++};

			if($line =~ /(.{30}\[.{30})/)
				{my $remaining = $1; print "\nwarning: $remaining\n"};
			print OUT $line;

			};
		};

print "Total, supports processed:$count_obervatins\n";

my $av_SH = $sum_SH / $count_obervatins;
print "
terminal_branch_lengths_removed:$terminal_branch_lengths_removed

node scores read:$count_obervatins
	average SH:$av_SH

";


	close OUT;close OUT7;
	exit;
	};


################################



if($basic_process == 1)
	{
	print "\njust doing basic processing!\n";
	open(IN, $file) || die "\n\nerror, cant open input file ($file). quitting\n\n";
	open(OUT, ">$file.procd") || die "\nerror 6.\n";
	while (my $line = <IN>)
		{
		if($line =~ /..+/)
			{
			$line =~ s/\[I[0-9]+\]//g;

#			in error:
#			if($print_support == 1)
#				{
#				$line =~ s/\[\&support=([\d\.]+)\]/$1/g; # sumtrees format
#				}else{
				$line =~ s/\[\&support=[\d\.]+\]//g; # sumtrees format
#				};

			$line =~ s/QUERY___//g;

			# noticed by Cong Xu, perhaps added to query ID's in a newer version of raxml?
			$line =~ s/___1([^a-z0-9])/$1/gi;


# _charteceus:0.257339378823[&support=0.430107526882][&length_mean=0.257339378823,length_median=0.230783784668,length_sd=0.214478539524,length_hpd95={1.34750946316e-06,0.898686520689},length_quant_5_95={0.0133875124881,0.529611694176},length_range={1.34750946316e-06,0.929595863936}],(Manti

			$line =~ s/\[\&[^\]]+\]//g;


			print OUT $line;
			}else{
			print OUT $line;
			}
		}
	close IN;close OUT;
	print "\ndone\n";exit;
	};


print "

	This script processes outtree of raxml epa.
	Terminals are changed from 
		'terminal_ID' 'branchlength' '[raxmlID]'
		to 
		'raxmlID' 'branchlength'
	A key file is made (XXX.3) giving terminal ID and raxml ID
	Also internal node labels are changed from 
		'branchlength' '[raxmlID]'
		to 
		'raxmlID' 'branchlength'
	Processed newick file is named XXX.2
\n";


open(IN, $file) || die "\n\nerror, cant open input file ($file). quitting\n\n";
open(OUT, ">$file.2") || die "\nerror 6.\n";
open(OUT3, ">$file.3") || die "\nerror 7.\n";


# perl /home/douglas/usr_scripts/newick/process_raxmlEPA_outtree.pl RAxML_originalLabelledTree.porter2014.raxmlEPA.G-0.05

my $processed=0;

while (my $line = <IN>)
	{
	if($line =~ /..+/)
		{
		# 898249[I17690],Ps
	#	$line =~ s/\[I\d+\]//g;


		if($terminal_format == 1) # 1 = $raxID\t$termnial; 2 = 
		{
		# replace terminal IDs with raxmls ids:
		while($line =~ s/([a-z]+\_[a-z]+)(\:[0-9\.]+)\[(I[0-9]+)\]/$3$2/i)
			{
			my $termnial= $1;my $branchlength=$2; my $raxID = $3;$processed++;
		#	print "termnial:$termnial raxID:$raxID branchlength:$branchlength\n";
			print OUT3 "$raxID\t$termnial\n";
			if($processed =~ /000$/){print "$processed\n"}
			};
		};
		if($terminal_format == 2) # 1 = $raxID\t$termnial; 2 = 
		{
		# replace terminal IDs 
		while($line =~ s/([a-z]+\_[a-z]+)(\:[0-9\.]+)\[(I[0-9]+)\]/$1.$3$2/i)
			{
			my $termnial= $1;my $branchlength=$2; my $raxID = $3;$processed++;
		#	print "termnial:$termnial raxID:$raxID branchlength:$branchlength\n";
			print OUT3 "$raxID\t$termnial\n";
			if($processed =~ /000$/){print "$processed\n"}
			};
		};


		# format node labels so can be read by tree veiwer
		$line =~ s/(\:[0-9\.]+)\[(I[0-9]+)\]/$2$1/g;
		$line =~ s/QUERY___//g;

		print OUT $line;
		
		}else{
		print OUT $line;
		
		}

	}


close IN;
close OUT;
close OUT3;


