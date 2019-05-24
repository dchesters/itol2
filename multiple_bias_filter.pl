




# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 



$ITOL_version = 2;# /COIevol=1, insect_TOL_analysis=2




#################################################################################################################
#
#
#
#################################################################################################################




if($ITOL_version ==1)
	{

# phylobayes:
# alignments=(insectaNUCL.*.RM1)
my @file_list = split /\n/ , `ls insectaNUCL.*.fa.RM1`;

foreach my $file(@file_list)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
	print "";
	if($file =~ /insectaNUCL.(WP_\d+).clo_pruned.fa/)
		{
		my $id = $1;$pb_orthoogs{$id}++;
		}else{
		die "\nerror 52, cannot extract ortholog name from:$file\n"
		};
	};


# blH:
# alignments=(insectaNUCL*.clo_pruned.fa.RM2)
my @file_list2 = split /\n/ , `ls insectaNUCL.*.fa.RM2`;
foreach my $file(@file_list2)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;$bl_orthoogs{$id}++;
		}else{
		die "\nerror 67, cant extrcat orth name from $file\n"
		};
	};


$mare = "/home/douglas/scripted_analyses/COIevol/corrections/transcriptomes/results/insectaNUCL.charsetfile8_reduced";
# charset gene619_insectaNUCL_WP_412469_clo_pruned_fa = 390 - 817;
# charset gene827_insectaNUCL_WP_412681_clo_pruned_fa = 818 - 945;
# charset gene1328_insectaNUCL_WP_413184_clo_pruned_fa = 946 - 1365;

open(IN, $mare) || die "";
while(my $line = <IN>)
	{
	if($line =~ /insectaNUCL_(WP_\d+)_clo/)
		{$mare_orthologs{$1}++};
	};
close IN;



my @file_list = split /\n/ , `ls insectaNUCL.*.fa`;

foreach my $file(@file_list)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
#	print "file:$file\n";
	if($file =~ /insectaNUCL.(WP_\d+).clo_pruned.fa/)
		{
		my $id = $1;
		if( exists($pb_orthoogs{$id}) && exists($bl_orthoogs{$id}) && exists($mare_orthologs{$1})  	)
			{
			$present++;
			#print "\npresent all 3\n"
			system("cp $file $file.RMmulti")
			};

		}else{
		die "\nerror 102, cant extract ortholog name from:$file\n";
		};
	};


print "
presnt all three:$present
";

}; # if($ITOL_version ==1)





#################################################################################################################
#
#
#
#################################################################################################################





if($ITOL_version == 2)
{
print "\nITOL_version == 2\n";



# phylobayes:
# alignments=(insectaNUCL.*.RM1)
my @file_list = split /\n/ , `ls insectaNUCL.*.fa.RM1`;

foreach my $file(@file_list)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
	print "";
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;$pb_orthoogs{$id}++;
		}else{
		die "\nerror 146, cant exatrct ortholog name from:$file\n"
		};
	};


my @pb_keys = keys %pb_orthoogs;
print "pb_keys:$#pb_keys\n";


# blH:
# alignments=(insectaNUCL*.clo_pruned.fa.RM2)
my @file_list2 = split /\n/ , `ls insectaNUCL.*.fa.RM2`;
foreach my $file(@file_list2)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;$bl_orthoogs{$id}++;
		}else{
		die "\nerror 167 can t extract orth id from $file\n"
		};
	};

my @bl_keys = keys %bl_orthoogs;
print "BL_keys:$#bl_keys\n";



$mare = "/home/douglas/scripted_analyses/insect_TOL_analysis/transcriptomes/insectaNUCL.charsetfile8_reduced";
# charset gene619_insectaNUCL_WP_412469_clo_pruned_fa = 390 - 817;
# charset gene827_insectaNUCL_WP_412681_clo_pruned_fa = 818 - 945;
# charset gene1328_insectaNUCL_WP_413184_clo_pruned_fa = 946 - 1365;

open(IN, $mare) || die "\nerror cant open $mare\n";
while(my $line = <IN>)
	{
	if($line =~ /insectaNUCL_([A-Z0-9]+)_clo/)
		{$mare_orthologs{$1}++};
	};
close IN;

my @mr_keys = keys %mare_orthologs;
print "MR_keys:$#mr_keys\n";
unless($#mr_keys >= 1){die "\nerror 191, failed to parse any mare orths\n"};


# alignemtn qual

my @file_list = split /\n/ , `ls insectaNUCL.*.fa.RM4`;
foreach my $file(@file_list)
	{
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;$al_orthoogs{$id}++;
		}else{
		die "\nerror 197, cant extract ortholog name from:$file\n"
		};
	};
my @al_keys = keys %al_orthoogs;
print "ALN_keys:$#al_keys\n";



# bootstrap 

my @file_list = split /\n/ , `ls insectaNUCL.*.fa.RM5`;
foreach my $file(@file_list)
	{
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;$boot_orthoogs{$id}++;
		}else{die "\nerror cant exatrct orth name from:$file\n"};
	};
my @boot_keys = keys %boot_orthoogs;
print "BOOT_keys:$#boot_keys\n";




# finally find orth presnt in all cases

my @file_list = split /\n/ , `ls insectaNUCL.*.fa`;

# print file for generation of venn digram
open(VENN , ">venn_input.txt") || die"";


foreach my $file(@file_list)
	{
# insectaNUCL.WP_413425.clo_pruned.fa 
#	print "file:$file\n";
	if($file =~ /insectaNUCL.([A-Z0-9]+).clo_pruned.fa/)
		{
		my $id = $1;
		if	( 
			exists($pb_orthoogs{$id}) && exists($bl_orthoogs{$id}) && exists($mare_orthologs{$1})
			&& exists($al_orthoogs{$id}) &&  exists($boot_orthoogs{$id})  	
			)
			{
			$present++;
	#		#print "\npresent all 3\n"
			system("cp $file $file.RMmulti")
			};

		unless(exists($pb_orthoogs{$id})){print VENN "$id\tCompositional Heterogeneity\n"};
		unless(exists($bl_orthoogs{$id})){print VENN "$id\tBranchlength Heterogeneity\n"};
		unless(exists($mare_orthologs{$id})){print VENN "$id\tInformation Content\n"};
		unless(exists($al_orthoogs{$id})){print VENN "$id\tAlignment Quality\n"};
		unless(exists($boot_orthoogs{$id})){print VENN "$id\tAverage Bootstrap\n"};

		}else{
		die "\nerror 252, cant exatrct ortholog name from:$file\n"};
	};


print "
out of $#file_list,
presnt all three:$present
";




close VENN;







}; # if($ITOL_version == 2)



#################################################################################################################
#
#
#
#################################################################################################################









