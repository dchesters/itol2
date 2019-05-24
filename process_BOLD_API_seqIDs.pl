



# 
# 
# 
# 2018-05-30: a few different accession formats encountered for plant barcodes.
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
#############################################################################################3



$in = $ARGV[0];
$out = $ARGV[1];

unless($in =~ /\w/ && $out =~ /\w/ ){die "\ncommand error.\n"};

open(IN, $in) || die "\nerror\n";
open(OUT, ">$out") || die "";


#>ASHYC4602-10|chalJanzen01 Janzen36|COI-5P
#>BKSTO512-11|Alloperla serrata|COI-5P
#>BLPAF695-07|Kakopoda stygia|COI-5P|JQ571001
#>BLPCK609-08|Hemiceras vecina|COI-5P|JQ556588
#>BLPEF264-12|Lepidoptera|COI-5P

$fake_accession_number = 0;


while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;	

	my $current_accession = "BD00$fake_accession_number";
	my $BIN; my $tax;

	if($line =~ /^>([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)/)
		{
		$BIN = $1; $tax = $2; $accession = $4;$accession =~ s/^(NC)_(\d+)$/$1$2/;
	#	print "line:$line\n";
	#	print "\tBIN:$BIN tax:$tax accession:$accession\n";

		$accession =~ s/AY038650AY010958/AY010958/; # jeez
		$accession =~ s/^([A-Z]{2}[0-9]{6})[A-Z]{2}[0-9]{6}$/$1/;

		if($accession =~ /^[A-Z]{1,2}\d{5,6}$/)
			{
			# great, has regular ncbi accession
			$current_accession = $accession;
			}elsif($accession =~ /^[A-Z]{1,2}\d{5,6}\.\d$/) # KM003292.1
			{
			$current_accession = $accession;
			}elsif($accession =~ /^[A-Z]{1,2}_\d{5,6}$/) # AC_000188
			{
			$current_accession = $accession;
			}elsif($accession =~ /Pending|WITHDRAWN|SUPPRESSED/)
			{
			# GU690135-SUPPRESSED
			$fake_accession_number++;
			}else{



			die "\nunexpected accession:$accession ... quitting.\n"


			};


		}elsif($line =~ /^>([^\|]+)\|([^\|]+)\|([^\|]+)/)
		{
		$BIN = $1; $tax = $2;	$fake_accession_number++;	
		}elsif($line =~ /^>/){die "\nhuh.\n"};

	$BIN =~ s/[^\w\d]//g;
	my $new_species_ID = "NA";
	if($line =~ /^>/)
		{

		$tax =~ s/ sp\. | nr\. | cf\. / /;
		$tax =~ s/ BOLD\:*(\S+)/ $1/; # CecidInt38 BOLD:ACC5872

		$tax =~ s/^([A-Z][a-z]+) (\d) (\S)/$1$2$3/;  # Muscina 1 AKR
		$tax =~ s/\-//g;

		if($tax =~ /^([A-Z][a-z]+) ([a-z]+) [a-z]+$/)
			{
			$new_species_ID = $1 . "_" . $2;
			
			}elsif($tax =~ /^([A-Z][a-z]+) ([a-z]+)$/)
			{
			$new_species_ID = $1 . "_" . $2;
			}elsif($tax =~ /^[A-Z][a-z]+$/)
			{	
			$new_species_ID = $tax . "_" . $BIN;
			}elsif($tax =~ /^([A-Z][a-z]+) ([a-z]+) \S+$/)
			{
			$new_species_ID = $1 . "_" . $2;

			}elsif($tax =~ /^([A-Z][a-z]+) (.+)/)
			{
			my $a1 = $1; my $rest = $2;$rest =~ s/[^\w\d]//g;
			$new_species_ID = $a1 . "_" . $rest;
			}
			else{
			$tax =~ s/[^\w\d]//g;
			$new_species_ID = $tax;
		#	print "crap tax:$tax\n";
			};
		
		print OUT ">" . $new_species_ID . "_" . $current_accession . "\n";

		$entires++;
		}else{   #  if($line =~ /^>/)

		print OUT "$line\n";
		};




	};


close IN;
close OUT;

print "\n$entires processed\nFIN\n\n";


