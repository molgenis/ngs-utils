#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Text::CSV;
use Spreadsheet::ParseExcel;
use Log::Log4perl qw(:easy);

my %log_levels = (
	'ALL'   => $ALL,
	'TRACE' => $TRACE,
	'DEBUG' => $DEBUG,
	'INFO'  => $INFO,
	'WARN'  => $WARN,
	'ERROR' => $ERROR,
	'FATAL' => $FATAL,
	'OFF'   => $OFF,
);

#
# Get options.
#
my %opts;
Getopt::Std::getopts('d:o:t:i:l:u:f:', \%opts);

my $dropsence	= $opts{'d'};
my $output		= $opts{'o'};
my $taqman		= $opts{'t'};
my $iplex		= $opts{'i'};
my $log_level	= $opts{'l'};
my $inputFolder = $opts{'f'};
my $UGT1A1		= $opts{'u'};

$log_level = (defined($log_level) ? $log_level : 'INFO');

# Reset log level to default if user specified illegal log level.
$log_level = (
	defined($log_levels{$log_level})
	? $log_levels{$log_level}
	: $log_levels{'INFO'});

#Log::Log4perl->init('log4perl.properties');
Log::Log4perl->easy_init(
	{
		level  => $log_level,
		file   => "STDOUT",
		layout => '%d L:%L %p> %m%n'
	},
);
my $logger = Log::Log4perl::get_logger();


#
# Global vars
#
my @verrichtingen = ('CYP1A2','CYP2B6','CYP2C19','CYP2C9','CYP2D6','CYP3A4','CYP3A5','DPYD','F5','HLA-B','MTHFR','SLCO1B1','TPMT','UGT1A1','VKORC1');
my $dropsenceError = "DNA isolation failed";
my $iPLEXUnexpected = "DNA onderzoek niet conclusief";
my $IPLEXEmpty = 'Unknown haplotype'; 
my $iPLEXManualCheck = "Graag handmatig naar kijken.";
my $PGX = "PGX";

#
# Check user input.
#
unless (defined($output)) {
	_Usage();
	exit(1);
}
my $outputDir = dirname($output);

#
# main 
#

#
# Always check Dropsence input. And store failed is in global array.
#

my @verrichtingList = ('CYP1A2','CYP2B6','CYP2C19','CYP2C9','CYP3A4','CYP3A5','DPYD','F5','HLA-B','MTHFR','SLCO1B1','TPMT','VKORC1');
my @verrichtingUGT1A1 = ('UGT1A1');
my @verrichtingCYP2D6 = ('CYP2D6');
my %samplesPassed;

#inputFolder only
if (defined($inputFolder) && defined($output)){
	$logger->warn('Parsing input folder nog niet geimplementeerd.');
	exit(1);
	$iplex  = $inputFolder . 'PGxReport_combined_voorbeeld.csv';
	$taqman = $inputFolder . 'Taqman_voorbeeld.xls';
	$UGT1A1 = $inputFolder . 'UGT1A1_voorbeeld.txt';
	$dropsence = $inputFolder . 'Dropsense_voorbeeld.xls';
	
	&_convertExcelToCSV($taqman);
	&_convertExcelToCSV($dropsence);
	%samplesPassed = %{&_checkDropsence($dropsence)};
	
	#voeg extra verrichting toe
	push(@verrichtingList , @verrichtingCYP2D6);
	push(@verrichtingList , @verrichtingUGT1A1);
	&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
	exit(0);
} 
# check only iplex data, without verrichtingen Tagman and UGT1A1
if (!defined($dropsence) && defined($iplex) && !defined($taqman)&& !defined($UGT1A1)){
	# situatie met alleen een iplexfile
	#_Usage();
	#exit(1);
	$logger->info('alleen iplex');
	&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}
# check only UGT1A1 file, without iplex verrichtingen
if (!defined($dropsence) && defined($UGT1A1) && !defined($iplex)){
	_Usage();
	exit(1);
	#$logger->info('alleen UGT1A1');
	#&_parseUGT1A1($UGT1A1,$output);
}

# check iplex including Taqman file, without UGT1A1 verrichting
if (!defined($dropsence) && defined($iplex) && defined($taqman)&& !defined($UGT1A1)){
	_Usage();
	exit(1);
	#$logger->info('alleen iplex en taqman');
	#add CYP2D6 verrichting for Taqman
	#&_convertExcelToCSV($taqman);
	#my $taqmanCSV = basename($taqman,  ".xls") . ".csv";
	#push(@verrichtingList, @verrichtingCYP2D6);
	#&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}

# check iplex + UGT1A1 file, without Tagman verrichting
if (!defined($dropsence) && defined($iplex) && defined($UGT1A1) && !defined($taqman)){
	_Usage();
	exit(1);
	#voeg extra verrichting toe
	#$logger->info('alleen iplex en UGT1A1');
	#push(@verrichtingList , @verrichtingUGT1A1);
	#&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}


# check iplex + Tagman and UGT1A1 file 
if (!defined($dropsence) && defined($iplex) && defined($UGT1A1) && defined($taqman)){
	$logger->info('Parsing Iplex, Taqman and UGT1A1 files.');
	&_convertExcelToCSV($taqman);
	#voeg extra verrichting toe
	push(@verrichtingList , @verrichtingCYP2D6);
	push(@verrichtingList , @verrichtingUGT1A1);
	&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
	exit(0);
}

# only dropsence test
if (defined( -f $dropsence) && !defined($iplex) && !defined($UGT1A1) && !defined($taqman)) {
	&_convertExcelToCSV($dropsence);
	&_checkDropsence($dropsence);
	exit(0);
}


#
## Functions
#

# gets dropsenceCSV path, and checks if 'Concentration (ng/ul)', 'A260/A280', 'A260/A230' are within boundaries
sub _checkDropsence{
my $dropsenceFile = shift;
my $dropsenceCSV = basename($dropsence,  ".xls") . ".csv";
unlink($output);

$logger->info('Checking samples in Dropsense file.');

#
# Retrieve header from CSV file
#
my $input_fh;
open($input_fh, '<', $dropsenceCSV) or die "Can't open Dropsense file.";

my $csv = Text::CSV->new();
my @column_labels;
while (<$input_fh>) {
    next unless ($. == 1); # header only
    if ($csv->parse($_)) {
        @column_labels = $csv->fields();
    } 
 }

close($input_fh);
my %header = (map { $column_labels[$_] => $_ } (0 .. $#column_labels)); 


# extract these columns:
my @dropsence_column_names = ('Sample name','Concentration (ng/ul)','A260/A280','A260/A230');
my ($sampleID,$concentration,$OD280,$OD230);
my %experiment_failed;
my %experiment_pass;
my $fail =0;

open($input_fh, '<', $dropsenceCSV);	

while (<$input_fh>) {	
	if ($csv->parse($_)) {
		
		my @fields = $csv->fields();
		$sampleID			= $fields[$header{'Sample name'}];
		$concentration		= $fields[$header{'Concentration (ng/ul)'}];
		$OD280				= $fields[$header{'A260/A280'}];
		$OD230				= $fields[$header{'A260/A230'}];

		#skip header
		if ($. == 1){
			next;
		}
		
		#Concentration should be above 10 ng/ul
		#The OD 260/280 ratio should be between 1.8 and 2.0
		#The OD 260/230 ratio should be above 1.5
		if (($concentration < 10) || ($OD280 < 1.8 || $OD280 > 2) || ($OD230 < 1.5)){
			$fail++;
			#add failed samples
			$experiment_failed{$sampleID} = 1;
			
			# error wat uitbreiden
			my $log_message = sprintf("$sampleID failed QC: Concentration (ng/ul): $concentration, A260/A280: $OD280, A260/A230: $OD230,\r\n" );
			$logger->warn($log_message);
			
			# error foreach possible verrichting;
			foreach my $verrichting (@verrichtingen){
				&_printGLIMP($sampleID,$verrichting,$dropsenceError,$output);
			}
		} else{
			$experiment_pass{$sampleID} = 1;
		}
		
	}
}

if ($fail == 0){
		$logger->info("All samples passed Dropsence test");
	}	
close($input_fh);
return \%experiment_pass;
}

#
# adds a verrichting row to outputfile for given sampleid 
sub _printGLIMP{
	my ($sampleID, $verrichting, $result, $output) = @_;
	my $outputDir = dirname($output);
	my $fileName = "/input_" . basename($output);
	my $output_fh; 
	my $outputFile = $outputDir . $fileName; 
	
	my $delimiter = ";";
	open($output_fh, '>>', $outputFile) or die "Can't open outputfile file.";;
	my $line = $sampleID . "$delimiter$verrichting$delimiter$delimiter" . $result ."$delimiter". "$delimiter". "$delimiter$PGX\n" ; 
	$logger->debug("RESULTLINE: $line");
	if (!$result eq ''){
		print($output_fh $line);
	} else{
		my $log_message = "$sampleID, $verrichting is empty";
		$logger->warn($log_message);
	}
	close($output_fh);
} 

#
# 
sub _createGlimpFile{
my ($iplex, $taqman, $UGT1A1, $output ,@verrichtingList) = @_;

#
# Retrieve header from CSV file
#
my $input_fh;
open($input_fh, '<', $iplex) or die "Can't open Iplex file.";;

my $csv = Text::CSV->new();
my @column_labels;
while (<$input_fh>) {
    next unless ($. == 1); # header only
    if ($csv->parse($_)) {
        @column_labels = $csv->fields();
    } 
 }

close($input_fh);
my %header = (map { $column_labels[$_] => $_ } (0 .. $#column_labels)); 

# extract these columns:
my ($sampleID,$gene);
open($input_fh, '<', $iplex)  or die "Can't open Iplex file.";	

while (<$input_fh>) {
	if ($csv->parse($_)) {
		if ($. == 1) {
				next;
			}
		my @fields= $csv->fields();
		$sampleID	= $fields[$header{'Sample'}];
		
		
		#select only A sampleIds, and merge B and C later.	
		if ($sampleID =~ /A$/) {
			
			my $sampleShort = $sampleID;
			chop($sampleShort);
			
			# Voor elke rij worden de aanwezige verrichting geparsed, geconverteerd en weggeschreven naar de outputfile.			
			foreach my $verrichting(@verrichtingList) {
				if ($verrichting eq 'CYP2D6') {
					
					my $BesteVerrichting = &_ABC_check($sampleID,$verrichting,$iplex);
					$logger->debug("$sampleID best match voor verrichting: $verrichting = $BesteVerrichting \n");
												
					$sampleID	= $fields[$header{'Sample'}];
					# get new CYP2D6 string;
					$gene		= &_convertCYP2D6( $sampleID, $BesteVerrichting, &_parseTaqman($taqman,$sampleShort));
										
					&_printGLIMP($sampleShort, $verrichting, $gene, $output);
			
				}elsif ($verrichting eq 'DPYD') {
					my $BesteVerrichting = &_ABC_check($sampleID,$verrichting,$iplex);
					my $DPYD_rs56038477 = &_ABC_check($sampleID,'DPYD (rs56038477)',$iplex);
					my $DPYD_rs67376798 = &_ABC_check($sampleID,'DPYD (rs67376798)',$iplex);
					$logger->debug("$sampleID best match voor verrichting: $verrichting zijn $BesteVerrichting , $DPYD_rs56038477 , $DPYD_rs67376798 \n");
					$sampleID	= $fields[$header{'Sample'}];
					#get correct DPYD;
					#$gene		= &_convertDPYD($fields[$header{$verrichting}],$fields[$header{'DPYD (rs56038477)'}],$fields[$header{'DPYD (rs67376798)'}]);
					$gene		= &_convertDPYD($BesteVerrichting,$DPYD_rs56038477,$DPYD_rs67376798);
					
					&_printGLIMP($sampleShort, $verrichting, $gene, $output);
							
				}elsif ($verrichting eq 'UGT1A1') {
				
					$sampleID	= $fields[$header{'Sample'}];
					
					#get correct UGT1A1;
					$gene = &_convertUGT1A1($UGT1A1,$sampleShort);
				
					&_printGLIMP($sampleShort, $verrichting, $gene, $output);
							
				}elsif ($verrichting eq 'TPMT') {
					my $BesteVerrichting = &_ABC_check($sampleID,$verrichting,$iplex);
					$logger->debug("$sampleID best match voor verrichting: $verrichting = $BesteVerrichting \n");
					$sampleID	= $fields[$header{'Sample'}];
					#get correct TPMT;
					$gene = &_convertTPMT($BesteVerrichting,$sampleID);
					#$gene = &_convertTPMT($fields[$header{$verrichting}],$sampleID);
					
					&_printGLIMP($sampleShort, $verrichting, $gene, $output);
							
				}
				else {
					my $BesteVerrichting = &_ABC_check($sampleID,$verrichting,$iplex);
					$logger->debug("$sampleID best match voor verrichting: $verrichting = $BesteVerrichting \n");
					$sampleID	= $fields[$header{'Sample'}];
					$gene		= $BesteVerrichting;
					$logger->debug("$sampleShort, $verrichting, $gene, $output \n");
				
					if ($gene eq 'Unknown haplotype') {
						$gene = $iPLEXUnexpected;
					}
					&_printGLIMP($sampleShort, $verrichting, $gene, $output);
					
					}		
				}
			} # close A if loop
			else {$logger->debug("SKIP: $sampleID\n")};
		}
	}
}

#
## functions.
#

#
# Gets 3 DPYD values, and returns converted results back.
sub _convertDPYD{
my $DPYD1 = shift;
my $DPYD2 = shift;
my $DPYD3 = shift;
my $newDPYD = '';
my $REF = '';
my $ALT = '';


if ($DPYD1 eq 'Unknown haplotype') {
	$newDPYD=$iPLEXUnexpected;
	}else{
	
	my ($GEN1, $GEN2) = split /\s*\/\s*/, $DPYD1;

	# Both wildtype? return $DPYD1
	if ($DPYD2 eq 'WT/WT' && $DPYD3 eq 'WT/WT'){
		$newDPYD = $DPYD1;
	}elsif ($GEN1 ne '*1' && $GEN2 ne '*1' ) {
		$newDPYD=$DPYD1;
	}else {
		my ($rs56038477, $RSWT1) = split /\s*\/\s*/, $DPYD2;
		my ($rs67376798, $RS2WT) = split /\s*\/\s*/, $DPYD3;
	
		if ($rs56038477 eq 'WT' || $rs67376798 eq 'WT') {
		
			if($GEN1 eq '*1') {
				$REF = $GEN2;
			}else {
				$REF = $GEN1;
			}
		
			if($rs56038477 eq 'WT') {
				$ALT = $rs67376798;
			}else {
				$ALT = $rs56038477;	
			}
			$newDPYD=$REF.'/'.$ALT;
		
		}else{
			if ($rs56038477 ne 'WT' || $rs67376798 ne 'WT'){
				if($GEN1 eq '*1'){
					$REF = $GEN2;
				}else {
					$REF = $GEN1;
				}
				$newDPYD=$REF.'/'.$rs56038477.' OR '. $REF.'/'.$rs67376798;
			}else{
				my $log_message = "UNKNOWN: $DPYD1 , $DPYD2 ,$DPYD3 " ;
				$logger->error($log_message);
				$newDPYD=$iPLEXUnexpected;
			}
		} 
	}
}
return $newDPYD;
}


sub _convertCYP2D6{
my $sampleID =shift;
my $CYP2D6 = shift;
my $taqmanCN = shift;
my $newCYP2D6 = '';

# Taqman opzoek hash
my %Exon_9 =(
	"*1" => 1,
	"*3" => 1,
	"*4" => 1,
	"*4N" => 0,
	"*4M" => 1,
	"*5" => 0,
	"*6" => 1,
	"*7" => 1,
	"*8" => 1,
	"*9" => 1,
	"*10" => 1,
	"*36" => 0,
	"*37" => 1,
	"*47" => 1,
	"*49" => 1,
	"*52" => 1,
	"*54" => 1,
	"*56B" => 1,
	"*57" => 0,
	"*65" => 1,
	"*72" => 1,
	"*87" => 1,
	"*94" => 1,
	"*95" => 1,
	"*100" => 1,
	"*101" => 1,
	"*12" => 1,
	"*14A" => 1,
	"*14B" => 1,
	"*17" => 1,
	"*40" => 1,
	"*58" => 1,
	"*29" => 1,
	"*70" => 1,
	"*41" => 1,
	"*91" => 1,
	"*64" => 1,
	"*69" => 1,
	"*82" => 1,
	"*109" => 1,
	"*80" => 0,
);
# Taqman opzoek hash
my %Intron_6_2 =(
	"*1" => 1,
	"*3" => 1,
	"*4" => 1,
	"*4N" => 1,
	"*4M" => 1,
	"*5" => 0,
	"*6" => 1,
	"*7" => 1,
	"*8" => 1,
	"*9" => 1,
	"*10" => 1,
	"*36" => 1,
	"*37" => 1,
	"*47" => 1,
	"*49" => 1,
	"*52" => 1,
	"*54" => 1,
	"*56B" => 1,
	"*57" => 1,
	"*65" => 1,
	"*72" => 1,
	"*87" => 1,
	"*94" => 1,
	"*95" => 1,
	"*100" => 1,
	"*101" => 1,
	"*12" => 1,
	"*14A" => 1,
	"*14B" => 1,
	"*17" => 1,
	"*40" => 1,
	"*58" => 1,
	"*29" => 1,
	"*70" => 1,
	"*41" => 1,
	"*91" => 1,
	"*64" => 1,
	"*69" => 1,
	"*82" => 1,
	"*109" => 1,
	"*80" => 0,
);

my %alternatieven =(
	"*4" => "*4N",
	);


if ($taqmanCN eq $iPLEXUnexpected | $CYP2D6 eq $IPLEXEmpty ){
	my $log_message = "$sampleID not in Taqman file, or CopyNumber is missing: $taqmanCN." ;
	$logger->warn($log_message);
	return "TaqmanError:$taqmanCN";
} else{
	
	my @options = split / OR /, $CYP2D6;
	foreach my $option (@options){
		  	
		my ($var1,$var2) = split /\//, $option;

		$logger->debug("CYP2D6 sample: $sampleID: $option and tagman: $taqmanCN\n");
		my @allel1 = split /;/, $var1;
		my @allel2 = split /;/, $var2;
	
		# add alternatives
		foreach my $key(@allel1){
			if (exists $alternatieven{$key}){
				push (@allel1, $alternatieven{$key});
			}
		}
		foreach my $key(@allel2){
			if (exists $alternatieven{$key}){
				push (@allel2, $alternatieven{$key});
			}
		}
	
		#
		# for all possible combinations off allel 1 and 2, calculate the predicted copyNumbers. That compare with the Taqman results and report the correct options. 
		foreach my $allelOne (@allel1) {
			foreach my $allelTwo (@allel2){
				my $voorspeld = (($Exon_9{$allelOne} + $Exon_9{$allelTwo}) . ',' . ($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo}) . ',' . ($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo}));
			
				if ( $voorspeld eq $taqmanCN){
					my $log_message	= "sample: $sampleID predicted allels eqal: $allelOne/$allelTwo: ";
					$log_message	.= ($Exon_9{$allelOne} + $Exon_9{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo})." Taqman called $taqmanCN";
					$logger->info($log_message);	
					$newCYP2D6 .= $allelOne.'/'.$allelTwo . ' OR ';
				}
				else{
					my $log_message	= "sample: $sampleID predicted allels differ: $allelOne/$allelTwo: ";
					$log_message	.= ($Exon_9{$allelOne} + $Exon_9{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo})." Taqman called $taqmanCN";
					$logger->warn($log_message);
					}
				
				}
			
			}
		}
	}

# remove extra OR at the back of string.
 $newCYP2D6 = substr($newCYP2D6, 0, -4);
 
 if ($newCYP2D6 eq ''){
	 $newCYP2D6 = $iPLEXManualCheck;
 }

return  $newCYP2D6;

}

sub _convertUGT1A1{
my $UGT1A1 = shift;
my $sampleID = shift;
my $UGT1A1line = '';
my $CorrectionBase = 253; # basis voor afronding.  $CorrectionBase(253) - Controlewaard = afronding§($correctionFactor)
my $correctionFactor = '';

# De UGT1A1 tandem repeat analyse bevat een vertaalslag van fragmentlengte naar haplotype zoals hieronder weergegeven.
# (TA)5: 251 bp (haplotype *36)
# (TA)6: 253 bp (wild type or haplotype *1)
# (TA)7: 255 bp (haplotype *28)
# (TA)8: 257 bp (haplotype *37)

my %UGT1A1Table = (
    "251"  => "*36",
    "253" => "*1",
    "255"  => "*28",
    "257"  => "*37",
);


my $fh;
open($fh, '<', $UGT1A1) or die "Can't open UGT1A1 file.";

# pakt de controle sample en berekend correctieFactor.
while (<$fh>) {
	if( my ( $id, $var1, $var2 ) = $_ =~ m/^(CONTROLE)\t([0-9.]+)\t([0-9.]+).+$/ ) {
		$correctionFactor = $CorrectionBase - $var1;
	}
}
close($fh);

open($fh, '<', $UGT1A1) or die "Can't open UGT1A1 file.";

# pakt beide UGT1A1 lengtes, rond af, haalt de correctieFactor eraf, en return 2 UGT1A1 allele.
while (<$fh>) {
	if( my ( $id, $var1, $var2 ) = $_ =~ m/^([0-9a-zA-z.]*)\t([0-9.]*)\t([0-9.]*)$/ ) {
				$logger->all("TEST:UGT1A1: sampleID: $id, " . &_roundup($var1) . ", " . &_roundup($var2) . ", " . 'correctie factor +'.$correctionFactor);
		if ($sampleID eq $id){
			#check for non existing values.
			if (! exists $UGT1A1Table{(&_roundup($var2)+$correctionFactor)} || ! exists $UGT1A1Table{(&_roundup($var1)+$correctionFactor)}){
				$UGT1A1line = $iPLEXUnexpected;
				$logger->warn("$sampleID: $var1 or $var2 does not exist in %UGT1A1Table");
				#next;
			}
			if ($UGT1A1Table{(&_roundup($var2)+$correctionFactor)} eq "*1"){
				$UGT1A1line = $UGT1A1Table{(&_roundup($var2)+$correctionFactor)}.'/'.$UGT1A1Table{(&_roundup($var1)+$correctionFactor)};
				$logger->info("UGT1A1: sampleID: $id, " . &_roundup($var2) . ", " . &_roundup($var1) . ", " . 'correctie factor +'.$correctionFactor);
			}else{
				$UGT1A1line = $UGT1A1Table{(&_roundup($var1)+$correctionFactor)}.'/'.$UGT1A1Table{(&_roundup($var2)+$correctionFactor)};	
				$logger->info("UGT1A1: sampleID: $id, " . &_roundup($var1) . ", " . &_roundup($var2) . ", " . 'correctie factor +'.$correctionFactor);
			}
		}else{
			next;
		}
		
		
	}
}	
	close($fh);
	return $UGT1A1line;
}

###

###


sub _convertTPMT{
my $TPMT = shift;
my $sampleID = shift;
my $newTPMT = '';

if ($TPMT =~ m/.1..3A.OR..3B..3C/){
	$newTPMT =	'*1/*3A';
} elsif($TPMT eq 'Unknown haplotype'){
	$newTPMT = $iPLEXUnexpected;
}
else{
	$newTPMT =	$TPMT;
}
	return $newTPMT;
}


# Leest Tagman file in en return de CopyNummers voor 3 metingen in 1 string. Indien 1 meting leeg is wordt een error gereturnd.
sub _parseTaqman{
my $taqmanFile = shift;
my $sampleID = shift;
my %Hs00010001_cn = (); # Hs00010001_cn (=exon 9)
my %Hs04502391_cn = (); # Hs04502391_cn (=intron 6)
my %Hs04083572_cn = (); # Hs04083572_cn (=intron 2)
my $taqmanFileCSV = basename($taqmanFile,  ".xls") . ".csv"; #BUG: schrijft output in scriptfolder

# watch out the encoding!
open(my $tagman_fh, '<:utf8', $taqmanFileCSV)
    or die "Can't open $taqmanFileCSV: $!";
    
# skip to first header af sampleblock
my $header = '';
while (<$tagman_fh>) {
    if (/^Sample Name,/) {
        $header = $_;
        last;
    }
}

my $csv = Text::CSV->new
    or die "Text::CSV error: " . Text::CSV->error_diag;

# define column names    
$csv->parse($header);
$csv->column_names([$csv->fields]);

# parse de Kollom CN Predicted tot aan sample ID 'NTC'
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^NTC/){
		last;
	} else{
		$Hs04083572_cn{$row->{'Sample Name'}} = $row->{'CN Predicted'};
    }
}

# skip to next sampleblock
while (<$tagman_fh>) {
    if (/^Sample Name,/) {
         last;
    }
}
# parse de Kollom CN Predicted tot aan sample ID 'NTC'
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^NTC/){
		last;
	} else{
		$Hs00010001_cn{$row->{'Sample Name'}} = $row->{'CN Predicted'};
    }
}
# skip to next sampleblock
while (<$tagman_fh>) {
    if (/^Sample Name,/) {
         last;
    }
}
# parse de Kollom CN Predicted tot aan sample ID 'NTC'
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^NTC/){
		last;
	} else{
		$Hs04502391_cn{$row->{'Sample Name'}} = $row->{'CN Predicted'};
    }
}

$csv->eof or $csv->error_diag;
close $tagman_fh;

#return de 3 metingen als die er  zijn, of andere een error $iPLEXUnexpected.
if (! $Hs00010001_cn{$sampleID} eq '' && ! $Hs04502391_cn{$sampleID} eq ''  && ! $Hs04083572_cn{$sampleID} eq '' ){
	$logger->debug("Tagman return for sample $sampleID: $Hs00010001_cn{$sampleID},$Hs04502391_cn{$sampleID},$Hs04083572_cn{$sampleID}");
	return ("$Hs00010001_cn{$sampleID},$Hs04502391_cn{$sampleID},$Hs04083572_cn{$sampleID}");
} else{
	$logger->debug("Tagman return for sample $sampleID: CN Predicted fields empty. Response: $iPLEXUnexpected ");
	return $iPLEXUnexpected;
	}
}

# geeft per ABC sample rij de meest gemeten waarde terug.
sub _ABC_check{
my $sampleID = shift;
my $verrichting = shift;
my $iplexFile = shift;

chop($sampleID);
my $sampleA = $sampleID . 'A';
my $sampleB = $sampleID . 'B';
my $sampleC = $sampleID . 'C';
my $celA = '';
my $celB = '';
my $celC = '';

$logger->debug("chopped ID:$sampleID, $sampleA $sampleB $sampleC");

# watch out the encoding!
open(my $iplex_fh, '<:utf8', $iplexFile)
    or die "Can't open $iplexFile: $!";
    
# skip to first header af sampleblock
my $header = '';
while (<$iplex_fh>) {
    if (/^Sample,/) {
        $header = $_;
        last;
    }
}

my $csv = Text::CSV->new
    or die "Text::CSV error: " . Text::CSV->error_diag;

# define column names    
$csv->parse($header);
$csv->column_names([$csv->fields]);

# parse the rest
while (my $row = $csv->getline_hr($iplex_fh)) {
	if ($row->{'Sample'} =~ m/^$sampleA/){
		$logger->debug("we hebben een $row->{'Sample'} match: $sampleA. Verrichting $verrichting is: $row->{$verrichting}");
		$celA = $row->{$verrichting};
	} elsif($row->{'Sample'} =~ m/^$sampleB/){
		$logger->debug("we hebben een $row->{'Sample'} match: $sampleB. Verrichting $verrichting is: $row->{$verrichting}");
		$celB = $row->{$verrichting};
	} elsif($row->{'Sample'} =~ m/^$sampleC/){
		$logger->debug("we hebben een $row->{'Sample'} match: $sampleC. Verrichting $verrichting is: $row->{$verrichting}");
		$celC = $row->{$verrichting};
	}else{
		#next;
	}
}

# matches 3 kollomen, en return 1 van 2 matchende
	if ($celA eq $celB | $celA eq $celC){
        $logger->debug("celA $celA match met $celB of $celC");
        return $celA;
    }elsif ($celB eq $celC){
		$logger->debug("celB $celB match met $celC");
		return $celB;
	} else{
		$logger->warn("$celA  $celB $celC matchen niet. return: $iPLEXUnexpected ");
		return $iPLEXUnexpected;
	}
}

#converts .xls to .csv, and outputs in same folder als input .xls
sub _convertExcelToCSV{

my $sourcename = shift;
my $source_excel = new Spreadsheet::ParseExcel;
my $source_book = $source_excel->Parse($sourcename)
    or die "Could not open source Excel file $sourcename: $!";
my $storage_book;

my $filename = basename($sourcename,  ".xls");
my $CSV_fh;
my $filenameOut = $outputDir.'/'.$filename . '.csv';
open($CSV_fh, '>', $filenameOut);
print($CSV_fh '');
open($CSV_fh, '>>', $filenameOut);

foreach my $source_sheet_number (0 .. $source_book->{SheetCount}-1) {
	my $source_sheet = $source_book->{Worksheet}[$source_sheet_number];

	next unless defined $source_sheet->{MaxRow};
	next unless $source_sheet->{MinRow} <= $source_sheet->{MaxRow};
	next unless defined $source_sheet->{MaxCol};
	next unless $source_sheet->{MinCol} <= $source_sheet->{MaxCol};

	foreach my $row_index ($source_sheet->{MinRow} .. $source_sheet->{MaxRow}) {
		foreach my $col_index ($source_sheet->{MinCol} .. $source_sheet->{MaxCol}) {
			my $source_cell = $source_sheet->{Cells}[$row_index][$col_index];
			if ($source_cell && $source_cell->Value) {
				my $regel = $source_cell->Value.",";
				print($CSV_fh $regel);
			} else {
				my $regel = ",";
				print($CSV_fh $regel);
			}
		}
		my $regel = "\n";
		print($CSV_fh $regel);
	}
}
close($CSV_fh);
} 

sub _roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 0.5))
}
 
 
sub _Usage {
	
print <<EOF;
#########################################################################################################
# This script creates a GLIMPs output file using 4 inputfile (Dropsense | Taqman | Iplex | UGT1A1 ).    #
#                                                                                                       #
#########################################################################################################
Usage: perl FGX.pl [options]

Options:
-o outputfile 	 (';' separated .csv/txt)
-d dropsensefile (.xls)
-t taqman	 (.xls)
-i iplex file	 (.csv)
-u UGT1A1 file	 (.txt)
-l log_level	 DEBUG or INFO(default)

#########################################################################################################
EOF

}
 
 
 
 
 
 
 
 
 
 
 
 
 
