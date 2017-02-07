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
Getopt::Std::getopts('d:o:t:i:l:u:', \%opts);

my $dropsence	= $opts{'d'};
my $output		= $opts{'o'};
my $taqman		= $opts{'t'};
my $iplex		=  $opts{'i'};
my $log_level	= $opts{'l'};
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
my $iPLEXUnexpected = "DNA analyse niet conclusief";
my $iPLEXManualCheck = "Graag handmatig naar kijken.";
my $PGX = "PGX";

#
# Check user input.
#
unless (defined($dropsence) && defined($output)) {
	_Usage();
	exit(1);
}
if ($dropsence =~ /^$/ || $output =~ /^$/) {
	_Usage();
	exit(1);
}
if ($dropsence eq $output) {
	_Usage();
	$logger->fatal('Output file is the same as the input file. Select a different file for the output.');
	exit(1);
}
unless (-f $dropsence && -r $dropsence) {
	$logger->fatal('Cannot read/access file: ' . $dropsence);
	exit(1);
}


#
# main 
#

#
# Always check Dropsence input. And store failed is in global array.
#
my %samplesFailed;
&_convertExcelToCSV($dropsence);
%samplesFailed = %{&_checkDropsence($dropsence)};

#if (defined( -f $dropsence) && !defined($iplex) && !defined($UGT1A1)) {
	#&_convertExcelToCSV($dropsence);
	#&_checkDropsence($dropsence);
#}

my @verrichtingList = ('CYP1A2','CYP2B6','CYP2C19','CYP2C9','CYP3A4','CYP3A5','DPYD','F5','HLA-B','MTHFR (rs1801133)','SLCO1B1','TPMT','VKORC1');
my @verrichtingUGT1A1 = ('UGT1A1');
my @verrichtingCYP2D6 = ('CYP2D6');

# check only iplex data, without verrichtingen Tagman and UGT1A1
if (defined($dropsence) && defined($iplex) && !defined($taqman)&& !defined($UGT1A1)){
	$logger->info('alleen iplex');
	# situatie met alleen een iplexfile
	_Usage();
	exit(1);
	#&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}

# check only UGT1A1 file, without iplex verrichtingen
if (defined($dropsence) && defined($UGT1A1) && !defined($iplex)){
	_Usage();
	exit(1);
	#$logger->info('alleen UGT1A1');
	#&_parseUGT1A1($UGT1A1,$output);
}

# check iplex including Taqman file, without UGT1A1 verrichting
if (defined($dropsence) && defined($iplex) && defined($taqman)&& !defined($UGT1A1)){
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
if (defined($dropsence) && defined($iplex) && defined($UGT1A1) && !defined($taqman)){
	_Usage();
	exit(1);
	#voeg extra verrichting toe
	#$logger->info('alleen iplex en UGT1A1');
	#push(@verrichtingList , @verrichtingUGT1A1);
	#&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}

# check iplex + Tagman and UGT1A1 file 
if (defined($dropsence) && defined($iplex) && defined($UGT1A1) && defined($taqman)){
	$logger->info('alle files zijn er');
	&_convertExcelToCSV($taqman);
	#voeg extra verrichting toe
	push(@verrichtingList , @verrichtingCYP2D6);
	push(@verrichtingList , @verrichtingUGT1A1);
	&_createGlimpFile($iplex, $taqman, $UGT1A1, $output ,@verrichtingList);
}

#
# Functions
#

# gets dropsenceCSV path, and checks if 'Concentration (ng/ul)', 'A260/A280', 'A260/A230' are within boundaries
sub _checkDropsence{
my $dropsenceFile = shift;
my $dropsenceCSV = basename($dropsence,  ".xls") . ".csv";
unlink($output);

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
			my $log_message = sprintf($concentration, $concentration, $OD280, $OD230,"\r\n" );
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
return \%experiment_failed;
}

#
# adds a verrichting row to outputfile for given sampleid 
sub _printGLIMP{
	my ($sampleID, $verrichting, $result, $output) = @_;
	my $output_fh;
	open($output_fh, '>>', $output) or die "Can't open outputfile file.";;
	my $line = $sampleID . ";$verrichting;;" . $result .";". ";". ";$PGX\n" ; 
	if (!$result eq ''){
		print($output_fh $line);
	} else{
		my $log_message = "$sampleID, $verrichting is empty";
		$logger->info($log_message);
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
		#
		## skip samples that fail QC. Maybe only good IDs?
		#
		my @fields= $csv->fields();
		$sampleID	= $fields[$header{'Sample'}];
		if (exists $samplesFailed{$sampleID}){
			$logger->info("$sampleID failed Dropsense QC.");
			next;
		}
		
		foreach my $verrichting(@verrichtingList) {
			if ($verrichting eq 'CYP2D6') {
							
				$sampleID	= $fields[$header{'Sample'}];
				# get new CYP2D6 string;
				$gene		= &_convertCYP2D6( $sampleID, $fields[$header{$verrichting}], &_parseTaqman($taqman,$sampleID));
				&_printGLIMP($sampleID, $verrichting, $gene, $output);
			
			}elsif ($verrichting eq 'DPYD') {
				
				$sampleID	= $fields[$header{'Sample'}];
				#get correct DPYD;
				$gene		= &_convertDPYD($fields[$header{$verrichting}],$fields[$header{'DPYD (rs56038477)'}],$fields[$header{'DPYD (rs67376798)'}]);
				
				#print "$sampleID, $verrichting, $gene, $output \n" ;
				&_printGLIMP($sampleID, $verrichting, $gene, $output);
							
			}elsif ($verrichting eq 'UGT1A1') {
				
				$sampleID	= $fields[$header{'Sample'}];
				
				#get correct UGT1A1;
				$gene = &_convertUGT1A1($UGT1A1,$sampleID);
				
				#print "$sampleID, $verrichting, $gene, $output \n" ;
				&_printGLIMP($sampleID, $verrichting, $gene, $output);
							
			}
			else {
			
				$sampleID	= $fields[$header{'Sample'}];
				$gene		= $fields[$header{$verrichting}];
				#print "$sampleID, $verrichting, $gene, $output \n" ;
				
				if ($gene eq 'unexpected') {
					$gene = $iPLEXUnexpected;
				}
				&_printGLIMP($sampleID, $verrichting, $gene, $output);
			}
		}
	}
}
} 

# TODO: Wat als er homozygote snp wordt gecalled?

#
# Gets 3 DPYD values, and returns converted results back.
sub _convertDPYD{
my $DPYD1 = shift;
my $DPYD2 = shift;
my $DPYD3 = shift;
my $newDPYD = '';
my $REF = '';
my $ALT = '';

my $log_message = "DPYD1 = $DPYD1 , DPYD2 = $DPYD2 , DPYD3 = $DPYD3.";
$logger->info($log_message);

my ($GEN1, $GEN2) = split /\s*\/\s*/, $DPYD1;

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


if ($taqman eq $iPLEXUnexpected){
	my $log_message = "$sampleID not in Taqman file: $taqmanCN." ;
	$logger->warn($log_message);
	return $taqmanCN;
} else{
	
	my ($var1,$var2) = split /\//, $CYP2D6;

	print $var1 . ' en ' . $var2."\n\n";
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
				
				#print "gelijk $voorspeld eq $taqmanCN \n\n";
				$newCYP2D6 .= $allelOne.'/'.$allelTwo . ' OR ';
			}
			else{
				my $log_message	= "sample: $sampleID predicted allels: $allelOne/$allelTwo: ";
				$log_message	.= ($Exon_9{$allelOne} + $Exon_9{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo}).','.($Intron_6_2{$allelOne} + $Intron_6_2{$allelTwo})." Taqman called $taqmanCN";
				$logger->warn($log_message);
				
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

while (<$fh>) {
	if( my ( $id, $var1, $var2 ) = $_ =~ m/^([0-9]*)\t([0-9]*)\t([0-9]*)$/ ) {
				#print "TEST: $id, $var1, $var2";
		if ($sampleID eq $id){
			#check for non existing values.
			if (! exists $UGT1A1Table{$var2} || ! exists $UGT1A1Table{$var1}){
				$UGT1A1line = $iPLEXUnexpected;
				$logger->warn("$sampleID: $var1 or $var2 does not exist");
				#next;
			}
			if ($UGT1A1Table{$var2} eq "*1"){
				$UGT1A1line = $UGT1A1Table{$var2}.'/'.$UGT1A1Table{$var1};
				$logger->info("$id, $var1, $var2");
			}else{
				print "\n\n bestaat? $sampleID\n\n";
				$UGT1A1line = $UGT1A1Table{$var1}.'/'.$UGT1A1Table{$var2};	
				$logger->info("$id, $var1, $var2");
			}
		}
	}
}	
	close($fh);
	return $UGT1A1line;
}

sub _parseUGT1A1{
my $UGT1A1 = shift;
my $output = shift;
unlink($output);

my $fh;
open($fh, '<', $UGT1A1) or die "Can't open UGT1A1 file.";

while (<$fh>) {
	
	if( my ( $id, $var1, $var2 ) = $_ =~ m/^([0-9]+)\t([0-9]+)\t([0-9]+).+$/ ) {
			
			&_printGLIMP($id,'UGT1A1', &_convertUGT1A1($UGT1A1,$id), $output); 
	}
}	
	close($fh);
}

sub _parseTaqman{
my $taqmanFile = shift;
my $sampleID = shift;
my %Hs00010001_cn = (); # Hs00010001_cn (=exon 9)
my %Hs04502391_cn = (); # Hs04502391_cn (=intron 6)
my %Hs04083572_cn = (); # Hs04083572_cn (=intron 2)
my $taqmanFileCSV = basename($taqmanFile,  ".xls") . ".csv";
print "filenaam: $taqmanFileCSV\n";

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

# parse the rest
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^Blanco/){
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
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^Blanco/){
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
while (my $row = $csv->getline_hr($tagman_fh)) {
	if ($row->{'Sample Name'} =~ m/^Blanco/){
		last;
	} else{
		$Hs04502391_cn{$row->{'Sample Name'}} = $row->{'CN Predicted'};
    }
}

$csv->eof or $csv->error_diag;
close $tagman_fh;

if (exists $Hs00010001_cn{$sampleID}){
	return ("$Hs00010001_cn{$sampleID},$Hs04502391_cn{$sampleID},$Hs04083572_cn{$sampleID}");
} else{
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
my $filenameOut = $filename . '.csv';
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
 
 
sub _Usage {
	
print <<EOF;
#########################################################################################################
# This script creates a GLIMPs output file using 4 inputfile (Dropsense | Taqman | Iplex | UGT1A1 ).    #
#                                                                                                       #
#########################################################################################################
Usage: perl FGX.pl [options]

Options:
-d dropsensefile (.xls)
-o outputfile 	 (';' seperated .csv/txt)
-t taqman		 (.xls)
-i iplex file	 (.csv)
-u UGT1A1 file	 (.txt)
-l log_level

#########################################################################################################
EOF

}
 
 
 
 
 
 
 
 
 
 
 
 
 
