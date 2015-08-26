#!/usr/bin/env perl

# This scripts retrieves protein coding canonical exon sequences from EnsEMBL and creates a BED file representing canonical exon ranges.
# Coded by A.T Miedema (a [dot] t [dot] miedema [at] st [dot] hanze [dot] nl)
#
# usage "CanonicalSeqRetriever.pl [output file]"

my $file= "";
if (@a){
	$file = $ARGV[0];
}else{
	print("This script retrieves protein coding canonical exon sequences from EnsEMBL and creates a BED file representing canonical exon ranges.\n");	
	print("usage \"CanonicalSeqRetriever.pl [output file]\".\n");
	exit 0;
}

# make things easier

use strict;
use warnings;

# import modules
use Bio::EnsEMBL::Registry;
use Data::Dumper;

# setup registry
my $registry = 'Bio::EnsEMBL::Registry';

# connect to EnsEMBL
$registry->load_registry_from_db(
				-host => "ensembldb.ensembl.org",
                -user => "anonymous",
				-db_version  => 75
				);

# Handy function for testing purposes.
# Prepares a nice toS0tring like function
sub feature2string
{
	my $feature = shift;

	my $stable_id  = $feature->stable_id();
	my $seq_region = $feature->slice->seq_region_name();
	my $start      = $feature->start();
	my $end        = $feature->end();
	my $strand     = $feature->strand();

	return sprintf( "%s: %s:%d-%d (%+d)",
		$stable_id, $seq_region, $start, $end, $strand );
}

# get gene adaptor object from registry for human core
my $gene_adaptor = $registry->get_adaptor("Human", "Core", "Gene");


my %geneSymbolCounts = ();
my %geneSymbolToBed = ();

my $count = 0;
# gene symbols with more then one transcripts
my $geneSymbolMTOTranscripts = 0;
my $defined_count = 0;
my $undefined_count = 0;

# get list of gene stable IDs
my $gene_ids = $gene_adaptor->list_stable_ids();
#print "Processing " . scalar(@{$gene_ids}) . " gene IDs...\n";
# traverse gene IDs
while (my $gene_id = shift(@{$gene_ids})) {
	
    # let user know current counts
    local $| = 1;
    print "[$count/" . scalar(@{$gene_ids}) . "]\r";
	
    # get gene object
    my $gene = $gene_adaptor->fetch_by_stable_id($gene_id);
	
    # get canonical transcript
    my $canonical_transcript = $gene->canonical_transcript();
	my $geneSymbol = $gene->external_name;
	# check if canonical transcript is defined
    if (defined $canonical_transcript && defined $geneSymbol) {

		# check if transcript is protein coding
		if(defined $canonical_transcript->translate()){
			# print("gene name: ",$gene->display_xref()->display_id() , "\n");
			my $isCorrectExonLocations = 1;
			my $currChrom = "";
			
			# retrieve all exons
			foreach my $exon ( @{$canonical_transcript->get_all_Exons() } ) {
				
				my $estring = feature2string($exon);
				my @splittedExon = split(":", $estring);
				$splittedExon[1] =~ s/^\s+|\s+$//g;
						
				# Check if the chromosome is not located on a patch contig
				if ($splittedExon[1] =~ /^[\dXxYyMTmt]+$/){
					#print($estring, "\n");
					# BED starts are zero-based and BED ends are one-based.
					# Hence we need to substract 1 from the exon start to cover the complete region in the BED file.
					# The other 5 bases is used to increase the window to capture splice variants aswel 
					my $chrom = $exon->slice->seq_region_name();
					my $start = $exon->start()-6;
					my $end = $exon->end()+5;
					my $name  = "HGNC_Symbol=". $geneSymbol . "|ENSMBL_id=". $exon->stable_id();
					my $score = "0";
					my $strand = $exon->strand();
					my $line = join("\t", $chrom, $start, $end, $name, $score, $strand);
					
					if(!exists($geneSymbolToBed{$geneSymbol})){

						$geneSymbolToBed{$geneSymbol} = [$line] ;
						
					}else{
						
						push ($geneSymbolToBed{$geneSymbol}, $line);
					}
					
					# checks if all exons are on the same chromosome
					if ($currChrom eq ""){
						$currChrom = $splittedExon[1];
					}
					if(!($currChrom eq $splittedExon[1])){
						$isCorrectExonLocations = 0;
					}
				}
				
			}
			if ($isCorrectExonLocations == 1){
				$defined_count++;
				# get gene name by cross referencing and add it to correctly validated ids
				if(!exists($geneSymbolCounts{$geneSymbol})){
					$geneSymbolCounts{$geneSymbol} = 1 ;
				}else{
					$geneSymbolCounts{$geneSymbol}++; 
				}
			}
					
			undef $canonical_transcript;			
		}
		else{
			$geneSymbolMTOTranscripts++;
		}
		
	}
    else {
        $undefined_count++;
    }
    $count++;

    # undef the transcript
    $canonical_transcript = undef;
	$gene = undef;
	$geneSymbol = undef;
}




open (my $fh, ">>", $file)
  or die("Error opening $file: $!");

my $badGeneCounts = 0;
my $goodGeneCounts = 0;
my $totalGeneSymbols = 0;
foreach my $key (keys %geneSymbolCounts) {
	
	if($geneSymbolCounts{$key} > 1){
		$badGeneCounts++;
		$totalGeneSymbols++;
	}
	else{
		if(exists($geneSymbolToBed{$key})){
			foreach my $bedFileLine (@{$geneSymbolToBed{$key}}){
				print $fh $bedFileLine."\n";
			#print($geneSymbolToBed{$key}, "\n");
			
			}
			$goodGeneCounts++;
			$totalGeneSymbols++;
		}
	}
}
close $fh;


#my $bed_seq			= $ex_seq;
#my $bed_start		= $ex_start - 1;
#my $bed_end			= $ex_end;
#my $bed_annotation	= join('|', 'HGNC_Symbol=' . $hgnc_symbol, 'REFSEQ_ID=' . $hgnc_id, $ex_id, $tr_id, $gene_id);


# Show the user some counts of filter steps
print("Total number of genes without protein coding conanical sequence: ", $geneSymbolMTOTranscripts, "\n");					
print("Total counts of gene symbols with protein coding canonical sequence: ", $totalGeneSymbols, "\n");
print("Number of gene symbols which have more than one canonical sequence in ensembl: ", $badGeneCounts, "\n");
print("Number of passed gene names: ", $goodGeneCounts, "\n");
# let the user know
print "$defined_count remained \& $undefined_count filtered in $count.\n";
print "...done!\n";
print "Results saved to: ", $file;