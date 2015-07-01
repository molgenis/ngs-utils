#!/usr/bin/perl
BEGIN{
	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
}
use warnings;
use strict;
use Vcf;
use Data::Dumper;
# input test
#$ARGV[0]='/home/terpstramm/Downloads/s12-193.annotated.vcf';
# use
my $use = <<"END";
	$0 in.vcf out.vcf
	ads Z scores based on the AD fields per sample
END
die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
AddZScoresBasedOnAdVals($ARGV[0],$ARGV[1]);

sub AddZScoresBasedOnAdVals{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn 'warning: no output vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	
	
	my $vcf = Vcf->new(file=>$vcfin);
	$vcf->parse_header();
	$vcf->add_header_line({key=>'FORMAT', ID=>'F',Number=>-1,Type=>'Float',Description=>'F as estimation of allele frequency. Derived from AD as F=$AD[1]/sum(@AD). if defaults to 0 when $AD[1] not present or when sum $AD[0] + $AD[1] <= 0 '});
	$vcf->add_header_line({key=>'FORMAT', ID=>'Z',Number=>-1,Type=>'Float',Description=>'Z score calcuated from AD vals: Z = 2 * ($AD[1] / ($AD[1] + $AD[0]) - 0.5) * sqrt($AD[1] + $AD[0]). Not calculated for multiallelic variants.'});
	$vcf->add_header_line({key=>'FORMAT', ID=>'Y',Number=>-1,Type=>'Float',Description=>'Y deived from Y=1-Z'});
	print {$out} $vcf->format_header();
	#$vcf->recalc_ac_an(0);
	while (my $x=$vcf->next_data_hash()){
		#warn $x->{'INFO'}{'VariantType'};
		if(not($x->{'INFO'}{'VariantType'} =~ m/^MULTIALLELIC/)){
			$vcf->add_format_field($x,'Z'); 
			$vcf->add_format_field($x,'Y'); 
			my @samples = keys(%{$x->{gtypes}});
			#warn Dumper(@samples); 
			for my $sample (@samples){
				my @AD = split(',',$x->{'gtypes'}{$sample}{'AD'})if(defined($x->{'gtypes'}{$sample}{'AD'})&& not($x->{'gtypes'}{$sample}{'AD'} =~ m/^\.,\./));
				my $Zscore;
				if( defined($AD[0]) && defined($AD[1]) && scalar(@AD)== 2 && $AD[0] + $AD[1] > 0 ){
					$Zscore = 0.5 * ($AD[1] / ($AD[1] + $AD[0]) - 0.5) * sqrt($AD[1] + $AD[0]);
				}else{
					$Zscore = 0;
				}
				$x->{'gtypes'}{$sample}{'Z'}=$Zscore;
				$x->{'gtypes'}{$sample}{'Y'}=(1-$Zscore);
				
				my $F;
				if( defined($AD[0]) && defined($AD[1]) && scalar(@AD)== 2 && $AD[0] + $AD[1] > 0 && $AD[1] > 0 ){
					
					$F = ($AD[1] / sum(@AD));
					
				}else{
					
					$F = '0';
				}
				
				$x->{'gtypes'}{$sample}{'F'}=$F;
			}
		}
		#}else{
		#	warn Dumper($x);die;
		#}
		
		#warn Dumper($vcf);
		#warn Dumper($x);
		
		print {$out} $vcf->format_line($x);
		#die;
		#die if ($. > 50);
	}
	$vcf->close();
}
