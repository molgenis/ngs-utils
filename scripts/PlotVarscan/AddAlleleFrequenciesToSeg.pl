#!/usr/bin/perl
BEGIN{
	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
}
use strict;
use warnings;
use Vcf;
use List::Util qw(sum min max);
use List::MoreUtils qw(uniq);
use Data::Dumper;

#errors: somewhere an empty record is defined, counsing an undef printout;

main();

sub main {
	#my $vcfFile='/home/terpstramm/Documents/projects/privsmeta/TestData/PriVsMeta_S12_204/vcf/PriVsMeta_S12_204.ug.vcf';
	#my $segFile='/home/terpstramm/Documents/projects/privsmeta/TestData/PriVsMeta_S12_204/cnv/S12_204_5.seg';
	#my $normalSample="S12_204_13";
	my $vcfFile=$ARGV[0];
	my $segFile=$ARGV[1];
	my $normalSample=$ARGV[2];
	my $hetFreqCutoff=0.45;
	my $vcfDPbySampleCutoff=20;
	
	my $vcf = Vcf->new(file=>$vcfFile) or die "cannot open vcf file $vcfFile\n";
	$vcf->parse_header();
	
	my $segdata = SegLoad($segFile);
	#die Dumper($segdata);
	$segdata = SegAnnotateWithVcfSampleSpecific($segdata, $vcf, $hetFreqCutoff, $vcfDPbySampleCutoff,$normalSample);
	
	print SegAsSegString($segdata);
	
}
sub SegAnnotateWithVcf{
	#($segdata, $vcf);
	my $seg = shift @_;
	my $vcf = shift @_;
	my $hetFreqCutoff = shift @_;
	my $vcfDPbySampleCutoff = shift @_;
	my $norm = shift @_;
	my $newSeg = $seg;
	#@{$$newSeg{"HEADER"}}=@{$$seg{"HEADER"}};
	
	while (my $x=$vcf->next_data_hash()){
		my @samples = keys(%{$x->{gtypes}});
		#warn Dumper(@samples); 
		my $f;
		for my $sample (@samples){
			next if(not($x->{'gtypes'}{$sample}{'DP'}) || $x->{'gtypes'}{$sample}{'DP'} < $vcfDPbySampleCutoff);
			next if(defined($x->{'gtypes'}{$sample}{'AD'}) && $x->{'gtypes'}{$sample}{'AD'} =~ m/\./);
			my @AD = DataHashGetADBySample($x, $sample);
			next if(not(sum(@AD)) || sum(@AD) < $vcfDPbySampleCutoff);
			
			my $g = DataHashGetGBySample($x, $sample);
			my $gw=$g*sum(@AD);
			#die "g=$g; f=$f; sumAD=".sum(@AD);
			#error if g=
			#my $zg=abs($f-0.5)*sqrt(sum(@AD));
			my $gnorm =DataHashGetGBySample($x, $norm);			
			next if($gnorm > $hetFreqCutoff);
			
			for (my $i = 0; $i < scalar((@{$$seg{"DATA"}})); $i++) {
				if(SegIntersectsVcfDatahash(${$$seg{"DATA"}}[$i],$x)){
					push(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}},$gw);
					
					${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"}+=sum(@AD);
					#push(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"ZG"}},$zg);
					#push @{@{$$newSeg{"DATA"}}});
					#die Dumper($x,$newSeg);
				}
			}
			#die Dumper($x,$seg);

		}	
	}
	
	for my $sample (keys(%{${${$$seg{"DATA"}}[0]}{'genotypes'}})){
		
			push(@{$$newSeg{"HEADER"}}, $sample.'.g')if(scalar(uniq(@{$$newSeg{"HEADER"}})) != scalar(uniq((@{$$newSeg{"HEADER"}},$sample.'.g'))));
			push(@{$$newSeg{"HEADER"}}, $sample.'.dp')if(scalar(uniq(@{$$newSeg{"HEADER"}})) != scalar(uniq((@{$$newSeg{"HEADER"}},$sample.'.dp'))));
		
		for (my $i = 0; $i < scalar((@{$$seg{"DATA"}})); $i++) {			
#			warn "DEBUG:What value is this??:'".scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})."'";
			#print Dumper(%{${$$seg{"DATA"}}[$i]});
			if(exists(${${$$seg{"DATA"}}[$i]}{"genotypes"}) && scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}}) > 0 ){
				#${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}=sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})/scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}});
				#warn "dpsum=".${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"}."Gwsum=".sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}});
				${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}=sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})/${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"};
				${${$$newSeg{"DATA"}}[$i]}{$sample.'.dp'}=${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"};
				
			}else{
				${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}='NA';
				${${$$newSeg{"DATA"}}[$i]}{$sample.'.dp'}='NA';
			}
			#push @{@{$$newSeg{"DATA"}}});
			#die Dumper($newSeg);
		}
	}
	return $newSeg;
}
sub SegAnnotateWithVcfSampleSpecific{
	#($segdata, $vcf);
	my $seg = shift @_;
	my $vcf = shift @_;
	my $hetFreqCutoff = shift @_;
	my $vcfDPbySampleCutoff = shift @_;
	my $norm = shift @_;
	my $newSeg = $seg;
	#@{$$newSeg{"HEADER"}}=@{$$seg{"HEADER"}};
	
	while (my $x=$vcf->next_data_hash()){
		my @samples = (${${$$seg{"DATA"}}[0]}{"ID"});
		#warn Dumper(@samples); 
		my $f;
		for my $sample (@samples){
			next if(not($x->{'gtypes'}{$sample}{'DP'}) || $x->{'gtypes'}{$sample}{'DP'} < $vcfDPbySampleCutoff);
			next if(defined($x->{'gtypes'}{$sample}{'AD'}) && $x->{'gtypes'}{$sample}{'AD'} =~ m/\./);
			my @AD = DataHashGetADBySample($x, $sample);
			next if(not(sum(@AD)) || sum(@AD) < $vcfDPbySampleCutoff);
			
			my $g = DataHashGetGBySample($x, $sample);
			my $gw=$g*sum(@AD);
			#die "g=$g; f=$f; sumAD=".sum(@AD);
			#my $zg=abs($f-0.5)*sqrt(sum(@AD));
			my $gnorm =DataHashGetGBySample($x, $norm);			
			next if($gnorm > $hetFreqCutoff);
			
			for (my $i = 0; $i < scalar((@{$$seg{"DATA"}})); $i++) {
				if(SegIntersectsVcfDatahash(${$$seg{"DATA"}}[$i],$x)){
					push(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}},$gw);
					
					${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"}+=sum(@AD);
					#push(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"ZG"}},$zg);
					#push @{@{$$newSeg{"DATA"}}});
					#die Dumper($x,$newSeg);
				}
			}
			#die Dumper($x,$seg);

		}	
	}
	
	my $sample = ${${$$seg{"DATA"}}[0]}{"ID"};
	push(@{$$newSeg{"HEADER"}}, $sample.'.g')if(scalar(uniq(@{$$newSeg{"HEADER"}})) != scalar(uniq((@{$$newSeg{"HEADER"}},$sample.'.g'))));
	push(@{$$newSeg{"HEADER"}}, $sample.'.dp')if(scalar(uniq(@{$$newSeg{"HEADER"}})) != scalar(uniq((@{$$newSeg{"HEADER"}},$sample.'.dp'))));
		
	for (my $i = 0; $i < scalar((@{$$seg{"DATA"}})); $i++) {
		
		#print Dumper((%{${${$$seg{"DATA"}}[$i]}{"genotypes"}}));
		#warn "DEBUG:What value is this??:'".scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})."'";
		if(exists(${${$$seg{"DATA"}}[$i]}{"genotypes"}) && scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}}) > 0 ){
			#${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}=sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})/scalar(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}});
			#warn "dpsum=".${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"}."Gwsum=".sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}});
			${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}=sum(@{${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"G"}})/${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"};
			${${$$newSeg{"DATA"}}[$i]}{$sample.'.dp'}=${${$$seg{"DATA"}}[$i]}{"genotypes"}{$sample}{"dpSum"};
			
		}else{
			${${$$newSeg{"DATA"}}[$i]}{$sample.'.g'}='NA';
			${${$$newSeg{"DATA"}}[$i]}{$sample.'.dp'}='NA';
		}
		#push @{@{$$newSeg{"DATA"}}});
		#die Dumper($newSeg);
		
	}
	return $newSeg;
}

sub DataHashGetADBySample{
	my $x = shift @_;
	my $sample = shift @_;
	die "sample '$sample' not found" if(not(defined($x->{'gtypes'}{$sample})));
	if(defined($x->{'gtypes'}{$sample}{'AD'}) && not($x->{'gtypes'}{$sample}{'AD'} =~ m/^\.,\./)){
		my $AD;
		@{$AD}=split(',',$x->{'gtypes'}{$sample}{'AD'});
		return @{$AD};
	};
}
sub DataHashGetGBySample{
	my $x = shift @_;
	my $sample = shift @_;
	die "sample '$sample' not found" if(not(defined($x->{'gtypes'}{$sample})));
	
	if(defined($x->{'gtypes'}{$sample}{'AD'}) && not($x->{'gtypes'}{$sample}{'AD'} =~ m/^\.,\./)){
		my $AD;
		@{$AD}=DataHashGetADBySample($x,$sample);
		my $F;
		my $i=0;
		while($i < scalar(@{$AD})){
			if( defined(${$AD}[0]) && defined(${$AD}[1]) && sum(@{$AD}) > 0 && ${$AD}[$i] > 0 ){
				${$F}[$i] = (${$AD}[$i] / sum(@{$AD}));
			}else{
				${$F}[$i] = '0';
			}
			$i++;
		}
		my $f=max(@{$F});
		
		my $g=abs($f-0.5);
		return $g;
	};
	
}
sub SegLoad {
	my $segFile = shift(@_);
	
	my $newSeg;
	
	open(my $segHandle,'<',$segFile) or die "Cannot open .seg file: '$segFile'";
	
	my $h=<$segHandle>;
	chomp $h;
	$h =~ s/"|'//g;
	my @h= split("\t",$h);
	@{$$newSeg{"HEADER"}} =  @h;
	
	while(<$segHandle>){
		chomp;
		s/"|'//g;
		my @segline = split("\t");
		my $i=0;
		my %record;
		map{s/"|'//g;$record{${$$newSeg{"HEADER"}}[$i]}=$_;$i++;}(@segline);
		push(@{$$newSeg{"DATA"}}, \%record);
	}
	return $newSeg;
}
sub SegAsSegString {
	my $seg=shift(@_);
	my $string=join("\t",@{$$seg{"HEADER"}})."\n";
	for my $record (@{$$seg{"DATA"}}){
		my @r;
		
		for my $header (@{$$seg{"HEADER"}}){
			if(defined($$record{$header})){
				push(@r,$$record{$header});
			}else{
				push(@r,"NA");
			}
			
		}
		
		$string=$string.join("\t",@r)."\n";
		#die Dumper( \$header, \@{$$seg{"HEADER"}}. \@r) if(not(defined($r[10])));
	}
	return $string;
}
sub SegAsSegStringDefault {
	my $seg=shift(@_);
	my $string=join("\t",@{$$seg{"HEADER"}})."\n";
	for my $record (@{$$seg{"DATA"}}){
		my @r;
		
		for my $header (("ID","chrom","loc.start","loc.end","num.mark","seg.mean")){
			if(defined($$record{$header})){
				push(@r,$$record{$header});
			}else{
				push(@r,"NA");
			}
			
		}
		
		$string=$string.join("\t",@r)."\n";
		#die Dumper( \$header, \@{$$seg{"HEADER"}}. \@r) if(not(defined($r[10])));
	}
	return $string;
}
sub SegAsBedString {
	my $seg=shift(@_);
	my $bed;
	#=join("\t",@{$$seg{"HEADER"}})."\n";
	for my $record (@{$$seg{"DATA"}}){
		my @r;
		
		for my $header (("chrom","loc.start","loc.end","ID")){
			if(defined($$record{$header})){
				push(@r,$$record{$header});
			}else{
				push(@r,"NA");
			}
			
		}
		
		$bed=$bed.join("\t",@r)."\n";
		#die Dumper( \$header, \@{$$seg{"HEADER"}}. \@r) if(not(defined($r[10])));
	}
	return $bed;
}
sub SegIntersectsVcfDatahash {
	my $recordSeg = shift @_;
	my $vcfHash = shift @_;
	#			${$x}{'CHROM'}
	#		${$x}{'POS'}
	if($$recordSeg{"chrom"} && $$recordSeg{"loc.start"} && $$recordSeg{"loc.end"} && ${$vcfHash}{'CHROM'} && ${$vcfHash}{'POS'}){
		if( $$recordSeg{"chrom"} eq ${$vcfHash}{'CHROM'} && $$recordSeg{"loc.start"} < ${$vcfHash}{'POS'} && $$recordSeg{"loc.end"} > ${$vcfHash}{'POS'}){
			return 1;
		}else{
			return 0;
		}
	}else{
		die "invalid seg/vcf record???"
	}
}
