#!/usr/bin/perl
use strict;
use warnings;
my $use = <<"END";
	calculates filtering flags for vcf files:
		per Patient specified by Gerard te Meerman and stored in the general info fields:
			4*sumAD-((maxAD/sumAD)-0.5)**2
			IsHetNormal
			IsNovelInTumor (not normal)
			MeanDepth
	use perl $0 NomalSampleID VCFFile.vcf
		NomalSampleID =  normal tissue other samples are assumed to be tumorSamples
		VCFFile.vcf = unannotated VCF file	
END

my $Normal = $ARGV[0];
my $VcfFile = $ARGV[1];

open(my $in,'<',$VcfFile)or die "cannot open $VcfFile\n$!";
my $line;
	
#general stuff
my %headToIndex;
my @samples;
#Parse header / define %validChroms
my $header = 1;

while($header == 1){
	$line = <$in>;
	$line =~ s/\r$|\n$|\r\n$//g;
	if($line =~ /##/){
		#SafetyParsing
		#warn "WARN: headerline has invalid contig/formatfield: skipped header line:'$line'\n" if(not($line =~ /^##contig=<ID=\d,length=|^##contig=<ID=|^##fileformat|^##fileDate|^##reference|^##INFO=<ID=[A-EG-Z]{2}|^##FILTER|^##FORMAT/));
		#next if(not($line =~ /^##contig=<ID=\d,length=|^##contig=<ID=|^##fileformat|^##fileDate|^##reference|^##INFO=<ID=[A-EG-Z]{2}|^##FILTER|^##FORMAT/));
		print $line."\n";
	}else{
		$header = 0;
	}
}
#print "%%%%%%%%%".$line."\n";
#Print Additional info to Header;
my $customINFO = <<"END";
##TumorGenotypeScriptUse="perl $0 Normal in.vcf > out.vcf"
##TumorGenotypeScriptCmd="perl $0 $Normal $VcfFile"
##INFO=<ID=IsHetrozygousInNormal,Number=0,Type=Flag,Description="Has two different alleles in normal sample (Is hetrogenous)">
##INFO=<ID=IsChangeInTumor,Number=0,Type=Flag,Description="Is changed in a tumor sample relative to Normal sample or other Tumor samples">
##INFO=<ID=IsChangeInTumorGQnormal,Number=1,Type=Float,Description="IsChangeInTumor normal GQ max for deciding IsChangeInTumor flag might not always be present.">
##INFO=<ID=IsChangeInTumorGQtumor,Number=1,Type=Float,Description="IsChangeInTumor tumor GQ max for deciding IsChangeInTumor flag">
##INFO=<ID=IsChangeInTumorGQtumorcritical,Number=0,Type=Flag,Description="Is changed in a tumor sample relative to Normal sample or other Tumor samples">
##INFO=<ID=Somatic,Number=0,Type=Flag,Description="Is changed in a tumor sample relative to Normal sample, the normal does not contain the tumor variant at an allele frequency of 0.05 and tumor gt is not ref (somewhat conservative)">
##INFO=<ID=TeMeermanAlleleBias,Number=1,Type=Float,Description="Allele bias as Specified by Gerard te Meerman '4*sumAD-((maxAD/sumAD)-0.5)**2' (only calculated for Het sites for the specified normal sample => $Normal )">
##INFO=<ID=MeanDP,Number=1,Type=Float,Description="Mean of sample depths (meanDP)">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum of sample depths (meanDP)">
##FORMAT=<ID=F,Number=1,Type=Float,Description="major allele frequency per sample for hetrozygous sites calculated by F=(AD[1]/sumAD) because depth in the VCF file is not always sumAD">
##FORMAT=<ID=Z,Number=1,Type=Float,Description="Z score calculated from AD vals: Z = (\$AD[1] / (\$AD[1] + \$AD[0]) - 0.5) * sqrt(\$AD[1] + \$AD[0]). Not calculated for multiallelic variants.'">
END
print $customINFO;

#
##Parse VCF header / define %headToIndex
#

my $head = substr($line,1);

#print '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'."\n";
my @head = split("\t",$head);
my $i = 0;
map{my $tmp = $_; $headToIndex{$tmp}=$i;push(@samples,$tmp)if($i>8);$i++;}(@head);
#print $headToIndex{"FORMAT"}."\n";

#addaditional tumor normal index
my $normalIndex = SetNormalIndex(\$Normal,\@samples);
my @tumorSamples = @samples;
splice(@tumorSamples,$normalIndex,1);#warn $_
#warn '"'.join('#',@tumorSamples).'"'.join('#',@samples).'"'."\n";
#print header
print '#'.join("\t",(@head))."\n";

#
##do variant parsing
#

my @tabSplit;	
while($line=<$in>){

	$line =~ s/\R$//g;
	@tabSplit = split("\t",$line);	
	
	#print $headToIndex{"CHROM"}."\n";
	#die "exit on line :'$.'\n$!" if ($. == 10000);#debugging 

	my %info = infoReader($tabSplit[$headToIndex{"INFO"}]);
	#print "'".join("#",@tabSplit[@headToIndex{@samples}])."'\n";
	my $a = $tabSplit[$headToIndex{"FORMAT"}];
	my @b = @tabSplit[@headToIndex{@samples}];
	my %formatInfo = formatReader(\$tabSplit[$headToIndex{"FORMAT"}],\@b,\@samples);
	#
	###
	###########youre able to recalculate/ delete tags freely here
	###
	#
	#delete $info{"AF"};
	#calulcate IsHetrozygousInNormal and IsChangeInTumor flags: allele based check for gain of allele or loss of allele normal to tumor.
#	if($formatInfo{$Normal}{'GT'} ne './.'){#assume detected known normal genotype and skip tumor sample if not detected genotype
#		my @normalGenotypeIndexes = split('/',$formatInfo{$Normal}{'GT'});
#		die "Invalid triallelic site \n'$line'\n$!" if(scalar(@normalGenotypeIndexes)>2);#never happens just to check if the BUG doesn't spit out triallelic genotypes like 0/0/1
#		if($normalGenotypeIndexes[0] ne $normalGenotypeIndexes[1]){
#			$info{"IsHetrozygousInNormal"} = 'true';
#		}
#		##isNovelInTumor
#		my %tumorGTs=();
#		for my $sampleField (@tumorSamples){
#			if($formatInfo{$sampleField}{'GT'} ne './.'){
#				map{$tumorGTs{$sampleField} = 1;}(split('/',$formatInfo{$sampleField}{'GT'}));
#			}
#		}
#		my %normalGTs;
#		map {$normalGTs{$_}=1;} (@normalGenotypeIndexes);
#	
#		map {$info{"IsChangeInTumor"} = 'true' if(not($normalGTs{$_}))}(keys(%tumorGTs));
#		map {$info{"IsChangeInTumor"}  = 'true' if(not( $tumorGTs{$_}))}(keys(%normalGTs))
#
#	}else{#assume undetected unknown normal genotype (changes in tumorsamples then the so only the IsNovelInTumor is assinged)
#		my $oldGt = '';
#		for my $tumorGt ($formatInfo{@tumorSamples}{'GT'}){
#			$info{"IsChangeInTumor"} = 'true' if($oldGt ne '' && $tumorGt ne './.' && $oldGt ne $tumorGt);
#			$oldGt = $tumorGt if($tumorGt ne './.');
#		}
#	}
###########################################################
####################################################################################################################
	
	if($formatInfo{$Normal}{'GT'} ne './.'){
		####################################################################################################################
		##########################Test for is het in normal
		my @normalGenotypeIndexes = split('/',$formatInfo{$Normal}{'GT'});
		#my $sumnormalAD=sum(split)
		if($normalGenotypeIndexes[0] ne $normalGenotypeIndexes[1]){
			$info{"IsHetrozygousInNormal"} = 'true';
		}
		####################################################################################################################
		##########################test for Somatic and test for Alteration of tumor genotypes in respect to the normal genotype
		die "Invalid triallelic in normal sample site \n'$line'\n$!" if(scalar(@normalGenotypeIndexes)>2);#never happens just to check if the BUG doesn't spit out triallelic genotypes like 0/0/1
		#$formatInfo{$Normal}{'GT'}
		for my $sampleField (@tumorSamples){
			my @tumorGenotypeIndexes=split('/',$formatInfo{$sampleField}{'GT'});
			die "Invalid triallelic in normal sample site \n'$line'\n$!" if(scalar(@tumorGenotypeIndexes)>2);#never happens just to check if the BUG doesn't spit out triallelic genotypes like 0/0/1
			
			if($formatInfo{$sampleField}{'GT'} ne './.' && join('/',sort(split('/',$formatInfo{$Normal}{'GT'}))) ne join('/',sort(split('/',$formatInfo{$sampleField}{'GT'}))) ){
				my @indexesTumor;
				@indexesTumor=getDifferentIndexesB($formatInfo{$Normal}{'GT'},$formatInfo{$sampleField}{'GT'});
				#warn "return value of formatInfo.Normal.'AD.indexesTumor':" . join(",",$formatInfo{$Normal}{'ADi'}{@indexesTumor})."\n";
				$info{"Somatic"} = 'true' if($normalGenotypeIndexes[0] eq $normalGenotypeIndexes[1] && $formatInfo{$sampleField}{'GT'} ne '0/0' &&(max($formatInfo{$Normal}{'ADi'}{@indexesTumor})==0|| max($formatInfo{$Normal}{'ADi'}{@indexesTumor})/Sum(split(',',$formatInfo{$Normal}{'AD'})) < 0.05  ));
				$info{"IsChangeInTumor"} = 'true';
				$info{"IsChangeInTumorGQnormal"} = $formatInfo{$Normal}{'GQ'};
				$info{"IsChangeInTumorGQtumor"} = $formatInfo{$sampleField}{'GQ'} if(not(defined($info{"IsChangeInTumorGQtumor"}))||$info{"IsChangeInTumorGQtumor"} < $formatInfo{$sampleField}{'GQ'});
			}
		}		
	}else{
		####################################################################################################################
		##########################test for Alteration of tumor genotypes between each other
		for my $sample1 (@tumorSamples){
			if($formatInfo{$sample1}{'GT'} ne './.'){
				for my $sample2 (@tumorSamples){
					if($formatInfo{$sample2}{'GT'} ne './.'){
						if($formatInfo{$sample1}{'GT'} ne $formatInfo{$sample2}{'GT'} ){
							$info{"IsChangeInTumor"} = 'true';
							$info{"IsChangeInTumorGQtumor"} = $formatInfo{$sample1}{'GQ'} if(not(defined($info{"IsChangeInTumorGQtumor"}))||$info{"IsChangeInTumorGQtumor"} < $formatInfo{$sample1}{'GQ'});
							$info{"IsChangeInTumorGQtumor"} = $formatInfo{$sample2}{'GQ'} if(not(defined($info{"IsChangeInTumorGQtumor"}))||$info{"IsChangeInTumorGQtumor"} < $formatInfo{$sample2}{'GQ'});
						}
					}
				}
			}
		}	
	}
	####################################################################################################################
	#zscore calc
	$tabSplit[$headToIndex{"FORMAT"}] = $tabSplit[$headToIndex{"FORMAT"}].':Z';
	for my $sample (@samples){
		my @AD = split(',',$formatInfo{$sample}{'AD'})if(defined($formatInfo{$sample}{'AD'})&& not($formatInfo{$sample}{'GT'} =~ m/^\.,\./));
		my $Zscore;
		if( defined($AD[0]) && defined($AD[1]) && scalar(@AD)== 2 && $AD[0] + $AD[1] > 0 ){
			$Zscore = ($AD[1] / ($AD[1] + $AD[0]) - 0.5) * sqrt($AD[1] + $AD[0]);
			#$Zscore = 0.5 * ($AD[1] / ($AD[1] + $AD[0]) - 0.5) * sqrt($AD[1] + $AD[0]);
		}else{
			$Zscore = 0;
		}
		$formatInfo{$sample}{'Z'}=$Zscore;
	}
	
	
	##ishetNormal
	#my @normalGenotypeIndexes = split('/',$formatInfo{$Normal}{'GT'});
	#die "Invalid triallelic site \n'$line'\n$!" if(scalar(@normalGenotypeIndexes)>2);
	#if($normalGenotypeIndexes[0] ne $normalGenotypeIndexes[1]){
	#	$info{"IsHetrozygousInNormal"} = 'true';
	#}
	##isNovelInTumor
	#my %tumorGTs=();
	#for my $sampleField (@tumorSamples){
	#	if($formatInfo{$sampleField}{'GT'} ne './.'){
	#		map{$tumorGTs{$sampleField} = 1;}(split('/',$formatInfo{$sampleField}{'GT'}));
	#	}
	#}
	#my $oldGT;
	#for my $GT (@$formatInfo{@tumorSamples}{'GT'}){
	#	
	#}
	#my %normalGTs;
	#if($formatInfo{$Normal}{'GT'} ne './.'){
	#	map {$normalGTs{$_}=1;} (@normalGenotypeIndexes);
	#}
	
	#map {$info{"IsNovelInTumor"} = 'true' if(not($normalGTs{$_}))}(keys(%tumorGTs));
	#map {$info{"IsLostInTumor"}  = 'true' if(not( $tumorGTs{$_}))}(keys(%normalGTs));
	##
	
	
	#major allele frequency = F as format field
	$tabSplit[$headToIndex{"FORMAT"}] = $tabSplit[$headToIndex{"FORMAT"}].':F';
	for my $sample (@samples){
		my @sampleAD =  split(',',$formatInfo{$sample}{'AD'}) if $formatInfo{$sample}{'GT'} ne './.';
		
		my @Gt = split('/',$formatInfo{$sample}{'GT'});
		if(($formatInfo{$sample}{'GT'} ne './.') && (Sum(@sampleAD) > 0 ) ){
			#die "Sanity error check the regex in the above line" if( substr($tmpGt,0,1) eq substr($tmpGt,1,1));
			
			$formatInfo{$sample}{'F'} = ($sampleAD[1]/Sum(@sampleAD));
			

		}elsif($formatInfo{$sample}{'GT'} ne './.'){
			$formatInfo{$sample}{'F'} = '.';
		}
	}
	#allele depth counts in array
	my @AD = split(',',$formatInfo{$Normal}{'AD'}) if $formatInfo{$Normal}{'GT'} ne './.';
	#TeMeermanAlleleBias
	if($formatInfo{$Normal}{'GT'} ne './.' && $info{"IsHetNormal"} && $info{"IsHetNormal"} eq 'true' && Sum(@AD) > 0){
		my @AD = split(',',$formatInfo{$Normal}{'AD'});
		#Sum(@AD), Mean(@AD), Min($formatInfo{@samples}{'DP'});
		$info{"TeMeermanAlleleBias"} = 4*Sum(@AD)-((max(@AD)/Sum(@AD))-0.5)**2;
	}
	
	#minDepth/meanDP
	my @DPs;
	map{push(@DPs, $formatInfo{$_}{'DP'}) if(defined($formatInfo{$_}{'DP'}) && $formatInfo{$_}{'DP'} ne "")}(@samples);
	#warn @DPs;
	$info{"MeanDP"} = Mean(@DPs);
	
	$info{"MinDP"}  = min(@DPs);
	
	#$formatInfo{"AD"}=join(",",calcADvals($formatInfo{"DP"},$formatInfo{"EC"}));
	#delete $formatInfo{"EC"};
	print join("\t",($tabSplit[$headToIndex{"CHROM"}],$tabSplit[$headToIndex{"POS"}],$tabSplit[$headToIndex{"ID"}],$tabSplit[$headToIndex{"REF"}],$tabSplit[$headToIndex{"ALT"}],$tabSplit[$headToIndex{"QUAL"}],$tabSplit[$headToIndex{"FILTER"}],&infoPrinter(%info),&formatPrinter(\$tabSplit[$headToIndex{"FORMAT"}],\%formatInfo)))."\n";
}

sub SetNormalIndex {
	my $index = 0;
	for(@{$_[1]}){
		if($_ eq ${$_[0]}){
			return $index;
		}
		$index++;
	}
	die "Name of sample not matched(did you make a typoo?)\n$!";
}

sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}
sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min;
}
sub Sum {
	my $sum =0;
	for (@_){
		$sum += $_;
	}
	return $sum;
}
sub Mean {
	my $count = scalar(@_);
	my $sum = 0;
	for (@_){
		$sum += $_ if($_);
	}
	my $mean=$sum/$count;
	return $mean;
}
sub infoReader {
	my $info = shift(@_);
	my @info = split(';', $info); 
	my %info;
	map{my $tmp = $_;($info{(split('=',$tmp))[0]}=(split('=',$tmp))[1])if($tmp =~ /=/);$info{$tmp} = 'true' if(not($tmp =~ /=/))}(@info);
	if($info{'AD'}){
		my $tmp = 0;
		map{$info{'ADi'}{$tmp} = $_;$tmp++}(split(',',$info{'AD'}));
		
	}
	return %info;
}
sub infoPrinter {
	my %info = @_ or warn "Internal error has occured at filehandle line '$.'\n$!\n";
	my @info = keys(%info);

	#GATK strigency >> don't ask dont tell..
	my @usedInfoOrder;
	my @gatkInfoOrder= ('GT','AD','DP','GQ','PL');#req order(if present): GT:AD:DP:GQ:PL
	for my $gatkInfokey (@gatkInfoOrder){
		for my $infokey (@info){
			if($infokey eq $gatkInfokey){
				push(@usedInfoOrder,$infokey);
			}
		}
	}
	for my $infokey (@info){
		my $isInGatkFormat = 0;
		for my $gatkInfokey (@gatkInfoOrder){
			if($infokey eq $gatkInfokey){
				$isInGatkFormat = 1;
			}
		}
		if($isInGatkFormat == 0){
			push(@usedInfoOrder,$infokey);
		}
	}
	my $text = "";
	for my $tmp (@usedInfoOrder){
		if($info{$tmp} eq 'true' && $text eq ""){
			$text = $tmp;
		}elsif($info{$tmp} eq 'true'){
			$text = $text.';'.$tmp;
		}elsif($info{$tmp} ne 'true' && $text eq ""){
			$text = $tmp.'='.$info{$tmp};
		}elsif($info{$tmp} ne 'true'){
			$text = $text.';'.$tmp.'='.$info{$tmp};
		}	
	}
	#map{my $tmp = $_;$text = $text.';'.$tmp$info if(not($info{$tmp} eq 'true'));$text = $text.$tmp if($info{$tmp} eq 'true')}(@info);
	return $text;
}
sub formatReader {
	my @format = split(':', ${$_[0]});
	my @sampleInfos = @{$_[1]};
	my @samples = @{$_[2]};
	my $sample_count = 0;
	my %perSampleFormatData;
	for my $sample (@samples){
		$perSampleFormatData{'sampleNames'}{$sample_count} = $sample;
		$sample_count++;
		my @formatInfo = split(":",shift(@sampleInfos));
		for my $formatfield (@format){
			$perSampleFormatData{$sample}{$formatfield} = shift(@formatInfo);
		}
		my $tmp = 0;
		map{$perSampleFormatData{$sample}{'ADi'}{$tmp}=$_;$tmp++;}(split(",",$perSampleFormatData{$sample}{'AD'}))if($perSampleFormatData{$sample}{'GT'} ne './.');
	}
	$perSampleFormatData{'sampleCount'} = $sample_count;
	##while(@_){
	#	my @sampleData = split(':',shift @_);
	#	my $i=0;
	#	#map{$perSampleFormatData{$sample_count}{$format[$i]}=$_;$i++}(@sampleData);
	#	map{$formatData{$format[$i]}=$_;$i++}(@sampleData);
	#	#$sample_count++;
	##}
	#die "still sampledata left while handling file line '$.':'".join("",@_)."'\nperl error'$!'\n" if(join("",@_) ne "");
	return %perSampleFormatData;
}
sub formatPrinter {
	my @format = split(':', ${$_[0]});
	my %formatData = %{$_[1]};
	#my $sample_count = 0;
	my @sampletext = ();
	my $formatText = "";
	#GATK strigency >> don't ask dont tell..
	my @usedFormatOrder;
	my @gatkFormatOrder= ('GT','AD','DP','GQ','PL');#req order(if present): GT:AD:DP:GQ:PL
	for my $gatkFormatkey (@gatkFormatOrder){
		for my $formatkey (@format){
			if($formatkey eq $gatkFormatkey){
				push(@usedFormatOrder,$formatkey);
			}
		}
	}
	for my $formatkey (@format){
		my $isInGatkFormat = 0;
		for my $gatkFormatkey (@gatkFormatOrder){
			if($formatkey eq $gatkFormatkey){
				$isInGatkFormat = 1;
			}
		}
		if($isInGatkFormat == 0){
			push(@usedFormatOrder,$formatkey);
		}
	}
	##endof GATK stringency
	my $count = 0;
	while($count < $formatData{'sampleCount'}){
		my $perSampleTekst="";

		if($formatData{ $formatData{'sampleNames'}{$count} }{'GT'} ne './.'){
			$formatText = "";#or else nonsensical format tekst
			for my $tmp (@usedFormatOrder){
				warn "tmp:".$tmp if(not(defined($tmp)));
				warn "formatData:".$formatData{ $formatData{'sampleNames'}{$count} }{$tmp} if(not(defined($formatData{ $formatData{'sampleNames'}{$count} }{$tmp})));
				if($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} eq 'true' && $formatText eq ""){
					$formatText = $tmp;
					$perSampleTekst = $tmp;
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} eq 'true'){
					$formatText = $formatText.':'.$tmp;
					$perSampleTekst = $perSampleTekst.':'.$tmp;
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} ne 'true' && $formatText eq ""){
					$formatText = $tmp;
					$perSampleTekst = $formatData{ $formatData{'sampleNames'}{$count} }{$tmp};
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} ne 'true'){
					$formatText = $formatText.':'.$tmp;
					$perSampleTekst = $perSampleTekst.':'.$formatData{ $formatData{'sampleNames'}{$count} }{$tmp};
				}
			}
		}else{
			$perSampleTekst='./.';
		}
		push(@sampletext,$perSampleTekst);
		$count++;
	}
	
	return ($formatText, join("\t",@sampletext));
}
#sub arrayFraction {
#	my $DP = shift(@_);
#	my @division;
#	for my $EC (@_){
#		push(@division,$EC/$DP);
#	}
#	return @division;
#}
sub calcADvals {
	my $DP = shift(@_);
	my @EC = split(",",shift(@_));
	my @ADvals;
	$ADvals[0] = $DP;
	for my $EC (@EC){
		push(@ADvals,$EC);
		$ADvals[0] =$ADvals[0] - $EC;
	}
	return @ADvals;
}
sub getDifferentIndexesB {
	my $A = shift @_;
	my $B = shift @_;
	my @A = split '/',$A;
	my %B;
	for(split '/',$A){
		$B{$_}=0;
	};
	
	for my $alleleA (@A) {
		for my $alleleB (keys(%B)) {
			$B{$alleleB}++ if($alleleB eq $alleleA);
		}
	}
	my @Balt;
	for my $alleleB (keys(%B)) {
		push(@Balt,$alleleB) if($B{$alleleB} eq 0);
	}
	return @Balt;
}
