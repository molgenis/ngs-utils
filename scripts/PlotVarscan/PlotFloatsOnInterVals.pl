#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opts;
my $use =<<"END";
use:
perl $0 -d <human.dict> <varscan.called> #default also looks for <varscan.called>.homdels so enable the homdels output
perl $0 -d <human.dict> -r 1,1,1000000 <varscan.called> #by region
perl $0 -d <human.dict> <varscan.called> <vcf.table> #also plot f/z vals
#note: the 
my actual use:
terpstramm\@Biolinux:~/Documents/gccclusterdata/data/projects/exomeseq/privsmeta/privsmeta/CNV\$ 
for i in $(ls S12_193_*.called);
	do perl /home/terpstramm/workspace/PlotGenome/PlotFloatsOnInterVals.pl -v ../../forGerard/s12-193.refilter.snps.Z.vcf.table \$i > S12_193_14.called.Rscript;
	Rscript \$i.Rscript ;
done
END
#R needs to be installed
$opts{'d'} = '/home/terpstramm/Documents/gccclusterdata/data/resources/gatkbundle2.5/b37/human_g1k_v37.dict';
#$ARGV[0] = '/home/terpstramm/Documents/gccclusterdata/data/projects/exomeseq/privsmeta/privsmeta/CNV/S12_193_3.called';
getopts('r:d:v:', \%opts);

my $Dict = $opts{'d'};

my %ChrLenData = ReadDict($Dict);
my $interval = $opts{'r'};

if(not($interval)&&$opts{'v'}&&-e $opts{'v'}&&  -e $ARGV[0]){
	MakeVarscanDataTable($ARGV[0], %ChrLenData);
	MakeVCFDataTable($opts{'v'}, %ChrLenData);
	MakePDFrscriptWithVCF($ARGV[0],$opts{'v'});
	
}elsif(not($interval)&& -e $ARGV[0]){
	MakeVarscanDataTable($ARGV[0], %ChrLenData);
	MakePDFrscript($ARGV[0]);
}elsif(-e $ARGV[0]){
	MakeVarscanDataTable($ARGV[0], %ChrLenData);
	MakeVarscanDataTableByInterval($ARGV[0],$interval);
	MakePDFrscriptByInterval($ARGV[0],$interval);
}else{
	$use && die;
}

sub ReadDict {
	open(my $DictHandle,'<',$_[0]) or die "cannot read dictionary file '$_[0]' ";
	my %ChrLenData;
	my $RefIndex = 0;
	my $OffSet = 1;
	while(my $line = <$DictHandle>){
		if($line =~ m/^\@SQ/){
			my @tabdelim = split("\t", $line);
			my @SN = split(':',$tabdelim[1]);
			my @LN = split(':',$tabdelim[2]);
			$ChrLenData{'ChrLen'}{$SN[1]} = $LN[1];
			
			
			$ChrLenData{'ChrOff'}{$SN[1]} = $OffSet;
			$OffSet=$OffSet+$LN[1];
			
			$ChrLenData{'order'}{$RefIndex} = $SN[1];
			$RefIndex++;
		}
		
	}
	$ChrLenData{'REFcount'}=$RefIndex;
	close($DictHandle);
	return %ChrLenData;
}
sub MakeVarscanDataTable {
	open(my $Varscan2CopyCallerHandle,'<',$_[0]) or die "cannot read Varscan2CopyCaller file '$_[0]' ";
	open(my $Varscan2CopyCallerHomsHandle,'<',$_[0].'.homdels') or die "cannot read <Varscan2CopyCaller>.homdels file '$_[0]' ";
	open(my $Varscan2CopyCallerOutHandle,'>',$_[0].'.table');
	shift @_;
	my %ChrLenData = @_;
	my $header = <$Varscan2CopyCallerHandle>;
	print $Varscan2CopyCallerOutHandle 'offset'."\t".$header;
	#die Dumper(\%ChrLenData);
	while(my @tabdelim = split("\t",<$Varscan2CopyCallerHandle>)){
		#warn "$tabdelim[0]\t$tabdelim[1]\t$ChrLenData{'ChrOff'}{$tabdelim[0]}\n";
		print $Varscan2CopyCallerOutHandle join("\t",(($ChrLenData{'ChrOff'}{$tabdelim[0]}+$tabdelim[1]),@tabdelim));
	}
	<$Varscan2CopyCallerHomsHandle>;#skips header
	while(my @tabdelim = split("\t",<$Varscan2CopyCallerHomsHandle>)){
		print $Varscan2CopyCallerOutHandle join("\t",($ChrLenData{'ChrOff'}{$tabdelim[0]}+$tabdelim[1],@tabdelim));
	}
}
sub MakePDFrscript {
	my $dataTable = shift @_;
	my $shortname = $dataTable;
	$shortname =~ s/.called//g;
	
	my $rscript=<<"END";
#load data
VarscanDat <- read.table("$dataTable.table", sep="\\t", header=TRUE)

#median centering
median <- quantile(VarscanDat\$adjusted_log_ratio, .5)
VarscanDat\$adjusted_log_ratio <- VarscanDat\$adjusted_log_ratio - median
#truncate CNV values
truncate.CNV <- function(x){ if(x < -3) { return(-3) } else if(x > 2) {return(2)} else {return(x) }}
VarscanDat\$adjusted_log_ratio<-apply(X=t(VarscanDat\$adjusted_log_ratio),MARGIN=2,FUN=truncate.CNV)
########################################^why the t() funtion ?????????????!.... should i learn R explicitly instead of grabbing codes from the webs????
#pdf open
pdf(file="$dataTable.pdf", width=20, height=6)
#plot
plot(VarscanDat\$offset, VarscanDat\$adjusted_log_ratio, pch='.', lty=3, xlab="", ylab="log2(tumor/normal)", xaxt="n", main="$shortname" , ylim = c(-3,2))
abline(a=log(x=2,base=2),b=0, col="red")
abline(a=log(x=1.5,base=2),b=0, col="orange")
abline(a=log(x=1,base=2),b=0, col="green")
abline(a=log(x=0.5,base=2),b=0, col="blue")
abline(a=log(x=0.25,base=2),b=0, col="black")
legend("bottomright", title="copynumber relative to normal",c("4n","3n","2n","1n","0n"), fill=c("red","orange","green","blue","black"))
END
	my $i = 0;
	my @at; my @labels;
	while($i < $ChrLenData{'REFcount'}){
		push(@labels,$ChrLenData{'order'}{$i});
		push(@at,$ChrLenData{'ChrOff'}{ $ChrLenData{'order'}{$i} });
		$i++;	
	}
	$rscript=$rscript.'axis(1, at=c("'.join('","',@at).'"),labels=c("'.join('","',@labels).'"), las=2)';
	$rscript=$rscript.<<"END";

#pdf close
dev.off()
END
	print $rscript;
}
sub MakeVarscanDataTableByInterval {
	open(my $Varscan2CopyCallerHandle,'<',$_[0]) or die "cannot read Varscan2CopyCaller file '$_[0]' ";
	open(my $Varscan2CopyCallerHomsHandle,'<',$_[0].'.homdels') or die "cannot read Varscan2CopyCaller file '$_[0]' ";
	open(my $Varscan2CopyCallerOutIntervalHandle,'>',$_[0].'.interval.table');
	shift @_;
	my ($chrom, $start,$end )= split(',',shift(@_));
	#my %ChrLenData = @_;
	my $header = <$Varscan2CopyCallerHandle>;
	print $Varscan2CopyCallerOutIntervalHandle 'offset'."\t".$header;
	#die Dumper(\%ChrLenData);
	while(my @tabdelim = split("\t",<$Varscan2CopyCallerHandle>)){
		#warn "$tabdelim[0]\t$tabdelim[1]\t$ChrLenData{'ChrOff'}{$tabdelim[0]}\n";
		if($tabdelim[0] eq $chrom &&$tabdelim[1] > $start&& $tabdelim[1] < $end){
			print $Varscan2CopyCallerOutIntervalHandle join("\t",(($tabdelim[1]-$start),@tabdelim));
		}
	}
	<$Varscan2CopyCallerHomsHandle>;#skips header
	while(my @tabdelim = split("\t",<$Varscan2CopyCallerHomsHandle>)){
		if($tabdelim[0] eq $chrom &&$tabdelim[1] > $start&& $tabdelim[1] < $end){
			print $Varscan2CopyCallerOutIntervalHandle join("\t",(($tabdelim[1]-$start),@tabdelim));
		}
	}
}
sub MakePDFrscriptByInterval {
	my $dataTable = shift @_;
	my $shortname = $dataTable;
	$shortname =~ s/.called//g;
	my ($chrom, $start,$end )= split(',',shift(@_));
	my $rscript=<<"END";
#load data
VarscanDat <- read.table("$dataTable.table", sep="\\t", header=TRUE)
VarscanIntervalDat <- read.table("$dataTable.interval.table", sep="\\t", header=TRUE)

#median centering
median <- quantile(VarscanDat\$adjusted_log_ratio, .5)
VarscanIntervalDat\$adjusted_log_ratio <- VarscanIntervalDat\$adjusted_log_ratio - median
#pdf open
pdf(file="$dataTable.pdf", width=20, height=6)
plot(VarscanIntervalDat\$offset, VarscanIntervalDat\$adjusted_log_ratio, pch='.', lty=3, xlab="", ylab="log2(tumor/normal)", xaxt="n", main="$shortname")
abline(a=log(x=2,base=2),b=0, col="red")
abline(a=log(x=1.5,base=2),b=0, col="orange")
abline(a=log(x=1,base=2),b=0, col="green")
abline(a=log(x=0.5,base=2),b=0, col="blue")
abline(a=log(x=0.25,base=2),b=0, col="black")
legend("bottomright", title="copynumber relative to normal",c("4n","3n","2n","1n","0n"), fill=c("red","orange","green","blue","black"))
END
	my $i = 0;
	my @at; my @labels;
	$rscript=$rscript.'axis(1, at=c("'.join('","',(1,sprintf ('%.0f',($end-$start)/10*1),sprintf ('%.0f',($end-$start)/10*2),sprintf ('%.0f',($end-$start)/10*3),sprintf ('%.0f',($end-$start)/10*4),sprintf ('%.0f',($end-$start)/10*5),sprintf ('%.0f',($end-$start)/10*6),sprintf ('%.0f',($end-$start)/10*7),sprintf ('%.0f',($end-$start)/10*8),sprintf ('%.0f',($end-$start)/10*9), $end-$start ))
	.'"),labels=c("'.join('","',($start,sprintf ('%.0f',$start+  ($end-$start)/10*1),sprintf ('%.0f',$start+($end-$start)/10*2),sprintf ('%.0f',$start+($end-$start)/10*3),sprintf ('%.0f',$start+($end-$start)/10*4),sprintf ('%.0f',$start+($end-$start)/10*5),sprintf ('%.0f',$start+($end-$start)/10*6),sprintf ('%.0f',$start+($end-$start)/10*7),sprintf ('%.0f',$start+($end-$start)/10*8),sprintf ('%.0f',$start+($end-$start)/10*9),$end)).'"), las=2)'."\n";
	$rscript=$rscript.<<"END";

#pdf close
dev.off()
END
	print $rscript;
}
sub MakeVCFDataTable {
	open(my $VcfTableHandle,'<',$_[0]) or die "cannot read Snv table file '$_[0]' ";
	#open(my $VcfTableHandle,'<',$_[0].'.indels.table') or die "cannot read indels file '$_[0]' ";
	open(my $VcfTableHandleOutHandle,'>',$_[0].'.table');
	shift @_;
	my %ChrLenData = @_;
	my $header = <$VcfTableHandle>;
	print $VcfTableHandleOutHandle 'Offsets'."\t".$header;
	#die Dumper(\%ChrLenData);
	while(my @tabdelim = split("\t",<$VcfTableHandle>)){
		#warn "$tabdelim[0]\t$tabdelim[1]\t$ChrLenData{'ChrOff'}{$tabdelim[0]}\n";
		print $VcfTableHandleOutHandle join("\t",(($ChrLenData{'ChrOff'}{$tabdelim[0]}+$tabdelim[1]),@tabdelim));
	}
	#<$VcfTableHandle>;#skips header
	#while(my @tabdelim = split("\t",<$VcfTableHandle>)){
	#	print $VcfTableHandleOutHandle join("\t",($ChrLenData{'ChrOff'}{$tabdelim[0]}+$tabdelim[1],@tabdelim));
	#}
}
sub MakePDFrscriptWithVCF {
	my $dataTable = shift @_;
	my $shortname = $dataTable;
	$shortname =~ s/.called//g;
	my $variantTable = shift @_;
	
	
	my $rscript=<<"END";
#load data
VarscanDat <- read.table("$dataTable.table", sep="\\t", header=TRUE)
VariantDat <- read.table("$variantTable.table", sep="\\t", header=TRUE)

#median centering
median <- quantile(VarscanDat\$adjusted_log_ratio, .5)
VarscanDat\$adjusted_log_ratio <- VarscanDat\$adjusted_log_ratio - median
#pdf open
pdf(file="$dataTable.pdf", width=20, height=8)
#plot
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE),
    heights=c(2,1,1))
plot(VarscanDat\$offset, VarscanDat\$adjusted_log_ratio, pch='.', lty=3, xlab="", ylab="log2(tumor/normal)", xaxt="n", main="$shortname" , ylim = c(-3,2))
abline(a=log(x=2,base=2),b=0, col="red")
abline(a=log(x=1.5,base=2),b=0, col="orange")
abline(a=log(x=1,base=2),b=0, col="green")
abline(a=log(x=0.5,base=2),b=0, col="blue")
abline(a=log(x=0.25,base=2),b=0, col="black")
legend("bottomright", title="copynumber relative to normal",c("4n","3n","2n","1n","0n"), fill=c("red","orange","green","blue","black"))

END
	my $i = 0;
	my @at; my @labels;
	while($i < $ChrLenData{'REFcount'}){
		push(@labels,$ChrLenData{'order'}{$i});
		push(@at,$ChrLenData{'ChrOff'}{ $ChrLenData{'order'}{$i} });
		$i++;	
	}
	$rscript=$rscript.'axis(1, at=c("'.join('","',@at).'"),labels=c("'.join('","',@labels).'"), las=2)';
	$rscript=$rscript.<<"END";

plot(VariantDat\$Offsets,VariantDat\$${shortname}.F, pch='.', lty=3, xlab="", ylab="log2(tumor/normal)", xaxt="n", main="F vals")
abline(a=1,b=0, col="red")
abline(a=0,b=0, col="red")
abline(a=.5,b=0, col="green")
plot(VariantDat\$Offsets,VariantDat\$${shortname}.Z, pch='.', lty=3, xlab="", ylab="log2(tumor/normal)", xaxt="n", main="Z vals")
abline(a=0,b=0, col="green")

#pdf close
dev.off()
END
	print $rscript;
}