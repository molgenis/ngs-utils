use warnings;
use strict;
use Data::Dumper;

die("Usage: cat dosage | $0 <name.lst>\n") if (@ARGV == 0);
my $fn = shift(@ARGV);
my %hash = ();
open(FILE, $fn) || die;
while (<FILE>){
	chomp;
	$hash{$_} = 1;
};
warn Dumper(\%hash);
close(FILE);
my $flag = 0;
while (<>) {
	my $rsID=substr $_,0,index($_,"\t");
	if(defined($hash{$rsID})){ 
		print; 
	};
}
