#!/usr/bin/perl

#
# Copy data from the (Groningen GCC) cluster to the (bbmri.nl VO on) Grid. 
#

#
# Initialize evironment
#
use strict;
use Getopt::Std;
use Log::Log4perl qw(:easy);

my %log_levels = (
	'ALL'	=> $ALL,
	'TRACE'	=> $TRACE,
	'DEBUG'	=> $DEBUG,
	'INFO'	=> $INFO,
	'WARN'	=> $WARN,
	'ERROR'	=> $ERROR,
	'FATAL'	=> $FATAL,
	'OFF'	=> $OFF,
);

#
# Get options.
#
my %opts;
Getopt::Std::getopts('t:f:d:p:c:l:', \%opts);

my $type			= $opts{'t'};
my $list_file		= $opts{'f'};
my $destination_dir = $opts{'d'}; # srm://srm.grid.sara.nl/pnfs/grid.sara.nl/data/bbmri.nl/....
my $dcache_proxy	= $opts{'p'};
my $certificate_dir	= $opts{'c'};
my $log_level		= $opts{'l'};

#
# Configure logging.
#
# Provides default if user did not specify log level:
$log_level = (defined($log_level) ? $log_level : 'INFO');
# Reset log level to default if user specified illegal log level.
$log_level = (defined($log_levels{$log_level}) ? $log_levels{$log_level} : $log_levels{'INFO'});
#Log::Log4perl->init('log4perl.properties');
Log::Log4perl->easy_init(
	#{ level    => $log_level,
	#  file     => ">>copy2grid.log",
	#  layout   => '%F{1}-%L-%M: %m%n' },
	{ level    => $log_level,
	  file     => "STDOUT",
	  layout   => '%d L:%L %p> %m%n' },
);
my $logger = Log::Log4perl::get_logger();

#
# Check options.
#
if ($list_file =~ /^$/) {
	_Usage();
	exit;
}
unless (-f $list_file && -r $list_file) {
	_Usage();
	exit;
}
unless ($destination_dir =~ /\/$/) {
	$destination_dir .= '/';
}
unless (defined($type)) {
	if ($destination_dir =~ m/srm:\/\//) {
		$type = 'c2g';
	} elsif ($destination_dir =~ m/\//) {
		$type = 'g2c';
	} else {
		$logger->fatal('Destination dir in unrecognised format...');
		exit;
	}
}

#
# Start job.
#
$logger->info('Starting...');

#
# Create filehandles.
#
my $list_fh;

eval {
	open($list_fh, "<$list_file");
};
if ($@) {
	$logger->fatal('Cannot create filehandle: ' . $@);
	exit;
}

my $found_dcache_srmcp = 0;
my $dcache_srmcp = `which srmcp`;
if ($dcache_srmcp =~ m/\/srmcp$/) {
	chomp($dcache_srmcp);
	$logger->trace('Found srmcp command in: ' . $dcache_srmcp);
	$found_dcache_srmcp = 1;
} else {
	$logger->trace('Could not find srmcp command in $PATH.');
	$logger->trace('Will try to find srmcp command in the "usual suspect" places of the GCC cluster.');
	my $gcc_home;
	if (defined($ENV{'GCC_HOME'})) {
		$gcc_home = $ENV{'GCC_HOME'};
	} else {
		$logger->fatal('${GCC_HOME} environment variable is not defined.');
		exit 1;
	}
	$dcache_srmcp = $gcc_home. '/tools/dcache/opt/d-cache/srm/bin/srmcp';
	if (-f $dcache_srmcp && -r $dcache_srmcp && -X $dcache_srmcp) {
		$logger->trace('Found srmcp command in: ' . $dcache_srmcp);
		$found_dcache_srmcp = 1;
	} else {
		$logger->trace('Could not find srmcp command in the "usual suspect" places on the GCC cluster either.');
	}
}

unless ($found_dcache_srmcp) {
	$logger->fatal('Could not find srmcp command!');
	exit;
}

my $cmd_base = $dcache_srmcp;
$cmd_base .= ' -x509_user_proxy=' . $dcache_proxy;
$cmd_base .= ' -x509_user_trusted_certificates=' . $certificate_dir;
$cmd_base .= ' -server_mode=passive';
# Disabled SRM V2: SE srm.grid.sara.nl supports SRM protocol version 2, but many of the other SEs do not.
#$cmd_base .= ' -2';

my $total_file_count = 0;

#
# Parse file that lists the files that need to be copied.
#
while (my $line = <$list_fh>) {
	my $file_path = $line;
	$file_path =~ s/[\n\r\f]+$//;
	my $file_name;
	my $cmd = $cmd_base;
	if ($file_path =~ /\/([^\/]+)$/) {
		$file_name = $1;
	} else {
		$file_name = $file_path;
		#$logger->fatal('Cannot parse file name from file path');
		#exit;
	}
	if ($type eq 'g2c') {
		$cmd .= ' ' . $file_path;		
		$cmd .= ' file:///' . $destination_dir . $file_name;
	} elsif ($type eq'c2g') {
		$cmd .= ' file:///' . $file_path;
		$cmd .= ' ' . $destination_dir . $file_name;
	}
	$logger->debug('Executing command: ' . $cmd);
	my $result = `$cmd 2>&1`;
	$logger->info($result);
	$total_file_count++;
}

#
# Close filehandles.
#
close($list_fh);

$logger->info('Processed ' . $total_file_count . ' files.');
$logger->info('Finished!');

#
##
### Subs.
##
#

sub _Usage {

	print "\n";
	print 'copy2grid.pl - Copies data from the (Groningen GCC) cluster to the (bbmri.nl VO on) Grid'."\n";
	print '               Assumes DCache is either in your $PATH or in the $GCC_HOME/tools/dcache dir on the Groningen GCC cluster'."\n";
	print "\n";
	print "Usage:\n";
	print "\n";
	print "   copy.grid.pl options\n";
	print "\n";
	print "Available options are:\n";
	print "\n";
	print "   -t [g2c|c2g]       Type or tranfser mode:\n";
	print "                       * c2g - To transfer from a Cluster to the Grid.\n";
	print "                       * g2c - To transfer from the Grid to a Cluster.\n";
	print "   -f [file]          Input file listing files that need to be transferred.\n";
	print "   -d [path]          Directory the files need to be tansfered to.\n";
	print "                      If the mode is cluster2grid a grid dir on an SRM.\n";
	print "                      For example srm://srm.grid.sara.nl/pnfs/grid.sara.nl/data/bbmri.nl/batch_x\n";
	print "                       * If the transfer mode is c2g, a grid dir on an SRM.\n";
	print "                       * If the transfer mode is g2c, an absolute dir on a local mount.\n";
	print "   -p [file]          Proxy certificate (your cert).\n";
	print "   -c [dir]           Certificate directory containing all the grid-security certs (certs from other machines in the grid).\n";
	print "   -l [LEVEL]         Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.\n";
	print "\n";
	exit;

}
