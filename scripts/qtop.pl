#!/usr/bin/env perl

# ============================================================================
#  Original Copyright 2002, 2003, 2004 University of Southern California
#  Later modifications Copyright 2014 University Medical Center Groningen
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  Original version of this software may be found at:
#      http://www-rcf.usc.edu/~garrick/perl-PBS
#      Please send comments to garrick@usc.edu.
#      Site appears to be down...
# ============================================================================

#
# Configurable defaults.
#
# Don't change the defaults here.
# Make a config file instead; see qtop(1).
#
my $columns         = 46;       # Columns in grid
my $sleeptime       = 30;       # Seconds between refreshes
my $colorize        = 1;        # 1 or 0
my $color_by        = 'job';    # Assign color per job, user or queue.
my $show_summary    = 1;        # 1 or 0
my $compact_summary = 1;        # 1 or 0
my $show_grid       = 1;        # 1 or 0
my $sorted_grid     = 1;        # 1 or 0
my $show_queue      = 1;        # 1 or 0
my $show_qqueue     = 1;        # 1 or 0
my $show_jobs       = 1;        # 1 or 0
my $show_user       = "all";    # Show only jobs for specific users
my $show_node       = "all";    # Show only jobs for specific nodes jobs
my @show_cpu        = ("0");    # List of cpu numbers. No longer used and hardcoded to CPU0.
my @host            = ();       # Hostname of your scheduler
                                # This is not necessarily the same as the user interface (UI) host
                                # on which you login to monitor/submit jobs. When empty we will try
                                # $ENV{"PBS_DEFAULT"}, SERVERHOST from showconfig or localhost.
my $maxrows         = 3000;     # maximum number of rows
my $maxcolumns      = 145;      # maximum number of columns
my $maxnodegrid     = 2;        # maximum number of CPUs on a node before it gets its own grid.

my $qmgr     = "/usr/local/pbs/bin/qmgr";
my $qstat    = "/usr/local/pbs/bin/qstat";
my $pbsnodes = "/usr/local/pbs/bin/pbsnodes";

#########################################################
#   Nothing else to adjust below here
#########################################################

# enable warnings if running under testing.
if (-d $INC[0]) {
	$^W = 1;
}

use strict;
use vars qw/$VERSION/;
use Curses;

$VERSION = "5.1";

#
# Initialize global vars.
#
my %Job_of_letter;
my @Colors = ();
my %searchobject  = ();

#
# Define which characters can be used to show jobs or node state in grid view.
# Characters used to indicate state must not be used as masterletters for jobs.
#
my %state_characters = (
	unknown => '?',
	busy    => '@',
	down    => 'X',
	idle    => '.',
	offline => '0',
	other   => '!'
);
#my $masterletters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
my $masterletters = "ABCDEHIKLNOSTUVYZabcdeiklnopqsuvwxyz-+=<>|/\_~";
my $letters       = $masterletters;
my $underline     = 0;


my ($y, $x, $Y, $X, $py, $px, $ly, $lx, $subY, $subX) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

readrc("/etc/qtoprc");
readrc("$ENV{HOME}/.qtoprc");

#
# Configure and find commandline utils.
#
-e $qmgr     or chomp($qmgr     = `which qmgr 2>/dev/null`);
-e $qmgr     or die "qmgr: Command not found\n";
-e $qstat    or chomp($qstat    = `which wqstat 2>/dev/null`);
-e $qstat    or chomp($qstat    = `which qstat 2>/dev/null`);
-e $qstat    or die "qstat: Command not found\n";
-e $pbsnodes or chomp($pbsnodes = `which pbsnodes 2>/dev/null`);
-e $pbsnodes or die "pbsnodes: Command not found\n";

# argument processing
my @argvhosts = ();
while (my $arg = shift @ARGV) {
	if ($arg eq '-c') {
		$columns = shift @ARGV;
		$columns =~ /^\d+$/ or $ARGV[0] = '-h';
		$columns > 0 or $ARGV[0] = '-h';
	} elsif ($arg eq '-s') {
		$sleeptime = shift @ARGV;
		$sleeptime =~ /^\d+$/ or $ARGV[0] = '-h';
		$sleeptime > 0 or $ARGV[0] = '-h';
	} elsif ($arg eq '-m') {
		$maxnodegrid = shift @ARGV;
		$maxnodegrid =~ /^\d+$/ or $ARGV[0] = '-h';
		$maxnodegrid > 0 or $ARGV[0] = '-h';
	} elsif ($arg eq '-C') {
		$colorize = !$colorize;
	} elsif ($arg eq '-S') {
		$show_summary = !$show_summary;
	} elsif ($arg eq '-G') {
		$show_grid = !$show_grid;
	} elsif ($arg eq '-Q') {
		$show_queue = !$show_queue;
	} elsif ($arg eq '-t') {
		$show_qqueue = !$show_qqueue;
	} elsif ($arg eq '-J') {
		$show_jobs = !$show_jobs;
	} elsif ($arg eq '-u') {
		if (defined $ARGV[0] and $ARGV[0] =~ /^([^-]+)/) {
			$show_user = join(' ', split(',', shift));
		} else {
			$ARGV[0] = '-h';
		}
	} elsif ($arg eq '-n') {
		if (defined $ARGV[0] and $ARGV[0] =~ /^([^-]+)/) {
			$show_node = join(' ', split(',', shift));
		} else {
			$ARGV[0] = '-h';
		}
	} elsif ($arg =~ /^-(\d+)$/) {
		@show_cpu = split(//, $1);
	} elsif ($arg =~ /^@(.*)/) {
		push(@argvhosts, $1);
	} elsif ($arg eq '-V') {
		print "qtop (p.k.a. pbstop) $VERSION\nCopyright 2002, 2003, 2004 University of Southern California\n";
		exit(0);
	} else {
		print "Usage:  qtop [-c columns] [-s seconds] [-m numcpus] [options] [\@host ...]\n";
		print "   Version: $VERSION\n";
		print "   Copyright 2002, 2003, 2004 University of Southern California\n";
		print "   garrick\@usc.edu http://www-rcf.usc.edu/~garrick/pbstop\n";
		print "   grep FIXME `which qtop` if you want to help out\n\n";
		print "   -s  seconds between refreshes\n";
		print "   -c  number of columns to display in the grid\n";
		print "   -m  max number of cpus in a node before it gets its own grid\n";
		print "   -u  show only a user's jobs\n";
		print "   -n  show only a node's jobs\n";
		print "   -C  toggle colorization\n";
		print "   -S  toggle state summary display\n";
		print "   -G  toggle grid display\n";
		print "   -Q  toggle queue display\n";
		print "   -t  toggle showing queued jobs in queue display\n";

		#print "   -[0-9]...  cpu numbers for grid display\n";
		print "   -J  toggle jobs in grid display\n";
		print "   -V  print version and exit\n";
		exit(1);
	}
}
if (scalar @argvhosts > 0) {
	@host = @argvhosts;
}
undef @argvhosts;

#
# If the scheduler was not specified,
# try to get it's name from the evironment, defaults, usual suspects, etc.
#
unless (defined($host[0])) {
	if (defined($ENV{"PBS_DEFAULT"})) {
		$host[0] = $ENV{"PBS_DEFAULT"};
	} elsif (`showconfig | fgrep SERVERHOST` =~ m/SERVERHOST\s+([^\s]+)/) {
		$host[0] = $1;
	} else {
		$host[0] = `hostname`;
	}
}
chomp(@host);

if ($show_user eq "all") {
	$show_user = 0;
} elsif ($show_user =~ /\bme\b/) {
	$show_user =~ s/\bme\b/$ENV{USER}/;
}
if ($show_node eq "all") {
	$show_node = 0;
}

use vars qw/$SIGWINCH/;
$SIGWINCH = 0;
$SIG{'WINCH'} = sub { $SIGWINCH = 1; };
$SIG{'INT'}  = sub { endwin; exit(0); };
$SIG{'TERM'} = sub { endwin; exit(0); };

#FIXME# Can someone tell me how to use filter() correctly?
-t STDOUT or filter();

# Is this portable?
my $CTRL_B = chr(ord("B") - ord("@"));
my $CTRL_F = chr(ord("F") - ord("@"));
my $CTRL_L = chr(ord("L") - ord("@"));
my $CTRL_G = chr(ord("G") - ord("@"));
my $CTRL_H = chr(ord("H") - ord("@"));

initscr;
cbreak;
noecho;
getmaxyx($Y, $X);
start_color;
$colorize = $colorize && has_colors();

my $pr = 0;

#
# This is every possible color combo list below.
# Over time, I've commented out color pairs that don't look very good.
# If your eyes disagree with my eyes, you are free to play around with this list.
# But don't forget... only the first $COLOR_PAIRS uncommented combos apply.
# $COLOR_PAIRS is set by your curses implementation.
# qtop's help screen (hit 'h' in qtop) will tell you the value of $COLOR_PAIRS.
#
init_pair(++$pr, COLOR_RED,     COLOR_BLACK);
init_pair(++$pr, COLOR_GREEN,   COLOR_BLACK);
init_pair(++$pr, COLOR_YELLOW,  COLOR_BLACK);
init_pair(++$pr, COLOR_BLUE,    COLOR_BLACK);
init_pair(++$pr, COLOR_MAGENTA, COLOR_BLACK);
init_pair(++$pr, COLOR_CYAN,    COLOR_BLACK);
init_pair(++$pr, COLOR_WHITE,   COLOR_BLACK);

#init_pair( ++$pr, COLOR_BLACK,   COLOR_BLACK );
#init_pair( ++$pr, COLOR_RED,     COLOR_WHITE );
#init_pair( ++$pr, COLOR_GREEN,   COLOR_WHITE );
#init_pair( ++$pr, COLOR_YELLOW,  COLOR_WHITE );
#init_pair( ++$pr, COLOR_BLUE,    COLOR_WHITE );
#init_pair( ++$pr, COLOR_MAGENTA, COLOR_WHITE );
#init_pair( ++$pr, COLOR_CYAN,    COLOR_WHITE );
#init_pair( ++$pr, COLOR_WHITE,   COLOR_WHITE );
init_pair(++$pr, COLOR_BLACK, COLOR_WHITE);

#init_pair( ++$pr, COLOR_RED,     COLOR_YELLOW );
#init_pair( ++$pr, COLOR_GREEN,   COLOR_YELLOW );
#init_pair( ++$pr, COLOR_YELLOW,  COLOR_YELLOW );
#init_pair( ++$pr, COLOR_BLUE,    COLOR_YELLOW );
init_pair(++$pr, COLOR_MAGENTA, COLOR_YELLOW);

#init_pair( ++$pr, COLOR_CYAN,    COLOR_YELLOW );
#init_pair( ++$pr, COLOR_WHITE,   COLOR_YELLOW );
#init_pair( ++$pr, COLOR_BLACK,   COLOR_YELLOW );
init_pair(++$pr, COLOR_RED, COLOR_CYAN);

#init_pair( ++$pr, COLOR_GREEN,   COLOR_CYAN );
init_pair(++$pr, COLOR_YELLOW, COLOR_CYAN);

#init_pair( ++$pr, COLOR_BLUE,    COLOR_CYAN );
init_pair(++$pr, COLOR_MAGENTA, COLOR_CYAN);

#init_pair( ++$pr, COLOR_CYAN,    COLOR_CYAN );
init_pair(++$pr, COLOR_WHITE,  COLOR_CYAN);
init_pair(++$pr, COLOR_BLACK,  COLOR_CYAN);
init_pair(++$pr, COLOR_RED,    COLOR_MAGENTA);
init_pair(++$pr, COLOR_GREEN,  COLOR_MAGENTA);    # current 16th
init_pair(++$pr, COLOR_YELLOW, COLOR_MAGENTA);
init_pair(++$pr, COLOR_BLUE,   COLOR_MAGENTA);

#init_pair( ++$pr, COLOR_MAGENTA, COLOR_MAGENTA );
init_pair(++$pr, COLOR_CYAN,  COLOR_MAGENTA);
init_pair(++$pr, COLOR_WHITE, COLOR_MAGENTA);

#init_pair( ++$pr, COLOR_BLACK,   COLOR_MAGENTA );
#init_pair( ++$pr, COLOR_RED,     COLOR_RED );
init_pair(++$pr, COLOR_GREEN,   COLOR_RED);
init_pair(++$pr, COLOR_YELLOW,  COLOR_RED);
init_pair(++$pr, COLOR_BLUE,    COLOR_RED);
init_pair(++$pr, COLOR_MAGENTA, COLOR_RED);
init_pair(++$pr, COLOR_CYAN,    COLOR_RED);
init_pair(++$pr, COLOR_WHITE,   COLOR_RED);
init_pair(++$pr, COLOR_BLACK,   COLOR_RED);
init_pair(++$pr, COLOR_RED,     COLOR_GREEN);

#init_pair( ++$pr, COLOR_GREEN,   COLOR_GREEN );
#init_pair( ++$pr, COLOR_YELLOW,  COLOR_GREEN );
init_pair(++$pr, COLOR_BLUE, COLOR_GREEN);

#init_pair( ++$pr, COLOR_MAGENTA, COLOR_GREEN );
init_pair(++$pr, COLOR_CYAN,  COLOR_GREEN);
init_pair(++$pr, COLOR_WHITE, COLOR_GREEN);
init_pair(++$pr, COLOR_BLACK, COLOR_GREEN);

sub init_colors {
	return (1 .. ($COLOR_PAIRS - 1 > $pr ? $pr : $COLOR_PAIRS - 1));
}

my $pad = newpad($maxrows, $maxcolumns);
my $cmdwin = newwin(1, $X - 1, $Y - 1, 0);
keypad($cmdwin, 1);
my $subwin = 0;

main_loop(\@host);

#   The original color set that I actually spent some time planning
#        "\033[07;34m",    "\033[07;35m",    "\033[07;36m",
#        "\033[07;37m",    "\033[01;37m",    "\033[35m",
#        "\033[36m",       "\033[37m",       "\033[34m",
#        "\033[33m",       "\033[32m",       "\033[01;36;45m",
#        "\033[01;30;47m", "\033[01;30;46m", "\033[36;45m",
#        "\033[30;47m",    "\033[30;46m",    "\033[01;33m",
#        "\033[01;34m",    "\033[01;35m",    "\033[01;36m",
#        "\033[01;31m",    "\033[01;32m",

###############################################################
## All subroutines below here
###############################################################

# main_loop() will 1) gather all of our data from pbs_server, 2) prep it a bit
# in letterize() and colorize(), 3) call update_display to draw our pretty
# grids and stuff, and finally calls 4) top_sleep which is where we spend most
# of our time.

# 1) We gather data by parsing the output of tools like qmgr and qstat.
#    Since this data is kept between cycles around the main loop, we take some care to remove old data.
#    The result is two large structures, one for Jobs and one for Nodes, which are used in the rest of the program.

# 2) letterize() and colorize() are fairly unexciting,
#    but they do assign letters and colors to each running job.
#    This info is stored in the large Job structure.
#    letterize() has probably the only original code left.

# 3) update_display(), by itself, is pretty coring.
#    It calls the functions responsible for the summary, colorful grid, and the job listing at the bottom.
#    show_grid() is pretty exciting; it first finds every node over $maxnodegrid,
#    calls them "timesharing" and shoves them aside, draws a big colorful grid for what's left,
#    and draws another colorful grid for the timesharing nodes.

# 4) top_sleep() is a big giant mess. It is far too monolithic.
#    If anyone wants chop it up a bit, feel free to send patches!
#    Anyways, it loops around on user input until the time expires and it is time to return back up to main_loop().
#    In the meantime, it does everything the user requests,
#    including griding through the main data structures looking for stuff.
#    All of the code responsible for moving around the giant pad is here.

sub main_loop {
	my $host = shift;
	my $maxlen;
	my %Nodes;
	my %Jobs;
	my %Users;
	my %Queues;
	my %State_count;

	# Main event loop.
	while (1) {

		%State_count          = ();
		$State_count{_nodes}  = 0;
		$State_count{_anodes} = 0;
		$State_count{_procs}  = 0;
		$State_count{_aprocs} = 0;
		$State_count{_mprocs} = 0;
		$State_count{_rjobs}  = 0;
		$State_count{_njobs}  = 0;

		foreach my $server (@$host) {

			get_info_cmdline($server, \%Nodes, \%Jobs, \%State_count);

			# trim out old nodes that are no longer seen
			foreach my $node (keys(%{ $Nodes{$server} })) {
				if ($Nodes{$server}{$node}{seen} != 1) {
					delete $Nodes{$server}{$node};
				} else {
					$Nodes{$server}{$node}{seen} = 0;
				}
			}
		}

		# trim out old jobs that are no longer seen
		foreach my $job (keys(%Jobs)) {
			if (   !exists $Jobs{$job}{seen}
				or !defined $Jobs{$job}{seen}
				or $Jobs{$job}{seen} != 1)
			{
				delete $Jobs{$job};
			} else {
				$Jobs{$job}{seen} = 0;
			}
		}

		#$maxlen |= getmaxkeylen(\%Nodes);
		$maxlen = getmaxkeylen(\%Nodes);
		if ($maxlen < 4) {
			$maxlen = 4;
		}

		letterize(\%Jobs);
		colorize(\%Jobs, \%Users, \%Queues);

		update_display(\%State_count, $State_count{"_mprocs"}, $maxlen, \%Nodes, \%Jobs, \%Users, \%Queues);
		-t STDOUT or do { endwin; exit; };
		top_sleep(\%State_count, $State_count{"_mprocs"}, $maxlen, \%Nodes, \%Jobs, \%Users, \%Queues);
	}
}

#sub get_info_modPBS {
#	my ( $server, $Nodes, $Jobs, $State_count ) = @_;
#
#	my $con = pbs_connect($server);
#	if ( $con <= 0 ) {
#		printwarning("Connect to $server failed: $PBS::pbs_errno\n");
#		return;
#	}
#	my $qmgr = pbs_statnode( $con, undef, undef, undef );
#	my $qstat = pbs_statjob( $con, undef, undef, undef );
#	pbs_disconnect($con);
#
#	my $jobs;
#	my $status;
#	my $statuses;
#	my $node;
#	my $job;
#	my $name;
#	my $value;
#
#	my $node = "";
#
#	foreach ( @{$qmgr} ) {
#		$node = $_->{name};
#		delete $Nodes->{$server}{$node};
#		$Nodes->{$server}{$node}{seen} = 1;
#		$State_count->{_nodes}++;
#
#		foreach ( @{ $_->{attribs} } ) {
#			$name  = $_->{name};
#			$value = $_->{value};
#
#			$name eq $PBS::ATTR_NODE_np and do {
#				$Nodes->{$server}{$node}{np} = $value;
#				$State_count->{_procs} += $value;
#				$State_count->{_mprocs} =
#				    $State_count->{_mprocs} > $value
#				  ? $State_count->{_mprocs}
#				  : $value;
#			};
#			$name eq $PBS::ATTR_NODE_state and do {
#				$Nodes->{$server}{$node}{state} = $value;
#				$State_count->{$value}++;
#			};
#			$name eq $PBS::ATTR_NODE_properties
#			  and $Nodes->{$server}{$node}{properties} = $value;
#			$name eq $PBS::ATTR_NODE_ntype
#			  and $Nodes->{$server}{$node}{ntype} = $value;
#			$name eq $PBS::ATTR_NODE_jobs and do {
#				$State_count->{"_anodes"}++;
#				foreach my $job ( split( /, /, $value ) ) {
#					if ( $job =~ m{(\d+)/(\d+)} ) {
#						$Nodes->{$server}{$node}{job}{$1} = $2;
#						$State_count->{"_aprocs"}++;
#					}
#				}
#			};
#			$name eq $PBS::ATTR_NODE_status and do {
#				foreach my $status ( split( /,/, $value ) ) {
#					if ( $status =~ m{(.+)=(.+)} ) {
#						$Nodes->{$server}{$node}{status}{$1} = $2;
#					}
#				}
#			};
#		}
#	}
#
#	my $job;
#	foreach ( @{$qstat} ) {
#		( $job = $_->{name} ) =~ s/\..*//;
#		$Jobs->{$job}{seen} = 1;
#		foreach ( @{ $_->{attribs} } ) {
#			$_->{name} eq $PBS::ATTR_server
#			  and $Jobs->{$job}{server} = $_->{value};
#			$_->{name} eq $PBS::ATTR_owner
#			  and ( $Jobs->{$job}{user} = $_->{value} ) =~ s/\@.*//;
#			$_->{name} eq $PBS::ATTR_queue
#			  and $Jobs->{$job}{queue} = $_->{value};
#			$_->{name} eq $PBS::ATTR_N and $Jobs->{$job}{jname} = $_->{value};
#			$_->{name} eq "$PBS::ATTR_l.nodect"
#			  and $Jobs->{$job}{ncount} = $_->{value};
#			$_->{name} eq "$PBS::ATTR_l.walltime"
#			  and $Jobs->{$job}{reqt} = $_->{value};
#			$_->{name} eq $PBS::ATTR_state
#			  and $Jobs->{$job}{state} = $_->{value};
#			$_->{name} eq "$PBS::ATTR_used.walltime"
#			  and $Jobs->{$job}{elpt} = $_->{value};
#		}
#		if ( exists $Jobs->{$job}{state} and $Jobs->{$job}{state} eq "R" ) {
#			$State_count->{"_rjobs"}++;
#		}
#		$State_count->{"_njobs"}++;
#	}
#}

sub get_info_cmdline {

	my ($server, $Nodes, $Jobs, $State_count) = @_;

	#
	# Use qmgr and qstat commands to get all data.
	#
	my @qmgr = `$qmgr -c 'l n \@$server' $server 2>/dev/null`;
	$? and do { printwarning("Connection to $server failed.") };

	my @qstat_f = `$qstat -f \@$server 2>/dev/null`;
	$? and do { printwarning("Connection to $server failed.") };

	my $jobs;
	my $eatingjobs = 0;
	my $status;
	my $statuses;
	my $eatingstatus = 0;
	my $node         = "";

	#
	# Parse qmgr output.
	#
	foreach (@qmgr) {
		chomp;

		if (/^Node /) {
			$node = $';
			$node =~ s/^targetgcc//;
			$node =~ s/-mgmt$//;
			$node =~ s/\s//g;
			delete $Nodes->{$server}{$node};
			$Nodes->{$server}{$node}{seen} = 1;
			$State_count->{_nodes}++;
			$eatingjobs   = 0;
			$eatingstatus = 0;
		} elsif (/\s+np = (.*)/) {
			$Nodes->{$server}{$node}{np} = $1;
			$State_count->{_procs} += $1;
			$State_count->{_mprocs} =
			    $State_count->{_mprocs} > $1
			  ? $State_count->{_mprocs}
			  : $1;
			$eatingstatus = 0;
			$eatingjobs   = 0;
		} elsif (/\s+properties = (.*)/) {
			$Nodes->{$server}{$node}{properties} = $1;
			$eatingstatus                        = 0;
			$eatingjobs                          = 0;
		} elsif (/\s+ntype = (.*)/) {
			$Nodes->{$server}{$node}{ntype} = $1;
			$eatingstatus                   = 0;
			$eatingjobs                     = 0;
		} elsif (/\s+state = (.*)/) {
			$Nodes->{$server}{$node}{state} = $1;
			$State_count->{$1}++;
			$eatingstatus = 0;
			$eatingjobs   = 0;
		} elsif (/\s+jobs = (.*)/) {
			$eatingjobs   = 1;
			$eatingstatus = 0;
			$jobs         = $1;
			$State_count->{"_anodes"}++;

			#
			# Reset list of jobs on this node.
			#
			$Nodes->{$server}{$node}{jobs} = {};
			foreach my $job_slot (split(/, /, $jobs)) {
				if ($job_slot =~ m{(\d+)/(\d+)}) {
					my $slot   = $1;
					my $job_id = $2;
					$Nodes->{$server}{$node}{job}{$slot} = $job_id;
					$State_count->{"_aprocs"}++;
					if (defined(${$Nodes}{$server}{$node}{jobs}{$job_id})) {
						${$Nodes}{$server}{$node}{jobs}{$job_id}++;
					} else {
						${$Nodes}{$server}{$node}{jobs}{$job_id} = 1;
					}
				} else {
					printwarning("ERROR: Cannot parse job ID from job slot $job_slot.");
				}
			}
		} elsif (/\s+status = (.*)/) {
			$eatingstatus = 1;
			$eatingjobs   = 0;
			$statuses     = $1;
			foreach my $status (split(/,/, $statuses)) {
				if ($status =~ m{(.+)=(.+)}) {
					$Nodes->{$server}{$node}{status}{$1} = $2;
				}
			}
		} elsif ($eatingjobs) {
			if ($_ =~ /\w/) {
				/^\s+(.*)$/;
				$jobs = $1;
				foreach my $job_slot (split(/, /, $jobs)) {
					if ($job_slot =~ m{(\d+)/(\d+)}) {
						my $slot   = $1;
						my $job_id = $2;
						$Nodes->{$server}{$node}{job}{$slot} = $job_id;
						$State_count->{"_aprocs"}++;
						if (defined(${$Nodes}{$server}{$node}{jobs}{$job_id})) {
							${$Nodes}{$server}{$node}{jobs}{$job_id}++;
						} else {
							${$Nodes}{$server}{$node}{jobs}{$job_id} = 1;
						}
					} else {
						printwarning("ERROR: Cannot parse job ID from job slot $job_slot.");
					}
				}
			} else {
				$eatingjobs = 0;
			}
		} elsif ($eatingstatus) {
			if ($_ =~ /\w/) {
				/^\s+(.*)$/;
				$statuses = $1;
				foreach my $status (split(/,/, $statuses)) {
					if ($status =~ m{(.+)=(.+)}) {
						$Nodes->{$server}{$node}{status}{$1} = $2;
					}
				}
			} else {
				$eatingstatus = 0;
			}
		}
	}

	#
	# Parse qstat output.
	#
	if (scalar(@qstat_f) > 0) {
		my $job;
		foreach (@qstat_f) {
			chomp;
			if (/^\s*Job Id: ([0-9]+)\.([a-z0-9A-Z-._]+)/) {
				$job = $1;
				my $scheduler = $2;
				if ($server eq $scheduler) {
					$Jobs->{$job}{seen}   = 1;
					$Jobs->{$job}{server} = $server;
					$State_count->{"_njobs"}++;
				} else {
					die "FATAL: server reported in qstat output ($scheduler) does not match server for which info was requested ($server).";
				}
			} elsif (/^\s*Job_Name = (.+)/i) {
				$Jobs->{$job}{jname} = $1;
			} elsif (/^\s*Job_Owner = ([^@]+)/i) {
				$Jobs->{$job}{user} = $1;
			} elsif (/^\s*job_state = ([CEHQRSTW])/i) {
				my $state = $1;
				$Jobs->{$job}{state} = $state;
				if ($state eq 'R') {
					$State_count->{"_rjobs"}++;
				}
			} elsif (/^\s*queue = (.+)/i) {
				$Jobs->{$job}{queue} = $1;
			} elsif (/^\s*Resource_List.nodes = ([0-9]+):ppn=([0-9]+)/i) {
				$Jobs->{$job}{ncount}  = $1;
				$Jobs->{$job}{cores_r} = $2;
			} elsif (/^\s*Resource_List.mem = ([0-9]+)([a-z]+)$/i) {
				my $mem_value = $1;
				my $mem_unit  = $2;
				if ($mem_unit eq 'b') {
					$mem_value = $mem_value / 1073741824;
				} elsif ($mem_unit eq 'kb') {
					$mem_value = $mem_value / 1048576;
				} elsif ($mem_unit eq 'mb') {
					$mem_value = $mem_value / 1024;
				} elsif ($mem_unit eq 'gb') {
					$mem_value = $mem_value;
				} else {
					die "FATAL: Unsupported mem unit $mem_unit in $qstat -f output.";
				}
				$Jobs->{$job}{mem_r} = $mem_value;
			} elsif (/^\s*Resource_List.walltime = ([0-9:]+)/i) {
				$Jobs->{$job}{walltime_r} = $1;
			} elsif (/^\s*resources_used.cput = (.+)$/i) {
				$Jobs->{$job}{cputime_u} = $1;
			} elsif (/^\s*resources_used.mem = ([0-9]+)([a-z]+)$/i) {
				my $mem_value = $1;
				my $mem_unit  = $2;
				if ($mem_unit eq 'b') {
					$mem_value = $mem_value / 1073741824;
				} elsif ($mem_unit eq 'kb') {
					$mem_value = $mem_value / 1048576;
				} elsif ($mem_unit eq 'mb') {
					$mem_value = $mem_value / 1024;
				} elsif ($mem_unit eq 'gb') {
					$mem_value = $mem_value;
				} else {
					die "FATAL: Unsupported mem unit $mem_unit in $qstat -f output.";
				}
				$Jobs->{$job}{mem_u} = $mem_value;
			} elsif (/^\s*resources_used.walltime = (.+)/i) {
				$Jobs->{$job}{walltime_u} = $1;
			}
		}
	}

	#
	# Revisit info about all nodes and copy some details to the data structure indexed by job.
	#
	foreach my $server (keys(%{$Nodes})) {
		foreach my $node (keys(%{ ${$Nodes}{$server} })) {
			foreach my $job (keys(%{ ${$Nodes}{$server}{$node}{jobs} })) {
				if (defined(${$Jobs}{$job}) && ${$Jobs}{$job}{server} eq $server) {
					$Jobs->{$job}{node}    = $node;
					$Jobs->{$job}{cores_u} = ${$Nodes}{$server}{$node}{jobs}{$job};
				} else {
					printwarning("ERROR: Cannot find job ID ${job} in list of jobs or the scheduler does not match (nodes server = $server | jobs server = ${$Jobs}{$job}{server}).");
				}
			}
		}
	}

}

sub update_display {

	my ($state_count, $mprocs, $maxlen, $nodes, $jobs, $users, $queues) = @_;

	my $foo;
	move($pad, 0, 0);
	getmaxyx($Y, $X);

	$y = 0, $x = 0;

	$show_summary and show_state_summary($state_count);
	if ($show_grid) {
		show_grid($jobs, $nodes, $users, $queues, $maxlen, $mprocs);
		my $legend = "[$state_characters{'unknown'}] unknown";
		$legend .= "  [$state_characters{'busy'}] busy";
		$legend .= "  [$state_characters{'down'}] down";
		$legend .= "  [$state_characters{'idle'}] idle";
		$legend .= "  [$state_characters{'offline'}] offline";
		$legend .= "  [$state_characters{'other'}] other";
		addstr($pad, $y, 0, $legend);
		clrtoeol($pad);
		move($pad, ++$y, $x = 0);
	}
	$show_queue and show_queue($jobs, $nodes, $users, $queues);
	
	getyx($pad, $ly, $foo);
	clrtobot($pad);
	
	pnoutrefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
	mvwin($cmdwin, $Y - 1, 0);
	refresh($cmdwin);
	update_subwin(@_);
	doupdate();
}

sub show_state_summary {
	my $t  = 1;
	my $t2 = 1;

	addstr $pad, sprintf("Usage Totals: %d/%d %s, %d/%d %s, %d/%d %s", ${ $_[0] }{_aprocs}, ${ $_[0] }{_procs}, "Procs", ${ $_[0] }{_anodes}, ${ $_[0] }{_nodes}, "Nodes", ${ $_[0] }{_rjobs}, ${ $_[0] }{_njobs}, "Jobs Running");

	my ($y1, $x1);
	getyx($pad, $y1, $x1);
	addstr $pad, " " x ($X - $x1 - 8);

	# Asbed asked for the time
	addstr $pad, sprintf("%02d:%02d:%02d", (localtime())[ 2, 1, 0 ]);

	my $line;
	my @states = sort grep !/^_/, keys %{ $_[0] };

	if ($compact_summary) {
		move($pad, ++$y, $x = 0);
		addstr($pad, 1, 0, "Node States:");
		for (my $i = 0 ; defined $states[$i] ; $i++) {
			$line = " " . ${ $_[0] }{ $states[$i] } . " " . $states[$i];
			$line .= defined $states[ $i + 1 ] ? "," : "";
			getyx($pad, $y1, $x1);
			if ($X - $x1 - 1 < length($line)) {
				clrtoeol($pad);
				move($pad, ++$y, $x = 0);
				addstr($pad, " " x 12);
			}
			addstr($pad, $line);
		}
		move($pad, ++$y, $x = 0);
		clrtoeol($pad);

	} else {

		for (my $i = 0 ; defined $states[$i] ; $i++) {
			move($pad, ++$y, $x = 0);
			$line = " " x 14;
			$line .= sprintf("%4s %-20s", ${ $_[0] }{ $states[$i] }, $states[$i]);

			$i++;
			$line .= sprintf("%4s %-20s", ${ $_[0] }{ $states[$i] }, $states[$i])
			  if defined $states[$i];
			$line .= " " x ($X - (length($line) + 14));
			clrtoeol($pad);
			addstr $pad, $line;
		}
		addstr($pad, 1, 0, "Node States:");
		move($pad, ++$y, $x = 0);
	}
}

sub show_grid {
	my ($jobs, $allnodes, $users, $queues, $maxlen, $maxprocs) = @_;
	my ($foo, $tmpx);
	$lx = 0;

	clrtoeol($pad);
	move($pad, ++$y, $x = 0);

	if (!scalar @show_cpu) {

		#addstr $pad, "  No CPUs selected!";
		#clrtoeol($pad);
		#return;
		$show_cpu[0] = "0";
	}

	#printvcpuline($maxlen);

	foreach my $server (keys %$allnodes) {

		my (@cluster, @ts);
		foreach my $node (sort keys %{ $allnodes->{$server} }) {
			$$allnodes{$server}{$node}{np} > $maxnodegrid
			  ? push(@ts,      $node)
			  : push(@cluster, $node);
		}

		# loop through each node, in lines and columns
		my $col = 0;

		#		if ( defined $cluster[0] ) {
		#			my $headerspaces = scalar @show_cpu;
		#			printnumberline( $maxlen, $headerspaces, $columns );
		#			printdashline( $maxlen, $headerspaces, $columns );
		#
		#			foreach my $node (@cluster) {
		#				$col = 0 if $col >= $columns;
		#				( addstr $pad, sprintf "  %${maxlen}s ", $node ) if $col == 0;
		#
		#				addstr $pad, "  " if ( $col != 0 and $col % 10 == 0 );
		#
		#				foreach my $this_cpu (@show_cpu) {
		#					if ( $this_cpu > $$allnodes{$server}{$node}{np} - 1 ) {
		#						addstr $pad, " ";
		#
		#					}
		#					else {
		#						my $state = $$allnodes{$server}{$node}{state};
		#						if ( exists $$allnodes{$server}{$node}{job}{$this_cpu} )
		#						{
		#							my $job =
		#							  $$allnodes{$server}{$node}{job}{$this_cpu};
		#							if (
		#								    exists $$jobs{$job}
		#								and exists $$jobs{$job}{letter}
		#								and ( $show_user
		#									? $show_user =~ /\b$$jobs{$job}{user}\b/
		#									: 1 )
		#							  )
		#							{
		#								my $letter    = $$jobs{$job}{letter};
		#								my $color     = $$jobs{$job}{color};
		#								my $underline = $jobs->{$job}{underline};
		#								printcpustate( $job, $letter, $state, $color,
		#									$underline );
		#							}
		#							else {
		#
		#						   # The job was deleted in between the time we captured
		#						   # the node info and the job info from pbs.
		#								printcpustate( 0, 0, $state, 0, 0 );
		#							}
		#						}
		#						else {
		#							printcpustate( 0, 0, $state, 0, 0 );
		#						}
		#
		#					}
		#				}
		#
		#				addstr $pad, " ";
		#				$col++;
		#				getyx( $pad, $foo, $tmpx );
		#				$lx = $lx > $tmpx ? $lx : $tmpx;
		#
		#				clrtoeol($pad);
		#				move( $pad, ++$y, $x = 0 ) if $col >= $columns;
		#
		#				#return if $y >= $Y;
		#			}
		#			clrtoeol($pad);
		#			move( $pad, ++$y, $x = 0 ) if $col != $columns;
		#			printdashline( $maxlen, $headerspaces, $columns );
		#			clrtoeol($pad);
		#			move( $pad, ++$y, $x = 0 );
		#		}

		my $headerspaces = 1;
		printnumberline($maxlen, $headerspaces, $columns);

		foreach my $node (@ts) {

			unless ($show_node ? $show_node =~ /\b$node\b/ : 1) {
				next;
			}

			printdashline($maxlen, $headerspaces, $columns);
			my $col = 0;
			my $load =
			  defined($$allnodes{$server}{$node}{status}{'loadave'})
			  ? $$allnodes{$server}{$node}{status}{'loadave'}
			  : '?';
			my $node_state      = $$allnodes{$server}{$node}{state};
			my $available_cores = $$allnodes{$server}{$node}{np};
			my $used_cores      = 0;

			if ($sorted_grid == 1) {

				#
				# Grid display sorted by job.
				#
				foreach my $job (sort(keys(%{ ${$allnodes}{$server}{$node}{jobs} }))) {

					my $letter    = 0;
					my $color     = 0;
					my $underline = 0;
					my $job_id    = 0;
					my $cpu_state = 'busy';

					if (
						    exists $$jobs{$job}
						and exists $$jobs{$job}{letter}
						and exists $$jobs{$job}{node}
						and (
							$show_user ? $show_user =~ /\b$$jobs{$job}{user}\b/
							: 1
						)
						and (
							$show_node ? $show_node =~ /\b$$jobs{$job}{node}\b/
							: 1
						)
					  )
					{
						$job_id    = $job;
						$letter    = $jobs->{$job}{letter};
						$color     = $jobs->{$job}{color};       # By default color by job.
						$underline = $jobs->{$job}{underline};
						my $user  = $jobs->{$job}{user};
						my $queue = $jobs->{$job}{queue};

						if ($color_by eq 'user') {
							$color = $users->{$user}{color};
						} elsif ($color_by eq 'queue') {
							$color = $queues->{$queue}{color};
						}
					}

					for (my $slots = 1 ; $slots <= ${$allnodes}{$server}{$node}{jobs}{$job} ; $slots++) {
						if ($col >= $columns) {
							move($pad, ++$y, $x = 0);
							$col = 0;
						}
						if ($col == 0) {
							addstr($pad, sprintf "  %${maxlen}s ", $node);
						} elsif ($col != 0 and $col % 10 == 0) {
							addstr $pad, "  ";
						}
						printcpustate($job_id, $letter, $cpu_state, $color, $underline);
						addstr $pad, " ";
						getyx($pad, $foo, $tmpx);
						$lx = $lx > $tmpx ? $lx : $tmpx;
						clrtoeol($pad);
						$col++;
						$used_cores++;
					}
				}

				while ($used_cores < $$allnodes{$server}{$node}{np}) {
					if ($col >= $columns) {
						move($pad, ++$y, $x = 0);
						$col = 0;
					}
					if ($col == 0) {
						addstr($pad, sprintf "  %${maxlen}s ", $node);
					} elsif ($col != 0 and $col % 10 == 0) {
						addstr $pad, "  ";
					}
					printcpustate(0, 0, $node_state, 0, 0);
					addstr $pad, " ";
					$col++;
					$used_cores++;
					getyx($pad, $foo, $tmpx);
					$lx = $lx > $tmpx ? $lx : $tmpx;
					clrtoeol($pad);
				}

				while ($col < $columns) {
					addstr $pad, "  " if ($col != 0 and $col % 10 == 0);
					addstr $pad, "  ";
					$col++;
					getyx($pad, $foo, $tmpx);
					$lx = $lx > $tmpx ? $lx : $tmpx;
					clrtoeol($pad);
				}
				addstr($pad, sprintf "%6s ", $load);
				move($pad, ++$y, $x = 0);

			} else {

				#
				# Grid display unsorted (disply in CPU/core order).
				#
				foreach my $this_cpu (0 .. $$allnodes{$server}{$node}{np} - 1) {

					$col = 0 if $col >= $columns;

					if ($col == 0) {
						addstr($pad, sprintf "  %${maxlen}s ", $node);
					}

					addstr $pad, "  " if ($col != 0 and $col % 10 == 0);

					#my $state = $$allnodes{$server}{$node}{state};

					if (exists $$allnodes{$server}{$node}{job}{$this_cpu}) {

						my $job = $$allnodes{$server}{$node}{job}{$this_cpu};
						if (
							    exists $$jobs{$job}
							and exists $$jobs{$job}{letter}
							and exists $$jobs{$job}{node}
							and (
								$show_user ? $show_user =~ /\b$$jobs{$job}{user}\b/
								: 1
							)
							and (
								$show_node ? $show_node =~ /\b$$jobs{$job}{node}\b/
								: 1
							)
						  )
						{
							my $letter    = $jobs->{$job}{letter};
							my $color     = $jobs->{$job}{color};
							my $underline = $jobs->{$job}{underline};
							my $cpu_state = 'busy';
							my $user      = $jobs->{$job}{user};
							my $queue     = $jobs->{$job}{queue};

							if ($color_by eq 'user') {
								$color = $users->{$user}{color};
							} elsif ($color_by eq 'queue') {
								$color = $queues->{$queue}{color};
							}

							printcpustate($job, $letter, $cpu_state, $color, $underline);
						} else {
							printcpustate(0, 0, $node_state, 0, 0);
						}
					} else {
						printcpustate(0, 0, $node_state, 0, 0);
					}

					addstr $pad, " ";
					$col++;
					getyx($pad, $foo, $tmpx);
					$lx = $lx > $tmpx ? $lx : $tmpx;

					clrtoeol($pad);

					if ($col >= $columns) {
						addstr($pad, sprintf "%6s ", $load);
						move($pad, ++$y, $x = 0);
					}
				}
			}

			clrtoeol($pad);
			move($pad, ++$y, $x = 0) if $col != $columns;
		}

		printdashline($maxlen, 1, $columns);
		clrtoeol($pad);
		move($pad, ++$y, $x = 0);

	}

	clrtoeol($pad);

}

# Print out the job queue
sub show_queue {
	my ($jobs, $allnodes, $users, $queues) = @_;
	my $format = "%-6s %-15.15s  %-12.12s  %-40.40s  %1s  %7.7s %7.7s  %8.8s %8.8s  %9.9s %9.9s";
	move($pad, ++$y, $x = 0);

	#return if $y >= $Y;
	attron($pad, A_BOLD);
	addstr($pad, '      ');
	addstr($pad, sprintf($format, 'JobID', 'Username', 'Queue', 'Jobname', 'S', '~CPU(%)', 'CPU(%)', 'Mem(GiB)', 'Mem(GiB)', 'Walltime', 'Walltime'));
	clrtoeol($pad);
	move($pad, ++$y, $x = 0);
	addstr($pad, '      ');
	addstr($pad, sprintf($format, '', '', '', '', '', 'Used', 'Req', 'Used', 'Req', 'Used', 'Req'));
	clrtoeol($pad);
	attroff($pad, A_BOLD);
	move($pad, ++$y, $x = 0);

	# Note: we never print our the server name, it is expected that the user will
	# recognize jobs by their jobid or queue name.

	foreach my $job (

		# Sort first by server, then by queue, last by jobid
		sort { defined $jobs->{$a}{server} && defined $jobs->{$b}{server} && $jobs->{$a}{server} cmp $jobs->{$b}{server} or defined $jobs->{$a}{queue} && defined $jobs->{$b}{queue} && $jobs->{$a}{queue} cmp $jobs->{$b}{queue} or $a <=> $b } keys %{$jobs}
	  )
	{

		my $l         = $jobs->{$job}{letter};
		my $color     = $jobs->{$job}{color};
		my $underline = $jobs->{$job}{underline};
		my $user      = $jobs->{$job}{user};
		my $queue     = $jobs->{$job}{queue};
		my $node      = 'NAmissingUndefined';
		my $cores_u;
		my $mem_u;

		if (defined($jobs->{$job}{node})) {
			$node = $jobs->{$job}{node};
		}
		if ($color_by eq 'user') {
			$color = $users->{$user}{color};
		} elsif ($color_by eq 'queue') {
			$color = $queues->{$queue}{color};
		}

		next
		  unless $show_user ? $show_user =~ /\b$jobs->{$job}{user}\b/ : 1;

		next
		  unless $show_node ? $show_node =~ /\b$node\b/ : 1;

		if (!($jobs->{$job}{state} =~ /Q|H/ and !$show_qqueue)) {
			addstr $pad, "  ";
			if (defined $l and $jobs->{$job}{state} eq "R") {
				print_colored_letter($l, $color, $underline);
				addstr $pad, " = ";
			} else {
				addstr $pad, "    ";
			}

			# The job is erroring, grab an eyeball
			if ($jobs->{$job}{state} eq "E") {
				attron($pad, A_REVERSE);
			}

			if (defined($jobs->{$job}{cputime_u}) && defined($jobs->{$job}{walltime_u})) {
				my @cputime = reverse(split(':', $jobs->{$job}{cputime_u}));
				my $cputime_sec = $cputime[0] + ($cputime[1] * 60) + ($cputime[2] * 3600);
				my @walltime = reverse(split(':', $jobs->{$job}{walltime_u}));
				my $walltime_sec = $walltime[0] + ($walltime[1] * 60) + ($walltime[2] * 3600);
				if ($walltime_sec > 0) {
					$cores_u = sprintf("%u", ($cputime_sec / $walltime_sec) * 100);
				}
				#printwarning("DEBUG: walltime_sec = $walltime[0] | $walltime[1] | $walltime[2]");
				#printwarning("DEBUG: cputime = $jobs->{$job}{cputime_u} | walltime = $jobs->{$job}{walltime_u} | cputime_sec = $cputime_sec | walltime_sec = $walltime_sec");
			}

			if (exists($jobs->{$job}{mem_u})) {
				$mem_u = sprintf("%.1f", $jobs->{$job}{mem_u});
			}

			addstr $pad, sprintf($format, $job, $jobs->{$job}{user}, $jobs->{$job}{queue}, $jobs->{$job}{jname}, $jobs->{$job}{state}, defined $cores_u ? $cores_u : "", exists $jobs->{$job}{cores_r} ? $jobs->{$job}{cores_r} * 100 : "", defined $mem_u ? $mem_u : "", exists $jobs->{$job}{mem_r} ? $jobs->{$job}{mem_r} : "", exists $jobs->{$job}{walltime_u} ? $jobs->{$job}{walltime_u} : "", exists $jobs->{$job}{walltime_r} ? $jobs->{$job}{walltime_r} : "");

			if ($jobs->{$job}{state} eq "E") {
				attroff($pad, A_REVERSE);
			}
			clrtoeol($pad);
			move($pad, ++$y, $x = 0);
		}
	}

	clrtoeol($pad);
}

# Pass in a reference to all jobs, and we'll assign a letter to each one.
sub letterize {
	my $Jobs = shift;

	# The original qtop only used one letter per job, and this info was held in %Job_of_letter.
	# Now that we've thrown out that limitation, it is only used to note the fact that _someone_ is using that letter.
	# If a job gets a letter that is already assigned, the second one will be noted in %Job_of_letter,
	# but that's ok because we don't care _who_ has that letter.

	# remove info about old jobs and jobs already assigned a letter
	foreach my $l (keys %Job_of_letter) {
		delete $Job_of_letter{$l}
		  if (exists $Jobs->{ $Job_of_letter{$l} }{user});
	}

	# pick a letter if not already choosen
	foreach my $job (keys %{$Jobs}) {
		next if !defined $Jobs->{$job}{state};
		next if $Jobs->{$job}{state} eq "Q";
		my $user = $Jobs->{$job}{user};
		if (!exists $Jobs->{$job}{letter}) {

			# find a letter that isn't already taken
			my $l = substr($user, 0, 1);
			if (exists $Job_of_letter{$l}) {
				$l = uc($l);
				if (exists $Job_of_letter{$l}) {
					if (length $letters <= 0) {

						# replenish our supply of letters
						$colorize
						  or printwarning("Reusing letters on B&W terminal.");
						$letters   = $masterletters;
						$underline = !$underline;
					}
					$letters =~ s/(.)//;
					$l = $1;
				}
			}
			$Job_of_letter{$l}       = $job;
			$Jobs->{$job}{letter}    = $l;
			$Jobs->{$job}{underline} = $underline;
			if ($l =~ m/[^a-zA-Z0-9_]/) {
				$l = '\\' . "$l";
			}
			$letters =~ s/$l//;
		}
	}
}

#
# Assign a color to each job, user and queue.
#
sub colorize {
	my ($Jobs, $Users, $Queues) = @_ or return;

	foreach my $job (keys %{$Jobs}) {
		next if defined $Jobs->{$job}{color};

		scalar @Colors == 0 and @Colors = init_colors();
		my $color = shift @Colors;
		$Jobs->{$job}{color} = $color;

		my $user  = $Jobs->{$job}{user};
		my $queue = $Jobs->{$job}{queue};
		$Users->{$user}{color}   = $color unless defined($Users->{$user}{color});
		$Queues->{$queue}{color} = $color unless defined($Queues->{$queue}{color});
	}
}

# This sucks, I wanted to seperate printing from colors,
# but I can't just pass back color escape strings.
# I'm forced to combine them here.
sub print_colored_letter {
	my ($letter, $color, $underline) = @_;
	$colorize or do { addstr $pad, $letter; return };

	attron($pad, A_BOLD | COLOR_PAIR($color) | ($underline && A_UNDERLINE));

	#printwarning("DEBUG: color = $color");
	addstr $pad, $letter;
	attroff($pad, A_BOLD | COLOR_PAIR($color) | ($underline && A_UNDERLINE));
}

# Used to find the longest hostname to align the left side of the node grid
sub getmaxkeylen {
	my $maxlen = 0;
	my $server;
	foreach $server (keys(%{ $_[0] })) {
		foreach (keys(%{ ${ $_[0] }{$server} })) {
			$maxlen =
			    length($_) > $maxlen
			  ? length($_)
			  : $maxlen;
		}
	}
	return $maxlen;
}

#
# Prints each character representing each core in the node grid.
#
sub printcpustate {

	my ($job, $letter, $state, $color, $underline) = @_;

	if ($job) {
		printwarning("WARN: $job has no letter.") unless $letter;
		$letter = '&' unless $letter;
	}

	#
	#FIXME# I'm ignoring job-sharing here because I've never seen it used.
	#FIXME# It would be more correct to use the constants defined in pbs_ifl.h,
	#FIXME# but that would break compatibility with non-perl-PBS environments.
	#
	if ($state =~ /down/ and $job) {
		print_colored_letter($letter, $color, $underline);
	} elsif ($job and $show_jobs) {
		print_colored_letter($letter, $color, $underline);
	} elsif ($state =~ /down/) {
		#print_colored_letter($state_characters{'down'}, 1, 0);    # Color pair 1 is red on black
		addch($pad, $state_characters{'down'});
	} elsif ($state =~ /offline/) {
		addch($pad, $state_characters{'offline'});
	} elsif ($state =~ /job-exclusive/) {
		addch($pad, $state_characters{'busy'});
	} elsif ($state =~ /busy/) {
		addch($pad, $state_characters{'busy'});
	} elsif ($state =~ /reserve/) {
		addch($pad, $state_characters{'busy'});
	} elsif ($state =~ /unknown/i) {
		addch($pad, $state_characters{'unknown'});
	} elsif ($state =~ /free/) {
		addch($pad, $state_characters{'idle'});
	} else {
		addch($pad, $state_characters{'other'});
	}
}

# Print the list of visible CPUs above the node grid.
sub printvcpuline {
	my $maxlen = shift;

	# inform the user of the visible CPUs.
	if (scalar @show_cpu == 1) {
		addstr $pad, "  CPU $show_cpu[0]" . " " x ($maxlen - 4);
	}

	#elsif ( scalar @show_cpu == $maxprocs ) {
	#addstr "   " . " " x ($maxlen);
	#}
	else {
		addstr $pad, " " x ($X - 1);
		move($pad, $y, 0);
		addstr $pad, " visible CPUs: " . join(",", @show_cpu);
	}
	clrtoeol($pad);
	move($pad, ++$y, $x = 0);
}

# Print the line of dashes above and below the node grid.
sub printdashline {
	my $maxlen  = shift;
	my $spaces  = shift;
	my $columns = shift;

	my $line = "   " . " " x $maxlen;
	for (my $i = 0 ; $i < $columns ; $i++) {
		$line .= "--" if ($i != 0 and $i % 10 == 0);
		$line .= "-" . "-" x $spaces;
	}
	$line =~ s/-$//;    # oops, we printed one extra, erase it
	addstr $pad, $line;
	clrtoeol($pad);
	move($pad, ++$y, $x = 0);
}

# Print the repetitive line of numbers along the top of the node grid.
sub printnumberline {
	my $maxlen  = shift;
	my $spaces  = shift;
	my $columns = shift;

	my $line = '  ' . 'node' . ' ' x ($maxlen - 4) . ' ';
	for (my $i = 0, my $j = 0 ; $i < $columns ; $i++, $j++) {
		if ($i != 0 and $i % 10 == 0) {
			$line .= "  ";
			$j = 0;
		}
		$line .= (($j + 1) % 10) . " " x $spaces;
	}
	addstr($pad, $line);
	addstr($pad, sprintf "%6s ", 'load');
	clrtoeol($pad);
	move($pad, ++$y, $x = 0);

}

# This is used in top_sleep to annoy the user.
sub printwarning {
	attron($cmdwin, A_REVERSE);
	addstr $cmdwin, 0, 0, join(" ", @_);
	attroff($cmdwin, A_REVERSE);
	refresh($cmdwin);
}

# This is used in top_sleep to solicit the user.
sub getstring {
	my $input = "";
	my $ch;
	my $x = 0;
	addstr $cmdwin, 0, 0, join(" ", @_);
	$x = length join(" ", @_);
	clrtoeol($cmdwin);
	refresh($cmdwin);
	echo;
	nodelay($cmdwin, 0);

	#getstr( $cmdwin,  $input );
	$ch = getch($cmdwin);
	while (1) {
		$ch = getch($cmdwin);
		if ($ch eq ERR) {    # ERR returned on timeout
			next;

			# why is this so freaking complicated??
		} elsif ($ch eq KEY_BACKSPACE or $ch eq KEY_DC or $ch eq $CTRL_H) {
			if (length $input) {
				$x--;
				move($cmdwin, 0, $x);
				delch($cmdwin);
				refresh($cmdwin);
				$input =~ s/.$//;
			} elsif ($ch eq KEY_BACKSPACE) {
				move($cmdwin, 0, $x);
			}
		} elsif ($ch eq $CTRL_G) {    # user abort
			$input = "";
			last;
		} elsif ($ch eq "\n") {
			last;
		} else {
			$x++;
			$input .= $ch;
		}
	}
	noecho;
	move($cmdwin, 0, 0);
	clrtoeol($cmdwin);
	refresh($cmdwin);
	return $input;
}

sub print_serverstatus_window {
	my $server = shift;
	destroy_subwin();
	printwarning("TODO: This information is not available without perl-PBS");
	return;

	my $con = pbs_connect($server);
	if ($con <= 0) {
		destroy_subwin();
		printwarning("Connect to $server failed: $PBS::pbs_errno\n");
		return;
	}
	my $ref = pbs_statserver($con, undef, undef);
	pbs_disconnect($con);
	print_status_window("$server", $ref->[0]->{attribs});
}

# since we don't store enough info about jobs in %Jobs, go ahead and get
# it from the server
sub print_jobstatus_window {
	my $job    = shift;
	my $server = shift;
	destroy_subwin();
	printwarning("TODO: This information is not available without perl-PBS");
	return;

	my $con = pbs_connect($server);
	if ($con <= 0) {
		destroy_subwin();
		printwarning("Connect to $server failed: $PBS::pbs_errno\n");
		return;
	}

	my $ref = pbs_statjob($con, "$job.$server", undef, undef);
	pbs_disconnect($con);

	print_status_window("$job.$server", $ref->[0]->{attribs}, "'l' for node load report");

}

sub print_jobloadstatus_window {
	my $job      = shift;
	my $server   = shift;
	my $allnodes = shift;

	my @loads;
	my $value;
	my $freephys;
	my $sessions;
	foreach $server (sort keys %{$allnodes}) {
	  NODE:
		foreach my $node (sort keys %{ $allnodes->{$server} }) {
			foreach my $this_cpu (0 .. $allnodes->{$server}{$node}{np} - 1) {
				if (    exists $allnodes->{$server}{$node}{job}{$this_cpu}
					and exists $allnodes->{$server}{$node}{status}{loadave})
				{
					if ($job eq $allnodes->{$server}{$node}{job}{$this_cpu}) {
						$value = "load: " . $allnodes->{$server}{$node}{status}{loadave};
						{    # yes, I'm TOTALLY cheating here.  Sue me.
							local $^W = 0;
							$freephys = int(($allnodes->{$server}{$node}{status}{physmem} - ($allnodes->{$server}{$node}{status}{totmem} - $allnodes->{$server}{$node}{status}{availmem})) / 1024);
							$sessions =
							  ($allnodes->{$server}{$node}{status}{nsessions} =~ /^\?/)
							  ? 0
							  : $allnodes->{$server}{$node}{status}{nsessions};

							#: scalar (split / /,$allnodes->{$server}{$node}{status}{sessions});
							$value .= "  physmem: " . int($allnodes->{$server}{$node}{status}{physmem} / 1024) . "MB";
							$value .= " avail: ${freephys}MB";
							$value .= "  sessions: $sessions";

						}
						push(
							@loads,
							{
								name  => $node,
								value => $value
							}
						);
						next NODE;
					}
				}
			}
		}
	}

	if (scalar @loads < 1) {
		push(
			@loads,
			{
				name  => "load report",
				value => "The job has ended, or the server is too old"
			}
		);
	}
	print_status_window("$job Load Report", \@loads, "'l' for job details");

}

# We already have everything we need to know about nodes in our big Nodes struct,
# so just pull info from there.
sub print_nodestatus_window {
	my $nodename = shift or return;
	my $ref      = shift or return;

	my (@attrs, $name, $value);
	while (($name, $value) = each %{ $ref->{status} }) {
		push(@attrs, { name => $name, value => $value });
	}
	push(@attrs, { name => "ntype", value => $ref->{ntype} });
	push(@attrs, { name => "state", value => $ref->{state} });
	if (exists $ref->{properties}) {
		push(@attrs, { name => "properties", value => $ref->{properties} });
	}

	#FIXME# only offer the 'j' option when only one job running, and not just on CPU0
	my $multiplejob = 0;
	foreach my $cpu (
		sort { $multiplejob ||= $ref->{job}{$a} != $ref->{job}{$b}; $a <=> $b }
		keys %{ $ref->{job} }
	  )
	{
		push(@attrs, { name => "CPU$cpu: job#", value => $ref->{job}{$cpu} });
	}

	#push(@attrs, { name => "multiple", value => $multiplejob });
	print_status_window("$nodename", \@attrs, "'j' for job details on CPU0");

}

#
# This is used by the subs above to actually paint the subpad.
#
sub print_status_window {
	my $title    = shift;
	my $ref      = shift;
	my $epilogue = shift;

	my $line = 1;
	my $string;
	my $maxlinelen;
	my $pat;
	my $indent;

	my $name;
	my $value;

	#
	# Subwin's width is 8 fewer than the main win, and with 2 chars padding inside,
	# each line will be 10 chars less wide than the main window.
	# Select the smallest sizes from default and size currently in use
	# to make sure the subpad is smaller than its parent pad.
	# When the dimensions for the subpad are larger than those of the parent,
	# subpad creating will fail!
	#
	my $subpad_width;
	my $subpad_heigth;
	if ($maxcolumns < $X) {
		$subpad_width = $maxcolumns - 8;
	} else {
		$subpad_width = $X - 8;
	}
	if ($maxrows < $Y) {
		$subpad_heigth = $maxrows - 5;
	} else {
		$subpad_heigth = $Y - 5;
	}

	#
	#FIXME# do we need to properly destroy the subpads?
	#
	$subwin = subpad($pad, $subpad_heigth, $subpad_width, 2, 5);

	$subwin or die "FATAL: Failed to print_status_window. X = $X | Y = $Y | maxcolumns = $maxcolumns | maxrows = $maxrows";
	move($subwin, $line, 0);
	clrtoeol($subwin);
	move($subwin, $line, 2);

	foreach my $attr (@{$ref}) {
		$name       = $attr->{name};
		$value      = $attr->{value};
		$indent     = length($name) + 2 + 3;     # 2 for padding, 3 for " = "
		$maxlinelen = $subpad_width - $indent;
		$pat        = ".{1,$maxlinelen}";

		addstr($subwin, $name . " = ");
		$string = $value;

		if ((length($string) + $indent) > $maxlinelen) {
			while (length($string)) {

				move($subwin, $line, $indent);
				$string =~ s/($pat)// or die "FATAL: Cannot truncate '$string' using regex '$pat'";
				addstr $subwin, "$1";
				clrtoeol($subwin);
				move($subwin, ++$line, 0);
				clrtoeol($subwin);
				move($subwin, $line, 2);

				$indent = 7;

			}
			move($subwin, $line, 0);
			clrtoeol($subwin);
			move($subwin, $line, 2);
		} else {
			addstr $subwin, $string;
			clrtoeol($subwin);
			move($subwin, ++$line, 0);
			clrtoeol($subwin);
			move($subwin, $line, 2);
		}
	}

	move($subwin, ++$line, 0);
	clrtoeol($subwin);
	move($subwin, $line, 2);
	if (defined $epilogue) {
		addstr $subwin, "'q' to exit  $epilogue";
	} else {
		addstr $subwin, "'q' to exit this window";
	}
	move($subwin, ++$line, 0);
	clrtoeol($subwin);

	# make a nice box window border for our output
	resize($subwin, $line + 2, $X - 10);
	attron($subwin, COLOR_PAIR(1));
	box($subwin, &ACS_VLINE, &ACS_HLINE);
	move($subwin, 0, 3);
	addch($subwin, &ACS_RTEE);
	addstr $subwin, "$title";
	addch($subwin, &ACS_LTEE);
	attroff($subwin, COLOR_PAIR(1));

	# currently not used (since the subwin is always smaller than the terminal)
	$subY = $line + 4;
	$subX = $X - 10;
}

sub update_subwin {
	$subwin or return;
	if ($searchobject{TYPE} eq "SERVER") {
		print_serverstatus_window($searchobject{VALUE});
	} elsif ($searchobject{TYPE} eq "NODE") {
		print_nodestatus_window($searchobject{VALUE}, $_[1]->{ $searchobject{SERVER} }{ $searchobject{VALUE} });
	} elsif ($searchobject{TYPE} eq "JOB") {
		print_jobstatus_window($searchobject{VALUE}, $searchobject{SERVER});
	} elsif ($searchobject{TYPE} eq "JOBLOAD") {
		print_jobloadstatus_window($searchobject{VALUE}, $searchobject{SERVER}, $_[1]);
	} else {
		printwarning("oddly, I'm on line " . __LINE__);
	}
	pnoutrefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
}

sub destroy_subwin {
	$subwin       = 0;
	%searchobject = ();
	$subY         = 0;
	$subX         = 0;
	if ($py >= $ly - $Y + 2) {
		$py = $ly - $Y + 3;
		pnoutrefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
	}
}

#FIXME#  top_sleep() is a kludge, I know it... It just keeps growing
#FIXME#  as I add new commands.  *shrug*
sub top_sleep {

	my $targettime = time() + $sleeptime;

	while (time() < $targettime) {
		halfdelay(1);
		my $input = getch($cmdwin);
		if ($SIGWINCH) {
			$SIGWINCH = 0;
			endwin;
			refresh();
			update_display(@_);
		}
		if (defined $input) {
			if ($input eq "q") {
				if ($subwin) {
					destroy_subwin();
					update_display(@_);
				} else {
					endwin;
					exit(0);
				}
			}

			# why doesn't curses do this automatically??
			elsif ($input eq $CTRL_L) {
				clear($pad);
				clrtoeol($cmdwin);
				update_display(@_);
			}

			#FIXME# $helpwin should be a scrollable pad
			elsif ($input eq "h" || $input eq "?") {
				my $helpwin = newwin(0, 0, 0, 0);
				attron($helpwin, A_REVERSE | COLOR_PAIR(6));
				addstr $helpwin, "qtop v$VERSION";
				attroff($helpwin, A_REVERSE | COLOR_PAIR(6));
				move($helpwin, 2, 0);
				addstr $helpwin, "Seconds Refresh ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, "$sleeptime";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nGrid Columns ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, "$columns";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nColorization ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $colorize ? "on" : "off";
				attroff($helpwin, A_BOLD);

				if ($colorize) {
					addstr $helpwin, "\nColor by ";
					attron($helpwin, A_BOLD);
					addstr $helpwin, $color_by;
					attroff($helpwin, A_BOLD);
				}
				addstr $helpwin, "\nState Summary Display ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $show_summary ? "on" : "off";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nGrid Display ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $show_grid ? "on" : "off";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nGrid Job Display ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $show_jobs ? "on" : "off";
				attroff($helpwin, A_BOLD);

				#addstr $helpwin, "\nShow CPU Number ";
				#attron( $helpwin, A_BOLD );
				#addstr $helpwin, join( " ", @show_cpu );
				#attroff( $helpwin, A_BOLD );
				addstr $helpwin, "\nShowing jobs from: ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, join(" ", @host);
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nQueue Display ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $show_queue ? "on" : "off";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nShow Queued Jobs ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $show_qqueue ? "on" : "off";
				attroff($helpwin, A_BOLD);
				addstr $helpwin, "\nNumber of possible colors ";
				attron($helpwin, A_BOLD);
				addstr $helpwin, $COLOR_PAIRS;
				attroff($helpwin, A_BOLD);

				if ($show_user) {
					addstr $helpwin, "\nLimiting job view to user(s) ";
					attron($helpwin, A_BOLD);
					addstr $helpwin, $show_user;
					attroff($helpwin, A_BOLD);
				}
				if ($show_node) {
					addstr $helpwin, "\nLimiting job view to node(s) ";
					attron($helpwin, A_BOLD);
					addstr $helpwin, $show_node;
					attroff($helpwin, A_BOLD);
				}

				#				addstr $helpwin, <<"__EOHELP__";
				#
				#
				#Interactive commands are:
				#
				# space   Update Display
				# /       Search for a server, node, or job and display details
				# q       Quit
				# h       Print this help
				# c       Grid Columns
				# u       Limit view to specific users' jobs
				# s       Seconds to refresh,
				#            accepts math operators (ie: 2*60)
				# C       Toggle Colorization
				# S       Toggle State Summary
				# G       Toggle Grid Display
				# Q       Toggle Queue Display
				# t       Toggle Queued Jobs in Queue Display
				# J       Toggle Show Jobs in Grid
				# 0-9     CPU number to display
				# l       Node load report
				#
				#Press any key to continue...
				#__EOHELP__
				addstr $helpwin, <<"__EOHELP__";


Interactive commands are:

 space   Update Display.
 /       Search for a server, node, or job and display details.
 q       Quit.
 h       Print this help.
 c       Grid Columns.
 u       Limit view to specific users' jobs.
 n       Limit view to specific nodes.
 s       Seconds to refresh.
            Accepts math operators (ie: 2*60)
 C       Toggle Colorization.
 S       Toggle State Summary.
 G       Toggle Grid Display.
 Q       Toggle Queue Display.
 t       Toggle Queued Jobs in Queue Display.
            When view is limited to specific nodes, 
            only jobs running on those nodes are displayed, 
            queued jobs are not displayed and t has no effect.
 J       Toggle Show Jobs in Grid.
 l       Node load report.

Press any key to continue...
__EOHELP__
				refresh($helpwin);

				# wait for the user to hit the any key.
				cbreak;
				nodelay($helpwin, 0);
				getch($helpwin);
				halfdelay(1);
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
				move($cmdwin, 0, 0);
				clrtoeol($cmdwin);
				refresh($cmdwin);

				delwin($helpwin);
			}

			# change the visible CPUs
			#FIXME# Don't allow the user to display CPUs that don't exist on any node
			elsif ($input =~ /^\d$/) {
				if (grep /^$input$/, @show_cpu) {
					@show_cpu = grep !/$input$/, @show_cpu;
				} else {
					my %seen = ();
					foreach (@show_cpu, $input) {
						$seen{$_} = 1;
					}
					@show_cpu = sort keys %seen;
				}

				update_display(@_);
			}

			elsif ($input eq "s") {
				$input = getstring("Number of seconds for refresh[$sleeptime]? ");

				if ($input) {
					my $tmp;

					# *grin* I love this use of eval
					if ($tmp = eval $input and $tmp > 0) {
						$sleeptime  = $tmp;
						$targettime = time() + $sleeptime;
					} else {
						printwarning("Invalid number!");
					}
				}

			}

			elsif ($input eq "c") {
				$input = getstring("Number of columns[$columns]? ");

				if ($input) {
					if ($input =~ /^\d+$/ and $input > 0) {
						$columns = $input;
						update_display(@_);
					} else {
						printwarning("Invalid number!");
					}
				}

			} elsif ($input eq "u") {

				$input = getstring("Limit view to \"all\", \"me\", a or username? ");

				if ($input) {
					if ($input eq "all" or $input eq "a") {
						$show_user = 0;
					} elsif ($input eq "me" or $input eq "m") {
						$show_user = $ENV{USER};
					} elsif ($input =~ /^\+(.*)/) {
						if ($1 eq "me" or $1 eq "m") {
							$show_user .= " " . $ENV{USER};
						} else {
							$show_user .= " $1";
						}
						$show_user =~ s/^ / /g;
					} elsif ($input =~ /^-(.*)/) {
						if ($1 eq "me" or $1 eq "m") {
							$show_user =~ s/\b$ENV{USER}\b//;
						} else {
							$show_user =~ s/\b$1\b//;
						}
						$show_user =~ s/  / /g;
						$show_user =~ s/^ / /g;
						$show_user =~ s/ $/ /g;
					} else {
						$show_user = $input;
					}
					update_display(@_);
				}

			} elsif ($input eq "n") {

				$input = getstring("Limit view to \"all\", a or nodename? ");

				if ($input) {
					if ($input eq "all" or $input eq "a") {
						$show_node = 0;
					} elsif ($input =~ m/^\+(.*)/) {
						$show_node .= " $1";
						$show_node =~ s/^ / /g;
					} elsif ($input =~ m/^-(.*)/) {
						$show_node =~ s/\b$1\b//;
						$show_node =~ s/  / /g;
						$show_node =~ s/^ / /g;
						$show_node =~ s/ $/ /g;
					} elsif ($input =~ m/^([0-9]+)-([0-9]+)$/) {
						$show_node = join(" ", $1 .. $2);
					} else {
						$show_node = $input;
					}
					update_display(@_);
				}
			}

			#			elsif ( $input eq "n" ) {
			#				$input = getstring("Node name? ");
			#
			#				if ($input) {
			#					foreach my $server ( %{ $_[1] } ) {
			#						if ( exists $_[1]->{$server}{$input} ) {
			#							print_nodestatus_window( $input,
			#								$_[1]->{$server}{$input} );
			#						}
			#					}
			#				}
			#			}
			elsif ( $input eq "j"
				and $subwin
				and $searchobject{TYPE} eq "NODE")
			{

				my $jobid = $_[1]->{ $searchobject{SERVER} }{ $searchobject{VALUE} }{job}{"0"};
				$searchobject{TYPE}   = "JOB";
				$searchobject{VALUE}  = $jobid;
				$searchobject{SERVER} = $_[2]->{$jobid}{server};
				$py                   = 0;
				$px                   = 0;
				update_subwin(@_);
				doupdate();
			} elsif ($input eq "l") {

				if ($subwin and $searchobject{TYPE} eq "JOBLOAD") {

					# We are currently looking a jobload report detail, immediately switch
					# this to a job detail window
					$searchobject{TYPE} = "JOB";
					$py                 = 0;
					$px                 = 0;
					update_subwin(@_);
					doupdate();

				} elsif ($subwin and $searchobject{TYPE} eq "JOB") {

					# We are currently looking at a job detail, immediately switch this
					# to a jobload report
					$searchobject{TYPE} = "JOBLOAD";
					$py                 = 0;
					$px                 = 0;
					update_subwin(@_);
					doupdate();

				} else {
					$input = getstring("Job ID Number? ");
					if ($input) {
						my @objects      = ();
						my $searchtype   = "";
						my $searchserver = "";

						# what other information can we extract?
						if ($input =~ /^\d+$/) {
							$searchtype ||= "job";
						} elsif ($input =~ /^\d+\./) {
							$searchtype ||= "job";
							($input, $searchserver) = split(/\./, $input, 2);
						}

						# we know everything we can, now go find stuff
						if ($searchtype eq "job") {
							foreach my $job (%{ $_[2] }) {
								if (
									$input eq $job
									and ( !$searchserver
										or $_[2]->{$job}{server} =~ m/^$searchserver/)
								  )
								{
									push(@objects, $job);
									$searchobject{TYPE}   = "JOBLOAD";
									$searchobject{SERVER} = $_[2]->{$job}{server};
									$searchobject{VALUE}  = "$job";
								}
							}
						}

						# if we have anything useful, go display it
						if (scalar @objects > 1) {
							printwarning("Multiple objects found. Please narrow your search.");
						} elsif (scalar @objects < 1) {
							printwarning("No objects found.");
						} elsif (exists $searchobject{TYPE} && defined $searchobject{TYPE}) {
							$subwin = 1;
							$py     = 0;
							$px     = 0;
							update_subwin(@_);
							doupdate();
						}
					}
				}
			}

			# Just about all of this should be moved out of here
			#FIXME# need to unify these searches somehow
			elsif ($input eq "/") {
				$input = getstring("Search string? ");

				if ($input) {
					my @objects      = ();
					my $searchtype   = "";
					my $searchserver = "";

					# did the user specify a pattern?
					if ($input =~ s/^~(.)\s?//) {
						if ($1 eq "s") {
							$searchtype = "server";
						} elsif ($1 eq "j") {
							$searchtype = "job";
						} elsif ($1 eq "n") {
							$searchtype = "node";
						} else {
							printwarning("Invalid search pattern");
							next;
						}
					}

					# what other information can we extract?
					if ($input =~ /^\d+$/) {
						$searchtype ||= "job";
					} elsif ($input =~ /^\d+\./) {
						$searchtype ||= "job";
						($input, $searchserver) = split(/\./, $input, 2);
					} elsif ($input =~ /^[a-z]+\s+\w/) {
						$searchtype ||= "node";
						($input, $searchserver) = split(/\s+/, $input, 2);
					}

					# we know everything we can, now go find stuff
					if ($searchtype eq "job") {
						foreach my $job (%{ $_[2] }) {
							if (
								$input eq $job
								and ( !$searchserver
									or $_[2]->{$job}{server} =~ /^$searchserver/)
							  )
							{
								push(@objects, $job);
								$searchobject{TYPE}   = "JOB";
								$searchobject{SERVER} = $_[2]->{$job}{server};
								$searchobject{VALUE}  = "$job";
							}
						}
					} else {
						foreach my $server (%{ $_[1] }) {
							if ($server =~ /^$input/
								and (!$searchtype or $searchtype eq "server"))
							{
								push(@objects, $input);
								$searchobject{TYPE}   = "SERVER";
								$searchobject{SERVER} = "$server";
								$searchobject{VALUE}  = "$server";
							}
							if (
								    exists $_[1]->{$server}{$input}
								and (!$searchtype or $searchtype eq "node")
								and ( !$searchserver
									or $server =~ /^$searchserver/)
							  )
							{
								push(@objects, $input);
								$searchobject{TYPE}   = "NODE";
								$searchobject{SERVER} = "$server";
								$searchobject{VALUE}  = "$input";
							}
						}
					}

					# if we have anything useful, go display it
					if (scalar @objects > 1) {
						printwarning("Multiple objects found.  Please narrow your search.");
					} elsif (scalar @objects < 1) {
						printwarning("no objects found matching $input");
					} elsif (exists $searchobject{TYPE}
						&& defined $searchobject{TYPE})
					{
						$subwin = 1;
						$py     = 0;
						$px     = 0;
						update_subwin(@_);
						doupdate();
					}
				}
			}

			#FIXME# had to remove the offline/clear options when I started
			#FIXME# to support multiple servers.  I haven't figured out
			#FIXME# a decent replacement interface yet
			#elsif ( $input eq "o" ) {
			#$input=getstring("Node to mark offline? ");
			#
			#if ($input) {
			#if ( exists $_[1]->{$input} ) {
			#system( $pbsnodes, "-o", "$input", "-s", "$host" );
			#}
			#else {
			#printwarning("Invalid nodename!");
			#}
			#}
			#}
			#elsif ( $input eq "r" ) {
			#$input=getstring("Node to clear offline? ");
			#
			#if ($input) {
			#if ( exists $_[1]->{$input} ) {
			#system( $pbsnodes, "-r", "$input", "-s", "$host" );
			#}
			#else {
			#printwarning("Invalid nodename!");
			#}
			#}
			#}

			# all of these toggles should be self-explanatory
			elsif ($input eq "C") {
				$colorize = !$colorize;
				if ($colorize && !has_colors()) {
					printwarning("Terminal doesn't support colors");
					$colorize = 0;
				}
				update_display(@_);
			} elsif ($input eq "b") {
				if (!has_colors()) {
					printwarning("Terminal doesn't support colors");
					$colorize = 0;
				} else {
					if (!$colorize) {

						# Toggle colorize if color was disabled.
						$colorize = !$colorize;
					}
					$input = getstring("Color jobs by \"j\" (job), \"u\" (user) or \"q\" (queue)? ");
					if ($input) {
						if ($input eq 'job' or $input eq 'j') {
							$color_by = 'job';
						} elsif ($input eq 'user' or $input eq 'u') {
							$color_by = 'user';
						} elsif ($input eq 'queue' or $input eq 'q') {
							$color_by = 'queue';
						} else {
							printwarning("Unsupported color by value.");
						}
						update_display(@_);
					}
				}
			} elsif ($input eq "G") {
				$show_grid = !$show_grid;
				update_display(@_);
			} elsif ($input eq "S") {
				$show_summary = !$show_summary;
				update_display(@_);
			} elsif ($input eq "Q") {
				$show_queue = !$show_queue;
				update_display(@_);
			} elsif ($input eq "t") {
				$show_qqueue = !$show_qqueue;
				update_display(@_);
			} elsif ($input eq "J") {
				$show_jobs = !$show_jobs;
				update_display(@_);
			}

			#FIXME# my home keyboard sends FIND and SELECT instead of HOME and END, weird?
			elsif ($input eq KEY_HOME
				or $input eq KEY_SHOME
				or $input eq KEY_FIND)
			{
				$py = 0;
				$px = 0;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq KEY_END
				or $input eq KEY_SEND
				or $input eq KEY_SELECT)
			{
				$py = ($ly > $subY ? $ly : $subY) + 2 - $Y;
				$px = 0;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq KEY_PPAGE
				or $input eq KEY_SPREVIOUS
				or $input eq $CTRL_B)
			{
				$py -= $Y - 2;
				$py <= 0 and $py = 0;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq KEY_NPAGE
				or $input eq KEY_SNEXT
				or $input eq $CTRL_F)
			{
				$py += $Y - 2;
				$py >= ($ly > $subY ? $ly : $subY) - $Y + 2
				  and $py = ($ly > $subY ? $ly : $subY) + 2 - $Y;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq "k" or $input eq KEY_UP) {
				$py <= 0 and $py = 0, next;
				$py--;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq "j" or $input eq KEY_DOWN) {
				$py >= ($ly > $subY ? $ly : $subY) - $Y + 2 and next;
				$py++;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq "h" or $input eq KEY_LEFT) {
				$px <= 0 and $px = 0, next;
				$px -= 2;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			} elsif ($input eq "l" or $input eq KEY_RIGHT) {
				$px >= $lx - $X + 1 and next;
				$px += 2;
				prefresh($pad, $py, $px, 0, 0, $Y - 2, $X - 1);
			}

			elsif ($input eq " ") {
				addstr $cmdwin, 0, 0, "Updating...";
				clrtoeol($cmdwin);
				refresh($cmdwin);
				move($cmdwin, 0, 0);
				clrtoeol($cmdwin);
				return;
			}

			#else {
			#addstr $cmdwin, 0, 0, ord($input);
			#}

		}

	}
}

sub readrc {
	my ($f) = shift or return;
	return unless (-f $f);

	open(F, $f) or die "$f: $!\n";
	while ($_ = <F>) {
		chomp;
		s/#.*//;
		next unless $_;
		my ($name, $value) = split(/=/);
		$name  =~ s/\s//g;
		$value =~ s/\s//g;
		next if (length($name) <= 0 or length($value) <= 0);

		if ($name eq "host" or $name eq "show_cpu") {

			# /em cringes
			eval "\@$name=(\"" . join('","', split(',', $value)) . "\")";
		} elsif ($name eq "show_user") {
			eval "\$$name=join(' ', split(',', \"$value\"))";
		} else {
			eval "\$$name=$value";
		}
	}
	close F;
}

__END__

=head1 NAME

qtop - monitoring utility for OpenPBS or Torque

=head1 SYNOPSIS

qtop [OPTION]... [@hostname]...

=head1 DESCRIPTION

Draws a full-terminal display of your nodes and jobs.  The default grid
shows each node's 1st CPU as a single character.  The specific character
denotes the state of the node or identifies the job running on that CPU.  The
job listing shows the job name, queue name, state, etc. and, on the far left,
the character used to identify nodes in the upper grid.  Pressing a number key
will toggle the display of that CPU on all of the nodes.

This program runs best if the C<perl-PBS> module is installed.  While there are
currently no loss of features if it isn't installed, it will run much faster
with it.  If you are unsure if PBS is installed, run this program, hit C<h>, and
look for the B<Backend> information at the top right.

=head1 COMMAND-LINE OPTIONS

=over 4

=item   B<-s> num

seconds between refreshes

=item   B<-c> num

number of columns to display in the grid

=item   B<-m> num

max number of cpus in a node before it gets its own grid

=item   B<-u>

usernames for limiting the view of the grid and job list.  Can be a
comma-seperated list of usernames or C<all>.  C<me> is a pseudonym for the
username running qtop(1).

=item   B<-C>

toggle colorization

=item   B<-S>

toggle state summary display

=item   B<-G>

toggle grid display

=item   B<-Q>

toggle queue display

=item   B<-t>

toggle showing queued jobs in queue display

=item   B<-[0-9]...>

cpu numbers for grid display

=item   B<-J>

toggle jobs in grid display

=item   B<-V>

print version and exit

=back

=head1 INTERACTIVE COMMANDS

Several single-key commands are recognized while qtop(1) is running.  The
arrow keys, PageUp, and PageDown keys will scroll the display if it doesn't fit
in your terminal.

When prompted to type something, ctrl-g can be used to cancel the command.

=over 4

=item   B<space>

Immediately update display

=item B<q>

Quit qtop(1)

=item B<h>

Display help screen, version, and current settings

=item B<c>

Prompts for the number of columns to display the node grid

=item B<s>

Prompts for the number of seconds to wait between display updates

=item B<u>

Prompts for a username.  The grid and job listing will be limited to the named
user.  Input C<all> will remove all limitations (the default), and C<me> will
limit to the current username running qtop(1).  If the username or C<me> is
prefixed with a C<+> or C<->, the username will be added or removed from the
list of usernames to be limited.  C<a> and C<m> are shortcuts for C<all> and
C<me>.

=item B</>

Prompts the user for a search string, for displaying the details of.  The
search can optionally begin with one of the following pattern specifiers
(think: mutt): C<~s> for a server, C<~n> for a node, or C<~j> for a job number.
If no pattern specifier is found, qtop will attempt to find the object that
best matches the search string. The string can be a server name, nodename, or a
job number.  Nodenames can optionally be followed by a space and the server
name.  Job numbers may optionally be followed by a dot and the server name.

If an object is found, a subwindow will be opened displaying details.  Hit C<q>
to exit the window.

When viewing a job detail subwindow, pressing C<l> is a shortcut for jumping
directly to the associated job's node load subwindow.

(Mnemonic: like using / to search for text in vi or less)

=item B<l>

Prompts the user for a job id.  A B<node load report> subwindow will be
displayed for the given jobid.  This subwindow shows the current load average,
the physical and available memory, and the number of sessions.  Available
physical memory will be negative in the event of swapping.  If the number of
sessions is 0, that might indicate a problem on that node.

Pressing C<l> in this subwindow jumps you directly to the associated job detail
subwindow; as if the user typed C</jobid>.

(Mnemonic: load average)

=item B<C>

Toggle the use of the colors in the display

=item B<S>

Toggle the display of the state summary

=item B<G>

Toggle the display of the node grid

=item B<Q>

Toggle the display of the job queue

=item B<t>

Toggle the display of currently queued (not running) jobs in the display.  This
can reduce the size of the queue display considerably in some environments.

(Mnemonic: I don't know, toggle?  C<Q> was already used for something more important)

=item B<J>

Toggle the display of job letters in the node grid.  This handy because you can
see the node state "hidden" behind the job letter.  For example, use this to
see which nodes are not yet "busy" that have jobs.

=for comment
# SMP stuff disabled in current version:
#= item B<Any single number (0-9)>
#
#Toggle display of that CPU number in the display.  This is confusing at first,
#but useful in SMP environments (See SMP section below).
#

=back

=head1 STARTUP

qtop(1) has many configuration variables that can set on the command line,
interactively, or from configuration files.  When qtop(1) starts, it first
initializes these variables with built-in defaults, then reads in
F</etc/qtoprc>, the reads F<~/.qtoprc>, and finally parses the command line
arguments.  Note that several of the command line arguments and interactive
commands are toggles, they don't directly set the value of the configuration.
In contrast, the configuration files are not toggles.

The configuration files may contain following name=value pairs:

=over 4

=item B<columns>

Number of columns in the node grid, positive integer

=item B<sleeptime>

Number of seconds to pause between display updates, positive integer

=item B<colorize>

Use colors in the display, 1 or 0

=item B<show_summary>

Display the summary at the top of the display, 1 or 0

=item B<compact_summary>

Show node state summary on one line, 1 or 0

=item B<show_grid>

Show the node grid, 1 or 0

=item B<show_queue>

Show the job queue, 1 or 0

=item B<show_qqueue>

Show queued (not running) jobs in the queue display, 1 or 0

=item B<show_jobs>

Show job and color information in the node grid, 1 or 0

=for comment
# SMP stuff disabled in current version:
#=item B<show_cpu>
#
#Comma seperated list of CPU numbers to display
#
=item B<show_user>

Usernames to limit the view in the grid and job list.  Can be a comma-seperated
list of users, C<all>, or C<me>.

It might be reasonable for a site to have C<show_user=me> in F</etc/qtoprc>
and for admin users to have C<show_user=all> in their own F<~/.qtoprc>.

Members of a group might want all of their teammate's usernames in their own
F<~/.qtoprc>.

=item B<host>

Comma seperated list of hostnames running pbs_server

=item B<maxrows>

Number of rows in the large scrollable panel

=item B<maxcolums>

Number of columns in the large scollable panel

=item B<maxnodegrid>

Nodes with more than this number of CPUs will be represented by a seperate grid

=back

A sample configuration file:

    # I'm grumpy and don't like color
    colorize=0

    # my 6 CPU machine should get a seperate grid
    maxnodegrid=5

    # all of my Torque servers
    host=teraserver,bigbird,testhpc

=head1 SMP ENVIRONMENTS

qtop(1) was developed with three specific clusters in mind, these are a 1000
node cluster of dual SMP machines, a 64 proc SMP with 16 single node machines,
and a 21 node cluster of single procs without nicely numbered hostnames.  With
this kind of pedigree, qtop(1) is fairly flexible.

By default, the node grid will show the state of the first CPU of 30 nodes in
each row.  The number of columns in the grid can be shrunk or expanded on the
command line with C<-C>, or interactively with C<c>.  Additional CPUs can be
displayed by pressing the appropriate number key.  Using the number keys is
confusing at first, but if you try it a few times it will became natural.
By default, nodes with 8 or more CPUs are displayed in a seperate grid.

The first two clusters mentioned above display well with the defaults.  The
third is typically displayed with the number of columns set to "1".

=head1 FILES

=over 4

=item F</etc/qtoprc>

The global configuration file  

=item F<~/.qtoprc>

The personal configuration file.

=back

=head1 ENVIRONMENTAL VARIABLES

=over 4

=item PBS_DEFAULT  

The server's hostname (same as most PBS client commands)

=back

=head1 SEE ALSO

=over 4

=item PBS(3pm), qstat(1B)

=back

=head1 BUGS

The large Job structure uses the servername supplied by the user, the Job
structure uses the servername returned by the server... so they don't match up
(this makes the jobloadreport imprecise).  
The curses code is very ineffecient, the screen flickers too much.
grep FIXME from qtop for more!

=head1 AUTHOR

qtop(1) was originally written by Garrick Staples E<lt>garrick@usc.eduE<gt>.
The node grid and lettering concept is from Dennis Smith.  Thanks to Egan Ford
and the xCAT mailing list for testing and feedback.

