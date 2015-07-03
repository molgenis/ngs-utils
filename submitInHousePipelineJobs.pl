#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Cwd;

use diagnostics;

####
# submitInHousePipelineJobs.pl
# WBKoetsier 20120928
####
# This script assesses in- and outputs of each of the given jobscripts and
# does or does not submit it to the cluster queue. 
# input: a list of jobscripts and a dependency list (might need a more elegant solution)
# output: log (no verbosity options) #### suggestion: be more verbose in general than I am now
# help: usage and exit.
# This script assumes it is run in the same directory as the jobscripts! (will report and exit if not)
####
# input jobscript list: one file per line
# example input jobscript list: 
#   $> ls s*_*_*.sh > input.list
# input dependency list: 'jobscript\tdependencyjobscript1:dependencyjobscript2:...'
# example input dependency list: 
#   $> awk '/qsub -N/
#    > {
#    >   sub( /)/, "" );
#    >   if ( $4=="-W" )
#    >   {
#    >     sub( /depend=afterok:/, "" );
#    >     gsub( /\$/, "" );
#    >     print $NF "\t" $5
#    >   } else {
#    >     print $NF 
#    >   }
#    > }' submit.sh > dependency.list 
 
#### Options
# get options
my $help;
my ($jobscript_list, $dependency_list);
my $logfile;
GetOptions(
  'h|help'	=> \$help,
  'j=s'		=> \$jobscript_list,
  'd=s'		=> \$dependency_list,
  'o=s'		=> \$logfile
);

# parse options
usage() and exit(0) if $help;
usage() and exit(0) unless $jobscript_list;
usage() and exit(0) unless $dependency_list;
$logfile = "submit.log" unless $logfile; 
my $log = "submitInHousePipelineJobs.pl run by " . getpwuid( $< ) . " at " . localtime() . "\n";
$log .= "Running on node: " . `hostname` . "\n";

#### Input
# load jobscripts list
my @jobscripts;
open ( JOBSCR, "<$jobscript_list" ) or die "Can't open $jobscript_list [$!]";
chomp ( @jobscripts = <JOBSCR> );
close ( JOBSCR );

# before continuing: check if the jobscripts are present in the cwd (or: is this the right project?)
my $cwd = cwd();
foreach ( @jobscripts ) {
  print "$_ does not exist in $cwd. Please check your input (or working directory). [$!]\n" and exit 1 if ! -e $_;
}

# load dependencies list
my %dependencies;
open ( DEPLST, "<$dependency_list" ) or die "Can't open $dependency_list [$!]";
while ( <DEPLST> ) {
  chomp;
  my ( $sc, $dp ) = split "\t";
  $dp = "" if !defined $dp;
  $dependencies{ $sc } = $dp;
}


#### Submit what has to be submitted
# get (unique) list of steps from jobscripts list
# step key is sNN
my %pipeline_steps;
foreach ( @jobscripts ) {
  my $st = ( split /\_/ )[0];
  push @{ $pipeline_steps{ $st } }, $_;
}
my @keys_pipeline_steps = sort keys %pipeline_steps;

# continue per step
my $abort = ":";
my %jobIDs;
$log .= "\n";
foreach my $stp ( @keys_pipeline_steps ) {

  $log .= "Step: $stp\n";

  # multiple jobscripts per step
  my $dep_list;
  my $ps_vals_ref = $pipeline_steps{ $stp };
  PSV: foreach my $psv ( @$ps_vals_ref ) {

    $log .= "Jobscript: $psv\n";

    # continue to the next jobscript if the previous one aborted
    # lazy: I did not link scriptnames <=s15 (1 to num bc/lanes) to >s15 (samplename)
    my $smpl = "$2" if $psv =~ /^($stp)_.+?(_.+?)\.sh$/;
    if ( "$abort" =~ /$smpl/ ) {
      $psv =~ s/\..+?$//;
      $jobIDs{ $psv } = "no job ID";
      next PSV;
    }

    # what dependencies does this jobscript have?
    $dep_list = $dependencies{ $psv };
    $dep_list = "" if !defined $dep_list;

    # grab inputs from jobscript (if any)
    my $inputs = extract_inputs_outputs ( "$cwd/$psv", "inputs" );

    # fetch dependency outputs and check existence
    foreach my $d ( split /:/, $dep_list ) {
      my $dep_outputs = extract_inputs_outputs ( "$cwd/$d.sh", "alloutputsexist" );
      foreach my $file ( @$dep_outputs ) {
        # only evaluate outputs that serve as input for $psv
        if ( grep /$file/, @$inputs ) {
          # allthough... I expect that we never see this, unless someone moved some data.
          unless ( -e $file ) {
            $log .= "Aborting submit of jobscript $psv: output\n";
            $log .= "\t$file\n";
            $log .= "of jobscript $d.sh (dependency of $psv) does not exist.\n";
            # and what of the other outputs of this dep? Want to report if they're missing?
            $abort .= "$2:" if $psv =~ /^($stp)_.+?(_.+?)$/;
            $psv =~ s/\..+?$//;
            $jobIDs{ $psv } = "no job ID";
            next PSV;
          }
          @$inputs = grep !/$file/, @$inputs;
        }
      }
    }

    # evaluate the other inputs of $psv
    foreach ( @$inputs ) {
      unless ( -e $_ ) {
        $log .= "Aborting submit of jobscript $psv: input\n";
        $log .= "\t$_\n";
        $log .= "of jobscript $psv does not exist.\n";
        $abort .= "$2:" if $psv =~ /^($stp)_.+?(_.+?)\.sh$/;
        $psv =~ s/\..+?$//;
        $jobIDs{ $psv } = "no job ID";
        next PSV;
      }
    }

    # evaluate the outputs of $psv
    my $psv_outputs = extract_inputs_outputs ( "$cwd/$psv", "alloutputsexist" );
    my $num_of_existing_outputs = 0;
    my $rm_str = "rm ";
    foreach ( @$psv_outputs ) {
      if ( -e $_ ) {
        $num_of_existing_outputs++;
        $rm_str .= "$_ ";
      }
    }

    if ( $num_of_existing_outputs == scalar @$psv_outputs ) {
      $log .= "Skipping submit of jobscript $psv: all outputs exist\n";
      # and prepend 'skipped' to .err and .out
      foreach ( (".out", ".err") ) {
        $psv =~ s/\..+?$/$_/;
        if ( -e $psv ) {
          # prepend
          #my $date = date(); # format like BASH date
          #my $msg = system ( "sed -i '1i skipped at $date' $cwd/$psv" );
          #$log .= "\"sed -i '1i skipped at $date' $cwd/$psv\" reported:\n[$msg]\n";
        } else {
          # new file
          #open ( OE, ">$cwd/$psv" ) or die "Could not open $cwd/$psv for writing: [$!]\n";
          #print OE "skipped\n";
          #close OE;
        }
      }
      $psv =~ s/\..+?$//;
      $jobIDs{ $psv } = "no job ID";
      next PSV;
    } else {
      $log .= "Some or no outputs of jobscript $psv exist. Removing remaining outputs.\n";
      # check string length just in case
      #my $msg = "";
      #$msg = system ( $rm_str ) if length $rm_str > 3;
      #$log .= "\"$rm_str\" reported:\n[$msg]\n";
    }

    # does dependency have a job ID? Add as -W
    my $afterok = "";
    foreach my $d ( split /:/, $dep_list ) {
      my $jid = $jobIDs{ $d };
      unless ( $jid eq "no job ID" ) {
        $afterok .= ":$jid";
      }
    }
    $afterok = "-W depend=afterok" . $afterok if length $afterok > 0;

    $psv =~ s/\.sh//;
    my $qsub = "qsub -N $psv $afterok $psv\.sh";

    # feature:
    # if the user already is at queue limits, wait for a slot before submitting or start using another queue
    # just run this script in the background and wait as you would wait for jobs to finish
    # (nooo wait, qsub this script and have it monitor node resources and run jobscripts directly on the node :P)

    # submit!
    $log .= "Submitting $qsub\n";
    #my $jobID = system ( $qsub );
    my $jobID = "12345";
    $jobIDs{ $psv } = $jobID;
    $log .= "JobID is: $jobID\n";


  } # foreach ps_vals_ref

  $log .= "\n\n";

} # foreach pipeline step


#### Write output
open ( OUTPUT, ">$logfile" ) or die "Can't open $logfile for writing [$!]";
print OUTPUT "$log\n"; 
close ( OUTPUT);


#### Sub: extract_outputs
# extract alloutputsexist from jobscript
# input: a valid jobscript filename (with full path)
# input: either 'inputs' or 'alloutputsexist'
# return: list of outputs defined in this jobscript
sub extract_inputs_outputs {

  my ( $scrpt, $type ) = @_;
  my $other_type = "alloutputsexist";
  $other_type = "inputs" if $type =~ /alloutputsexist/;

  # see if $type was defined in this script
  open ( SCRPT, "<$scrpt" ) or die "Could not open $scrpt for reading: [$!]\n";
  my $io_files_str = "";
  while ( <SCRPT> ) {
    chomp;
    if ( /$type/ && /[^#]/ ) {
      $io_files_str .= "$_ " if /\//;
      {
        my $nextline = <SCRPT>;
        chomp $nextline;
        if ( $nextline =~ /\// && $nextline !~ /$other_type/) {
          $io_files_str .= "$nextline ";
          redo;
        }
      } 
    }
  }
  close SCRPT;

  # we only want the files and their paths
  my @io_files = ();
  $io_files_str =~ s/$type//g;
  $io_files_str =~ s/" "/\n/g;
  $io_files_str =~ s/\\//g;
  $io_files_str =~ s/"//g;
  $io_files_str =~ s/\s+/ /g;
  $io_files_str =~ s/^\s*//g;
  @io_files = split /\s+/, $io_files_str;

  return \@io_files;

}


#### Sub: usage
# prints info message and usage
sub usage {
        print <<EOF;
#######################################################################################################
# submitInHousePipelineJobs.pl 
#######################################################################################################
# This script assesses in- and outputs of each of the given jobscripts and
# does or does not submit it to the cluster queue. 
# input: a list of jobscripts and a dependency list
# output: log (no verbosity options)
# help: usage and exit.
# This script assumes it is run in the same directory as the jobscripts!
#######################################################################################################
# Usage: perl createSubmitScript.pl -j jobscripts.list -d dependency.list -o submit.log
# \t -j \t list of jobscripts to submit
# \t -d \t list of dependencies that goes with the jobscript list
# \t -o \t logfile (defaults to submit.log)
# \t -h \t this message
#######################################################################################################
# input jobscript list: one file per line
# example input jobscript list: 
#   \$> ls s*_*_*.sh > input.list
# input dependency list: 'jobscript\tdependencyjobscript1:dependencyjobscript2:...'
# example input dependency list: 
#   \$> awk '/qsub -N/ { sub( /)/, "" ); if ( \$4=="-W" ) { sub( /depend=afterok:/, "" );
#    >   gsub( /\$/, "" ); print \$NF "\t" \$5 } else { print \$NF } }' submit.sh > dependency.list 
#######################################################################################################
EOF
} 
