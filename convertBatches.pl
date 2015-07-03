#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

my $reencode = $ARGV[0];
my $batches = $ARGV[1];

open (REENCODE, "< $reencode") or die $!;
open (BATCHES, "< $batches") or die $!;

my $reencodeline;
my %reencode;
while ($reencodeline = <REENCODE>){
    chomp $reencodeline;
    #print "$reencodeline\n";
    my @array = split('\t', $reencodeline);
    my $ll_id = $array[1];
    my $id = $array[3];
    $reencode{ $ll_id } = $id;
}

my $batchesline;
my %batches;
while ($batchesline = <BATCHES>){
    chomp $batchesline;
    my @array = split('\t', $batchesline);
    my $batch_id = $array[0];
    my $ll_id = $array[1];
    $batches{ $ll_id } = $batch_id;
}

foreach my $key (%batches){
    if (exists($reencode{$key})){
        print $batches{$key} . "\t" . $reencode{$key} . "\n";
    }  
}
