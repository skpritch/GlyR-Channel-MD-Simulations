#!/usr/bin/perl -w

use strict;

# This script processes all .mdp files in the MDP directory and generates modified .mdp files
# with 'pull_coord1_init' values adjusted according to indices from -48 to 48.

# Input directory containing template .mdp files
my $input_dir = 'cg_output/MDP';

# Output directory where the generated .mdp files will be stored
my $output_dir = 'david_MDP_CLA';
# my $output_dir = 'MDP_POT';
# my $output_dir = 'MDP_SOD';

# Check if the output directory exists; if not, create it
unless (-d $output_dir) {
    mkdir $output_dir or die "Cannot create directory $output_dir: $!";
}

# Define the specific .mdp files to process
my @files = ('NPT.mdp', 'PROD.mdp');

foreach my $file (@files) {
    my $key = "$input_dir/$file";
    my @temp = split('\.', $file);
    my $base = $temp[0];

    open(IN, "<$key") or die "Cannot open input file $key: $!";
    my @in = <IN>;
    close(IN);

    for (my $i = -45; $i < 46; $i++) {
        my $filename = "${base}_${i}.mdp";
        my $filepath = "$output_dir/$filename";
        open(OUT, ">$filepath") or die "Cannot open output file $filepath: $!";
        foreach $_ (@in) {
            unless ($_ =~ /^pull_coord1_init\s*=/) {
                print OUT $_;
            }
            if ($_ =~ /^pull_coord1_init\s*=/) {
                printf OUT "%s %f \n", $&, $i / 10;
            }
        }
        close(OUT);
    }
}

exit;
