#!/usr/bin/env perl

# Run a test on the desired test set, and print the results.
# ./test.pl
#
# If the script is executed with the argument "verbose", each
# sequence will be shown with the corresponding prediction

use strict;
use warnings;

use HMM;

# Select a HMM. You can find all the available models in
# MODELS.pm (see the end of the file)
HMM::select_model("U1_all");

my $verbose = 0;
if (scalar @ARGV > 0 && shift @ARGV eq "verbose") {
  $verbose = 1;
}

# Run the test on the sequences found in testset_full.txt
# (one sequence per line), taking into account that the
# donor G is in position 51 of ALL the sequences
#
print "Running test. This can take several minutes...\n";
HMM::test("tests/testset_full.txt", 51, $verbose);
