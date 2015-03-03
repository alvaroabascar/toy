#!/usr/bin/env perl

use strict;
use warnings;
use HMM;

# Sample a sequence of n nucleotides. This script must be called as:
# ./sampling.pl n

# Given a state, sample the next state, according to the probabilities of the HMM.
sub sampleNextState {
  my $currentstate = shift @_;
  my @probabilities = HMM::getProbabilitiesOfTransition($currentstate);
  my $index = getRandomIndex(\@probabilities);

  return HMM::getStateByIndex($index);
}

# Given a state, sample an observation, according to the probabilities of the HMM.
sub sampleNextObservation {
  my $currentstate = shift @_;
  my @probabilities = HMM::getProbabilitiesOfState($currentstate);
  my $index = getRandomIndex(\@probabilities);
  return HMM::getSymbolByIndex($index);
}

# Given an array of probabilities, return an integer "i" in the range
# [ 0, length(array) ), such that p(i) = array[i].
#
# For example, if this function is called with the array (0.1, 0.2, 0.7),
# it will return 0 the 10% of the times, 1 the 20% of them, and 2 the 70%.
sub getRandomIndex {
  my @probabilities = @{shift @_};
  my $random = rand();
  my $sum_probabilities = 0;

  for (0..$#probabilities) {
    $sum_probabilities += $probabilities[$_];
    if ($random < $sum_probabilities) {
      return $_;
    }
  }
}

# Sample a sequence of n observations.
sub sampling {
  my $n = shift @_;
  my $currentstate = 'begin';
  my $sequence;
  my $observation;
  for (1..$n) {
    $currentstate = sampleNextState($currentstate);
    $observation = sampleNextObservation($currentstate);
    $sequence .= $observation;
  }
  return $sequence;
}

# Select a model, by default, the most basic one.
HMM::select_model("TIA2");

# 
print sampling(shift @ARGV), "\n";
