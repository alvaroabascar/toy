#!/usr/bin/env perl

use strict;
use warnings;
use MODELS;
use AUX;

package HMM;

# minusinfinite
our $minusinf = -9999999;

# model chosen by default
our %hmm = %{$MODELS::models{toy}};

# call this function to choose another model
sub select_model {
  my $model = shift @_;
  # if the name is not recognized, issue an error and show available models
  if (not exists($MODELS::models{$model})) {
  print "Model $model is not available. Please choose among the following:\n";
  foreach my $key (keys %MODELS::models) {
    print "$key ";
  }
  print "\n";
  die();
  }
  # if the name is recognized, select the model
  else {
  %hmm = %{$MODELS::models{$model}};
  }
}

# given a state name, return its index
sub getIndexOfState {
  my $state = shift @_;
  my $all_states = $hmm{states};
  return AUX::getIndex($all_states, $state);
}

# given a state index, return its name
sub getStateByIndex {
  my $index = shift;
  return $hmm{states}->[$index];
}

# given a symbol index, return its name
sub getSymbolByIndex {
  my $index = shift;
  return $hmm{symbols}->[$index];
}

# given a state name, return the array of EMISSION probabilities
# of this state
sub getProbabilitiesOfState {
  my $currentstate = shift @_;
  my $index = AUX::getIndex($hmm{states}, $currentstate);
  my $probabilities = $hmm{emp}->[$index];

  return @{$probabilities};
}

# given a state name, return the array of TRANSITION probabilities
# from this state
sub getProbabilitiesOfTransition {
  my $currentstate = shift @_;
  my $index = AUX::getIndex($hmm{states}, $currentstate);
  my $probabilities = $hmm{trp}->[$index];

  return @{$probabilities};
}

# given a state name and a symbol name, return the emission probability
# for this symbol, from this state
sub getEmissionProbability {
  my $currentstate = shift @_;
  my $symbol = shift @_;
  my $symbolindex = AUX::getIndex ($hmm{symbols}, $symbol);
  my $currentindex = AUX::getIndex ($hmm{states}, $currentstate);
  return $hmm{emp}[$currentindex][$symbolindex];

}

# Given a matrix of pointers and a matrix of log probabilities built by
# viterbi, return the most probable path (as an array of state names)
sub getPath {
  my $dic = shift @_;
  my $matrix = $dic->{matrix};
  my $pointers = $dic->{pointers};

  my $i = @{$matrix} - 1;
  my $laststate = getIndexOfState("end");
  my @path;
  for (my $state = $laststate; $i >= 0; $i--){
  unshift @path, getStateByIndex($state);
  $state = $pointers -> [$i][$state];
  }
  return @path;
}

# Given a sequence, use the Viterbi algorithm to fill a matrix of probabilities
# and pointers, then call getPath to obtain the most probable path and return
# this path as an array of state names.
sub viterbi {
  my $sequence = shift @_;
  my @path;
  my @matrix;
  my $emission;
  my $transition;
  my @pointers;

  #initialization step
  $matrix[0] = [];
  $pointers[0] = [];

  # fill the column for the "begin" position
  # (all log probs are -inf except for the begin state, which has 0)
  foreach my $j (0..scalar(@{$hmm{states}})-1) {
    # for the pointers, set an invalid number
    $pointers[0][$j] = -1;

    if (getStateByIndex($j) eq $hmm{initstate}){
      $matrix[0][$j] = 0 ;
    }
    else {
      $matrix[0][$j] = $minusinf;
    }
  }

  # Recursive step
  my $symbol;
  my $state;
  # iterate through the sequence...
  foreach my $i (1..length($sequence)+1) {
    # iterate through the states...
    foreach my $j (0..scalar(@{$hmm{states}})-1) {
      $symbol = substr($sequence, $i-1, 1);
      $state = getStateByIndex($j);;
      $emission = getEmissionProbability($state, $symbol);

      # can't emmit this nucleotide? then the log prob of this cell
      # is -inf
      if ($emission == 0) {
        $matrix[$i][$j] = $minusinf;
        $pointers[$i][$j] = -1;
        next;
      }

      # what is the previous state which gives the maximum probability
      # transition * emission?
      $transition = bestTransition($matrix[$i-1], $state);

      # all the probabilities of transition are 0?
      if ($transition->{probability} <= $minusinf) {
        $matrix[$i][$j] = $minusinf;
        $pointers[$i][$j] = -1;
        next;
      # we got the best feasible transition
      } else {
        # log prob of this cell is the product of the transition probability
        # by the emissio probability (or sum, if working with logs)
        $matrix[$i][$j] = AUX::milog($emission) + $transition->{probability} ;
        $pointers[$i][$j] = $transition->{pointer};
      }
    }
  }

  # Termination step
  # we are in the "end" (fictional) position of the sequence
  my $i = length($sequence) + 1;
  # iterate through all the states...
  foreach my $j (0..scalar(@{$hmm{states}})-1) {
    $state = getStateByIndex($j);
    # get the best transition to $state from previous column
    $transition = bestTransition($matrix[$i-1], $state);

    # all transitions are impossible?
    if ($transition->{probability} <= $minusinf) {
    $matrix[$i][$j] = $minusinf;
    $pointers[$i][$j] = -1;
    next;
     } else {
    $matrix[$i][$j] = $transition->{probability} ;
    $pointers[$i][$j] = $transition->{pointer};
    }
    }

  return getPath {matrix => \@matrix, pointers => \@pointers};
}

# bestTransition:
#   @probabilities, $state (string)
#   returns:
#    {  probability => imum probability product:
#              PROB_k(i-1) * a_kl
#       pointer => k
#    }
#
# Given the probabilities of the previous column "i-1", and the current state
# $currentstate, find the state "k" which minimizes the product
# PROB_k(i-1) * a_kl (transition from state "k" to current state "l")
#
sub bestTransition {
  my @probabilities = @{shift @_};
  my $currentstate = shift @_;
  my $currentindex;
  my $transition;
  my @total_probs = ();
  # create an array with the product log(probability) of a cell +
  # log(prob of transition to current state)
  foreach my $j (0..scalar(@probabilities)-1) {
  $currentindex = getIndexOfState($currentstate);
  $transition = $hmm{trp}[$j][$currentindex];
  $total_probs[$j] = $probabilities[$j] + AUX::milog($transition);
  }
  # get the imum of these probabilities (by index)
  my $max_index = AUX::maxIndex(\@total_probs);
  return {pointer => $max_index, probability => $total_probs[$max_index]};
}

# Like calling viterbi, but prints the sequence and a string representing the
# path (each letter is the first letter of the state)
sub predict {
  my $seq = shift;
  my @path = viterbi $seq;
  print $seq, "\n";
  for my $i (1..length $seq) {
  print substr $path[$i], 0, 1;
  }
  print "\n";
  return @path;
}

# Given a file with test sequences (one sequence per file), and
# an integer indicating the position of the donor site in this
# sequences, run the HMM on all of them and return a dictionary
# {SN => sensitivity, SP => specificity,
#  TP => true positives, FP => false positives};
sub test {
  my $file = shift;
  my $donor_position = shift;
  my $verbose = shift;

  open FILE, "<$file";
  my $TP=0;
  my $FP=0;
  my $SN;
  my $SP;

  # if called with verbose=1, use the function predict, which shows
  # the prediction for each function. Else use the function viterbi,
  # which does the same but silently
  my $predict_fn;
  if ($verbose) {
    $predict_fn = \&predict;
  } else {
    $predict_fn = \&viterbi;
  }


  # for each sequence in the file...
  while (<FILE>){
    chomp $_;
    # get the path
    my @path = &$predict_fn($_);

    # did the model predict the donor position correctly?
    if ($path[$donor_position] eq "Gdonor"){
      $TP++;
    }
    else{
      $FP++;
    }
  }

  # calculate sensitivity and specificity (the same, in this case)
  $SN = $SP = $TP/($TP + $FP);

  print "Sensitivity: $SN\nSpecificity: $SP\nTrue positive: $TP\nFalse Positive/False negative: $FP\n";
  return {SN => $SN, SP => $SP, TP => $TP, FP => $FP};
}

1;
