#!/usr/bin/env perl

# Some auxiliary functions, aiming to ease the development, but unrelated
# to the concept of Hidden Markov Models

package AUX;

our $minusinf = -9999999;

# Given a numerical array, return the index of the biggest element in it.
sub maxIndex {
  my @array = @{shift @_};
  my $max = $array[0];
  my $max_index = 0;
  for my $i (0..$#array) {
    if ($array[$i] > $max) {
      $max = $array[$i];
      $max_index = $i;
    }
  }
  return $max_index;
}

# Given an array and an scalar, return the index of this scalar in the
# array. If this scalar appears more than once, return the index
# corresponding to its first appearance.
sub getIndex {
  my @array = @{shift @_};
  my $item = shift @_;
  my $i;
  for ($i = 0; $i <= $#array; $i++) {
    if ($item eq $array[$i]) {
      return $i;
    }
  }
}

# (for debugging) Just take a matrix and print it
sub printMatrix {
    my @matrix = @{shift @_};
    my $i;
    my $j;
    for ($j = 0; $j < scalar(@{$matrix[0]}); $j++) {
        for ($i = 0; $i < scalar(@matrix); $i++) {
            printf "%.3e ", $matrix[$i][$j];
         }
         print "\n";
    }
    print "\n";
}

# As log, but log(0) = -99999..
sub milog {
  my $x = shift;
  return $x == 0 ? $minusinf: log($x); 
}

1;
