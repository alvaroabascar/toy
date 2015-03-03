#!/usr/bin/env perl

# A simple script used to calculate the Jensen - Shanon
# divergency among the distribution of nucleotides in the case
# of real and false donors, for each position, as well as the
# frequencies.

use strict;
use warnings;

# Given a file of sequences, calculate the frequency of each
# nucleotide at each position (using pseudocounts). Return an array
# of dictionaries {nucleotide -> frequency}
sub countfrequencies {

  my $file = shift;

  # an array of dicts, one per position in the sequences
  my @array;

  open(FILE, "<$file") or die($!);

  # initialize all dicts with a count of 0 for each nucletoide, then
  # add the counts for the first line.
  # this way we cope can cope with files with different sequence lengths.
  my $line = <FILE>;
  chomp $line;
  foreach my $i (0..length($line)-1){
    $array[$i] = {A => 0, T=>0, C=>0, G=>0};
    $array[$i]{substr($line,$i,1)}++;
  }

  my $numlines = 1;

  # for each sequence...
  while (<FILE>){
    chomp;
    # for each position in the sequence...
    foreach my $i (0..length($_)-1){
      # increase the count for this nucleotide
      $array[$i]{substr($_,$i,1)}++;
    }
    $numlines++;
  }

  # for each position...
  foreach my $i (0..$#array){
    # for each nucleotide...
    foreach my $key (keys %{$array[$i]}){
      # turn number of occurrences into frequencies (with pseudocounts)
      $array[$i]{$key} = ($array[$i]{$key} + 1)/($numlines + 4);
    }
  }
  return @array;
}

# Given two probability distributions (each represented by a dictionary
# nucleotide -> probability), calculate the Kullback-Leibler divergence.
sub KL_divergence {
  my %distrib1 = %{shift @_};
  my %distrib2 = %{shift @_};
  my $divergence = 0;
  foreach my $key (keys %distrib1) {
    $divergence += $distrib1{$key} * log($distrib1{$key} / $distrib2{$key});
  }
  return $divergence;
}

# Given two probability distributions (each represented by a dictionary
# nucleotide -> probability), calculate the Jensen-Shannon divergence.
sub JS_divergence {
  my @real = @{shift @_};
  my @false = @{shift @_};

  my @diverg;

  foreach my $position (0..$#real) {
    $diverg[$position] = 0.5 *
                         (KL_divergence($real[$position], $false[$position]) +
                         KL_divergence($false[$position], $real[$position]));
  }
  return @diverg;
}

# Take a reference to an array of floats and print them.
sub print_array {
  my $arref = shift @_;
  foreach my $item (@{$arref}) {
  printf "%.4f  ", $item;
  }
  print "\n";
}

# Pretty print for the frequencies at each position
sub print_freqs {
  my @freqs = @{shift @_};
  foreach my $key (keys %{$freqs[0]}) {
  print "$key: ";
  foreach my $i (0..$#freqs) {
    printf "%.4f  ", $freqs[$i]{$key};
  }
  print "\n";
  }
}

# Get the nucleotide frequencies per position in each file
my @real = countfrequencies("false_human_donors.txt");
my @false = countfrequencies("real_human_donors");
my @divergencies = JS_divergence(\@real, \@false);

# Print the frequencies
print "Frequencies in real_donors:\n";
print_freqs(\@real);
print "\nFrequencies in false_donors:\n";
print_freqs(\@false);
# Print the divergence
print "\nJensen-Shannon divergence per position:\n   ";
print_array(\@divergencies);





