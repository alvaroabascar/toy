#!/usr/bin/env perl

use strict;
use warnings;

package MODELS;

# This module contains all the models implemented along this project.
# Each model is composed by a set of states, a set of symbols, and a
# set of transition {emp} and emission {trp} probabilities.
#
# At the end of this file you will find a dictionary %models which
# gives to each model a representative name. New models must be
# added to this dictionary using a new name.

##################################################################
#
# toy model, considering intron, exon and only the G of the donor
#
##################################################################
our %hmm_toy = ( states    => ["begin", "exon", "Gdonor", "intron", "end"],
            symbols => ["A", "C", "G", "T"],
            emp       => [],
            trp       => [],
            initstate =>  "begin"  );

$hmm_toy{emp} = [ [0.00, 0.00, 0.00, 0.00], # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],     # exon nuc are equiprobable
              [0.05, 0.00, 0.95, 0.00],     # donor sites are mostly G
              [0.40, 0.10, 0.10, 0.40],     # introns are A-T rich
              [0.0, 0.0, 0.0, 0.0] ];       # end state

$hmm_toy{trp} = [ [0.0, 1.0, 0.0, 0.0, 0.0],   # outgoing from begin
                  [0.0, 0.9, 0.1, 0.0, 0.0],   # outgoing from exon
                  [0.0, 0.0, 0.0, 1.0, 0.0],   # outgoing from donor
                  [0.0, 0.0, 0.0, 0.9, 0.1], # outgoing from intron
                  [0.0, 0.0, 0.0, 0.0, 0.0] ]; # end state

##################################################################
#
# model considering the most relevant positions of the U1 snRNP
# binding region (positions 3, 4, 5, 6 of real_human_donors)
#
##################################################################
our %hmm_U1 = ( states    => ["begin", "exon", "1donor", "Gdonor", "3donor", "4donor","5donor","6donor", "intron", "end"],
                symbols => ["A", "C", "G", "T"],
                emp       => [],
                trp       => [],
                initstate =>  "begin"  );

$hmm_U1{emp} = [ [0.00, 0.00, 0.00, 0.00],          # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],             # exon nuc are equiprobable
              [0.09456, 0.03039, 0.80668, 0.06837], # donor1 sites are mostly G
              [0.00, 0.00, 1, 0.00],                # this site is ALWAYS G
              [0.00, 0.00, 0.00, 1],                # this site is ALWAYS T
              [0.58916, 0.02529, 0.36006, 0.02549], # donor4
              [0.70532, 0.07147, 0.11455, 0.10866], # donor5 sites are mostly A
              [0.08357, 0.05418, 0.79268, 0.06957], # donor6 sites are mostly G
              [0.25, 0.25, 0.25, 0.25],             # introns are A-T rich
              [0.0, 0.0, 0.0, 0.0] ];               # end state

$hmm_U1{trp} = [ 
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from begin
              [0.0, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from exon
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor1
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donorG
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor3
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],   # outgoing from donor4
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],   # outgoing from donor5
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],   # outgoing from donor6
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1], # outgoing from intron
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ]; # end state

######################################################################
#
# model considering ALL the positions of the U1 snRNP binding region
#
#####################################################################
our %hmm_U1_all = ( states    => ["begin", "exon", "1donor","2donor","3donor", "Gdonor", "5donor", "6donor","7donor","8donor","9donor", "intron", "end"],
            symbols => ["A", "C", "G", "T"],
            emp       => [],
            trp       => [],
            initstate =>  "begin"  );

$hmm_U1_all{emp} = [ [0.00, 0.00, 0.00, 0.00],      # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],             # exon nuc are equiprobable
              [0.33906, 0.36066, 0.18862, 0.11166], # donor1 sites are mostly G
              [0.62795, 0.11555, 0.11635, 0.14014], # donor2 sites are mostly G
              [0.09456, 0.03039, 0.80668, 0.06837], # donor3 sites are mostly G
              [0.00, 0.00, 1, 0.00],                # this site is ALWAYS G
              [0.00, 0.00, 0.00, 1],                # this site is ALWAYS T
              [0.58916, 0.02529, 0.36006, 0.02549], # donor6 
              [0.70532, 0.07147, 0.11455, 0.10866], # donor7 sites are mostly A
              [0.08357, 0.05418, 0.79268, 0.06957], # donor8 sites are mostly G
              [0.17283, 0.15424, 0.19172, 0.48121], # donor9 sites are mostly G
              [0.25, 0.25, 0.25, 0.25],             # introns are A-T rich
              [0.0, 0.0, 0.0, 0.0] ];               # end

$hmm_U1_all{trp} = [ 
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from begin
              [0.0, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from exon
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor1
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donorG
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor3
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor4
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor5
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],       # outgoing from donor6
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],       # outgoing from donor7
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],       # outgoing from donor8
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1],     # outgoing from intron
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ];     # end


##################################################################
#
# model considering ALL the positions of the U1 snRNP binding region
# and the splice regulator
#
##################################################################

our %hmm_TIA = ( states    => ["begin", "exon", "1donor","2donor","3donor", "Gdonor", "5donor", "6donor","7donor","8donor","9donor", "stia", "intron", "end"],
            symbols => ["A", "C", "G", "T"],
            emp       => [],
            trp       => [],
            initstate =>  "begin"  );

$hmm_TIA{emp} = [ [0.00, 0.00, 0.00, 0.00],           # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],               # exon nuc are equiprobable
              [0.33906, 0.36066, 0.18862, 0.11166],   # donor1 sites are mostly G
              [0.62795, 0.11555, 0.11635, 0.14014],   # donor2 sites are mostly G
              [0.09456, 0.03039, 0.80668, 0.06837],   # donor3 sites are mostly G
              [0.00, 0.00, 1, 0.00],                  # this site is ALWAYS G
              [0.00, 0.00, 0.00, 1],                  # this site is ALWAYS T
              [0.58916, 0.02529, 0.36006, 0.02549],   # donor6 
              [0.70532, 0.07147, 0.11455, 0.10866],   # donor7 sites are mostly A
              [0.08357, 0.05418, 0.79268, 0.06957],   # donor8 sites are mostly G
              [0.17283, 0.15424, 0.19172, 0.48121],   # donor9 sites are mostly G
              [0.06, 0.06, 0.06, 0.82],               # tia emits mostly T
              [0.25, 0.25, 0.25, 0.25],               # introns
              [0, 0, 0, 0] ];                         # end state

$hmm_TIA{trp} = [
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from begin
              [0.0, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from exon
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor1
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor2
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor3
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor5
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor6
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor7
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],       # outgoing from donor8
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0],       # outgoing from donor9
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0],       # outgoing from tia
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.79, 0.1],     # outgoing from intron
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ];     # end


##################################################################
#
# model considering ALL the positions of the U1 snRNP binding region
# and the splice regulator
#
##################################################################

our %hmm_TIA2 = (states => ["begin", "exon", "1donor","2donor","3donor", "Gdonor", "5donor", "6donor","7donor","8donor","9donor", "stia", "stia2", "intron", "end"],
            symbols => ["A", "C", "G", "T"],
            emp       => [],
            trp       => [],
            initstate =>  "begin"  );

$hmm_TIA2{emp} = [ [0.00, 0.00, 0.00, 0.00],   # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],   # exon nuc are equiprobable
              [0.33906, 0.36066, 0.18862, 0.11166],   # donor1 sites are mostly G
              [0.62795, 0.11555, 0.11635, 0.14014],   # donor2 sites are mostly G
              [0.09456, 0.03039, 0.80668, 0.06837],   # donor3 sites are mostly G
              [0.00, 0.00, 1, 0.00],   # donor sites are mostly G---> NO PSEUDOCOUNTS!
              [0.00, 0.00, 0.00, 1],   # donor5 sites are mostly T
              [0.58916, 0.02529, 0.36006, 0.02549],   # donor6 
              [0.70532, 0.07147, 0.11455, 0.10866],   # donor7 sites are mostly A
              [0.08357, 0.05418, 0.79268, 0.06957],   # donor8 sites are mostly G
              [0.17283, 0.15424, 0.19172, 0.48121],   # donor9 sites are mostly G
              [0.06, 0.06, 0.06, 0.82],   # tia emits mostly T
              [0.06, 0.06, 0.06, 0.82],   # tia emits mostly T
              [0.25, 0.25, 0.25, 0.25],   # introns are A-T rich
              [0, 0, 0, 0] ];

$hmm_TIA2{trp} = [ 
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from begin
              [0.0, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from exon
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor1
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor2
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor3
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor5
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor6
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor7
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],       # outgoing from donor8
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0],       # donor9 can go to tia or intron
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0],       # outgoing from tia
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0],       # outgoing from tia
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.89, 0.01],     # intron can go to tia, itself or end
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ];     # outgoing from end

##################################################################
#
# model considering ALL the positions of the U1 snRNP binding region
# and the splice regulator tia 3
#
##################################################################

our %hmm_TIA3 = ( states    => ["begin", "exon", "1donor","2donor","3donor", "Gdonor", "5donor", "6donor","7donor","8donor","9donor", "intron", "stia", "end"],
                  symbols => ["A", "C", "G", "T"],
                  emp       => [],
                  trp       => [],
                  initstate =>  "begin"  );

$hmm_TIA3{emp} = [ [0.00, 0.00, 0.00, 0.00],        # begin is a silent state
              [0.25, 0.25, 0.25, 0.25],             # exon nuc are equiprobable
              [0.33906, 0.36066, 0.18862, 0.11166], # donor1 sites are mostly G
              [0.62795, 0.11555, 0.11635, 0.14014], # donor2 sites are mostly G
              [0.09456, 0.03039, 0.80668, 0.06837], # donor3 sites are mostly G
              [0.00, 0.00, 1, 0.00],                # this site is ALWAYS G
              [0.00, 0.00, 0.00, 1],                # this site is ALWAYS T
              [0.58916, 0.02529, 0.36006, 0.02549], # donor6
              [0.70532, 0.07147, 0.11455, 0.10866], # donor7 sites are mostly A
              [0.08357, 0.05418, 0.79268, 0.06957], # donor8 sites are mostly G
              [0.17283, 0.15424, 0.19172, 0.48121], # donor9 sites are mostly G
              [0.25, 0.25, 0.25, 0.25],             # introns are A-T rich
              [0.06, 0.06, 0.06, 0.82],             # tia emits mostly T
              [0, 0, 0, 0] ];                       # end state

$hmm_TIA3{trp} = [ 
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from begin
              [0.0, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from exon
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor1
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor2
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor3
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donorG
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor5
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor6
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],   # outgoing from donor7
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],   # outgoing from donor8
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0],   # outgoing from donor9
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.89, 0.1, 0.01], # outgoing from intron
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.89, 0.01], # outgoing from tia
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ]; # outgoing from end

# A dictionary linking each model names to the model itself
our %models = (toy => \%hmm_toy, U1 => \%hmm_U1, U1_all => \%hmm_U1_all, TIA => \%hmm_TIA, TIA2 => \%hmm_TIA2, TIA3 => \%hmm_TIA3);

1;
