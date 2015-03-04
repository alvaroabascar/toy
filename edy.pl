#!/usr/bin/env perl

use strict;
use warnings;

use HMM;

HMM::select_model("toy");

HMM::predict("CTTCATGTGAAAGCAGACGTAAGTCA");
