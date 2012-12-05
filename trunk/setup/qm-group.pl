#!/usr/bin/perl -w

use strict;

my @qm = (11..19, 22..44);
my @pseudo = (8, 45);
my @zero_charge = (5, 6, 7, 9, 10, 46, 47, 48);

my @pdb = ();

my $qm_pdb_file = "SAM-LYH.zpdb";
open(FIN, "<$qm_pdb_file") || die "$0: Failed to open $qm_pdb_file\n";

while(<FIN>) {
    push(@pdb, $_);
}

close(FIN);

my @qm_list = ();
for my $q (@qm) {
    print substr($pdb[$q-1], 7, 5), "\n";
}

