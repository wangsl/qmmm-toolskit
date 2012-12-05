#!/bin/env perl

use warnings;
use strict;

$| = 1;

my @key_words = ("TINKER-PARAMETER", "ZPDB-INPUT", "ZPDB-OUTPUT", "TINKER-XYZ",
		 "BOND-SEARCH-EXE", "QM-ATOMS", "PSEUDO-ATOMS", "ZERO-CHARGE-ATOMS");

# data from input

my $tinker_parameter_file = undef;
my $zpdb_input = undef;
my $bond_construct_exe = undef;

my $zpdb_output = "qmmm.zpdb";
my $tinker_xyz = "qmmm.xyz";

my @qm_atoms = ();
my @pseudo_atoms = ();
my @zero_charge_atoms = ();

# global variables

my @atom_type_parameters = ();
my @zpdb = ();
my @qmmm_zpdb = ();

my @qm_list = ();
my @pseudo_list = ();
my @zero_charge_list = ();
my @mm_list = ();

# data for 
my @atom_name = ();
my @qmmm_xyz = ();
my @tinker_type = ();
my @atomic_number = ();
my @atom_bond_number = ();

my $atom_bond_list = [];

# real work
{
    my ($input) = @ARGV;
    die "$0: input_file\n" unless defined $input;
    
    &read_input_file($input);
    
    @qm_atoms = sort {$a <=> $b} &uniq(@qm_atoms);
    @pseudo_atoms = sort {$a <=> $b} &uniq(@pseudo_atoms);
    @zero_charge_atoms = sort {$a <=> $b} &uniq(@zero_charge_atoms);
    
    &print_input_data();
    
    &read_zpdb();
    
    &generate_qmmm_atom_group_list();
    &reorder_zpdb_to_get_qmmm_zpdb();
    &write_qmmm_zpdb_to_file();
    
    &read_tinker_parameter_file();
    
    &generate_data_for_file_conversion();
    &construct_bond_list();
    &write_tinker_qmmm_xyz_to_file();

    exit;
}


# same radius as in VMD 1.8.7 PeriodicTable.C

sub radius {
    my ($i_ele) = @_;
    die "$0: element index should > 0" unless ($i_ele > 0);
    if($i_ele == 1) {
	return 1.00;
    } elsif ($i_ele == 6) {
	return 1.50;
    } elsif ($i_ele == 7) {
	return 1.40;
    } elsif ($i_ele == 8) {
	return 1.30;
    } elsif ($i_ele == 9) {	
	return 1.20;
    } elsif ($i_ele == 15) {
	return 1.50;
    } elsif ($i_ele == 16) {
	return 1.90;
    } else {
	return 1.50;
    }
}

sub read_tinker_parameter_file {
    print " Read Tinker force field parameter from file: $tinker_parameter_file\n";
    open(FIN, "<$tinker_parameter_file") || die "$0: Failed to open $tinker_parameter_file\n";
    while(<FIN>) {
	next if(/^s*$/);
	my @tmp = &mysplit($_);
	next unless (uc($tmp[0]) eq "ATOM");
	my $itype = $tmp[1];
	my $iele = $tmp[-3];
	my $nbonds = $tmp[-1];
	$atom_type_parameters[$itype][0] = $iele;
	$atom_type_parameters[$itype][1] = $nbonds;
    }
    close(FIN);
}

sub mysplit {
    my ($input) = @_;
    my @output = ();
    foreach (split(/\s+/, $input)) {
        if(/\S+/) {
            push(@output, $_);
        }
    }
    return @output;
}

sub read_zpdb {
    print " Read zpdb from file: $zpdb_input\n";
    @zpdb = ();
    open(PDBIN, "<$zpdb_input") || die "$0: Failed to open $zpdb_input\n";
    while(<PDBIN>) {
	chomp;
	my $line = $_;
	next if(substr($line, 0, 4) ne "ATOM" && substr($line, 0, 6) ne "HETATM");
	push(@zpdb, $line);
    }
    close(PDBIN);
}

sub uniq {
    my (@a) = @_;
    my %hash = map { $_ => 1 } @a;
    return keys %hash;
}

sub generate_qmmm_atom_group_list {
    print " Generate QMMM atom list\n";
    @qm_list = ();
    @pseudo_list = ();
    @zero_charge_list = ();
    @mm_list = ();
    for(my $i = 0; $i < @zpdb; $i++) {
	my $atom_serial_number = substr($zpdb[$i], 6, 5);
	if(&is_qm_atom($atom_serial_number)) {
	    push(@qm_list, $i);
	} elsif (&is_pseudo_atom($atom_serial_number)) {
	    push(@pseudo_list, $i);
	} elsif (&is_zero_charge_atom($atom_serial_number)) {
	    push(@zero_charge_list, $i);
	} else {
	    push(@mm_list, $i);
	}
    }
}

sub is_qm_atom {
    my ($atom_index) = @_;
    foreach (@qm_atoms) {
	return 1 if ($_ == $atom_index);
    }
    return 0;
}

sub is_pseudo_atom {
    my ($atom_index) = @_;
    foreach (@pseudo_atoms) {
	return 1 if ($_ == $atom_index);
    }
    return 0;
}

sub is_zero_charge_atom {
    my ($atom_index) = @_;
    foreach (@zero_charge_atoms) {
	return 1 if ($_ == $atom_index);
    }
    return 0;
}

sub reorder_zpdb_to_get_qmmm_zpdb {
    print " Reorder zpdb for QMMM\n";
    @qmmm_zpdb = ();
    foreach(@qm_list, @pseudo_list, @zero_charge_list, @mm_list) {
	my $i = $_;
	push(@qmmm_zpdb, $zpdb[$i]);
    }
}

sub write_qmmm_zpdb_to_file {
    print " QMMM zpdb will be written to file $zpdb_output\n";
    open(ZPDBOUT, ">$zpdb_output") || die "$0: Failed to open $zpdb_output\n";
    print ZPDBOUT (" QMMM-atom-list: ", @qm_list+0, "  ", @pseudo_list+0, "  ",
		   @zero_charge_list+0, "  ", @mm_list+0, "\n");
    print ZPDBOUT join("\n", @qmmm_zpdb), "\n";
    close(ZPDBOUT);
}

sub generate_data_for_file_conversion {
    print " Generate data for zpdb to Tinker xyz conversion\n";
    @atom_name = ();
    @qmmm_xyz = ();
    @tinker_type = ();
    @atomic_number = ();
    @atom_bond_number = ();
    
    for(my $i = 0; $i < @qmmm_zpdb; $i++) {
	my $line = $qmmm_zpdb[$i];
	$atom_name[$i] = substr($line, 12, 4);
	$qmmm_xyz[$i][0] = substr($line, 30, 8);
	$qmmm_xyz[$i][1] = substr($line, 38, 8);
	$qmmm_xyz[$i][2] = substr($line, 46, 8);
	my $itype = substr($line, 54, 6);
	$tinker_type[$i] = $itype;
	$atomic_number[$i] = $atom_type_parameters[$itype][0];
	$atom_bond_number[$i] = $atom_type_parameters[$itype][1];
    }
}

sub generate_bond_list_with_perl {
    print " Generate bond list\n";
    $atom_bond_list = [];
    my $n_atoms = @qmmm_zpdb;
    for(my $i = 0; $i < $n_atoms; $i++) {
	print " $i";
	print "\n" if(($i+1)%10 == 0);
	my $atomic_number_i = $atomic_number[$i];
	my $rad_i = &radius($atomic_number_i);
	my $coord_i = $qmmm_xyz[$i];
	for(my $j = $i+1; $j < $n_atoms; $j++) {
	    my $atomic_number_j = $atomic_number[$j];
	    next if ($atomic_number_i == 1 && $atomic_number_j == 1);
	    my $rad_j = &radius($atomic_number_j);
	    my $coord_j = $qmmm_xyz[$j];
	    my $cutoff = 0.60*($rad_i + $rad_j);
	    my $rij2 = &dist2($coord_i, $coord_j);
	    if($rij2 < $cutoff*$cutoff) {
		push(@{$atom_bond_list->[$i]}, $j+1);
		push(@{$atom_bond_list->[$j]}, $i+1);
	    }
	}
    }
    print "\n";
    &check_bond_number();
}

sub check_bond_number {
    print " Check bond number\n";
    for(my $i = 0; $i < @{$atom_bond_list}; $i++) {
	if(defined $atom_bond_list->[$i]) {
	    @{$atom_bond_list->[$i]} = sort {$a <=> $b} &uniq(@{$atom_bond_list->[$i]});
	    if(@{$atom_bond_list->[$i]} != $atom_bond_number[$i]) {
		print " Bond number error for: ", $i+1, "\n";
		print $qmmm_zpdb[$i], "\n";
		print " Found ", @{$atom_bond_list->[$i]}+0, " bonds";
		print " Suppose to be ", $atom_bond_number[$i], "\n";
		die "$0: bond number error";
	    }
	}
    }
}

sub dist2 {
    my ($r1, $r2) = @_;
    return 
	($r1->[0] - $r2->[0])**2 + 
	($r1->[1] - $r2->[1])**2 + 
	($r1->[2] - $r2->[2])**2;
}

sub write_tinker_qmmm_xyz_to_file {
    print " Write tinker xyz to file $tinker_xyz\n";
    open(FOUT, ">$tinker_xyz") || die "$0: Failed to open $tinker_xyz";
    
    my $n_atoms = @qmmm_xyz;
    my $n_qm_atoms = @qm_list;
    my $n_pseudo_atoms = @pseudo_list;
    my $n_zero_charge_atoms = @zero_charge_list;
    my $n_mm_atoms = @mm_list;
    
    die "$0: QMMM atom group number error\n" unless ($n_atoms == $n_qm_atoms + $n_pseudo_atoms +
						     $n_zero_charge_atoms + $n_mm_atoms);
    
    printf FOUT ("%6i%6i%6i%6i%6i\n", $n_atoms, $n_qm_atoms, $n_pseudo_atoms,
		 $n_zero_charge_atoms, $n_mm_atoms);
    for(my $i = 0; $i < $n_atoms; $i++) {
	printf FOUT ("%6i  %4s %8.3f %8.3f %8.3f %6i ", $i+1, $atom_name[$i],
		     $qmmm_xyz[$i][0], $qmmm_xyz[$i][1], $qmmm_xyz[$i][2], $tinker_type[$i]);
	if(defined $atom_bond_list->[$i]) {
	    for(my $j = 0; $j < @{$atom_bond_list->[$i]}; $j++) {
		printf FOUT "%6i", $atom_bond_list->[$i]->[$j];
	    }
	}
	printf FOUT "\n";
    }
    close(FOUT);
}

sub generate_bond_list_with_exe {
    print " Generate bond list with executable file: $bond_construct_exe\n";
    my $tmp_input = "bond_list-$$.input";
    my $tmp_output = "bond_list-$$.output";
    open(FOUT, ">$tmp_input") || die "$0: Failed to open $tmp_input\n";
    my $n_atoms = @qmmm_xyz;
    printf FOUT "%6i\n", $n_atoms;
    for(my $i = 0; $i < $n_atoms; $i++) {
	my $i_ele = $atomic_number[$i];
	printf FOUT ("%6i%4i %8.3f %8.3f %8.3f  %6.3f\n", $i+1, $i_ele,
		     $qmmm_xyz[$i][0], $qmmm_xyz[$i][1], $qmmm_xyz[$i][2], 
		     &radius($i_ele));
    }
    close(FOUT);
    
    die "$0: File $tmp_input does not exist\n" unless (-e $tmp_input);
    unlink $tmp_output || die "$0: Failed to unlink $tmp_output\n" if (-e $tmp_output);
    my $cmd = "$bond_construct_exe $tmp_input $tmp_output";
    system "$cmd" || die "$0: Failed to run '$cmd'\n";
    unlink $tmp_input || die "$0: Failed to remove file $tmp_input\n" if (-e $tmp_input);
    die "$0: File $tmp_output does not exist\n" unless (-e $tmp_output);

    print " Reade bond list from file $tmp_output\n";
    $atom_bond_list = [];

    open(FIN, "<$tmp_output") || die "$0: Failed to open $tmp_output\n";
    die "$0: first line error" unless defined ($_ = <FIN>);
    my @tmp = split;
    die "$0: atomic number error\n" unless $tmp[0] == @qmmm_xyz;
    for(my $i = 0; $i < @qmmm_xyz; $i++) {
	die "$0: file reading error" unless defined ($_ = <FIN>);
	my @tmp = split;
	die "$0: line index error\n" unless ($tmp[0] == $i+1);
	if(@tmp > 1) {
	    push(@{$atom_bond_list->[$i]}, @tmp[1..$#tmp]);
	}
    }
    close(FIN);

    unlink $tmp_output || die "$0: Failed to unlink $tmp_output\n" if (-e $tmp_output);

    &check_bond_number();
}

sub construct_bond_list {
    if(defined $bond_construct_exe) {
	&generate_bond_list_with_exe();
    } else {
	&generate_bond_list_with_perl();
    }
}

sub read_input_file {
    my ($input_file) = @_;
    print " Read input data from file: ", $input_file, "\n";
    open(FIN, "<$input_file") || die "$0: Failed to open $input_file\n";
    while(<FIN>) {
	next if (/\s*\#/ || /^\s*$/);
	my @tmp = split;
	my $key_word = uc($tmp[0]);
	if($key_word eq "TINKER-PARAMETER") {
	    $tinker_parameter_file = $tmp[1];
	} elsif ($key_word eq "ZPDB-INPUT") {
	    $zpdb_input = $tmp[1];
	} elsif ($key_word eq "ZPDB-OUTPUT") {
	    $zpdb_output = $tmp[1];
	} elsif ($key_word eq "TINKER-XYZ") {
	    $tinker_xyz = $tmp[1];
	} elsif ($key_word eq "BOND-SEARCH-EXE") {
	    $bond_construct_exe = $tmp[1];
	} elsif ($key_word eq "QM-ATOMS") {
	    push(@qm_atoms, @tmp[1..$#tmp]);
	} elsif ($key_word eq "PSEUDO-ATOMS") {
	    push(@pseudo_atoms, @tmp[1..$#tmp]);
	} elsif ($key_word eq "ZERO-CHARGE-ATOMS") {
	    push(@zero_charge_atoms, @tmp[1..$#tmp]);
	} else {
	    print " Valid key words: ", join(", ", @key_words), "\n";
	    die "$0: key word '$key_word' error\n";
	}
    }
    close(FIN);
    
    die "$0: TINKER-PARAMETER undefined\n" unless defined  $tinker_parameter_file;
    die "$0: ZPDB-INPUT undefined\n" unless defined $zpdb_input;
}

sub print_input_data {
    print "\n";
    print " Tinker parameter file: ", $tinker_parameter_file, "\n";
    print " ZPDB input file: ", $zpdb_input, "\n";
    print " QMMM ZPDB output file: ", $zpdb_output, "\n";
    print " QMMM Tinker xyz file: ",  $tinker_xyz, "\n";
    print " QM atoms: ", join(", ", @qm_atoms), "\n";
    print " Pseudo atoms: ", join(", ", @pseudo_atoms), "\n";
    print " Zero charge atoms: ", join(", ", @zero_charge_atoms), "\n";
    print " Bond search executable file: ";
    print $bond_construct_exe if(defined $bond_construct_exe);
    print "\n\n";
}
