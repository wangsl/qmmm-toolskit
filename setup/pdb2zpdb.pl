#!/bin/env perl 

use warnings;
use strict;

$| = 1;

my @key_words = ("PDB-INPUT", "ZPDB-OUTPUT", "AMBER-TINKRT-TABLE", 
		 "N-TERMINAL-RESID", "C-TERMINAL-RESID");

my ($input_file) = @ARGV;
die "$0: input-file\n" unless defined $input_file;

#input data

my $pdb_input = undef;
my $zpdb_output = undef;
my @N_terminal = ();
my @C_terminal = ();
my $amber_tinker_file = undef;

# global variables

my %parameters = ();
my @pdb = ();
my @zpdb = ();
my $ZPDBOUT = undef;

# main part
{
    &read_input_file($input_file);
    @N_terminal = sort {$a <=> $b} &uniq(@N_terminal);
    @C_terminal = sort {$a <=> $b} &uniq(@C_terminal);
    
    &print_input_data();
    
    &create_parameters_table();
    &read_input_pdb();
    
    &assgin_N_terminal_residues();
    &assgin_C_terminal_residues();
    &assign_non_terminal_residues();
    
    &open_zpdb_out();
    &print_zpdb();
    &close_zpdb_out();
    
    exit;
}

sub create_parameters_table {
    open(PIN, "<$amber_tinker_file") || die "$0: Failed to open $amber_tinker_file\n";
    while(<PIN>) {
	chomp;
	next if(/^\s*#/ || /^\s*$/);
	my @tmp = split;
	my $key = uc("$tmp[0]-$tmp[2]");
	die "'$_' exist\n" if(exists $parameters{$key});
	$parameters{$key} = $tmp[1];
    }
    close(PIN);
}

sub print_parameter_table {
    print " Amber-Tinker force field parameters table\n";
    foreach my $key (sort keys %parameters) {
	printf("%-10s%10i\n", $key, $parameters{$key});
    }
}

sub read_input_pdb {
    print " Read pdb structure from file $pdb_input\n";
    @pdb = ();
    open(PDBIN, "<$pdb_input") || die "$0: Failed to open $pdb_input\n";
    while(<PDBIN>) {
	chomp;
	my $line = $_;
        next if(substr($line, 0, 4) ne "ATOM" && substr($line, 0, 6) ne "HETATM");
	push(@pdb, $line);
    }
    close(PDBIN);
}

sub assgin_N_terminal_residues {
    print " Assign Tinker parameter for N-terminal residues\n";

    for(my $i = 0; $i < @pdb; $i++) {
	my $resid = substr($pdb[$i], 22, 4);
	next if(! &is_N_terminal($resid));
	my $resname = uc(substr($pdb[$i], 17, 3));
	my $atom_name = undef;
	if(substr($pdb[$i], 12, 4) =~ /\s*(\S+)\s*/) {
	    $atom_name = $1;
	} else {
	    die "$0: atom name error: $atom_name\n";
	}
	my $key = uc("$resname-$atom_name");
	my $N_key = uc("N$resname-$atom_name");
	my $line = $pdb[$i];
	my $tinker_type = undef;
	if(exists $parameters{$N_key}) {
	    $tinker_type = $parameters{$N_key};
	} elsif (exists $parameters{$key}) {
	    $tinker_type = $parameters{$key};
	} else {
	    die "$0: $key and $N_key does not exist\n";
	}

	substr($line, 54, 6) = sprintf("  %4i", $tinker_type);
	$zpdb[$i] = substr($line, 0, 60);
    }
}

sub is_N_terminal {
    my ($resid) = @_;
    foreach my $id(@N_terminal) {
	return 1 if($resid == $id);
    }
    return 0;
}

sub assgin_C_terminal_residues {
    print " Assign Tinker parameter for C-terminal residues\n";

    for(my $i = 0; $i < @pdb; $i++) {
	my $resid = substr($pdb[$i], 22, 4);
	next if(! &is_C_terminal($resid));
	my $resname = uc(substr($pdb[$i], 17, 3));
	my $atom_name = undef;
	if(substr($pdb[$i], 12, 4) =~ /\s*(\S+)\s*/) {
	    $atom_name = $1;
	} else {
	    die "$0: atom name error: $atom_name\n";
	}
	my $key = uc("$resname-$atom_name");
	my $C_key = uc("C$resname-$atom_name");
	my $line = $pdb[$i];
	my $tinker_type = undef;
	if(exists $parameters{$C_key}) {
	    $tinker_type = $parameters{$C_key};
	} elsif (exists $parameters{$key}) {
	    $tinker_type = $parameters{$key};
	} else {
	    die "$0: $key and $C_key does not exist\n";
	}
	
	substr($line, 54, 6) = sprintf("  %4i", $tinker_type);
	$zpdb[$i] = substr($line, 0, 60);
    }
}

sub is_C_terminal {
    my ($resid) = @_;
    foreach my $id(@C_terminal) {
	return 1 if($resid == $id);
    }
    return 0;
}

sub assign_non_terminal_residues {
    print " Assign Tinker parameter for non-terminal residues\n";
    for(my $i = 0; $i < @pdb; $i++) {
	next if(defined $zpdb[$i]);
	my $resid = substr($pdb[$i], 22, 4);
	my $resname = uc(substr($pdb[$i], 17, 3));
	my $atom_name = undef;
	if(substr($pdb[$i], 12, 4) =~ /\s*(\S+)\s*/) {
	    $atom_name = $1;
	} else {
	    die "$0: atom name error: $atom_name\n";
	}
	my $key = uc("$resname-$atom_name");
	my $line = $pdb[$i];
	my $tinker_type = undef;
	if(exists $parameters{$key}) {
	    $tinker_type = $parameters{$key};
	} else {
	    die "$0: $key does not exist\n";
	}
	
	substr($line, 54, 6) = sprintf("  %4i", $tinker_type);
	$zpdb[$i] = substr($line, 0, 60);
    }
}

sub print_zpdb {
    my $i = 1;
    foreach(@zpdb) {
	print $ZPDBOUT $_, "\n";
	$i++;
    }
}

sub read_input_file {
    my ($input_file) = @_;
    open(FIN, "<$input_file") || die "$0: Failed to open input file '$input_file'\n";
    while(<FIN>) {
	next if(/\s*\#/ || /^\s*$/);
	my @tmp = split;
	my $key_word = uc($tmp[0]);
	if($key_word eq "PDB-INPUT") {
	    $pdb_input = $tmp[1];
	} elsif ($key_word eq "ZPDB-OUTPUT") {
	    $zpdb_output = $tmp[1];
	} elsif ($key_word eq "AMBER-TINKER-TABLE") {
	    $amber_tinker_file = $tmp[1];
	} elsif ($key_word eq "N-TERMINAL-RESID") {
	    push(@N_terminal, @tmp[1..$#tmp]);
	} elsif ($key_word eq "C-TERMINAL-RESID") {
	    push(@C_terminal, @tmp[1..$#tmp]);
	} else {
	    print "Valid key words: ", join(", ", @key_words), "\n";
	    die "$0: keyword $key_word error\n";
	}
    }
    close(FIN);
}

sub print_input_data {
    print " PDB input file: ", $pdb_input, "\n";
    print " ZPDB output file: ", $zpdb_output, "\n";
    print " Amber-Tinker force field parameter file: ", $amber_tinker_file, "\n";
    print " N-terminal residue id: ", join(", ", @N_terminal), "\n";
    print " C-terminal residue id: ", join(", ", @C_terminal), "\n";
}

sub uniq {
    my (@a) = @_;
    my %hash = map { $_ => 1 } @a;
    return keys %hash;
}

sub open_zpdb_out {
    open(ZPDBOUT, ">$zpdb_output") || die "$0: Failed to open $zpdb_output\n";
    $ZPDBOUT = \*ZPDBOUT;
}

sub close_zpdb_out {
    close($ZPDBOUT) || die "$0: Failed to close ZPDBOUT\n" if(defined $ZPDBOUT);
}
