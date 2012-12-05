#!/bin/env perl 

# $Id: cut-sphere.pl 6 2012-12-05 04:48:45Z wangsl2001@gmail.com $ 

use warnings;
use strict;
use Getopt::Std;

# add Tue Dec  4 23:42:58 EST 2012
# this script can work with any number of water molecule systems now

# This script can not be used if the water molecule resid number is greater than 9999
# for this case, you have to modify the subroutine generate_water_index,
# which is used to determine the begin and end index of water molecules

# global variables

my ($input_file) = @ARGV;

my @key_words = ("PDB-INPUT", "PDB-OUTPUT", "CENTER-COORDINATE",
		 "SPHERE-RADIUS", "BOX-DIMENSION", "IONS-SAVE",
		 "WATER-NAME", "IONS-NAME");

my $pdb_input_file = undef;
my $pdb_output_file = undef;
my @center_coordinates = undef;
my $cutoff = 27.0; 
my @box_dimensions = (2000.0, 2000.0, 2000.0); 
my @ions_save_name = ();

my @water_name = ("WAT", "HOH", "CRY", "SOL");
my @ions_name = ("NA+", "Cl-");

### You don't need to modify the following part

my @protein = ();
my @water = ();
my @ion = ();

my @water_begin_index = ();
my @water_end_index = ();
my @water_save = ();

my @ion_save = ();

my $cutoff2 = undef;

my $PDBOUT = undef;

my $n_ion_save = 0;

my $usage = <<EOF;
 usage: $0;
EOF

# main part    
{    
    $| = 1;

    &read_input($input_file);

    &upcase(\@water_name);
    @water_name = sort &uniq(@water_name);
    
    &upcase(\@ions_name);
    @ions_name = sort &uniq(@ions_name);

    &upcase(\@ions_save_name);
    @ions_save_name = sort &uniq(@ions_save_name);

    &print_input_data();

    # cut sphere real work
    
    $cutoff2 = $cutoff * $cutoff;
    
    &read_pdb();
    &generate_water_index();
    &determine_water_inside_sphere();
    &determine_ions_inside_sphere();
    
    print "\n";

    &open_pdb_out();
    &write_protein();
    &write_ion();
    &write_water();
    &close_pdb_out();
    
    print "\n Job finished\n";
}

exit;

sub read_pdb {
    print " Reading PDB from file $pdb_input_file\n";
    open(PDBIN, "<$pdb_input_file") || die "$0: Failed to open file $pdb_input_file\n";
    while(<PDBIN>) {
	my $line = $_;
	next if(substr($line, 0, 4) ne "ATOM" && substr($line, 0, 6) ne "HETATM");
	
	my $resname = uc(substr($line, 17, 3));
	my $is_water_ = &is_water($resname);
	my $is_ion_ = 0;
	if(!$is_water_) { $is_ion_ = &is_ion($resname); }
	
	if($is_water_) {
	    push(@water, $line);
	} elsif ($is_ion_) {
	    push(@ion, $line);
	} else {
	    push(@protein, $line);
	}
    }
    
    close(PDBIN);
    
    print " There are ", @protein+0, " protein atoms\n";
    print " There are ", @water/3, " water molecules\n";
    print " There are ", @ion+0, " ions\n";
}

sub generate_water_index {
    print " Generate water index\n";

    my @begin = ();
    
    my $iwater = -1;
    for(my $i = 0; $i < @water; $i++) {
	my $id = substr($water[$i], 22, 4);
	$begin[$id] = 0 if(!$begin[$id]);
	if($begin[$id]%3 == 0) {
	    $iwater++;
	    $water_begin_index[$iwater] = $i;
	    $water_end_index[$iwater] = $i;
	} else {
	    $water_end_index[$iwater] = $i;
	}
	
	$begin[$id]++;
    }

    die "$0: water number error\n" unless (@water/3 == @water_begin_index && 
					   @water/3 == @water_end_index);
}

sub determine_water_inside_sphere {
    print " Determine water inside sphere\n";
    
    my @water_xyz = ();
    for(my $i = 0; $i < @water; $i++) {
	$water_xyz[$i][0] = substr($water[$i], 30, 8) - $center_coordinates[0];
	$water_xyz[$i][1] = substr($water[$i], 38, 8) - $center_coordinates[1];
	$water_xyz[$i][2] = substr($water[$i], 46, 8) - $center_coordinates[2];
	$water_save[$i] = 0;
    }

    my $n_water_save = 0;

    for(my $i = 0; $i < @water_begin_index; $i++) {
	my $ibegin = $water_begin_index[$i];
	my $iend = $water_end_index[$i];

	my ($ic, $jc, $kc) = (0, 0, 0);
	my $inside = 0;
	for(my $iwater = $ibegin; $iwater <= $iend; $iwater++) {
	    ($ic, $jc, $kc, $inside) = &inside_sphere_with_neighbour_cell(@{$water_xyz[$iwater]});
	    last if($inside);
	}
	
	if($inside) {
	    $n_water_save++;
	    for(my $iwater = $ibegin; $iwater <= $iend; $iwater++) {
		my $x = substr($water[$iwater], 30, 8) + $ic*$box_dimensions[0];
		my $y = substr($water[$iwater], 38, 8) + $jc*$box_dimensions[1];
		my $z = substr($water[$iwater], 46, 8) + $kc*$box_dimensions[2];
		substr($water[$iwater], 30, 8) = sprintf("%8.3f", $x);
		substr($water[$iwater], 38, 8) = sprintf("%8.3f", $y);
		substr($water[$iwater], 46, 8) = sprintf("%8.3f", $z);
		$water_save[$iwater] = 1;
	    }
	}
    }
    
    print " There are ", $n_water_save, " water molecules will be kept\n";
}

sub inside_sphere_with_neighbour_cell {
    my (@r) = @_;
    my $r2 = $r[0]*$r[0] + $r[1]*$r[1] + $r[2]*$r[2];
    if($r2 <= $cutoff2) {
	return (0, 0, 0, 1);
    }
    
    for(my $i = -1; $i <= 1; $i++) {
	my $x = $r[0] + $i*$box_dimensions[0];
	for(my $j = -1; $j <= 1; $j++) {
	    my $y = $r[1] + $j*$box_dimensions[1];
	    for(my $k = -1; $k <= 1; $k++) {
		my $z = $r[2] + $k*$box_dimensions[2];
		$r2 = $x*$x + $y*$y + $z*$z;
		if($r2 <= $cutoff2) {
		    return ($i, $j, $k, 1);
		}
	    }
	}
    }

    return (0, 0, 0, 0);
}

sub write_protein {
    print " Write protein\n" if(@protein);
    for(my $i = 0; $i < @protein; $i++) {
	print $PDBOUT $protein[$i];
    }
}

sub write_water {
    print " Write water\n" if(@water);
    for(my $i = 0; $i < @water; $i++) {
	if($water_save[$i]) {
	    print $PDBOUT $water[$i];
	}
    }
}

sub write_ion {
    print " Write ions\n" if($n_ion_save);
    for(my $i = 0; $i < @ion; $i++) {
	my $line = $ion[$i];
	print $PDBOUT $line if($ion_save[$i]);
	#my $resname = uc(substr($line, 17, 3));
	#if(&is_save_ion($resname)) {
	#    print $PDBOUT $line;
	#}
    }
}

sub is_water {
    my ($resname) = @_;
    foreach(@water_name) {
	return 1 if($_ eq $resname);
    }
    return 0;
}

sub is_ion {
    my ($resname) = @_;
    foreach(@ions_name) {
	return 1 if($_ eq $resname);
    }
    return 0;
}

sub is_save_ion {
    return 0 if(!@ions_save_name);
    my ($resname) = @_;
    foreach(@ions_save_name) {
	return 1 if($_ eq $resname);
    }
    return 0;
}

sub open_pdb_out {
    open(PDBOUT, ">$pdb_output_file") || die "$0: Failed to open $pdb_output_file\n";
    $PDBOUT = \*PDBOUT;
}

sub close_pdb_out {
    close ($PDBOUT) || die "$0: Failed to close PDBOUT\n" if(defined $PDBOUT);
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

sub read_input {
    my ($input_file) = @_;
    open(FIN, "<$input_file") || die "$0: Failed to open $input_file\n";
    while(<FIN>) {
	next if(/\s*\#/ || /^s*$/);
	my @tmp = split;
	my $key_word = uc($tmp[0]);
	if($key_word eq "PDB-INPUT") {
	    $pdb_input_file = $tmp[1];
	} elsif ($key_word eq "PDB-OUTPUT") {
	    $pdb_output_file = $tmp[1];
	} elsif ($key_word eq "CENTER-COORDINATE") {
	    $center_coordinates[0] = $tmp[1];
	    $center_coordinates[1] = $tmp[2];
	    $center_coordinates[2] = $tmp[3];
	} elsif ($key_word eq "SPHERE-RADIUS") {
	    $cutoff = $tmp[1];
	} elsif ($key_word eq "BOX-DIMENSION") {
	    $box_dimensions[0] = $tmp[1];
	    $box_dimensions[1] = $tmp[2];
	    $box_dimensions[2] = $tmp[3];
	} elsif ($key_word eq "IONS-NAME") {
	    push(@ions_name, @tmp[1..$#tmp]);
	} elsif ($key_word eq "WATER-NAME") {
	    push(@water_name, @tmp[1..$#tmp]);
	} elsif ($key_word eq "IONS-SAVE") {
	    push(@ions_save_name, @tmp[1..$#tmp]);
	} else {
	    print "Valid key words: ", join(", ", @key_words), "\n";
	    die "$0: keyword '$key_word' error\n";
	}
    }

    close(FIN);
}

sub uniq {
    my (@a) = @_;
    my %hash = map { $_ => 1 } @a;
    return keys %hash;
}

sub print_input_data {
    print " PDB input file: ", $pdb_input_file, "\n";
    print " PDB output file: ", $pdb_output_file, "\n";
    print " Center coordinate: ", join(", ", @center_coordinates), "\n";
    print " Sphere radius: ", $cutoff, "\n";
    print " Box dimensions: ", join(", ", @box_dimensions), "\n";
    print " Water residue: ", join(", ", @water_name), "\n";
    print " Ion residue: ", join(", ", @ions_name), "\n";
    print " Ion to be saved: ", join(", ", @ions_save_name), "\n";
    print "\n";
}

sub upcase {
    my ($a) = @_;
    for(my $i = 0; $i < @{$a}; $i++) {
	$a->[$i] = uc($a->[$i]);
    }
}

sub determine_ions_inside_sphere {
    print " Determine ion inside sphere\n";
    my @ion_xyz = ();
    for(my $i = 0; $i < @ion; $i++) {
	$ion_xyz[$i][0] = substr($ion[$i], 30, 8) - $center_coordinates[0];
	$ion_xyz[$i][1] = substr($ion[$i], 38, 8) - $center_coordinates[1];
	$ion_xyz[$i][2] = substr($ion[$i], 46, 8) - $center_coordinates[2];
	$ion_save[$i] = 0;
    }
    
    $n_ion_save = 0;
    for(my $i = 0; $i < @ion; $i++) {
	my $r2 = 
	    $ion_xyz[$i][0]*$ion_xyz[$i][0] + 
	    $ion_xyz[$i][1]*$ion_xyz[$i][1] +
	    $ion_xyz[$i][2]*$ion_xyz[$i][2];
	if($r2 <= $cutoff2) {
	    $ion_save[$i] = 1;
	    $n_ion_save++;
	}
    }
    
    print " There are ", $n_ion_save, " ions inside sphere\n";
}
