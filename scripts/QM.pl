#!/usr/bin/perl -w

use strict;

my @elements = qw(LI BE NE MG AL SI CL AR FE CO NI CU ZN);

for my $xyz(@ARGV) {
    &read_tinker_xyz($xyz);
}

exit;

sub read_tinker_xyz
{
    my ($tinker_xyz) = @_;
    die "$0: argument error\n" unless defined $tinker_xyz;
    die "$0: $tinker_xyz file does not exist\n" unless (-e $tinker_xyz);
    
    open(FIN, "<$tinker_xyz") || die "$0: Failed to open $tinker_xyz\n";
    
    for(my $n_mol = 0; ; $n_mol++) {
	last unless defined ($_ = <FIN>);
	my ($n_atoms, $n_qm_atoms, $n_pseudo_atoms) = split;
	last unless (defined $n_atoms && defined $n_qm_atoms &&
		     defined $n_pseudo_atoms);
	
	printf " %d\n", $n_qm_atoms+$n_pseudo_atoms;
	printf " %s, index = %d\n", $tinker_xyz, $n_mol+1;
	
	for(my $i = 0; $i < $n_qm_atoms+$n_pseudo_atoms; $i++) {
	    my $line = $i+2;
	    die "$0: $tinker_xyz $line errorn\n" unless defined ($_ = <FIN>);
	    my ($index, $sym_tinker, $x, $y, $z) = split;
	    my $j = $i+1;
	    die "$0: atomic index error: $j\n" unless ($index == $j);
	    my $sym = undef;
	    if($sym_tinker =~ /([A-Za-z]+)/) {
		$sym = $1;
	    } else {
		die "$0: element symbol error\n";
	    }
	    
	    my $e = undef;
	    if(length $sym >= 2) {
		my $OK = undef;
		my $E2 = uc substr($sym, 0, 2);
		foreach (@elements) {
		    if($E2 eq $_) {
			$OK = 1;
			last;
		    }
		}
		if($OK) {
		    $e = ucfirst substr($sym, 0, 2);
		} else {
		    $e = ucfirst substr($sym, 0, 1);
		}
	    } else {
		$e = ucfirst substr($sym, 0, 1);
	    }
	    die "$0:\n:" unless (defined $e);
	    printf "  %-2s%12.6f%12.6f%12.6f\n", $e, $x, $y, $z;
	}
	
	for(my $i = $n_qm_atoms+$n_pseudo_atoms; $i < $n_atoms; $i++) {
	    die "$0: $tinker_xyz line errorn\n" unless defined ($_ = <FIN>);
	}
    }
    close(FIN);
}

