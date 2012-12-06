#!/usr/bin/perl -w

use strict;

my @ftran_files = @ARGV;

foreach (@ftran_files) {
    &find_include($_);
}

exit;

sub find_include 
{
    my ($ftran_source) = @_;
    die "$0: please input fortran source code name\n" unless defined $ftran_source;
    die "$0: $ftran_source does not exis\n" unless (-e $ftran_source);
    $_ = $ftran_source;
    my $base_name = undef;
    if(/(\S+).f\b/ || /(\S+).F\b/ || /(\S+).for\b/ || /(\S+).f90\b/) {
	$base_name = $1;
    } else {
	die "$0: '$ftran_source' is not a fotran file, must be *.F *.for *.f90\n";
    }
    
    my %include_i;
    
    my $grep_cmd = "grep -i include $ftran_source";
    open(F, "$grep_cmd |") || die "$0: Failed in $grep_cmd";
    while(<F>) {
	if(/^\s+include\s+\'(\S+).i\'/i) {
	    $include_i{$1} = 1;
	}
    }
    close(F);
    
    my @keys = keys %include_i;

    return if(@keys == 0);
    
    print "\$(O)/$base_name.o:";
    foreach my $key(@keys) {
	print " $key.i";
    }
    print "\n";
}



