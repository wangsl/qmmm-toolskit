#!/usr/bin/perl -w
 
#use strict;
use Getopt::Std;

my $usage = <<EOF;
usage: pathenergy [-f] [-b] [-l <rc_min>] [-r <rc_max>] [-n <n_digit>] [-d <drc>]
EOF

$opt_f = 0;
$opt_b = 0;
$opt_d = 0.01;
$opt_l = -20.0;
$opt_r = 20.0;
$opt_n = 2;

die $usage unless getopts("fbd:l:r:n:");

die "$0: can not assign -f -b together\n" if($opt_f && $opt_b);

=head1
print $opt_d, "\n";
print $opt_l, "  ", $opt_r, "\n";
print $opt_n, "\n";
=cut

my $rc_min = $opt_l;
my $rc_max = $opt_r;
my $n_digit = $opt_n;
my $drc = $opt_d;

my $base_name = "RCF";
$base_name = "RCB" if($opt_b);

my $format = "%.$n_digit"."f";

my @potent = [];
my $ii = 0;
for(my $rc = $rc_min; $rc <= $rc_max; $rc += $drc) {
    my $index = sprintf($format, $rc);
    my $fname = "$base_name.$index";
    my $i = 1;
    my $file_name = "$fname-$i.log";
    next if (! -e $file_name);
    for(; ; ) {
	$i++;
	my $file_name_next = "$fname-$i.log";
	if(! -e $file_name_next) {
	    last;
	} else {
	    $file_name = $file_name_next;
	}
    }
    my $cmd = "grep \"Optimization step = 1,\" $file_name";
    open(FIN, "$cmd |") || die "$0: Failed to $cmd\n";
    my $line = undef;
    my $j = 0;
    while(<FIN>) {
	$line = $_;
	$j++;
    }
    close(FIN);
    next if ($j < 2);
    #die "$0: error\n" unless defined $line;
    next unless defined $line;
    $_ = $line;
    if(/Optimization step = 1, energy = (\S+) Hartree/) {
	$potent[$ii][0] = $index;
	$potent[$ii][1] = $1;
	$ii++;
    } else {
	die "$0: $_ error\n";
    }
}

my $len = @potent + 0;

exit if (! $len);

my $e0 = $potent[0][1];
my $Hartree = 627.5095;

for(my $i = 0; $i < $len; $i++) {
    printf("%6.3f%20.10f\n", $potent[$i][0], ($potent[$i][1]-$e0)*$Hartree);
}
