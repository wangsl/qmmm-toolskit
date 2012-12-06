#!/usr/bin/perl -w

use strict;

my @exes = qw(alchemy  analyze  anneal  archive  correlate  crystal 
	      diffuse  distgeom  document  dynamic  gda  intedit  
	      intxyz  minimize  minirot  minrigid  monte  newton  newtrot 
	      nucleic  optimize  optirot  optrigid  path  pdbxyz  polarize
	      prmedit  protein  pss  pssrigid  pssrot  radial  saddle  scan
	      sniffer  spacefill  spectrum  superpose  sybylxyz  testgrad 
	      testhess  testlight  testrot  timer  timerot  vibrate  vibrot
	      xtalfit  xtalmin  xyzedit  xyzint  xyzpdb  xyzsybyl
	      );

foreach (@exes) {
    print "\$(E)/$_ : \$(O)/$_.o \$(TINKERLIB) \$(SIMLIB)\n";
    print "\t\$(LINK) \$(LINKFLAGS) -o \$(E)/$_ \$(O)/$_.o \$(LIBS)\n\n";
}

