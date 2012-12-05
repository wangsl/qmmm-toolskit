#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
# $Id: dep,v 1.3 2000/01/13 19:07:25 hstern Exp $

$a = '$(O)/';

while(<>) {
  s/\/usr\S*\.h//g; 
  s/\w*\.[Cc]//g;
  if (/\\$/) {
    s/\s*\\$//;
    chop;
  }
  s/ +/ /g;
  s/(\S+\.o)/$a$1/g;
  print;
} 

$allfiles = join " ", grep !/rcsid.C/, glob( "*.[hcCf]" );
print "rcsid.C: $allfiles\n";
print "\t perl rcsid.pl > rcsid.C\n";
