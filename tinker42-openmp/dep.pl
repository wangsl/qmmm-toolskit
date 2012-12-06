#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
# $Id: dep.pl 9 2008-08-13 15:25:52Z wangsl $

$a = '';

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

#$allfiles = join " ", grep !/rcsid.C/, glob( "*.[hcCfF]" );
#print "rcsid.C: $allfiles\n";
#print "\t perl rcsid.pl > rcsid.C\n";
