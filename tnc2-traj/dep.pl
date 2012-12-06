#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
# $Id: dep.pl,v 1.1 2007/09/02 21:11:37 sw77 Exp $

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
