#!/usr/bin/perl -w
require 5.002;
# $Id: io,v 1.15 2002/04/03 00:38:19 hstern Exp $

$usage = <<EOF;
usage: mkio [-o][-D<name>] <.h files>
-o print to standard output instead of io.C
-D define name

EOF

$opt_o = 0;
die $usage unless @ARGV;
while (defined($_ = shift)) {
  if (/\-D(.*)/) {
    $def{$1} = 1;
  } elsif (/\-o/) {
    $opt_o = 1
  } elsif (/\.h$/) {
    push(@f,$_);
  }
}

$io = '\/\/ (in|out|io)';
$io2 = '\/\* (in|out|io) \*\/';

foreach (@f) {
  die "mkio: not a .h file" unless /(.*)\.h/;
  $fname = $1 . "io.C";
  open(F, $_) || die "mkio: cannot open $_";
  unless ($opt_o) {
    open(G, ">$fname") || die "mkio: cannot open $fname";
    select(G);
  }
  print <<EOF;
#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include \"indent.h\"
#include \"$_\"
#include \"die.h\"
#include \"out.h\"

EOF
  undef $super_read;
  undef $super_write;
  undef $show_self;
  undef @ifdef;
 loop: while(<F>) {
  if (/^\#ifdef (\S+)/) {
	 push(@ifdef,$1);
    } elsif (/^\#endif/) {
      pop(@ifdef);
    } elsif (/^\#define (\S+)/) {
      $def{$1} = 1;
    } elsif (/^struct (\S+)/) {
      $class = $1;
      $super_read = $1 if (/public (\S+)\s+\/\/ (in|io)/);
      $super_write = $1 if (/public (\S+)\s+\/\/ (out|io)/);
    } elsif (/^class (\S+)/) {
      $class = $1;
      $super_read = $1 if (/public (\S+)\s+\/\/ (in|io)/);
      $super_write = $1 if (/public (\S+)\s+\/\/ (out|io)/);
    } elsif (/$io/ || /$io2/) {
      foreach $t (@ifdef) {
	next loop unless $def{$t};
      }
      my($which) = $1;
      s/\bconst\b//;
      s/\bstatic\b//;
      s/\bvirtual\b//;
      s/\<.*\>//;
      s/[\&\<\>\;]//g;
      die "io: cannot parse io line: $_"
	unless /^\s*(\S+)\s+(\S.*)$io/ || /^\s*(\S+)\s+(\S.*)$io2/;
      my($x) = $1;
      die unless defined $x;
      $_ = $2;
      my(@nm);
      if (/\(/) {
	  push(@nm,$_);
      } else {
	  @nm = split(',',$_);
      }
      $show_self = $1 if /^show_self$/;
      foreach (@nm) {
        push(@type, $x);
	push(@fields, $_);
	push(@which, $which);
      }
    } elsif (/^\}\;/) {
      undef $anyin;
      undef $anyout;
      foreach (@which) {
	die unless (/in|io|out/);
	$anyin = 1 if (/in|io/);
	$anyout = 1 if (/out|io/);
      }
      $anyin = 1 if $super_read;
      $anyout = 1 if $super_write;
      &do_write if $anyout && (@fields || $super_write || $super_read);
      &do_read if $anyin && (@fields || $super_write || $super_read);
      undef @fields;
      undef @which;
      undef @type;
      undef $super_read;
      undef $super_write;
      undef $show_self;
    } elsif (/void show_self\(\)/) {
      $show_self = 1;
    }
  }
}

sub do_write {
    #cout << *this << "\\n" << flush;
  print <<"EOF" if $show_self;
void $class\:\:show_self() 
{ 
  Out() << *this << "\\n" << flush;
}

EOF
  print <<"EOF";
ostream & operator<<(ostream &s, const $class &c)
{
  s << \"{\\n\"\;
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << \"}\\n\";
}

void $class\:\:write_fields(ostream &s) const
{
EOF
print "  $super_write\:\:write_fields(s);\n" if defined $super_write;
$i = -1;
foreach (@fields) {
  $i++;
  local($pointer,$key);
  $key = $_;
  $key =~ s/\s//g;
  next if $which[$i] =~ /in/;
  next if /\(/;
  if ($key =~ /\*/) {
    $pointer = 1;
    $key =~ s/\*//;
  }
  if (defined ($pointer)) {
      print "  if ($key)\n";
      if($type[$i] eq "char") {
	  print "    s << Indent() << \"$key \" << $key << \"\\n\";\n";
      } else {
	  print "    s << Indent() << \"$key \" << *$key << \"\\n\";\n";
      }
	  
  } elsif ($type[$i] eq "Str") {
    print "  if (strlen($key) > 0)\n";
    print "    s << Indent() << \"$key \" << $key << \"\\n\";\n";
  } else {
    print "  s << Indent() << \"$key \" << $key << \"\\n\";\n";
  }
}
print "}\n\n";
}

sub do_read {
  print <<EOF;
istream & operator>>(istream &s, $class &c)
{
  char buf[256], bracket;
  s >> bracket;
  if (bracket != '{')
    die("error reading class '$class': expected '{', found something else\\n");
  while (s) {
    s >> bracket;
    if (bracket == '}') {
      return s;
    } else {
      s.putback(bracket);
    }
    s >> buf;
    if (c.read_field(s, buf) == 0) {
      Out() << "error: unknown field name: '" << buf << "'\\n"
               "while reading class '$class'.\\n"
               "Must be one of:\\n";
      c.write_field_names();
      die("");
    }
  }
  die("EOF or error reading class '$class'\\n");
  return s;
}

int $class\:\:read_field(istream &s, const char *buf)
{
  if (!s)
    die("EOF or error reading class '$class'\\n");
EOF
  print <<EOF if $show_self && $anyout;
  else if (!strcmp(buf,"show_self"))
    show_self();
EOF
  undef $i;
  undef @readfields;
  foreach (@fields) {
    next if $which[$i++] =~ /out/;
    $nm = $_;
    $nm =~ s/\(.*//g;
    $nm =~ s/[\*\s]//g;
    push(@readfields,"\'$nm\'");
    print "  else if (!strcmp(buf,\"$nm\"))\n";
    if (/\(([^\)]*)\)/) {
	my($args) = $1;
	if ($args =~ /istream/) {
	    s/\(.*//;
	    s/\s//g;
	    print "    $_(s);\n";
	} elsif ($args =~ /\S/) {
	    @args = split(",",$args);
	    print "  {\n";
	    my(@v);
	    foreach (@args) {
		my($t,$v) = /(\S+)\s+(\S+)/;
		print "    $t $v;\n";
		print "    s >> $v;\n";
		push(@v,$v);
	    }
	    print "    if (!s)\n";
	    print "      die(\"error reading arguments to $nm, class $class\\n\"\n";
	    print "          \"Expecting $args\");\n";
	    print "    $nm(", join(",",@v),");\n";
	    print "  }\n";
	} else {
	    s/\(.*//;
	    s/\s//g;
	    print "    $_();\n";
	}
    } elsif (/\*/) {
	print "    s >> *$nm;\n";
    } else {
	s/\s//g;
	print "    s >> $_;\n";
    }
  }
  $fldlst = join(", ", @readfields);
  print "  else\n";
  if (defined $super_read) {
    print "    return $super_read\:\:read_field(s, buf);\n";
  } else {
    print "    return 0;\n";
  }
  print "  return 1;\n";
  print "}\n\n";
  print "void $class\:\:write_field_names() const\n";
  print "{\n";
  print "  Out() << \"$fldlst\\n\";\n";
  print "  $super_read\:\:write_field_names();\n" if defined $super_read;
  print "}\n\n";
}
