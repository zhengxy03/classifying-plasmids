#!/user/bin/perl
use strict;
use warnings;
use autodie;

use lib "/home/zxy0303/perl5/lib/perl5/";
use Path::Tidy;

my $file = shift;

if( !$file ){
  die "you don't provide a .gbff file.\n";
}elsif( ! -e $file){
  die "[$file] doesn't exist.\n";
}

$file = path($file);

my $content = $file->slurp;

my @gbs = grep {/\S+/} split(/^\/\//m, $content);

printf STDERR "There are [%d] sequences.\n", scalar @gbs;

  
