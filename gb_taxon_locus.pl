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

my $count;
for my $gb (@gbs){
  my ($locus, $taxon);

  if ($gb =~ /LOCUS\s+(\w+)\s+/){
    $locus = $1;
  }else{
    warn "Can't get locus\n";
    next;
  }

  if ($gb =~ /db_xref\=\"taxon\:(\d+)/){
    $taxon = $1;
  }else{
    warn "Can't get taxon\n";
    next;
  }

  printf "%s,%s\n", $locus, $taxon;
  $count++;
}

printf STDERR "There are [%d] valid sequences.\n", $count;
__END__

=head1 NAME
gb_taxon_locus.pl -scan a multi-sequence .gb file

=head1 SYNOPSIS

  wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/organelle/plastid/plastid.1.genomic.gbff.gz
  gzip -d plastid.1.genomic.gbff.gz
  perl gb_taxon_locus.pl plastid.1.genomic.gbff.gz > refseq_id_seq.csv
  
=cut
