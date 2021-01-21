#!/usr/bin/env perl

# Copyright [2017-2021] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 NAME

 transcripts2fasta.pl

=head1 DESCRIPTION

  A script to dump transcript sequences from an Ensembl core database into a FASTA file.
  It dumps the exon sequences concatenated with flanking sequences of the given length.
  The transcript stable ID is the header for each record.

=head1 OPTIONS
  -dbhost
  -dbport
  -dbuser
  -dbname
  -dnahost
  -dnaport
  -dnauser
  -dnadbname
  -output_file
  -biotype
  -flanking_length

=head1 EXAMPLE

  perl transcripts2fasta.pl \
      --dbhost=MYDBHOST \
      --dbuser=MYREADONLYUSER  \
      --dbname=MYCOREDBNAME \
      --port=MYPORT
      --output_file=tpp_transcripts.fa
      --biotype=transcribed_processed_pseudogene
      --flanking_length=50

=cut

use warnings;
use strict;
use feature 'say';
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $dbname;
my $dbuser;
my $dbhost;
my $dbport;
my $dbpass;

my $dnadbname;
my $dnauser;
my $dnahost;
my $dnaport;
my $dnapass;

my $output_file;

my $biotype = "protein_coding";
my $flanking_length = 0;

my $result = GetOptions ("user|dbuser|u=s" => \$dbuser,
                         "host|dbhost|h=s" => \$dbhost,
                         "port|dbport|P=i" => \$dbport,
                         "dbname|db|D=s"   => \$dbname,
                         "dbpass|pass|p=s" => \$dbpass,
                         "dnauser=s"   => \$dnauser,
                         "dnahost=s"   => \$dnahost,
                         "dnaport=i"   => \$dnaport,
                         "dnadbname=s" => \$dnadbname,
                         "dnadbpass=s"  => \$dnapass,
                         "output_file=s" => \$output_file,
                         "biotype=s" => \$biotype,
                         "flanking_length=i" => \$flanking_length);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dbport,
  -user    => $dbuser,
  -host    => $dbhost,
  -dbname  => $dbname,
  -pass    => $dbpass);

my $dnadb;
if ($dnadbname and $dnauser and $dnahost and $dnaport) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dnaport,
  -user    => $dnauser,
  -host    => $dnahost,
  -dbname  => $dnadbname,
  -pass   => $dnapass);
  $db->dnadb($dnadb);
}

my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

if ($output_file) {
  open(OUT,">".$output_file);
}

foreach my $transcript (@{$transcript_adaptor->fetch_all_by_biotype($biotype)}) {
  
  my $seq_region_start = $transcript->seq_region_start();
  my $seq_region_end = $transcript->seq_region_end();
  my $seq_region_strand = $transcript->seq_region_strand();
  my $seq_region_name = $transcript->seq_region_name();

  my $left_flanking_slice = $slice_adaptor->fetch_by_region('toplevel',
                                                               $seq_region_name,
                                                               $seq_region_start-1-$flanking_length+1,
                                                               $seq_region_start-1,
                                                               $seq_region_strand);
  my $right_flanking_slice = $slice_adaptor->fetch_by_region('toplevel',
                                                               $seq_region_name,
                                                               $seq_region_end+1,
                                                               $seq_region_end+1+$flanking_length-1,
                                                               $seq_region_strand);;

  my $transcript_seq = $transcript->spliced_seq(); # concatenated exon sequences

  if ($seq_region_strand) {
    $transcript_seq = $left_flanking_slice->seq().$transcript_seq.$right_flanking_slice->seq();
  } else {
    my $reverse_left_flanking_seq = Bio::PrimarySeq->new(-seq => $left_flanking_slice->seq())->revcom()->seq();
    my $reverse_right_flanking_seq = Bio::PrimarySeq->new(-seq => $right_flanking_slice->seq())->revcom()->seq();
    $transcript_seq = $reverse_left_flanking_seq.$transcript_seq.$reverse_right_flanking_seq;
  }

  $transcript_seq =~ s/(\w{60})/$1\n/g; # format into 60-base lines

  my $fasta_record_str = $transcript->stable_id_version()."\n".$transcript_seq;
  if ($output_file) {
    say OUT ">".$fasta_record_str;
  } else {
    say ">".$fasta_record_str;
  }
}

if ($output_file) {
  close OUT;
}

;
