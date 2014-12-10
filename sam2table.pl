#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through);

my $seqFileName = "";
my $refseq = "";

my $showRef = 0;

GetOptions('reffasta|f=s' => \$seqFileName,
          'showref!' => \$showRef);

if(!$seqFileName){
  print(STDERR "Error: reference FASTA file is required (-f <fname>)\n");
  exit(1);
}

open(my $seqFile, "<", $seqFileName);
while(<$seqFile>){
  chomp;
  if(/^[^>]/){
    $refseq .= $_;
  }
}
close($seqFile);

my %delPoss = ();
my %insPoss = ();
my %poss = ();

my %sampleIDs = ();
my %sampVariants = ();

my $processed = 0;

print(STDERR "Processing input file (assuming SAM format)...");

while(<>){
  chomp;
  my @F = split();
  my $id = $F[0];
  my $pos = $F[3];
  my $cigar = $F[5];
  my $altseq = $F[9];
  if(exists($sampleIDs{$id})){ # don't process things that have been seen before
    if($sampleIDs{$id} == 1){
      print(STDERR "D");
    }
    $sampleIDs{$id} = -1;
    next;
  }
  if($pos != 1){
    print(STDERR "X"); # ignore incomplete matches
    $sampleIDs{$id} = -2;
    next;
  }
  $sampleIDs{$id} = 1;
  my $refpos = 0;
  my $altpos = 0;
  while($cigar =~ s/^([0-9]+)(.)//){
    my $length = $1;
    my $type = $2;
    if($type eq "M"){
      for(my $i = 0; $i < $length; $i++){
        my $rss = substr($refseq,$refpos+$i,1);
        my $ass = substr($altseq,$altpos+$i,1);
        if(($rss ne $ass) && ($ass ne "N") && ($rss ne "N")){
          $sampVariants{$refpos+$i}{$id} = $ass;
        }
      }
      $refpos += $length;
      $altpos += $length;
    }
    if($type eq "I"){
      my $sseq = substr($altseq,$altpos,$length);
      if($sseq !~ /N/){ # inserted Ns are not informative
        if((!$insPoss{$refpos}) || ($insPoss{$refpos} < $length)){
          $insPoss{$refpos} = $length;
          $poss{$refpos+0.5} = 1;
        }
        $sampVariants{$refpos+0.5}{$id} = $sseq;
      }
      $altpos += $length;
    }
    if($type eq "D"){
      for(my $i = 0; $i < $length; $i++){
        if(substr($refseq, $refpos + $i, 1) ne "N"){ # exclude positions with N reference
          $delPoss{$refpos+$i} = "-";
          $poss{$refpos+$i} = 1;
          $sampVariants{$refpos+$i}{$id} = "-";
        }
      }
      $refpos += $length;
    }
  } # end while($cigar...)
  if($processed++ % 100 == 0){
    print(STDERR ".");
  }
}
print(STDERR "\n");

my @sPoss = sort {$a <=> $b} (keys(%poss));
my %refGT = ();

printf("%s","sampleID");
foreach my $pos (@sPoss){
  printf("\t%g", $pos);
  if($pos =~ /\.5$/){
    $refGT{$pos} = "-";
  } else {
    $refGT{$pos} = substr($refseq,$pos,1);
  }
}
print("\n");
foreach my $id (sort(keys(%sampleIDs))){
  if($sampleIDs{$id} < 0){ # ignore samples with multiple mappings, or any bad mappings
    next;
  }
  printf("%s",$id);
  foreach my $pos (@sPoss){
    printf("\t%s", $sampVariants{$pos}{$id}?$sampVariants{$pos}{$id}:
             ($showRef?$refGT{$pos}:"."));
  }
  print("\n");
}
