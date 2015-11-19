#!/usr/bin/perl

my $pos = -1;
my $seqName = "";
my $bestFlags = "";
my $bestID = "";
my $bestSeq = "";
my $bestQual = "";

sub printSeq {
    my ($id, $seq, $qual) = @_;
    if($id){
	printf("@%s\n%s\n+\n%s\n", $id, $seq, $qual);
    }
}

while(<>){
    chomp;
    my @F = split(/\t/);
    if(($F[2] ne $seqName) || ($F[3] != $pos) || (length($bestSeq) < length($F[9]))){
	if(($F[2] ne $seqName) || ($F[3] != $pos)){
	    printSeq($bestID, $bestSeq, $bestQual);
	}
	$seqName = $F[2];
	$pos = $F[3];
	$bestID = $F[0];
	$bestFlags = $F[1];
	$bestSeq = $F[9];
	$bestQual = $F[10];
    }
}
printSeq($bestID, $bestSeq, $bestQual);
