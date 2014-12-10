#!/usr/bin/perl
use warnings;
use strict;

open(my $treeFile, "<", "mtDNA_haplogroups.csv");
my %header = ();
my %lookup = ();
my $maxArgs = 0;
while(<$treeFile>){
    chomp;
    $_ =~ tr/\"//d;
    my @F = split(/,/, $_);
    if(!%header){
        $maxArgs = scalar(@F);
        for(my $i = 0; $i < scalar(@F); $i++){
            $header{$F[$i]} = $i;
        }
    } else {
        my $mutString = $F[$header{"normalised"}];
        $mutString = "" unless defined($mutString);
        my @mutArr = ($mutString)?split(/ /, $mutString):();
        $lookup{"mutations"}{$F[$header{"name"}]} = \@mutArr;
        if($F[$header{"parent"}] ne "NA"){
            $lookup{"parent"}{$F[$header{"name"}]} = $F[$header{"parent"}];
        }
        $lookup{"level"}{$F[$header{"name"}]} = $F[$header{"level"}];
    }
}
close($treeFile);

sub getLoc($){
    my $sig = shift;
    $sig =~ s/;.*$//;
    return($sig);
}

sub getAllMutations($){
    my $name = shift;
    if(!defined($lookup{"mutations"}{$name})){
        return(());
    }
    if(defined($lookup{"allMutations"}{$name})){
        return(@{$lookup{"allMutations"}{$name}});
    }
    my @out = @{$lookup{"mutations"}{$name}};
    if(defined($lookup{"parent"}{$name})){
        @out = (&getAllMutations($lookup{"parent"}{$name}),@out);
    }
    my %addedSNPs = ();
    my %addedINDELs = ();
    while(@out){
        my $sig = shift(@out);
        my ($loc, $ref, $alt) = split(/;/, $sig);
        if(length($ref) == length($alt)){
            if($addedSNPs{$loc}){
                my $sigP = $addedSNPs{$loc};
                my ($locP, $refP, $altP) = split(/;/, $sigP);
                $sig = ($refP eq $alt)?"":sprintf("%d;%s;%s",$loc,$refP,$alt);
            }
            $addedSNPs{$loc} = $sig;
        } else {
            if($addedINDELs{$loc}){
                my $sigP = $addedINDELs{$loc};
                my ($locP, $refP, $altP) = split(/;/, $sigP);
                my $lr = length($refP.$ref);
                my $la = length($altP.$alt);
                my $ldiff = $la - $lr;
                $ref = ($ldiff < 0)?("X" x -$ldiff):""; # deletion
                $alt = ($ldiff > 0)?("X" x $ldiff):""; # insertion
                $sig = ($ldiff==0)?"":sprintf("%d;%s;%s",$loc,$ref,$alt);
            }
            $addedINDELs{$loc} = $sig;
        }
    }
    @out = grep {$_} (values(%addedSNPs),values(%addedINDELs));
    @out = sort {getLoc($a) <=> getLoc($b)} @out;
    $lookup{"allMutations"}{$name} = \@out;
    return(@out);
}

sub normaliseVCF($){
    my $fileName = shift;
    my @allSigs = ();
    open(my $vcfFile, "<", $fileName);
    while(<$vcfFile>){
        if(substr($_,0,1) ne "#"){
            chomp;
            my @F = split(/\t/, $_);
            my ($loc, $ref, $alt) = @F[(1,3,4)];
            if($alt =~ /,/){
                foreach my $tAlt (split(/,/, $alt)){
                    push(@allSigs, sprintf("%d;%s;%s",$loc,$ref,$tAlt));
                }
            } else {
                push(@allSigs, sprintf("%d;%s;%s",$loc,$ref,$alt));
            }
        }
    }
    close($vcfFile);
    my @out = ();
    while(@allSigs){
        my $sig = shift(@allSigs);
        my ($loc, $ref, $alt) = split(/;/, $sig);
        if(length($ref) == length($alt)){
            push(@out, $sig); # SNP, nothing more to do
        } else { # INDEL
            my $lr = length($ref);
            my $la = length($alt);
            my $ldiff = $la - $lr;
            my $lminadj = (($la<$lr)?$la:$lr) + 1;
            $ref = ($ldiff < 0)?("X" x -$ldiff):""; # deletion
            $alt = ($ldiff > 0)?("X" x $ldiff):""; # insertion
            for(my $i = 0; $i < $lminadj; $i++){
                push(@out, sprintf("%d;%s;%s",$loc+$i,$ref,$alt));
            }
        }
    }
    @out = sort {getLoc($a) <=> getLoc($b)} @out;
    return(@out);
}

sub getProportion($$){
    my ($typeGroupRef, $testGroupRef) = @_;
    my @typeGroup = @{$typeGroupRef};
    my @testGroup = @{$testGroupRef};
    my %typeSet = ();
    foreach (@typeGroup){
        $typeSet{$_} = 1;
    }
    my $matchCount = 0;
    foreach (@testGroup){
        if($typeSet{$_}){
            $matchCount++;
            # stop double counting for duplicated variants in the testGroup
            delete $typeSet{$_};
        }
    }
    if($matchCount == scalar(@typeGroup)){
        return(1);
    } else {
        return($matchCount / scalar(@typeGroup));
    }
}

sub getMissing($$){
    my ($typeGroupRef, $testGroupRef) = @_;
    my @typeGroup = @{$typeGroupRef};
    my @testGroup = @{$testGroupRef};
    my %testSet = ();
    foreach (@testGroup){
        $testSet{$_} = 1;
    }
    my @missing = ();
    foreach (@typeGroup){
        if(!$testSet{$_}){
            push(@missing, $_);
        }
    }
    return(@missing);
}

my $maxNameLength = 0;
my %vcfMuts = ();
foreach my $vcfFileName (@ARGV){
    if(length($vcfFileName) > $maxNameLength){
	$maxNameLength = length($vcfFileName);
    }
    my @normalisedMuts = normaliseVCF($vcfFileName);
    $vcfMuts{$vcfFileName} = \@normalisedMuts;
}

my %levelMaxes = ();
foreach my $name (keys($lookup{"mutations"})){
    my @groupMuts = getAllMutations($name);
    foreach my $vcfFileName (keys(%vcfMuts)){
        my @testMuts = @{$vcfMuts{$vcfFileName}};
        my $prop = getProportion(\@groupMuts, \@testMuts);
        my $level = $lookup{"level"}{$name};
        if(!$levelMaxes{$vcfFileName}{$level}{"prop"} || 
           ($levelMaxes{$vcfFileName}{$level}{"prop"} < $prop)){
            $levelMaxes{$vcfFileName}{$level}{"prop"} = $prop;
            $levelMaxes{$vcfFileName}{$level}{"names"} = ();
        }
        if($levelMaxes{$vcfFileName}{$level}{"prop"} == $prop){
            push(@{$levelMaxes{$vcfFileName}{$level}{"names"}}, $name);
        }
    }
}

printf("%s,%s,%s,%s,%s\n", "Name","Level","Match Proportion","Best Haplotype","Other Information");
foreach my $vcfFileName (sort(keys(%vcfMuts))){
    my %vcfMaxes = %{$levelMaxes{$vcfFileName}};
    foreach my $level (sort {$b <=> $a} keys(%vcfMaxes)){
        my $prop = $levelMaxes{$vcfFileName}{$level}{"prop"};
        my @names = @{$levelMaxes{$vcfFileName}{$level}{"names"}};
        my $firstName = shift(@names);
        my $otherString = "";
        if(($prop >= 0.9) && ($prop < 1) && (!@names)){
            my @groupMuts = getAllMutations($firstName);
            my @testMuts = @{$vcfMuts{$vcfFileName}};
            $otherString = 
                sprintf(" [missing: %s]", 
                        join(" ",getMissing(\@groupMuts, \@testMuts)));
        }
        if(@names){
            $otherString = sprintf(" [also: %s]",
                                   join(",", @names));
        }
	$otherString =~ s/^ /,/;
        printf("%s,%d,%0.4f,%s%s\n", $vcfFileName,
               $level, $levelMaxes{$vcfFileName}{$level}{"prop"}, 
               $firstName, $otherString);
    }
}

