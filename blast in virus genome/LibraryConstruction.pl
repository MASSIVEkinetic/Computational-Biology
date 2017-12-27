#!/usr/bin/perl
use strict;
use warnings;
use Storable;

#E-value is not required.

my %DB;
my $ALLSEQ = "";
my $ALLSEQLEN;
#my $SpeciesName;
my $WordLength = 11;
my $counter;
my $criteria = 0;

open DATABASE, "/1_disk/public_resources/hg19.fa";
#open MYDB, ">LIBRARY";

my $line;
while($line=<DATABASE>){
    if($line=~/>chr22/){
        $criteria = 0;
        last;
    }
    if($criteria == 1){
        $line =~ s/\s//g;
        $line =~ tr/[a-z]/[A-Z]/;
        $ALLSEQ = $ALLSEQ.$line;
    }
    if($line=~/>chr21/){
        $criteria = 1;
    }
}
=pod
my $first = shift @line;
if($first=~/>/){
    $SpeciesName = "";
}
=cut

#This block is not thought to be efficient
=pod
for my $elements (@line){
    $elements =~ s/\s//g;
}
=cut

$ALLSEQLEN = length($ALLSEQ);
print $ALLSEQLEN."\n";

my $words;
for($counter=0;$counter<($ALLSEQLEN+1-$WordLength);$counter++){
    $words = substr($ALLSEQ,$counter,$WordLength);
    if($DB{$words}){
        $DB{$words} = $DB{$words}."+$counter";
    }else{
        $DB{$words} = "$counter";
    }
}
print "completed!\n";
store \%DB,"LIBRARY";