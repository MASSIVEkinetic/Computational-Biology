#!/usr/bin/perl
use strict;
use warnings;

my %DB;
my $ALLSEQ = "";
my $ALLSEQLEN;
my $WordLength = 8;
my $counter;

open DATABASE, "lambda_virus.fa";
open MYDB, ">LIBRARY";


my @line = <DATABASE>;

$ALLSEQ = join("",@line);
$ALLSEQ =~ s/\s//g;
$ALLSEQ =~ tr/[a-z]/[A-Z]/;
$ALLSEQLEN = length($ALLSEQ);

my $words;
for($counter=0;$counter<($ALLSEQLEN+1-$WordLength);$counter++){
    $words = substr($ALLSEQ,$counter,$WordLength);
    if($DB{$words}){
        $DB{$words} = $DB{$words}."+$counter";
    }else{
        $DB{$words} = "$counter";
    }
}

my @WordList = sort keys %DB;
for my $word (@WordList){
    print MYDB $word.">".$DB{$word}."\n";
}
