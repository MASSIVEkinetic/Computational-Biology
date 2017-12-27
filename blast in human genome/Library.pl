#!/usr/bin/perl
use strict;
use warnings;
use Storable;

my %DB;
my $ALLSEQ = "";
my $ALLSEQLEN;
my $WordLength = 11;
my $counter;
my $line;
my $words;
my $criteria = 0;
my @ChrName = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"M","X","Y");

for my $NO (@ChrName){
    open DATABASE, "/1_disk/public_resources/hg19.fa";
    #the file path can be modified to fit your own needs
    while($line=<DATABASE>){
        $line =~ s/\s//g;
        if($line=~/>chr(\w*)/){
            if($1 eq $NO){
                $criteria = 1;
                next;
            }else{
                $criteria = 0;
                next;
            }
        }
        if($criteria == 1){
            $line =~ tr/[a-z]/[A-Z]/;
            $ALLSEQ = $ALLSEQ.$line;
        }
    }

    $ALLSEQLEN = length($ALLSEQ);
    print "chr".$NO.": ".$ALLSEQLEN."\n";

    for($counter=0;$counter<($ALLSEQLEN+1-$WordLength);$counter++){
       $words = substr($ALLSEQ,$counter,$WordLength);
        if($words=~/N/g){
            next;
        }
        if(exists $DB{$words}){
            $DB{$words} = $DB{$words}.pack("L",$counter);
        }else{
            $DB{$words} = pack("L",$counter);
        }
    }
    store \$ALLSEQ, "ALLSEQ$NO";
    undef $ALLSEQ;

    if($NO ne "Y"){
        for $words (keys %DB){
            $DB{$words} = $DB{$words}.pack("L",0);
        }
    }
    
    print "chr".$NO." finished!\n";
    close DATABASE;
}

store \%DB,"LIBRARY";