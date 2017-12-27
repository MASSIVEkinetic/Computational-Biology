#!/usr/bin/perl
use strict;
use warnings;
use Storable;
use threads ('yield',
             'stack_size' => 64*4096,
             'exit' => 'threads_only',
             'stringify');
use Thread qw(:DEFAULT async yield);

my $startTime = time();

#The defined variables
my $QueryName = $ARGV[0];
my $WordLength = $ARGV[1];
my $threshold = $ARGV[2];

my $gapp = -2;
my $match = 2;
my $mismatch = -1;
#E-value is not required.

my $DB;
my $DBName;
my @ALLSEQ;

my $counter;
my $subcounter;
my $chrcounter;
my $words;
my $Query;
my @temporary;
my $result;
my @results;
my @resultsQ;
my @resultsD;
my @scores;
my @chrs;
my @locations;

my @location01;
my @location02;
my @location03;
my @location04;
my @location05;
my @location06;
my @location07;
my @location08;
my @location09;
my @location10;
my @location11;
my @location12;
my @location13;
my @location14;
my @location15;
my @location16;
my @location17;
my @location18;
my @location19;
my @location20;
my @location21;
my @location22;
my @locationM;
my @locationX;
my @locationY;

#my @KillRepeats;

my %locations = (
    1 => [@location01],
    2 => [@location02],
    3 => [@location03],
    4 => [@location04],
    5 => [@location05],
    6 => [@location06],
    7 => [@location07],
    8 => [@location08],
    9 => [@location09],
    10 => [@location10],
    11 => [@location11],
    12 => [@location12],
    13 => [@location13],
    14 => [@location14],
    15 => [@location15],
    16 => [@location16],
    17 => [@location17],
    18 => [@location18],
    19 => [@location19],
    20 => [@location20],
    21 => [@location21],
    22 => [@location22],
    23 => [@locationM],
    24 => [@locationX],
    25 => [@locationY],
);

my $StartPoint;

my $LOQ;
my $AGOQ;
my $AGOD;
my $AGMax;
my $AGMin;
my $DBPiece;

my @matrix;
my $x;
my $y;
my $top_bottom = 0;
my $up_down = 0;
my $left_right =0;
my $alignedX = " ";
my $alignedY = " ";
my $score;
my $matched;
my $totalAA;

#parameter related to the score filter
#the filter aims to predict the alignment score
#AMax and AMin record number of As

if(not defined $threshold){
    $threshold = 0.8;
    print "No threshold, default threshold = $threshold\n";
}
if(not defined $WordLength){
    $WordLength = 11;
    print "No word length, default length = $WordLength\n";
}

#load the whole sequence
my ($thr1) = threads->create(
    sub{
        $_[0] = \${retrieve("ALLSEQ1")} or die("Cannot open the database!\n");
        $_[1] = \${retrieve("ALLSEQ2")} or die("Cannot open the database!\n");
        $_[2] = \${retrieve("ALLSEQ3")} or die("Cannot open the database!\n");
        $_[3] = \${retrieve("ALLSEQ4")} or die("Cannot open the database!\n");
        return(@_);
    }
);
my ($thr2) = threads->create(
    sub{
        $_[0] =  \${retrieve("ALLSEQ5")} or die("Cannot open the database!\n");
        $_[1] =  \${retrieve("ALLSEQ6")} or die("Cannot open the database!\n");
        $_[2] =  \${retrieve("ALLSEQ7")} or die("Cannot open the database!\n");
        $_[3] =  \${retrieve("ALLSEQ8")} or die("Cannot open the database!\n");
        $_[4] =  \${retrieve("ALLSEQ9")} or die("Cannot open the database!\n");
        return(@_);
    }
);
my ($thr3) = threads->create(
    sub{
        $_[0] =  \${retrieve("ALLSEQ10")} or die("Cannot open the database!\n");
        $_[1] =  \${retrieve("ALLSEQ11")} or die("Cannot open the database!\n");
        $_[2] =  \${retrieve("ALLSEQ12")} or die("Cannot open the database!\n");
        $_[3] =  \${retrieve("ALLSEQ13")} or die("Cannot open the database!\n");
        $_[4] =  \${retrieve("ALLSEQ14")} or die("Cannot open the database!\n");
        $_[5] =  \${retrieve("ALLSEQ15")} or die("Cannot open the database!\n");
        $_[6] =  \${retrieve("ALLSEQ16")} or die("Cannot open the database!\n");
        return(@_);
    }
);
my ($thr4) = threads->create(
    sub{
        $_[0] =  \${retrieve("ALLSEQ17")} or die("Cannot open the database!\n");
        $_[1] =  \${retrieve("ALLSEQ18")} or die("Cannot open the database!\n");
        $_[2] =  \${retrieve("ALLSEQ19")} or die("Cannot open the database!\n");
        $_[3] =  \${retrieve("ALLSEQ20")} or die("Cannot open the database!\n");
        $_[4] =  \${retrieve("ALLSEQ21")} or die("Cannot open the database!\n");
        $_[5] =  \${retrieve("ALLSEQ22")} or die("Cannot open the database!\n");
        $_[6] =  \${retrieve("ALLSEQM")} or die("Cannot open the database!\n");
        $_[7] =  \${retrieve("ALLSEQX")} or die("Cannot open the database!\n");
        $_[8] =  \${retrieve("ALLSEQY")} or die("Cannot open the database!\n");
        return(@_);
    }
);

push(@ALLSEQ,$thr1->join());
push(@ALLSEQ,$thr2->join());
push(@ALLSEQ,$thr3->join());
push(@ALLSEQ,$thr4->join());

#open LOGS,">LOGS";
print "BLAST is ready, at ".(time()-$startTime)."seconds";

#read the query
my $QueryCounter = 0;
open QUERYHANDLE, "$QueryName" or die("Cannot open the query!\n");
while($Query=<QUERYHANDLE>){
    next if $Query =~ />/;
    $QueryCounter++;
    $Query =~ s/\s//g;
    $Query =~ tr/[a-z]/[A-Z]/;
    $LOQ = length($Query);
    $matrix[0][0] = 0;
    for($x=0;$x<$LOQ;$x++){
        $matrix[$x+1][0] = ($x+1)*(-2);
    }
    for($y=0;$y<$LOQ;$y++){
        $matrix[0][$y+1] = ($y+1)*(-2);
    }

    #process the query and do the allocation
    $scores[0] = $threshold;

    for($counter=0;$counter<($LOQ+1-$WordLength);$counter=$counter+8){
        $words = substr($Query,$counter,$WordLength);
        $DBName = substr($words,0,6);
        $DB = \%{retrieve("lib/$DBName")} or die("Cannot open the database!\n");
        if(exists $$DB{$words}){
            @temporary = unpack("L*",$$DB{$words});
            $chrcounter = 1;
            for my $location (@temporary){
                if($location == 0){
                    $chrcounter++;
                }else{
                    $StartPoint = $location-$counter;
                    if(grep(/$StartPoint/,@{$locations{$chrcounter}})){
                        next;
                    }else{
                        $DBPiece = substr(${$ALLSEQ[$chrcounter-1]},($location-$counter),$LOQ);
                        for($y=0;$y<$LOQ;$y++){
                            for($x=0;$x<$LOQ;$x++){
                                if (substr($Query,$y,1) eq substr($DBPiece,$x,1)){
                                    $top_bottom = $matrix[$x][$y] + $match;
                                }else{
                                    $top_bottom = $matrix[$x][$y] + $mismatch;
                                }
                                $left_right = $matrix[$x][$y+1]+$gapp;
                                $up_down = $matrix[$x+1][$y]+$gapp;
                                if($top_bottom>=$left_right){
                                    if($top_bottom>=$up_down){
                                        $matrix[$x+1][$y+1] = $top_bottom;
                                    }else{
                                        $matrix[$x+1][$y+1] = $up_down;
                                    }
                                }else{
                                    if($left_right>=$up_down){
                                        $matrix[$x+1][$y+1] = $left_right;
                                    }else{
                                        $matrix[$x+1][$y+1] = $up_down;
                                    }
                                }
                            }
                        }

                        #calculate the number of total matches, and find the aligned chains
                        $y = $LOQ-1;
                        $x = $LOQ-1;
                        $matched = 0;
                        $totalAA = $LOQ;
                        do{
                            if(substr($Query,$y,1) eq substr($DBPiece,$x,1)){
                                $matched++;
                                $alignedX = substr($DBPiece,$x,1).$alignedX;
                                $alignedY = substr($Query,$y,1).$alignedY;
                                $y--;
                                $x--;
                            }elsif(($matrix[$x][$y]>=$matrix[$x+1][$y]) && ($matrix[$x][$y]>=$matrix[$x][$y+1])){
                                $alignedX = substr($DBPiece,$x,1).$alignedX;
                                $alignedY = substr($Query,$y,1).$alignedY;
                                $y--;
                                $x--;
                            }elsif($matrix[$x][$y+1]>=$matrix[$x+1][$y]){
                                $alignedX = substr($DBPiece,$x,1).$alignedX;
                                $alignedY = "-".$alignedY;
                                $totalAA++;
                                $x--;
                            }else{
                                $alignedX = "-".$alignedX;
                                $alignedY = substr($Query,$y,1).$alignedY;
                                $y--;
                            }
                        }while($x>=0 && $y>=0);
                        $score = $matched/$totalAA;
                        if($score>$threshold){
                            for($subcounter=0;$subcounter<scalar(@scores);$subcounter++){
                                if($score>=$scores[$subcounter]){
                                    splice(@scores,$subcounter,0,$score);
                                    splice(@resultsD,$subcounter,0,$alignedX);
                                    splice(@resultsQ,$subcounter,0,$alignedY);
                                    splice(@locations,$subcounter,0,$StartPoint." at chromosome ".$chrcounter);
                                    push(@{$locations{$chrcounter}},$StartPoint);
                                    last;
                                }
                            }
                        }

                        $alignedX = " ";
                        $alignedY = " ";
                    }
                }
            }
        }
    }

    #print the result in proper order
    print "\nQuery $QueryCounter: $#scores alignment results in total.\n";
    for($subcounter=0;$subcounter<$#scores;$subcounter++){
        printf "Score: %.4f%s\nPosition: %s\n",100*$scores[$subcounter],"%",$locations[$subcounter];        printf "%-10s%-s\n","Query:", $resultsQ[$subcounter];
        printf "%-10s%-s\n","DataBase:", $resultsD[$subcounter];
    }

    undef @scores;
    undef @resultsD;
    undef @resultsQ;
    undef @locations;
    undef %locations;
    undef @matrix;

    print "Finished the $QueryCounter query, at ".(time()-$startTime)."seconds\n";

}