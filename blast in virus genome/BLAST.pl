#!/usr/bin/perl
use strict;
use warnings;

#The defined variables
my $DBName = $ARGV[0];
my $QueryName = $ARGV[1];
my $WordLength = $ARGV[2];
my $threshold = $ARGV[3];

my $gapp = -2;
my $GOP = -1;
my $match = 2;
my $mismatch = -1;

#E-value is not required.
my $counter;
my $subcounter;
if(not defined $threshold){
    $threshold = 0.8;
    print "No threshold, default threshold = $threshold\n";
}
if(not defined $WordLength){
    $WordLength = 11;
    print "No word length, default length = $WordLength\n";
}

#data that related with the database
my %DB;
my $ALLSEQ = "";
my $ALLSEQLEN;
my $SpeciesName;

#data that related with the query sequence
#my $Query = "CATGAGGCTGCCCCGTATTCAGTGTCGCGGATTTGTACTGTCTG";
#my $Query = "CATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTG";

#build the library
open(DATABASE, "$DBName") or die("Cannot open the database!\n");
#open MYDB, ">>LIBRARY";


my @line = <DATABASE>;
my $first = shift @line;
if($first=~/>/){
    $SpeciesName = "";
}
$ALLSEQ = join("",@line);
$ALLSEQ =~ s/\s//g;
$ALLSEQLEN = length($ALLSEQ);

#read the query
open QUERYHANDLE, "$QueryName" or die("Cannot open the query!\n");
my $Query = <QUERYHANDLE>;
$Query =~ s/\s//g;
my $LOQ = length($Query);
my $LOD = $WordLength+2*int(($LOQ-$WordLength)/(2*$threshold)+2);
my @QueryWord;
my @QueryWordLocation;

#save the library to a hash
my $words;
for($counter=0;$counter<($ALLSEQLEN+1-$WordLength);$counter++){
    $words = substr($ALLSEQ,$counter,$WordLength);
    if(exists $DB{$words}){
        $DB{$words} = $DB{$words}."+$counter";
    }else{
        $DB{$words} = "$counter";
    }
}

=pod
#print out the library
my @WordList = sort keys %DB;
for my $word (@WordList){
    print MYDB $word.">".$DB{$word}."\n";
}
=cut

#process the query and do the allocation
my @temporary;
my $DBOS;
my @result;
my @resultsQ;
my @resultsD;
my @scores;
my @locations;
my $StartPoint;
$scores[0] = $threshold;
for($counter=0;$counter<($LOQ+1-$WordLength);$counter++){
    $words = substr($Query,$counter,$WordLength);
    push(@QueryWord,$words);
    if(exists $DB{$words}){
        $QueryWordLocation[$counter] = $DB{$words};
        @temporary = split(/\+/,$QueryWordLocation[$counter]);
        for my $location (@temporary){
            $StartPoint = $location-$counter;
            if(grep(/$StartPoint/,@locations)){
                next;
            }else{
                $DBOS = 2*int($counter/(2*$threshold)+1);
                @result = alignment($Query,substr($ALLSEQ,($location-$DBOS),$LOD));
                if($result[0]!=0){
                    for($subcounter=0;$subcounter<scalar(@scores);$subcounter++){
                        if($result[0]>=$scores[$subcounter]){
                            splice(@scores,$subcounter,0,$result[0]);
                            splice(@resultsD,$subcounter,0,$result[1]);
                            splice(@resultsQ,$subcounter,0,$result[2]);
                            splice(@locations,$subcounter,0,$StartPoint);
                            last;
                        }
                    }
                }
            }
        }         
    }else{
        $QueryWordLocation[$counter] = 0;
    }
}

#print the result in proper order
print "$#scores alignment results in total.\n";
for($subcounter=0;$subcounter<$#scores;$subcounter++){
    printf "Score: %.4f%s\nPosition: %s in the database.\n",100*$scores[$subcounter],"%",$locations[$subcounter];
    printf "%-16s%-s\n","Query:", $resultsQ[$subcounter];
    printf "%-16s%-s\n\n","DataBase:", $resultsD[$subcounter];
}

#0->up_down:from xi,y(i-1) to xi,yi
#1->left_right:from x(i-1),yi to xi,yi
#2->top_bottom:from x(i-1),y(i-1) to xi,yi
#sub max return the maximum of those three values and in which direction it comes from
#the diretion has its meaning if gap open penalty takes into count
sub max{
    my $max = pop(@_);
    my $direction = 0;
    for (my $counter = 0 ; $counter < scalar @_ ; $counter++) {
       if ($_[$counter]>$max){
            $max = $_[$counter];
            $direction = $counter+1;
       }
    }
    return $max,$direction;
}

#modified version
#align two sequences and return the score
sub alignment{
    my @chainY =  split(//,$_[0]);
    my @chainX =  split(//,$_[1]);

    my $counter;
    my @direction;
#initialization
    for ($counter = 0; $counter <= scalar(@chainX); $counter++){
        push(@direction,2);
    }

    my $alignedX;
    my $alignedY;
    my $score;

    my $matched = 0;
    my $totalAA = scalar @chainY;

    my $top_bottom = 0;
    my $up_down = 0;
    my $left_right =0;

    my $x;
    my $y;

    my @matrix;
#initiate score matrix
    $matrix[0][0] = 0;
    for($x=0;$x<scalar(@chainX);$x++){
        $matrix[$x+1][0] = ($x+1)*(-2);
    }
    for($y=0;$y<scalar(@chainY);$y++){
        $matrix[0][$y+1] = ($y+1)*(-2);
    }
    
    $x = 1;
    $y = 1;

#assign values to the matrix
    for my $Yaa (@chainY){
        for my $Xaa (@chainX){
            if ($Yaa eq $Xaa){
                $top_bottom = $matrix[$x-1][$y-1] + $match;
            }else{
                $top_bottom = $matrix[$x-1][$y-1] + $mismatch;
            }
            $left_right = $matrix[$x-1][$y]+$gapp+$GOP;
            $left_right = $left_right - $GOP if $direction[$x-1] == 1;
            $up_down = $matrix[$x][$y-1]+$gapp+$GOP;
            $up_down = $up_down - $GOP if $direction[$x] == 0;
            my @result = max($up_down,$left_right,$top_bottom);
            $matrix[$x][$y] = $result[0];
            $direction[$x] = $result[1];
            $x++;
        }
        $y++;
        $x=1;
    }

#calculate the number of total matches, and find the aligned chains
    $y = scalar(@chainY);
    $x = scalar(@chainX);
    do{
        if(($matrix[$x-1][$y-1]>=$matrix[$x][$y-1] && $matrix[$x-1][$y-1]>=$matrix[$x-1][$y]) || ($chainY[$y-1] eq $chainX[$x-1])){
            $matched++ if $chainY[$y-1] eq $chainX[$x-1];
            $alignedX = $alignedX.$chainX[$x-1];
            $alignedY = $alignedY.$chainY[$y-1];
            $y--;
            $x--;
        }elsif($matrix[$x-1][$y]>$matrix[$x-1][$y-1] && $matrix[$x-1][$y]>$matrix[$x][$y-1]){
            if($y==scalar(@chainY)){
                $x--;
            }else{
                $alignedX = $alignedX.$chainX[$x-1];
                $alignedY = $alignedY."-";
                $totalAA++;
                $x--;
            }
        }else{
            $alignedX = $alignedX."-";
            $alignedY = $alignedY.$chainY[$y-1];
            $y--;
        }
    }while($x>0 && $y>0);

    $score = $matched/$totalAA;

    if($score<$threshold){
        return  0;
    }else{
        my $temporary =$alignedY;
        my $diff;
        my $AAtemp;
        my $Matchedtemp;
        my $Scoretemp;
        my $switch = 1;
        while(0==0){
            $Matchedtemp =0;
            if($temporary =~ s/-+//){
                $AAtemp = length($temporary);
                $diff = length($alignedX)-$AAtemp;
                for($counter = $diff; $counter<length($alignedX); $counter++){
                    $Matchedtemp++ if (substr($alignedX,$counter,1) eq substr($temporary, $counter-$diff,1));
                }
                $Scoretemp = $Matchedtemp/$AAtemp;
                if($Scoretemp>=$score){
                    $alignedY = $temporary;
                    $alignedX = substr($alignedX,$diff);
                    $score = $Scoretemp;
                }else{
                    last;
                }
            }else{
                last;
            }
        }
        $alignedY = reverse $alignedY;
        $alignedX = reverse $alignedX;
        return  $score,$alignedX,$alignedY;
    }
    

#the offset may be useful
=pod
    my $temporary = $alignedX;
    $temporary =~ s/-//g;
    $offset = length($_[1])-$offset-length($temporary);
=cut

=pod
#This is another way to calculate identity score
    for(my $counter = 0; $counter <= $#alignedX; $counter++){
        $matched++ if $alignedX[$counter] eq $alignedY[$counter];
        $totalAA = scalar(@alignedX);
    }
=cut
}