#!/usr/bin/perl
use strict;
use warnings;

my $gapp = -2;
my $GOP = -1;
my $match = 2;
my $mismatch = -1;

print "Enter file name(.fasta or other file that contains amino acid sequences):\n";
my $name = <STDIN>;
chomp $name;

open FILE,"$name";
my @fasta = <FILE>;

my $ChainCounter = -1;
my @Chains;
my @SpeciesName;

for my $line (@fasta){
    if($line=~/>/){
        $ChainCounter++;
        $Chains[$ChainCounter] = "";
        $SpeciesName[$ChainCounter] = substr($1,0,1).".".$2 if $line =~ /\[(\w*)\s(\w*)\]/g;
    }else{
        $line =~ s/\s//g;
        $Chains[$ChainCounter] = $Chains[$ChainCounter].$line;
    }
}

arrangement(@Chains);

#up_down:from xi,y(i-1) to xi,yi
#left_right:from x(i-1),yi to xi,yi
#top_bottom:from x(i-1),y(i-1) to xi,yi
#sub max return the maximum of those three values and in which direction it comes from
#the diretion has its meaning if gap open penalty takes into count
sub max{
    my $max = pop(@_);
    my $temporary = 0;
    my $direction;
    for (my $counter = 0 ; $counter < scalar @_ ; $counter++) {
       if ($_[$counter]>$max){
            $max = $_[$counter];
            $temporary = $counter+1;
       }
    }
    if($temporary == 2){
        $direction = "top_bottom";
    }elsif($temporary == 1){
        $direction = "left_right";
    }else{
        $direction = "up_down";
    }
    return $max,$direction;
}

sub arrangement{
    my @chains = @_;
    
    my @ScoreMatrix;
    my $x;
    my $y;

    my $total = scalar(@chains);

    #set up the final matrix which shows diferent scores
    $ScoreMatrix[0][0] = "species";
    for($x=0;$x<$total;$x++){
        if(@SpeciesName){
            $ScoreMatrix[$x+1][0] = $SpeciesName[$x];
        }else{
            $ScoreMatrix[$x+1][0] = ($x+1);
        }
        $ScoreMatrix[0][$x+1] = $ScoreMatrix[$x+1][0];
    }

=pod
these lines are useless, they are used to fill the matrix with zeros;     
    for($y=0;$y<$total;$y++){
        for($x=0;$x<$total;$x++){
            $ScoreMatrix[$x+1][$y+1] = 0;
        }
    }
=cut

    for($y=1;$y<$total;$y++){
        my $DefaultChain = shift @chains;
        $x = $y+1;
        for my $chain (@chains){
            my @result = alignment($DefaultChain,$chain);
            $ScoreMatrix[$x][$y] = (sprintf "%.4f",100*$result[0])."%";
            print "\nThe alignment of polypeptide chain from $ScoreMatrix[0][$y] and $ScoreMatrix[0][$x]:\n";
            printf "%-20s%-s\n","$ScoreMatrix[0][$y]: ",join("",@{$result[2]});
            printf "%-20s%-s\n\n","$ScoreMatrix[0][$x]: ",join("",@{$result[1]});
            $x++;
        }
    }

    for($y=0;$y<=$total;$y++){
        for($x=0;$x<=$total;$x++){
            if($x==0 && $y==0){
                printf '%-20s', $ScoreMatrix[0][0];
            }elsif($x<$y){
                printf '%-20s', $ScoreMatrix[$y][$x];
            }elsif($x==$y){
                printf '%-20s', "100.0000%";
            }else{
                printf '%-20s', $ScoreMatrix[$x][$y];
            }
        }
        print "\n";
    }
}

sub alignment{
    my @chainY =  split(//,$_[0]);
    my @chainX =  split(//,$_[1]);

    my @direction;
#initialization
    for (my $counter = 0; $counter <= scalar(@chainX); $counter++){
        push(@direction,"top_bottom");
    }

    my @alignedX;
    my @alignedY;

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
            $left_right = $left_right - $GOP if $direction[$x-1] eq "left_right";
            $up_down = $matrix[$x][$y-1]+$gapp+$GOP;
            $up_down = $up_down - $GOP if $direction[$x] eq "up_down";
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
            unshift(@alignedX, $chainX[$x-1]);
            unshift(@alignedY, $chainY[$y-1]);
            $y--;
            $x--;
        }elsif($matrix[$x-1][$y]>$matrix[$x-1][$y-1] && $matrix[$x-1][$y]>$matrix[$x][$y-1]){
            unshift(@alignedX, $chainX[$x-1]);
            unshift(@alignedY, "-");
            $totalAA++;
            $x--;
        }else{
            unshift(@alignedX, "-");
            unshift(@alignedY, $chainY[$y-1]);
            $y--;
        }
    }while($x>0 && $y>0);

=pod
#This is another way to calculate identity score
    for(my $counter = 0; $counter <= $#alignedX; $counter++){
        $matched++ if $alignedX[$counter] eq $alignedY[$counter];
        $totalAA = scalar(@alignedX);
    }
=cut

    return  $matched/$totalAA,\@alignedX,\@alignedY;
}