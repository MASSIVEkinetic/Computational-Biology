#!/usr/bin/perl
use strict;
use warnings;
use Storable;

my %DB = %{retrieve('LIBRARY')} or die("Cannot open the library!\n");
print "\nFinish loading, keys: ".(scalar keys %DB);

my @entries = sort keys %DB;
my @temporary;
my %group;
my %GroupData;
my $entry;
for $entry (@entries){
    if($group{substr($entry,0,6)}){
        $group{substr($entry,0,6)} = $group{substr($entry,0,6)}."+".$entry;
    }else{
        $group{substr($entry,0,6)} = $entry;
    }
}

for my $name (keys %group){
    @temporary = split(/\+/,$group{$name});
    for $entry (@temporary){
        $GroupData{$entry} = $DB{$entry};
    }
    store \%GroupData,"lib/$name";
    undef %GroupData;
}


