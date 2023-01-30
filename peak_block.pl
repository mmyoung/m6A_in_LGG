#!/usr/bin/perl

use strict;
use warnings;
my ($in_file,$out_file)=@ARGV;
our %peak;
our %block;

open ME,"<$in_file" or die $!; 
while (my $aline=<ME>) {
        chomp($aline);
        my @tmp=split/\t/,$aline;
	my $peak_id=$tmp[0]."\t".$tmp[1]."\t".$tmp[2];
        push @{$block{$peak_id}},$tmp[13]; 
		if(! defined $peak{$peak_id}){
                $peak{$peak_id}=$tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[4];
        }
}
close ME;
open OUT,">$out_file" or die $!;
foreach my $id (keys %peak){
	my $len=0;
	foreach my $block_size (@{$block{$id}}){
		$len+=$block_size;
	}
        #my $blockSize = join(',',@{$block{$id}});    
        print OUT "$peak{$id}\t$len\n";
}

