#!/usr/bin/perl

use strict;
use warnings;
my ($peak_count,$peak_len,$trans_count,$out_file)=@ARGV;
our %PeakLen;
our %block;
our %PeakGen;
our %trans_cov;

open TRANS, "<$trans_count" or die $!;
my $aline=<TRANS>;
chomp $aline;
my @s_name;
my @arr=split/\t/, $aline;
foreach (@arr){
	if($_=~/.*hisat2_res\.(\S+)_Q20_sorted_nonrRNA.bam/){
		push(@s_name,$1);
	}
}
#print join(",",@s_name);
while (my $aline=<TRANS>) {
	chomp $aline;
	my @arr=split/\t/, $aline;
#	my @s_name;
#	if($.==1){
#		foreach (@arr){
#			if($_=~/.*hisat2_res\.(\S+)_Q20_sorted_nonrRNA.bam/){
#				push(@s_name,$1);
#			}
#		print join(",",@s_name);
#		}
#	}
#	else{
		foreach my $idx(2..$#arr){
			$trans_cov{$arr[0]}{$s_name[$idx-2]}=$arr[$idx];		
		}
#	}
}
close TRANS;

open LEN, "<$peak_len" or die $!;
while (my $bline=<LEN>) {
	chomp $bline;
	my @brr=split/\t/, $bline;
	$PeakLen{$brr[0]."\t".$brr[1]."\t".$brr[2]}=$brr[5];	
	$PeakGen{$brr[0]."\t".$brr[1]."\t".$brr[2]}=$brr[4];
}
close LEN;

open PEAK,"<$peak_count" or die $!;
open OUT, ">$out_file" or die $!;

print OUT "seqname\tstart\tend\tgene_id";

my $pline=<PEAK>;
chomp $pline;
my @prr=split/\t/, $pline;
my %s_idx;
foreach (3..(@prr-1)){
	$prr[$_]=~/\'(\S+)_Q20_sorted_nonrRNA.bam/;
	$s_idx{$_}=$1;
	print OUT "\t$1";
}
print OUT "\n";

while (my $pline=<PEAK>){
	chomp $pline;
	my @prr=split/\t/, $pline;
#	my %s_idx;
#	if($.==1){
#                foreach (3..(@prr-1)){
#                        $prr[$_]=~/\'(\S+)_Q20_sorted_nonrRNA.bam/;
#			$s_idx{$_}=$1;
#			print OUT "\t$1";
#                }
#		print OUT "\n";
#	}
#	else{
		my $peak_id=$prr[0]."\t".$prr[1]."\t".$prr[2];
		print OUT "$peak_id\t$PeakGen{$peak_id}";
		foreach (3..$#prr){
			my $sample_id=$s_idx{$_};
			if(defined $PeakLen{$peak_id} & defined $trans_cov{$PeakGen{$peak_id}}{$sample_id} & $PeakLen{$peak_id} > 0 & $trans_cov{$PeakGen{$peak_id}}{$sample_id} > 0){
#			if($PeakLen{$peak_id} > 0 & $trans_cov{$PeakGen{$peak_id}}{$sample_id} > 0){
			#	print "$PeakLen{$peak_id}\t$trans_cov{$PeakGen{$peak_id}}{$sample_id}\n";
				my $peak_score=$prr[$_]/($PeakLen{$peak_id} * $trans_cov{$PeakGen{$peak_id}}{$sample_id});
				print OUT "\t$peak_score";
			}
		}
		print OUT "\n";
#	}
}
