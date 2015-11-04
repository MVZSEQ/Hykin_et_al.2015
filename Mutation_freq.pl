#!/usr/bin/perl

#read_mutfreq.pl

use warnings;
use strict;
use Getopt::Std;

my %opts = (n=>100);
getopts('n:', \%opts);
die (qq/
Usage: plot mutation frequencies along reads

read_mutfreq.pl [options] <in.sam>

Options: -n INT   number of positions along read to examine

Note: This script only uses reads with bitwise flag 99 and 147

\n/) unless @ARGV;

#set up hashes to count mutations, keys are reference_alternate
my %forward = (refmatch=>0, AC=>0, AG=>0, AT=>0, GA=>0, GC=>0, 
GT=>0, CA=>0, CG=>0, CT=>0, TA=>0, TG=>0, TC=>0);
foreach (keys %forward) {
	my $setup;
	for (my $i = 0; $i <= $opts{n}-1; $i++) {
		$setup->[$i] = 0;
	}
	$forward{$_} = $setup;
}
my %reverse = (refmatch=>0, AC=>0, AG=>0, AT=>0, GA=>0, GC=>0, 
GT=>0, CA=>0, CG=>0, CT=>0, TA=>0, TG=>0, TC=>0);
foreach (keys %reverse) {
	my $setup;
	for (my $i = 0; $i <= $opts{n}-1; $i++) {
		$setup->[$i] = 0;
	}
	$reverse{$_} = $setup;
}
# starting processing SAM file
open(SAM,'<',$ARGV[0]);
my ($read, $ref, $type, $length, $elem);
my $process = 0; # for tracking progress
while (<SAM>) {
	next if $_ =~ /^\@HD\s+|^\@PG\s+|^\@SQ\s+/i;
	my @s = split /\t/, $_;
	next unless ($s[1] == 99 || $s[1] == 147); # try adding 83 and 163 bitwise flags later
	my @seq = split "", $s[9];
	my @md = $1 =~ /(\^[A-Za-z]+|0[A-Za-z]|[A-Za-z]|\d+)/g if $_ =~ /MD:Z:([^\s+]+)\s+/i;
	if ($s[1] == 147) { # reverse CIGAR, SEQ, MD:Z: 
		my @revcig = $s[5] =~ /(\d+[A-Za-z])/ig;
		@revcig = @revcig[reverse(0 .. $#revcig)];
		$s[5] = join "", @revcig;
		@seq = @seq[reverse(0 .. $#seq)];
		@md = @md[reverse(0 .. $#md)]; 
	}
	my $pos = 0;
	my %insert;
# adjust for insertions
	if ($s[5] =~ /I/) {
		my @mis;
		if ($s[5] =~ /(^\d+H\d+S|^\d+H|^\d+S)/) {
			@mis = $' =~ /(\d+[A-Za-z])/g;
		} else {
			@mis = $s[5] =~ /(\d+[A-Za-z])/g;
		}
		my $adj = 0;
		my @adj;
		push @adj, $adj;
		my $ipos;
		while (my $k = shift @mis) {
			if ($k =~ /(\d+)[^D]/) {
				push @adj, $adj;
				$adj = $adj + $1; 
			if ($k =~ /(\d+)I/) {
				$ipos = (pop @adj);
				$insert{$ipos} = $1;
				}	
			}
		}
	}	
#subtract soft-clipped positions from SEQ
	if ($s[5] =~ /(\d+)S/) {
		for (my $j = 0; $j < $1; $j++) {
			shift @seq;
		}
	}
# read1 forward/read1 reverse (flipped) for 5-3' along frag
	if ($s[1] == 99) {
		foreach (@md) {		
			if ($_ =~ /(^0[A-Za-z]|^[A-Za-z])/) {
				$pos += $insert{$pos} if $insert{$pos}; #adjust for insert
				if ($seq[$pos] =~ /N/i) {
					$pos++; last if $pos > $opts{n}-1;
					next;
				}
				my $base = $1; $base =~ s/0//;
				$type = $base . $seq[$pos];
				${$forward{$type}}[$pos]++;
				$pos++; last if $pos > $opts{n}-1
			} elsif ($_ =~ /^(\d+)$/) {
				next if $1 == 0;
				my $start = $pos;
				for (my $i = $start; $i < $start + $1; $i++) {
					$pos += $insert{$pos} if $insert{$pos}; #adjust for insert
					${$forward{refmatch}}[$pos]++;
					$pos++; last if $pos > $opts{n}-1;
				}
			} elsif ($_ =~ /\^[A-Za-z]+/) {
				next;
			}
		}
	$process++;
	print STDERR "$process reads analyzed...\n" if ($process % 5000 == 0);
	} elsif ($s[1] == 147) { #read2 reverse/read2 forwards (flipped) for 3-5' along frag
		foreach (@md) {		
			if ($_ =~ /(^0[A-Za-z]|^[A-Za-z])/) {
				$pos += $insert{$pos} if $insert{$pos}; #adjust for insert
				if ($seq[$pos] =~ /N/i) {
					$pos++; last if $pos > $opts{n}-1;
					next;
				}
				my $base = $1; $base =~ s/0//;
				$type = $base . $seq[$pos];
				${$reverse{$type}}[$pos]++;
				$pos++; last if $pos > $opts{n}-1;
			} elsif ($_ =~ /^(\d+)$/) {
				next if $1 == 0;
				my $start = $pos;
				for (my $i = $start; $i < $start + $1; $i++) {
					$pos += $insert{$pos} if $insert{$pos}; #adjust for insert
					${$reverse{refmatch}}[$pos]++;
					$pos++; last if $pos > $opts{n}-1;
				}
			} elsif ($_ =~ /\^[A-Za-z]+/) {
				next;
			}
		}
	$process++;
	print STDERR "$process reads analyzed...\n" if ($process % 5000 == 0);
	}
}
print STDERR "$process total reads analyzed...finished!\n";
close SAM;
# mutation type (and ref match) frequency at each site
for(my $i=0; $i <= $opts{n} - 1; $i++) {
	my $total = 0;
	foreach (keys %forward) {
		$total += ${$forward{$_}}[$i] if (${$forward{$_}}[$i]);
	}
	foreach (keys %forward) {
		${$forward{$_}}[$i] /= $total if (${$forward{$_}}[$i]);
	}
}
for(my $i=0; $i <= $opts{n} - 1; $i++) {
	my $total = 0;
	foreach (keys %reverse) {
		$total += ${$reverse{$_}}[$i] if (${$reverse{$_}}[$i]);
	}
	foreach (keys %reverse) {
		${$reverse{$_}}[$i] /= $total if (${$reverse{$_}}[$i]);
	}
}
#print results: 1 for read1-forward (5'-3') & 2 for read2-reverse (3'-5')
print "ref1\tAC1\tAG1\tAT1\tGA1\tGC1\tGT1\tCA1\tCG1\tCT1\tTA1\tTC1\tTG1\n";
for (my $i=0; $i <= $opts{n} - 1; $i++) {
	print "${$forward{refmatch}}[$i]\t";
	print "${$forward{AC}}[$i]\t";
	print "${$forward{AG}}[$i]\t";
	print "${$forward{AT}}[$i]\t";
	print "${$forward{GA}}[$i]\t";
	print "${$forward{GC}}[$i]\t";
	print "${$forward{GT}}[$i]\t";
	print "${$forward{CA}}[$i]\t";
	print "${$forward{CG}}[$i]\t";
	print "${$forward{CT}}[$i]\t";
	print "${$forward{TA}}[$i]\t";
	print "${$forward{TC}}[$i]\t";
	print "${$forward{TG}}[$i]\n";
}
print "\n\n"; #space between forward and reverse
print "ref2\tAC2\tAG2\tAT2\tGA2\tGC2\tGT2\tCA2\tCG2\tCT2\tTA2\tTC2\tTG2\n";
for (my $i=0; $i <= $opts{n} - 1; $i++) {
	print "${$reverse{refmatch}}[$i]\t";
	print "${$reverse{AC}}[$i]\t";
	print "${$reverse{AG}}[$i]\t";
	print "${$reverse{AT}}[$i]\t";
	print "${$reverse{GA}}[$i]\t";
	print "${$reverse{GC}}[$i]\t";
	print "${$reverse{GT}}[$i]\t";
	print "${$reverse{CA}}[$i]\t";
	print "${$reverse{CG}}[$i]\t";
	print "${$reverse{CT}}[$i]\t";
	print "${$reverse{TA}}[$i]\t";
	print "${$reverse{TC}}[$i]\t";
	print "${$reverse{TG}}[$i]\n";
}
