#!/usr/bin/env perl -w
use strict;

my @dec2nuc;
my %nuc2dec;
my @nuc = qw(A C G T);
for my $n1 (@nuc) {
	for my $n2 (@nuc) {
		for my $n3 (@nuc) {
			for my $n4 (@nuc) {
				my $mer = $n1.$n2.$n3.$n4;
				push @dec2nuc, $mer;
				$nuc2dec{$mer} = @dec2nuc - 1;
			}
		}
	}
}

my @res;
for(my $i=0;$i<256;$i++) {
	my $n = $dec2nuc[$i];
	my $rc = revcomp($n);
	print join(",",$i,$n,$rc,$nuc2dec{$rc}),"\n";
	push @res, $nuc2dec{$rc};
}
for(my $i=0;$i<16;$i++) {
	for(my $j=0;$j<16;$j++) {
		my $k = $i*16 + $j;
		print "$res[$k],";
	}
	print "\n";
}

sub revcomp {
	my $str = shift;
	print "str  $str\n";
	
	$str =~ tr/ACGT/TGCA/;
	print "flip $str\n";
	return reverse $str;
}
