#!/usr/bin/env perl -w
use strict;

my @res;
for(my $i=0;$i<256;$i++) {
	my $b = dec2bin($i);
	my $rc = revcomp($b);
	print join(",",$i,$b,$rc,bin2dec($rc)),"\n";
	push @res, bin2dec($rc);
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
	$str =~ tr/01/10/;
	return reverse $str;
}
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/0{24}//;   # otherwise you'll get leading zeros
    return $str;
}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
