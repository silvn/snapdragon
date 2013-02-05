#!/usr/bin/env perl -w
use strict;

use IO::Compress::Gzip qw(gzip $GzipError);
my $max_open_fh = 200;

my ($gi2tax_file,$nodes_file,$outdir) = @ARGV;
$outdir or die "usage: gzip -cd nt.gz | $0 gi_taxid_nucl.dmp(.gz) nodes.dmp <tax tree root>\n";
open(my $nodes_fh, "<", $nodes_file) or die "failed to open $nodes_file for reading: $!\n";
my %tax2parent;
while (<$nodes_fh>) {
	my @x = split /\t/, $_;
	$tax2parent{$x[0]} = $x[2];
}
close $nodes_fh;

print STDERR "finished filling %tax2parent from $nodes_file\n";

my $gi2tax_fh;
if ($gi2tax_file =~ m/\.gz$/) {
	open($gi2tax_fh,"gzip -cd $gi2tax_file |") or die "failed to 'gzip -cd $gi2tax_file |' : $!\n";
}
else {
	open($gi2tax_fh, "<", $gi2tax_file) or die "failed to open $gi2tax_file for reading: $!\n";
}
my %gi2tax;
my %tax2path;
my %active_fh;
my $tax;
while (<>) { # usage: gzip -cd nt.gz | split_nt.pl gi_taxid_nucl.dmp nodes.dmp outdir
	if (/^>(.+)/) {
		# parse the header line
		my @gis = ($_ =~ m/gi\|(\d+)\|/g);
		my $min_gi = $gis[0];
		for(my $i=1;$i<@gis;$i++) {
			$min_gi = $gis[$i] if ($gis[$i] < $min_gi);
		}
		# use $min_gi to lookup tax_id
		while (not exists $gi2tax{$min_gi}) {
			# empty the hash because all the keys are < $min_gi
			undef %gi2tax;
			# read 1000 lines
			for(my $i=0;$i<1000;$i++) {
				my $line = <$gi2tax_fh>;
				$line or last;
				chomp $line;
				my ($gi,$tax) = split /\t/, $line;
				$gi2tax{$gi} = $tax;
			}
		}
		$tax = $gi2tax{$min_gi};
		if (not exists $active_fh{$tax}) {
			if (not exists $tax2path{$tax}) {
				# construct the path
				my @path;
				while($tax != 1) {
					unshift @path, $tax;
					$tax = $tax2parent{$tax};
				}
				$tax2path{$tax} = join("/",$outdir,@path,"split.fa.gz");
			}
			# check if there are too many files open
			if (keys %active_fh == $max_open_fh) {
				# close something
				my ($t,$z) = each %active_fh;
				close $z;
				delete $active_fh{$t};
			}
			# open tax2path{$tax} in append mode and store in $active{$tax}
			$active_fh{$tax} = new IO::Compress::Gzip $tax2path{$tax}, Merge => 1 or die "failed to open $tax2path{$tax} in Merge mode: $!\n";
		}
	}
	# write the sequence to $active{$tax}
	$active_fh{$tax}->print($_);
}

close $gi2tax_fh;
for my $tax (keys %active_fh) {
	close $active_fh{$tax};
}


exit;

