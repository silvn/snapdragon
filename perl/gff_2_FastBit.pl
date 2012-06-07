#!/usr/bin/perl -w
use strict;
use FastBit::Loader;

my $outdir = shift @ARGV;
my $set    = shift @ARGV;

$set or die "usage: $0 <output directory> <name of data set> gff (on stdin)\n";

# GFF files are assumed to be grouped by chromosome
# It would help the start and end indexes if the GFF is sorted by start position.
# The outdir will have a subdirectory for each chromosome
# and each type of feature will be stored in it's own table (a subdirectory of $outdir)
# each -part.txt file will have a corresponding meta data for the chromosome
# metaTags = chr="$chromosome"
# columns will be defined as follows:
# source - key
# start - int
# end - int
# score - float
# strand - key
# phase - key
# the type of each attribute will be automatically detected - null values will
# be assigned as needed
#

my %col_type = (
    source  => 'CATEGORY',
    start   => 'UINT',
    end     => 'UINT',
    score   => 'FLOAT',
    strand  => 'BYTE',
    phase   => 'CATEGORY',
    ID      => 'TEXT',
    biotype => 'CATEGORY',
    Name    => 'TEXT',
    Parent  => 'TEXT'
);

my %data;

my $chr = '';
while (<>) {
    next if (/^#/);
    chomp;
    my (
        $seqname, $source, $type,  $start, $end,
        $score,   $strand, $phase, $attributes
    ) = split /\t/, $_;
    $start--;
    if ( $strand eq '+' ) {
        $strand = 1;
    }
    elsif ( $strand eq '-' ) {
        $strand = -1;
    }
    else {
        $strand = 0;
    }
    if ( $seqname ne $chr ) {

        if ($chr) {
            write_fastbit_data($chr);
            flush_data();
        }
        $chr = $seqname;
    }
    push @{ $data{$type}{source} }, $source;
    push @{ $data{$type}{start} },  $start;
    push @{ $data{$type}{end} },    $end;
    push @{ $data{$type}{score} }, ( $score eq '.' ) ? 0 : $score;
    push @{ $data{$type}{strand} }, $strand;
    push @{ $data{$type}{phase} },  $phase;

    for my $attribute ( split /\s*;\s*/, $attributes ) {
        my ( $k, $v ) = $attribute =~ m/^(\w+)=?\s?"?(.*?)"?$/;

#    next unless $k and $v; #### sometimes the attribute is present but no value is present
        push @{ $data{$type}{$k} }, $v;
    }
}
if ($chr) {
    write_fastbit_data($chr);
}
exit;

sub flush_data {
    for my $type ( keys %data ) {
        delete $data{$type};
    }
}

sub write_fastbit_data {
    for my $type ( keys %data ) {
        my @columns;
        my $nrows = @{$data{$type}{start}};
        system("mkdir -p $outdir/$type/$chr");
        my $part =
          FastBit::Loader->new( "$outdir/$type/$chr",
            { metaTags => "FBchr = $chr, FBfeature = $type, FBset = $set" } );
		my @sort_idx = sort {
			$data{$type}{start}[$a] <=> $data{$type}{start}[$b]
			or
			$data{$type}{end}[$a] <=> $data{$type}{end}[$b]
			} (0..$nrows-1);
        for my $col ( keys %{ $data{$type} } ) {
			if (@{$data{$type}{$col}} != $nrows) {
				print STDERR "WARNING: $type column $col doesn't contain $nrows rows\n";
				next;
			}
            $part->new_column( $col,
                exists $col_type{$col} ? $col_type{$col} : "TEXT",
				$col eq "start" ? {sorted => 'true'} : {} );
            $part->add_column_data( $col, $data{$type}{$col}, \@sort_idx );
        }
    }
}

sub sort_by {
	my @arrays = @_;
	my @idx = sort {$arrays[0][$a] <=> $arrays[0][$b]} (0..$arrays[0]); 
}
