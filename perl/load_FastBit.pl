#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FastBit::Loader;
my ( $man, $help );
my (
    $part_dir, $key_column, $val_column, $vtype, $data_file,
    $format,   $delim,      $meta,       $buffer
);
GetOptions(
    'help|?'   => \$help,
    'man'      => \$man,
    'part=s'   => \$part_dir,
    'key=s'    => \$key_column,
    'val=s'    => \$val_column,
    'vtype=s'  => \$vtype,
    'data=s'   => \$data_file,
    'format=s' => \$format,
    'delim=s'  => \$delim,
    'meta=s'   => \$meta,
    'buffer=s' => \$buffer
);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

( $part_dir and $data_file or ( $format or ( $val_column and $vtype ) ) )
  or pod2usage(1);

$delim ||= "\t";

my $part =
  $meta
  ? FastBit::Loader->new( $part_dir, parse_meta($meta) )
  : FastBit::Loader->new($part_dir);

if ($format) {
    my $column_types = parse_format($format);
    $part->load_from_file( $data_file, $column_types, $delim, $buffer );
}
elsif ($val_column) {
    $part->has_column($val_column)
      and die "can't append data to existing column '$val_column'";
    $part->new_column( $val_column, $vtype );

    if ($key_column) {
        $part->has_column($key_column)
          or die "key column '$key_column' does not exist\n";
        $part->primary_key($key_column);
        $part->add_column_kv_data( $key_column, $val_column,
            read_kv_data($data_file) );
    }
    else {
		if($data_file) {
			$part->add_column_data( $val_column, read_data($data_file) );
		}
		else {
			$part->write_part_file();
		}
    }
}
exit;

sub read_data {
    my $infile = shift;
    my $fh     = can_opener($infile);
    my @data;
    while (<$fh>) {
        chomp;
        push @data, $_;
    }
    close $fh;
    return \@data;
}

sub read_kv_data {
    my $infile = shift;
    my $fh     = can_opener($infile);
    my %data;
    while (<$fh>) {
        chomp;
        my ( $k, $v ) = split /$delim/, $_;
        $k and $v or next;
        $data{$k} = $v;
    }
    close $fh;
    return \%data;
}

sub can_opener {
    my $filename = shift;
    my $cmd;
    if ( $filename =~ m/\.gz$/ ) {
        $cmd = "gzip -cd $filename";
    }
    elsif ( $filename =~ m/\.bz2$/ ) {
        $cmd = "bzip2 -cd $filename";
    }
    else {
        open( my $fh, "<", $filename ) or die "failed to open $filename : $!\n";
        return $fh;
    }
    open( my $fh, "-|", $cmd ) or die "failed to $cmd";
    return $fh;
}

sub parse_format {
    my $format = shift;
    my @cols = split /,/, $format;
    for ( my $i = 0 ; $i < @cols ; $i++ ) {
        my ( $name, $type, @extras ) = split /:/, $cols[$i];
        $cols[$i] = [ $name, $type ];
        if (@extras) {
            my %x = @extras;
            $cols[$i][2] = \%x;
        }
    }
    return \@cols;
}

sub parse_meta {
    my $meta = shift;
    my %hsh;
    for my $kv ( split /,/, $meta ) {
        my ( $k, $v ) = split /:/, $kv;
		$v =~ s/_comma_/, /g;
        $hsh{$k} = $v;
    }
    return \%hsh;
}

=head1 NAME

load_FastBit.pl - adds data to a FastBit partition

=head1 SYNOPSIS

	load_FastBit.pl -part <path to part dir> -val <column name> -vtype <data type> -data <data file>
	load_FastBit.pl -part <path to part dir> -val <column name> -vtype <data type> -data <data file> -key <key colum name>
	load_FastBit.pl -part <path to part dir> -data <data file> -format col1:FLOAT,col2:UINT:sorted:true,col3:CATEGORY

=head1 OPTIONS

		-help            brief help message
		-man             full documentation
		-part            FastBit partition
		-key             key column
		-val             value column
		-vtype           data type of value column
		-data            tab delimited data
		-format          comma separated list of name:type:extras for each column
		-meta            comma separated list of key:value pairs for the header

=head1 DESCRIPTION

This program adds new columns of data to a FastBit partition.  If it doesn't
exist, a new partition is created.  It can be used to load a tab (or other)
delimited file or just one column of data at a time.  If the partition has a
column that acts as a primary key, you can add a new column of data based on
key-value pairs.  For keys that are not already in the partition, new empty
values will be appended to every other column.

=cut
