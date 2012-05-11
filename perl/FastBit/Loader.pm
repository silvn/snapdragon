
=head1 NAME

FastBit::Loader - A module for loading data into a FastBit partition

=head1 SYNOPSIS

	use FastBit::Loader;
	my $partition = FastBit::Loader->new( "part_dir" );
	$partition->new_column( "col1", "FLOAT" );
	$partition->add_column_data( "col1", \@col1_data );
	$partition->load_from_file( $filename,
		[[ 'col1', 'FLOAT' ],
		 [ 'col2', 'UINT', {sorted=>'true'} ],
		 [ 'col3', 'CATEGORY' ]
		], ",");


The valid FastBit data types are:

UBYTE, BYTE, USHORT, SHORT, UINT, INT, FLOAT, DOUBLE, CATEGORY, TEXT

If the column type doesn't match one of these words, the guess_type()
function will examine the data and try to choose the most compact
data type required to hold the data.

=head1 DESCRIPTION

This module lets you create an empty FastBit partition and add data
to it, or add data to an existing partition.  Since FastBit is a column
based storage engine, it is trivial to add additional columns to
a partition.  If a partition has a primary key, it can be used to add
a column of data provided as key value pairs.  Any new keys are added
at the end of the primary key column and other data columns are padded
with nulls.

=head2 METHODS

=over 12

=item C<new($path,\%meta)>

Given a path to a directory, it will load an existing partition into memory
or create an empty partition at the specified location.  An optional hash
reference can be used to set parameters in the header section of the --part.txt
file.  For example, {metaTags => "FBchr = 1, FBfeature = exon"}

=item C<new_column($name,$type)>

Creates a new column in the partition called $name of data type $type.
Valid types include: UBYTE,BYTE,USHORT,SHORT,UINT,INT,FLOAT,DOUBLE,TEXT,CATEGORY

=item C<load_from_file($filename,\@columns,$delim)>

Populates the partition with $delim delimited data in $filename.  Column
[names, types, and extras] are passed as a list of array refs.  The default
delimiter is a tab.  Automatically sniffs for .bz2 or .gz suffixes.


=item C<add_column_data($name,\@data)>

Packs the array of data into the appropriate format for the named column
and writes the packed representation to disk.

=item C<add_column_kv_data($kcol,$vcol,\%data)>

The column $kcol should be a primary key and $vcol is a new column that
will hold the associated data from the hash.  The order of the data in $kcol
has to be preserved, so new keys are appended at the end, and other columns
get padded with nulls (of the appropriate type).

=item C<primary_key($name)>

This is the getter/setter function for the primary key of a partition.

=back

=item C<has_column($name)>

Checks if the named column exists in the partition

=item C<columns()>

Returns the list of column names

=back

=head1 LICENSE

This is release under whatever license the Ware lab uses for its code.

=head1 AUTHOR

Andrew Olson

=cut

package FastBit::Loader;
use strict;
use feature "switch";

use constant MAXUBYTE    => 255;
use constant MAXBYTE     => 127;
use constant MINBYTE     => -128;
use constant MAXUSHORT   => 65535;
use constant MAXSHORT    => 32767;
use constant MINSHORT    => -32768;
use constant MAXCATEGORY => 256;

my %pack_type = (
    UBYTE    => 'C',
    BYTE     => 'c',
    USHORT   => 'S',
    SHORT    => 's',
    UINT     => 'I',
    INT      => 'i',
    FLOAT    => 'f',
    DOUBLE   => 'd',
    CATEGORY => '(Z*)',
    TEXT     => '(Z*)'
);

my %doublecheck = ( UBYTE => 1, BYTE => 1, USHORT => 1, SHORT => 1, UINT => 1 );

sub new {
    my $name     = shift;
    my $part_dir = shift;
    my $extras   = shift;
    $part_dir or die "FastBit::Loader:: missing partition directory\n";
    my $class = ref($name) || $name;
    my $this = {};
    bless $this, $class;
    $this->init( $part_dir, $extras );
    return $this;
}

sub init {
    my ( $this, $part_dir, $extras ) = @_;

    $this->{_part_dir} = $part_dir;

    if ( -f "$part_dir/-part.txt" ) {
        $this->load_partition($part_dir);
    }
    else {
        ( $this->{Name} ) = $part_dir =~ m/(\w+)\/?$/;
        $this->{_columns}          = {};
        $this->{number_of_rows}    = 0;
        $this->{number_of_columns} = 0;
        system("mkdir -p $part_dir")
          and die "failed to mkdir -p $part_dir: $!\n";
    }
    if ($extras) {
        while ( my ( $k, $v ) = each %$extras ) {
            $this->{$k} = $v;
        }
    }
}

sub load_partition {
    my ( $this, $part_dir ) = @_;
    open( my $fh, "<", "$part_dir/-part.txt" )
      or die "failed to open $part_dir/-part.txt for reading : $!\n";
    my $header = parse_section( $fh, 'header' );
    while ( my ( $k, $v ) = each %$header ) {
        $this->{$k} = $v;
    }
    while ( my $column = parse_section( $fh, 'column' ) ) {
        $this->{_columns}{ $column->{name} } = $column;
        if ( $column->{_primary_key} ) {
            if ( $this->{_primary_key} ) {
                die "partition $part_dir should only have one primary key\n";
            }
            $this->{_primary_key} = $column->{name};
        }
    }
    close $fh;
}

sub parse_section {
    my $fh      = shift;
    my $section = shift;
    my %meta;
    while (<$fh>) { last if /^BEGIN\s+$section/i; }
    while (<$fh>) {
        last if /^END\s$section/i;
        my ( $k, $v ) = /(\w+)\s*=\s*(.+)/;
        $v =~ s/\s+$//;
        if ( substr( $v, 0, 1 ) eq '"' ) {
            $v =~ s/"//g;
        }
        $meta{ lc $k } = $v;
    }
    return keys %meta ? \%meta : undef;
}

sub write_part_file {
    my $this = shift;
    open( my $fh, ">", $this->{_part_dir} . "/-part.txt" )
      or die "failed to open -part.txt for writing : $!\n";
    write_section( $fh, "HEADER", $this );
    while ( my ( $col_name, $col_meta ) = each %{ $this->{_columns} } ) {
        write_section( $fh, "Column", $col_meta );
    }
    close $fh;
}

sub write_section {
    my ( $fh, $section, $hsh ) = @_;
    print $fh "BEGIN $section\n";
    while ( my ( $k, $v ) = each %$hsh ) {
        next if ref($v);
        if ( $k ne "metaTags" and $v =~ m/\s/ ) {
            $v = '"' . $v . '"';
        }
        print $fh "$k = $v\n";
    }
    print $fh "END $section\n\n";
}

sub has_column {
    my ( $this, $col_name ) = @_;
    return exists $this->{_columns}{$col_name};
}

sub new_column {
    my ( $this, $col_name, $col_type, $extras ) = @_;
    my %column = (
        name      => $col_name,
        data_type => $col_type
    );
    $this->{_columns}{$col_name} = \%column;
    $this->{number_of_columns}++;
    if ($extras) {
        while ( my ( $k, $v ) = each %$extras ) {
            $this->{_columns}{$col_name}{$k} = $v;
        }
    }
}

sub columns {
    my $this = shift;
    return keys %{ $this->{_columns} };
}

# getter/setter of primary key
sub primary_key {
    my $this     = shift;
    my $col_name = shift;
    if ($col_name) {
        if ( $this->{_primary_key} and $this->{_primary_key} ne $col_name ) {
            die
"partition alreay has a primary key : $this->{_primary_key} ne $col_name\n";
        }
        $this->{_primary_key} = $col_name;
    }
    else {
        return $this->{_primary_key} ? $this->{_primary_key} : undef;
    }
}

sub add_column_data {
    my ( $this, $col_name, $data, $idx ) = @_;
    my $nrows = @$data;
    $this->{number_of_rows} == 0
      or $this->{number_of_rows} == $nrows
      or die
"add_column_data($col_name) failed: number of rows doesn't match partition\n";

    $this->{number_of_rows} = $nrows;

    my $col = $this->{_columns}{$col_name};
    $col
      or die
      "add_column_data($col_name) failed: column $col_name doesn't exist\n";

    if ( not exists $pack_type{ $col->{data_type} } ) {
        my $data_type = guess_type($data);
        print STDERR
          "updated data_type from '$col->{data_type}' to '$data_type'\n";
        $col->{data_type} = $data_type;
    }
    my $format = $pack_type{ $col->{data_type} };

	if ($idx) {
		$data = reorder($data,$idx);
	}
    $this->write_column( $col_name, $format, $data );
    $this->write_part_file();
}

sub reorder {
	my ($in,$idx) = @_;
	my @out;
	for my $i (@$idx) {
		push @out, $in->[$i];
	}
	return \@out;
}

sub add_column_kv_data {
    my ( $this, $kcol, $vcol, $data ) = @_;

    # check if $kcol is primary key
    # read $kcol values into an array
    # iterate over the array and mark old data
    # append new keys to key column
    # recheck data type of key column
    # write column data for vcol
    # rewrite column data for kcol
    # append null data for all other columns

    if ( $kcol ne $this->primary_key() ) {
        die "$kcol is not the primary key\n";
    }
    my $kdata = $this->read_column($kcol);
    my ( @vdata, %seen );
    my $v_column = $this->{_columns}{$vcol};
    my $v_format = $pack_type{ $v_column->{data_type} };
    my $vnull    = $v_format eq '(Z*)' ? '' : 0;
    for my $k (@$kdata) {
        if ( exists $data->{$k} ) {
            push @vdata, $data->{$k};
            $seen{$k} = 1;
        }
        else {
            push @vdata, $vnull;
        }
    }
    my $new = 0;
    while ( my ( $k, $v ) = each %$data ) {
        if ( not $seen{$k} ) {
            push @$kdata, $k;
            push @vdata,  $v;
            $new++;
        }
    }
    if ($new) {
        $this->{number_of_rows} += $new;
        my $k_column = $this->{_columns}{$kcol};
        if ( exists $doublecheck{ $k_column->{data_type} } ) {
            $k_column->{data_type} = guess_type($kdata);
        }
        $this->write_column( $kcol, $pack_type{ $k_column->{data_type} },
            $kdata );

        while ( my ( $col_name, $col ) = each %{ $this->{_columns} } ) {
            next if ( $col_name eq $kcol );
            next if ( $col_name eq $vcol );

            # append $new nulls to $col_name
            $this->append_nulls( $col_name, $new );
        }

    }
    if ( exists $doublecheck{ $v_column->{data_type} } ) {
        $v_column->{data_type} = guess_type( \@vdata );
    }
    $this->write_column( $vcol, $pack_type{ $v_column->{data_type} }, \@vdata );
    $this->write_part_file();
}

sub load_from_file {
    my ( $this, $filename, $column_types, $delim, $buffer ) = @_;
    $delim ||= "\t";
    my $fh;
    if ( $filename =~ m/\.gz$/ ) {
        open( $fh, "-|", "gzip -cd $filename" )
          or die "failed to gzip -cd $filename : $!\n";
    }
    elsif ( $filename =~ m/\.bz2$/ ) {
        open( $fh, "-|", "bzip2 -cd $filename" )
          or die "failed to bzip2 -cd $filename : $!\n";
    }
    else {
        open( $fh, "<", $filename ) or die "failed to open $filename : $!\n";
    }
    my $n_cols = @$column_types;
    my @data;
    my $nrows = 0;
    while (<$fh>) {
        next if /^#/;
        chomp;
        my @cols = split /$delim/, $_;
        @cols == $n_cols
          or die "number of columns in line ("
          . scalar @cols
          . ") != $n_cols\n";
        for ( my $i = 0 ; $i < $n_cols ; $i++ ) {
            push @{ $data[$i] }, $cols[$i];
        }
        $nrows++;
        if ( $buffer and @{$data[0]} == $buffer ) {
            for ( my $i = 0 ; $i < $n_cols ; $i++ ) {
                if ( $nrows == @{$data[$i]} ) {
                    $this->new_column( @{ $column_types->[$i] } );
                    $this->add_column_data( $column_types->[$i][0], $data[$i] );
                }
                else {
                    $this->write_column(
                        $column_types->[$i][0],
                        $pack_type{
                            $this->{_columns}{ $column_types->[$i][0] }
                              ->{data_type}
                          },
                        $data[$i],">>"
                    );
                }
				$data[$i] = [];
            }
        }
    }
    close $fh;
    for ( my $i = 0 ; $i < $n_cols ; $i++ ) {
        if ( $nrows == @{$data[$i]} ) {
            $this->new_column( @{ $column_types->[$i] } );
            $this->add_column_data( $column_types->[$i][0], $data[$i] );
        }
        else {
            $this->write_column(
                $column_types->[$i][0],
                $pack_type{ $this->{_columns}{ $column_types->[$i][0] }
                      ->{data_type} },
                $data[$i],">>"
            );
        }
    }
	$this->{number_of_rows} = $nrows;
	$this->write_part_file()
}

sub append_nulls {
    my ( $this, $col_name, $n ) = @_;
    open( my $colfh, ">>", $this->{_part_dir} . "/$col_name" )
      or die
      "append_nulls($col_name, $n) failed to open data file for writing: $!\n";

    my $format = $pack_type{ $this->{_columns}{$col_name}{data_type} };
    my $null   = $format eq '(Z*)' ? '' : 0;
    my @nulls  = ($null) x $n;
    print $colfh pack( "$format$n", @nulls );
    close $colfh;
}

sub read_column {
    my ( $this, $col_name ) = @_;
    open( my $colfh, "<", $this->{_part_dir} . "/$col_name" )
      or die "failed to open $col_name: $!\n";
    my @stat = stat $colfh;
    read $colfh, my $buffer, $stat[7];
    close $colfh;
    my $format = $pack_type{ $this->{_columns}{$col_name}{data_type} };
    my @values = unpack( $format . '*', $buffer );
    return \@values;
}

sub write_column {
    my ( $this, $col_name, $format, $data, $mode ) = @_;
	$mode ||= ">";
    open( my $colfh, $mode, $this->{_part_dir} . "/$col_name" )
      or die
      "write_column($col_name) failed to open data file for writing: $!\n";

    print $colfh pack( $format . '*', @$data );
    close $colfh;
}

sub guess_type {
    my $data = shift;
    my $type = "UINT";
    my ( $min, $max ) = ( 'NaN', 'NaN' );
    for my $v (@$data) {
        given ($type) {
            when ("UINT") {
                given ($v) {
                    when (/^\d+$/) {
                        $max = $v if ( $max eq 'NaN' or $v > $max );
                    }
                    when (/^-\d+$/) {
                        $type = "INT";
                        $min = $v if ( $min eq 'NaN' or $v < $min );
                    }
                    when (/^-?\d?\.\d+$/) { $type = "FLOAT"; }
                    default { $type = "str"; last; }
                }
            }
            when ("INT") {
                given ($v) {
                    when (/^-?\d+$/) {
                        $min = $v if ( $min eq 'NaN' or $v < $min );
                        $max = $v if ( $max eq 'NaN' or $v > $max );
                    }
                    when (/^-?\d?\.\d+$/) { $type = "FLOAT"; }
                    default { $type = "str"; last; }
                }
            }
            when ("FLOAT") {
                if ( $v !~ m/^-?\d?\.\d+$/ ) {
                    $type = "str";
                    last;
                }
            }
        }
    }
    return refine_type( $type, $max, $min, scalar @$data );
}

sub refine_type {
    my ( $type, $max, $min, $nvalues ) = @_;
    if ( $type eq 'UINT' ) {
        if ( $max <= MAXUBYTE ) {
            $type = 'UBYTE';
        }
        elsif ( $max <= MAXUSHORT ) {
            $type = 'USHORT';
        }
    }
    elsif ( $type eq 'INT' ) {
        if ( $max <= MAXBYTE and $min >= MINBYTE ) {
            $type = 'BYTE';
        }
        elsif ( $max <= MAXSHORT and $min >= MINSHORT ) {
            $type = 'SHORT';
        }
    }
    if ( $type eq "str" ) {
        $type = $nvalues > MAXCATEGORY ? "TEXT" : "CATEGORY";
    }
    return $type;
}

1;
