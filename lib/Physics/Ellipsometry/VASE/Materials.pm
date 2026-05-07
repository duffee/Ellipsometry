package Physics::Ellipsometry::VASE::Materials;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Exporter 'import';

our @EXPORT_OK = qw(load_material interpolate_material);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::Materials - Optical constants file loader

=head1 DESCRIPTION

Loads point-by-point (PBP) optical constant files (.mat format) used by
WVASE and other ellipsometry software. Handles automatic eV↔nm conversion
and provides interpolation to arbitrary wavelength grids.

=head1 SUPPORTED FORMATS

=over

=item Woollam .mat format

3-line header:
  Line 1: material name
  Line 2: units (nm or eV)
  Line 3: number of data points
  Data: wavelength/energy  n  k

=item Generic 3-column format

Tab or space-separated: wavelength(nm)  n  k
No header, or lines starting with # are comments.

=back

=cut

sub load_material {
    my ($filepath) = @_;
    open my $fh, '<', $filepath or die "Cannot open material file $filepath: $!";
    my @lines = <$fh>;
    close $fh;

    # Strip Windows CR from all lines
    s/\r//g for @lines;
    chomp @lines;

    my ($name, $units, $npts);
    my @data;

    # Detect format by checking header
    if (@lines >= 3 && $lines[1] =~ /^\s*(nm|eV)\s*$/i) {
        # Woollam .mat format
        $name  = $lines[0];
        $units = lc($lines[1]);
        $units =~ s/\s+//g;
        $npts  = $lines[2] + 0 if $lines[2] =~ /^\d+/;

        for my $i (3 .. $#lines) {
            next if $lines[$i] =~ /^\s*$/;
            my @fields = split /\s+/, $lines[$i];
            next unless @fields >= 3 && $fields[0] =~ /^[-+]?\d/;
            push @data, [@fields[0..2]];
        }
    } else {
        # Generic 3-column format
        $name  = $filepath;
        $units = 'nm';
        for my $line (@lines) {
            next if $line =~ /^\s*#/;
            next if $line =~ /^\s*$/;
            my @fields = split /\s+/, $line;
            next unless @fields >= 3 && $fields[0] =~ /^[-+]?\d/;
            push @data, [@fields[0..2]];
        }
    }

    die "No data found in $filepath" unless @data;

    my $arr = pdl \@data;
    my $wav = $arr->(0,:)->flat->sever;
    my $n   = $arr->(1,:)->flat->sever;
    my $k   = $arr->(2,:)->flat->sever;

    # Convert eV to nm if needed
    if ($units eq 'ev') {
        $wav = 1239.842 / $wav;
        # Reverse if now in descending order
        if ($wav->at(0) > $wav->at(-1)) {
            $wav = $wav->(-1:0)->sever;
            $n   = $n->(-1:0)->sever;
            $k   = $k->(-1:0)->sever;
        }
    }

    return {
        name       => $name,
        wavelength => $wav,
        n          => $n,
        k          => $k,
        npts       => $wav->nelem,
        wav_min    => $wav->min->sclr,
        wav_max    => $wav->max->sclr,
    };
}

# Interpolate material optical constants to a given wavelength grid
sub interpolate_material {
    my ($material, $lambda_nm) = @_;
    my $n = $lambda_nm->interpol($material->{wavelength}, $material->{n});
    my $k = $lambda_nm->interpol($material->{wavelength}, $material->{k});
    return ($n, $k);
}

1;
