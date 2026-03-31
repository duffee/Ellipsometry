package Physics::Ellipsometry::VASE;
use strict;
use warnings;
use PDL;
use PDL::Fit::LM;
use PDL::NiceSlice;

our $VERSION = '0.01';

# Constructor
sub new {
    my ($class, %args) = @_;
    my $self = {
        layers => $args{layers} || 1,
        model  => $args{model}  || sub { die "Model not implemented" },
        data   => undef,
    };
    bless $self, $class;
    return $self;
}

# Load data from file (auto-detects Woollam VASE format)
sub load_data {
    my ($self, $filename) = @_;
    open my $fh, '<', $filename or die "Cannot open $filename: $!";
    my @lines = <$fh>;
    close $fh;

    # Detect Woollam VASE format: line 2 starts with VASEmethod[
    my $is_woollam = (@lines >= 4 && $lines[1] =~ /^VASEmethod\[/);

    my @data;
    if ($is_woollam) {
        # Parse Woollam header
        chomp(my $sample_name = $lines[0]);
        $self->{sample_name} = $sample_name;

        if ($lines[1] =~ /^VASEmethod\[(.+)\]/) {
            $self->{vase_method} = $1;
        }

        if ($lines[2] =~ /^Original\[(.+)\]/) {
            $self->{original_file} = $1;
        }

        chomp(my $units = $lines[3]);
        $self->{units} = $units;

        # Data starts at line 5 (index 4)
        for my $i (4 .. $#lines) {
            my $line = $lines[$i];
            chomp $line;
            next if $line =~ /^\s*$/;
            my @fields = split ' ', $line;
            # Skip lines whose first field is not numeric (e.g. Mueller labels)
            next unless $fields[0] =~ /^[-+]?\d*\.?\d+([eE][-+]?\d+)?$/;
            push @data, \@fields;
        }
    } else {
        # Original format: skip #comments and blank lines
        for my $line (@lines) {
            next if $line =~ /^\s*#/;   # skip comment lines
            next if $line =~ /^\s*$/;    # skip blank lines
            my @fields = split ' ', $line;
            push @data, \@fields;
        }
    }

    my $data = pdl \@data;

    # Woollam files have 6 columns; extract sigma and keep standard 4
    if ($is_woollam && $data->getdim(0) >= 6) {
        $self->{sigma} = $data->(4:5,:);
        $data = $data->(0:3,:)->sever;
    }

    $self->{data} = $data;
    return $data;
}

# Set model function
sub set_model {
    my ($self, $model) = @_;
    $self->{model} = $model;
}

# Fit data to model
sub fit {
    my ($self, $initial_params) = @_;
    my $data = $self->{data};
    my $model = $self->{model};

    # data shape is (nfields, npts); dim0=fields, dim1=data points
    # Extract x as (npts, 2) so $x->(:,0)=wavelength, $x->(:,1)=angle
    my $x_data = $data->(0:1,:)->xchg(0,1);

    # Build y to match model output order: [psi_all, delta_all]
    my $y_data = $data->((2),:)->flat->append($data->((3),:)->flat);

    # Use measured uncertainties from Woollam data when available
    my $sigma;
    if (defined $self->{sigma}) {
        $sigma = $self->{sigma}->((0),:)->flat->append($self->{sigma}->((1),:)->flat);
    } else {
        $sigma = ones($y_data->nelem);
    }

    # lmfit sizes $dyda from $x->getdim(0), so pass a dummy x whose
    # first dimension equals the number of y values (2*npts)
    my $x_fit = sequence($y_data->nelem);

    # Wrapper: adapts user model to lmfit interface ($x, $par, $ym, $dyda)
    my $fit_func = sub {
        my ($x, $par, $ym, $dyda) = @_;

        my $y_model = &$model($par, $x_data);
        $ym .= $y_model;

        # Numerical partial derivatives via finite differences
        my $np  = $par->nelem;
        my $eps = 1e-7;
        for my $i (0 .. $np - 1) {
            my $par_h = $par->copy;
            $par_h->slice("($i)") += $eps;
            $dyda->slice(",($i)") .= (&$model($par_h, $x_data) - $y_model) / $eps;
        }
    };

    my ($ym, $finalp, $covar, $iters) = lmfit(
        $x_fit, $y_data, $sigma, $fit_func, $initial_params,
        {Maxiter => 300, Eps => 1e-7}
    );

    return $finalp;
}

# Plot raw data with model fit overlay
sub plot {
    my ($self, $fit_params, %opts) = @_;
    require PDL::Graphics::Gnuplot;

    my $data  = $self->{data};
    my $model = $self->{model};

    my $wavelength = $data->((0),:)->flat;
    my $psi_data   = $data->((2),:)->flat;
    my $delta_data = $data->((3),:)->flat;

    # Evaluate model at fitted parameters
    my $x_data  = $data->(0:1,:)->xchg(0,1);
    my $y_model = &$model($fit_params, $x_data);
    my $npts    = $wavelength->nelem;
    my $psi_fit   = $y_model->slice("0:" . ($npts - 1));
    my $delta_fit = $y_model->slice("$npts:" . (2 * $npts - 1));

    my $output = $opts{output};
    my $title  = $opts{title} // 'VASE Fit Results';

    # Select terminal and construct gpwin
    my $gp;
    if ($output) {
        my ($term, @topts);
        if    ($output =~ /\.png$/i) { $term = "pngcairo"; @topts = (size => [900,700,"px"]) }
        elsif ($output =~ /\.pdf$/i) { $term = "pdfcairo"; @topts = (size => [7,5.5,"in"]) }
        elsif ($output =~ /\.svg$/i) { $term = "svg";      @topts = (size => [900,700,"px"]) }
        elsif ($output =~ /\.eps$/i) { $term = "epscairo" }
        else                         { $term = "pngcairo"; @topts = (size => [900,700,"px"]) }
        $gp = PDL::Graphics::Gnuplot::gpwin($term, output => $output, enhanced => 1, @topts);
    } else {
        $gp = PDL::Graphics::Gnuplot::gpwin(enhanced => 1);
    }

    # Multiplot: Psi on top, Delta on bottom (rows, cols)
    $gp->multiplot(layout => [1, 2], title => $title);

    # --- Psi panel ---
    $gp->plot(
        { title  => '{/Symbol Y} (Psi)',
          xlabel => '',
          ylabel => '{/Symbol Y} (deg)',
          key    => 'top right box' },
        with => 'points', legend => 'Data',
            pt => 7, ps => 1.2, lc => 'rgb "#0072B2"',
            $wavelength, $psi_data,
        with => 'lines', legend => 'Fit',
            lw => 2, dt => 1, lc => 'rgb "#D55E00"',
            $wavelength, $psi_fit,
    );

    # --- Delta panel ---
    $gp->plot(
        { title  => '{/Symbol D} (Delta)',
          xlabel => 'Wavelength (nm)',
          ylabel => '{/Symbol D} (deg)',
          key    => 'top left box' },
        with => 'points', legend => 'Data',
            pt => 7, ps => 1.2, lc => 'rgb "#0072B2"',
            $wavelength, $delta_data,
        with => 'lines', legend => 'Fit',
            lw => 2, dt => 1, lc => 'rgb "#D55E00"',
            $wavelength, $delta_fit,
    );

    $gp->end_multi;

    if ($output) {
        $gp->close;
        print "Plot saved to $output\n";
    }

    return $gp;
}

1;

__END__

=head1 NAME

Physics::Ellipsometry::VASE - Variable Angle Spectroscopic Ellipsometry data fitting

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

    use PDL;
    use PDL::NiceSlice;
    use Physics::Ellipsometry::VASE;

    my $vase = Physics::Ellipsometry::VASE->new(layers => 1);
    $vase->load_data('data.dat');

    sub my_model {
        my ($params, $x) = @_;
        my $a = $params->(0);
        my $b = $params->(1);
        my $c = $params->(2);
        my $d = $params->(3);
        my $wavelength = $x->(:,0);

        my $psi   = $a - $b * $wavelength;
        my $delta = $c + $d * $wavelength;

        return cat($psi, $delta)->flat;
    }

    $vase->set_model(\&my_model);

    my $fit_params = $vase->fit(pdl [65, 0.05, 80, 0.1]);
    print "Fitted: $fit_params\n";

=head1 DESCRIPTION

Physics::Ellipsometry::VASE provides Levenberg-Marquardt fitting of
user-defined optical models to variable angle spectroscopic ellipsometry
(VASE) data.  It wraps L<PDL::Fit::LM> and handles the bookkeeping of
mapping between the simple user model interface and the lmfit calling
convention, including automatic numerical computation of the Jacobian
via finite differences.

Ellipsometry measures the change in polarization state of light reflected
from a surface.  The two measured quantities are B<Psi> (amplitude ratio)
and B<Delta> (phase difference).  This module fits models that predict
both Psi and Delta simultaneously as a function of wavelength and
angle of incidence.

=head1 DATA FORMAT

Input files should contain whitespace-separated columns:

    # Wavelength(nm)  Angle(deg)  Psi(deg)  Delta(deg)
    400  70  45.0  120.0
    410  70  44.5  121.0

Lines beginning with C<#> and blank lines are skipped.

=head1 MODEL FUNCTIONS

A model function receives two arguments and must return a single
flattened PDL piddle:

    sub model {
        my ($params, $x) = @_;

        # $params: PDL piddle of fit parameters
        # $x:      PDL piddle of shape (npoints, 2)
        #          $x->(:,0) = wavelength (nm)
        #          $x->(:,1) = angle of incidence (deg)

        my $psi   = ...;   # compute Psi  (npoints)
        my $delta = ...;   # compute Delta (npoints)

        return cat($psi, $delta)->flat;
    }

The Jacobian (partial derivatives with respect to parameters) is
computed automatically via numerical finite differences.

=head1 METHODS

=head2 new

    my $vase = Physics::Ellipsometry::VASE->new(%args);

Constructor.  Accepts:

=over 4

=item layers

Number of layers in the optical model (default: 1).

=item model

Optional code reference to a model function.

=back

=head2 load_data

    $vase->load_data($filename);

Reads ellipsometry data from a whitespace-delimited file.
Returns the loaded data as a PDL piddle.

=head2 set_model

    $vase->set_model(\&model_func);

Sets the model function used for fitting.

=head2 fit

    my $fitted_params = $vase->fit($initial_params);

Performs Levenberg-Marquardt fitting.  C<$initial_params> is a PDL
piddle of initial guesses.  Returns a PDL piddle of fitted parameters.

=head2 plot

    $vase->plot($fit_params);
    $vase->plot($fit_params, output => 'fit.png');
    $vase->plot($fit_params, output => 'fit.pdf', title => 'My Fit');

Plots raw data points with model fit overlay in a two-panel layout
(Psi on top, Delta on bottom).  Requires L<PDL::Graphics::Gnuplot>.

Options:

=over 4

=item output

File path for saving the plot.  Format is inferred from the extension
(C<.png>, C<.pdf>, C<.svg>, C<.eps>).  If omitted, displays an
interactive window.

=item title

Overall plot title (default: C<VASE Fit Results>).

=back

=head1 DEPENDENCIES

L<PDL>, L<PDL::Fit::LM>, L<PDL::NiceSlice>

L<PDL::Graphics::Gnuplot> is required only for the C<plot> method.

=head1 AUTHOR

jtrujil1

=head1 LICENSE AND COPYRIGHT

This software is copyright (c) 2026.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
