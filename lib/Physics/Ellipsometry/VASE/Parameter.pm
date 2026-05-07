package Physics::Ellipsometry::VASE::Parameter;
use strict;
use warnings;
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(param params_to_pdl pdl_to_params make_fit_model get_values);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::Parameter - Parameter bounds and vary/fix control

=head1 DESCRIPTION

Provides a parameter object for ellipsometry model fitting with:
- Named parameters with values
- Min/max bounds (enforced via transformation)
- Vary/fix control (fixed parameters excluded from fit)
- Scaling for numerical stability

Parameters with bounds are internally transformed using a logit-like mapping
to enforce constraints without discontinuous derivatives.

=cut

sub param {
    my (%args) = @_;
    return {
        name   => $args{name}  // 'unnamed',
        value  => $args{value} // 0.0,
        min    => $args{min},           # undef = unbounded below
        max    => $args{max},           # undef = unbounded above
        vary   => $args{vary}  // 1,    # 1=vary, 0=fixed
        scale  => $args{scale} // 1.0,  # internal scaling factor
    };
}

# Convert parameter list to internal PDL for fitting (only varying params)
sub params_to_pdl {
    my ($params) = @_;  # arrayref of param hashes
    my @values;
    for my $p (@$params) {
        next unless $p->{vary};
        my $internal = _value_to_internal($p);
        push @values, $internal;
    }
    return pdl(\@values);
}

# Update parameter list from fitted PDL values
sub pdl_to_params {
    my ($params, $fitted_pdl) = @_;
    my $idx = 0;
    for my $p (@$params) {
        next unless $p->{vary};
        my $internal = $fitted_pdl->at($idx);
        $p->{value} = _internal_to_value($p, $internal);
        $idx++;
    }
    return $params;
}

# Get all parameter values as a list (both fixed and varying)
sub get_values {
    my ($params) = @_;
    return map { $_->{value} } @$params;
}

# Create a model wrapper that handles parameter transformation
# Returns a closure compatible with VASE fit()
sub make_fit_model {
    my ($params, $full_model) = @_;
    # $full_model receives ($all_values_pdl, $x_data)
    # Returns wrapper that receives ($vary_pdl, $x_data)

    return sub {
        my ($vary_pdl, $x_data) = @_;
        # Reconstruct full parameter vector
        my @all_values;
        my $vary_idx = 0;
        for my $p (@$params) {
            if ($p->{vary}) {
                my $internal = $vary_pdl->at($vary_idx);
                push @all_values, _internal_to_value($p, $internal);
                $vary_idx++;
            } else {
                push @all_values, $p->{value};
            }
        }
        my $all_pdl = pdl(\@all_values);
        return &$full_model($all_pdl, $x_data);
    };
}

# Transform value to internal (unbounded) space for fitting
sub _value_to_internal {
    my ($p) = @_;
    my $v = $p->{value};

    if (defined $p->{min} && defined $p->{max}) {
        # Bounded: use logit transform
        my $range = $p->{max} - $p->{min};
        my $norm = ($v - $p->{min}) / $range;
        # Clamp to avoid log(0)
        $norm = 0.001 if $norm < 0.001;
        $norm = 0.999 if $norm > 0.999;
        return log($norm / (1.0 - $norm)) * $p->{scale};
    } elsif (defined $p->{min}) {
        # Lower bounded: use log transform
        my $shifted = $v - $p->{min};
        $shifted = 0.001 if $shifted < 0.001;
        return log($shifted) * $p->{scale};
    } elsif (defined $p->{max}) {
        # Upper bounded: use negative log transform
        my $shifted = $p->{max} - $v;
        $shifted = 0.001 if $shifted < 0.001;
        return -log($shifted) * $p->{scale};
    } else {
        # Unbounded
        return $v * $p->{scale};
    }
}

# Transform from internal (unbounded) space back to actual value
sub _internal_to_value {
    my ($p, $internal) = @_;
    my $x = $internal / ($p->{scale} || 1.0);

    if (defined $p->{min} && defined $p->{max}) {
        # Inverse logit (sigmoid)
        my $range = $p->{max} - $p->{min};
        my $sigmoid = 1.0 / (1.0 + exp(-$x));
        return $p->{min} + $range * $sigmoid;
    } elsif (defined $p->{min}) {
        # Inverse log
        return $p->{min} + exp($x);
    } elsif (defined $p->{max}) {
        # Inverse negative log
        return $p->{max} - exp(-$x);
    } else {
        # Unbounded
        return $x;
    }
}

1;
