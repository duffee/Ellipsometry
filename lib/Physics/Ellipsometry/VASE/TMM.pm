package Physics::Ellipsometry::VASE::TMM;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use Exporter 'import';

our @EXPORT_OK = qw(tmm_reflect psi_delta);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::TMM - Transfer Matrix Method for multilayer optics

=head1 DESCRIPTION

Implements the 2×2 transfer matrix method for calculating reflection
coefficients of multilayer thin film stacks. Uses the physics sign
convention (e^{-iωt} time dependence → e^{+2iβ} phase propagation).

Supports arbitrary layer stacks with complex refractive indices.

=head1 CONVENTIONS

=over

=item Phase: Physics convention (e^{+2iβ} for forward propagation)

=item Fresnel rp: Verdet convention: rp = (Nf·cosθi - Ni·cosθf)/(Nf·cosθi + Ni·cosθf)

=item Delta: -arg(ρ) mapped to [0°, 360°) matching WVASE/refellips

=back

=cut

# Calculate reflection coefficients for a multilayer stack
# Input:
#   $lambda_nm - wavelength piddle (npts)
#   $theta_deg - angle of incidence piddle (npts), in degrees
#   $N_layers  - arrayref of complex refractive index piddles [N0, N1, ..., Ns]
#                N0 = ambient, Ns = substrate
#   $d_nm      - arrayref of layer thicknesses in nm [d1, d2, ..., d_{s-1}]
#                (no thickness for ambient or substrate)
# Returns: ($rp, $rs) complex reflection coefficients
sub tmm_reflect {
    my ($lambda_nm, $theta_deg, $N_layers, $d_nm) = @_;

    my $n_media = scalar @$N_layers;  # number of media (ambient + layers + substrate)
    die "Need at least 2 media (ambient + substrate)" if $n_media < 2;
    die "Need exactly N-2 thicknesses" if scalar @$d_nm != $n_media - 2;

    my $theta_rad = $theta_deg * (PI / 180.0);
    my $cos_t0 = cos($theta_rad);
    my $sin_t0 = sin($theta_rad);
    my $N0 = $N_layers->[0];

    # Calculate cos(θ) in each layer via Snell's law
    my @cos_t;
    push @cos_t, $cos_t0;
    for my $j (1 .. $n_media - 1) {
        my $Nj = $N_layers->[$j];
        $cos_t[$j] = sqrt(1.0 - ($N0 * $sin_t0)**2 / $Nj**2);
    }

    # Build system by nested Airy formula (back-to-front)
    # Start from the deepest interface and propagate backward
    my $last = $n_media - 1;

    # Fresnel coefficients at last interface (layer n-2 → substrate)
    my $rs = _fresnel_s($N_layers->[$last-1], $N_layers->[$last],
                        $cos_t[$last-1], $cos_t[$last]);
    my $rp = _fresnel_p($N_layers->[$last-1], $N_layers->[$last],
                        $cos_t[$last-1], $cos_t[$last]);

    # Propagate backward through each layer
    for (my $j = $last - 2; $j >= 0; $j--) {
        # Phase thickness of layer j+1
        my $d = $d_nm->[$j];  # thickness of layer j+1 (0-indexed in d_nm)
        my $beta = (2 * PI / $lambda_nm) * $N_layers->[$j+1] * $d * $cos_t[$j+1];

        # Fresnel at interface j → j+1
        my $r_s_ij = _fresnel_s($N_layers->[$j], $N_layers->[$j+1],
                                $cos_t[$j], $cos_t[$j+1]);
        my $r_p_ij = _fresnel_p($N_layers->[$j], $N_layers->[$j+1],
                                $cos_t[$j], $cos_t[$j+1]);

        # Airy formula: r = (r_ij + r_below·e^{+2iβ}) / (1 + r_ij·r_below·e^{+2iβ})
        my $phase = exp(2.0 * i() * $beta);
        $rs = ($r_s_ij + $rs * $phase) / (1.0 + $r_s_ij * $rs * $phase);
        $rp = ($r_p_ij + $rp * $phase) / (1.0 + $r_p_ij * $rp * $phase);
    }

    return ($rp, $rs);
}

# Calculate Psi and Delta from a layer stack
# Returns: ($psi_deg, $delta_deg) with Delta in [0, 360)
sub psi_delta {
    my ($lambda_nm, $theta_deg, $N_layers, $d_nm, %opts) = @_;

    my ($rp, $rs) = tmm_reflect($lambda_nm, $theta_deg, $N_layers, $d_nm);

    my $rho = $rp / $rs;
    my $psi = atan(abs($rho)) * (180.0 / PI);

    # Delta = -arg(ρ) mapped to [0°, 360°)
    my $delta_rad = -carg($rho)->re;
    $delta_rad += 2*PI * ($delta_rad < 0);
    my $delta = $delta_rad * (180.0 / PI);

    # Optionally align to reference data (avoid 0/360 wrap)
    if (my $delta_ref = $opts{delta_ref}) {
        my $diff = $delta - $delta_ref;
        $delta -= 360.0 * rint($diff / 360.0);
    }

    return ($psi->re->double, $delta->double);
}

# Fresnel s-polarization: rs = (Ni·cosθi - Nf·cosθf) / (Ni·cosθi + Nf·cosθf)
sub _fresnel_s {
    my ($Ni, $Nf, $cos_ti, $cos_tf) = @_;
    return ($Ni*$cos_ti - $Nf*$cos_tf) / ($Ni*$cos_ti + $Nf*$cos_tf);
}

# Fresnel p-polarization (Verdet): rp = (Nf·cosθi - Ni·cosθf) / (Nf·cosθi + Ni·cosθf)
sub _fresnel_p {
    my ($Ni, $Nf, $cos_ti, $cos_tf) = @_;
    return ($Nf*$cos_ti - $Ni*$cos_tf) / ($Nf*$cos_ti + $Ni*$cos_tf);
}

1;
