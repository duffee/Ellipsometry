package Physics::Ellipsometry::VASE::Dispersion;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use Exporter 'import';

our @EXPORT_OK = qw(cauchy_nk sellmeier_nk tauc_lorentz_nk drude_nk genosc_nk);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::Dispersion - Optical dispersion models

=head1 DESCRIPTION

Provides standard dispersion models for thin film optical constants.
All functions return (n, k) piddles as a function of wavelength (nm).

Note: Functions are called directly (not via closures) to avoid
PDL::NiceSlice source filter conflicts with $ref->() syntax.

=cut

# Cauchy dispersion: n(λ) = A + B/λ² + C/λ⁴, k=0 (or Urbach tail)
# $lambda in nm, B in µm², C in µm⁴
sub cauchy_nk {
    my ($lambda_nm, $A, $B, $C, %opts) = @_;
    $A //= 2.0;
    $B //= 0.01;
    $C //= 0.0;
    my $k_amp = $opts{k_amp} // 0;
    my $k_exp = $opts{k_exp} // 0;

    my $lam_um = $lambda_nm / 1000.0;
    my $n = $A + $B / $lam_um**2 + $C / $lam_um**4;
    my $k;
    if ($k_amp > 0) {
        $k = $k_amp * exp($k_exp * (1.0/$lam_um - 1.0/0.4));
    } else {
        $k = zeroes($lambda_nm);
    }
    return ($n, $k);
}

# Sellmeier dispersion: n²(λ) = 1 + Σ Bᵢλ²/(λ² - Cᵢ)
# $B_terms, $C_terms are arrayrefs; C in µm²
sub sellmeier_nk {
    my ($lambda_nm, $B_terms, $C_terms) = @_;
    $B_terms //= [1.0];
    $C_terms //= [0.01];

    my $lam_um_sq = ($lambda_nm / 1000.0)**2;
    my $n_sq = ones($lambda_nm) + 0.0;
    for my $i (0 .. $#$B_terms) {
        $n_sq += $B_terms->[$i] * $lam_um_sq / ($lam_um_sq - $C_terms->[$i]);
    }
    $n_sq = $n_sq->clip(0.01, 1e6);
    my $n = sqrt($n_sq);
    my $k = zeroes($lambda_nm);
    return ($n, $k);
}

# Tauc-Lorentz dispersion (Jellison & Modine, 1996)
sub tauc_lorentz_nk {
    my ($lambda_nm, $A, $E0, $Gamma, $Eg, $eps_inf) = @_;
    $A       //= 100;
    $E0      //= 4.0;
    $Gamma   //= 1.0;
    $Eg      //= 3.5;
    $eps_inf //= 1.0;

    my $E = 1239.842 / $lambda_nm;

    # ε₂(E): Tauc-Lorentz imaginary part
    my $eps2 = zeroes($lambda_nm);
    my $above_gap = which($E > $Eg);
    if ($above_gap->nelem > 0) {
        my $Ea = $E->index($above_gap);
        my $numer = $A * $E0 * $Gamma * ($Ea - $Eg)**2;
        my $denom = (($Ea**2 - $E0**2)**2 + $Gamma**2 * $Ea**2) * $Ea;
        $eps2->index($above_gap) .= $numer / $denom;
    }

    # ε₁(E): numerical Kramers-Kronig
    my $eps1 = _kk_transform($E, $eps2) + $eps_inf;

    my $eps = $eps1 + i() * $eps2;
    my $N = sqrt($eps);
    return ($N->re, $N->im->abs);
}

# Drude model for metals: ε(E) = ε_inf - ωp²/(E² + iEΓ)
sub drude_nk {
    my ($lambda_nm, $eps_inf, $omega_p, $gamma) = @_;
    $eps_inf //= 1.0;
    $omega_p //= 10.0;
    $gamma   //= 0.1;

    my $E = 1239.842 / $lambda_nm;
    my $eps = $eps_inf - $omega_p**2 / ($E**2 + i() * $E * $gamma);
    my $N = sqrt($eps);
    return ($N->re, $N->im->abs);
}

# General oscillator: sum of Lorentz oscillators
# $oscillators: arrayref of [A, E0, Gamma] triplets
sub genosc_nk {
    my ($lambda_nm, $oscillators, $eps_inf) = @_;
    $oscillators //= [];
    $eps_inf //= 1.0;

    my $E = 1239.842 / $lambda_nm;
    my $eps = ones($lambda_nm) * $eps_inf + i() * zeroes($lambda_nm);

    for my $osc (@$oscillators) {
        my ($Ai, $E0i, $Gi) = @$osc;
        $eps += $Ai * $E0i**2 / ($E0i**2 - $E**2 - i() * $Gi * $E);
    }

    my $N = sqrt($eps);
    return ($N->re, $N->im->abs);
}

# Numerical Kramers-Kronig transform
sub _kk_transform {
    my ($E, $eps2) = @_;
    my $npts = $E->nelem;
    my $eps1 = zeroes($npts);

    for my $j (0 .. $npts - 1) {
        my $Ej = $E->at($j);
        my $mask = abs($E - $Ej) > 0.001;
        my $idx = which($mask);
        next if $idx->nelem < 2;
        my $E_prime = $E->index($idx);
        my $eps2_p  = $eps2->index($idx);
        my $integrand = $E_prime * $eps2_p / ($E_prime**2 - $Ej**2);
        my $dE = $E_prime->(1:) - $E_prime->(0:-2);
        my $avg = ($integrand->(1:) + $integrand->(0:-2)) / 2;
        $eps1->set($j, (2.0 / PI) * sum($dE * $avg)->sclr);
    }
    return $eps1;
}

1;
