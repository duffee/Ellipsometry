package Physics::Ellipsometry::VASE::EMA;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Exporter 'import';

our @EXPORT_OK = qw(ema_linear ema_bruggeman ema_maxwell_garnett);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::EMA - Effective Medium Approximation models

=head1 DESCRIPTION

Provides EMA mixing rules for calculating effective optical constants
of composite layers (e.g., porous films, interfacial layers, surface roughness).

All functions take dielectric functions (ε = N²) of constituent materials
and a volume fraction, returning the effective dielectric function.

=cut

# Linear (volume-weighted) mixing: ε_eff = (1-f)ε_a + f·ε_b
sub ema_linear {
    my ($eps_a, $eps_b, $vf) = @_;
    return (1.0 - $vf) * $eps_a + $vf * $eps_b;
}

# Bruggeman EMA: f_a·(ε_a - ε_eff)/(ε_a + 2ε_eff) + f_b·(ε_b - ε_eff)/(ε_b + 2ε_eff) = 0
# Solved analytically for 2-component mixtures:
#   ε_eff = (b ± sqrt(b² + 8ε_aε_b)) / 4
#   where b = (3f_b - 1)ε_b + (3f_a - 1)ε_a = (3vf - 1)ε_b + (2 - 3vf)ε_a
sub ema_bruggeman {
    my ($eps_a, $eps_b, $vf) = @_;
    my $fa = 1.0 - $vf;

    # Quadratic solution for 2-component Bruggeman
    my $b = (3*$vf - 1) * $eps_b + (3*$fa - 1) * $eps_a;
    my $discriminant = $b**2 + 8 * $eps_a * $eps_b;

    # Take the root with positive real part
    my $sqrt_disc = sqrt($discriminant + i()*0);
    my $eps_eff = ($b + $sqrt_disc) / 4.0;

    # Ensure physical result (positive real part of ε)
    my $eps_eff2 = ($b - $sqrt_disc) / 4.0;
    my $use_alt = ($eps_eff->re < 0) & ($eps_eff2->re > 0);
    if (ref $use_alt && $use_alt->any) {
        my $idx = which($use_alt);
        $eps_eff->index($idx) .= $eps_eff2->index($idx);
    }

    return $eps_eff;
}

# Maxwell-Garnett EMA: inclusion (b) in host matrix (a)
# ε_eff = ε_a · (ε_b + 2ε_a + 2f(ε_b - ε_a)) / (ε_b + 2ε_a - f(ε_b - ε_a))
sub ema_maxwell_garnett {
    my ($eps_a, $eps_b, $vf) = @_;

    my $numer = $eps_b + 2*$eps_a + 2*$vf*($eps_b - $eps_a);
    my $denom = $eps_b + 2*$eps_a -   $vf*($eps_b - $eps_a);

    return $eps_a * $numer / $denom;
}

1;
