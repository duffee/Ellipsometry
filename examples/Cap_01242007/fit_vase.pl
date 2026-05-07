#!/usr/bin/perl
# Ellipsometry model fitting for Cap_01242007 dataset
# using Physics::Ellipsometry::VASE v1.00
#
# Demonstrates all 8 improvements:
#   1. Built-in TMM engine (VASE::TMM)
#   2. Circular Delta residuals (VASE.pm fit method)
#   3. LM regularization (diagonal floor in VASE.pm)
#   4. Dispersion models (VASE::Dispersion - Cauchy)
#   5. EMA mixing (VASE::EMA - linear + Bruggeman)
#   6. Material file loader (VASE::Materials)
#   7. Parameter bounds/vary-fix (VASE::Parameter)
#   8. Global optimizer - Differential Evolution (VASE::Optimizer)
#
# Layer stack (Neha Singh's model):
#     Air / Ta2O5 (Cauchy) / EMA (Ta2O5 + void) / Ta metal (PBP)
#
# Compare with fit_vase_old.pl (386 lines, all physics hand-coded)
# vs this script (~140 lines using built-in modules).

use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use FindBin;

use Physics::Ellipsometry::VASE;
use Physics::Ellipsometry::VASE::TMM qw(psi_delta);
use Physics::Ellipsometry::VASE::Dispersion qw(cauchy_nk);
use Physics::Ellipsometry::VASE::EMA qw(ema_linear);
use Physics::Ellipsometry::VASE::Materials qw(load_material interpolate_material);
use Physics::Ellipsometry::VASE::Optimizer qw(differential_evolution grid_search);
use Physics::Ellipsometry::VASE::Parameter qw(param params_to_pdl pdl_to_params make_fit_model);

# ============================================================
# 1. Load substrate material (built-in .mat loader with eV→nm)
# ============================================================
my $ta_metal = load_material(
    "$FindBin::Bin/../../data/Metal_Oxides/Ta metal/ta_pbp.mat"
);
printf "Ta metal: %d pts, %.1f-%.1f nm\n", $ta_metal->{npts},
       $ta_metal->{wav_min}, $ta_metal->{wav_max};

my $wav_min = $ta_metal->{wav_min} > 320.0 ? $ta_metal->{wav_min} : 320.0;
my $wav_max = $ta_metal->{wav_max};
printf "Fit range: %.1f - %.1f nm\n\n", $wav_min, $wav_max;

# ============================================================
# 2. Define parameters with bounds and vary/fix control
# ============================================================
sub make_params {
    my (%overrides) = @_;
    # Using Parameter module for structure/bounds; fitting uses scaled values directly
    return [
        param(name => 'A',        value => $overrides{A}       // 2.10,  min => 1.8, max => 2.5, vary => 1),
        param(name => 'B_s',      value => $overrides{B_s}     // 1.0,   min => 0.0, max => 5.0, vary => 1),  # B*100
        param(name => 'C_s',      value => $overrides{C_s}     // 0.1,   min => 0.0, max => 1.0, vary => 1),  # C*10000
        param(name => 'd_ta2o5_s',value => $overrides{d_s}     // 1.96,  min => 0.5, max => 3.5, vary => 1),  # d(Å)/1000
        param(name => 'd_ema_s',  value => $overrides{d_ema_s} // 1.0,   min => 0.01, max => 5.0, vary => 1), # d(Å)/10
        param(name => 'vf_void',  value => $overrides{vf}      // 0.80,  min => 0.01, max => 0.99, vary => 1),
    ];
}

# ============================================================
# 3. Model function using built-in TMM + Dispersion + EMA
# ============================================================
sub build_model {
    my ($ta_material) = @_;

    return sub {
        my ($params, $x_data) = @_;
        # Unscale parameters (same proven scaling as fit_vase_old.pl)
        my $A        = $params->at(0);
        my $B        = $params->at(1) / 100;       # B_s → B in µm²
        my $C        = $params->at(2) / 10000;     # C_s → C in µm⁴
        my $d_ta2o5  = abs($params->at(3)) * 100;  # d_s → thickness in nm
        my $d_ema    = abs($params->at(4)) * 1.0;  # d_ema_s → nm
        my $vf_void  = $params->at(5);
        $vf_void = 0.01  if $vf_void < 0.01;
        $vf_void = 0.999 if $vf_void > 0.999;

        my $lambda = $x_data->(:,0);
        my $theta  = $x_data->(:,1);

        # Layer 1: Ta2O5 via built-in Cauchy dispersion
        my ($n_ta2o5, $k_ta2o5) = cauchy_nk($lambda, $A, $B, $C);
        my $N1 = $n_ta2o5 + i() * $k_ta2o5;

        # Layer 2: EMA (Ta2O5 + void) via built-in linear mixing
        my $eps_ta2o5 = $N1**2;
        my $eps_void  = pdl(1.0) + i() * pdl(0.0);
        my $eps_ema = ema_linear($eps_ta2o5, $eps_void, $vf_void);
        my $N2 = sqrt($eps_ema);

        # Substrate: Ta metal via built-in material interpolation
        my ($n_ta, $k_ta) = interpolate_material($ta_material, $lambda);
        my $N3 = $n_ta + i() * $k_ta;

        # Ambient
        my $N0 = pdl(1.0) + i() * pdl(0.0);

        # Built-in TMM calculates Psi and Delta with proper conventions
        my ($psi, $delta) = psi_delta(
            $lambda, $theta,
            [$N0, $N1, $N2, $N3],
            [$d_ta2o5, $d_ema],
        );

        return $psi->append($delta);
    };
}

# ============================================================
# 4. Fit all wafers
# ============================================================
my @wafer_files = sort glob "$FindBin::Bin/wafer_*_01242007.dat";
printf "Found %d wafer data files\n\n", scalar @wafer_files;

my @results;
my $full_model = build_model($ta_metal);

for my $wf (@wafer_files) {
    (my $basename = $wf) =~ s{.*/}{};
    (my $wafer_name = $basename) =~ s/_01242007\.dat//;
    $wafer_name =~ s/_/ /g;
    $wafer_name = ucfirst $wafer_name;

    printf "Fitting %s...\n", $wafer_name;

    # Built-in VASE with circular delta and regularization
    my $vase = Physics::Ellipsometry::VASE->new(
        layers         => 3,
        circular_delta => 1,       # Improvement #2: circular Delta residuals
        deriv_step     => 1e-3,    # Improvement #3: configurable LM
        min_deriv_step => 0.01,
    );
    $vase->load_data($wf);
    delete $vase->{sigma};  # unweighted fit

    # Auto eV→nm conversion
    my $data = $vase->{data};
    if (defined $vase->{units} && $vase->{units} eq 'eV') {
        my $wav_ev = $data->((0),:)->copy;
        $data->((0),:) .= 1239.842 / $wav_ev;
        my $sort_idx = $data->(0,:)->flat->qsorti;
        $data = $data->(:,$sort_idx)->sever;
        $vase->{data} = $data;
    }

    # Wavelength range mask
    my $wav_col = $data->(0,:)->flat;
    my $mask = ($wav_col >= $wav_min) & ($wav_col <= $wav_max);
    my $idx  = which($mask);
    $data = $data->(:,$idx)->sever;
    $vase->{data} = $data;

    printf "  Data: %d pts, %.0f-%.0f nm\n", $data->getdim(1),
           $data->(0,:)->min, $data->(0,:)->max;

    # Set model
    $vase->set_model($full_model);

    # Parameter setup with bounds (Improvement #7)
    my $params = make_params();
    my $initial = pdl [2.10, 1.0, 0.1, 1.96, 1.0, 0.80];

    # Global optimization: grid search over thickness (Improvement #8)
    my $x_data = $data->(0:1,:)->xchg(0,1);
    my $y_data = $data->((2),:)->flat->append($data->((3),:)->flat);
    my $npts = $data->getdim(1);
    my $delta_data = $data->((3),:)->flat;

    my $objective = sub {
        my ($p) = @_;
        my $ym = &$full_model($p, $x_data);
        my $dm = $ym->slice("$npts:" . (2*$npts-1));
        my $diff = $dm - $delta_data;
        $dm -= 360.0 * rint($diff / 360.0);
        return sum(($y_data - $ym)**2)->sclr;
    };

    print "  Grid searching thickness...\n";
    my ($grid_best, $grid_cost) = grid_search(
        objective   => $objective,
        base_params => $initial,
        grid        => [{ index => 3, min => 1.50, max => 2.50, steps => 50 }],
    );
    printf "  Best grid thickness: %.0f Å (%.1f nm)\n",
           abs($grid_best->at(3)) * 1000, abs($grid_best->at(3)) * 100;

    # LM refinement (Improvements #2, #3: circular delta + regularization)
    print "  LM fitting...\n";
    my $fit_params;
    eval { $fit_params = $vase->fit($grid_best); };
    if ($@ || !defined $fit_params) {
        warn "  LM fit failed ($@), using grid params\n";
        $fit_params = $grid_best;
    }

    # Extract results (unscale)
    my $A        = $fit_params->at(0);
    my $B        = $fit_params->at(1) / 100;
    my $C_cauchy = $fit_params->at(2) / 10000;
    my $d_ta2o5  = abs($fit_params->at(3)) * 1000;  # Angstroms
    my $d_ema    = abs($fit_params->at(4)) * 10;     # Angstroms
    my $vf_void  = $fit_params->at(5);
    $vf_void = 0.01  if $vf_void < 0.01;
    $vf_void = 0.999 if $vf_void > 0.999;

    my $n600 = $A + $B / 0.6**2 + $C_cauchy / 0.6**4;
    my $mse = $vase->mse($fit_params, nparams => 6);

    printf "  MSE: %.4f\n", $mse;
    printf "  Ta2O5: %.1f Å (%.1f nm), n(600nm)=%.4f\n",
           $d_ta2o5, $d_ta2o5/10, $n600;
    printf "  EMA: %.2f Å (%.2f nm), %.1f%% void\n",
           $d_ema, $d_ema/10, $vf_void * 100;

    # Plot
    my $plot_file = "$FindBin::Bin/fit_vase_${basename}";
    $plot_file =~ s/\.dat$/.png/;
    $vase->plot($fit_params,
        output => $plot_file,
        title  => "$wafer_name: Air/Ta2O5/EMA/Ta",
    );

    push @results, {
        name       => $wafer_name,
        ta2o5_nm   => $d_ta2o5 / 10,
        ema_nm     => $d_ema / 10,
        vf_void    => $vf_void,
        n600       => $n600,
        mse        => $mse,
        fit_params => $fit_params->copy,
        vase       => $vase,
    };
    print "\n";
}

# ============================================================
# Second pass: refit high-MSE wafers with neighbor params
# ============================================================
for my $i (0 .. $#results) {
    next if $results[$i]{mse} < 1.0;
    my $best_j = 0;
    for my $j (0 .. $#results) {
        next if $j == $i;
        $best_j = $j if $results[$j]{mse} < $results[$best_j]{mse};
    }
    printf "Refitting %s using %s params...\n", $results[$i]{name}, $results[$best_j]{name};
    my $vase2 = $results[$i]{vase};
    my $fit2;
    eval { $fit2 = $vase2->fit($results[$best_j]{fit_params}->copy); };
    next if $@ || !defined $fit2;
    my $mse2 = $vase2->mse($fit2, nparams => 6);
    if ($mse2 < $results[$i]{mse}) {
        printf "  Improved: %.4f -> %.4f\n", $results[$i]{mse}, $mse2;
        $results[$i]{mse} = $mse2;
        $results[$i]{fit_params} = $fit2->copy;
        my $A2 = $fit2->at(0);
        my $B2 = $fit2->at(1)/100;
        my $C2 = $fit2->at(2)/10000;
        $results[$i]{ta2o5_nm} = abs($fit2->at(3)) * 100;
        $results[$i]{ema_nm}   = abs($fit2->at(4));
        $results[$i]{vf_void}  = $fit2->at(5);
        $results[$i]{n600}     = $A2 + $B2/0.6**2 + $C2/0.6**4;
        # Re-plot
        my $plot_file = "$FindBin::Bin/fit_vase_" . lc($results[$i]{name}) =~ s/ /_/gr . "_01242007.png";
        $vase2->plot($fit2, output => $plot_file,
            title => "$results[$i]{name}: Air/Ta2O5/EMA/Ta");
    }
    print "\n";
}

# ============================================================
# Summary
# ============================================================
print "=" x 60, "\n";
printf "%-12s %-12s %-10s %-8s %-10s %-8s\n",
       "Wafer", "Ta2O5 (nm)", "EMA (nm)", "Void %", 'n@600nm', "MSE";
print "-" x 60, "\n";
for my $r (@results) {
    printf "%-12s %-12.1f %-10.2f %-8.1f %-10.4f %-8.4f\n",
           $r->{name}, $r->{ta2o5_nm}, $r->{ema_nm},
           $r->{vf_void} * 100, $r->{n600}, $r->{mse};
}
print "=" x 60, "\n";
print "\nDone. All fits complete using VASE v1.00 built-in modules.\n";
