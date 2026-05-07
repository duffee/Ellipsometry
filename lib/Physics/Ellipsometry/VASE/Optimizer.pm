package Physics::Ellipsometry::VASE::Optimizer;
use strict;
use warnings;
use PDL;
use Exporter 'import';

our @EXPORT_OK = qw(differential_evolution grid_search);

our $VERSION = '0.01';

=head1 NAME

Physics::Ellipsometry::VASE::Optimizer - Global optimization algorithms

=head1 DESCRIPTION

Provides global optimization methods for finding initial parameter estimates
before Levenberg-Marquardt refinement. Currently implements:

=over

=item Differential Evolution (DE/rand/1/bin)

Population-based stochastic optimizer. Robust against local minima.
Parameters: population size, mutation factor (F), crossover probability (CR).

=item Grid Search

Systematic search over a parameter grid. Good for 1-2 parameters with
known ranges; impractical for high dimensions.

=back

=cut

# Differential Evolution (DE/rand/1/bin)
# $objective->($params_pdl) returns scalar cost (e.g., chi²)
# $bounds: arrayref of [min, max] pairs for each parameter
# Returns: best parameter PDL
sub differential_evolution {
    my (%args) = @_;
    die "Need objective function" unless defined $args{objective};
    die "Need parameter bounds"   unless defined $args{bounds};
    my $objective = $args{objective};
    my $bounds    = $args{bounds};
    my $np        = $args{pop_size}  // 30;
    my $F         = $args{F}         // 0.7;
    my $CR        = $args{CR}        // 0.9;
    my $maxiter   = $args{maxiter}   // 200;
    my $tol       = $args{tol}       // 1e-6;
    my $seed      = $args{seed};
    my $verbose   = $args{verbose}   // 0;

    srand($seed) if defined $seed;

    my $ndim = scalar @$bounds;
    $np = $ndim * 10 if $np < $ndim * 5;  # ensure adequate population

    # Initialize population randomly within bounds
    my @population;
    my @costs;
    for my $i (0 .. $np - 1) {
        my @individual;
        for my $d (0 .. $ndim - 1) {
            my ($lo, $hi) = @{$bounds->[$d]};
            push @individual, $lo + rand() * ($hi - $lo);
        }
        my $ind_pdl = pdl(\@individual);
        push @population, $ind_pdl;
        push @costs, $objective->($ind_pdl);
    }

    # Find initial best
    my $best_idx = 0;
    for my $i (1 .. $#costs) {
        $best_idx = $i if $costs[$i] < $costs[$best_idx];
    }
    my $best_cost = $costs[$best_idx];
    my $best = $population[$best_idx]->copy;

    printf "  DE: initial best cost = %.4f\n", $best_cost if $verbose;

    # Evolution loop
    for my $gen (1 .. $maxiter) {
        my $improved = 0;

        for my $i (0 .. $np - 1) {
            # Select 3 distinct random indices ≠ i
            my @r;
            while (@r < 3) {
                my $idx = int(rand($np));
                next if $idx == $i || grep { $_ == $idx } @r;
                push @r, $idx;
            }

            # Mutation: v = x_r0 + F*(x_r1 - x_r2)
            my $v = $population[$r[0]] + $F * ($population[$r[1]] - $population[$r[2]]);

            # Clip to bounds
            for my $d (0 .. $ndim - 1) {
                my ($lo, $hi) = @{$bounds->[$d]};
                my $val = $v->at($d);
                $val = $lo if $val < $lo;
                $val = $hi if $val > $hi;
                $v->set($d, $val);
            }

            # Crossover: binomial
            my $trial = $population[$i]->copy;
            my $j_rand = int(rand($ndim));
            for my $d (0 .. $ndim - 1) {
                if (rand() < $CR || $d == $j_rand) {
                    $trial->set($d, $v->at($d));
                }
            }

            # Selection
            my $trial_cost = $objective->($trial);
            if ($trial_cost < $costs[$i]) {
                $population[$i] = $trial;
                $costs[$i] = $trial_cost;
                $improved++;

                if ($trial_cost < $best_cost) {
                    $best_cost = $trial_cost;
                    $best = $trial->copy;
                }
            }
        }

        if ($verbose && $gen % 20 == 0) {
            printf "  DE gen %d: best=%.4f, improved=%d/%d\n",
                   $gen, $best_cost, $improved, $np;
        }

        # Convergence check: population diversity
        if ($gen > 10) {
            my $spread = 0;
            for my $d (0 .. $ndim - 1) {
                my @vals = map { $_->at($d) } @population;
                my $min_v = (sort { $a <=> $b } @vals)[0];
                my $max_v = (sort { $b <=> $a } @vals)[0];
                my ($lo, $hi) = @{$bounds->[$d]};
                $spread += ($max_v - $min_v) / (($hi - $lo) || 1);
            }
            $spread /= $ndim;
            last if $spread < $tol;
        }
    }

    printf "  DE: final best cost = %.4f\n", $best_cost if $verbose;
    return ($best, $best_cost);
}

# Grid search over specified parameter dimensions
# $objective->($params_pdl) returns scalar cost
# $base_params: PDL with default parameter values
# $grid_spec: arrayref of {index => param_idx, min => val, max => val, steps => N}
sub grid_search {
    my (%args) = @_;
    die "Need objective function" unless defined $args{objective};
    die "Need base_params PDL"   unless defined $args{base_params};
    die "Need grid specification" unless defined $args{grid};
    my $objective   = $args{objective};
    my $base_params = $args{base_params};
    my $grid_spec   = $args{grid};
    my $verbose     = $args{verbose} // 0;

    my $best_cost   = 1e30;
    my $best_params = $base_params->copy;

    # For 1D or 2D grid search
    my @axes;
    for my $spec (@$grid_spec) {
        my $step = ($spec->{max} - $spec->{min}) / ($spec->{steps} - 1);
        my @values;
        for my $i (0 .. $spec->{steps} - 1) {
            push @values, $spec->{min} + $i * $step;
        }
        push @axes, { index => $spec->{index}, values => \@values };
    }

    # Recursive grid evaluation
    _grid_recurse(\@axes, 0, $base_params->copy, $objective,
                  \$best_cost, \$best_params);

    printf "  Grid: best cost = %.4f\n", $best_cost if $verbose;
    return ($best_params, $best_cost);
}

sub _grid_recurse {
    my ($axes, $depth, $params, $objective, $best_cost_ref, $best_params_ref) = @_;

    if ($depth >= scalar @$axes) {
        my $cost = $objective->($params);
        if ($cost < $$best_cost_ref) {
            $$best_cost_ref = $cost;
            $$best_params_ref = $params->copy;
        }
        return;
    }

    my $axis = $axes->[$depth];
    for my $val (@{$axis->{values}}) {
        my $p = $params->copy;
        $p->set($axis->{index}, $val);
        _grid_recurse($axes, $depth + 1, $p, $objective, $best_cost_ref, $best_params_ref);
    }
}

1;
