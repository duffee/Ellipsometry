use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use FindBin;
use Physics::Ellipsometry::VASE;

# Create VASE object with 1 layer
my $vase = Physics::Ellipsometry::VASE->new(layers => 1);

# Load sample data
$vase->load_data("$FindBin::Bin/sample.dat");

# Define model function (linear model example)
sub model {
    my ($params, $x) = @_;
    
    # Unpack parameters
    my $a = $params->(0);
    my $b = $params->(1);
    my $c = $params->(2);
    my $d = $params->(3);
    
    # Compute linear model (using only wavelength)
    my $wavelength = $x->(:,0);   # first column: wavelength

    my $psi = $a - $b * $wavelength;
    my $delta = $c + $d * $wavelength;
    
    return cat($psi, $delta)->flat;
}

$vase->set_model(\&model);

# Initial parameters: [a, b, c, d] for linear model
my $initial_params = pdl [65, 0.05, 80, 0.1];

# Perform fit
my $fit_params = $vase->fit($initial_params);

# Extract results
my ($a, $b, $c, $d) = list $fit_params;
print "Fit results:\n";
printf "  a = %.6f\n", $a;
printf "  b = %.6f\n", $b;
printf "  c = %.6f\n", $c;
printf "  d = %.6f\n", $d;
printf "  MSE = %.6f\n", $vase->mse($fit_params, nparams => 4);
printf "  Iterations: %d\n", $vase->{iters};
