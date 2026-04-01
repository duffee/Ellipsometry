use Physics::Ellipsometry::VASE;
use PDL;
use PDL::NiceSlice;

my $vase = Physics::Ellipsometry::VASE->new(Layers => 1);
my $data = $vase->load_data('../data/Metal_Oxides/tantalum oxide/Cap_11012006/w1_11012006.dat');

print "Points: ", $data->getdim(1), "\n";
print "Sample: ", $vase->{sample_name}, "\n";
print "Has sigma: ", defined $vase->{sigma} ? "yes" : "no", "\n";

# Access sigma piddle: column 0 = sigma_psi, column 1 = sigma_delta
my $sigma = $vase->{sigma};

# Extract individual columns (using NiceSlice syntax in your model)
my $sigma_psi = $sigma((0),:)->flat;
my $sigma_delta = $sigma((1),:)->flat;

# The fit() method already uses these automatically as weights for the Levenberg-Marquardt fit when they're present - points with larger uncertainties get less weight. 
