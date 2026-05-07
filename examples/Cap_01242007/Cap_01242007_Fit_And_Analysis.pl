#!perl
# Data fitting and analysis for Cap_01242007 data set
# Jovan Trujillo
# Advanced Electronics and Photonics Core
# Arizona State University
# 4/27/2026

use PDL;
use PDL::NiceSlice;
use Physics::Ellipsometry::VASE;


my $vase = Physics::Ellipsometry::VASE->new(layers => 1);
my $data = $vase->load_data('wafer_1_01242007.dat');

print "Sample: ", $vase->{sample_name}, "\n";
print "Method: ", $vase->{vase_method}, "\n";
print "Original: ", $vase->{original_file}, "\n";
print "Units: ", $vase->{units}, "\n";

# Access measurement uncertainties (automatically used as fit weights)
if (defined $vase->{sigma}) {
	my $sigma_psi = $vase->{sigma}->(0,:); # psi uncertainties
	my $sigma_delta = $vase->{sigma}->(1,:); #delta uncertainties
	print "Has sigma columns for weighted fitting\n";
}
