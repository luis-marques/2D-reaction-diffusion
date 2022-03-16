/* Program:
 * 
 * Tested: 
 * 
 * Author: CID-01715222
 */

#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>

// To reduce typing
using namespace std;
namespace po = boost::program_options;

#include "ReactionDiffusion.h"

int main(int argc, char* argv[]) {
    // Specify the options we want to make available to the user
    po::options_description opts(
        "Sorts a list of random numbers using the insertion sort algorithm.");
    opts.add_options()
        ("dt", po::value<double>(),
                 "Time-step to use.")
        ("T",  po::value<int>(),
                 "Total integration time.")
        ("Nx",  po::value<int>()->default_value(101),
                 "Number of grid points in x.")
        ("Ny",  po::value<int>()->default_value(101),
                 "Number of grid points in y.")
        ("a",  po::value<double>()->default_value(0.75),
                 "Value of parameter a.")
        ("b",  po::value<double>()->default_value(0.06),
                 "Value of parameter b.")
        ("mu1",  po::value<double>()->default_value(5.0),
                 "Value of parameter mu1.")
        ("mu2",  po::value<double>()->default_value(0.0),
                 "Value of parameter mu2.")
        ("eps",  po::value<double>()->default_value(50.0),
                 "Value of parameter epsilon.")
        ("descending", "Indicate the array should be reversed.")
        ("help",       "Print help message.");

    // Tell Boost to parse the command-line arguments using the list of
    // possible options and generate a map (vm) containing the options and
    // values actually specified by the user.
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    // Check if the user gave the "--help" option and print the usage.
    if (vm.count("help")) {
        cout << "Performs an insertion sort algorithm on an array of "
             << "random numbers." << endl;
        cout << opts << endl;
        return 0;
    }

    // Extract the values given to other parameters using the appropriate
    // data type.
    const double dt = vm["dt"].as<double>();
    const int T = vm["T"].as<int>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();
    const double a = vm["a"].as<double>();
    const double b = vm["b"].as<double>();
    const double mu1 = vm["mu1"].as<double>();
    const double mu2 = vm["mu2"].as<double>();
    const double eps = vm["eps"].as<double>();
    
    // Initializing Reaction Diffusion class
    ReactionDiffusion react_diff;
    
    // Setting the problem parameters as specific in the command-line
    react_diff.SetParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);
    
    react_diff.SetInitialConditions();
    
    react_diff.TimeIntegrate();
    
    react_diff.Terminate();
    
    return 0;
}
