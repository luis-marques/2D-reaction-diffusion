/* Program: This program numerically solves a 2D reaction-diffusion problem
 * using parallel computing. 
 * 
 * Tested on: Typhoon (& Spitfire), Debian GNU/Linux 11 (bullseye).
 * 
 * Author: CID-01715222.
 */

#include <iostream>

#include <boost/program_options.hpp>
#include <omp.h>

#include "ReactionDiffusion.h"

// To improve code readibility
using namespace std;
namespace po = boost::program_options;


int main(int argc, char* argv[]) {
    
    // 1 <= c <= 8 ==> Valid c = {1, 2,  3,  4,  5,  6,  7,   8}
    // C = c^2     ==> Valid C = {1, 4,  9, 16, 25, 36, 49,  64} (number of Cores)
    // T = C * 2   ==> Valid T = {2, 8, 16, 32, 50, 72, 98, 128} (number of threads = number of Cores * 2)
    
    // Getting number of threads (as set in OMP_NUM_THREADS env variable) and
    // checking that satisfies the requirement set in Section 2.4 of the handout.
    int max_nr_threads = omp_get_max_threads();
    if  ( (max_nr_threads != 1) && (max_nr_threads != 4) && (max_nr_threads != 9) &&
          (max_nr_threads != 16) && (max_nr_threads != 25) && (max_nr_threads != 36) &&
          (max_nr_threads != 49) && (max_nr_threads != 64) ) {
          
            cout << "Error: Program is being run with " << max_nr_threads << " threads." << endl;
            cout << "This program should only be run with 1, 4, 9, 16, 25, 36, 49 or 64 threads." << endl;
            cout << "Please change the number of threads by altering the 'OMP_NUM_THREADS' env variable." << endl;
            return 0;
        }
    else {
        cout << "Running program with " << max_nr_threads << " thread(s)." << endl;
    }
    
    // Specify the options we want to make available to the user
    po::options_description opts("Allowed options");
    // dt and T are only required parameters since they do not appear in Table 1.
    opts.add_options()
        ("dt", po::value<double>()->required(),
                 "Time-step to use. (double)")
        ("T",  po::value<int>()->required(),
                 "Total integration time. (int)")
        ("Nx",  po::value<int>()->default_value(101),
                 "Number of grid points in x. (int)")
        ("Ny",  po::value<int>()->default_value(101),
                 "Number of grid points in y. (int)")
        ("a",  po::value<double>()->default_value(0.75, "0.75"),
                 "Value of parameter a. (double)")
        ("b",  po::value<double>()->default_value(0.06, "0.06"),
                 "Value of parameter b. (double)")
        ("mu1",  po::value<double>()->default_value(5.0, "5.0"),
                 "Value of parameter mu1. (double)")
        ("mu2",  po::value<double>()->default_value(0.0, "0.0"),
                 "Value of parameter mu2. (double)")
        ("eps",  po::value<double>()->default_value(50.0, "50.0"),
                 "Value of parameter epsilon. (double)")
        ("help", "Prints the help message.")

        ;

    // Tell Boost to parse the command-line arguments using the list of
    // possible options and generate a map (vm) containing the options and
    // values actually specified by the user.
    po::variables_map vm;
    
    // Handles errors with the passed arguments and provides helpful error messages.
    try {
        po::store(po::parse_command_line(argc, argv, opts), vm);
        
        // If the user asked for help, or no arguments were supplied, provide
        // help message and lists the possible arguments.
        if (vm.count("help") || argc == 1) {
            cout << "Please provide a value for all the required options (optional options are indicated with round brackets)." << endl
                 << opts << endl;
            return 0;
        }
        
        po::notify(vm);

    } catch (exception& e) {
        cout << "Error: " << e.what() << endl;
        cout << opts << endl;
        return 0;
    }
    
    // Extracts values given to parameters using the appropriate datatype.
    const double dt = vm["dt"].as<double>();
    const int T = vm["T"].as<int>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();
    const double a = vm["a"].as<double>();
    const double b = vm["b"].as<double>();
    const double mu1 = vm["mu1"].as<double>();
    const double mu2 = vm["mu2"].as<double>();
    const double eps = vm["eps"].as<double>();
 
 //{
    // Initializing Reaction Diffusion class
    ReactionDiffusion react_diff;
    //ReactionDiffusion* react_diff = new ReactionDiffusion();
    
    // Setting the problem parameters as specific in the command-line
    react_diff.SetParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);
    
    //react_diff->SetParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);
    //react_diff->SetInitialConditions();
    //react_diff->TimeIntegrate();
    //react_diff->Terminate();
    react_diff.SetInitialConditions();
    
    //delete react_diff;
    react_diff.TimeIntegrate();
    

    react_diff.Terminate();
 //}
    return 0;
}
