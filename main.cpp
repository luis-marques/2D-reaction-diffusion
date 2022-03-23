/* Purpose: This program numerically solves a 2D reaction-diffusion problem
 * (Barkley model) using central finite-differences and Explicit (Forward)
 * Euler. It supports parallel programming via the shared-memory paradigm, which
 * is implemented using the OpenMP API. 
 * 
 * Tested on: Spitfire (& Typhoon), Debian GNU/Linux 11 (bullseye).
 * 
 * Author: CID 01715222.
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
    // C = c^2     ==> Valid C = {1, 4,  9, 16, 25, 36, 49,  64} (Number of cores to be supported)
    // While one could constrain the number of threads to solely those present in the C set above,
    // it was opted instead to simply ensure that the number of threads would not exceed 64. While for testing,
    // only {1, 4, 9, 16, 25, 36, 49, 64} threads were used, technically if running with e.g. 30 threads, the 
    // program would be necessarily running in less than 64 cores. 
    // Note: The 4 commented-out lines bellow implement the case where only a number of threads present in the
    // C set is valid.
    
    // Getting the number of threads, as set in OMP_NUM_THREADS env variable.
    int max_nr_threads = omp_get_max_threads();
    
//    if  ( (max_nr_threads != 1) && (max_nr_threads != 4) && (max_nr_threads != 9) &&
//          (max_nr_threads != 16) && (max_nr_threads != 25) && (max_nr_threads != 36) &&
//          (max_nr_threads != 49) && (max_nr_threads != 64) ) {
    if (max_nr_threads > 64) {
            cout << "Error: Program is being run with " << max_nr_threads << " threads." << endl;
//            cout << "This program should only be run with 1, 4, 9, 16, 25, 36, 49 or 64 threads." << endl;
            cout << "The program must be run with 64 or less threads." << endl;
            cout << "Please change the number of threads by altering the 'OMP_NUM_THREADS' env variable." << endl;
            return 0;
    }
    else {  //
        cout << "Running program with " << max_nr_threads << " thread(s)." << endl;
    }
    
    // Specifying the options we want to make available to the user.
    po::options_description opts("Allowed options");
    opts.add_options()
        ("dt", po::value<double>()->required(),             // dt and T are the only required arguments
                 "Time-step to use. (double)")              // since these are the only parameters not
        ("T",  po::value<int>()->required(),                // present in Table 2 of the handout.
                 "Total integration time. (int)")           // for all other arguments default values are
        ("Nx",  po::value<int>()->default_value(101),       // those specified in Table 2.
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

    // Parsing command-line arguments and generating map containg the options and values.
    po::variables_map vm;
    
    // Handling errors with the passed arguments and provides helpful error messages.
    try {
        po::store(po::parse_command_line(argc, argv, opts), vm);
        
        // If the user asked for help, or no arguments were supplied, provide
        // help message and list the possible arguments.
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
    
    // Extracts the parameter values given to varibles using the appropriate datatype.
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
    
    // Passing the parsed parameter values to the class.
    react_diff.SetParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);
    
    // Setting the initial conditions of the 'u','v' solution fields.
    react_diff.SetInitialConditions();
    
    // Perform time-integration from t=0 to t=T with timestep = dt.
    react_diff.TimeIntegrate();
    
    // Saving result of simulation to file called 'output.txt'.
    react_diff.SaveToFile();

    return 0;
}
