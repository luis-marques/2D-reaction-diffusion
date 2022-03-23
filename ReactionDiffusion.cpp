#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

using namespace std;

/**
 * @brief Takes in the command-line arguments as parsed by boost_program_options and stores
 * their values in the class members. Allocates memory to solutions fields. Pre-calculates any
 * variables that will be repeatedly used in time integration.
 * @param arg_dt Integration time-step.
 * @param arg_T Total integration time.
 * @param arg_Nx Number of grid-points/nodes along x direction.   
 * @param arg_Ny Number of grid-points/nodes along y direction
 * @param arg_a Parameter 'a' of Barkley model problem.
 * @param arg_b Parameter 'b' of Barkley model problem.
 * @param arg_mu1 Diffusion coefficient for 'u'.
 * @param arg_mu2 Diffusion coefficient for 'v'.
 * @param arg_eps Parameter 'eps' of Barkley model problem.
 */
void ReactionDiffusion::SetParameters(
                    const double& arg_dt, const int& arg_T,
                    const int& arg_Nx, const int& arg_Ny,
                    const double& arg_a, const double& arg_b,
                    const double& arg_mu1, const double& arg_mu2,
                    const double& arg_eps) {
                        
                        
        // Storing arguments in class variables.
        dt  = arg_dt;
        T   = arg_T;
        Nx  = arg_Nx;
        Ny  = arg_Ny;
        a   = arg_a;
        b   = arg_b;
        mu1 = arg_mu1;
        mu2 = arg_mu2;
        eps = arg_eps;
        
        // Pre-computing helpful variables, as specified in header file.
        recip_a     = 1.0 / a; 
        u_grad_coef = mu1 * dt;
        v_grad_coef = mu2 * dt;
        dt_eps      = dt * eps;
        b_over_a    = b / a;
        half_a      = a / 2.0;
        Lx_index    = (Nx-1);
        Ly_index    = Nx * (Ny-1);
        
        // Allocating memomry for the solution fields.
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        u_next = new double[Nx*Ny];
        v_next = new double[Nx*Ny];
        

        // Add check in case T is not properly divisible by dt (check remainder)!!!!!!!
        nr_timesteps = T / dt;
        cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << endl; 
        
        // Printing values of all the parameters to the terminal.
        cout << endl;
        cout << "Parameters of PDE problem to solve:"        << endl;
        cout << left << setw(18) << "    - dt"        << dt  << endl;
        cout << left << setw(18) << "    - T:"        << T   << endl;
        cout << left << setw(18) << "    - Nx:"       << Nx  << endl;
        cout << left << setw(18) << "    - Ny:"       << Ny  << endl;
        cout << left << setw(18) << "    - a:"        << a   << endl;
        cout << left << setw(18) << "    - b:"        << b   << endl;
        cout << left << setw(18) << "    - mu1:"      << mu1 << endl;
        cout << left << setw(18) << "    - mu2:"      << mu2 << endl;
        cout << left << setw(18) << "    - epsilon:"  << eps << endl;
        cout << endl;
        
}

/**
 * @brief Populates the u,v arrays with their initial conditions. These are stores in
 * column-major fomat.
 */
void ReactionDiffusion::SetInitialConditions() {
    
        
    // when u allocate on heap others are 0 automatically
    
    // Data 
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            // If y < Ly/2 (note that Ly = (Ny-1)*dy, where we are taking dy=1.
            if (j > (Ny-1)/2.0) {
                u[i + j*Nx] = 1.0;
            }
            else{
                u[i + j*Nx] = 0.0;
            } 

            // If x < Lx/2 (note that Lx = (Nx-1)*dx, where we are taking dx=1.
            if (i < (Nx-1)/2.0) {
                v[i + j*Nx] = half_a;
            }
            else{
                v[i + j*Nx] = 0.0;
            }        
        
        }
    }
}


/**
 * @brief Performs integration over time using the Explicit (Forward)
 * Euler numerical scheme.
 * 
 * A discrete update of the u,v arrays is run \f$ n=T/dt \f$ times
 * where \f$ dt \f$ is timestep used for integration and
 * \f$ T \f$ the total integration time.
 */
void ReactionDiffusion::TimeIntegrate() {
    cout << "Starting numerical solving of PDE.\n" << endl;
    
    // Each iteration of this for-loop performs one update of the 'u','v' solution fields.
    // i.e. 
    for (int timestep = 0; timestep < nr_timesteps; ++timestep) {
        
        // Single omp parallel, so there is only one forking/joining of threads per iteration.
        #pragma omp parallel
        { 
            
            // We are calculating first the update to the 'u' solution field and only later the update to the 'v' solution field
            // to try to exploit cache locality.
            
            // Corner (0, 0) (i=0, j=0)
            #pragma omp single nowait
            u_next[0] = u[0] + u_grad_coef * (u[1] + u[Nx] - 2*u[0]) 
                      + dt_eps * u[0] * (1.0 - u[0]) * (u[0] - v[0] * recip_a - b_over_a);
            
            // Corner (Lx, 0) (i=Nx-1, j=0) Recall: Lx_index = (Nx-1)
            #pragma omp single nowait
            u_next[Lx_index] = u[Lx_index] + u_grad_coef*(u[Lx_index - 1] + u[Lx_index + Nx] - 2*u[Lx_index])
                             + dt_eps * u[Lx_index] * (1.0 - u[Lx_index]) * (u[Lx_index] - v[Lx_index] * recip_a - b_over_a);
            
            // Corner (0, Ly) (i=0, j=Ny-1) Recall: Ly_index = Nx*(Ny-1)
            #pragma omp single nowait
            u_next[Ly_index] = u[Ly_index] + u_grad_coef*(u[Ly_index + 1] + u[Ly_index - Nx] - 2*u[Ly_index])
                             + dt_eps * u[Ly_index] * (1.0 - u[Ly_index]) * (u[Ly_index] - v[Ly_index] * recip_a - b_over_a);
           
            // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
            #pragma omp single nowait
            u_next[Lx_index + Ly_index] = u[Lx_index + Ly_index] + u_grad_coef*(u[Lx_index + Ly_index - 1] + u[Lx_index + Ly_index - Nx] - 2*u[Lx_index + Ly_index]) 
                                        + dt_eps * u[Lx_index+Ly_index] * (1.0 - u[Lx_index+Ly_index]) * (u[Lx_index+Ly_index] - v[Lx_index+Ly_index] * recip_a - b_over_a);
            
            // Edge BCs for fixed v values.
            #pragma omp for nowait
            for (int i = 1; i < (Nx-1); i++) {
                // Along (y==0) (0<i<Nx, j=0)
                u_next[i] = u[i] + u_grad_coef*(u[i+1] + u[i-1] + u[i + Nx] - 3*u[i]) 
                          + dt_eps * u[i] * (1.0 - u[i]) * (u[i] - v[i] * recip_a - b_over_a);
                          
                          
                // Along (y==Ly) (0<i<Nx, j=Ny-1)
                u_next[i + Ly_index] = u[i + Ly_index] + u_grad_coef*(u[i+1 + Ly_index] + u[i-1 + Ly_index] + u[i - Nx + Ly_index] - 3*u[i + Ly_index])
                                 + dt_eps * u[i+Ly_index] * (1.0 - u[i+Ly_index]) * (u[i+Ly_index] - v[i+Ly_index] * recip_a - b_over_a);
            }
            

            // Central points (0<i<Nx-1, 0<j<Ny-1)
            #pragma omp for nowait 
            for (int j = 1; j < (Ny-1); ++j) {
                for (int i = 1; i < (Nx-1); ++i) {
                    u_next[i+j*Nx] = u[i+j*Nx] + u_grad_coef*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 4*u[i+j*Nx])
                                   + dt_eps * u[i+j*Nx] * (1.0 - u[i+j*Nx]) * (u[i+j*Nx] - v[i+j*Nx] * recip_a - b_over_a);
                    
                }
            }
            
            // Edge BCs for fixed x values.
            #pragma omp for nowait
            for (int j = 1; j < (Ny-1); j++) {
                
                // Along (x==0) (i=0, 1<j<Ny-1)
                u_next[j*Nx] = u[j*Nx] + u_grad_coef*(u[1 + j*Nx] + u[(j+1)*Nx] + u[(j-1)*Nx] - 3*u[j*Nx])
                             + dt_eps * u[j*Nx] * (1.0 - u[j*Nx]) * (u[j*Nx] - v[j*Nx] * recip_a - b_over_a);
                
                // Along (x==Lx) (i=Nx-1, 1<j<Ny-1)
                u_next[Lx_index + j*Nx] = u[Lx_index + j*Nx] + u_grad_coef*(u[Lx_index - 1 + j*Nx] + u[Lx_index + (j+1)*Nx] + u[Lx_index + (j-1)*Nx] - 3*u[Lx_index + j*Nx])
                                        + dt_eps * u[Lx_index + j*Nx] * (1.0 - u[Lx_index + j*Nx]) * (u[Lx_index + j*Nx] - v[Lx_index + j*Nx] * recip_a - b_over_a);

            
            }
            
            
            // Now calculating the updates for the 'v' field.
            
            
             // Corner (0, 0) (i=0, j=0)
            #pragma omp single nowait
            v_next[0] = v[0] + v_grad_coef * (v[1] + v[Nx] - 2*v[0])
                      + dt * (u[0] * u[0] * u[0] - v[0]);
                      
            // Corner (Lx, 0) (i=Nx-1, j=0) Lx_index = (Nx-1)
            #pragma omp single nowait
            v_next[Lx_index] = v[Lx_index] + v_grad_coef*(v[Lx_index - 1] + v[Lx_index + Nx] - 2*v[Lx_index])
                             + dt * (u[Lx_index] * u[Lx_index] * u[Lx_index] - v[Lx_index]);
                             
            // Corner (0, Ly) (i=0, j=Ny-1) Ly_index = Nx*(Ny-1)
            #pragma omp single nowait
            v_next[Ly_index] = v[Ly_index] + v_grad_coef*(u[Ly_index + 1] + v[Ly_index - Nx] - 2*v[Ly_index])
                            + dt * (u[Ly_index] * u[Ly_index] * u[Ly_index] - v[Ly_index]);
                            
            // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
            #pragma omp single nowait
            v_next[Lx_index + Ly_index] = v[Lx_index + Ly_index] + v_grad_coef*(v[Lx_index + Ly_index - 1] + v[Lx_index + Ly_index - Nx] - 2*v[Lx_index + Ly_index])
                                        + dt * (u[Lx_index+Ly_index] * u[Lx_index+Ly_index] * u[Lx_index+Ly_index] - v[Lx_index+Ly_index]);
            
            // Edge BCs for fixed y values.
            #pragma omp for nowait
            for (int i = 1; i < (Nx-1); i++) {
                
                // Along (y==0) (0<i<Nx, j=0)
                v_next[i] = v[i] + v_grad_coef*(v[i+1] + v[i-1] + v[i + Nx] - 3*v[i])
                          + dt * (u[i] * u[i] * u[i] - v[i]);
                          
                // Along (y==Ly) (0<i<Nx, j=Ny-1)
                v_next[i + Ly_index] = v[i+Ly_index] + v_grad_coef*(v[i+1 + Ly_index] + v[i-1 + Ly_index] + v[i - Nx + Ly_index] - 3*v[i + Ly_index])
                                     + dt * (u[i+Ly_index] * u[i+Ly_index] * u[i+Ly_index] - v[i+Ly_index]);

            }
            
            // Central points (0<i<Nx-1, 0<j<Ny-1)
            #pragma omp for nowait 
            for (int j = 1; j < (Ny-1); ++j) {
                for (int i = 1; i < (Nx-1); ++i) {     
                    v_next[i+j*Nx] = v[i+j*Nx] + v_grad_coef*(v[i+1 + j*Nx] + v[i-1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 4*v[i+j*Nx])
                                   + dt * (u[i+j*Nx] * u[i+j*Nx] * u[i+j*Nx] - v[i+j*Nx]);
                }
            }

        
            // Edge BCs for fixed x values.
            #pragma omp for nowait
            for (int j = 1; j < (Ny-1); j++) {
                
                // Along (x==Lx) (i=Nx-1, 0<j<Ny-1)
                v_next[Lx_index + j*Nx] = v[Lx_index + j*Nx] + v_grad_coef*(v[Lx_index - 1 + j*Nx] + v[Lx_index + (j+1)*Nx] + v[Lx_index + (j-1)*Nx] - 3*v[Lx_index + j*Nx])
                                        + dt * (u[Lx_index+j*Nx] * u[Lx_index+j*Nx] * u[Lx_index+j*Nx] - v[Lx_index+j*Nx]);
            
                // Along (x==0) (i=0, 0<j<Ny-1)
                v_next[j*Nx] = v[j*Nx] + v_grad_coef*(v[1 + j*Nx] + v[(j+1)*Nx] + v[(j-1)*Nx] - 3*v[j*Nx])
                             + dt * (u[j*Nx] * u[j*Nx] * u[j*Nx] - v[j*Nx]);
            
            }
        
        }   // End of pragma parallel region (there is an implicit barrier)
        
        // All the threads have finished running due to implicit barrier above, thus can now store the
        // calculated u^{n+1} = u_next and v^{n+1} = v_next into 'u','v' and prepare for following time-step.
        swap(u, u_next);
        swap(v, v_next);
        
        
        // Providing feedback to the user of how fast the program is running during execution.
        // Has negligible performance impact.
        if (timestep % 10000 == 0) { 
            cout << "Time-step: " << right << setw(6) << timestep 
                 << "/" << nr_timesteps 
                 << " (" << right << setw(2) << 100*timestep/nr_timesteps << "%)" << endl;
        }
              
    }   // End of time-loop (i.e. finished integrating over time).
    
    cout << "\nFinished solving PDE (from t_i = 0 to t_f = T).\n" << endl;
}


/**
 * @brief Responsible for saving results of numerical simulation to a text
 * file called 'output.txt'. Also deals with memory de-allocation and
 * other procedures that should be run at the end of the program execution.
 * 
 * To visualize the result of the simulation one can run the following commands:
 * $ gnuplot
 * $ set pm3d at st
 * $ set view map
 * $ set cbrange[0:1]
 * 
 * and then
 * $ splot 'output.txt' using 1:2:3 w l palette
 * to plot 'u' over the domain, or instead
 * $ splot 'output.txt' using 1:2:4 w l palette
 * to plot 'v' over the domain.
 * 
 * xrange[0:150] for last plot
 * 
 */
void ReactionDiffusion::SaveToFile() {
    
    cout << "Writting output of simulation to file 'output.txt'." << endl;
    
    ofstream vOut("output.txt", ios::out | ios::trunc);
    
    // Checking that file opened successfuly.
    if (vOut.is_open()) {
        cout << "File opened successfully." << endl;
        
        // Writing solution row-by-row (x y u v)
        // Not super efficient since vectors are stored column-wise. 
        // One could instead store the vectors row-wise to improve performance here but
        // since this is only called once, it was deemed that the potential performance
        // increase (through cache-locality) was not worth it to change from the convention
        // of storing matrices column-wise.
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i){
                vOut << i << " " 
                     << j << " "
                     << u[i + j*Nx] << " "
                     << v[i + j*Nx] << endl;
            }
            vOut << endl;   // Empty line after each row of points
        }
    }
    else {
        cout << "Did not open vOut successfully!" << endl;
    }
    vOut.close();   // Closing file.
    cout << "Finished writing to file." << endl;
        
}

/**
 * @brief Destructor de-allocates all dynamically assigned memory.
 */
ReactionDiffusion::~ReactionDiffusion() {
        
    // De-allocating memory.
    delete[] u;
    delete[] v;
    delete[] u_next;
    delete[] v_next;
}