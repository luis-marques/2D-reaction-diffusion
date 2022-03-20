#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>

#include <omp.h>
using namespace std;

/**
 * @brief Saves
 * @param arg_dt
 * @param arg_T
 * @param arg_Nx
 * @param arg_Ny
 * @param arg_a
 * @param arg_b
 * @param arg_mu1
 * @param arg_mu2
 * @param arg_eps
 */
void ReactionDiffusion::SetParameters(
                    const double& arg_dt, const int& arg_T,
                    const int& arg_Nx, const int& arg_Ny,
                    const double& arg_a, const double& arg_b,
                    const double& arg_mu1, const double& arg_mu2,
                    const double& arg_eps) {
                        
        dt = arg_dt;
        T = arg_T;
        Nx = arg_Nx;
        Ny = arg_Ny;
        a = arg_a;
        b = arg_b;
        mu1 = arg_mu1;
        mu2 = arg_mu2;
        eps = arg_eps;
        
        // 
        recip_a = 1.0 / a; // reciprocal of "a"
        u_grad_coef = mu1 * dt;
        v_grad_coef = mu2 * dt;
        dt_eps = dt * eps;
        b_over_a = b / a;
        half_a = a / 2.0;
        
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        u_next = new double[Nx*Ny];
        v_next = new double[Nx*Ny];
        
        Ly_index = Nx*(Ny-1);
        Lx_index = (Nx-1);

        // Add check in case T is not properly divisible by dt (check remainder)
        nr_timesteps = T / dt;
        cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << endl; 
        
        // For Debugging
        cout << "Parameters of PDE problem to solve:" << endl;
        cout << "\t*Time-step (dt)" << right << setw(30) << setfill(' ') << dt << endl;
        cout << "\t*Integration time (T)" << right << setw(30) << setfill(' ') << T << endl;
        cout << "\t*Nx" << right << setw(30) << setfill(' ') << Nx << endl;
        cout << "\t*Ny" << right << setw(30) << setfill(' ') << Ny << endl;
        cout << "\t*a" << right << setw(30) << setfill(' ') << a << endl;
        cout << "\t*b" << right << setw(30) << setfill(' ') << b << endl;
        cout << "\t*mu1" << right << setw(30) << setfill(' ') << mu1 << endl;
        cout << "\t*mu2" << right << setw(30) << setfill(' ') << mu2 << endl;
        cout << "\t*epsilon" << right << setw(30) << setfill(' ') << eps << endl;
        
}

/**
 * @brief Populates the u,v arrays with their initial conditions.
 */
void ReactionDiffusion::SetInitialConditions() {
    
    //int Lx = (Nx - 1);
    //int Ly = (Ny - 1);
    
//    cout << "Lx=" << Lx << "; Ly=" << Ly << endl;
//    cout << "Ly/2 = " << Lx/2.0 << "; Lx/2 = " << Ly/2.0 << endl;
//    cout << "floor(Ny/2) = " << floor(Ny/2.0) << "; Ny = " << Ny << endl;
//    cout << "ceil(Nx/2) = " << ceil(Nx/2.0) << "; Nx = " << Nx << endl;
        
    // Storing vectors in columnwise (column-major format)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            if (j > (Ny-1)/2.0) {
                u[i + j*Nx] = 1.0;
            }
            else{
                u[i + j*Nx] = 0.0;
            }            
        }
    }

     for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            if (i < (Nx-1)/2.0) {
                v[i + j*Nx] = half_a; // make new variable = a/2 to pre-compute this!
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
    cout << "Starting numerical solving of PDE." << endl;
    
    
    //

   
    for (int timestep = 0; timestep < nr_timesteps; ++timestep) {
        
        
        #pragma omp parallel
        { 

                
                
                
                // U !!!!!!!!!

                    // Along same column j=0 (3 below)
                    // Corner (0, 0) (i=0, j=0)
            #pragma omp sections
            {
                #pragma omp section
                u_next[0] = u[0] + u_grad_coef * (u[1] + u[Nx] - 2*u[0]);
                
               
                // Along (y==0) (0<i<Nx, j=0)
                #pragma omp section
                for (int i = 1; i < (Nx-1); i++) {
                   // #pragma omp task firstprivate(i)
                   // {
                        u_next[i] = u[i] + u_grad_coef*(u[i+1] + u[i-1] + u[i + Nx] - 3*u[i]);
                   // }
                }
            
                
               // #pragma omp task
                //{
                    // Corner (Lx, 0) (i=Nx-1, j=0) Lx_index = (Nx-1)
                #pragma omp section
                u_next[Lx_index] = u[Lx_index] + u_grad_coef*(u[Lx_index - 1] + u[Lx_index + Nx] - 2*u[Lx_index]);
                //}
                
                
                
                // Central points (0<i<Nx-1, 0<j<Ny-1)
                #pragma omp section
                for (int j = 1; j < (Ny-1); ++j) {
                    for (int i = 1; i < (Nx-1); ++i) {
                        u_next[i+j*Nx] = u[i+j*Nx] + u_grad_coef*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 4*u[i+j*Nx]);
                    }
                }
                
                
                // Along (x==0) (i=0, 0<j<Ny-1)
                #pragma omp section
                for (int j = Nx; j < Nx*(Ny-1); j+=Nx) {
                    u_next[j] = u[j] + u_grad_coef*(u[1 + j] + u[j+Nx] + u[j-Nx] - 3*u[j]);
                }
                
            
                // Along (x==Lx) (i=Nx-1, 0<j<Ny-1)
                #pragma omp section
                for (int j = 1; j < (Ny-1); j++) {
                    u_next[Lx_index + j*Nx] = u[Lx_index + j*Nx] + u_grad_coef*(u[Lx_index - 1 + j*Nx] + u[Lx_index + (j+1)*Nx] + u[Lx_index + (j-1)*Nx] - 3*u[Lx_index + j*Nx]);
                }
                
                
                
                    // Along same column j=Ny-1 (3 below)
                    // Corner (0, Ly) (i=0, j=Ny-1) Ly_index = Nx*(Ny-1)
                #pragma omp section
                u_next[Ly_index] = u[Ly_index] + u_grad_coef*(u[Ly_index + 1] + u[Ly_index - Nx] - 2*u[Ly_index]);
                 
 
                // Along (y==Ly) (0<i<Nx, j=Ny-1)
                #pragma omp section
                for (int i = 1; i < (Nx-1); i++) {
                    u_next[i + Ly_index] = u[i + Ly_index] + u_grad_coef*(u[i+1 + Ly_index] + u[i-1 + Ly_index] + u[i - Nx + Ly_index] - 3*u[i + Ly_index]);
                }
              
                    // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
                #pragma omp section
                u_next[Lx_index + Ly_index] = u[Lx_index + Ly_index] + u_grad_coef*(u[Lx_index + Ly_index - 1] + u[Lx_index + Ly_index - Nx] - 2*u[Lx_index + Ly_index]);

                // V !!!!!!!!!!!!
                    // Along same column j=0 (3 below)
                    // Corner (0, 0) (i=0, j=0)
                #pragma omp section
                v_next[0] = v[0] + v_grad_coef * (v[1] + v[Nx] - 2*v[0]);
                 

                // Along (y==0) (0<i<Nx, j=0)
                #pragma omp section
                for (int i = 1; i < (Nx-1); i++) {
                    v_next[i] = v[i] + v_grad_coef*(v[i+1] + v[i-1] + v[i + Nx] - 3*v[i]);
                }
                
                  
                    // Corner (Lx, 0) (i=Nx-1, j=0) Lx_index = (Nx-1)
                #pragma omp section
                v_next[Lx_index] = v[Lx_index] + v_grad_coef*(v[Lx_index - 1] + v[Lx_index + Nx] - 2*v[Lx_index]);
                
          
                // Along (x==Lx) (i=Nx-1, 0<j<Ny-1)
                #pragma omp section
                for (int j = 1; j < (Ny-1); j++) {
                    v_next[Lx_index + j*Nx] = v[Lx_index + j*Nx] + v_grad_coef*(v[Lx_index - 1 + j*Nx] + v[Lx_index + (j+1)*Nx] + v[Lx_index + (j-1)*Nx] - 3*v[Lx_index + j*Nx]);
                }
               
                 
                 
                // Central points (0<i<Nx-1, 0<j<Ny-1)
                #pragma omp section
                for (int j = 1; j < (Ny-1); ++j) {
                    for (int i = 1; i < (Nx-1); ++i) {
                        v_next[i+j*Nx] = v[i+j*Nx] + v_grad_coef*(v[i+1 + j*Nx] + v[i-1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 4*v[i+j*Nx]);
                    }
                }


                // Along (x==0) (i=0, 0<j<Ny-1)
                #pragma omp section
                for (int j = 1; j < (Ny-1); j++) {
                    v_next[j*Nx] = v[j*Nx] + v_grad_coef*(v[1 + j*Nx] + v[(j+1)*Nx] + v[(j-1)*Nx] - 3*v[j*Nx]);
                }
                
                // Along same column j=Ny-1 (3 below)
                // Corner (0, Ly) (i=0, j=Ny-1) Ly_index = Nx*(Ny-1)
                #pragma omp section
                v_next[Ly_index] = v[Ly_index] + v_grad_coef*(u[Ly_index + 1] + v[Ly_index - Nx] - 2*v[Ly_index]);
                 
       
                // Along (y==Ly) (0<i<Nx, j=Ny-1)
                #pragma omp section
                for (int i = 1; i < (Nx-1); i++) {
                    v_next[i + Ly_index] = v[i+Ly_index] + v_grad_coef*(v[i+1 + Ly_index] + v[i-1 + Ly_index] + v[i - Nx + Ly_index] - 3*v[i + Ly_index]);
                }
            
                
                // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
                #pragma omp section
                v_next[Lx_index + Ly_index] = v[Lx_index + Ly_index] + v_grad_coef*(v[Lx_index + Ly_index - 1] + v[Lx_index + Ly_index - Nx] - 2*v[Lx_index + Ly_index]);
                
                        
                
                //#pragma omp parallel
        //        //#pragma omp parallel for collapse(2)
        //        for (int j = 0; j < Ny; ++j) {
        //            //#pragma omp for nowait schedule(static)
        //            for (int i = 0; i < Nx; ++i) {

//                #pragma omp section
//                if (omp_get_thread_num() == 0) {
//                    if (timestep % 10000 == 0) { 
//                        cout << "Timestep " << timestep << endl;
//                    }
//                }
                
                // dont copy entire array, do pointer arithmetic instead!!!!!!!!!!!
//                    temp = u;
//                    u = u_1;
//                    u_1 = temp;

                    #pragma omp wait
                    
                    #pragma omp section
                    for (int k = 0; k < Nx*Ny; ++k) {
                        //#pragma omp task firstprivate(k)
                       // {
                            
                            //u_squared = u_prev[k] * u_prev[k];
                            
                            u_next[k] = u_next[k] + dt_eps * u[k] * (1.0 - u[k]) * (u[k] - v[k] * recip_a - b_over_a);
                            v_next[k] = v_next[k] + dt * (u[k] * u[k] * u[k] - v[k]);
                        
                            u[k] = u_next[k];
                            v[k] = v_next[k];
                      //  }

                    }
                
            } // End of sections (there's an implicit barrier)
            
        //    } // End of single
//        
//        }    // End of parallel
        }   // End of parallel
    } // End time-integrate
    // delete[] dummy;
    cout << "Finished solving PDE (from t_i = 0 to t_f = T)." << endl;
    
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
 */
void ReactionDiffusion::Terminate() {
    
    cout << "Writting output of simulation to file 'output.txt'." << endl;
    
    ofstream vOut("output.txt", ios::out | ios::trunc);
    
    if (vOut.is_open()) {
        cout << "File opened successfully" << endl;
        
        // Writing solution row-by-row (x y u v)
        // Not super efficient since vectors are stored column-wise. 
        // One could instead store the vectors row-wise to improve performance here but
        // since this is only called once, it was deemed that the potential performance
        // increase (through cache-locality) was not worth it to chnage from the convention.
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
    vOut.close();
        
    // De-allocating memory
    delete[] u;
    delete[] v;
    delete[] u_next;
    delete[] v_next;
}