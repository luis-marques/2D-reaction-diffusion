#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>

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
        
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        u_prev = new double[Nx*Ny];
        v_prev = new double[Nx*Ny];
        
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
    
    int Lx = (Nx - 1);
    int Ly = (Ny - 1);
    
    cout << "Lx=" << Lx << "; Ly=" << Ly << endl;
    cout << "Ly/2 = " << Lx/2.0 << "; Lx/2 = " << Ly/2.0 << endl;
    cout << "floor(Ny/2) = " << floor(Ny/2.0) << "; Ny = " << Ny << endl;
    cout << "ceil(Nx/2) = " << ceil(Nx/2.0) << "; Nx = " << Nx << endl;
        
    //#pragma omp parallel
    for (int i = 0; i < Nx; ++i) {
      //  #pragma omp for nowait
        for (int j = 0; j < Ny; ++j) {
            
            if (j > Ly/2.0) {
                u[i + j*Nx] = 1.0;
            }
            else{
                u[i + j*Nx] = 0.0;
            }
            
            if (i < Lx/2.0) {
                v[i + j*Nx] = a/2.0; // make new variable = a/2 to pre-compute this!
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
        
        // Copying u,v to u_prev, v_prev
        //#pragma omp parallel for \
        //    schedule(static)
        for (int k = 0; k < Nx*Ny; ++k) {
            u_prev[k] = u[k];
            v_prev[k] = v[k];
            
            //u_squared = u_prev[k] * u_prev[k];
            
            u[k] += dt_eps * u_prev[k] * (1.0 - u_prev[k]) * (u_prev[k] - v_prev[k] * recip_a - b_over_a);
            v[k] += dt * (u_prev[k] * u_prev[k] * u_prev[k] - v_prev[k]);

        }

        // U !!!!!!!!!
        
        // Along same column j=0 (3 below)
        // Corner (0, 0) (i=0, j=0)
        u[0] += u_grad_coef * (u_prev[1] + u_prev[Nx] - 2*u_prev[0]);

        // Along (y==0) (0<i<Nx, j=0)
        for (int i = 1; i < (Nx-1); i++) {
            u[i] += u_grad_coef*(u_prev[i+1] + u_prev[i-1] + u_prev[i + Nx] - 3*u_prev[i]);
        }

        // Corner (Lx, 0) (i=Nx-1, j=0) Lx_index = (Nx-1)
        u[Lx_index] += u_grad_coef*(u_prev[Lx_index - 1] + u_prev[Lx_index + Nx] - 2*u_prev[Lx_index]);
        
        
        
        
        // Central points (0<i<Nx-1, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); ++j) {
            for (int i = 1; i < (Nx-1); ++i) {
                u[i+j*Nx] += u_grad_coef*(u_prev[i+1 + j*Nx] + u_prev[i-1 + j*Nx] + u_prev[i+(j+1)*Nx] + u_prev[i+(j-1)*Nx] - 4*u_prev[i+j*Nx]);
            }
        }
        
        // Along (x==0) (i=0, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); j++) {
            u[j*Nx] += u_grad_coef*(u_prev[1 + j*Nx] + u_prev[(j+1)*Nx] + u_prev[(j-1)*Nx] - 3*u_prev[j*Nx]);
        }
        
    
        // Along (x==Lx) (i=Nx-1, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); j++) {
            u[Lx_index + j*Nx] += u_grad_coef*(u_prev[Lx_index - 1 + j*Nx] + u_prev[Lx_index + (j+1)*Nx] + u_prev[Lx_index + (j-1)*Nx] - 3*u_prev[Lx_index + j*Nx]);
        }
        
        
        
        
        // Along same column j=Ny-1 (3 below)
        // Corner (0, Ly) (i=0, j=Ny-1) Ly_index = Nx*(Ny-1)
        u[Ly_index] += u_grad_coef*(u_prev[Ly_index + 1] + u_prev[Ly_index - Nx] - 2*u_prev[Ly_index]);
        
        // Along (y==Ly) (0<i<Nx, j=Ny-1)
        for (int i = 1; i < (Nx-1); i++) {
            u[i + Ly_index] += u_grad_coef*(u_prev[i+1 + Ly_index] + u_prev[i-1 + Ly_index] + u_prev[i - Nx + Ly_index] - 3*u_prev[i + Ly_index]);
        }
        
        // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
        u[Lx_index + Ly_index] += u_grad_coef*(u_prev[Lx_index + Ly_index - 1] + u_prev[Lx_index + Ly_index - Nx] - 2*u_prev[Lx_index + Ly_index]);


        // V !!!!!!!!!!!!
        
        // Along same column j=0 (3 below)
        // Corner (0, 0) (i=0, j=0)
        v[0] += v_grad_coef * (v_prev[1] + v_prev[Nx] - 2*v_prev[0]);
        
        // Along (y==0) (0<i<Nx, j=0)
        for (int i = 1; i < (Nx-1); i++) {
            v[i] += v_grad_coef*(v_prev[i+1] + v_prev[i-1] + v_prev[i + Nx] - 3*v_prev[i]);
        }
        
        // Corner (Lx, 0) (i=Nx-1, j=0) Lx_index = (Nx-1)
        v[Lx_index] += v_grad_coef*(v_prev[Lx_index - 1] + v_prev[Lx_index + Nx] - 2*v_prev[Lx_index]);


        // Along (x==Lx) (i=Nx-1, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); j++) {
            v[Lx_index + j*Nx] += v_grad_coef*(v_prev[Lx_index - 1 + j*Nx] + v_prev[Lx_index + (j+1)*Nx] + v_prev[Lx_index + (j-1)*Nx] - 3*v_prev[Lx_index + j*Nx]);
        }
        
        // Central points (0<i<Nx-1, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); ++j) {
            for (int i = 1; i < (Nx-1); ++i) {
                v[i+j*Nx] += v_grad_coef*(v_prev[i+1 + j*Nx] + v_prev[i-1 + j*Nx] + v_prev[i+(j+1)*Nx] + v_prev[i+(j-1)*Nx] - 4*v_prev[i+j*Nx]);
            }
        }

        // Along (x==0) (i=0, 0<j<Ny-1)
        for (int j = 1; j < (Ny-1); j++) {
            v[j*Nx] += v_grad_coef*(v_prev[1 + j*Nx] + v_prev[(j+1)*Nx] + v_prev[(j-1)*Nx] - 3*v_prev[j*Nx]);
        }
        
        // Along same column j=Ny-1 (3 below)
        // Corner (0, Ly) (i=0, j=Ny-1) Ly_index = Nx*(Ny-1)
        v[Ly_index] += v_grad_coef*(u_prev[Ly_index + 1] + v_prev[Ly_index - Nx] - 2*v_prev[Ly_index]);
        
        // Along (y==Ly) (0<i<Nx, j=Ny-1)
        for (int i = 1; i < (Nx-1); i++) {
            v[i + Ly_index] += v_grad_coef*(v_prev[i+1 + Ly_index] + v_prev[i-1 + Ly_index] + v_prev[i - Nx + Ly_index] - 3*v_prev[i + Ly_index]);
        }
        
        // Corner (Lx, Ly) (i=Nx-1, j=Ny-1)
        v[Lx_index + Ly_index] += v_grad_coef*(v_prev[Lx_index + Ly_index - 1] + v_prev[Lx_index + Ly_index - Nx] - 2*v_prev[Lx_index + Ly_index]);

        
                
        
        
        //#pragma omp parallel
//        //#pragma omp parallel for collapse(2)
//        for (int j = 0; j < Ny; ++j) {
//            //#pragma omp for nowait schedule(static)
//            for (int i = 0; i < Nx; ++i) {


        
        if (timestep % 10000 == 0) { 
            cout << "Timestep " << timestep << endl;
        }
    
    }
    
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
    delete[] u_prev;
    delete[] v_prev;
}