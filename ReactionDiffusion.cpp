#include "ReactionDiffusion.h"

#include <iostream>
#include <cmath>


void ReactionDiffusion::SetParameters(
                    const double arg_dt, const int arg_T,
                    const int arg_Nx, const int arg_Ny,
                    const double arg_a, const double arg_b,
                    const double arg_mu1, const double arg_mu2,
                    const double arg_eps) {
                        
        dt = arg_dt;
        T = arg_T;
        Nx = arg_Nx;
        Ny = arg_Ny;
        a = arg_a;
        b = arg_b;
        mu1 = arg_mu1;
        mu2 = arg_mu2;
        eps = arg_eps;
        
        
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        
        // A = nabla^2
        A = new double[Ny * Nx*Ny];
        B = new double[Ny * Nx*Ny];

        // Add check in case T is not properly divisible by dt (check remainder)
        nr_timesteps = T / dt;
        std::cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << std::endl; 
        
        std::cout << "ceil(Ny/2) = " << ceil(Ny/2) << "; Ny = " << Ny << std::endl; 
        
        f1 = new double[Nx*Ny];
        f2 = new double[Nx*Ny];
};
    
void ReactionDiffusion::SetInitialConditions() {
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            if (j > ceil(Ny/2)) {
                u[i + j*Ny] = 1;
            }
            else{
                u[i + j*Ny] = 0;
            }
            
            if (i > ceil(Nx/2)) {
                v[i + j*Ny] = a/2; // make new variable = a/2 to pre-compute this!
            }
            else{
                v[i + j*Ny] = 0;
            }
            
            
        }
    }
};
    
void ReactionDiffusion::TimeIntegrate (){
    std::cout << "Time integrate" << std::endl;
    
    for (int timestep = 0; timestep < nr_timesteps; ++timestep) {
    
    }
    
};