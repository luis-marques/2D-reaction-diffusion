#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>

using namespace std;

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
        
        
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        
        f1 = new double[Nx*Ny];
        f2 = new double[Nx*Ny];

        // Add check in case T is not properly divisible by dt (check remainder)
        nr_timesteps = T / dt;
        cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << endl; 
        
//        cout << "Ly/2 = " << (Ny-1)/2.0 << "; Lx/2 = " << (Nx-1)/2.0 << endl;
//        cout << "floor(Ny/2) = " << floor(Ny/2.0) << "; Ny = " << Ny << endl;
//        cout << "ceil(Nx/2) = " << ceil(Nx/2.0) << "; Nx = " << Nx << endl;
//        
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
        
};
    
void ReactionDiffusion::SetInitialConditions() {
    
    int Lx = (Nx - 1);
    int Ly = (Ny - 1);
    
    cout << "Lx=" << Lx << "; Ly=" << Ly << endl;
    cout << "Ly/2 = " << Lx/2.0 << "; Lx/2 = " << Ly/2.0 << endl;
    cout << "floor(Ny/2) = " << floor(Ny/2.0) << "; Ny = " << Ny << endl;
    cout << "ceil(Nx/2) = " << ceil(Nx/2.0) << "; Nx = " << Nx << endl;
        
    
    for (int i = 0; i < Nx; ++i) {
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
};

// dt * f1 directly & dt * f2 directly
void ReactionDiffusion::f_functions() {

    for (int i = 0; i < Nx*Ny; ++i) {
        f1[i] = eps * u[i] * (1.0 - u[i]) * (u[i] - (v[i] + b) / a);
        
        f2[i] = u[i] * u[i] * u[i] - v[i] ; // also test with pow() to see performance increase
    }
    
};


void ReactionDiffusion::TimeIntegrate() {
    cout << "Starting numerical solving of PDE." << endl;
    
    for (int timestep = 0; timestep < nr_timesteps; ++timestep) {
        
        // Get f vectors for time-step n
        f_functions();
        
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
            
                if (i==0){
                    
                    // Corner (0,0)
                    if (j==0) {
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i+(j+1)*Nx] - 2*u[i+j*Nx]) + dt *f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+(j+1)*Nx] - 2*v[i+j*Nx]) + dt *f2[i+j*Nx];
                    }
                    
                    // Other corner (0, Ly)
                    else if (j == (Ny-1)) {
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i+(j-1)*Nx] - 2*u[i+j*Nx]) + dt *f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+(j-1)*Nx] - 2*v[i+j*Nx]) + dt *f2[i+j*Nx];
                    }
                    
                    // along x == 0
                    else { 
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 3*u[i+j*Nx]) + dt *f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 3*v[i+j*Nx]) + dt *f2[i+j*Nx];
                    }
                    
                }
                
                else if (i==(Nx-1)){
                    
                    // Corner (Lx,0)
                    if (j==0) {
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i-1 + j*Nx] + u[i+(j+1)*Nx] - 2*u[i+j*Nx]) + dt *f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i-1 + j*Nx] + v[i+(j+1)*Nx] - 2*v[i+j*Nx]) + dt *f2[i+j*Nx];
                    }
                    
                    // Other corner (Lx, Ly)
                    else if (j == (Ny-1)) {
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i-1 + j*Nx] + u[i+(j-1)*Nx] - 2*u[i+j*Nx]) + dt *f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i-1 + j*Nx] + v[i+(j-1)*Nx] - 2*v[i+j*Nx]) + dt *f2[i+j*Nx];
                    }
                    
                    // along x == Lx
                    else { 
                        u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i-1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 3*u[i+j*Nx]) + dt * f1[i+j*Nx];
                        v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i-1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 3*v[i+j*Nx]) + dt * f2[i+j*Nx];
                    }
                    
                }
                
                // along y == 0
                else if (j==0){
                    u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j+1)*Nx] - 3*u[i+j*Nx]) + dt * f1[i+j*Nx];
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i-1 + j*Nx] + v[i+(j+1)*Nx] - 3*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }
                
                // along y == Ly
                else if (j==(Ny-1)){
                    u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j-1)*Nx] - 3*u[i+j*Nx]) + dt * f1[i+j*Nx];
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i-1 + j*Nx] + v[i+(j-1)*Nx] - 3*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }
                
                
                // Central points
                else {
                    u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 4*u[i+j*Nx]) + dt * f1[i+j*Nx];
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i-1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 4*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }

            }
        }
        
        if (timestep % 10000 == 0) { 
            cout << "Timestep " << timestep << endl;
        }
    
    }
    
    // delete[] dummy;
    cout << "Finished solving PDE (from t_i = 0 to t_f = T)." << endl;
    
};

/*
 * 
 * 
 * 
 * 
 * Commands to plot u,v on domain:
 * $ gnuplot
 * $ set pm3d at st
 * $ set view map
 * $ set cbrange[0:1]
 * 
 * $ splot 'output.txt' using 1:2:3 w l palette <- Plots 'u'
 * OR
 * $ splot 'output.txt' using 1:2:4 w l palette <- Plots 'v'
 * 
 * */
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
    delete[] f1;
    delete[] f2;
};