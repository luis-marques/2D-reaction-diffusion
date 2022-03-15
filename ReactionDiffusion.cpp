#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>

using namespace std;

void ReactionDiffusion::SetParameters(
                    const double& arg_dt, const int arg_T,
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
        B_u = new double[Ny * Nx*Ny];
        B_v = new double[Ny * Nx*Ny];

        // Add check in case T is not properly divisible by dt (check remainder)
        nr_timesteps = T / dt;
        cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << endl; 
        
        cout << "ceil(Ny/2) = " << ceil(Ny/2) << "; Ny = " << Ny << endl; 
        
        f1 = new double[Nx*Ny];
        f2 = new double[Nx*Ny];
        
        
        // For Debugging
        cout << "\tParameters of PDE problem to solve, from the ReactionDiffusion class:" << endl;
        cout << "\tTime-step (dt)" << right << setw(30) << setfill(' ') << dt << endl;
        cout << "\tIntegration time (T)" << right << setw(30) << setfill(' ') << T << endl;
        cout << "\tNx" << right << setw(30) << setfill(' ') << Nx << endl;
        cout << "\tNy" << right << setw(30) << setfill(' ') << Ny << endl;
        cout << "\ta" << right << setw(30) << setfill(' ') << a << endl;
        cout << "\tb" << right << setw(30) << setfill(' ') << b << endl;
        cout << "\tmu1" << right << setw(30) << setfill(' ') << mu1 << endl;
        cout << "\tmu2" << right << setw(30) << setfill(' ') << mu2 << endl;
        cout << "\tepsilon" << right << setw(30) << setfill(' ') << eps << endl;
        
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

void ReactionDiffusion::f_f1() {
    for (int i = 0; i < Nx*Ny; ++i) {
        f1[i] = eps * u[i] * (1 - u[i]) * (u[i] - (v[i] + b) / a);
    }
};


void ReactionDiffusion::f_f2() {
    for (int i = 0; i < Nx*Ny; ++i) {
        f2[i] = u[i] * u[i] * u[i] - v[i] ; // also test with pow() to see performance increase
    }
};


void ReactionDiffusion::TimeIntegrate() {
    cout << "Starting numerical solving of PDE." << endl;
    
    for (int timestep = 0; timestep < nr_timesteps; ++timestep) {
        
        // Get f vectors for time-step n
        f_f1();
        f_f2();
        
        if (timestep % 10000 == 0) { 
            cout << "Timestep " << timestep << endl;
        }
    
    }
    
    cout << "Finished solving PDE (from t_i = 0 to t_f = T)." << endl;
    
};


void ReactionDiffusion::Terminate() {
    
    cout << "Writting output of simulation to file 'output.txt'." << endl;
    
    ofstream vOut("output.txt", ios::out | ios::trunc);
    
    if (vOut.is_open()) {
        cout << "File opened successfully" << endl;
    
        for (int j = 0; j < Ny; ++j){
            for (int i = 0; i < Nx; ++i) {
                vOut << i << " " 
                     << j << " "
                     << u[i + j*Ny] << " "
                     << v[i + j*Ny] << endl;
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
    delete[] A;
    delete[] B_u;
    delete[] B_v;
    delete[] f1;
    delete[] f2;
};