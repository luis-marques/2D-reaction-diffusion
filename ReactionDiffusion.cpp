#include "ReactionDiffusion.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>

using namespace std;

// (Nx+1) rows and Nx*Ny columns - symmetric banded matrix
void ReactionDiffusion::PopulateA(double* MATRIX, const int& ROWS, const int& COLUMNS) {
    
    
    //cout << "Populating A" << endl;
    //cout << "Rows " << ROWS << "       COLUMNS " << COLUMNS << endl;
    
    // ROWS = Nx+1 [i is row index]  - Nx = (ROWS-1)
    // COLUMNS = Nx*Ny [j is column index]
    int aNX = ROWS-1;
    for (int j = 0; j < COLUMNS; ++j) {
       for (int i = 0; i < ROWS; ++i) {
           
           // 0th row all ones
           if (i == 0) {
               MATRIX[i + j*ROWS] = 1.0; // Top Upper diagonal
           }
           
           // Row above main diagonal
           else if (i == (ROWS-2)) {
               
                // If the ROWS-1=Nx column, then should be 0
                if ( (j % aNX) == 0) {
                   MATRIX[i + j*ROWS] = 0.0;
                }
                
                // Else is 1
                else{
                    MATRIX[i + j*ROWS] = 1.0; // 1st Upper diagonal
                }
            }
           
           // Main diagonal
           else if (i == aNX) {
               
                // First and last columns
                if ( (j == 0) || (j == COLUMNS-1) ) {
                    MATRIX[i + j*ROWS] = -2.0;
                }
                
                // Columns less than Nx, columns more than Nx*Nx - Nx (first and last row)
                // OR end and start of every row
                else if ( (j < aNX) || (j > COLUMNS - aNX) || (j % aNX == 0) || ( (j % aNX) == 1) ) {
                    MATRIX[i + j*ROWS] = -3.0;
                }
               
                else {
                    MATRIX[i + j*ROWS] = -4.0; // Main diagonal
                }
            }
           
           
           
           // All other rows
           else {
               MATRIX[i + j*ROWS] = 0.0; // (Nx-2) rows of zeros
           }
       }
    }
};


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
        
//        // A = nabla^2
//        myA = new double[(Nx+1) * Nx*Ny];
//        B_u = new double[(Nx+1) * Nx*Ny];
//        B_v = new double[(Nx+1) * Nx*Ny];
//        
//        ReactionDiffusion::PopulateA(myA, (Nx+1), Nx*Ny);
//        
//        ofstream vOut("matrixA.csv", ios::out | ios::trunc);
//        
//        cout << "File opened successfully" << endl;
//        
//        for (int i = 0; i < Nx+1; ++i) {
//            for (int j = 0; j < Nx*Ny; ++j) {
//                vOut << myA[i + j*(Nx+1)] << ",";
//            }
//            vOut << endl;
//        }
//        cout << "Closed myA matrix file fine." << endl;
//        vOut.close();
        
    
//        // Populate B_u and B_v
//        for (int j = 0; j < Nx*Ny; ++j) {
//           for (int i = 0; i < (Nx+1); ++i) {
//               
//                B_u[i + j*(Nx+1)] = mu1 * dt * myA[i + j*(Nx+1)]; // divide by h^2 = divide by 1
//                B_v[i + j*(Nx+1)] = mu2 * dt * myA[i + j*(Nx+1)];
//                
//                // Identity matrix -> add 1 to main diag
//                if (i == Nx) {
//                    B_u[i + j*(Nx+1)] += 1.0;
//                    B_v[i + j*(Nx+1)] += 1.0;
//                }
//            }
//        }


        // Add check in case T is not properly divisible by dt (check remainder)
        nr_timesteps = T / dt;
        cout << "nr_timesteps = " << nr_timesteps << "; T = " << T << "; dt = " << dt << endl; 
        
        cout << "floor(Ny/2) = " << floor(Ny/2.0) << "; Ny = " << Ny << endl;
        cout << "ceil(Nx/2) = " << ceil(Nx/2.0) << "; Nx = " << Nx << endl;
        
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
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            if (j > floor(Ny/2.0)) {
                u[i + j*Nx] = 1.0;
            }
            else{
                u[i + j*Nx] = 0.0;
            }
            
            if (i < ceil(Nx/2.0)) {
                v[i + j*Nx] = a/2.0; // make new variable = a/2 to pre-compute this!
            }
            else{
                v[i + j*Nx] = 0;
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
//    double* dummy = new double[Nx*Ny];
//    double* dummy2 = new double[Nx*Ny];

//    int problem_rank = Nx*Ny;
    
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
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+1 + j*Nx] + v[i+(j+1)*Nx] - 3*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }
                
                // along y == Ly
                else if (j==(Ny-1)){
                    u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j-1)*Nx] - 3*u[i+j*Nx]) + dt * f1[i+j*Nx];
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+1 + j*Nx] + v[i+(j-1)*Nx] - 3*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }
                
                
                // Central points
                else {
                    u[i+j*Nx] = u[i+j*Nx] + mu1*dt*(u[i+1 + j*Nx] + u[i-1 + j*Nx] + u[i+(j+1)*Nx] + u[i+(j-1)*Nx] - 4*u[i+j*Nx]) + dt * f1[i+j*Nx];
                    v[i+j*Nx] = v[i+j*Nx] + mu2*dt*(v[i+1 + j*Nx] + v[i+1 + j*Nx] + v[i+(j+1)*Nx] + v[i+(j-1)*Nx] - 4*v[i+j*Nx]) + dt * f2[i+j*Nx];
                }

            }
        }
        
//        // Get f vectors for time-step n
//        f_functions();
//        
//        // y <- alpha*A*x + beta*y (dsbmv)(UPLO, N, K, alpha, A, lda, x, incx, beta, y, incy)
//        F77NAME(dsbmv)('U', problem_rank, Nx, 1.0, B_u, (Nx+1), u, 1, 0.0, u, 1); // u^{n+1} = B_u*u^{n}
//        F77NAME(dsbmv)('U', problem_rank, Nx, 1.0, B_v, (Nx+1), v, 1, 0.0, v, 1); // v^{n+1} = B_v*v^{n}
//
////        F77NAME(dcopy)(Nx*Ny, dummy,1, u, 1);
////        F77NAME(dcopy)(Nx*Ny, dummy2,1, u, 1);
//
//        // y <- ax + y (daxpy)(N, x, incx, y, incy)
//        F77NAME(daxpy)(problem_rank, 1.0, f1, 1, u, 1); // u^{n+1} += f_1
//        F77NAME(daxpy)(problem_rank, 1.0, f2, 1, v, 1); // v^{n+1} += f_2
////        for (int i = 0; i < Nx*Ny; ++i) {
////            u[i] += f1[i];
////            v[i] += f2[i];
////        }
        
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
//    delete[] myA;
//    delete[] B_u;
//    delete[] B_v;
    delete[] f1;
    delete[] f2;
};