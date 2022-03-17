#ifndef REACT_DIFF
#define REACT_DIFF


//#include "cblas.h"
//
// Macro to append an underscore to the Fortran symbols (to mimic behaviour of some Fortran compilers)
#define F77NAME(x) x##_

// Prototypes for the BLAS routines. Using "C" naming to avoid name mangling
extern "C" {
    void F77NAME(dsbmv) (const char& uplo, const int& N, const int& BW,
                         const double& alpha, double* A, const int& lda,
                         double* x, const int& incx,
                         const double& beta, double* y, const int& incy);
    void F77NAME(daxpy) (const int& N, const double& alpha, 
                         const double* x, const int& incx,
                         const double* y, const int& incy);
    void F77NAME(dcopy) (const int& N, const double* x, const int& incx,
                         const double* y, const int& incy);
}


class ReactionDiffusion {
    
    // Should be private?
    //private:
    private:
    
        // Class variables (should be static const??)
        double dt;
        int T;
        int Nx;         // Number of grid-points/nodes along x direction
        int Ny;         // Number of grid-points/nodes along y direction
        double a;
        double b;
        double mu1;
        double mu2;
        double eps;
        
        int nr_timesteps;
        
        // Remember to de-allocate this!
        
        // Solutions
        double* u;
        double* v;
        
        // myA = nabla^2
        double* myA;
        double* B_u;
        double* B_v;

        double* f1;
        double* f2;
        
        
    public:
        // Constructor as want to set the variables as constant and that can't be done if they are 
        // initialized in SetParams
        void SetParameters(const double& arg_dt, const int& arg_T,
                           const int& arg_Nx, const int& arg_Ny,
                           const double& arg_a, const double& arg_b,
                           const double& arg_mu1, const double& arg_mu2,
                           const double& arg_eps);

        
        void SetInitialConditions();
        
        void PopulateA(double* MATRIX, const int& ROWS, const int& COLUMNS);
        
        void TimeIntegrate();
        
        void f_functions();
        
        void Terminate();
};

#endif