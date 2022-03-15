#ifndef REACT_DIFF
#define REACT_DIFF

class ReactionDiffusion {
    
    // Should be private?
    //private:
    public:
    
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
        
        // A = nabla^2
        double* A;
        double* B_u;
        double* B_v;

        double* f1;
        double* f2;
        
        // Constructor as want to set the variables as constant and that can't be done if they are 
        // initialized in SetParams
        void SetParameters(const double& arg_dt, const int arg_T,
                           const int arg_Nx, const int arg_Ny,
                           const double arg_a, const double arg_b,
                           const double arg_mu1, const double arg_mu2,
                           const double arg_eps);

        
        void SetInitialConditions();
        
        void TimeIntegrate();
        
        void f_f1();
        void f_f2();
        
        void Terminate();
};

#endif