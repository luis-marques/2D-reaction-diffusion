#ifndef REACT_DIFF
#define REACT_DIFF
// Good practice to include header guards so that the .h file is never
// include more than once in the same .cpp file.

/**
 * @class ReactionDiffusion
 * @author Luis Marques
 * @file ReactionDiffusion.h
 * @brief Encapsulates the reaction-diffusion problem, while following the
 * specifications set out in the briefing.
 */
class ReactionDiffusion {
    
    // Encapsulating solution fields u,v and other parameters/variables.
    // Ensures these variables can only be accessed and modified by class methods.
    private:
        
        // Problem parameters required by project briefing.
        double dt;              // Integration time-step.
        int    T;               // Total integration time.
        int    Nx;              // Number of grid-points/nodes along x direction.
        int    Ny;              // Number of grid-points/nodes along y direction
        double a;               // Parameter 'a' of Barkley model problem.
        double b;               // Parameter 'b' of Barkley model problem.
        double mu1;             // Diffusion coefficient for 'u'.
        double mu2;             // Diffusion coefficient for 'v'.
        double eps;             // Parameter 'eps' of Barkley model problem.
        
        // Solution fields.
        double* u;              // u^{n} -> Solution field for 'u' at time-step n.
        double* v;              // v^{n} -> Solution field for 'v' at time-step n.
        double* u_next;         // u^{n+1} -> Solution field for 'u' at time-step n+1.
        double* v_next;         // v^{n+1} -> Solution field for 'v' at time-step n+1.
        
        // Variables defined to improve code performance.
        double recip_a;         // Reciprocal of parameter 'a' (i.e. = 1/a)
        double u_grad_coef;     // (i.e. = mu1*dt / h^2, where we take h=1)
        double v_grad_coef;     // (i.e. = mu2*dt / h^2, where we take h=1)
        double dt_eps;          // (i.e. = dt*eps
        double b_over_a;        // (i.e. = b/a)
        double half_a;          // (i.e. = a/2.0)
        int    nr_timesteps;    // (i.e. = T/dt)
        

        // Recall 'u','v' are stored in column-major format, 'Nx' is the number
        // of rows and 'Ny' the number of columns. These 2 variables help with
        // code readability and slightly with performance.
        
        // Provides number of elements that must be traversed to go from
        // u_{0,j} to u_{Lx,j} (i.e. = Nx-1) 
        int    Lx_index;
        
        // Provides number of elements that must be traversed to go from
        // u_{i,0} to u_{i,Ly} (i.e. = Nx*(Ny-1))
        int    Ly_index; 
      
          


    // Class methods are public, as these can be called from outside the class.
    public:
    
        /// Takes parsed arguments, stores them and calculates various variables for performance improvements.
        void SetParameters(const double& arg_dt, const int& arg_T,
                           const int& arg_Nx, const int& arg_Ny,
                           const double& arg_a, const double& arg_b,
                           const double& arg_mu1, const double& arg_mu2,
                           const double& arg_eps);

        /// Sets the initial state of the 'u','v' solution fields. 
        void SetInitialConditions();
        
        /// Solves the problem using time integration (from t=0 to t=T).
        void TimeIntegrate();
        
        /// Saves the result of the simulation to a .txt file.
        void SaveToFile();
        
        /// Deconstructor de-allocates all dynamically allocated memory.
        ~ReactionDiffusion();   
        
};

#endif