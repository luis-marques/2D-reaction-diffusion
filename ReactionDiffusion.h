#ifndef REACT_DIFF
#define REACT_DIFF

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
        int T;                  // Total integration time.
        int Nx;                 // Number of grid-points/nodes along x direction.
        int Ny;                 // Number of grid-points/nodes along y direction
        double a;               // Parameter 'a' of Barkley model problem.
        double b;               // Parameter 'b' of Barkley model problem.
        double mu1;             // Diffusion coefficient for 'u'.
        double mu2;             // Diffusion coefficient for 'v'.
        double eps;             // Paramter 'eps' of Barkley model problem.
        
        // Variables defined to improve code performance.
        double recip_a;         // Reciprobal of parameter 'a' (i.e. =1/a)
        double u_grad_coef;     // (i.e. = mu1*dt / h^2, where we take h=1)
        double v_grad_coef;     // (i.e. = mu2*dt / h^2, where we take h=1)
        double dt_eps;          // (i.e. = dt*eps
        double b_over_a;        // (i.e. = b/a)
        double half_a;          // (i.e. = a/2.0)
        int nr_timesteps;       // (i.e. = T/dt)
        int Lx_index;           // ( = Nx-1, explained further in source file) 
        int Ly_index;           // ( = Nx*(Ny-1, explained further in source file) 
                
        // Solution fields.
        double* u;              // u^{n} -> Solution field for 'u' at time-step n.
        double* v;              // v^{n} -> Solution field for 'v' at time-step n.
        double* u_next;         // u^{n+1} -> Solution field for 'u' at time-step n+1.
        double* v_next;         // v^{n+1} -> Solution field for 'v' at time-step n+1.
        
        
        
    // Class methods are public, as these can be called outside from outside the class.
    public:
    
        /// Takes parsed arguments, stores them in class and calculates
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
        
        /// De-allocates all dynamically allocated memory.
        ~ReactionDiffusion();   
        
};

#endif