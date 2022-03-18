#ifndef REACT_DIFF
#define REACT_DIFF

/**
 * @class ReactionDiffusion
 * @author Luis Marques
 * @file ReactionDiffusion.h
 * @brief 
 */
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
        
        
        double recip_a;
        double u_grad_coef;
        double v_grad_coef;
        double u_squared;
        double dt_eps;
        double b_over_a;
        
        int nr_timesteps;
                
        // Solutions
        double* u;
        double* v;
        double* u_prev;
        double* v_prev;
        
        
    public:
        // Constructor as want to set the variables as constant and that can't be done if they are 
        // initialized in SetParams
        void SetParameters(const double& arg_dt, const int& arg_T,
                           const int& arg_Nx, const int& arg_Ny,
                           const double& arg_a, const double& arg_b,
                           const double& arg_mu1, const double& arg_mu2,
                           const double& arg_eps);

        
        void SetInitialConditions();
                
        void TimeIntegrate();
                
        void Terminate();
};

#endif