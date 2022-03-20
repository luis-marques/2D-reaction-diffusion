#ifndef REACT_DIFF
#define REACT_DIFF

#include <omp.h>

/**
 * @class ReactionDiffusion
 * @author Luis Marques
 * @file ReactionDiffusion.h
 * @brief 
 */
class ReactionDiffusion {
    
    // Encapsulating solution fields u,v and other problem parameters.
    // Ensures they can only be accessed and modified by class methods.
    private:
    
        // Problem parameters (should be static const??)
        double dt;
        int T;
        int Nx;         // Number of grid-points/nodes along x direction
        int Ny;         // Number of grid-points/nodes along y direction
        double a;
        double b;
        double mu1;
        double mu2;
        double eps;
        
        // Variables defined to improve performance.
        double recip_a;
        double u_grad_coef;
        double v_grad_coef;
        double u_squared;
        double dt_eps;
        double b_over_a;
        int nr_timesteps;
                
        // Solution fields.
        double* u;
        double* v;
        double* u_prev;
        double* v_prev;
        int Ly_index;
        int Lx_index;
        //double* u_next;
        //double* v_next;
        
        
    // Methods are public as these must be accessed outside.
    public:
        void SetParameters(const double& arg_dt, const int& arg_T,
                           const int& arg_Nx, const int& arg_Ny,
                           const double& arg_a, const double& arg_b,
                           const double& arg_mu1, const double& arg_mu2,
                           const double& arg_eps);

        /// Sets the initial state of the u,v solutions field. 
        void SetInitialConditions();
        
        /// Solves the problem using time integration.
        void TimeIntegrate();
        
        /// Saves
        void Terminate();
        
};

#endif