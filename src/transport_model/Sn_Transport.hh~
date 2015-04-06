#ifndef Lumped_Linear_Discontinuous_hh
#define Lumped_Linear_Discontinuous_hh

#include <vector>

using namespace std;

class Sn_Transport
{
public:

    // Creator
    Sn_Transport(vector<double> &internal_source,
                 vector<double> &boundary_sources,
                 vector<string> &boundary_conditions,
                 vector<double> &cell_length,
                 vector<double> &sigma_t,
                 vector<double> &sigma_s,
                 vector<double> &nu_sigma_f,
                 vector<double> &chi,
                 vector<double> &ordinates,
                 vector<double> &weights,
                 unsigned &number_of_cells,
                 unsigned &number_of_groups,
                 unsigned &number_of_ordinates,
                 unsigned &max_iterations,
                 double &tolerance);

    // implementation
    
    void lumped_linear_discontinuous(vector<double> &psi,
                                     vector<double> &leakage,
                                     unsigned &iterations);

    void psi_to_phi(vector<double> &phi,
                    vector<double> &psi);
    
protected:

    // vector operators
    
    void calculate_leakage(vector<double> &psi,
                           vector<double> &leakage);

    void calculate_source(vector<double> &q,
                          vector<double> &phi);

    void lumped_linear_discontinuous_sweep(vector<double> &psi,
                                           vector<double> &q);
    
    void check_convergence(bool &converged,
                           vector<double> &phi,
                           vector<double> &phi_old,
                           vector<double> &error_phi,
                           vector<double> &error_phi_old);
    
    
    // data
    
    vector<double> &internal_source;
    vector<double> &boundary_sources;
    vector<string> &boundary_conditions;
    vector<double> &cell_length;
    vector<double> &sigma_t;
    vector<double> &sigma_s;
    vector<double> &nu_sigma_f;
    vector<double> &chi;
    vector<double> &ordinates;
    vector<double> &weights;
    unsigned &number_of_cells;
    unsigned &number_of_groups;
    unsigned &number_of_ordinates;
    unsigned number_of_nodes;
    unsigned &max_iterations;
    double &tolerance;
};

#endif
