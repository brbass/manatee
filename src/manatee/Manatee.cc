#include <iostream>
#include <string>

//#include <mpi.h>

#include "Parser.hh"
#include "Transport_Model.hh"

void run_problem(string input_folder)
{
    using namespace std;
    using namespace monte_carlo_ns;
    using namespace parser_ns;
    using namespace transport_ns;
    
    Parser parser(input_folder);

    string transport_type;
    
    parser.parse_data(transport_type, "transport_type");

    unsigned
        number_of_cells,
        number_of_groups;
    double
        side_length;
    vector<double>
        internal_source,
        boundary_sources,
        sigma_t,
        sigma_s;
    vector<string>
        boundary_conditions;
        
    parser.parse_data(number_of_cells, "number_of_cells");
    parser.parse_data(number_of_groups, "number_of_groups");
    parser.parse_data(side_length, "side_length");
    parser.parse_data(internal_source, "internal_source");
    parser.parse_data(boundary_sources, "boundary_sources");
    parser.parse_data(sigma_t, "sigma_t");
    parser.parse_data(sigma_s, "sigma_s");
    parser.parse_data(boundary_conditions, "boundary_conditions");
    
    if (transport_type.compare("fem_transport") == 0)
    {
        parser.parse_data(number_of_scattering_moments, "number_of_scattering_moments");
        parser.parse_data(number_of_groups, "nu_sigma_f");
        parser.parse_data(side_length, "chi");
        parser.parse_data(side_length, "chi");
        
        FEM_Transport(number_of_cells,
                      number_of_groups,
                      number_of_scattering_moments,
                      side_length,
                      internal_source,
                      boundary_sources,
                      sigma_t,
                      sigma_s,
                      nu_sigma_f,
                      chi,
                      boundary_conditions,
                      max_iterations,
                      tolerance,
                      geometry);
        
    }
    else if (transport_type.compare("sn_transport") == 0)
    {
    parser.parse_data();
    parser.parse_data();
    parser.parse_data();
        
    }
    else if (transport_type.compare("sp1_transport") == 0)
    {
    parser.parse_data();
    parser.parse_data();
    parser.parse_data();
        
    }
    else
    {
        cout << "transport type not available" << endl;
    }

}

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace transport_ns;
    
    if (argc != 2)
    {
        cout << "usage: manatee [input_folder]" << endl;
        return 1;
    }

    //MPI_Init(&argc, &argv);
    
    string input_folder = argv[1];
    
    run_problem(input_folder);
    
    //MPI_Finalize();
}
