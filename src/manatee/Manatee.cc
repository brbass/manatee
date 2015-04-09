#include <iostream>
#include <string>

#include <mpi.h>

#include "Transport_Model.hh"

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace transport_ns;
    
    if (argc != 2)
    {
        cout << "usage: manatee [input_folder]" << endl;
        return 1;
    }

    MPI_Init(&argc, &argv);
    
    string input_folder = argv[1];
    
    Transport_Model transport_model(input_folder);

    MPI_Finalize();
}
