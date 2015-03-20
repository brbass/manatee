#include <iostream>
#include <string>

#include <mpi.h>

#include "Neutronics.hh"

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace neutronics_ns;

    if (argc != 2)
    {
        cout << "usage: manatee [input_folder]" << endl;
        return 1;
    }

    MPI_Init(&argc, &argv);
    
    string input_folder = argv[1];
    
    Neutronics neutronics(input_folder);

    MPI_Finalize();
}
