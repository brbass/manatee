#include "SP1_Transport.hh"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include <Amesos.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>

#include "Data.hh"
#include "Mesh.hh"

namespace transport_ns
{
    using namespace std;
    using namespace data_ns;
    using namespace mesh_ns;
    
    SP1_Transport::
    SP1_Transport(Data &data,
                  Mesh &mesh,
                  string problem_type):
        data_(data),
        mesh_(mesh),
        problem_type_(problem_type)
    {
        data_.check();
        mesh_.check();
        
        number_of_edges_ = mesh_.number_of_cells() + 1;
        
        num_global_elements_ = number_of_edges_ * data_.number_of_groups();

        comm_ = unique_ptr<Epetra_MpiComm> (new Epetra_MpiComm(MPI_COMM_WORLD));
        map_ = unique_ptr<Epetra_Map> (new Epetra_Map(num_global_elements_, index_base_, *comm_));
        
        initialize_d();
        initialize_matrix();
        initialize_lhs();
        initialize_rhs();
        initialize_problem();
        initialize_solver();
    }
    
    int SP1_Transport::
    solve()
    {
        solver_->SymbolicFactorization();
        
        solver_->NumericFactorization();
        
        solver_->Solve();
        
        // solver_->PrintStatus();

        return 0;
    }
    
    int SP1_Transport::
    initialize_matrix()
    {
        num_my_elements_ = map_->NumMyElements();
        my_global_elements_.assign(map_->MyGlobalElements(), map_->MyGlobalElements() + num_my_elements_);
        
        num_entries_per_row_.resize(num_my_elements_, 0);
        
        for (unsigned i = 0; i < number_of_edges_; ++i)
        {
            num_entries_per_row_[i + 0] = 3 * data_.number_of_groups();
        }
        
        matrix_ = unique_ptr<Epetra_FECrsMatrix> (new Epetra_FECrsMatrix (Copy, *map_, &num_entries_per_row_[0], true));
        
        unsigned n = 0;
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            unsigned dims = mesh_.number_of_nodes() * data_.number_of_groups();
            
            Epetra_SerialDenseMatrix fill_matrix (dims, dims);
            Epetra_IntSerialDenseVector fill_vector (dims);

            double dx2 = pow(mesh_.cell_length(i), 2);
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                {
                    unsigned k1 = o1 + mesh_.number_of_nodes() * g;
                    
                    for (unsigned o2 = 0; o2 < mesh_.number_of_nodes(); ++o2)
                    {
                        unsigned k2 = o2 + mesh_.number_of_nodes() * g;
                        
                        fill_matrix(k1, k2) += dx2 / 4 * data_.sigma_t(i, g) * mesh_.stiffness(i, o1, o2);
                    }
                }
            }

            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                {
                    unsigned k1 = o1 + mesh_.number_of_nodes() * gt;

                    for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                    {
                        for (unsigned o2 = 0; o2< mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;

                            if (problem_type_.compare("forward") == 0)
                            {
                                fill_matrix(k1, k2) -= dx2 / 4 * data_.sigma_s(i, gf, gt, n) * mesh_.stiffness(i, o1, o2);
                            }
                            else if (problem_type_.compare("adjoint") == 0)
                            {
                                fill_matrix(k1, k2) -= dx2 / 4 * data_.sigma_s(i, gt, gf, n) * mesh_.stiffness(i, o1, o2);
                            }
                            else
                            {
                                cerr << "problem type not available" << endl;
                            }
                        }
                    }
                }
            }
            
            for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
            {
                for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                {
                    for (unsigned o1 = 0; o1 < mesh_.number_of_nodes(); ++o1)
                    {
                        unsigned k1 = o1 + mesh_.number_of_nodes() * gt;
                        
                        for (unsigned o2 = 0; o2 < mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                            
                            fill_matrix(k1, k2) += d(i, gf, gt) * mesh_.stiffness(i, o1, o2, 1, 1);
                        }
                        
                        fill_vector(k1) = i + o1 + number_of_edges_ * gt;
                    }
                }
            }
            if (i == 0)
            {
                for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                { 
                    for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                    {
                        unsigned o1 = 0;
                            
                        unsigned k1 = o1 + mesh_.number_of_nodes() * gt;
                            
                        for (unsigned o2 = 0; o2 < mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                                    
                            fill_matrix(k1, k2) = 0;
                        }
                    }
                }
                if (data_.boundary_condition(0).compare("dirichlet") == 0)
                {
                    for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                    {
                        unsigned o = 0;
                            
                        unsigned k = o + mesh_.number_of_nodes() * g;
                            
                        fill_matrix(k, k) = 4;
                    }
                }
                else if (data_.boundary_condition(0).compare("marshak") == 0)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    { 
                        for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                        {
                            unsigned k1 = 0 + mesh_.number_of_nodes() * gt;
                            unsigned k2 = 0 + mesh_.number_of_nodes() * gf;
                            unsigned k3 = 1 + mesh_.number_of_nodes() * gf;
                            
                            fill_matrix(k1, k2) = 2 / mesh_.cell_length(i) * d(i, gf, gt);
                            fill_matrix(k1, k3) = - 2 / mesh_.cell_length(i) * d(i, gf, gt);
                        }
                    }
                }
                else if (data_.boundary_condition(0).compare("robin") == 0)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    { 
                        for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                        {
                            unsigned k1 = 0 + mesh_.number_of_nodes() * gt;
                            unsigned k2 = 0 + mesh_.number_of_nodes() * gf;
                            unsigned k3 = 1 + mesh_.number_of_nodes() * gf;
                            
                            fill_matrix(k1, k2) = 2 / mesh_.cell_length(i) * d(i, gf, gt);
                            fill_matrix(k1, k3) = - 2 / mesh_.cell_length(i) * d(i, gf, gt);

                            if (gt == gf)
                            {
                                fill_matrix(k1, k2) += 1;
                            }
                        }
                    }
                }
                else
                {
                    cout << "boundary condition " << data_.boundary_condition(0) << " not available" << endl;
                }
            }
            
            if (i == mesh_.number_of_cells() - 1)
            {
                for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                { 
                    for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                    {
                        unsigned o1 = 1;
                            
                        unsigned k1 = o1 + mesh_.number_of_nodes() * gt;
                            
                        for (unsigned o2 = 0; o2 < mesh_.number_of_nodes(); ++o2)
                        {
                            unsigned k2 = o2 + mesh_.number_of_nodes() * gf;
                                    
                            fill_matrix(k1, k2) = 0;
                        }
                    }
                }
                if (data_.boundary_condition(1).compare("dirichlet") == 0)
                {
                    for (unsigned g = 0; g < data_.number_of_groups(); ++g)
                    {
                        unsigned o = 1;

                        unsigned k = o + mesh_.number_of_nodes() * g;
                        
                        fill_matrix(k, k) = 4;
                    }
                }
                else if (data_.boundary_condition(1).compare("marshak") == 0)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    { 
                        for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                        {
                            unsigned k1 = 1 + mesh_.number_of_nodes() * gt;
                            unsigned k2 = 0 + mesh_.number_of_nodes() * gf;
                            unsigned k3 = 1 + mesh_.number_of_nodes() * gf;
                        
                            fill_matrix(k1, k2) = - 2 / mesh_.cell_length(i) * d(i, gf, gt);
                            fill_matrix(k1, k3) = 2 / mesh_.cell_length(i) * d(i, gf, gt);
                        }
                    }
                }
                else if (data_.boundary_condition(1).compare("robin") == 0)
                {
                    for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                    { 
                        for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
                        {
                            unsigned k1 = 1 + mesh_.number_of_nodes() * gt;
                            unsigned k2 = 0 + mesh_.number_of_nodes() * gf;
                            unsigned k3 = 1 + mesh_.number_of_nodes() * gf;
                        
                            fill_matrix(k1, k2) = - 2 / mesh_.cell_length(i) * d(i, gf, gt);
                            fill_matrix(k1, k3) = 2 / mesh_.cell_length(i) * d(i, gf, gt);

                            if (gt == gf)
                            {
                                fill_matrix(k1, k3) += 1;
                            }
                        }
                    }
                }
                else
                {
                    cout << "boundary condition " << data_.boundary_condition(1) << " not available" << endl;
                }
            }

            matrix_->InsertGlobalValues(fill_vector, fill_matrix);
        }
        
        matrix_->GlobalAssemble();
        
        return 0;
    }
    
    int SP1_Transport::
    initialize_lhs()
    {
        lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
        lhs_->PutScalar(1.0);
        
        return 0;
    }
    
    int SP1_Transport::
    initialize_rhs()
    {
        rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
        rhs_->PutScalar(0.0);
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            double dx2 = pow(mesh_.cell_length(i), 2);
            
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                for (unsigned o = 0; o < mesh_.number_of_nodes(); ++o)
                {
                    unsigned k = o + i + number_of_edges_ * g;
                    
                    (*rhs_)[k] += dx2 / 4 * data_.internal_source(i, g);
                }
            }
        }
        
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            unsigned i = 0;
            unsigned o = 0;
            unsigned k = o + i + number_of_edges_ * g;
            
            (*rhs_)[k] = 4 * data_.boundary_source(0, g);
        }
        for (unsigned g = 0; g < data_.number_of_groups(); ++g)
        {
            unsigned i = mesh_.number_of_cells() - 1;
            unsigned o = 1;
            unsigned k = o + i + number_of_edges_ * g;
            
            (*rhs_)[k] = 4 * data_.boundary_source(1, g);
        }
        
        return 0;
    }

    int SP1_Transport::
    initialize_problem()
    {
        problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
        
        return 0;
    }

    int SP1_Transport::
    initialize_solver()
    {
        solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *problem_));
        
        if (solver_ == NULL)
        {
            cerr << "Specified solver \"" << solver_type_ << "\" is not available." << endl;
        }
        
        // list_.set("PrintTiming", true);
        // list_.set("PrintStatus", true);
        
        solver_->SetParameters(list_);
        
        return 0;
    }

    int SP1_Transport::
    initialize_d()
    {
        unsigned n = 0;
        
        Epetra_SerialDenseSolver solver;
        
        d_.resize(pow(data_.number_of_groups(), 2) * mesh_.number_of_cells(), 0);
        
        for (unsigned i = 0; i < mesh_.number_of_cells(); ++i)
        {
            Epetra_SerialDenseMatrix invert_matrix (data_.number_of_groups(), data_.number_of_groups());
                
            for (unsigned g = 0; g < data_.number_of_groups(); ++g)
            {
                invert_matrix(g, g) += data_.sigma_t(i, g);
            }
                
            for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
            {
                for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                {
                    if (problem_type_.compare("forward") == 0)
                    {
                        invert_matrix(gf, gt) -= data_.sigma_s(i, gf, gt, n+1);
                    }
                    else if (problem_type_.compare("adjoint") == 0)
                    {
                        invert_matrix(gf, gt) -= data_.sigma_s(i, gt, gf, n+1);
                    }
                    else
                    {
                        cerr << "problem type not available" << endl;
                    }
                }
            }
            
            solver.SetMatrix(invert_matrix);
            solver.Invert();
                
            for (unsigned gf = 0; gf < data_.number_of_groups(); ++gf)
            {
                for (unsigned gt = 0; gt < data_.number_of_groups(); ++gt)
                {
                    unsigned k = gf + data_.number_of_groups() * (gt + data_.number_of_groups() * i);

                    d_[k] = pow(n+1, 2) / ((2*n+1) * (2*n+3)) * invert_matrix(gf, gt);
                }
            }
        }
        
        return 0;
    }
}
