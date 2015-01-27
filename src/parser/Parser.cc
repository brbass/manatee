#include "Parser.hh"

#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace parser_ns
{
    using namespace std;

    Parser::
    Parser(string &folder)
    {
        folder_ = folder;
        
        system(("mkdir -p " + folder_).c_str());
    }

    void Parser::
    set_folder(string folder)
    {
        folder_ = folder;

        system(("mkdir -p " + folder_).c_str());
    }
    
    unsigned Parser::
    parse_data(unsigned &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ifstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file >> data;

            return 0;
        }
        else
        {
            cout << "unable to load " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    parse_data(double &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ifstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file >> data;

            return 0;
        }
        else
        {
            cout << "unable to load " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    parse_data(vector<unsigned> &data, string data_filename)
    {
        data.resize(0);
        
        string data_path = get_filepath(data_filename);
        
        ifstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            unsigned temp;
            
            while (data_file >> temp)
            {
                data.push_back(temp);
            }
            
            return 0;
        }
        else
        {
            cout << "unable to load " << data_filename << endl;
            
            return 1;
        }
    }
    
    unsigned Parser::
    parse_data(vector<double> &data, string data_filename)
    {
        data.resize(0);
    
        string data_path = get_filepath(data_filename);
    
        ifstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            double temp;
        
            while (data_file >> temp)
            {
                data.push_back(temp);
            }

            return 0;
        }
        else
        {
            cout << "unable to load " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    parse_data(vector<string> &data, string data_filename)
    {
        data.resize(0);
        
        string data_path = get_filepath(data_filename);
        
        ifstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            string temp;
        
            while (data_file >> temp)
            {
                data.push_back(temp);
            }

            return 0;
        }
        else
        {
            cout << "unable to load " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    write_data(bool &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    write_data(int &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    write_data(unsigned &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    write_data(double &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;

            return 1;
        }
    }

    unsigned Parser::
    write_data(string &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;

            return 1;
        }
    }
    
    unsigned Parser::
    write_data(vector<bool> &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ofstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            for (unsigned i = 0; i < data.size(); ++i)
            {
                data_file << data[i] << endl;
            }
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;
            
            return 1;
        }
    }

    unsigned Parser::
    write_data(vector<unsigned> &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ofstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            for (unsigned i = 0; i < data.size(); ++i)
            {
                data_file << data[i] << endl;
            }
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;
            
            return 1;
        }
    }
    
    unsigned Parser::
    write_data(vector<double> &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ofstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            for (unsigned i = 0; i < data.size(); ++i)
            {
                data_file << data[i] << endl;
            }
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;
            
            return 1;
        }
    }

        unsigned Parser::
        write_data(vector<string> &data, string data_filename)
        {
        string data_path = get_filepath(data_filename);
    
        ofstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            for (unsigned i = 0; i < data.size(); ++i)
            {
                data_file << data[i] << endl;
            }
            
            return 0;
        }
        else
        {
            cout << "unable to save to " << data_filename << endl;
            
            return 1;
        }
    }
    
    string Parser::
    get_filepath(string &data_filename)
    {
        return folder_ + "/" + data_filename;
    }
}
