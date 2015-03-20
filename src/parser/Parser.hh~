#ifndef Parser_hh
#define Parser_hh

#include <string>
#include <vector>

namespace parser_ns
{
    using std::string;
    using std::vector;

    class Parser
    {
    private:

        string folder_;
    
    public:

        Parser(string &folder);
    
        unsigned parse_data(unsigned &data, string data_filename);
        unsigned parse_data(double &data, string data_filename);
        unsigned parse_data(vector<unsigned> &data, string data_filename);
        unsigned parse_data(vector<double> &data, string data_filename);
        unsigned parse_data(vector<string> &data, string data_filename);

        unsigned write_data(bool &data, string data_filename);
        unsigned write_data(int &data, string data_filename);
        unsigned write_data(unsigned &data, string data_filename);
        unsigned write_data(double &data, string data_filename);
        unsigned write_data(string &data, string data_filename);
        unsigned write_data(vector<bool> &data, string data_filename);
        unsigned write_data(vector<unsigned> &data, string data_filename);
        unsigned write_data(vector<double> &data, string data_filename);
        unsigned write_data(vector<string> &data, string data_filename);

        void set_folder(string folder);

        string get_filepath(string &value_filename);
    };
}
#endif
