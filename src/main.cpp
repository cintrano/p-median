/******************************************************************************
 * Main File
 * Author: Christian Cintrano
 * Date: 2018-05-27
 * Version: 2.0
 * Main file to solve the p-median problem
 * ****************************************************************************/

#include <iostream>
#include <iomanip>
#include <array>
#include <fstream>      /* reading a text file */
#include <string>
#include <cstring>
#include <chrono>       /* measure time */
#include <unistd.h>     /* sleep */
#include <algorithm>    /* find, sort */
#include <set>
#include <map>
#include <vector>
#include "utilities.h"
#include "structs.h"
#include "localsearches.h"
#include "algorithms.h"

using namespace std;

// Program params. Default values
std::map<std::string, int> params_int;
std::map<std::string, bool> params_bool;
std::map<std::string, double> params_double;
std::map<std::string, unsigned> params_unsigned;
std::map<std::string, string> params_string;

int num_inputs = 0;
string PATH = "";
string INTERCHANGE = "SHAKING";
string NEIGHBORHOOD_MODE = "NEAR";
string ALGO;
string LOCAL_SEARCH_1;
string LOCAL_SEARCH_2;
string GENERATE_NEW_SOLUTION;
bool WRITE_FILES = true;


// Functions
void delete_to_end();
void local_search(TSolution &x, unsigned mode);
TSolution initial_solution();
void default_params();
void run_algorithm(std::vector<TSolution*> &pop, TSolution* &best);
// Input
void read_conf(string filename);
void read_inputs(string filename);
void read_inputs(const std::string& f_customer, const std::string& f_facility);
void read_args(int size, char* args[]);
// Input of p-median problem
void read_customers(const std::string& f_customer);
void read_facilities(const std::string& f_customer);
void read_distance_matrix(const std::string& f_customer);
// Output
string get_file_name();
void write_results(string &filename, int time_init, int time_run, int iter, float best, float worst);
void write_population(string &filename, const std::vector<TSolution*> pop, int size, int ind_size);
void write_individual(string &filename, TSolution* ind, int ind_size);


int main(int argc, char* argv[]) {
    default_params();
    if (argc == 2)
    {
        read_conf(argv[1]);
    }
    else
    {
        read_args(argc, argv);
    }
    log("=== START ===");

    std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

    /* initialize random seed: */
    srand(params_int["seed"]);  // srand(time(NULL));

    log("SEED: " + std::to_string(params_int["seed"]));
    log("Inicializing distance matrix...", false);
    // The distance matrix must be computed
    if (num_inputs == 2)
    {    
        initialize_distance_matrix_two_points_array();
    }
    if (num_inputs == 1)
    {
        initialize_distance_matrix();
    }
    log("DONE");
    log("Inicializing neighborhood...", false);
    if (NEIGHBORHOOD_MODE == "NEAR")
    {
        log("mode near...", false);
        initialize_neighborhood();
    }
    if (NEIGHBORHOOD_MODE == "QUAD")
    {
        log("mode quad...", false);
        initialize_neighborhood_quadrants();
    }
    log("DONE");

    std::vector<TSolution*> pop(params_int["POP_SIZE"]);
    TSolution* best;

    std::chrono::steady_clock::time_point time_start= std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point time_end;

    // Algorithm
    log("Running algorithm...");
    run_algorithm(pop, best);

    time_end= std::chrono::steady_clock::now();

    
    log("Final population: ");
    if (DEBUG) print_population(pop, params_int["POP_SIZE"]);

    int time_init = std::chrono::duration_cast<std::chrono::milliseconds> (time_start - time_begin).count();
    int time_run = std::chrono::duration_cast<std::chrono::milliseconds> (time_end - time_start).count();

    if (WRITE_FILES)
    {
        string file_name = get_file_name();

        string opt_file_name = "opt";
        opt_file_name.append(file_name);
        string name = PATH + opt_file_name;
        write_results(name, time_init, time_run, 0, best->fitness, 0);

        string pop_file_name = "pop";
        pop_file_name.append(file_name);
        name = PATH + pop_file_name;
        write_population(name, pop, params_int["POP_SIZE"], P);

        string best_file_name = "best";
        best_file_name.append(file_name);
        name = PATH + best_file_name;
        write_individual(name, best, P);
    }


    std::cout << best->fitness;
    log("\n=== END ===");

    remove_structs();

    for (unsigned i = 0; i < pop.size(); ++i)
    {
        pop[i]->individual.clear();
        delete pop[i];
    }
    pop.clear();

    return 0;
}

void default_params()
{
    params_bool["DEBUG"] = false;

    params_int["seed"] = 13;
    params_int["MAX_TIME"] = std::numeric_limits<int>::max(); // Max time for execution in seconds
    params_int["POP_SIZE"] = 1; // Population size
    params_int["Gmax"] = 1; // MAX generations

    params_double["THETA"] =  0.25; // START_algorithm
    params_double["Laux"] =  10; 
    params_double["m"] =  0.2; 
    params_double["lambda"] =  1; 
    params_double["ALPHA"] =  10; // max{0.1p, 10} IS-RW 

    params_unsigned["kmax"]= 20;
    params_unsigned["Kmayus"] = 50;
}

void run_algorithm(std::vector<TSolution*> &pop, TSolution* &best)
{
    if (ALGO == "VNS") // same familly
    {
        double *next_opt_params = new double[2]{params_double["m"], params_double["lambda"]};
        best = VNS(params_string["init_sol"], params_double["init_sol"], params_int["Kmayus"], params_int["kmax"], params_int["MAX_TIME"], 
                   params_string["VNS_next_opt"], next_opt_params, params_string["shake"], params_string["ls1"], params_double["ls1"], 
                   params_string["ls2"], params_double["ls2"], params_string["accept"], params_double["accept"]);
        pop[0] = best;
    }
    if (ALGO == "SA") // same familly
    {
        double *next_opt_params = new double[2]{params_double["m"], params_double["lambda"]};
        best = SA(params_string["init_sol"], params_double["init_sol"], params_int["Kmayus"], params_int["kmax"], params_int["MAX_TIME"], 
                   params_string["VNS_next_opt"], next_opt_params, params_string["shake"], params_string["ls1"], params_double["ls1"], params_double["accept"]);
        pop[0] = best;
    }
    /*
    if (ALGO == "TS") // same familly
    {
        best = TS(params_int["Gmax"], params_string["init_sol"], params_double["init_sol"], params_string["ls1"], params_double["ls1"]);
        //pop = new TSolution[1];
        pop[0] = best;
    }
    */
    if (ALGO == "GA") // same familly
    {
        GA(params_int["Gmax"], params_int["MAX_TIME"], pop, params_int["POP_SIZE"], params_double["lambda"], params_string["init_sol"], params_double["init_sol"], 
            params_string["sel_mode"], params_string["cross_mode"], params_string["mut_mode"], params_double["mut_prob"], params_string["repl_mode"]);
        best = find_best(pop, params_int["POP_SIZE"]);
    }
}

/******************************************************************************
 * INPUT
 * ****************************************************************************/

vector<string> split(const string& str, int delimiter(int) = ::isspace)
{
  vector<string> result;
  auto e=str.end();
  auto i=str.begin();
  while(i!=e)
  {
    i = find_if_not(i,e, delimiter);
    if(i==e) break;
    auto j = find_if(i,e, delimiter);
    result.push_back(string(i,j));
    i=j;
  }
  return result;
}

void read_conf(string filename)
{
    string line;
    ifstream myfile(filename);
    if (myfile.is_open())
    {
        int i = 0;
        while ( getline(myfile,line) )
        {
            // line
            vector<string> result = split(line);
            string key = result[0];
            if (key == "P")
            {
                P = strtof((result[1]).c_str(),0);
            }
            if (key == "N")
            {
                N = strtof((result[1]).c_str(),0);
            }
            if (key == "F")
            {
                F = strtof((result[1]).c_str(),0);
            }
            if (key == "input")
            {
                read_inputs(result[1]);
            }
            if (key == "POP_SIZE")
            {
                params_int["POP_SIZE"] = strtof((result[1]).c_str(),0);
            }
            if (key == "ITER_MAX")
            {
                params_int["Gmax"] = strtof((result[1]).c_str(),0);
            }
            if (key == "K")
            {
                k = strtof((result[1]).c_str(),0);
            }
            i++;
        }
        myfile.close();
    }
    else cout << "Unable to open file";
}

void read_inputs(string filename)
{
    points = new double*[N]; // dynamic array (size 10) of pointers to int
    for (int i = 0; i < N; ++i) {
      points[i] = new double[3];
      // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
    }

    string line;
    ifstream myfile(filename);
    if (myfile.is_open())
    {
        int i = 0;
        while ( getline(myfile,line) )
        {
            // line
            vector<string> result = split(line);
            for (int j = 0; j < 3; ++j)
            {
                points[i][j] = strtof((result[j]).c_str(),0); 
            }
            i++;
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
}

void read_customers(const std::string& f_customer)
{
    log("Reading: " + f_customer);
    //points = new float[N][3];
    customer_points = new double*[N]; // dynamic array (size 10) of pointers to int
    for (int i = 0; i < N; ++i) {
      customer_points[i] = new double[3];
      // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
    }
    string line;
    ifstream myfile(f_customer);
    if (myfile.is_open())
    {
        int i = 0;
        while ( getline(myfile,line) )
        {
            // line
            vector<string> result = split(line);
            for (int j = 1; j < 3; ++j)
            {
                customer_points[i][j] = strtof((result[j-1]).c_str(),0); 
            }
            customer_points[i][0] = strtof((result[2]).c_str(),0); // weights
            i++;
        }
        myfile.close();
    }
    else std::cout << "Unable to open file\n";
}

void read_distance_matrix(const std::string& f_matrix)
{
    log("Reading: " + f_matrix);
    if (N * F < MAX_DISTANCE_SIZE)
    {
        double* temp = new double[N * F];
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            D[i] = (temp + i * F);
        }
    }
    else
    {
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            D[i] = new double[F];
        }
    }

    string line;
    ifstream myfile(f_matrix);
    if (myfile.is_open())
    {
        for (int i = 0; i < N; ++i)
        {
            getline(myfile,line);
            // line
            vector<string> result = split(line);
            for (int j = 0; j < F; ++j)
            {
                D[i][j] = strtof((result[j]).c_str(),0);
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file\n";
}

void read_facilities(const std::string& f_facility)
{
    log("Reading: " + f_facility);
    facility_points = new double*[F]; // dynamic array (size 10) of pointers to int
    for (int i = 0; i < F; ++i) {
      facility_points[i] = new double[3];
      // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
    }
    string line;
    ifstream myfile2(f_facility);
    if (myfile2.is_open())
    {
        int i = 0;
        while ( getline(myfile2, line) )
        {
            // line
            vector<string> result = split(line);
            for (int j = 0; j < 3; ++j)
            {
                facility_points[i][j] = strtof((result[j]).c_str(),0); 
            }
            i++;
        }
        myfile2.close();
    }
    else std::cout << "Unable to open file\n";
}

void read_inputs(const std::string& f_customer, const std::string& f_facility)
{
    read_customers(f_customer);
    read_facilities(f_facility);    
}

void read_inputs(const std::string& f_matrix, const std::string& f_customer, const std::string& f_facility)
{
    read_customers(f_customer);
    read_facilities(f_facility);  
    read_distance_matrix(f_matrix);
    for (int j = 0; j < P; ++j)
    {
        individual_index.push_back(j);
    }
}

void read_args(int size, char* args[])
{
    P = std::stoi(args[1]);
    N = std::stoi(args[2]);
    F = std::stoi(args[3]);
    int i = 4;
    try
    {
        while(i < size) {
            if (std::strcmp(args[i], "--input") == 0)
            {
                i++;
                num_inputs = std::stoi(args[i]);
                switch(num_inputs) {
                    case 1 : read_inputs(args[i]);
                             break;
                    case 2 : 
                        {
                            i++;
                            string f_customer(args[i]);
                            i++;
                            string f_facility(args[i]);
                            read_inputs(f_customer, f_facility);
                            break;
                        }
                    case 3 : 
                        {
                            i++;
                            string f_matrix(args[i]);
                            i++;
                            string f_customers(args[i]);
                            i++;
                            string f_facilities(args[i]);
                            read_inputs(f_matrix, f_customers, f_facilities);
                            break;
                        }
                }
            }
            if (std::strcmp(args[i], "--output") == 0)
            {
                i++;
                PATH = args[i];
            }
            if (std::strcmp(args[i], "--no-output") == 0)
            {
                WRITE_FILES = false;
            }
            if (std::strcmp(args[i], "--pop") == 0)
            {
                i++;
                params_int["POP_SIZE"] = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--timer") == 0)
            {
                i++;
                params_int["MAX_TIME"] = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--iter") == 0)
            {
                i++;
                if (std::strcmp(args[i], "Np/5") == 0)
                {
                    params_int["Gmax"] = N * P / 5;
                }
                else if (std::strcmp(args[i], "100p") == 0)
                {
                    params_int["Gmax"] = 100 * P;
                }
                else if (std::strcmp(args[i], "MAX2N-100") == 0) // Default TS in Rolland1997
                {
                    params_int["Gmax"] = std::max(2* N, 100);
                }
                else
                {
                    params_int["Gmax"] = std::stoi(args[i]);
                }
            }
            if (std::strcmp(args[i], "--neighborhood") == 0)
            {
                i++;
                NEIGHBORHOOD_MODE = args[i];
            }
            if (std::strcmp(args[i], "--k") == 0)
            {
                i++;
                k = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--kmax") == 0)
            {
                i++;
                params_int["kmax"] = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--Kmayus") == 0)
            {
                i++;
                params_int["Kmayus"] = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--lambda") == 0)
            {
                i++;
                params_double["lambda"] = std::stod(args[i]);
            }
            if (std::strcmp(args[i], "--m") == 0)
            {
                i++;
                params_double["m"] = std::stod(args[i]);
            }
            if (std::strcmp(args[i], "--accept") == 0)
            {
                i++;
                params_string["accept"] = args[i];
                if (params_string["accept"] == "SA")
                {
                    i++;
                    params_double["accept"] = std::stod(args[i]);
                }
            }
            if (std::strcmp(args[i], "--next") == 0)
            {
                i++;
                params_string["VNS_next_opt"] = args[i];
            }
            if (std::strcmp(args[i], "--debug") == 0)
            {
                DEBUG = true;
            }
            if (std::strcmp(args[i], "--seed") == 0)
            {
                i++;
                params_int["seed"] = std::stoi(args[i]);
            }
            if (std::strcmp(args[i], "--algo") == 0)
            {
                i++;
                ALGO = args[i];
            }
            if (std::strcmp(args[i], "--gen") == 0) // Generate new solution mode
            {
                i++;
                params_string["init_sol"] = args[i];

                if (params_string["init_sol"] == "START")
                {
                    i++;
                    params_double["init_sol"] = std::stod(args[i]);
                }
            }
            if (std::strcmp(args[i], "--change") == 0)
            {
                i++;
                INTERCHANGE = args[i];
            }
            if (std::strcmp(args[i], "--shake") == 0)
            {
                i++;
                params_string["shake"] = args[i];
            }
            if (std::strcmp(args[i], "--ls1") == 0)
            {
                i++;
                params_string["ls1"] = args[i];

                if (params_string["ls1"] == "IALT" || params_string["ls1"] == "ISRW" || params_string["ls1"] == "IMP")
                {
                    i++;
                    params_double["ls1"] = std::stod(args[i]);
                }
            }
            if (std::strcmp(args[i], "--ls2") == 0)
            {
                i++;
                params_string["ls2"] = args[i];

                if (params_string["ls2"] == "IALT" || params_string["ls2"] == "ISRW"|| params_string["ls2"] == "IMP")
                {
                    i++;
                    params_double["ls2"] = std::stod(args[i]);
                }
            }
            if (std::strcmp(args[i], "--ga-sel") == 0)
            {
                i++;
                params_string["sel_mode"] = args[i];
            }
            if (std::strcmp(args[i], "--ga-cross") == 0)
            {
                i++;
                params_string["cross_mode"] = args[i];
            }
            if (std::strcmp(args[i], "--ga-mut") == 0)
            {
                i++;
                params_string["mut_mode"] = args[i];
                i++;
                params_double["mut_prob"] = std::stod(args[i]);
            }
            if (std::strcmp(args[i], "--ga-repl") == 0)
            {
                i++;
                params_string["repl_mode"] = args[i];
            }
            i++;
        }
    }
    catch (std::exception const &e)
    {
        // This could not be parsed into a number so an exception is thrown.
        // atoi() would return 0, which is less helpful if it could be a valid value.
        cerr << "Incorrect format.\n";
    }
}

/******************************************************************************
 * OUTPUT
 * ****************************************************************************/

string get_file_name()
{
    string file_name = "_";
    file_name.append(ALGO).append("_")
             .append(to_string(N)).append("_")
             .append(to_string(F)).append("_")
             .append(to_string(P)).append("_");
    file_name.append("S-" + to_string(params_int["seed"]) + "_");
    // GENERATE_NEW_SOLUTION
    if (GENERATE_NEW_SOLUTION == "RAND" || GENERATE_NEW_SOLUTION == "100IMP" || GENERATE_NEW_SOLUTION == "100RAND") file_name.append(GENERATE_NEW_SOLUTION);
    if (GENERATE_NEW_SOLUTION == "START") file_name.append(GENERATE_NEW_SOLUTION).append("-").append(to_string(params_double["THETA"]));
    file_name.append("_");
    // LS 1
    if (LOCAL_SEARCH_1 == "IALT") file_name.append(LOCAL_SEARCH_1).append("-").append(to_string(params_double["Laux"]));
    if (LOCAL_SEARCH_1 == "ISRW") file_name.append(LOCAL_SEARCH_1).append("-").append(to_string(params_double["ALPHA"]));
    if (LOCAL_SEARCH_1 == "IMP" || LOCAL_SEARCH_1 == "FI") file_name.append(LOCAL_SEARCH_1);
    file_name.append("_");
    // LS 2
    if (LOCAL_SEARCH_2 == "IALT") file_name.append(LOCAL_SEARCH_2).append("-").append(to_string(params_double["Laux"]));
    if (LOCAL_SEARCH_2 == "ISRW") file_name.append(LOCAL_SEARCH_2).append("-").append(to_string(params_double["ALPHA"]));
    if (LOCAL_SEARCH_2 == "IMP" || LOCAL_SEARCH_2 == "FI") file_name.append(LOCAL_SEARCH_2);
    file_name.append(".ssv");
    return file_name;
}

void write_results(string &filename, int time_init, int time_run, int iter, float best, float worst) {
    ofstream myfile (filename);
    myfile.precision(std::numeric_limits<double>::digits10);
    myfile.precision(std::numeric_limits<float>::digits10);
    if (myfile.is_open())
    {
        myfile << time_init << " " << time_run << " " << iter << " " << best << " " << worst;
        myfile.close();
    }
    else std::cout << "Unable to open file";
}

void write_population(string &filename, const std::vector<TSolution*> pop, int size, int ind_size) {
    ofstream myfile (filename);
    myfile.precision(std::numeric_limits<double>::digits10);
    myfile.precision(std::numeric_limits<float>::digits10);
    if (myfile.is_open())
    {
        for (int i = 0; i < size; ++i)
        {
            myfile << pop[i]->fitness << " ";
            for (int j = 0; j < ind_size; ++j)
            {
                myfile << pop[i]->individual[j] << " ";
            }
            myfile << '\n';
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
}

void write_individual(string &filename, TSolution* ind, int ind_size) {
    ofstream myfile (filename);
    myfile.precision(std::numeric_limits<double>::digits10);
    myfile.precision(std::numeric_limits<float>::digits10);
    if (myfile.is_open())
    {
        myfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << ind->fitness << " ";
        for (int j = 0; j < ind_size; ++j)
        {
            myfile << ind->individual[j] << " ";
        }
        myfile << '\n';
        myfile.close();
    }
    else std::cout << "Unable to open file";
}
