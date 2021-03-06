#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <string>
#include <list>
#include <vector>
#include "structs.h"


TSolution* VNS(std::string gen_mode, double gen_param, int Kmayus, int kmax, int max_time, std::string next_opt, double *next_opt_param, 
              std::string shake_opt, std::string ls1, double ls1_param, std::string ls2, double ls2_param, std::string accept_mode, double accept_param, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol);
//TSolution VNS(int max_iterations, int max_time, Parameters &param);
TSolution* SA(std::string gen_mode, double gen_param, int kmax, int max_time, std::string next_opt, double *next_opt_param, 
              std::string shake_opt, std::string cooling_opt, double cooling_param, double t0, std::string ls1, double ls1_param, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol);
TSolution* ILS(std::string gen_mode, double gen_param, int kmax, int max_time, 
              std::string shake_opt, std::string ls1, double ls1_param, int n_perturbations, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol);
TSolution* TS(int max_iterations, std::string gen_mode, double gen_param, std::string ls, double ls_param);
//TSolution TS(int max_iterations, int max_time, Parameters &param);
TSolution* MM();
void GA(int max_iterations, int max_time, std::vector<TSolution*> &pop, int pop_size, int lambda, std::string gen_mode, double gen_param, 
    std::string sel_mode, std::string cross_mode, std::string mut_mode, double mut_prob, std::string repl_mode, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol);
void PSO(int max_iterations, int max_time, std::vector<TSolution*> &pop, int pop_size, std::string gen_mode, double gen_param, 
         float omega, float phip, float phig, int lb, int ub, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol);


#endif
