#ifndef GENERATE_H
#define GENERATE_H

#include <string>
#include "structs.h"


TSolution* initial_solution(std::string mode, double param);
TSolution* new_rand_solution();
TSolution* START_algorithm(double **delta, double theta = 0.25);

#endif