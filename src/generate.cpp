#include <limits> 
#include <string>
#include "generate.h"
#include "structs.h"
#include "utilities.h"
#include "localsearches.h"

TSolution* initial_solution(std::string mode, double param)
{
    log("-- initialize sol..." + mode + "...", false);
    TSolution* x;
    if (mode == "RAND")
    {
        x = new_rand_solution();
    }
    else if (mode == "START")
    {
        double **delta;
        delta = new double*[F];
        for (int delt_index = 0; delt_index < F; ++delt_index)
        {
            delta[delt_index] = new double[F];
        }
        x = START_algorithm(delta, param); // params_double["THETA"]
        for (int i = 0; i < F; ++i)
        {
            delete delta[i];
        }
        delete delta;
    }
    else if (mode == "100RAND" || mode == "100IMP")
    {
        int population_size = 100;
        std::vector<TSolution*> pop(population_size);
        float f_min, f_max, f_current;
        pop[0] = new_rand_solution();
        if (mode == "100IMP")
        {
            IMP(pop[0], param);
            //log("IMP ", false);
            //if (DEBUG) print_solution(pop[0]);
        }
        f_min = pop[0]->fitness;
        f_max = pop[0]->fitness;
        int best = 0;
        for (int i = 1; i < population_size; ++i)
        {
            pop[i] = new_rand_solution();
            if (mode == "100IMP")
            {
                IMP(pop[i], param);
                //log("IMP ", false);
                //if (DEBUG) print_solution(pop[i]);
            }
            f_current = pop[i]->fitness;
            if (f_current > f_max) {
                f_max = f_current;
            } else if (f_current < f_min)
            {
                f_min = f_current;
                best = i;
            }
        }
        x = new TSolution;
        for (int i = 0; i < P; ++i)
        {
            x->individual.push_back(pop[best]->individual[i]);
        }
        x->fitness = pop[best]->fitness;
        //log("--->");
        //if (DEBUG) print_solution(x);

        for (int i = 0; i < population_size; ++i)
        {
            delete pop[i];
        }
        pop.clear();
    }
    else { // DEFAULT is RAND
        x = new_rand_solution();
    }
    if (DEBUG) print_solution(x);
    log("END");
    return x;
}

TSolution* new_rand_solution()
{
    TSolution* sol = new TSolution;
    //generate random numbers:
    for (int i=0;i<P;i++) {
        bool check; //variable to check or number is already used
        int n = -1; //variable to store the number in
        do {
            n=rand()%F;
            //check or number is already used:
            check=true;
            for (int j=0;j<i;j++)
                if (n == sol->individual[j]) { //if number is already used
                    check=false; //set check to false
                    break; //no need to check the other elements of value[]
                }
        } while (!check); //loop until new, unique number is found
        sol->individual.push_back(n); //store the generated number in the array
    }
    sol->fitness = fitness(sol);
    return sol;
}

TSolution* START_algorithm(double **delta, double theta)
{
    log("Generate by START",false);
    std::vector<double> v(F);
    for (int i = 0; i < F; ++i) // TODO Cahange to weights values
    {
        v[i] = 1.0;
    }
    int removed = 0;
    double u;
    // Step 0
    for (int i = 0; i < F; ++i)
    {
        delta[i][i] = 1;
        for (int j = i+1; j < F; ++j)
        {
            u = (double) rand() / RAND_MAX;
            double dist = distanceCalculate(facility_points[i][1], facility_points[i][2], facility_points[j][1], facility_points[j][2]);
            delta[i][j] = (v[i]*v[j]) * (dist) * (theta * u) / (v[i]+v[j]);
        }
    }

    while(F - removed > P)
    {
        log(".",false);
        // 1. Find the pair i<j for which Delta_{ij} is minimized.
        int imin = 0, jmin = 0;
        double delta_min = std::numeric_limits<double>::infinity();
        double current;
        for (int i = 0; i < F; ++i)
        {
            if (delta[i][i])
            {
                for (int j = i+1; j < F; ++j)
                {
                    if (delta[j][j])
                    {
                        current = delta[i][j];
                        if (current < delta_min)
                        {
                            delta_min = current;
                            imin = i;
                            jmin = j;
                        }
                    }
                }
            }
        }
        // 2. Create a location for a new facility at (ciXi+vjXj)/vi+vj with a weight vi+vj.
        // TODO falta la relocalizacion
        double weight = v[imin] + v[jmin];
        // 3.1 Remove current facilities i and j, and label the index of the new facility and its weight with i. 
        delta[jmin][jmin] = 0;
        removed++;
        v[imin] = weight;
        // 3.2 Calculate for all r != i: Delta_{ri} for r<i and Delta_{ir} for r>i by (2). The number of facilities is reduced by one.
        // Column
        for (int r = 0; r < imin; ++r)
        {
            u = (double) rand() / RAND_MAX;
            double dist = distanceCalculate(facility_points[r][1], facility_points[r][2], facility_points[imin][1], facility_points[imin][2]);
            delta[r][imin] = (v[r]*v[imin]) * (dist) * (theta * u) / (v[r]+v[imin]);
        }
        // Row
        for (int r = imin + 1; r < F; ++r)
        {
            u = (double) rand() / RAND_MAX;
            double dist = distanceCalculate(facility_points[imin][1], facility_points[imin][2], facility_points[r][1], facility_points[r][2]);
            delta[imin][r] = (v[imin]*v[r]) * (dist) * (theta * u) / (v[imin]+v[r]);
        }
    }

    TSolution *sol = new TSolution;
    for (int i = 0; i < F; ++i)
    {
        if (delta[i][i])
        {
            sol->individual.push_back(i);
        }
    }
    sol->fitness = fitness(sol);
    return sol;
}
