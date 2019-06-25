#include <iostream>
#include <math.h>       /* sin, cos, M_PI, pow */
#include <chrono>       /* measure time */
#include <random>
#include <algorithm>
#include <set>
#include <map>
#include "algorithms.h"
#include "structs.h"
#include "utilities.h"
#include "localsearches.h"
#include "generate.h"


/******************************************************************************
 * Variable Neighborhood Search
 * Author: Christian Cintrano
 * Date: 2018-04-16
 * Updated: 2018-04-16
 * Implementation of the VNS to solve the p-median problem
 *
 * Original paper:
 * Mladenović, N., & Hansen, P. (1997). 
 * Variable neighborhood search. 
 * Computers & Operations Research, 24(11), 1097–1100. 
 * http://doi.org/10.1016/S0305-0548(97)00031-2
 * ****************************************************************************/
int next_k(std::string mode, int index, double *params, int kmax);
double psi(double x, double m, double lambda); // VNS: DVNS next_param
void shake(std::string mode, TSolution* original, TSolution* &x, int i); // VNS: shakes
bool acceptation(std::string mode, double param, TSolution* &x_old, TSolution *x_new);
TSolution* shake_rand(TSolution* x, int i);
TSolution* shake_rand_neighborhood(TSolution* x, int i);

TSolution* VNS(std::string gen_mode, double gen_param, int Kmayus, int kmax, int max_time, std::string next_opt, double *next_opt_param,
              std::string shake_opt, std::string ls1, double ls1_param, std::string ls2, double ls2_param, std::string accept_mode, double accept_param, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol)
{
    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;

    TSolution* xprime;
    xprime = initial_solution(gen_mode, gen_param);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

    t_start = std::chrono::steady_clock::now();
    t_current = std::chrono::steady_clock::now();
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    local_search(xprime, ls1, ls1_param, max_time - current_time);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

    int counter = 0;
    log("---- G " + std::to_string(counter) + " Max time: " + std::to_string(max_time));

    bool restart = true; 

    int kindex;

    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    // RUN

    TSolution* x;
    int index, j;
    log(std::to_string(current_time));
    while (restart && (current_time < max_time)) // General loop
    { 
        restart = false;
        if (DEBUG)
        {
            counter++;
            if (counter % 50 == 0)
            {
                log("---- G " + std::to_string(counter) + " " + std::to_string(xprime->fitness));
            }
        }

        j = 1;
        while(!restart && j <= Kmayus) // Number of no new optimal found
        {
            index = 1;
            while(!restart && index <= kmax) { // kindex neighbourhood

                kindex = next_k(next_opt, index, next_opt_param, kmax);
                shake(shake_opt, xprime, x, kindex);
                t_current= std::chrono::steady_clock::now();
                current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
                local_search(x, ls2, ls2_param, max_time - current_time);
                if (DEBUG) print_solution(x);
                restart = acceptation(accept_mode, accept_param, xprime, x); // true if is new optimal
                if (SAVING)
                {
                    std::cerr << "\t F " << xprime->fitness << "\n";
                    evo_fitness.push_back(xprime->fitness);
                    std::vector<int> v(xprime->individual);
                    evo_sol.push_back(v);
                }
                index++;
            }
            j++;
        }
        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    }

    return xprime;
}

int next_k(std::string mode, int index, double *params, int kmax)
{
    int kindex;
    if (mode == "DVNS")
    {
        double u = (double) rand() / RAND_MAX;
        kindex = round(psi(u, params[0], params[1]) * kmax) + 1; //next_opt_param[0] = m, next_opt_param[1] = lambda
    }
    else
    {
        kindex = index;
    }
    return kindex;
}

double psi(double x, double m, double lambda)
{
    // m and lambda are params of this function
    double division = ((3 * lambda) / (lambda - 1)) * m * m;
    double term1 = x / (1 - 3 * m + division) ;
    double term2 = (x * x) - 3 * m * x + division;
    return term1 * term2;
}

void shake(std::string mode, TSolution* original, TSolution* &x, int i)
{
    if (mode == "SHAKING")
    {
        x = shake_rand(original, i);
    }
    else if (mode == "NRAND")
    {
        x = shake_rand_neighborhood(original, i);
    }
    else
    {
        x = new TSolution;
        for (int j = 0; j < P; ++j)
        {
            x->individual.push_back(original->individual[j]);
        }
        x->fitness = original->fitness;
        log("Shake NONDE ");
        if (DEBUG) print_solution(x);
    }
}

bool acceptation(std::string mode, double param, TSolution* &x_old, TSolution *x_new)
{
    if (mode == "ELITIST")
    {
        if (x_new->fitness < x_old->fitness)
        {
            x_old->individual.clear();
            delete x_old;
            x_old = x_new;
            return true;
        }
        else
        {
            x_new->individual.clear();
            delete x_new;
            return false;
        }
    }
    if (mode == "WALK")
    {
        x_old->individual.clear();
        delete x_old;
        x_old = x_new;
        return true;
    }
    if (mode == "SA")
    {
        if (x_new->fitness < x_old->fitness)
        {
            x_old->individual.clear();
            delete x_old;
            x_old = x_new;
            return true;
        }
        else
        {
            if (rand()/RAND_MAX <= param)
            {
                x_old->individual.clear();
                delete x_old;
                x_old = x_new;
                return true;
            }
            else
            {
                x_new->individual.clear();
                delete x_new;
                return false;   
            }
        }
    }
    return false;
}

TSolution* shake_rand(TSolution* x, int i)
{
    TSolution* out = new TSolution;
    for (int j = 0; j < P; ++j)
    {
        out->individual.push_back(x->individual[j]);
    }
    
    std::set<int> checked;
    while(checked.size() < (unsigned) (i * 2))
    {
        int v = rand() % P;
        int dest = rand() % F;
        while(checked.find(dest) != checked.end() || ind_contains(out, dest))
        {
            dest = rand() % F;
        }
        checked.insert(out->individual[v]);
        checked.insert(dest);
        out->individual[v] = dest;
    }

    out->fitness = fitness(out);
    checked.clear();
    return out;  
}

TSolution* shake_rand_neighborhood(TSolution* x, int i)
{
    log("\tshake_rand_neighborhood...");
    TSolution* out = new TSolution;
    for (int j = 0; j < P; ++j)
    {
        out->individual.push_back(x->individual[j]);
    }
    log("");

    std::vector<int> neighbors_indexes(k);
    for (int j = 0; j < k; ++j)
    {
        neighbors_indexes[j] = j;
    }

    std::random_shuffle(individual_index.begin(), individual_index.end());
    std::set<int> checked;
    int position = 0;
    int v, neighbors_index, new_point;
    bool checked_previously;
    while(position < (i-1) && position < P)
    {
        v = out->individual[individual_index[position]]; // size of the neighborhood between facilities
        checked.insert(v);
        std::random_shuffle(neighbors_indexes.begin(), neighbors_indexes.end());

        neighbors_index = 0;
        new_point = Neighborhood[v][neighbors_indexes[neighbors_index]];
        checked_previously = checked.find(new_point) != checked.end();
        checked.insert(new_point);
        neighbors_index++;

        while((neighbors_index < k && ind_contains(out, new_point)) || (neighbors_index < k && checked_previously))
        {
            //std::cout << "NP "<< new_point;
            new_point = Neighborhood[v][neighbors_indexes[neighbors_index]];
            //std::cout << " " << new_point << "  ___ ";
            checked_previously = checked.find(new_point) != checked.end();
            checked.insert(new_point);
            neighbors_index++;
        }

        if (neighbors_index < k)
        {
            //std::cout << "]   ["<< new_point;
            out->individual[individual_index[position]] = new_point;
        }
        position++;
    }
    out->fitness = fitness(out);
    neighbors_indexes.clear();
    checked.clear();
    log("\tend");
    return out;   
}

/******************************************************************************
 * Basic Simulated Annealing Algorithm
 * Author: Christian Cintrano
 * Date: 2018-09-20
 * Updated: 2018-09-20
 * Version: 1.0
 * Implementation of a SA to solve the p-median problem
 * ****************************************************************************/
double cooling(std::string mode, double param, double t0, double t);
void potential(TSolution* &x_new, TSolution* &x_old, double t0, double t);

TSolution* SA(std::string gen_mode, double gen_param, int kmax, int max_time, std::string next_opt, double *next_opt_param,
              std::string shake_opt, std::string cooling_opt, double cooling_param, double t0, std::string ls1, double ls1_param, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol)
{
    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;
    t_start= std::chrono::steady_clock::now();   

    TSolution* xprime; 
    xprime = initial_solution(gen_mode, gen_param);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

                t_current= std::chrono::steady_clock::now();
                current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    local_search(xprime, ls1, ls1_param, max_time - current_time);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

    int counter = 0;
    log("---- G " + std::to_string(counter) + " Max time: " + std::to_string(max_time));

    int kindex;

    //t_start= std::chrono::steady_clock::now();
    t_start= std::chrono::steady_clock::now();   
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    // RUN

    TSolution* x;
    log(std::to_string(current_time));
    int index = 1;
    double t;
    while (index <= kmax && (current_time < max_time)) // General loop
    { 
        if (DEBUG)
        {
            counter++;
            if (counter % 50 == 0)
            {
                log("---- G " + std::to_string(counter) + " " + std::to_string(xprime->fitness));
            }
        }
        log("--- cooling");
        t = cooling(cooling_opt, cooling_param, t0, index);
        log("--- next");
        kindex = next_k(next_opt, index, next_opt_param, kmax);
        log("--- perturbation");
        shake(shake_opt, xprime, x, kindex); // xprime = original, x = new
        //local_search(x, ls2, ls2_param);
        if (DEBUG) print_solution(xprime);
        if (DEBUG) print_solution(x);
        log("--- potential");
        potential(x, xprime, t0, t); // x = new, xprime = original
        
        if (SAVING)
        {
            evo_fitness.push_back(xprime->fitness);
            std::vector<int> v(xprime->individual);
            evo_sol.push_back(v);
        }
        log("--- .");
        index++;

        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    }

    return xprime;
}

double cooling(std::string mode, double param, double t0, double t)
{
    if (mode == "EXPONENTIAL")
    {
        double Tc = t0 * pow(param, t);
        return Tc;
    }
    else if (mode == "LINEAL")
    {
        double Tc = t0 - (param * t);
        return Tc;
    }
    else
    {
        return t;
    }
}

void potential(TSolution* &x_new, TSolution* &x_old, double t0, double t)
{
    if (x_new->fitness < x_old->fitness)
    {
        x_old->individual.clear();
        delete x_old;
        x_old = x_new;
    }
    else
    {
        double tnew = t0 * ((double) rand() / (RAND_MAX));
        log(std::to_string(t0) + " " + std::to_string(tnew) + " " + std::to_string(t));
        if (tnew <= t)
        {
            x_old->individual.clear();
            delete x_old;
            x_old = x_new;
        }
        else
        {
            x_new->individual.clear();
            delete x_new; 
        }
    }
}

/******************************************************************************
 * Basic iterated Local Search
 * Author: Christian Cintrano
 * Date: 2018-09-25
 * Updated: 2018-09-25
 * Version: 1.0
 * Implementation of a ILS to solve the p-median problem
 * ****************************************************************************/
TSolution* ILS(std::string gen_mode, double gen_param, int kmax, int max_time, 
              std::string shake_opt, std::string ls1, double ls1_param, int n_perturbations, std::list<double> &evo_fitness, std::list<std::vector<int>> &evo_sol)
{
    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;

    TSolution* xprime;
    xprime = initial_solution(gen_mode, gen_param);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    local_search(xprime, ls1, ls1_param, max_time - current_time);

    log("Initial individual: ", false);
    if (DEBUG) print_solution(xprime);

    int counter = 0;
    log("---- G " + std::to_string(counter) + " Max time: " + std::to_string(max_time));


    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    // RUN

    TSolution* x;
    log(std::to_string(current_time));
    int index = 1;
    while (index <= kmax && (current_time < max_time)) // General loop
    { 
        if (DEBUG)
        {
            counter++;
            if (counter % 50 == 0)
            {
                log("---- G " + std::to_string(counter) + " " + std::to_string(xprime->fitness));
            }
        }
        
        shake(shake_opt, xprime, x, n_perturbations);

        t_current= std::chrono::steady_clock::now();
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
        local_search(x, ls1, ls1_param, max_time - current_time);
        acceptation("ELITIST", 0, xprime, x); // true if is new optimal

        if (SAVING)
        {
            evo_fitness.push_back(xprime->fitness);
            std::vector<int> v(xprime->individual);
            evo_sol.push_back(v);
        }
        index++;

        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    }

    return xprime;
}

/******************************************************************************
 * Basic Genetic Algorithm
 * Author: Christian Cintrano
 * Date: 2018-03-16
 * Updated: 2018-04-10
 * Version: 1.6
 * Implementation of a GA to solve the p-median problem
 * ****************************************************************************/
void generate_population(std::vector<TSolution*> &pop, int pop_size, std::string gen_mode, double gen_param, TSolution* &best, TSolution* &worst);
void selection(std::string mode, std::vector<TSolution*> &pop, int pop_size, TSolution* &parent1, TSolution* &parent2);
TSolution* crossover(std::string mode, TSolution* parent1, TSolution* parent2);
void mutation(std::string mode, double prob, TSolution* &ind);
int replacement(std::string mode, TSolution* offspring, std::vector<TSolution*> &pop, int pop_size, TSolution* &best, int worst_idx);
void replacement(std::string mode, std::vector<TSolution*> &childrens, int lambda, std::vector<TSolution*> &pop, int mu);
// Crossover options
TSolution* merging(TSolution* p1, TSolution* p2);
TSolution* cross_1point(TSolution* p1, TSolution* p2);
TSolution* cupcap(TSolution* p1, TSolution* p2);


void GA(int max_iterations, int max_time, std::vector<TSolution*> &pop, int pop_size, int lambda, std::string gen_mode, double gen_param, 
    std::string sel_mode, std::string cross_mode, std::string mut_mode, double mut_prob, std::string repl_mode)
{
    TSolution* best;
    TSolution* worst;
    generate_population(pop, pop_size, gen_mode, gen_param, best, worst);
    log("After population");
    int worst_idx = find_idx_of(worst, pop, pop_size);

    log("Initial population: ", false);
    if (DEBUG) print_population(pop, pop_size);

    int g = 0;
    int counter = 0;
    log("---- G " + std::to_string(g));
    TSolution* parent1;
    TSolution* parent2;
    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;
    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    int l;
    std::vector<TSolution*> childrens;
    TSolution* offspring;
    log("A");
    while(g < max_iterations && (current_time < max_time))
    {        
        log("B");
        g++;
        if (DEBUG)
        {
            counter++;
            if (counter % 50 == 0)
            {
                log("---- G " + std::to_string(counter));
            }
        }
        log("C");
        l = 0;
        while(l < lambda) {
        log("D");
            // Random Selection
            log("sel ");
            selection(sel_mode, pop, pop_size, parent1, parent2);
            if (DEBUG) print_solution(parent1);
            if (DEBUG) print_solution(parent2);
            // Crossover
            log("cross ");
            offspring = crossover(cross_mode, parent1, parent2);
            // Mutation
            log("mut ");
            mutation(mut_mode, mut_prob, offspring);
            log("fit ");
            if(DEBUG)
            {
                for (int i = 0; i < P; ++i)
                {
                    std::cout << offspring->individual[i] << " ";
                }
                std::cout << std::endl;
            }
            offspring->fitness = fitness(offspring);
            if (lambda == 1)
            {
            log("repl ");
                worst_idx = replacement(repl_mode, offspring, pop, pop_size, best, worst_idx);
            }
            else
            {
                log("#" + std::to_string(l));
                if (DEBUG) print_solution(offspring);
                childrens.push_back(offspring);
            }
            l++;
        }
        if (lambda != 1)
        {
            log("repl ");
            replacement(repl_mode, childrens, lambda, pop, pop_size);
        }
        childrens.clear();
        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();
    }
}

void generate_population(std::vector<TSolution*> &pop, int pop_size, std::string gen_mode, double gen_param, TSolution* &best, TSolution* &worst)
{
    double f_min, f_max, f_current;
    pop[0] = initial_solution(gen_mode, gen_param);
    if (DEBUG) print_solution(pop[0]);
    log("-- ind: 0\t" + std::to_string(pop[0]->fitness));
    best = pop[0];
    worst = pop[0];
    f_min = best->fitness;
    f_max = worst->fitness;
    for (int i = 1; i < pop_size; ++i)
    {
        pop[i] = initial_solution(gen_mode, gen_param);
        log("-- ind: " + std::to_string(i) + "\t" + std::to_string(pop[0]->fitness));
        f_current = pop[i]->fitness;
        if (f_current > f_max) {
            f_max = f_current;
            worst = pop[i];
        } else if (f_current < f_min)
        {
            f_min = f_current;
            best = pop[i];
        }
    }
}

void selection(std::string mode, std::vector<TSolution*> &pop, int pop_size, TSolution* &parent1, TSolution* &parent2)
{
    if (mode == "RAND")
    {
        int id1, id2;
        id1 = rand() % pop_size;
        id2 = rand() % pop_size;
        while(id1 == id2)
        {
            id2 = rand() % pop_size;
        }
        log("---->" + std::to_string(id1));
        if (DEBUG) print_solution(pop[id1]);
        log("---->" + std::to_string(id2));
        if (DEBUG) print_solution(pop[id2]);
        parent1 = pop[id1];
        parent2 = pop[id2];
    }
    else if (mode == "BETTERS")
    {
        int b1, b2;
        int f1, f2, fc;
        b1 = b2 = 0;
        f1 = pop[b1]->fitness;
        f2 = f1;
        for (int i = 0; i < pop_size; ++i)
        {
            fc = pop[i]->fitness;
            if (fc < f1)
            {
                b2 = b1;
                f2 = f1;
                b1 = i;
                f1 = fc;
            }
            else if (fc <= f2)
            {
                b2 = i;
                f2 = fc;
            }
        }
        parent1 = pop[b1];
        parent2 = pop[b2];
    }
    else // if (mode == "WORSE")
    {
        int w1, w2;
        int f1, f2, fc;
        w1 = w2 = 0;
        f1 = pop[w1]->fitness;
        f2 = f1;
        for (int i = 0; i < pop_size; ++i)
        {
            fc = pop[i]->fitness;
            if (fc > f1)
            {
                w2 = w1;
                f2 = f1;
                w1 = i;
                f1 = fc;
            }
            else if (fc >= f2)
            {
                w2 = i;
                f2 = fc;
            }
        }
        parent1 = pop[w1];
        parent2 = pop[w2];
    }
}

TSolution* crossover(std::string mode, TSolution* parent1, TSolution* parent2)
{
    if (DEBUG) print_solution(parent1);
    if (DEBUG) print_solution(parent2);
    if (mode == "MERGING")
    {
        return merging(parent1, parent2);
    }
    else if (mode == "1POINT")
    {
        return cross_1point(parent1, parent2);
    }
    else if (mode == "CUPCAP")
    {
        return cupcap(parent1, parent2);
    }
    else
    {
        TSolution* sol = new TSolution;
        if (0 == (rand() % 2))
        {
            //return parent1;
            for (int i = 0; i < P; ++i)
            {
                sol->individual.push_back(parent1->individual[i]);
            }
        }
        else
        {
            //return parent2;
            for (int i = 0; i < P; ++i)
            {
                sol->individual.push_back(parent2->individual[i]);
            }
        }
        sol->fitness = fitness(sol);
        return sol;
    }
}

void mutation(std::string mode, double prob, TSolution* &ind)
{
    double calc = (double) ((double)rand()/(double)RAND_MAX);
    if (calc <= prob)
    {
        if (mode == "SHAKING")
        {
            int i = rand() % P;
            i++;
            TSolution* temp = shake_rand(ind, i);
            ind->individual.clear();
            delete ind;
            ind = temp;
        }
        else if (mode == "NRAND")
        {
            int i = rand() % P;
            i++;
            if (DEBUG) print_solution(ind);
            TSolution* temp = shake_rand_neighborhood(ind, i);
            ind->individual.clear();
            delete ind;
            ind = temp;
        }
    }
}

int replacement(std::string mode, TSolution* offspring, std::vector<TSolution*> &pop, int pop_size, TSolution* &best, int worst_idx)
{
    if (mode == "ELITIST")
    {
        if (offspring->fitness < pop[worst_idx]->fitness && !find_in_pop(offspring, pop, pop_size))
        {
            pop[worst_idx]->individual.clear();
            delete pop[worst_idx];
            pop[worst_idx] = offspring;
            for (int i = 0; i < pop_size; ++i)
            {
                if (pop[i]->fitness > pop[worst_idx]->fitness)
                {
                    worst_idx = i;
                }
            }

            if (offspring->fitness < best->fitness)
            {
                best = offspring;
            }
        }
    }
    return worst_idx;
}

void replacement(std::string mode, std::vector<TSolution*> &childrens, int lambda, std::vector<TSolution*> &pop, int mu)
{
    if (DEBUG) print_population(pop, mu);
    if (DEBUG) print_population(childrens, lambda);

    if (mode == "mu,lambda")
    {
        std::vector<TSolution*> total;

        for (int i = 0; i < mu; ++i)
        {
            total.push_back(pop[i]);
        }
        sort(total.begin(),total.end(), [](const TSolution* x, const TSolution* y){ return (x->fitness < y->fitness);});
        for (int i = 0; i < lambda; ++i)
        {
            total[i + (mu - lambda)]->individual.clear();
            delete total[i + (mu - lambda)];
            total[i + (mu - lambda)] = childrens[i];
        }
        for (int i = 0; i < mu; ++i)
        {
            pop[i] = total[i];
        }
    }
    if (mode == "mu+lambda")
    {
        std::vector<TSolution*> total;

        for (int i = 0; i < lambda; ++i)
        {
            total.push_back(childrens[i]);
        }
        for (int i = 0; i < mu; ++i)
        {
            total.push_back(pop[i]);
        }
        sort(total.begin(),total.end(), [](const TSolution* x, const TSolution* y){ return (x->fitness < y->fitness);});
        for (int i = 0; i < mu; ++i)
        {
            pop[i] = total[i];
        }
        for (int i = 0; i < lambda; ++i)
        {
            total[i + mu]->individual.clear();
            delete total[i + mu];
        }
    }
    if (DEBUG) print_population(pop, mu);
}

TSolution* cross_1point(TSolution* p1, TSolution* p2)
{
    TSolution* sol = new TSolution;
    int point = rand() % P;
    log("cross_1point " + std::to_string(point));
    for (int i = 0; i < point; ++i)
    {
        sol->individual.push_back(p1->individual[i]);
    }
    for (int i = point; i < P; ++i)
    {
        if (sol->individual.begin() + i == std::find(sol->individual.begin(), sol->individual.begin() + i, p2->individual[i]))
        {
            sol->individual.push_back(p2->individual[i]);
        }
    }
    int idx = point;
    while ((sol->individual.size() != (unsigned) P) && (idx < P))
    {
        if (sol->individual.end() == std::find(sol->individual.begin(), sol->individual.end(), p1->individual[idx]))
        {
            sol->individual.push_back(p1->individual[idx]);
        }
        idx++;
    }
    idx = 0;
    while ((sol->individual.size() != (unsigned) P) && (idx < point))
    {
        if (sol->individual.end() == std::find(sol->individual.begin(), sol->individual.end(), p2->individual[idx]))
        {
            sol->individual.push_back(p2->individual[idx]);
        }
        idx++;
    }
    return sol;
}

TSolution* merging(TSolution *p1, TSolution *p2)
{
    // 1. Randomly generate theta \in [0;2pi].
    double theta=((double)rand()/(double)(2.0 * M_PI));
    // 2. Calculate vi =xi cos theta + yi sin theta for the locations of each parent.
    //    This rotates the axes and facilities by theta and vi is the projection of facility i on the rotated x-axis.
    std::vector<double> v1(P);
    std::vector<double> v2(P);
    int idx1;
    int idx2;
    for (int i = 0; i < P; ++i)
    {
        idx1 = p1->individual[i];
        v1[i] = facility_points[idx1][1] * cos(theta) + facility_points[idx1][2] * sin(theta);

        idx2 = p2->individual[i];
        v2[i] = facility_points[idx2][1] * cos(theta) + facility_points[idx2][2] * sin(theta);
    }
    // 3. Sort the vector {vi} for i=1,...,p for each parent.
    // 4. Select the smallest p2= [p/2] facilities from the first parent and the largest p-p2 facilities from the second parent.
    std::vector<double> lowestP2(P/2);
    kLowest(v2, P, (P/2), lowestP2);
    std::vector<double> largestP1(P - (P/2));
    kLargest(v1, P, P - (P/2), largestP1);
    TSolution* sol = new TSolution;
    for (int i = 0; i < P/2; ++i)
    {
        sol->individual.push_back(find_in_array(lowestP2[i], v2, p2->individual));
    }
    for (int i = 0; i < P - (P/2); ++i)
    {
        sol->individual.push_back(find_in_array(largestP1[i], v1, p1->individual));
    }
    sol->fitness = fitness(sol);
    return sol;
}

TSolution* cupcap(TSolution* p1, TSolution* p2)
{
    //log("CUPCAP");
    //if (DEBUG) print_solution(p1);
    //if (DEBUG) print_solution(p2);
    std::vector<int> fixG = intersection(p1->individual, p2->individual);
    std::vector<int> freeG = unionarray(p1->individual, p2->individual);

    //if (DEBUG) std::cout << fixG.size() << "\n";
    //if (DEBUG) std::cout << freeG.size() << "\n";
    minus(freeG, fixG);
    //if (DEBUG) std::cout << fixG.size() << "\n";
    //if (DEBUG) std::cout << freeG.size() << "\n";
    std::vector<int> sprime1;
    std::vector<int> sprime2;
    for (unsigned i = 0; i < fixG.size(); ++i)
    {
        sprime1.push_back(fixG[i]);
        sprime2.push_back(fixG[i]);
    }
    int s;
    //log("cp1");
    //if (DEBUG) std::cout << fixG.size() << "\n";
    //if (DEBUG) std::cout << freeG.size() << "\n";
    while(freeG.size() > 0) {
        if (DEBUG) std::cout << freeG.size() << "\n";
        s = rand() % freeG.size(); // return position
        sprime1.push_back(freeG[s]);
        freeG.erase(freeG.begin()+s);

        s = rand() % freeG.size(); // return position
        sprime2.push_back(freeG[s]);
        freeG.erase(freeG.begin()+s);
    }
    //log("cp2");
    TSolution *sol = new TSolution;
    for (unsigned i = 0; i < sprime1.size(); ++i)
    {
        sol->individual.push_back(sprime1[i]);
    }
    /*
    for (int i = 0; i < sprime2.size(); ++i)
    {
        sol->individual.push_back(sprime2[i]);
    }
    */
    //log("cp3");
    return sol;
}


/******************************************************************************
 * Test code to the Tabu Search
 * Author: Christian Cintrano
 * Date: 2018-05-29
 * Updated: 2018-05-29
 * Version: 1.0
 * Implementation of a TS to solve the p-median problem
 *
 * Rolland, E., Schilling, D. A., & Current, J. R. (1997). 
 * An efficient tabu search procedure for the p-Median Problem. 
 * European Journal of Operational Research, 96(2), 329–342. 
 * http://doi.org/10.1016/S0377-2217(96)00141-5
 * ****************************************************************************/
/*std::map<int, int> initialize_map_values(int value);
void initialize(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time, std::string mode = "moveADD", double param = 0.25); // same param than START
void choose_move(std::set<int> &S, int slack, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time, std::string mode = "randomADD-DROP");
void save(std::set<int> &S, TSolution &sol); //TODO no asegura mejora ya que se tienen en cuenta dos fitness distintos.
double fitness(std::set<int> &S, std::map<int, int> freq, double kindex);
void ADD_best(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time);
void DROP_best(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time);
bool is_tabu(int i, std::map<int, int> &add_time, int current_time, int tabu_time);
void update_tabu_list(int elem, std::map<int, int> &add_time, int tabu_time);
int set_tabu_time();

TSolution TS(int max_iterations, std::string gen_mode, double gen_param, std::string ls, double ls_param)// TODO: diferenciar N y F
{
    // Initialize:
    //Stable lterations : = 0.2 * Max Iteration, 
    double stable_iterations = 0.2 * max_iterations;
    // Iteration : = 1, 
    int iteration = 1;
    // Best_Solution : = Infinity, 
    double best_solution = std::numeric_limits<double>::infinity();
    TSolution best;
    best.individual = new int[P];
    // Slack : = O,
    int slack = 0;
    // Add_Time(v) := - infinity and Freq(v) := 0 for all nodes v.
    int tabu_time = 1;
    std::map<int, int> add_time = initialize_map_values(std::numeric_limits<int>::min());
    std::map<int, int> freq = initialize_map_values(0);
    // S:= {Empty}
    std::set<int> S;
    // k := Max_{ij}d_{ij}
    double kindex = 0;// TODO He supuesto que es la distancia máxima
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < F; ++j)
        {
            if (D[i][j] > kindex) kindex = D[i][j];
        }
    }
    double new_solution = std::numeric_limits<double>::infinity();
    // Step 0. Generate a feasible starting solution: 
    // While |S| < p, perform procedure ADD. 
    initialize(S, add_time, freq, kindex, tabu_time, iteration, gen_mode, gen_param);
    // WHILE Iteration < Max_Iteration do the following
    while(iteration < max_iterations) {
        // Step 1. Choose a candidate for moving using procedure CHOOSE_MOVE. 
        choose_move(S, slack, add_time, freq, kindex, tabu_time, iteration);
        new_solution = fitness(S, freq, kindex);
        // Step 2. Iteration := Iteration + 1 
        iteration++;
        bool improvement = false;
            // If |S| = p and NewSolution <Best_Solution then
        if (S.size() == (int) P && new_solution < best_solution)
        {
            improvement = true;
            // set Best_Solution:= NewSolution, 
            best_solution = new_solution;
            // set Slack := O, and 
            slack = 0;
            // save the best solution configuration. 
            save(S, best);
        }
            // otherwise Go to Step 3.
        // Step 3. If no improvement has been found in a number of iterations = Stable_Iterations X 2 
        if (!improvement && iteration == stable_iterations * 2)
        {
            // then Slack:= Slack + 1.
            slack++;
        }
        // Step 4. If no improvement has been found in a number of iterations = Stable_Iterations/2 
        if (!improvement && iteration == stable_iterations * 2)
        {
            // then Set Tabu_time:= a uniformly distributed random number in the range [1,p + 1]
            tabu_time = set_tabu_time();
        }
        // Step 5. If a feasible solution has been obtained, but no improvement in the best feasible solution has been found in a number of iterations = Stable__Iterations then 
        if (!improvement && S.size() == (int) P && iteration == stable_iterations)
        {
            // Iteration : = Max__Iteration 
            iteration = max_iterations;
        }
        else
        {
            // otherwise Iteration: = Iteration + 1.
            log(std::to_string(iteration));
            iteration++;
        }
    // END WHILE
    }
    // Step 6. To ensure that a local optima is found, run the NS local search heuristic using BestSolution as a starting solution. 
    local_search(best, ls, 1, ls_param); //TODO: Use NS.
    // END.
    return best;
}

double fitness(std::set<int> &S, std::map<int, int> freq, double kindex)
{
    double fitness = 0.0f;
    double current = 0.0f;
    double penalty = 0.0f;
    for (int i = 0; i < N; ++i)
    {
        double min = std::numeric_limits<double>::infinity();
        for(auto elem : S)
        {
            penalty = kindex * static_cast<double>(freq[elem]);
            current = D[i][elem] + penalty;
            if (current < min)
            {
               min = current;
            }
        }
        fitness = fitness + min;
    }
    return fitness;
}

// To ADD procedure
double fitness(std::set<int> &S, int facility, std::map<int, int> freq, double kindex)
{
    double fitness = 0.0f;
    double current = 0.0f;
    double penalty = 0.0f;
    double min;
    for (int i = 0; i < N; ++i)
    {
        penalty = kindex * static_cast<double>(freq[facility]);
        min = D[i][facility] + penalty;
        for(auto elem : S)
        {
            penalty = kindex * static_cast<double>(freq[elem]);
            current = D[i][elem] + penalty;
            if (current < min)
            {
               min = current;
            }
        }
        fitness = fitness + min;
    }
    return fitness;
}

// To DROP procedure
double fitness_drop(std::set<int> &S, int facility, std::map<int, int> freq, double kindex)
{
    double fitness = 0.0f;
    double current = 0.0f;
    double penalty = 0.0f;
    for (int i = 0; i < N; ++i)
    {
        double min = std::numeric_limits<double>::infinity();
        for(auto elem : S)
        {
            if(elem != facility)
            {
                penalty = kindex * static_cast<double>(freq[elem]);
                current = D[i][elem] + penalty;
                if (current < min)
                {
                   min = current;
                }
            }
        }
        fitness = fitness + min;
    }
    return fitness;
}

void save(std::set<int> &S, TSolution &sol)
{
    int i = 0;
    for(auto elem : S)
    {
        sol.individual[i] = elem;
        i++;
    }
    sol.fitness = fitness(sol);
}

std::map<int, int> initialize_map_values(int value)
{
    std::map<int, int> mapping;
    for (int i = 0; i < F; ++i)
    {
        mapping[i] = value;
    }
    return mapping;
}

void initialize(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time, std::string mode, double param)
{
    log("...initialize...", false);
    if (mode == "moveADD")
    {
        while(S.size() < (int) P) {
            ADD_best(S, add_time, freq, kindex, tabu_time, current_time);
        }
    }
    else
    {
        TSolution x = initial_solution(mode, param);
        if (DEBUG) print_solution(x);
        for (int i = 0; i < P; ++i)
        {
            S.insert(x.individual[i]);
        }
        delete x.individual;
    }
    log("end");
}

void choose_move(std::set<int> &S, int slack, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time, std::string mode) // TODO: añadir modos a ADD y DROP para contemplar random y best option
{
    log("..moving..", false);
    if (mode == "randomADD-DROP ")
    {
        if (S.size() < (int) (P-slack))
        {
            ADD_best(S, add_time, freq, kindex, tabu_time, current_time);
        }
        else if (S.size() < (int) (P-slack))
        {
            DROP_best(S, add_time, freq, kindex, tabu_time, current_time);
        }
        else
        {
            int flip;
            flip = rand() % 2 + 1;// assign random numbers
            if (flip == 1 && S.size() > 0) // The flip was heads.
                DROP_best(S, add_time, freq, kindex, tabu_time, current_time);
            else // The flip was tails
                ADD_best(S, add_time, freq, kindex, tabu_time, current_time);
        }
    }
    log("end");
}

int set_tabu_time()
{
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,P+1);
    return distribution(generator);
}

bool is_tabu(int i, std::map<int, int> &add_time, int current_time, int tabu_time)
{
    return add_time[i] >= (current_time - tabu_time);
}

void update_tabu_list(int elem, std::map<int, int> &add_time, int current_time)
{
    add_time[elem] = current_time;
}

void ADD_best(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time)
{
    double best_fitness = fitness(S, freq, kindex);
    double current_fitness = 0; 
    int best = -1;
    for (int i = 0; i < F; ++i)
    {
        if (!is_tabu(i, add_time, current_time, tabu_time) && S.find(i) == S.end())
        {
            current_fitness = fitness(S, i, freq, kindex);
            if (current_fitness < best_fitness)
            {
                best_fitness = current_fitness;
                best = i;
            }
        }
    }
    if (best != -1)
    {
        S.insert(best);
        freq[best] = freq[best] + 1;
        update_tabu_list(best, add_time, current_time);
    }
}

void DROP_best(std::set<int> &S, std::map<int, int> &add_time, std::map<int, int> &freq, int kindex, int tabu_time, int current_time)
{
    double best_fitness = fitness(S, freq, kindex);
    double current_fitness = 0; 
    int best = -1;

    for (int i = 0; i < F; ++i)
    {
        if (S.find(i) == S.end() && is_tabu(i, add_time, current_time, tabu_time))
        {
            current_fitness = fitness_drop(S, i, freq, kindex);
            if (current_fitness < best_fitness)
            {
                best_fitness = current_fitness;
                best = i;
            }
        }
    }
    if (best == -1)
    {
        for (int i = 0; i < F; ++i)
        {
            if (S.find(i) == S.end() && !is_tabu(i, add_time, current_time, tabu_time))
            {
                current_fitness = fitness_drop(S, i, freq, kindex);
                if (current_fitness < best_fitness)
                {
                    best_fitness = current_fitness;
                    best = i;
                }
            }
        }
    }
    if (best != -1)
    {
        S.erase(best);
    }
}
*/