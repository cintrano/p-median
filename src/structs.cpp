
#include <iostream>
#include <vector>
#include <algorithm>
#include "utilities.h"
#include "structs.h"

const int MAX_DISTANCE_SIZE = 500000;

int P, N, F;

double **D; // Distance matrix [customers x facilities] NxF
int **Neighborhood; // Neighborhood matrix [facilities x facilities] Fxk
double **points;
double **customer_points;
double **facility_points;
int k = 2; // Neighborhood size parameter
int *Nearest; // Closest facility id of each facility 
std::vector<int> individual_index;
// Auxiliars
int* tempFk;
double * tempNF;

void print_population(std::vector<TSolution*> pop, int size)
{
    for (int i = 0; i < size; ++i)
        {
            std::cout << "Individual " << i << ": ";
            print_solution(pop[i]);
        }
}

void print_solution(TSolution* x)
{
    for (int i = 0; i < P; ++i)
    {
        std::cout << x->individual[i] << '\t';
    }
    std::cout << " | " << x->fitness << std::endl;
}

double fitness_without_weights(TSolution* x)
{
    double fitness = 0.0f;
    for (int i = 0; i < N; ++i)
    {
        double min = D[i][x->individual[0]];
        for (int j = 1; j < P; ++j)
        {
           if (D[i][x->individual[j]] < min)
           {
               min = D[i][x->individual[j]];
           }
        }
        fitness = fitness + min;
    }
    return fitness;
}

double fitness(TSolution* x)
{
    double fitness = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double min = D[i][x->individual[0]];
        for (int j = 1; j < P; ++j)
        {
           if (D[i][x->individual[j]] < min)
           {
               min = D[i][x->individual[j]];
           }
        }
        fitness = fitness + (min * customer_points[i][0]);
    }
    return fitness;
}

bool ind_contains(TSolution *x, int elem)
{
    for (int i = 0; i < P; ++i)
    {
        if (x->individual[i] == elem)
        {
            return true;
        }
    }
    return false;
}

TSolution* find_best(std::vector<TSolution*> pop, int size)
{
    TSolution *best = pop[0];
    for (int i = 1; i < size; ++i)
    {
        if (pop[i]->fitness < best->fitness)
        {
            best = pop[i];
        }
    }
    return best;
}

bool find_in_pop(TSolution *sol, std::vector<TSolution*> pop, int size) // WARNING: Improve this function check the solution values
{ 
    for (int i = 0; i < size; ++i) {
        if (sol->fitness == pop[i]->fitness) {
            return true;
        }
    }
    return false;
}

int find_idx_of(TSolution *sol, std::vector<TSolution*> pop, int size)
{
    for (int i = 0; i < size; ++i) {
        if (sol->fitness == pop[i]->fitness) {
            return i;
        }
    }
    return -1;
}

void initialize_distance_matrix()
{
    if (N * F < MAX_DISTANCE_SIZE)
    {
        tempNF = new double[N * F];
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            //D[i] = new int[F];
            D[i] = (tempNF + i * F);
            for (int j = 0; j < F; ++j)
            {
                D[i][j] = distanceCalculate(points[i][1], points[i][2], points[j][1], points[j][2]);
            }
        }
    }
    else
    {
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            //D[i] = new int[F];
            D[i] = new double[F];
            for (int j = 0; j < F; ++j)
            {
                D[i][j] = distanceCalculate(points[i][1], points[i][2], points[j][1], points[j][2]);
            }
        }
    }

    for (int j = 0; j < P; ++j)
    {
        individual_index.push_back(j);
    }
}

void initialize_distance_matrix_two_points_array()
{
    if (N * F < MAX_DISTANCE_SIZE)
    {
        tempNF = new double[N * F];
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            //D[i] = new int[F];
            D[i] = (tempNF + i * F);
            for (int j = 0; j < F; ++j)
            {
                D[i][j] = distanceCalculate(customer_points[i][1], customer_points[i][2], facility_points[j][1], facility_points[j][2]);
            }
        }
    }
    else
    {
        D = new double*[N];
        for (int i = 0; i < N; ++i)
        {
            //D[i] = new int[F];
            D[i] = new double[F];
            for (int j = 0; j < F; ++j)
            {
                D[i][j] = distanceCalculate(customer_points[i][1], customer_points[i][2], facility_points[j][1], facility_points[j][2]);
            }
        }
    }

    for (int j = 0; j < P; ++j)
    {
        individual_index.push_back(j);
    }
}

void initialize_neighborhood() // WARNING Solo es valido para el caso NxN
{
    unsigned neg_count = 0;
    tempFk = new int[F * k];
    Neighborhood = new int*[F]; // dynamic array (size 10) of pointers to int

    std::vector<double> lineD(F);
    //std::vector<double> lineN(F);
    for (int i = 0; i < F; ++i)
    {
        //Neighborhood[i] = new int[k];
        Neighborhood[i] = (tempFk + i * k);
        std::vector<double> neighbors(k+1);
        for (int line_i = 0; line_i < F; ++line_i)
        {
            lineD[line_i] = D[i][line_i];
        }
        kLowest(lineD, N, k+1, neighbors); // first is the same element dist=0
        for (int j = 0; j < k; ++j)
        {
            int index = 0;
            int idx = -1;
            while(idx == -1) {
                if (D[i][index] == neighbors[j+1])
                {
                    idx = index;
                }
                index++;
            }
            //Neighborhood[i][j] = neighbors[j+1]; // Not take into account the same element
            Neighborhood[i][j] = idx;
            if (Neighborhood[i][j] > N)
            {
                std::cout << "Nerror " << Neighborhood[i][j] << "\n";
            }
            if (Neighborhood[i][j] < 0)
            {
                neg_count++;
            }
        }
        //delete neighbors;
    }
    if (neg_count != 0)
    {
        std::cout << "Nerror number of -1 is " << neg_count << "\n";
    }
    initialize_nearest();
}


std::vector<double> computeDistances(int i)
{
    std::vector<double> distances(F);
    for (int j = 0; j < F; ++j)
    {
        distances[j] = distanceCalculate(facility_points[i][1], facility_points[i][2], facility_points[j][1], facility_points[j][2]);
    }
    return distances;
}

void initialize_neighborhood_quadrants() // WARNING Solo es valido para el caso NxN
{
    tempFk = new int[F * k];
    Neighborhood = new int*[F]; // dynamic array (size 10) of pointers to int

    for (int i = 0; i < F; ++i)
    {
        //Neighborhood[i] = new int[k];
        Neighborhood[i] = (tempFk + i * k);

        std::vector<double> x = computeDistances(i); // The same is infiny

        std::vector<int> y(x.size());
        std::size_t n(0);
        std::generate(std::begin(y), std::end(y), [&]{ return n++; });

        std::sort(std::begin(y), 
                  std::end(y),
                  [&](int i1, int i2) { return x[i1] < x[i2]; } );

        int ntr, nrb, nbl, nlt;
        ntr = nrb = nbl = nlt = 0;
        int n_quadrant = k/4;
        std::vector<bool> check(F);
        int neighbor;
        for (int j = 0; j < F; ++j)
        {
            neighbor = y[j];
            if (ntr < n_quadrant && (facility_points[i][1] < facility_points[neighbor][1] && facility_points[i][2] < facility_points[neighbor][2]))
            {
                ntr++;
                check[neighbor] = true;
            }
            else if (nrb < n_quadrant && (facility_points[i][1] > facility_points[neighbor][1] && facility_points[i][2] < facility_points[neighbor][2]))
            {
                nrb++;
                check[neighbor] = true;
            }
            else if (nbl < n_quadrant && (facility_points[i][1] > facility_points[neighbor][1] && facility_points[i][2] > facility_points[neighbor][2]))
            {
                nbl++;
                check[neighbor] = true;
            }
            else if (nlt < n_quadrant && (facility_points[i][1] < facility_points[neighbor][1] && facility_points[i][2] > facility_points[neighbor][2]))
            {
                nlt++;
                check[neighbor] = true;
            }
            else 
            {
                check[neighbor] = false;
            }
        }
        int nnearest = ntr + nrb + nbl + nlt;
        int j = 1;
        while(nnearest != k && j < F) 
        {
            neighbor = y[j];
            if (!check[neighbor])
            {
                check[neighbor] = true;
                nnearest++;
            }
            j++;
        }

        int idx = 0;
        int pos = 0;

        while(pos < k) {
            if (check[idx])
            {
                Neighborhood[i][pos] = idx;
                pos++;
            }
            idx++;
        }
    }
    initialize_nearest();
}

void initialize_nearest()
{
    Nearest = new int[F];
    int near, current;
    double d_near, d_current;
    for (int i = 0; i < F; ++i)
    {
        near = Neighborhood[i][0];
        d_near = distanceCalculate(facility_points[i][1], facility_points[i][2], facility_points[near][1], facility_points[near][2]);
        for (int j = 1; j < k; ++j)
        {
            current = Neighborhood[i][j];
            d_current = distanceCalculate(facility_points[i][1], facility_points[i][2], facility_points[current][1], facility_points[current][2]);
            if (d_current < d_near)
            {
                near = current;
                d_near = d_current;
            }
        }
        Nearest[i] = near;
    }
}

void print_distances()
{
    std::cout << "D: \n";
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < F; ++j)
        {
            std::cout << D[i][j] << '\t';
        }
        std::cout << '\n';
    }
}

void print_neighborhood()
{
    std::cout << "Neighborhoods: \n";
    for (int i = 0; i < F; ++i)
    {
        std::cout << i << " : ";
        for (int j = 0; j < k; ++j)
        {
            std::cout << Neighborhood[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void remove_structs()
{
    for (int i = 0; i < N; ++i)
    {
        if (N * F >= MAX_DISTANCE_SIZE)
            delete [] D[i];
        delete [] customer_points[i];
    }
    delete [] D;
    delete [] customer_points;

    for (int i = 0; i < F; ++i)
    {
        // delete [] Neighborhood[i];
        delete [] facility_points[i];
    }
    delete [] Neighborhood;
    delete [] facility_points;

    delete [] Nearest;
    if (N * F < MAX_DISTANCE_SIZE)
        delete [] tempNF;
    delete [] tempFk;
}
