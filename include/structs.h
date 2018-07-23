#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>

struct TSolution {
  std::vector<int> individual;
  double fitness;
};

double fitness(TSolution* x);

extern const int MAX_DISTANCE_SIZE;
extern int P, N, F;

extern double **D;
extern int **Neighborhood;
extern double **points;
extern double **customer_points;
extern double **facility_points;
extern int k;
extern int *Nearest;
extern std::vector<int> individual_index;

bool ind_contains(TSolution* x, int elem);
TSolution* find_best(std::vector<TSolution*> pop, int size);
bool find_in_pop(TSolution *sol, std::vector<TSolution*> pop, int size); // TODO: Improve this function check the solution values
int find_idx_of(TSolution *sol, std::vector<TSolution*> pop, int size);
void print_population(std::vector<TSolution*> pop, int size);
void print_solution(TSolution* x);

void initialize_distance_matrix();
void initialize_distance_matrix_two_points_array();
void initialize_neighborhood();
void initialize_neighborhood_quadrants();
void initialize_nearest();
void print_distances();
void print_neighborhood();

void remove_structs();

#endif