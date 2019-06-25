
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <limits>
#include <algorithm>    /* find, sort */
#include <chrono>       /* measure time */
#include <set>
#include <map>
#include <iterator>
#include "structs.h"
#include "localsearches.h"
#include "utilities.h"

using namespace std;

// Default values
double BSSS_EPSILON = 0.0000000001; 

int BSSS(std::vector<double> &Dist, double epsilon);
//void IMP(TSolution* X, double param);
//void IALT(TSolution* X, double L);
//void ISRW(TSolution* x, int alpha);
//void FI(TSolution* X);


void local_search(TSolution* x, string ls, double param, int max_time)
{
    log("-- local search " + ls, false);
    
    if (ls == "IALT")
    {
        IALT(x, param, max_time); // param = Laux
    }
    else if (ls == "IMP")
    {
        IMP(x, param, max_time); // param = epsilon
    }
    else if (ls == "ISRW")
    {
        ISRW(x, param, max_time); // param = ALPHA
    }
    else if (ls == "FI")
    {
        FI(x, max_time);
    }
    else if (ls == "NONE")
    {
        // No local search perform
    }
    else
    {
        log("ERROR: Invalid local search");
    }
    log("...end");
}

//------------------------------------------------
// IMP
//------------------------------------------------

void printDis(double *Dist)
{
    for (int i = 0; i < N; ++i)
    {
        std::cout << Dist[i] << '\t';
    }
    std::cout << '\n';
}

void IMP(TSolution* X, double param, int max_time) {
    std::vector<double> Dist(N);
    int kindex;
    // 1. Obtain a starting solution.
    // 2. Set index=0.
    bool index = true;
    int facility_location;
    int bsss_index;
    bool in_sol;
    double previous_fitness;
    double new_fitness;
    double min;

    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;
    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();


    while(index && (current_time < max_time)) {
        index = false;
        // 3. Repeat the following for each facility k in random order.
        std::random_shuffle(individual_index.begin(), individual_index.end());
        for (int pos = 0; pos < P; pos++) {
            kindex = individual_index[pos];
            // (a) Calculate Di by (8) for i=1,…,n .
            for (int i = 0; i < N; ++i) { // distancias del customer
                min = std::numeric_limits<double>::infinity();
                for (int j = 0; j < P; ++j) { // j facility
                    if (j != kindex && D[i][X->individual[j]] < min)
                    {
                        min = D[i][X->individual[j]];
                    }
                }
                Dist[i] = min;
            }
            facility_location = X->individual[kindex];
            // (b) Relocate facility k by solving the LD1 problem.
            //log("---");
            bsss_index = BSSS(Dist, param);
            //log("__________");
            //log(std::to_string(kindex));
            //log(std::to_string(bsss_index));
            in_sol = ind_contains(X, bsss_index);
            previous_fitness = X->fitness;
            if (!in_sol)
            {
                X->individual[kindex] = bsss_index;
                new_fitness = fitness(X);
                if (new_fitness > previous_fitness)
                {
                    X->individual[kindex] = facility_location;
                }
                else
                {
                    X->fitness = new_fitness;
                }
            }
            // (c) If the location of facility k changed, set index=1.
            if (facility_location != X->individual[kindex])
            {
                index=true;
            }
        }
        // 4. If index=0, stop with the solution being the current locations of the p facilities.

        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    }
    //log("pre dist");
    Dist.clear();
    //log("pos dist");
    // 5. Otherwise, go to Step 2.
    X->fitness = fitness(X);
    if (DEBUG) print_solution(X);
    log("After fitness");
}


struct TSquare {
  //int individual[P];
  double top, right, bottom, left;
  int UB;
  double UB_fitness;
  double LB;
};

/*
inline bool operator <(const TSquare* &lhs, const TSquare* &rhs)
{
  // return lhs.UB < rhs.UB;
  //return lhs.UB < rhs.UB && lhs.top < rhs.top && lhs.right < rhs.right && lhs.bottom < rhs.bottom && lhs.left < rhs.left;
  return (lhs->UB < rhs->UB || lhs->LB < rhs->LB) !=0;
}
*/

int point_nearest(double lat, double lon) // TODO: Optimize this operation to be able to use Weber operations
{
    int nearest = 0;
    double dist = distanceCalculate(facility_points[nearest][1], facility_points[nearest][2], lat, lon);
    double new_dist;
    for (int i = 1; i < F; ++i)
    {
        new_dist = distanceCalculate(facility_points[i][1], facility_points[i][2], lat, lon);
        if (new_dist < dist)
        {
            dist = new_dist;
            nearest = i;
        }
    }
    return nearest;
}

bool inSquare(TSquare* &square, int i)
{
    return facility_points[i][2] >= square->left 
           && facility_points[i][2] <= square->right
           && facility_points[i][1] >= square->bottom
           && facility_points[i][1] <= square->top;
}

int point_nearest(TSquare* &square, double lat, double lon) // TODO: Optimize this operation to be able to use Weber operations
{
    int nearest = -1;
    double dist = std::numeric_limits<double>::infinity();
    double new_dist;
    for (int i = 0; i < F; ++i)
    {
        if (inSquare(square, i))
        {
            new_dist = distanceCalculate(facility_points[i][1], facility_points[i][2], lat, lon);
            if (new_dist < dist)
            {
                dist = new_dist;
                nearest = i;
            }
        }
    }
    return nearest;
}

double LD1Problem(int index, const std::vector<double> &Dist)
{
    double fitness = 0.0;
    for (int i = 0; i < N; ++i)
    {
        int min = D[i][index];
        if (Dist[i] < min)
        {
            min = Dist[i];
        }
        fitness = fitness + (double) min;
    }
    return fitness;
}

double calculateLB(TSquare* &square, const std::vector<double> &Dist)
{
    double fitnessTL = 0.0f;
    double fitnessTR = 0.0f;
    double fitnessBR = 0.0f;
    double fitnessBL = 0.0f;
    int tr = point_nearest(square, square->top, square->right);
    int br = point_nearest(square, square->bottom, square->right);
    int bl = point_nearest(square, square->bottom, square->left);
    int tl = point_nearest(square, square->top, square->left);

    if (tr != -1) fitnessTR = LD1Problem(tr, Dist);
    if (br != -1) fitnessBR = LD1Problem(br, Dist);
    if (bl != -1) fitnessBL = LD1Problem(bl, Dist);
    if (tl != -1) fitnessTL = LD1Problem(tl, Dist);
    double fitness = std::numeric_limits<double>::infinity();
    if (tr != -1 && fitness > fitnessTR)
    {
        fitness = fitnessTR;
    }
    if (br != -1 && fitness > fitnessBR)
    {
        fitness = fitnessBR;
    }
    if (bl != -1 && fitness > fitnessBL)
    {
        fitness = fitnessBL;
    }
    if (tl != -1 && fitness > fitnessTL)
    {
        fitness = fitnessTL;
    }
    if (fitness == 0.0f)
    {
        fitness = std::numeric_limits<double>::infinity();
    }
    return fitness;
}

void calculateUB(TSquare* &square, const std::vector<double> &Dist)
{
    double x, y;
    x = square->left + ((square->right - square->left) / 2.0);
    y = square->bottom + ((square->top - square->bottom) / 2.0);
    square->UB = point_nearest(square, y, x);
    square->UB_fitness = LD1Problem(square->UB, Dist);
}

TSquare* newSquare(double top, double right, double bottom, double left, const std::vector<double> &Dist)
{
    TSquare* square = new TSquare;
    square->top = top;
    square->right = right;
    square->bottom = bottom;
    square->left = left;
    calculateUB(square, Dist);
    square->LB = calculateLB(square, Dist);
    return square;
}

string printSquare(TSquare* &square)
{
    return to_string(square->UB) + "_" + to_string(square->LB);// + "[" + to_string(square.top) + "," + to_string(square.right) +  "," + to_string(square.bottom) +  "," + to_string(square.left) + "]";
}

void divideSquare(TSquare* &square, std::set<TSquare*> &set, const std::vector<double> &Dist)
{
    // |-------|
    // | 1 | 2 |
    // |---ub--|
    // | 3 | 4 |
    // |-------|

    TSquare* s1;
    TSquare* s2;
    TSquare* s3;
    TSquare* s4;
    s1 = newSquare(square->top, facility_points[square->UB][2], facility_points[square->UB][1], square->left, Dist);
    s2 = newSquare(square->top, square->right, facility_points[square->UB][1], facility_points[square->UB][2], Dist);
    s3 = newSquare(facility_points[square->UB][1], facility_points[square->UB][2], square->bottom, square->left, Dist);
    s4 = newSquare(facility_points[square->UB][1], square->right, square->bottom, facility_points[square->UB][2], Dist);
    set.insert(s1);
    set.insert(s2);
    set.insert(s3);
    set.insert(s4);
}

bool sameSquare(const TSquare* s1, const TSquare* s2)
{
    return s1->UB == s2->UB && s1->LB == s2->LB &&
           s1->top == s2->top && s1->right == s2->right &&
           s1->bottom == s2->bottom && s1->left == s2->left;
}

int BSSS(std::vector<double> &Dist, double epsilon) {
    // 1. The set of squares consists of a square enclosing the feasible area. The best upper bound UB∗ is the value of the objective function at the center of the square. A lower bound LB in the square is calculated.
    std::set<TSquare*> squares;
    std::set<TSquare*>::iterator it;
    int UBstar;
    double UBstart_fitness;

    double max_lat, max_lon, min_lat, min_lon, lat, lon;
    max_lat = min_lat = facility_points[0][1];
    max_lon = min_lon = facility_points[0][2];
    for (int i = 1; i < F; ++i)
    {
        lat = facility_points[i][1];
        lon = facility_points[i][2];
        if (lat > max_lat)
        {
            max_lat = lat;
        }
        if (lat < min_lat)
        {
            min_lat = lat;
        }
        if (lon > max_lon)
        {
            max_lon = lon;
        }
        if (lon < min_lon)
        {
            min_lon = lon;
        }
    }
    TSquare* first_square = new TSquare;
    first_square->top = max_lat;
    first_square->right = max_lon;
    first_square->bottom = min_lat;
    first_square->left = min_lon;
    calculateUB(first_square, Dist);
    first_square->LB = calculateLB(first_square, Dist);

    log(printSquare(first_square));

    UBstar = first_square->UB;
    UBstart_fitness = first_square->UB_fitness;
    squares.insert(first_square);
    double bound;
    while(squares.size() != 0) {//!squares.empty()) {
    // 2. The square with the smallest LB is selected and four small squares are constructed by two perpendicular lines through its center parallel to its sides.
        TSquare* square_lower_LB;

    //log(printSquare(square_lower_LB));
        double lower_LB = std::numeric_limits<double>::infinity();

        for (auto elem : squares)
        {
            if (elem->LB < lower_LB)
            {
               square_lower_LB = elem;
               lower_LB = elem->LB;
            }
        }
        divideSquare(square_lower_LB, squares, Dist); 
    // 3. For each of the small squares an upper bound UB (the value of the objective function at the center of the square) and a lower bound LB are calculated.
    // Done in the creation
 
    // 4. The best upper bound UB∗ is updated if necessary.
    // TODO: Change inneficient way to update the best
        for (auto elem : squares)
        {
            if (elem->UB_fitness < UBstart_fitness)
            {
                UBstart_fitness = elem->UB_fitness;
                UBstar = elem->UB;
            }   
        }
 
    // 5. The big square and all squares for which LB>=UB∗(1−epsilon) are discarded from the set of squares. All other squares remain or are added to the set of squares.
        bound = UBstart_fitness * (1.0f - epsilon );
        for (it=squares.begin(); it!=squares.end(); ++it)
        {
            if ((*it)->LB >= bound || sameSquare((*it), square_lower_LB))
            {
                squares.erase(it);
            }
        }
 
    // 6. If the set of squares is not empty, go to step 2.
    }
    // 7. Stop with UB∗ as the solution. UB∗ is within a relative accuracy epsilon of the optimum.
    squares.clear();
    return UBstar;
}

//------------------------------------------------
// IALT
//------------------------------------------------

/*
bool containsZero(int vec[], int size) {
    for (int i = 0; i < size; ++i)
    {
        if(vec[i] == 0) {
            return true;
        }
    }
    return false;
}
*/
bool containsZero(std::vector<int> &vec, int size) {
    for (int i = 0; i < size; ++i)
    {
        if(vec[i] == 0) {
            return true;
        }
    }
    return false;
}

void calculateS(const TSolution* x, std::map<int,std::set<int>> &S)
{
    S.clear();
    double min, current_min;
    int facility;
    for (int i = 0; i < N; ++i)
    {
        facility = x->individual[0];
        min = D[i][x->individual[0]];
        for (int j = 1; j < P; ++j)
        {
            current_min = D[i][x->individual[j]];
            if (current_min < min)
            {
                min = current_min;
                facility = x->individual[j];
            }
        }
        /*
        if (S.find(facility) != S.end()) {
            std::set< int>& s_ref = S[facility];
            s_ref.insert( i);
        } else {
            std::set< int> s;
            S.insert(std::make_pair(facility, s));
        }
        */
        S[facility].insert(i);
    }
}

void initializeA(int *a)
{
    int min, current_min, facility;
    for (int i = 0; i < N; ++i)
    {
        facility = 0;
        min = D[i][0];
        for (int j = 1; j < F; ++j)
        {
            current_min = D[i][j];
            if (current_min < min)
            {
                min = current_min;
                facility = j;
            }
        }
        a[i] = facility;
    }
}

void initializeA(int *a, std::map<int,std::set<int>> &map)
{
    for ( const auto &myPair : map )
    {
        for(auto n : myPair.second)
        {
            a[n] = myPair.first;
        }
    }
}

double* centroid(std::set<int> set)
{
    double x, y;
    x = y = 0;
    for (auto elem : set)
    {
        x += facility_points[elem][2];
        y += facility_points[elem][1];
    }
    // std::cerr << "\tx=" << x << " \t y=" << y << '\n';

    double* yx = new double[2];
    yx[0] = y / (double) set.size();
    yx[1] = x / (double) set.size();
    // std::cerr << "yx = " << yx[0] << " " << yx[1] << '\n';
    return yx;
}

int WeberSolver(std::map<int,std::set<int>> &S, int j)
{
    double* yx = centroid(S[j]);
                // std::cerr << j << "_______________ " << yx[0] << " " << yx[1] << '\n';
    int i = point_nearest(yx[0], yx[1]);
    //std::cerr << "_______________ " << j << " " << distanceCalculate(yx[0], yx[1], facility_points[j][1], facility_points[j][2])<< '\n';
    //std::cerr << "_______________ " << i << " " << distanceCalculate(yx[0], yx[1], facility_points[i][1], facility_points[i][2])<< '\n';
    delete yx;
    return i;
}

void calculateD(TSolution* X, int **Dist, int *a, int *b)
{
    // Customer i
    int first = 0, second = 0;
    int d_first, d_second;
    int d_current;
    for (int i = 0; i < N; ++i)
    {
        d_first = d_second = std::numeric_limits<int>::max();
        for (int j = 0; j < P; ++j)
        {
            d_current = D[i][X->individual[j]];
            if (d_current < d_first)
            {
                // first becomes second
                second = first;
                d_second = d_first;

                // new first
                first = j;
                d_first = d_current;
            }
            else if (d_current < d_second)
            {
                second = j;
                d_second = d_current;
            }
        }
        a[i] = first;
        b[i] = second;
        Dist[i][0] = d_first;
        Dist[i][1] = d_second;
    }
}

int find_index(int arr[], int len, int seek)
{
    for (int i = 0; i < len; ++i)
    {
        if (arr[i] == seek) return i;
    }
    return -1;
}

/**
 * param:
 * delta: distances array
 * indexes: facility indexes array for earch i-delta
 * sorted: sort array of index;
 */
void sortDeltas(int *delta, int *indexes, int *sorted)
{
    std::vector<int> delta_sorted(N);
    for (int i = 0; i < N; ++i)
    {
        delta_sorted[i] = delta[i];
    }
    sort(delta_sorted.begin(), delta_sorted.begin() + N);
    int idx;
    for (int i = 0; i < N; ++i)
    {
        idx = find_index(delta, N, delta_sorted[i]);
        if (idx == -1)
        {
            cerr << "ERROR IN SORT\n";
        }
        sorted[i] = indexes[idx];
    }
    //delete delta_sorted;
}

void printS(std::map<int,std::set<int>> myset) {
    std::map<int,std::set<int>>::iterator it;
    for(it = myset.begin(); it != myset.end(); ++it)
    {
        cout<<it->first<<": {"; //it->first gives you the key of the map.

        //it->second is the value -- the set. Iterate over it.
        for (set<int>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++)
            std::cerr<<*it2<<" ";
        std::cerr<<"}\n";
    }
}

void IALT(TSolution* X, double L, int max_time)
{

    TSolution* Xprime = new TSolution();

    std::map<int,std::set<int>> S;
    calculateS(X, S); // |S| = P, demand points 

    int *a= new int[N]; // Facility closest to the demand point i
    initializeA(a);
    std::vector<int> I(P);
    for (int i = 0; i < P; ++i)
    {
        I[i] = 0;
    }
    bool first_part_algorithm = true;


        // TIMER
        int current_time;
        std::chrono::steady_clock::time_point t_start, t_current;
        t_start= std::chrono::steady_clock::now();
        t_current= std::chrono::steady_clock::now();   
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();

    do { // Step 9        
        while(containsZero(I, P) || !first_part_algorithm) { // Step 5
            //printS(S);
            if (first_part_algorithm) // Condition to avoid these steps in the step 9 phase
            {
                // Step 1
                int j = rand() % P;
                while(I[j] == 1) {
                    j = rand() % P;
                }
                int previous = X->individual[j];
                int Xstar = WeberSolver(S, X->individual[j]);
                if (!ind_contains(X, Xstar))
                {
                    X->individual[j] = Xstar;
                    double new_fitness = fitness(X);
                    if(X->fitness < new_fitness)
                    {
                        X->individual[j] = previous;
                    } else
                    {
                        X->fitness = new_fitness;
                    }
                }
                //std::cerr << "\t\t" << Xstar << "\t"  << X->individual[j] << "\t" << X->fitness << "\t" << fitness(X) << "\n";
                I[j] = 1;
            }
            // Step 3
            int *aprime = new int[N]; // Facility closest to the demand point i
            initializeA(aprime, S);
            calculateS(X, S);
            // Step 4
            for (int i = 0; i < N; ++i)
            {
                if (aprime[i] != a[i])
                {
                    for (int kqux = 0; kqux < P; ++kqux)
                    {
                       if (I[kqux] == a[i] || I[kqux] == aprime[i])
                       {
                           I[kqux] = 0;
                       }
                    }
                    
                    a[i] = aprime[i];
                }
            }
        // Step 5
            first_part_algorithm = true;
            delete aprime;
        }
        // Step 6
        int **Dist;
        int* temp = new int[N * 2];
        Dist = new int*[N];
        for (int indez = 0; indez < N; ++indez)
        {
            Dist[indez] = (temp + indez * 2);
            
        }
        int *aclosest = new int[N];
        int *b = new int[N];
        calculateD(X, Dist, aclosest, b);
                
        int *delta = new int[N];
        for (int i = 0; i < N; ++i)
        {
            delta[i] = Dist[i][1] - Dist[i][0];
        }
        // Step 7
        int *index_sort = new int[N];
        sortDeltas(delta, aclosest, index_sort);
        int kindex = 0;
        // Step 8
        do {
            std::map<int,std::set<int>> SprimeA;
            std::map<int,std::set<int>> SprimeB;

            std::set<int> auxiliar_setA;
            for (auto elem : S[aclosest[index_sort[kindex]]])
            {
                if (elem != index_sort[kindex])
                {
                    auxiliar_setA.insert(elem);
                }
            }
            SprimeA[aclosest[index_sort[kindex]]] = auxiliar_setA;
            std::set<int> auxiliar_setB;
            auxiliar_setB.insert(index_sort[kindex]);
            for (auto elem : S[b[index_sort[kindex]]])
            {
                auxiliar_setB.insert(elem);
            }
            SprimeB[b[index_sort[kindex]]] = auxiliar_setB;

            int XprimeA = WeberSolver(SprimeA, aclosest[index_sort[kindex]]);
            int XprimeB = WeberSolver(SprimeB, b[index_sort[kindex]]);
            for (int i = 0; i < P; ++i)
            {
                Xprime->individual.push_back(X->individual[i]);
            }
            Xprime->individual[index_sort[kindex]] = XprimeA;
            Xprime->individual[index_sort[kindex]] = XprimeB;
            Xprime->fitness = fitness(Xprime);

            // Step 9
            if (Xprime->fitness < X->fitness)
            {
                first_part_algorithm = false;
                for (int i = 0; i < P; ++i)
                {
                    X->individual[i] = Xprime->individual[i];
                    X->fitness = Xprime->fitness;
                }
            }
            // Step 10
            if (Xprime->fitness >= X->fitness)
            {
                kindex = kindex + 1;
            }

            SprimeA.clear();
            SprimeB.clear();
            auxiliar_setA.clear();
            auxiliar_setB.clear();
        } while(Xprime->fitness >= X->fitness && kindex < L); // Condition of 10
        delete [] aclosest;
        delete [] b;
        delete [] delta;
        delete [] index_sort;
        delete [] Dist;

        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 

    } while(Xprime->fitness < X->fitness && (current_time < max_time)); // Condition of 9

    delete [] a;
    delete Xprime;
    S.clear();

    X->fitness = fitness(X);
}



//------------------------------------------------
// IS-RW
//------------------------------------------------

std::set<int> nearest(const TSolution* x, int index, int alpha)
{
    int *near = new int[alpha];
    int count = 1;
    near[0] = x->individual[0];
    // TODO con la estructura del k-d tree se podría optimizar mucho esta operacion
    for(int i = 1; i < P; ++i)
    {
        if (i != index)
        {
            if (count != alpha)
            {
                near[count] = x->individual[i];
                count++;
            }
            else
            {
                int max = D[x->individual[index]][near[0]];
                int index_max = 0;
                int current_dist;
                for (int j = 1; j < alpha; ++j)
                {
                    current_dist = D[x->individual[index]][near[j]];
                    if (current_dist > max)
                    {
                        max = current_dist;
                        index_max = j;
                    }
                }
                if (D[x->individual[index]][x->individual[i]] < max)
                {
                    near[index_max] = i;
                }
            }
        }
    }
    std::set<int> near_set;
    for (int i = 0; i < alpha; ++i)
    {
        near_set.insert(near[i]);
    }
    delete near;
    return near_set;
}

/**
 * The customers link to the facilities x_near in the solution x
 * x_near is a subset of x.individual
 */
std::set<int> demand_points(const TSolution* x, std::set<int> &x_near)
{
    std::set<int> customers;
    int dmin, index, current;
    bool is_in;
    for (int i = 0; i < N; ++i)
    {
        index = x->individual[0];
        dmin = D[i][index];
        for (int j = 1; j < P; ++j)
        {
            current = D[i][x->individual[j]];
            if (current < dmin)
            {
                dmin = current;
                index = x->individual[j];
            }
        }
        is_in = x_near.find(index) != x_near.end();
        if (is_in)
        {
            customers.insert(i);
        }
    }
    return customers;
}

std::set<int> demand_points(const TSolution* x, int point)
{
    std::set<int> customers;
    int dmin, index, current;
    for (int i = 0; i < N; ++i)
    {
        index = x->individual[0];
        dmin = D[i][index];
        for (int j = 1; j < P; ++j)
        {
            current = D[i][x->individual[j]];
            if (current < dmin)
            {
                dmin = current;
                index = x->individual[j];
            }
        }
        if (index == x->individual[point])
        {
            customers.insert(i);
        }
    }
    return customers;
}

int FindBestCustomer(int j, TSolution* &x, int alpha, int *a1, int *a2)
{
    int best = 0, loss;
    loss = std::numeric_limits<int>::max();
    std::set<int> x_near = nearest(x, j, alpha);
    std::set<int> c_near = demand_points(x, x_near);
    std::set<int> c_facj = demand_points(x, j);

    int w, g;
    for(auto i : c_facj)
    {
        if (i != x->individual[j])
        {
            w = 0;
            g = 0;
            for(auto z : c_facj)
            {
                if (D[z][i] < D[z][a1[z]])
                {
                    w = w + D[z][i] - D[z][a1[z]];
                }
                else
                {
                    if (D[z][i] < D[z][a2[z]])
                    {
                        g = g + D[z][i] - D[z][a1[z]];
                    }
                    else
                    {
                        g = g + D[z][a2[z]] - D[z][a1[z]];
                    }
                }
            }
            for(auto z : c_near)
            {
                if (D[z][i] < D[z][a1[z]])
                {
                    w = w + D[z][i] - D[z][a2[z]];
                }
            }
            w = w + g;
            if (w < loss)
            {
                best = i;
            }
        }
    }
    x_near.clear();
    c_near.clear();
    c_facj.clear();
    return best;
}

/**
 * in an array of customers the values are the nearest (and the second nearest) facility for each a_i customer
 * a1 first nearest
 * a2 second nearest
 */
void nearFacility2Customer(int *a_first, int *a_second)
{
    int first = 0, second = 0;
    double d_first, d_second;
    double d_current;
    for (int i = 0; i < N; ++i)
    {
        d_first = d_second = std::numeric_limits<int>::max();
        for (int j = 0; j < F; ++j)
        {
            d_current = D[i][j];
            if (d_current < d_first)
            {
                // first becomes second
                second = first;
                d_second = d_first;

                // new first
                first = j;
                d_first = d_current;
            }
            else if (d_current < d_second)
            {
                second = j;
                d_second = d_current;
            }
        }
        a_first[i] = first;
        a_second[i] = second;
    }
}

void ISRW(TSolution* x, int alpha, int max_time)
{
    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;
    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();

    int j = rand() % P;
    int *a1 = new int[N];
    int *a2 = new int[N];
    nearFacility2Customer(a1, a2);
    t_current= std::chrono::steady_clock::now();  
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    if (current_time < max_time)
        x->individual[j] = FindBestCustomer(j, x, alpha, a1, a2);
    delete a1;
    delete a2;
    x->fitness = fitness(x);
}


//------------------------------------------------
// Fast Interchange (FI)
// Whitaker R (1983) A fast algorithm for the greedy interchange for large-scale clustering and median location problems. INFOR 21:95–108
//------------------------------------------------

void FI(TSolution* x, int max_time) // TODO: Optimizar estructura de vecinos para tener los mas cercanos
{
    double fitness_old = x->fitness;

    int v, iter;
    int neighbor_best;
    double fitness_new;

    // TIMER
    int current_time;
    std::chrono::steady_clock::time_point t_start, t_current;
    t_start= std::chrono::steady_clock::now();
    t_current= std::chrono::steady_clock::now();   
    current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count();

    bool restart = true;
    while(restart && (current_time < max_time)) {
        restart = false;
 
        std::random_shuffle(individual_index.begin(), individual_index.end());

        iter = 0;
        while(!restart && iter < P)
        {
            v = x->individual[individual_index[iter]]; // size of the neighborhood between facilities
            neighbor_best = Nearest[v];

            if (!ind_contains(x, neighbor_best))
            {
                x->individual[individual_index[iter]] = neighbor_best;
                fitness_new = fitness(x);
                if (fitness_new < fitness_old)
                {
                    fitness_old = fitness_new;
                    restart = true;
                }
                else
                {
                   x->individual[individual_index[iter]] = v; 
                }
            }

            iter++;
        }

        t_current= std::chrono::steady_clock::now();  
        current_time = std::chrono::duration_cast<std::chrono::seconds> (t_current - t_start).count(); 
    }
    x->fitness = fitness_old;
}
