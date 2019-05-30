/******************************************************************************
 * Utilities functions
 * Author: Christian Cintrano
 * Date: 2018-03-26
 * Updated: 2019-03-26
 * Version: 1.0
 * 
 * ****************************************************************************/
#include <iostream>
#include <string>
#include <algorithm>    /* find */
#include <math.h>
#include <bits/stdc++.h>
#include "utilities.h"

using namespace std;

bool DEBUG = false;
bool SAVING = false;


void log(const std::string& text, bool break_line)
{
    if (DEBUG)
    {
        std::cerr << text;
        if (break_line)
        {
            std::cerr << std::endl;
        }
    }
}
void log(const char*& text, bool break_line)
{
    if (DEBUG)
    {
        std::cerr << text;
        if (break_line)
        {
            std::cerr << std::endl;
        }
    }
}

/*
void kLargest(double *array, int n, int k, double *p)
{
    double arr[n];
    for (int i = 0; i < n; ++i)
        arr[i] = array[i];
    // Sort the given array arr in reverse 
    // order.
    sort(arr, arr+n, greater<int>());
    //sort(arr, arr+n, less<int>());
 
    // Print the first kth largest elements
    for (int i = 0; i < k; i++)
         p[i] = arr[i];
}
*/
void kLargest(std::vector<double> &array, int n, int k, std::vector<double> &p)
{
    std::vector<double> arr(n);
    for (int i = 0; i < n; ++i)
        arr[i] = array[i];
    // Sort the given array arr in reverse 
    // order.
    std::sort(arr.begin(), arr.begin()+n, greater<int>());
    //sort(arr, arr.begin()+n, greater<int>());
    //sort(arr, arr+n, less<int>());
 
    // Print the first kth largest elements
    for (int i = 0; i < k; i++)
         p[i] = arr[i];
}

/**
 * params:
 * array: origin array
 * n: size of array n
 * k: number of elements to select
 * p: return array with the elements of array
 */
/*
void kLowest(double *array, int n, int k, double *p)
{
    double arr[n];
    for (int i = 0; i < n; ++i)
        arr[i] = array[i];
    // Sort the given array arr in reverse 
    // order.
    //sort(arr, arr+n, greater<int>());
    sort(arr, arr+n, less<int>());
    // Print the first kth largest elements
    for (int i = 0; i < k; i++)
         p[i] = arr[i];
}
*/

void kLowest(std::vector<double> &array, int n, int k, std::vector<double> &p)
{
    std::vector<double> arr(n);
    for (int i = 0; i < n; ++i)
        arr[i] = array[i];
    // Sort the given array arr in reverse 
    // order.
    //sort(arr, arr+n, greater<int>());
    std::sort(arr.begin(), arr.begin()+n, less<int>());
    // Print the first kth largest elements
    for (int i = 0; i < k; i++)
         p[i] = arr[i];
}

double distanceCalculate(double x1, double y1, double x2, double y2)
{
    double x = x1 - x2; //calculating number to square in next step
    double y = y1 - y2;
    double dist;

    dist = pow(x, 2) + pow(y, 2); // calculating Euclidean distance
    dist = sqrt(dist);                  

    return dist;
}

bool contains(int elem, int* individual)
{
    size_t myArraySize = sizeof(individual) / sizeof(int);
    int *end = individual + myArraySize;
    // find the value 0:
    int *result = std::find(individual, end, elem);
    if (result != end)
    {
      // found value at "result" pointer location...
        return true;
    }
    else
    {
        return false;
    }
}

// Search de element in search and return the same index of the vector values
/*
int find_in_array(double elem, double* search, int* values)
{
    int i = 0;
    while(elem != search[i])
    {
        i++;
    }
    return values[i];
}
*/

int find_in_array(double elem, std::vector<double> &search, std::vector<int> &values)
{
    int i = 0;
    while(elem != search[i])
    {
        i++;
    }
    return values[i];
}

std::vector<int> intersection(int arr1[], int arr2[], int m, int n)
{
    std::sort(arr1, arr1 + m);
    std::sort(arr2, arr2 + n);
    std::vector<int> out;
    int i = 0, j = 0;
    while (i < m && j < n)
    {
        if (arr1[i] < arr2[j])
            i++;
        else if (arr2[j] < arr1[i])
            j++;
        else /* if arr1[i] == arr2[j] */
        {
            out.push_back(arr2[j]);
            i++;
            j++;
        }
    }
    return out;
}

std::vector<int> intersection(std::vector<int> &arr1, std::vector<int> &arr2)
{
    std::sort(arr1.begin(), arr1.end());
    std::sort(arr2.begin(), arr2.end());
    std::vector<int> out;
    unsigned i = 0, j = 0;
    while (i < arr1.size() && j < arr2.size())
    {
        if (arr1[i] < arr2[j])
            i++;
        else if (arr2[j] < arr1[i])
            j++;
        else /* if arr1[i] == arr2[j] */
        {
            out.push_back(arr2[j]);
            i++;
            j++;
        }
    }
    return out;
}

std::vector<int> unionarray(int arr1[], int arr2[], int arr1_size, int arr2_size)
{
    std::vector<int> out;
    // Taking max element present in either array
    int m = arr1[arr1_size - 1];
    int n = arr2[arr2_size - 1];
     
    int ans = 0;
     
    if(m > n)
    {
        ans = m;
    }
    else
    ans = n;
     
    // Finding elements from 1st array
    // (non duplicates only). Using 
    // another array for storing union 
    // elements of both arrays
    // Assuming max element present 
    // in array is not more than 10^7
    int *newtable = new int[ans + 1];
     
    // First element is always 
    // present in final answer
    out.push_back(arr1[0]);
     
    // Incrementing the First element's count
    // in it's corresponding index in newtable
    ++newtable[arr1[0]];
     
    // Starting traversing the first 
    // array from 1st index till last
    for(int i = 1; i < arr1_size; i++)
    {
        // Checking whether current element 
        // is not equal to it's previous element
        if(arr1[i] != arr1[i - 1])
        {
            out.push_back(arr1[i]);
            ++newtable[arr1[i]];
        }
    }
     
    // Finding only non common
    // elements from 2nd array        
    for(int j = 0; j < arr2_size; j++)
    {
        // By checking whether it's already 
        // present in newtable or not
        if(newtable[arr2[j]] == 0)
        {
            out.push_back(arr2[j]);
            ++newtable[arr2[j]];
        }
    }
    return out;
}

std::vector<int> unionarray(std::vector<int> &arr1, std::vector<int> &arr2)
{
    std::vector<int> out;

    std::vector<int> v(arr1.size() + arr2.size());                      // 0  0  0  0  0  0  0  0  0  0
    std::vector<int>::iterator it;

    std::sort (arr1.begin(),arr1.end());     //  5 10 15 20 25
    std::sort (arr2.begin(),arr2.end());   // 10 20 30 40 50

    it=std::set_union (arr1.begin(),arr1.end(),arr2.begin(),arr2.end(), v.begin());
                                                   // 5 10 15 20 25 30 40 50  0  0
    v.resize(it-v.begin());                      // 5 10 15 20 25 30 40 50

    for (it=v.begin(); it!=v.end(); ++it)
        out.push_back(*it);

    return out;
}

void minus(std::vector<int> &a, std::vector<int> &b)
{
    sort( begin(b), end(b) );
    a.erase( remove_if( begin(a),end(a),
    [&](int x){return binary_search(begin(b),end(b),x);}), end(a) );
}
