#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>

extern bool DEBUG;
extern bool SAVING;


void log(const std::string& text, bool break_line = true);
void log(const char*& text, bool break_line = true);
//void kLargest(double *array, int n, int k, double *p);
void kLargest(std::vector<double> &array, int n, int k, std::vector<double> &p);
//void kLowest(double *array, int n, int k, double *p);
void kLowest(std::vector<double> &array, int n, int k, std::vector<double> &p);
double distanceCalculate(double x1, double y1, double x2, double y2);
bool contains(int elem, int* individual);
//int find_in_array(int elem, int* search, int* values);
//int find_in_array(double elem, double* search, int* values);
int find_in_array(double elem, std::vector<double> &search, std::vector<int> &values);
std::vector<int> intersection(int arr1[], int arr2[], int m, int n);
std::vector<int> intersection(std::vector<int> &arr1, std::vector<int> &arr2);
std::vector<int> unionarray(int arr1[], int arr2[], int arr1_size, int arr2_size);
std::vector<int> unionarray(std::vector<int> &arr1, std::vector<int> &arr2);
void minus(std::vector<int> &a, std::vector<int> &b);

#endif
