#ifndef LOCALSEARCHES_H
#define LOCALSEARCHES_H

#include<string>
#include "structs.h"

void local_search(TSolution* x, std::string ls, double param, int max_time = std::numeric_limits<int>::max());
void IMP(TSolution* X, double param, int max_time = std::numeric_limits<int>::max());
void IALT(TSolution* X, double L, int max_time = std::numeric_limits<int>::max());
void ISRW(TSolution* x, int alpha, int max_time = std::numeric_limits<int>::max());
void FI(TSolution* x, int max_time = std::numeric_limits<int>::max());

#endif