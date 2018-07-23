#ifndef LOCALSEARCHES_H
#define LOCALSEARCHES_H

#include<string>
#include "structs.h"

void local_search(TSolution* x, std::string ls, double param);
void IMP(TSolution* X, double param);
void IALT(TSolution* X, double L);
void ISRW(TSolution* x, int alpha);
void FI(TSolution* x);

#endif