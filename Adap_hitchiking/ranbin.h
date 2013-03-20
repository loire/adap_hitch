/*
 * ranbin.h
 *
 *  Created on: 11 févr. 2013
 *      Author: etienne
 */
#ifndef RANBIN_H
#define RANBIN_H

#include <cmath>
#include "MersenneTwister.h"

using namespace std;

extern MTRand rnd;

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gasdev();

#endif
