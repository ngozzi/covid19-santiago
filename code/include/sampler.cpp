//
// Created by Staff on 22/06/2020.
//

#include "sampler.h"
random_device device;
mt19937 g(device());

/**
 * Binomial distribution sampler
 */
double binomial(int n, double p)
{
    binomial_distribution<int> distr(n, p);
    return distr(g);
}

/**
 * Uniform random sampler (0, 1)
 */
double rnd()
{
    uniform_real_distribution<double> distr(0.0, 1.0);
    return distr(g);
}




