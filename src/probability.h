
#pragma once

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
// NB: execution time found to be ~69% of R::runif() method
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
// NB: execution time found to be ~69% of R::runif() method
double runif1(double a=0, double b=1.0);

//------------------------------------------------
// draw from Bernoulli(p) distribution
// NB: execution time found to be ~34% of R::rbinom() method
bool rbernoulli1(double p=0.5);

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum=1);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability. Works on positive or negative values of a or b, and works irrespective of which of a or b is larger.
int sample2(int b, int a=1);

//------------------------------------------------
// sample without replacement from vector x
std::vector<int> sample3(std::vector<int> x, int n);

//------------------------------------------------
// draw from binomial(n,p) distribution
// NB: execution time found to be ~64% of c++ binomial_distribution method
int rbinom1(int n, double p);

//------------------------------------------------
// draw from Poisson(rate) distribution
// NB: execution time found to be generally equal or faster than c++ poisson_distribution method (e.g. faster when rate~=5)
int rpois1(double rate);

//------------------------------------------------
// draw from exponential(rate) distribution (expectation = 1/rate)
// NB: execution time found to be ~25% of R::rexp() method
double rexp1(double rate);

//------------------------------------------------
// draw from gamma(shape,rate) distribution
// NB: execution time found to be ~22% of R::rgamma() method
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(alpha,beta) distribution
// NB: execution time found to be ~50% of cpp ratio of gamma draws method
double rbeta1(double alpha, double beta);

//------------------------------------------------
// draw from univariate normal distribution
// NB: execution time found to be ~78% of R::rnorm() method
double rnorm1(double mean=0, double sd=1);

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
// draw from logit-normal distribution, given mean and sd of underlying normal random variable. If hardLimits==true then cannot produce exactly 0 or 1.
double rlogitnorm1(double mean, double sd, bool hardLimits=true);

//------------------------------------------------
// draw from logit-normal distribution, centred on p and given sd of underlying normal random variable (useful in e.g. Metropolis-Hastings). If hardLimits==true then cannot produce exactly 0 or 1.
double rlogitnorm2(double p, double sd, bool hardLimits=true);

//------------------------------------------------
// time evaluation of code blocks
void functionTimer_cpp(int reps);
