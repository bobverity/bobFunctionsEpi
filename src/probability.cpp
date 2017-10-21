
#include <Rcpp.h>
#include <random>
#include "probability.h"
#include "misc.h"

using namespace std;

//-- set random seed --
random_device rd;
default_random_engine generator(rd());

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
// NB: execution time found to be ~69% of R::runif() method
double runif_0_1() {
    uniform_real_distribution<double> uniform_0_1(0.0,1.0);
    return uniform_0_1(generator);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
// NB: execution time found to be ~69% of R::runif() method
double runif1(double a, double b) {
    uniform_real_distribution<double> uniform_a_b(a,b);
    return(uniform_a_b(generator));
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
// NB: execution time found to be ~34% of R::rbinom() method
bool rbernoulli1(double p) {
    bernoulli_distribution dist_bernoulli(p);
    return dist_bernoulli(generator);
}

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(vector<double> &p, double pSum) {
    double rand = pSum*runif_0_1();
    double z = 0;
    for (int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z) {
            return i+1;
        }
    }
    return 0;
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability. Works on positive or negative values of a or b, and works irrespective of which of a or b is larger.
int sample2(int a, int b) {
    if (a<b) {
        return floor(runif1(a, b+1));
    } else {
        return floor(runif1(b, a+1));
    }
}

//------------------------------------------------
// sample without replacement from vector x
vector<int> sample3(vector<int> x, int n) {
    int rand1, tmp1;
    for (int i=0; i<n; i++) {
        rand1 = sample2(i,x.size()-1);
        tmp1 = x[rand1];
        x[rand1] = x[i];
        x[i] = tmp1;
    }
    vector<int> ret(x.begin(), x.begin()+n);
    return ret;
}

//------------------------------------------------
// draw from binomial(n,p) distribution
// NB: execution time found to be ~64% of c++ binomial_distribution method
int rbinom1(int n, double p) {
    return R::rbinom(n,p);
}

//------------------------------------------------
// draw from Poisson(rate) distribution
// NB: execution time found to be generally equal or faster than c++ poisson_distribution method (e.g. faster when rate~=5)
int rpois1B(double rate) {
    return R::rpois(rate);
}

//------------------------------------------------
// draw from exponential(rate) distribution (expectation = 1/rate)
// NB: execution time found to be ~25% of R::rexp() method
double rexp1(double rate) {
    exponential_distribution<double> dist_exponential(rate);
    return dist_exponential(generator);
}

//------------------------------------------------
// draw from gamma(shape,rate) distribution
// NB: execution time found to be ~22% of R::rgamma() method
double rgamma1(double shape, double rate) {
    gamma_distribution<double> rgamma(shape,1.0/rate);
    double x = rgamma(generator);
    
    // check for zero or infinite values (catches bug present in Visual Studio 2010)
    if (x<UNDERFLO) {
        x = UNDERFLO;
    }
    if (x>OVERFLO) {
        x = OVERFLO;
    }
    
    return x;
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
// NB: execution time found to be ~50% of cpp ratio of gamma draws method
double rbeta1(double alpha, double beta) {
    if (alpha==1 && beta==1) {
        return runif_0_1();
    }
    return R::rbeta(alpha,beta);
}

//------------------------------------------------
// draw from univariate normal distribution
// NB: execution time found to be ~78% of R::rnorm() method
double rnorm1(double mean, double sd) {
    normal_distribution<double> normal_dist(mean,sd);
    return normal_dist(generator);
}

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {
    
    // draw raw value relative to a
    double ret = rnorm1(mean, sd) - a;
    
    // reflect off boundries at 0 and (b-a)
    if (ret<0 || ret>(b-a)) {
        // use multiple reflections to bring into range [-(b-a), 2(b-a)]
        while (ret < -(b-a)) {
            ret += 2*(b-a);
        }
        while (ret > 2*(b-a)) {
            ret -= 2*(b-a);
        }
        
        // use one more reflection to bring into range [0, (b-a)]
        if (ret < 0) {
            ret = -ret;
        }
        if (ret > (b-a)) {
            ret = 2*(b-a) - ret;
        }
    }
    
    // no longer relative to a
    ret += a;
    
    // don't let ret equal exactly a or b
    if (ret==a) {
        ret += UNDERFLO;
    } else if (ret==b) {
        ret -= UNDERFLO;
    }
    
    return ret;
}

//------------------------------------------------
// draw from logit-normal distribution, given mean and sd of underlying normal random variable. If hardLimits==true then cannot produce exactly 0 or 1.
double rlogitnorm1(double mean, double sd, bool hardLimits) {
    double x = rnorm1(mean,sd);
    double ret = 1/(1+exp(-x));
    
    // apply hard limits
    if (hardLimits) {
        if (ret<UNDERFLO) {
            ret = UNDERFLO;
        }
        if (ret>(1.0-UNDERFLO)) {
            ret = 1.0-UNDERFLO;
        }
    }
    
    return ret;
}

//------------------------------------------------
// draw from logit-normal distribution, centred on p and given sd of underlying normal random variable (useful in e.g. Metropolis-Hastings). If hardLimits==true then cannot produce exactly 0 or 1.
double rlogitnorm2(double p, double sd, bool hardLimits) {
    double x = rnorm1(log(p/(1-p)) ,sd);
    double ret = 1/(1+exp(-x));
    
    // apply hard limits
    if (hardLimits) {
        if (ret<UNDERFLO) {
            ret = UNDERFLO;
        }
        if (ret>(1.0-UNDERFLO)) {
            ret = 1.0-UNDERFLO;
        }
    }
    
    return ret;
}

//------------------------------------------------
// time evaluation of code blocks
// [[Rcpp::export]]
void functionTimer_cpp(int reps) {
    
    // timer objects
    chrono::high_resolution_clock::time_point t1, t2;
    chrono::duration<double> diff1, diff2;
    
    // other variables
    double x;
    
    // first code block
    t1 = chrono::high_resolution_clock::now();
    for (int i=0; i<reps; i++) {
        
        x = rlogitnorm2(0.5, 1);
        
    }
    t2 = chrono::high_resolution_clock::now();
    diff1 = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print(diff1.count(), "seconds");
    
    
    // second code block
    t1 = chrono::high_resolution_clock::now();
    for (int i=0; i<reps; i++) {
        
        x = rlogitnorm1(0,1);
    }
    t2 = chrono::high_resolution_clock::now();
    diff2 = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print(diff2.count(), "seconds");
    
    // relative speed
    if (diff1<diff2) {
        print("first method faster: execution time", diff1/diff2*100, "% of second method");
    } else {
        print("second method faster: execution time", diff2/diff1*100, "% of first method");
    }
    
}


