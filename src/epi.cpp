
#include <Rcpp.h>
#include "misc.h"
#include "probability.h"
#include "epi.h"

using namespace std;

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    int I_start = Rcpp_to_int(args["I_start"]);
    int N = Rcpp_to_int(args["N"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> I_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    I_vec[0] = I_start;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*I*(1-I*N_inv);
        rate2 = r*I;
        rateTotal = rate1+rate2;
        
        // abort if reached stable state
        if (rateTotal==0) {
            break;
        }
        
        // draw new time
        t += rexp1(rateTotal);
        
        // draw event
        rand1 = runif_0_1();
        if (rand1<(rate1/rateTotal)) {
            I++;
        } else {
            I--;
        }
        
        // store values
        t_vec[i] = t;
        I_vec[i] = I;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec, Rcpp::Named("I")=I_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    int I_start = Rcpp_to_int(args["I_start"]);
    int N = Rcpp_to_int(args["N"]);
    vector<double> t_vec = Rcpp_to_vector_double(args["t_vec"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> I_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    I_vec[0] = I;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    int j=0;
    for (int i=0; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*I*(1-I*N_inv);
        rate2 = r*I;
        rateTotal = rate1+rate2;
        
        // draw new time
        t += rexp1(rateTotal);
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size) {
                break;
            }
            I_vec[j] = I;
            j++;
        }
        
        // draw event
        rand1 = unif_rand();
        if (rand1<(rate1/rateTotal)) {
            I++;
        } else {
            I--;
        }
        
        // abort on end of t_vec
        if (j==t_size) {
            break;
        }
        
        // report if maxIterations reached
        if (i==(maxIterations-1)) {
            print("maxIterations reached\n");
        }
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec, Rcpp::Named("I")=I_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    int I_start = Rcpp_to_int(args["I_start"]);
    int N = Rcpp_to_int(args["N"]);
    vector<double> t_vec = Rcpp_to_vector_double(args["t_vec"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> I_vec(t_size);
    I_vec[0] = I;
    
    // carry out simulation loop
    double rate1, rate2, prob1, prob2, delta_t;
    int rand1, rand2;
    for (int i=1; i<t_size; i++) {
        
        // calculate rates of all events
        rate1 = beta*I*N_inv;
        rate2 = r;
        
        // convert to probabilities
        delta_t = t_vec[i]-t_vec[i-1];
        prob1 = 1 - exp(-rate1*delta_t);
        prob2 = 1 - exp(-rate2*delta_t);
        
        // draw events
        rand1 = rbinom1(N-I, prob1); // new infections
        rand2 = rbinom1(I, prob2); // recoveries
        I += rand1 - rand2;
        
        // store values
        I_vec[i] = I;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec, Rcpp::Named("I")=I_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    double mu = Rcpp_to_double(args["mu"]);
    int I_init = Rcpp_to_int(args["I_init"]);
    int R_init = Rcpp_to_int(args["R_init"]);
    int N = Rcpp_to_int(args["N"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> S_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> I_vec(maxIterations, -1);
    vector<double> R_vec(maxIterations, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rateTotal;
    vector<double> rates(5);
    int rand1;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rates[0] = beta*S*I*N_inv; // infection
        rates[1] = r*I; // recovery
        rates[2] = mu*S; // natural death in S
        rates[3] = mu*I; // natural death in I
        rates[4] = mu*R; // natural death in R
        rateTotal = sum(rates);
        
        // abort if reached stable state
        if (rateTotal==0) {
            break;
        }
        
        // draw new time
        t += rexp1(rateTotal);
        t_vec[i] = t;
        
        // draw event
        rand1 = sample1(rates, rateTotal);
        switch(rand1) {
            case 1: // infection
                S --;
                I ++;
                break;
            case 2: // recovery
                I --;
                R ++;
                break;
            case 3: // natural death in S
                // (death exactly matched by new birth)
                break;
            case 4: // natural death in I
                I --;
                S ++;
                break;
            case 5: // natural death in R
                R --;
                S ++;
        }
        
        // store values
        S_vec[i] = S;
        I_vec[i] = I;
        R_vec[i] = R;
        
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    double mu = Rcpp_to_double(args["mu"]);
    int I_init = Rcpp_to_int(args["I_init"]);
    int R_init = Rcpp_to_int(args["R_init"]);
    int N = Rcpp_to_int(args["N"]);
    vector<double> t_vec = Rcpp_to_vector_double(args["t_vec"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rateTotal;
    vector<double> rates(5);
    int rand1, j=0;
    for (int i=0; i<maxIterations; i++) {
        
        // calculate rates of all events
        rates[0] = beta*S*I*N_inv; // infection
        rates[1] = r*I; // recovery
        rates[2] = mu*S; // natural death in S
        rates[3] = mu*I; // natural death in I
        rates[4] = mu*R; // natural death in R
        rateTotal = sum(rates);
        
        // draw new time
        t += rexp1(rateTotal);
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size) {
                break;
            }
            S_vec[j] = S;
            I_vec[j] = I;
            R_vec[j] = R;
            j++;
        }
        
        // draw event
        rand1 = sample1(rates, rateTotal);
        switch(rand1) {
            case 1: // infection
                S --;
                I ++;
                break;
            case 2: // recovery
                I --;
                R ++;
                break;
            case 3: // natural death in S
                // (death exactly matched by new birth)
                break;
            case 4: // natural death in I
                I --;
                S ++;
                break;
            case 5: // natural death in R
                R --;
                S ++;
        }
        
        // abort on end of t_vec
        if (j==t_size) {
            break;
        }
        
        // report if maxIterations reached
        if (i==(maxIterations-1)) {
            print("maxIterations reached");
        }
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double r = Rcpp_to_double(args["r"]);
    double mu = Rcpp_to_double(args["mu"]);
    int I_init = Rcpp_to_int(args["I_init"]);
    int R_init = Rcpp_to_int(args["R_init"]);
    int N = Rcpp_to_int(args["N"]);
    vector<double> t_vec = Rcpp_to_vector_double(args["t_vec"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1);
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double delta_t, rate_inf;
    vector<double> probs(5);
    vector<int> rands(6);
    for (int i=1; i<t_size; i++) {
        
        // calculate rates of all events
        rate_inf = beta*I*N_inv; // infection
        
        // convert to probabilities, allowing for competing hazards
        delta_t = t_vec[i]-t_vec[i-1];
        if (rate_inf==0 && mu==0) {
            probs[0] = 0;
            probs[1] = 0;
        } else {
            probs[0] = 1 - exp(-(rate_inf + mu)*delta_t); // infection or death in S
            probs[1] = rate_inf/(rate_inf+mu); // infection
        }
        if (r==0 && mu==0) {
            probs[2] = 0;
            probs[3] = 0;
        } else {
            probs[2] = 1 - exp(-(r + mu)*delta_t); // recovery or death in I
            probs[3] = r/(r+mu); // recovery
        }
        probs[4] = 1 - exp(-mu*delta_t); // death in R
        
        // draw events
        rands[0] = rbinom1(S, probs[0]); // infection or death in S
        rands[1] = rbinom1(rands[0], probs[1]); // infection
        rands[2] = rbinom1(I, probs[2]); // recovery or death in I
        rands[3] = rbinom1(rands[2], probs[3]); // recovery
        rands[4] = rands[2] - rands[3]; // death in I
        rands[5] = rbinom1(R, probs[4]); // death in R
        
        S += -rands[1] + rands[4] + rands[5];
        I += rands[1] - rands[3] - rands[4];
        R += rands[3] - rands[5];
        
        // store values
        S_vec[i] = S;
        I_vec[i] = I;
        R_vec[i] = R;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SLIR_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double dur_lag = Rcpp_to_double(args["dur_lag"]);
    double r = Rcpp_to_double(args["r"]);
    int I_init = Rcpp_to_int(args["I_init"]);
    int R_init = Rcpp_to_int(args["R_init"]);
    int N = Rcpp_to_int(args["N"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int L = 0;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> S_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> L_vec(maxIterations, -1);
    vector<double> I_vec(maxIterations, -1);
    vector<double> R_vec(maxIterations, -1);
    S_vec[0] = S;
    L_vec[0] = 0;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    vector<double> lagList;
    int lagList_size = 0;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*S*I*N_inv; // infection (move to lag state)
        rate2 = r*I; // recovery
        rateTotal = rate1 + rate2;
        
        // abort if reached stable state
        if (rateTotal==0) {
            break;
        }
        
        // draw new time
        rand1 = rexp1(rateTotal);
        if (lagList_size>0) {
            if ((t+rand1)<lagList[0]) {
                t += rand1;
            } else {
                t = lagList[0];
                lagList.erase (lagList.begin());
                lagList_size --;
                L--;
                I++;
                goto store;
            }
        } else {
            t += rand1;
        }
        
        // draw event
        rand1 = runif1();
        if (rand1<(rate1/rateTotal)) {
            S--;
            L++;
            lagList.push_back(t + dur_lag);
            lagList_size ++;
        } else {
            I--;
            R++;
        }
        
        // store values
        store:
        t_vec[i] = t;
        S_vec[i] = S;
        L_vec[i] = L;
        I_vec[i] = I;
        R_vec[i] = R;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("L")=L_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
}

//------------------------------------------------
// See R function of the same name (without _cpp) for details
// [[Rcpp::export]]
Rcpp::List SLIR_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp_to_double(args["beta"]);
    double dur_lag = Rcpp_to_double(args["dur_lag"]);
    double r = Rcpp_to_double(args["r"]);
    int I_init = Rcpp_to_int(args["I_init"]);
    int R_init = Rcpp_to_int(args["R_init"]);
    int N = Rcpp_to_int(args["N"]);
    vector<double> t_vec = Rcpp_to_vector_double(args["t_vec"]);
    int maxIterations = Rcpp_to_int(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int L = 0;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> L_vec(t_size, -1);
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    L_vec[0] = L;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    vector<double> lagList;
    int lagList_size=0, j=0;
    bool new_event;
    for (int i=0; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*S*I*N_inv; // infection (move to lag state)
        rate2 = r*I; // recovery
        rateTotal = rate1 + rate2;
        
        // draw new time
        rand1 = rexp1(rateTotal);
        new_event = true;
        if (lagList_size>0) {
            if ((t+rand1)<lagList[0]) {
                t += rand1;
            } else {
                t = lagList[0];
                lagList.erase (lagList.begin());
                lagList_size --;
                new_event = false;
            }
        } else {
            t += rand1;
        }
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size) {
                break;
            }
            S_vec[j] = S;
            L_vec[j] = L;
            I_vec[j] = I;
            R_vec[j] = R;
            j++;
        }
        
        // draw event
        if (new_event) {
            rand1 = runif1();
            if (rand1<(rate1/rateTotal)) {
                S--;
                L++;
                lagList.push_back(t + dur_lag);
                lagList_size ++;
            } else {
                I--;
                R++;
            }
        } else {
            L--;
            I++;
        }
        
        // abort on end of t_vec
        if (j==t_size) {
            break;
        }
        
        // report if maxIterations reached
        if (i==(maxIterations-1)) {
            print("maxIterations reached");
        }
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("L")=L_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
}

