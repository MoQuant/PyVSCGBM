#include <stdlib.h>
#include <math.h>
#include <time.h>

double dWT()
{
    int num = 35;
    double dw = (rand() % (2*num + 1)) - num;
    return dw / 10.0;
}

double simulation(double S, double drift, double vol, double t, int steps, int paths){
    srand(time(NULL));
    double total = 0.0;
    double dt = t / (double) steps;
    for(int i = 0; i < paths; ++i){
        double S0 = S;
        for(int j = 0; j < steps; ++j){
            S0 += drift*S0*dt + vol*S0*dWT();
        }
        total += S0;
    }
    return total / (double) paths;
}