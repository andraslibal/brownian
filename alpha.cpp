#include <iostream>
#include <cstdlib>
#include "random_numrec.c"

using namespace std;

// number of elements
int N;
// ids of the particles
int *ID;
// coordinate arrays
double *x;
double *y;
// forces acting on a paricle
double *fx;
double *fy;
// particle's color array
int *color;
// system size
double Sx, Sy;

int initParticles() {
    N = 2000;
    ID = new int(N);

    x = new double(N);
    y = new double(N);

    fx = new double(N);
    fy = new double(N);

    color =  new int(N);

    Sx = 20.0;
    Sy = 20.0;
}

int generateCoordinates() {
    for (int i = 0; i < N; i++) {
        ID[i] = i;
        // choose a random position in the system for the particle
        double tmpX = Sx * Rand();
        double tmpY = Sy * Rand();
        // checking if in this position there already is an element
        for (int j = 0; j < i; j++) {
            double diffX = x[j] - tmpX;
            double diffY = y[j] - tmpY;

            
        }

    }
}

int main() {
    
    return 0;
}