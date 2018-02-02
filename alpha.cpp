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
// particles color array
int *color;
// particles electric charge array
double *q;
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
    double Sx_2 = Sx / 2.0;
    double Sy_2 = Sy / 2.0;

    for (int i = 0; i < N; i++) {
        ID[i] = i;

        int j;
        double tmpX;
        double tmpY;
        do {
            // choose a random position in the system for the particle
            tmpX = Sx * Rand();
            tmpY = Sy * Rand();
            // checking if in this position there already is an element
            for (j = 0; j < i; j++) {
                double diffX = x[j] - tmpX;
                double diffY = y[j] - tmpY;

                if (diffX > Sx_2) diffX -= Sx;
                if (diffX < -Sx_2) diffX += Sx;
                if (diffX > Sy_2) diffX -= Sy;
                if (diffX < -Sy_2) diffX += Sy;

                double diffSquare = diffX * diffX + diffY * diffY;
                if (diffSquare < 0.2 * 0.2) // the distance between two particles has to 
                                            // be bigger then 0.2
                    break;
            }
        // if a managed to iterate over the whole bunch of particles which have their 
        // position that means there wasn't any overlap so the new particle's coordinates
        // are alright 
        // ------
        // ELSE the new element has overlap with one element from the system hence is 
        // needed to regenerate them
        } while (j != i);

        x[i] = tmpX;
        y[i] = tmpY;

        // setting color and charge to particle
        if (Rand() > 0.5) {
            color[i] = 1;
            q[i] = 1;
        }
        else {
            color[i] = 0;
            q[i] = -1;
        }

        // setting the initial forces
        fx[i] = 0;
        fy[i] = 0;
    }
}

int main() {
    return 0;
}