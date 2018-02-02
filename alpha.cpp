#include <iostream>

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

    for(int i = 0; i < N; i++) {

    }
}

int main() {
    cout << "Helloka" << endl;
}