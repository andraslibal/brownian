#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include "random_numrec.c"
#include <cmath>
#include <string.h>

using namespace std;

// number of elements
int N;
// ids of the particles
int *ID = NULL;
// coordinate arrays
double *x = NULL;
double *y = NULL;
// forces acting on a paricle
double *fx = NULL;
double *fy = NULL;
// particles color array
int *color = NULL;
// particles electric charge array
double *q = NULL;
// system size
double Sx, Sy;
double Sx_2, Sy_2;
// distance 
double r;
//time
double dt;
time_t beginning_time = 0;
time_t ending_time = 0;

char* moviefile;

void start_timing()
{
    time(&beginning_time);
}

void stop_timing()
{
    time(&ending_time);
    double difference = difftime(ending_time, beginning_time);
    tm* time_info = localtime(&beginning_time);
    cout << "Program started at:" << asctime(time_info) << endl;

    time_info = localtime(&ending_time);
    cout << "Program ended at:" << asctime(time_info) << endl;

    cout << "Program running time was " << difference << " seconds" << endl;
}

void readDataFromFile(const char* filename)
{
    ifstream f(filename);
    if (f.is_open()) {
        f >> moviefile;
        strcat(moviefile, ".mvi");
        f >> Sx;
        Sx_2 = Sx / 2;
        f >> Sy;
        Sy_2 = Sy / 2;
        f >> N;
    }

    f.close();
}

void initParticles() {
    N = 51200;
    ID = new int[N];

    x = new double[N];
    y = new double[N];

    fx = new double[N];
    fy = new double[N];

    color =  new int[N];
    q = new double[N];

    moviefile = new char[30];
    strcpy(moviefile, "first.mvi");

    Sx = 320.0;
    Sy = 320.0;

    r = 4.0;
    dt = 0.01;
}

void freeData() {
    delete[] ID;
    delete[] x;
    delete[] y;
    delete[] fx;
    delete[] fy;
    delete[] color;
    delete[] q;
    delete[] moviefile;
}

void generateCoordinates() {
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
            q[i] = 1.0;
        }
        else {
            color[i] = 0;
            q[i] = -1.0;
        }

        // setting the initial forces
        fx[i] = 0.0;
        fy[i] = 0.0;
    }
}

void writeToFile(char* filename) {
    ofstream f(filename);
    cout << filename << endl;
    f << setw(20) << "ID";
    f << setw(20) << "X";
    f << setw(20) << "Y";
    f << setw(20) << "FX";
    f << setw(20) << "FY";
    f << setw(20) << "COLOR";
    f << setw(20) << "Q" << endl;
    for(int i = 0; i < N; i++) {
        f << setw(20) << ID[i];
        f << setw(20) << x[i];
        f << setw(20) << y[i];
        f << setw(20) << fx[i];
        f << setw(20) << fy[i];
        f << setw(20) << color[i];
        f << setw(20) << q[i] << endl;
    }
    f.close();
}

void calculateForces() {
    double Sx_2 = Sx / 2.0;
    double Sy_2 = Sy / 2.0;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double diffX = x[i] - x[j];
            double diffY = y[i] - y[j];

            if (diffX < -Sx_2) diffX += Sx;
            if (diffX > Sx_2) diffX -=Sx;
            if (diffY < -Sy_2) diffY += Sy;
            if (diffY > Sy_2) diffY -= Sy;

            double distance = diffX * diffX + diffY * diffY;

            double f = 1 / (distance * distance) * exp(- distance / r);

            double dist_sqrt = sqrt(distance);

            fx[i] += f * diffX / dist_sqrt;
            fy[i] += f * diffY / dist_sqrt;
            fx[j] -= f * diffX / dist_sqrt;
            fy[j] -= f * diffY / dist_sqrt;
        }
    }
}

void calculateExternalForces() {
    for (int i = 0; i < N; i++) {
        fx[i] += q[i] * 0.5;
    }
}

void moveParticles() {
    for (int i = 0; i < N; i++) {
        x[i] += fx[i] * dt;
        y[i] += fy[i] * dt;

        if (x[i] < 0) x[i] += Sx;
        if (x[i] > Sx) x[i] -= Sx;
        if (y[i] < 0) y[i] += Sy;
        if (y[i] > Sy) y[i] -=Sy;

        fx[i] = 0.0;
        fy[i] = 0.0;
    }
}

void writeCmovie(char* file, int t)
{
    int i;
    float floatholder;
    int intholder;

    ofstream f(file, ios_base::app | ios_base::out | ios::binary);
    
    f << N;
    f << t;
    
    for (int i = 0; i < N; i++)
    {
        f << color[i] + 2;
        f << ID[i];
        f << x[i];
        f << y[i];
        f << 1.0;
    }
    f.close();
} 

void write_cmovie(FILE* moviefile, int t)
{
    int i;
    float floatholder;
    int intholder;
    
    intholder = N;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    intholder = t;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    for (int i = 0; i < N; i++)
    {
        intholder = color[i] + 2;
        fwrite(&intholder,sizeof(int),1,moviefile);
        intholder = ID[i];//ID
        fwrite(&intholder,sizeof(int),1,moviefile);
        floatholder = (float)x[i];
        fwrite(&floatholder,sizeof(float),1, moviefile);
        floatholder = (float)y[i];
        fwrite(&floatholder,sizeof(float),1,moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder,sizeof(float),1,moviefile);
    }
    
}

int main(int argc, char* argv[]) {
    initParticles();
    FILE* f;
    if (argc == 2) 
        readDataFromFile(argv[1]);
    f = fopen(moviefile, "wb");
    generateCoordinates();
    cout << "particles generated" << endl;
    start_timing();
    int time_echo = 100, total_time = 20000;
    for (int i = 0; i < total_time; i++) {
        calculateForces();
        calculateExternalForces();
        moveParticles();

        if (i % 100 == 0)
            write_cmovie(f, i);

        if (i % time_echo == 0)
		{
			double perc = (double) i / (double) total_time * 100;
			cout << i << "/" << total_time << " " << perc << "%" << endl;
		}
    }
    freeData();
    stop_timing();
    return 0;
}