#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <vector>
#include <string>
#include <cmath>

#include "random_numrec.c"

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
// distance cutoff for the interaction 
double r0;
// distance cutoff for the Verlet list 
double r_verlet;
//distance a particle needs to travel
//to potentially destroy the Verlet list
double r_travel_max;


//time
double dt;

//total runtime
int total_runtime;
int time_echo;

//verlet lists
vector<vector<int> > verlet;
int flag_rebuild_verlet;
double *x_so_far;
double *y_so_far;

//time
double dt = 0.0;
time_t beginning_time = 0;
time_t ending_time = 0;

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

void initParticles()
{
    N = 6400;
    ID = new int[N];

    x = new double[N];
    y = new double[N];

    fx = new double[N];
    fy = new double[N];

    color =  new int[N];
    q = new double[N];

    Sx = 80.0;
    Sy = 80.0;

    r0 = 4.0;
    r_verlet = 6.0;
    r_travel_max = r_verlet - r0;
    dt = 0.01;

    flag_rebuild_verlet = 1;
    x_so_far = new double[N];
    y_so_far = new double[N];
}

void freeData()
{
    delete[] ID;
    delete[] x;
    delete[] y;
    delete[] fx;
    delete[] fy;
    delete[] color;
    delete[] q;
    delete[] x_so_far;
    delete[] y_so_far;
}

void generateCoordinates()
{
    double Sx_2 = Sx / 2.0;
    double Sy_2 = Sy / 2.0;

    for (int i = 0; i < N; i++)
    {
        ID[i] = i;

        int j;
        double tmpX;
        double tmpY;
        do {
            // choose a random position in the system for the particle
            tmpX = Sx * Rand();
            tmpY = Sy * Rand();
            // checking if in this position there already is an element
            for (j = 0; j < i; j++)
            {
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
        if (Rand() > 0.5)
        {
            color[i] = 1;
            q[i] = 1.0;
        }
        else
        {
            color[i] = 0;
            q[i] = -1.0;
        }

        // setting the initial forces
        fx[i] = 0.0;
        fy[i] = 0.0;
        x_so_far[i] = 0.0;
        y_so_far[i] = 0.0;
    }
}

void calculateVerletList()
{

    double Sx_2 = Sx / 2.0;
    double Sy_2 = Sy / 2.0;

    for (int i = 0; i < verlet.size(); i++)
        verlet.at(i).clear();
    verlet.clear();

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            double diffX = x[i] - x[j];
            double diffY = y[i] - y[j];

            if (diffX < -Sx_2) diffX += Sx;
            if (diffX > Sx_2) diffX -= Sx;
            if (diffY < Sy_2) diffY += Sy;
            if (diffY > Sy_2) diffY -= Sy;

            double distance = diffX * diffX + diffY * diffY;

            if (distance < r_verlet * r_verlet)
            {
                vector<int> v;
                v.push_back(i);
                v.push_back(j);
                verlet.push_back(v);
            }
        }
        
    }

    //clear all accumulated distances
    for (int i = 0; i < N; i++)
    {
        x_so_far[i] = 0.0;
        y_so_far[i] = 0.0;
    }
    //clear rebuild flag
    flag_rebuild_verlet = 0;

    printf(".");fflush(stdout);
}

void colorverlet()
{
    for(int i = 0; i < N; i++)
    {
        if (q[i] == -1) color[i] = 0;        
        if (q[i] == 1) color[i] = 1;
    }

    for(int i = 0; i < verlet.size(); i++)
    {
        if (verlet.at(i).at(0) == 30) color[verlet.at(i).at(1)] = 2;
        if (verlet.at(i).at(1) == 30) color[verlet.at(i).at(0)] = 2;
    }
}


void writeToFile(char* filename)
{
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

    for (int it = 0; it < verlet.size(); it++)
    {

        int i = verlet.at(it).at(0);
        int j = verlet.at(it).at(1);
        
        double diffX = x[i] - x[j];
        double diffY = y[i] - y[j];

        if (diffX < -Sx_2) diffX += Sx;
        if (diffX > Sx_2) diffX -=Sx;
        if (diffY < -Sy_2) diffY += Sy;
        if (diffY > Sy_2) diffY -= Sy;

        double distance = diffX * diffX + diffY * diffY;
        double dist_sqrt = sqrt(distance);

        double f = exp(- dist_sqrt / r0) / distance;

        fx[i] += f * diffX / dist_sqrt;
        fy[i] += f * diffY / dist_sqrt;
        fx[j] -= f * diffX / dist_sqrt;
        fy[j] -= f * diffY / dist_sqrt;
    }
}

void calculateExternalForces()
{
    for (int i = 0; i < N; i++)
    {
        fx[i] += q[i] * 0.5;
    }
}

void moveParticles()
{
    double deltax, deltay;

    for (int i = 0; i < N; i++) 
    {
        deltax = fx[i] * dt;
        deltay = fy[i] * dt;

        x[i] += deltax;
        y[i] += deltay;

        x_so_far[i] += deltax;
        y_so_far[i] += deltay;

        if (x[i] < 0) x[i] += Sx;
        if (x[i] > Sx) x[i] -= Sx;
        if (y[i] < 0) y[i] += Sy;
        if (y[i] > Sy) y[i] -=Sy;

        //only check on unset flag
        if (flag_rebuild_verlet == 0)
            if (x_so_far[i] * x_so_far[i] + y_so_far[i] * y_so_far[i] >= r_travel_max * r_travel_max)
                flag_rebuild_verlet = 1;

        fx[i] = 0.0;
        fy[i] = 0.0;
    }
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

int main(int argc, char* argv[]) 
{
    printf("Generic BD simulation\n");
    printf("Simulating spontaneous lane formation\n");

    initParticles();
    generateCoordinates();
    const char* moviefile = new char[20];
    FILE* f;
    if (argc == 2) 
        moviefile = argv[1];
    else
        moviefile = "moviefile.mvi";
    f = fopen(moviefile, "wb");

    total_runtime = 20000;
    time_echo = 500;
    start_timing();
    for (int i = 0; i < total_runtime; i++) 
    {
        if (i % time_echo == 0)
        {
            double perc = (double) i / (double) total_runtime * 100;
            cout << i << "/" << total_runtime << " " << perc << "%" << endl;
        }

        if (flag_rebuild_verlet == 1) 
        {
            calculateVerletList();
            colorverlet();
        }

        calculateForces();
        calculateExternalForces();
        moveParticles();

        if (i % 100 == 0)
            write_cmovie(f, i);
    }
    freeData();
    start_timing();
    return 0;
}