#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <vector>
#include <string.h>
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
double Sx_2, Sy_2;
// distance cutoff for the interaction 
double r0;
// distance cutoff for the Verlet list 
double r_verlet;
//distance a particle needs to travel
//to potentially destroy the Verlet list
double r_travel_max;

//verlet lists
vector<vector<int> > verlet;
int rebuild_verlet_flag;
double *x_so_far;
double *y_so_far;

//f in each point depending on r
double *tabulated_force;
//measure of tabulation
int N_tabulate;
double tab_start;
double tab_measure;

// number of pinning sites
int N_pinning_sites;
double r_pinning_site;
double f_max_pinning_sites;
// coordinates of pinning sites
double *x_pinning_site;
double *y_pinning_site;

char* moviefile;

//time
double dt = 0.0;
time_t beginning_time = 0;
time_t ending_time = 0;

//total runtime
int time_echo;
int current_time = 0;
int total_time = 0;

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

void initParticles()
{
    N = 51200;
    ID = new int[N];

    moviefile = new char[30];

    x = new double[N];
    y = new double[N];

    fx = new double[N];
    fy = new double[N];

    color =  new int[N];
    q = new double[N];

    N_tabulate = 50000;
    tabulated_force = new double[N_tabulate];

    Sx = 320.0;
    Sx_2 = Sx / 2;
    Sy = 320.0;
    Sy_2 = Sy / 2;

    r0 = 4.0;
    r_verlet = 6.0;
    r_travel_max = r_verlet - r0;

    dt = 0.01;
    current_time = 0;
    total_time = 100000;
    time_echo = 500;

    rebuild_verlet_flag = 1;
    x_so_far = new double[N];
    y_so_far = new double[N];

    N_pinning_sites = 50;
    r_pinning_site = 1.0;
    f_max_pinning_sites = 2.0;
    x_pinning_site = new double[N_pinning_sites];
    y_pinning_site = new double[N_pinning_sites];
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
        // if it managed to iterate over the whole bunch of particles which have their 
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

void initPinningSites()
{
	for (int i = 0; i < N_pinning_sites; i++)
	{
		int j;
        double tmpX;
        double tmpY;
        do {
            // choose a random position in the system for the pinning site
            tmpX = Sx * Rand();
            tmpY = Sy * Rand();
            // checking if in this position there already is an element
            for (j = 0; j < i; j++)
            {
                double diffX = x_pinning_site[j] - tmpX;
                double diffY = y_pinning_site[j] - tmpY;

                if (diffX > Sx_2) diffX -= Sx;
                if (diffX < -Sx_2) diffX += Sx;
                if (diffX > Sy_2) diffX -= Sy;
                if (diffX < -Sy_2) diffX += Sy;

                double diffSquare = diffX * diffX + diffY * diffY;
                if (diffSquare < 2.2 * 2.2) // the distance between two particles has to 
                                            // be bigger then 2.2
                    break;
            }
        // if it managed to iterate over the whole bunch of particles which have their 
        // position that means there wasn't any overlap so the new particle's coordinates
        // are alright 
        // ------
        // ELSE the new element has overlap with one element from the system hence is 
        // needed to regenerate them
        } while (j != i);

        x_pinning_site[i] = tmpX;
        y_pinning_site[i] = tmpY;
	}
}

void calculateTabulatedForces() 
{
    double r = 0.1;
    double r2 = r * r;
    tab_measure = (r_verlet * r_verlet - r2) / (N_tabulate - 1.0);
    tab_start = r;
    for (int i = 0; i < N_tabulate; i++) {
        tabulated_force[i] = exp(- r / r0) / r2 * r;
        r2 += tab_measure;
        r = sqrt(r2);
    }
}

void calculateVerletList()
{
	// clearing previous data from verlet list
    for (int i = 0; i < verlet.size(); i++)
        verlet.at(i).clear();
    verlet.clear();

    // iterating through the particles
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

            // if they are closer than a specified value, that means that they are in interaction so I put the in the verlet list
            // one pair will be just once in that list
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
    rebuild_verlet_flag = 0;
}

void colorVerlet()
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

void calculateForces()
{
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
        int tab_index = (int)floor((distance-tab_start)/(tab_measure));
        double f = tabulated_force[tab_index];

        fx[i] += f * diffX;
        fy[i] += f * diffY;
        fx[j] -= f * diffX;
        fy[j] -= f * diffY;
    }
}

void calculateExternalForces()
{
    for (int i = 0; i < N; i++)
    {
        fx[i] += 5.0 * current_time / total_time;
    }
}

void calculatePinningSitesForce() 
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N_pinning_sites; j++)
		{
			double diffX = x[i] - x_pinning_site[j];
			double diffY = y[i] - y_pinning_site[j];

			if (diffX < -Sx_2) diffX += Sx;
	        if (diffX > Sx_2) diffX -=Sx;
	        if (diffY < -Sy_2) diffY += Sy;
	        if (diffY > Sy_2) diffY -= Sy;

        	double distance = diffX * diffX + diffY * diffY;
        	if (distance < r_pinning_site * r_pinning_site)
        	{
        		double f = f_max_pinning_sites / r_pinning_site;
        		fx[i] += -f * diffX;
        		fy[i] += -f * diffY;
        	}
		}
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
        if (rebuild_verlet_flag == 0)
            if (x_so_far[i] * x_so_far[i] + y_so_far[i] * y_so_far[i] >= r_travel_max * r_travel_max)
                rebuild_verlet_flag = 1;

        fx[i] = 0.0;
        fy[i] = 0.0;
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

void writeContourFile()
{
    ofstream f("../Plotter/contour.txt");

    f << N_pinning_sites << endl;

	for(int i = 0; i < N_pinning_sites; i++)
	{
		f << x_pinning_site[i] << endl;
		f << y_pinning_site[i] << endl;
		f << r_pinning_site << endl;
		f << r_pinning_site << endl;
		f << r_pinning_site << endl;
	}
    f.close();
}

void writeCmovie(FILE* moviefile, int t)
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
    cout << "Generic BD simulation" << endl;
    cout << "Simulating spontaneous lane formation" << endl;

    initParticles();
    FILE* f;
    if (argc == 2) 
        readDataFromFile(argv[1]);

    f = fopen(moviefile, "wb");
    cout << "Sx = " << Sx << endl;
    cout << "Sy = " << Sy << endl;
    cout << "N = " << N << endl;
    cout << "moviefile: " << moviefile << endl;

    generateCoordinates();
    initPinningSites();
    writeContourFile();
    calculateTabulatedForces();

    start_timing();
    for (current_time = 0; current_time < total_time; current_time++) 
    {
        if (current_time % time_echo == 0)
        {
            double perc = (double) current_time / (double) total_time * 100;
            cout << current_time << "/" << total_time << " " << perc << "%" << endl;
        }

        if (rebuild_verlet_flag == 1) 
        {
            calculateVerletList();
            colorVerlet();
        }

        calculateForces();
        calculateExternalForces();
        calculatePinningSitesForce();
        moveParticles();

        if (current_time % 100 == 0)
            writeCmovie(f, current_time);
    }
    freeData();
    stop_timing();
    return 0;
}