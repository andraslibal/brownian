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

#define PI 3.14159265358979323846264338327950288419716939937510

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
// radius of a pinning site
double r_pinning_site;
double half_length_pinning_site;
// pinning site's length
double pin_length;
// forces in pinning site
// in the ends of the pinning site
double K_max_pinning_site;
// in middle
double K_middle_pinning_site;
// coordinates of pinning sites (center)
double *x_pinning_site;
double *y_pinning_site;
// angle of the pinning site
double *cos_fi;
double *sin_fi;
// particles ID in pinning site
int* particle_ID;
// verteces near the pinning site
int** verteces_id_near_pinning_site;

// number of verteces
int N_vertex;
// coordinates of verteces
double* x_vertex;
double* y_vertex;
// vertex type
int* vertex_type;
// vertex color
int* vertex_color;
// number of neighbors for the vertex
int vertex_neighbor_number;
// vertex's z type
int *vertex_z_type;
// neighbor pinning site's id
int** vertex_neighbor_id;
int* vertex_chess_color;
int* particle_is_near_to_vertex;

// temperature
double T;

int	pinning_lattice_Nx, pinning_lattice_Ny;
double pinning_lattice_ax, pinning_lattice_ay;

char* moviefile;

//time
double dt = 0.0;
time_t beginning_time = 0;
time_t ending_time = 0;

//total runtime
int time_echo;
int current_time = 0;
int total_time = 0;
double multiplier;

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
		f >> Sy;
		f >> N;
	}

	f.close();
}

void initSize()
{
	strcpy(moviefile, "moviefile.mvi");
	Sx = 320.0;
	Sy = 320.0;
	N = 51200;
}

void initParticles()
{
	N = N_pinning_sites;
	ID = new int[N];

	x = new double[N];
	y = new double[N];

	fx = new double[N];
	fy = new double[N];

	color =  new int[N];
	q = new double[N];

    particle_is_near_to_vertex = new int[N];

	N_tabulate = 50000;
	tabulated_force = new double[N_tabulate];

	r0 = 4.0;
	r_verlet = 6.0;
	r_travel_max = r_verlet - r0;

	dt = 0.01;
	current_time = 0;
	total_time = 300000;
    multiplier = 0;
    time_echo = 500;

	rebuild_verlet_flag = 1;
	x_so_far = new double[N];
	y_so_far = new double[N];
	printf("Initialization complete\n");fflush(stdout);
}

void initArraysForPinningSites(int multiplier)
{
	int i;

	r_pinning_site = 0.2;
	half_length_pinning_site = 0.6;
	K_max_pinning_site = 5.0;
	K_middle_pinning_site = 0.1;
	T = 3.0;

	pinning_lattice_Nx = 10;
	pinning_lattice_Ny = 10;
	pinning_lattice_ax = 2.0;
	pinning_lattice_ay = 2.0;

	// calculating the system's size depending on the number of pinning sites
	pin_length = 2 * (half_length_pinning_site + r_pinning_site);

	if (multiplier == 2) 
	{
		Sx = pinning_lattice_Nx * pinning_lattice_ax;
		Sy = pinning_lattice_Ny * pinning_lattice_ay;
		N_vertex = pinning_lattice_Nx * pinning_lattice_Ny;
	    vertex_neighbor_number = 4;
	} else 
	if (multiplier == 6)
	{
		Sx = pinning_lattice_Nx * pinning_lattice_ax * 3;
		Sy = pinning_lattice_Ny * pinning_lattice_ay * sqrt(3);
		N_vertex = pinning_lattice_Nx * pinning_lattice_Ny * 4;
	    vertex_neighbor_number = 3;
	}
	Sx_2 = Sx / 2;
	Sy_2 = Sy / 2;

	// calculating number of pinning sites in system
	N_pinning_sites = pinning_lattice_Nx * pinning_lattice_Ny * multiplier;

	particle_ID = new int[N_pinning_sites];

	x_pinning_site = new double[N_pinning_sites];
	y_pinning_site = new double[N_pinning_sites];

	sin_fi = new double[N_pinning_sites];
	cos_fi = new double[N_pinning_sites];

	verteces_id_near_pinning_site = new int*[N_pinning_sites];
	for (i = 0; i < N_pinning_sites; i++)
		verteces_id_near_pinning_site[i] = new int[2];

	// init verteces array
	x_vertex = new double[N_vertex];
	y_vertex = new double[N_vertex];
	vertex_color = new int[N_vertex];
	vertex_type = new int[N_vertex];
    vertex_z_type = new int[N_vertex];
	vertex_neighbor_id = new int*[N_vertex];
    vertex_chess_color = new int[N_vertex];
}

void freeData()
{
	int i;

	delete[] ID;
	delete[] x;
	delete[] y;
	delete[] fx;
	delete[] fy;
	delete[] color;
	delete[] q;
	delete[] x_so_far;
	delete[] y_so_far;
    delete[] particle_is_near_to_vertex;

	delete[] x_pinning_site;
	delete[] y_pinning_site;
	delete[] cos_fi;
	delete[] sin_fi;
	delete[] particle_ID;

	for (i = 0; i < N_pinning_sites; i++)
		delete[] verteces_id_near_pinning_site[i];
	delete[] verteces_id_near_pinning_site;

	delete[] tabulated_force;

	delete[] x_vertex;
	delete[] y_vertex;
	delete[] vertex_color;
	delete[] vertex_type;
	delete[] vertex_z_type;

	for (i = 0; i < N_vertex; i++)
		delete[] vertex_neighbor_id[i];
	delete[] vertex_neighbor_id;
    delete[] vertex_chess_color;
}

void generateCoordinates()
{
	int i, j;
	double tmpX, tmpY;
	double diffX, diffY, diffSquare;
	for (i = 0; i < N; i++)
	{
		ID[i] = i;
		do {
			// choose a random position in the system for the particle
			tmpX = Sx * Rand();
			tmpY = Sy * Rand();
			// checking if in this position there already is an element
			for (j = 0; j < i; j++)
			{
				diffX = x[j] - tmpX;
				diffY = y[j] - tmpY;

				if (diffX > Sx_2) diffX -= Sx;
				if (diffX < -Sx_2) diffX += Sx;
				if (diffY > Sy_2) diffY -= Sy;
				if (diffY < -Sy_2) diffY += Sy;

				diffSquare = diffX * diffX + diffY * diffY;
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

	printf("Generated coordinates for the particles\n");fflush(stdout);
}

void putParticleInPinningSite() 
{
	int i;
	for (i = 0; i < N_pinning_sites; i++)
	{
		particle_ID[i] = ID[i] = i;
		
		if (Rand() > 0.5) {
			if (cos_fi[i] == 1.0) {
				x[i] = x_pinning_site[i] - half_length_pinning_site;
				y[i] = y_pinning_site[i];
			} else {
				if (cos_fi[i] == 0.0) {
					x[i] = x_pinning_site[i];
					y[i] = y_pinning_site[i] - half_length_pinning_site;
				} else {
					x[i] = x_pinning_site[i] - cos_fi[i] * half_length_pinning_site;
					y[i] = y_pinning_site[i] - sin_fi[i] * half_length_pinning_site;
				}
			}
		} else {
			if (cos_fi[i] == 1.0) {
				x[i] = x_pinning_site[i] + half_length_pinning_site;
				y[i] = y_pinning_site[i];
			} else {
				if (cos_fi[i] == 0.0) {
					x[i] = x_pinning_site[i];
					y[i] = y_pinning_site[i] + half_length_pinning_site;
				} else {
					x[i] = x_pinning_site[i] + cos_fi[i] * half_length_pinning_site;
					y[i] = y_pinning_site[i] + sin_fi[i] * half_length_pinning_site;
				}
			}
		}

		color[i] = 1;
		q[i] = 1.0;

		// setting the initial forces
		fx[i] = 0.0;
		fy[i] = 0.0;
		x_so_far[i] = 0.0;
		y_so_far[i] = 0.0;
        particle_is_near_to_vertex[i] = 0;
	}
}

void initSquareVerteces()
{
	int i, j, k;
	double x, y;

	x = y = 0.0;
	k = 0;

	for (i = 0; i < pinning_lattice_Ny; i++)
	{
		x = 0.0;
		for (j = 0; j < pinning_lattice_Nx; j++)
		{
            // 0 - pinning site from right
            // 1 - pinning site from bottom
            // 2 - pinning site from top
            // 3 - pinning site from left
			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_type[k] = 0;
			vertex_color[k] = 4;
            vertex_chess_color[k] = (i + j) % 2;
            vertex_z_type[k] = 4;

			vertex_neighbor_id[k] = new int[vertex_neighbor_number];
			if (i == 0) {
                if (j == 0) {
                    vertex_neighbor_id[k][1] = 2 * (pinning_lattice_Nx * pinning_lattice_Ny) - 1;
                } else {
                    vertex_neighbor_id[k][1] = 2 * (pinning_lattice_Nx * (pinning_lattice_Ny - 1) + j) - 1;
                }
			} else {
				if (j == 0) {
					vertex_neighbor_id[k][1] = i * pinning_lattice_Nx * 2 - 1;
				} else {
					vertex_neighbor_id[k][1] = 2 * ((i - 1) * pinning_lattice_Nx + j) - 1;
				}
			}
			vertex_neighbor_id[k][0] = 2 * k;
			if (j == 0) {
				vertex_neighbor_id[k][2] = pinning_lattice_Nx * 2 * (i + 1) - 1;
				vertex_neighbor_id[k][3] = pinning_lattice_Nx * 2 * (i + 1) - 2;
			} else  {
				vertex_neighbor_id[k][2] = 2 * (pinning_lattice_Nx * i + j - 1) + 1;
				vertex_neighbor_id[k][3] = 2 * (pinning_lattice_Nx * i + j - 1);
			}
			x += pinning_lattice_ax;
			k++;
		}
		y += pinning_lattice_ay;
	}

	for (i = 0; i < N_vertex; i++) {
		cout << i << ": " << vertex_neighbor_id[i][0] << " " << vertex_neighbor_id[i][1] << " " << vertex_neighbor_id[i][2] <<" " << vertex_neighbor_id[i][3] << endl;
	}
	
	printf("Init square verteces complet\n");
	fflush(stdout);
}

void initHexaVerteces()
{
	int i, j, k;
	double x, y, sqrt3;

	sqrt3 = sqrt(3);
	x = 0.0;
	y = pinning_lattice_ax * sqrt3 / 2;
	k = 0;

	for (i = 0; i < pinning_lattice_Nx; i++)
	{
		x = 0.0;
		for (j = 0; j < pinning_lattice_Ny; j++)
		{
            // be kell tenni a chess color meghatarozasar es annak a meghatarozasat, hogy melyik pinning siteok vannak hozza kozel
			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_type[k] = 0;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax;
			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_type[k] = 0;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax / 2;
			y -= pinning_lattice_ay * sqrt3 / 2;

			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_type[k] = 0;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax;
			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_type[k] = 0;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax / 2;
			y += pinning_lattice_ay * sqrt3 / 2;
		}
		y += pinning_lattice_ay * sqrt3;
	}
	
	printf("Init hexa verteces complet\n");
	fflush(stdout);
}

void initPinningSites()
{
	int i, j;
	double tmpX, tmpY;
	double diffX, diffY, diffSquare;
	initArraysForPinningSites(1);
	for (i = 0; i < N_pinning_sites; i++)
	{   
		do {
			// choose a random position in the system for the pinning site
			tmpX = Sx * Rand();
			tmpY = Sy * Rand();
			// checking if in this position there already is an element
			for (j = 0; j < i; j++)
			{
				diffX = x_pinning_site[j] - tmpX;
				diffY = y_pinning_site[j] - tmpY;

				if (diffX > Sx_2) diffX -= Sx;
				if (diffX < -Sx_2) diffX += Sx;
				if (diffX > Sy_2) diffX -= Sy;
				if (diffX < -Sy_2) diffX += Sy;

				diffSquare = diffX * diffX + diffY * diffY;
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

void initSquarePinningSites()
{
	double x, y, ax, ay;
	int i, j, k;

	x = y = 0.0;
	k = 0;

	initArraysForPinningSites(2);

	for (i = 0; i < pinning_lattice_Ny; i++)
	{
		x = pinning_lattice_ax / 2;
		for (j = 0; j < pinning_lattice_Nx; j++)
		{
			x_pinning_site[k] = x;
			y_pinning_site[k] = y;
			
			verteces_id_near_pinning_site[k][0] = i * pinning_lattice_Nx + j;
			verteces_id_near_pinning_site[k][1] = j == pinning_lattice_Nx - 1 ? i * pinning_lattice_Nx : i * pinning_lattice_Nx + j + 1;

			sin_fi[k] = 0.0;
			cos_fi[k++] = 1.0;

			x += pinning_lattice_ax / 2;
			y += pinning_lattice_ay / 2;
			x_pinning_site[k] = x;
			y_pinning_site[k] = y;
			
			verteces_id_near_pinning_site[k][0] = j == pinning_lattice_Nx - 1 ? i * pinning_lattice_Nx : i * pinning_lattice_Nx + j + 1;
			if (i == pinning_lattice_Ny - 1)
				if (j == pinning_lattice_Nx - 1)
					verteces_id_near_pinning_site[k][1] = 0;
				else
					verteces_id_near_pinning_site[k][1] = j + 1;
			else
				if (j == pinning_lattice_Nx - 1)
					verteces_id_near_pinning_site[k][1] = (i + 1) * pinning_lattice_Nx;
				else
					verteces_id_near_pinning_site[k][1] = (i + 1) * pinning_lattice_Nx + j + 1;

			sin_fi[k] = 1.0;
			cos_fi[k++] = 0.0;

			x += pinning_lattice_ax / 2;
			y -= pinning_lattice_ay / 2;
		}
		y += pinning_lattice_ay;
	}

  printf("Initialized square pinning array\n");fflush(stdout);  
}

void initHexaPinningSites()
{
	double x, y, sqrt3;
	int i, j, k;

	initArraysForPinningSites(6);

	sqrt3 = sqrt(3);
	x = 0.0;
	y = pinning_lattice_ax * sqrt3 / 2;
	k = 0;

	for (i = 0; i < pinning_lattice_Nx; i++)
	{
		x = pinning_lattice_ax / 2;
		for (j = 0; j < pinning_lattice_Ny; j++)
		{
			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = 0.0;
			cos_fi[k++] = 1.0;

			x += 3 * pinning_lattice_ax / 4;
			y += pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = sin(PI / 3);
			cos_fi[k++] = cos(PI / 3);

			y -= pinning_lattice_ay * sqrt3 / 2;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = sin(-PI / 3);
			cos_fi[k++] = cos(-PI / 3);

			x += 3 * pinning_lattice_ax / 4;
			y -= pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = 0.0;
			cos_fi[k++] = 1.0;

			x += 3 * pinning_lattice_ax / 4;
			y += pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = sin(PI / 3);
			cos_fi[k++] = cos(PI / 3);

			y += pinning_lattice_ay * sqrt3 / 2;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			sin_fi[k] = sin(-PI / 3);
			cos_fi[k++] = cos(-PI / 3);

			x += 3 * pinning_lattice_ax / 4;
			y -= pinning_lattice_ax  * sqrt3 / 4;
		}
		y += pinning_lattice_ay * sqrt3;
	}
}

void calculateGroudState(int i)
{
    int right_pinning_site,
        left_pinning_site,
        bottom_pinning_site,
        top_pinning_site;
    
    right_pinning_site = vertex_neighbor_id[i][0];
    left_pinning_site = vertex_neighbor_id[i][3];
    bottom_pinning_site = vertex_neighbor_id[i][1];
    top_pinning_site = vertex_neighbor_id[i][2];

    if (right_pinning_site != -1 && left_pinning_site != -1)
        if (particle_is_near_to_vertex[right_pinning_site] && particle_is_near_to_vertex[left_pinning_site])
            vertex_type[i] = 5 + vertex_chess_color[i];
        else
            if (bottom_pinning_site != -1 && top_pinning_site != -1)
                if (particle_is_near_to_vertex[bottom_pinning_site] && particle_is_near_to_vertex[top_pinning_site])
                    vertex_type[i] = 6 - vertex_chess_color[i];
}

void calculateVertexType()
{
	int i, j;
	int particle_id;
	double diffX, diffY;
    double x_rotated, y_rotated;

	// for each vertex is calculating it's type, which means iterating through it's neighbor 
	// pinning sites and defining whather it's particle is close enough to the vertex or not
	for (i = 0; i < N_vertex; i++)
	{
		vertex_type[i] = 0;
		for (j = 0; j < vertex_neighbor_number; j++)
		{
			particle_id = particle_ID[vertex_neighbor_id[i][j]];
			// calculating the particles distance, from the vertex
			diffX = x[particle_id] - x_vertex[i];
			diffY = y[particle_id] - y_vertex[i];

			if (diffX < -Sx_2) diffX += Sx;
			if (diffX > Sx_2) diffX -= Sx;
			if (diffY < -Sy_2) diffY += Sy;
			if (diffY > Sy_2) diffY -= Sy;

            x_rotated = diffX * cos_fi[particle_id] + diffY * sin_fi[particle_id];
		    y_rotated = -diffX * sin_fi[particle_id] + diffY * cos_fi[particle_id];

			// it won't work well if pinning_lattice_az is different from pinning_lattice_ay
            // if the particle is closet than the half of the pinning length, then:
			if (x_rotated * x_rotated + y_rotated * y_rotated <= pinning_lattice_ax * pinning_lattice_ax / 4.0)
            {
				vertex_type[i] ++;
                particle_is_near_to_vertex[particle_id] = 1;
            } else {
                particle_is_near_to_vertex[particle_id] = 0;
            }
		}
        if (vertex_type[i] == 2) calculateGroudState(i);
		vertex_color[i] = 2 + vertex_type[i];
	}
}

void calculateTabulatedForces() 
{
	int i;
	double r, r2;

	r = 0.1;
	r2 = r * r;

	tab_measure = (r_verlet * r_verlet - r2) / (N_tabulate - 1.0);
	tab_start = r;

	for (i = 0; i < N_tabulate; i++) {
		tabulated_force[i] = exp(- r / r0) / r2 / r;
		r2 += tab_measure;
		r = sqrt(r2);
	}
}

void calculateVerletList()
{
	int i, j;
	double diffX, diffY, distance;

	// clearing previous data from verlet list
	for (i = 0; i < verlet.size(); i++)
		verlet.at(i).clear();
	verlet.clear();

	// iterating through the particles
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			diffX = x[i] - x[j];
			diffY = y[i] - y[j];

			if (diffX < -Sx_2) diffX += Sx;
			if (diffX > Sx_2) diffX -= Sx;
			if (diffY < -Sy_2) diffY += Sy;
			if (diffY > Sy_2) diffY -= Sy;

			distance = diffX * diffX + diffY * diffY;

			// if they are closer than a specified value, that means that they are in interaction so I put them in the verlet list
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
	for (i = 0; i < N; i++)
	{
		x_so_far[i] = 0.0;
		y_so_far[i] = 0.0;
	}

	//clear rebuild flag
	rebuild_verlet_flag = 0;
	
	printf("Verlet list built\n");
	fflush(stdout);
}

void colorVerlet()
{
	int i;
	for(i = 0; i < N; i++)
	{
		if (q[i] == -1.0) color[i] = 0;        
		if (q[i] == 1.0) color[i] = 1;
	}

	for(i = 0; i < verlet.size(); i++)
	{
		if (verlet.at(i).at(0) == 30) color[verlet.at(i).at(1)] = 2;
		if (verlet.at(i).at(1) == 30) color[verlet.at(i).at(0)] = 2;
	}
}

void calculatePairwiseForces()
{
	int it, i, j, tab_index;
	double f;
	double diffX, diffY, distance;

	for (it = 0; it < verlet.size(); it++)
	{

		i = verlet.at(it).at(0);
		j = verlet.at(it).at(1);
		
		diffX = x[i] - x[j];
		diffY = y[i] - y[j];

		if (diffX < -Sx_2) diffX += Sx;
		if (diffX > Sx_2) diffX -= Sx;
		if (diffY < -Sy_2) diffY += Sy;
		if (diffY > Sy_2) diffY -= Sy;

		distance = diffX * diffX + diffY * diffY;
		if (distance < 0.1) tab_index = 0;
		else tab_index = (int) floor((distance - tab_start) / tab_measure);
  
		if (tab_index >= 50000) f = 0.0;
		else f = tabulated_force[tab_index] * multiplier;

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

void calculateThermalForce()
{
	int i;
	for (i = 0; i < N; i++)
	{
		fx[i] += gasdev() * T;
		fy[i] += gasdev() * T;
	}
}

void calculatePinningSitesForce() 
{
	int i, j;
	double diffX, diffY, distance, f;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N_pinning_sites; j++)
		{
			diffX = x[i] - x_pinning_site[j];
			diffY = y[i] - y_pinning_site[j];

			if (diffX < -Sx_2) diffX += Sx;
			if (diffX > Sx_2) diffX -= Sx;
			if (diffY < -Sy_2) diffY += Sy;
			if (diffY > Sy_2) diffY -= Sy;

			distance = diffX * diffX + diffY * diffY;
			if (distance < r_pinning_site * r_pinning_site)
			{
				f = K_max_pinning_site / r_pinning_site;
				fx[i] -= f * diffX;
				fy[i] -= f * diffY;
			}
		}
	}
}

void calculateModifiedPinningiteForces()
{
	int i, j;
	double diffX, diffY, distance;
	double x_rotated, y_rotated;
	double fx_rotated, fy_rotated;
	double f;
	for (i = 0; i < N_pinning_sites; i++)
	{
		j = particle_ID[i];

		diffX = x[j] - x_pinning_site[i];
		diffY = y[j] - y_pinning_site[i];

		if (diffX < -Sx_2) diffX += Sx;
		if (diffX > Sx_2) diffX -= Sx;
		if (diffY < -Sy_2) diffY += Sy;
		if (diffY > Sy_2) diffY -= Sy;

		x_rotated = diffX * cos_fi[i] + diffY * sin_fi[i];
		y_rotated = -diffX * sin_fi[i] + diffY * cos_fi[i];

		if (x_rotated > half_length_pinning_site) 
		{
			x_rotated = x_rotated - half_length_pinning_site;
			fx_rotated = -K_max_pinning_site * x_rotated;
			fy_rotated = -K_max_pinning_site * y_rotated;

		} else
		if (x_rotated < -half_length_pinning_site) 
		{
			x_rotated = x_rotated + half_length_pinning_site;
			fx_rotated = -K_max_pinning_site * x_rotated;
			fy_rotated = -K_max_pinning_site * y_rotated;
		} else {
			fx_rotated = K_middle_pinning_site * x_rotated;
			fy_rotated = -K_max_pinning_site * y_rotated;
		}

		fx[j] += fx_rotated * cos_fi[i] - fy_rotated * sin_fi[i];
		fy[j] += fx_rotated * sin_fi[i] + fy_rotated * cos_fi[i];
		
		// circle shaped
		//fx[j] += -K_max_pinning_site * diffX;
		//fy[j] += -K_max_pinning_site * diffY;
	}
}

void moveParticles()
{
	int i;
	double deltax, deltay;
	double maxfx = 0.0, maxfy = 0.0;

	for (i = 0; i < N; i++) 
	{
		maxfx = maxfx < fx[i] ? fx[i] : maxfx;
		maxfy = maxfy < fy[i] ? fy[i] : maxfy;

		deltax = fx[i] * dt;
		deltay = fy[i] * dt;

		x[i] += deltax;
		y[i] += deltay;

		x_so_far[i] += deltax;
		y_so_far[i] += deltay;

		if (x[i] < 0) x[i] += Sx;
		if (x[i] > Sx) x[i] -= Sx;
		if (y[i] < 0) y[i] += Sy;
		if (y[i] > Sy) y[i] -= Sy;

		//only check on unset flag
		if (rebuild_verlet_flag == 0)
			if (x_so_far[i] * x_so_far[i] + y_so_far[i] * y_so_far[i] >= r_travel_max * r_travel_max)
				rebuild_verlet_flag = 1;

		fx[i] = 0.0;
		fy[i] = 0.0;
	}
	// cout << maxfx << " " << maxfy << endl;
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
	double sqrt3;
    int i, j;
	
	sqrt3 = sqrt(3);
	f << N_pinning_sites * 3 << endl;

    // for(j = 0; j < vertex_neighbor_number; j++)
	// {
    //     i = vertex_neighbor_id[49][j];
	// 	f << x_pinning_site[i] << endl;
	// 	f << y_pinning_site[i] << endl;
	// 	f << r_pinning_site << endl;
	// 	f << r_pinning_site << endl;
	// 	f << r_pinning_site << endl;

	// 	if (cos_fi[i] == 1)
	// 	{
	// 		f << x_pinning_site[i] + half_length_pinning_site << endl;
	// 		f << y_pinning_site[i] << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;

	// 		f << x_pinning_site[i] - half_length_pinning_site << endl;
	// 		f << y_pinning_site[i] << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 	} else 
	// 	if (cos_fi[i] == cos(PI / 3) && sin_fi[i] == sin(PI / 3)) {
	// 		f << x_pinning_site[i] + half_length_pinning_site / 2 << endl;
	// 		f << y_pinning_site[i] + half_length_pinning_site * sqrt3 / 2 << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;

	// 		f << x_pinning_site[i] - half_length_pinning_site / 2 << endl;
	// 		f << y_pinning_site[i] - half_length_pinning_site * sqrt3 / 2 << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 	} else 
	// 	if (cos_fi[i] == 0) {
	// 		f << x_pinning_site[i] << endl;
	// 		f << y_pinning_site[i] + half_length_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;

	// 		f << x_pinning_site[i] << endl;
	// 		f << y_pinning_site[i] - half_length_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 	} else {
	// 		f << x_pinning_site[i] - half_length_pinning_site / 2<< endl;
	// 		f << y_pinning_site[i] + half_length_pinning_site * sqrt3 / 2 << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;

	// 		f << x_pinning_site[i] + half_length_pinning_site / 2 << endl;
	// 		f << y_pinning_site[i] - half_length_pinning_site * sqrt3 / 2 << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 		f << r_pinning_site << endl;
	// 	}
	// }

	for(i = 0; i < N_pinning_sites; i++)
	{
		f << x_pinning_site[i] << endl;
		f << y_pinning_site[i] << endl;
		f << r_pinning_site << endl;
		f << r_pinning_site << endl;
		f << r_pinning_site << endl;

		if (cos_fi[i] == 1)
		{
			f << x_pinning_site[i] + half_length_pinning_site << endl;
			f << y_pinning_site[i] << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;

			f << x_pinning_site[i] - half_length_pinning_site << endl;
			f << y_pinning_site[i] << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
		} else 
		if (cos_fi[i] == cos(PI / 3) && sin_fi[i] == sin(PI / 3)) {
			f << x_pinning_site[i] + half_length_pinning_site / 2 << endl;
			f << y_pinning_site[i] + half_length_pinning_site * sqrt3 / 2 << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;

			f << x_pinning_site[i] - half_length_pinning_site / 2 << endl;
			f << y_pinning_site[i] - half_length_pinning_site * sqrt3 / 2 << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
		} else 
		if (cos_fi[i] == 0) {
			f << x_pinning_site[i] << endl;
			f << y_pinning_site[i] + half_length_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;

			f << x_pinning_site[i] << endl;
			f << y_pinning_site[i] - half_length_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
		} else {
			f << x_pinning_site[i] - half_length_pinning_site / 2<< endl;
			f << y_pinning_site[i] + half_length_pinning_site * sqrt3 / 2 << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;

			f << x_pinning_site[i] + half_length_pinning_site / 2 << endl;
			f << y_pinning_site[i] - half_length_pinning_site * sqrt3 / 2 << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
			f << r_pinning_site << endl;
		}
	}
	f.close();
	
	printf("Written the contour file\n");
	fflush(stdout);
}

void writeCmovie(FILE* moviefile, int t)
{
	int i;
	float floatholder;
	int intholder;

	intholder = N + N_vertex;
	fwrite(&intholder, sizeof(int), 1, moviefile);

	intholder = t;
	fwrite(&intholder, sizeof(int) ,1, moviefile);
	
	for (i = 0; i < N; i++)
	{
		intholder = color[i] + 2;
		fwrite(&intholder, sizeof(int), 1, moviefile);
		intholder = ID[i]; //ID
		fwrite(&intholder, sizeof(int), 1, moviefile);
		floatholder = (float)x[i];
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = (float)y[i];
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = 1.0; //cum_disp, cmovie format
		fwrite(&floatholder, sizeof(float), 1, moviefile);
	}

	for (i = 0; i < N_vertex; i++)
	{
		intholder = vertex_color[i] + 2;
		fwrite(&intholder, sizeof(int), 1, moviefile);
		intholder = N + i; //ID
		fwrite(&intholder, sizeof(int), 1, moviefile);
		floatholder = (float)x_vertex[i];
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = (float)y_vertex[i];
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = 1.0; //cum_disp, cmovie format
		fwrite(&floatholder, sizeof(float), 1, moviefile);
	}
}

void writeGfile()
{
	ofstream f("../Plotter/gfile");

	f << "set xrange -1 " << (int) Sx + 1 << endl;
	f << "set yrange -1 " << (int) Sy + 1 << endl;
	f << "loadcolormap colors.txt" << endl;
	f << "loadcontour contour.txt" << endl;
	f << "cmovie" << endl;

	f.close();
}

void writeStatistics()
{

}

int main(int argc, char* argv[]) 
{
	cout << "Generic BD simulation" << endl;
	cout << "Simulating spontaneous lane formation" << endl;

	FILE* f;
	moviefile = new char[30];
	if (argc == 2) 
		readDataFromFile(argv[1]);
	else
		initSize();

	f = fopen(moviefile, "wb");

	initSquarePinningSites();
	initSquareVerteces();
	initParticles();
	writeGfile();
	cout << "Sx = " << Sx << endl;
	cout << "Sy = " << Sy << endl;
	cout << "N = " << N << endl;
	cout << "moviefile: " << moviefile << endl;

	// generateCoordinates();
	putParticleInPinningSite();
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
			// colorVerlet();
		}

		calculatePairwiseForces();
		// calculateExternalForces();
		calculateThermalForce();
		calculateModifiedPinningiteForces();

		moveParticles();

		if (current_time % 100 == 0)
        {
            calculateVertexType();
			writeCmovie(f, current_time);
            multiplier = current_time / total_time;
            writeStatistics();
        }
	}
	fclose(f);
	freeData();
	stop_timing();
	return 0;
}