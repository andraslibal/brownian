#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <regex.h>

#include "random_numrec.c"
#include "jsmn.h"

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
// pinning site's id, in which the given particle is
int* pinning_ID;
// spin flip count
int* spin_flip;
// where is the particle in the pinning site (left or right end)
// -1 - left
//  1 - right
int* position_in_pinning_site;

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
// decimate number of pinning sites
int decimate_number;
// marking if pinning site is removed
int* is_pinning_site_removed;

// number of verteces
int N_vertex;
int N_vertex_z_type_3;
int N_vertex_z_type_4;
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
int min_z_type;
// topological charge on vertex
int* vertex_q;

// temperature
double T;

int	pinning_lattice_Nx, pinning_lattice_Ny;
double pinning_lattice_ax, pinning_lattice_ay;

char* moviefile;
char* statfile;

int* vertex_average;
int vertex_type_number;
int charge_3_near_3;
int charge_3_near_4;
int charge_4_near_3;
int charge_4_near_4;
int N_z3_charged;
int N_z4_charged;
int stat_steps;

//time
double dt = 0.0;
time_t beginning_time = 0;
time_t ending_time = 0;

//total runtime
int time_echo;
int current_time = 0;
int total_time = 0;
double multiplier, max_multiplier;

// variables for JSON parsing
string JSON;
jsmntok_t* t = NULL;

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
	string tmp;
	if (f.is_open())
	{
		while (!f.eof())
		{
			f >> tmp;
			JSON += tmp;
		}
		// for this command (string.pop_back()) 11th standard of c++ is needed
		if (JSON[JSON.length() - 2] == '}' && JSON[JSON.length() - 1] == '}')
			JSON.pop_back();
		cout << JSON << endl;
	} else {
		JSON = "{}";
	}

	f.close();
}

bool checkFileName(const char* regex_str, string toTest)
{
	regex_t regexCompiled;
	bool result;
	regmatch_t pmatch[1];

	if (regcomp(&regexCompiled, regex_str, REG_EXTENDED))
	{
		cout << "\033[1;31mCould not create regex!\033[0m" << endl;
		return false;
	}

	result = !regexec(&regexCompiled, toTest.c_str(), 1, pmatch, 0) ? ((pmatch[0].rm_so == 0) && (pmatch[0].rm_eo == toTest.length())) : false;
	regfree(&regexCompiled);
	return result;
}

int getToken(int start, int end, const char* token)
{
	int l2 = strlen(token);
	for (int i = start; i < end; i++)
	{
		if (t[i].type == JSMN_STRING)
		{
			int l1 = t[i].end - t[i].start;
			if (l1 == l2)
			{
				char* tmp = new char[l1];
				strncpy(tmp, JSON.substr(t[i].start, l1).c_str(), l1);
				if (strncmp(tmp, token, l1) == 0)
				{
					delete[] tmp;
					return i;
				}
				delete[] tmp;
			}
		}
	}
	return -1;
}

void initData()
{
	//T es decimated_number kimaradt
	jsmn_parser p;
	int n, r, x, y, object_end;
	size_t len, regex_len;
	char* regex_str;

	regex_len = 36;

	len = (size_t) JSON.length();
	jsmn_init(&p);
	n = jsmn_parse(&p, JSON.c_str(), len, t, 0);
	if (n < 0) {
		cout << "\033[1;31mSome error occured!\033[0m" << endl;
		return;
	} else {
		t = new jsmntok_t[n];
	}
	jsmn_init(&p);
	r = jsmn_parse(&p, JSON.c_str(), len, t, n);

	if (r < 1 && !t[0].type == JSMN_OBJECT) {
		cout << "\033[1;31mJSON object not found!\033[0m" << endl;
		return;
	}
	if ((x = getToken(1, r, "moviefile")) >= 0) {
		regex_str = new char[regex_len];
		strncpy(regex_str, "^[a-z0-9A-Z_]+$", (size_t) 15);
		string toTest = JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start);
		if (checkFileName(regex_str, toTest)) {
			toTest += ".mvi";
		}
		regex_str[strlen(regex_str) - 1] = '\0';
		strncat(regex_str, "(\\.mvi)$", (size_t) 7);
		if (!checkFileName(regex_str, toTest)) {
			strncpy(moviefile, "result.mvi", 11);
		} else {
			strncpy(moviefile, toTest.c_str(), toTest.length() + 1);
		}
		delete[] regex_str;
	} else
		strncpy(moviefile, "result.mvi", 10);
	cout << "moviefile: " << moviefile << endl;

	if ((x = getToken(1, r, "statfile")) >= 0) {
		regex_str = new char[regex_len];
		strncpy(regex_str, "^[a-z0-9A-Z_]+$", (size_t) 15);
		string toTest = JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start);
		if (checkFileName(regex_str, toTest)) {
			toTest += ".txt";
		}
		regex_str[strlen(regex_str) - 1] = '\0';
		strncat(regex_str, "(\\.txt)$", (size_t) 7);
		if (!checkFileName(regex_str, toTest)) {
			strncpy(statfile, "stat.txt", 9);
		} else {
			strncpy(statfile, toTest.c_str(), toTest.length() + 1);
		}
	} else
			strncpy(statfile, "stat.txt", 8);
	cout << "statfile: " << statfile << endl;
	
	if ((x = getToken(1, r, "pinning_lattice")) >= 0) {
		if (t[x+1].type == JSMN_OBJECT) {
			object_end = x + 1 + 2*t[x+1].size;
			if ((y = getToken(x + 2, object_end, "pinning_lattice_Nx")) >= 0)
				pinning_lattice_Nx = stoi(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				pinning_lattice_Nx = 10;

			if ((y = getToken(x + 2, object_end, "pinning_lattice_Ny")) >= 0)
				pinning_lattice_Ny = stoi(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				pinning_lattice_Ny = 10;

			if ((y = getToken(x + 2, object_end, "pinning_lattice_ax")) >= 0)
				pinning_lattice_ax = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				pinning_lattice_ax = 2.0;

			if ((y = getToken(x + 2, object_end, "pinning_lattice_ay")) >= 0)
				pinning_lattice_ay = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				pinning_lattice_ay = 2.0;

			if ((x = getToken(1, r, "decimate_number")) >= 0)
				decimate_number = stoi(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
			else
				decimate_number = pinning_lattice_Nx * pinning_lattice_Ny / 4;
		}
	} else {
		pinning_lattice_Nx = 10;
		pinning_lattice_Ny = 10;
		pinning_lattice_ax = 2.0;
		pinning_lattice_ay = 2.0;
	}
	cout << "pinning_lattice_Nx: " << pinning_lattice_Nx << endl;
	cout << "pinning_lattice_Ny: " << pinning_lattice_Ny << endl;
	cout << "pinning_lattice_ax: " << pinning_lattice_ax << endl;
	cout << "pinning_lattice_ay: " << pinning_lattice_ay << endl;
	cout << "decimate_number: " << decimate_number << endl;

	if ((x = getToken(1, r, "pinning_site")) >= 0) {
		if (t[x+1].type == JSMN_OBJECT) {
			object_end = x + 1 + 2*t[x+1].size;
			if ((y = getToken(x + 2, object_end, "r_pinning_site")) >= 0)
				r_pinning_site = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				r_pinning_site = 0.2;

			if ((y = getToken(x + 2, object_end, "half_length_pinning_site")) >= 0)
				half_length_pinning_site = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				half_length_pinning_site = 0.6;

			if ((y = getToken(x + 2, object_end, "K_max_pinning_site")) >= 0)
				K_max_pinning_site = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				K_max_pinning_site = 2.0;

			if ((y = getToken(x + 2, object_end, "K_middle_pinning_site")) >= 0)
				K_middle_pinning_site = stod(JSON.substr(t[y + 1].start, t[y + 1].end - t[y + 1].start));
			else
				K_middle_pinning_site = 0.15;
		}
	} else {
		r_pinning_site = 0.2;
		half_length_pinning_site = 0.6;
		K_max_pinning_site = 2.0;
		K_middle_pinning_site = 0.15;
	}
	cout << "r_pinning_site: " << r_pinning_site << endl;
	cout << "half_length_pinning_site: " << half_length_pinning_site << endl;
	cout << "K_max_pinning_site: " << K_max_pinning_site << endl;
	cout << "K_middle_pinning_site: " << K_middle_pinning_site << endl;

	if ((x = getToken(1, r, "T")) >= 0)
		T = stod(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		T = 3.0;

	if ((x = getToken(1, r, "N_tabulate")) >= 0)
		N_tabulate = stoi(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		N_tabulate = 50000;
	
	if ((x = getToken(1, r, "r0")) >= 0)
		r0 = stod(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		r0 = 4.0;
	
	if ((x = getToken(1, r, "r_verlet")) >= 0)
		r_verlet = stod(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		r_verlet = 6.0;
	
	if ((x = getToken(1, r, "dt")) >= 0)
		dt = stod(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		dt = 0.002;
	
	if ((x = getToken(1, r, "total_time")) >= 0)
		total_time = stoi(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		total_time = 100000;
	
	if ((x = getToken(1, r, "time_echo")) >= 0)
		time_echo = stoi(JSON.substr(t[x + 1].start, t[x + 1].end - t[x + 1].start));
	else
		time_echo = 500;
	cout << "T: " << T << endl;
	cout << "N_tabulated: " << N_tabulate << endl;
	cout << "r0: " << r0 << endl;
	cout << "r_verlet: " << r_verlet << endl;
	cout << "dt: " << dt << endl;
	cout << "total_time: " << total_time << endl;
	cout << "time_echo: " << time_echo << endl;
	delete[] t;
}

void initParticles()
{
	N = N_pinning_sites - decimate_number;
	ID = new int[N];

	x = new double[N];
	y = new double[N];

	fx = new double[N];
	fy = new double[N];

	color =  new int[N];
	q = new double[N];

    particle_is_near_to_vertex = new int[N];

	tabulated_force = new double[N_tabulate];

	r_travel_max = r_verlet - r0;

	current_time = 0;
    multiplier = 0.0;
    max_multiplier = 1.0;

	rebuild_verlet_flag = 1;
	x_so_far = new double[N];
	y_so_far = new double[N];
	pinning_ID = new int[N];
	spin_flip = new int[N];
	position_in_pinning_site = new int[N];
	printf("Initialization complete\n");fflush(stdout);
}

void initArraysForPinningSites(int multiplier)
{
	int i;
	// calculating the system's size depending on the number of pinning sites
	pin_length = 2 * (half_length_pinning_site + r_pinning_site);

	if (multiplier == 2) 
	{
		Sx = pinning_lattice_Nx * pinning_lattice_ax;
		Sy = pinning_lattice_Ny * pinning_lattice_ay;
		N_vertex = pinning_lattice_Nx * pinning_lattice_Ny;
	    vertex_neighbor_number = 4;
		min_z_type = 3;
	} else 
	if (multiplier == 6)
	{
		Sx = pinning_lattice_Nx * pinning_lattice_ax * 3;
		Sy = pinning_lattice_Ny * pinning_lattice_ay * sqrt(3);
		N_vertex = pinning_lattice_Nx * pinning_lattice_Ny * 4;
	    vertex_neighbor_number = 3;
		min_z_type = 2;
	}
	Sx_2 = Sx / 2;
	Sy_2 = Sy / 2;

	// calculating number of pinning sites in system
	N_pinning_sites = pinning_lattice_Nx * pinning_lattice_Ny * multiplier;
	N_vertex_z_type_3 = decimate_number * 2; 
	N_vertex_z_type_4 = N_vertex - N_vertex_z_type_3;

	particle_ID = new int[N_pinning_sites];

	x_pinning_site = new double[N_pinning_sites];
	y_pinning_site = new double[N_pinning_sites];

	sin_fi = new double[N_pinning_sites];
	cos_fi = new double[N_pinning_sites];

	verteces_id_near_pinning_site = new int*[N_pinning_sites];
	for (i = 0; i < N_pinning_sites; i++)
		verteces_id_near_pinning_site[i] = new int[2];
	is_pinning_site_removed = new int[N_pinning_sites];

	// init verteces array
	x_vertex = new double[N_vertex];
	y_vertex = new double[N_vertex];
	vertex_color = new int[N_vertex];
	vertex_type = new int[N_vertex];
    vertex_z_type = new int[N_vertex];
	vertex_neighbor_id = new int*[N_vertex];
    vertex_chess_color = new int[N_vertex];
	vertex_q = new int[N_vertex];

	vertex_type_number = 7 + 8;
	vertex_average = new int[vertex_type_number];
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
	delete[] pinning_ID;
	delete[] spin_flip;
	delete[] position_in_pinning_site;
    delete[] particle_is_near_to_vertex;

	delete[] x_pinning_site;
	delete[] y_pinning_site;
	delete[] cos_fi;
	delete[] sin_fi;
	delete[] particle_ID;

	for (i = 0; i < N_pinning_sites; i++)
		delete[] verteces_id_near_pinning_site[i];
	delete[] verteces_id_near_pinning_site;
	delete[] is_pinning_site_removed;

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
	delete[] vertex_q;

	delete[] vertex_average;
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

void zeros(int* array, int n)
{
	int i;
	#pragma omp parallel for simd
	for (i = 0; i < n; i++) array[i] = 0;
}

void zeros(double* array, int n)
{
	int i;
	#pragma omp parallel for simd
	for (i = 0; i < n; i++) array[i] = 0.0;
}

void putParticleInPinningSite() 
{
	int i, k;
	k = 0;
	for (i = 0; i < N_pinning_sites; i++)
	{
		if (!is_pinning_site_removed[i])
		{
			particle_ID[i] = ID[k] = k;
			pinning_ID[k] = i;
			
			if (Rand() > 0.5) {
				if (cos_fi[i] == 1.0) {
					x[k] = x_pinning_site[i] - half_length_pinning_site;
					y[k] = y_pinning_site[i];
				} else {
					if (cos_fi[i] == 0.0) {
						x[k] = x_pinning_site[i];
						y[k] = y_pinning_site[i] - half_length_pinning_site;
					} else {
						x[k] = x_pinning_site[i] - cos_fi[i] * half_length_pinning_site;
						y[k] = y_pinning_site[i] - sin_fi[i] * half_length_pinning_site;
					}
				}
				position_in_pinning_site[k] = -1;
			} else {
				if (cos_fi[i] == 1.0) {
					x[k] = x_pinning_site[i] + half_length_pinning_site;
					y[k] = y_pinning_site[i];
				} else {
					if (cos_fi[i] == 0.0) {
						x[k] = x_pinning_site[i];
						y[k] = y_pinning_site[i] + half_length_pinning_site;
					} else {
						x[k] = x_pinning_site[i] + cos_fi[i] * half_length_pinning_site;
						y[k] = y_pinning_site[i] + sin_fi[i] * half_length_pinning_site;
					}
				}
				position_in_pinning_site[k] = 1;
			}

			color[k] = 1;
			q[k] = 1.0;
			k++;
		}
	}
	
	// setting the initial forces
	zeros(fx, N);
	zeros(fy, N);
	zeros(x_so_far, N);
	zeros(y_so_far, N);
	zeros(spin_flip, N);
	zeros(particle_is_near_to_vertex, N);
	printf("Particles placed in system");
	fflush(stdout);
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

	// for (i = 0; i < N_vertex; i++) {
	// 	cout << i << ": " << vertex_neighbor_id[i][0] << " " << vertex_neighbor_id[i][1] << " " << vertex_neighbor_id[i][2] <<" " << vertex_neighbor_id[i][3] << endl;
	// }
	
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
			vertex_color[k++] = 4;

			x += pinning_lattice_ax;
			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax / 2;
			y -= pinning_lattice_ay * sqrt3 / 2;

			x_vertex[k] = x;
			y_vertex[k] = y;
			vertex_color[k++] = 4;

			x += pinning_lattice_ax;
			x_vertex[k] = x;
			y_vertex[k] = y;
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
			is_pinning_site_removed[k] = 0;

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
			is_pinning_site_removed[k] = 0;

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

	// meg kell hatarozni a szomszedis vertex id-kat

	for (i = 0; i < pinning_lattice_Nx; i++)
	{
		x = pinning_lattice_ax / 2;
		for (j = 0; j < pinning_lattice_Ny; j++)
		{
			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = 0.0;
			cos_fi[k++] = 1.0;

			x += 3 * pinning_lattice_ax / 4;
			y += pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = sin(PI / 3);
			cos_fi[k++] = cos(PI / 3);

			y -= pinning_lattice_ay * sqrt3 / 2;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = sin(-PI / 3);
			cos_fi[k++] = cos(-PI / 3);

			x += 3 * pinning_lattice_ax / 4;
			y -= pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = 0.0;
			cos_fi[k++] = 1.0;

			x += 3 * pinning_lattice_ax / 4;
			y += pinning_lattice_ay * sqrt3 / 4;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = sin(PI / 3);
			cos_fi[k++] = cos(PI / 3);

			y += pinning_lattice_ay * sqrt3 / 2;

			x_pinning_site[k] = x;
			y_pinning_site[k] = y;

			is_pinning_site_removed[k] = 0;

			sin_fi[k] = sin(-PI / 3);
			cos_fi[k++] = cos(-PI / 3);

			x += 3 * pinning_lattice_ax / 4;
			y -= pinning_lattice_ax  * sqrt3 / 4;
		}
		y += pinning_lattice_ay * sqrt3;
	}
}

void removePinningSiteFromVertexNeighborhood(int vertexID, int removeID)
{
	int i;
	for (i = 0; i < vertex_neighbor_number; i++) {
		if (vertex_neighbor_id[vertexID][i] == removeID)
		{
			vertex_neighbor_id[vertexID][i] = -1;
			break;
		}
	}
}

void removePinningsite(int index, int vertex1, int vertex2)
{
	is_pinning_site_removed[index] = 1;
	vertex_z_type[vertex1] --;
	vertex_z_type[vertex2] --;
	
	removePinningSiteFromVertexNeighborhood(vertex1, index);
	removePinningSiteFromVertexNeighborhood(vertex2, index);
}

void decimating()
{
	int i, index, id1, id2;
	for (i = 0; i < decimate_number; i++)
	{
		do {
			index = (int) (Rand() * N_pinning_sites) % N_pinning_sites;
			id1 = verteces_id_near_pinning_site[index][0];
			id2 = verteces_id_near_pinning_site[index][1];
		} while (vertex_z_type[id1] <= min_z_type || vertex_z_type[id2] <= min_z_type);
		removePinningsite(index, id1, id2);
	}
	printf("Decimation complete!\n"); fflush(stdout);
}

void calculateTopologicalChargeOnVertex()
{
	int i;
	#pragma omp parallel for simd
	for (i = 0; i < N_vertex; i++)
	{
		if (vertex_z_type[i] == 3 && vertex_type[i] == 4) 
			cout << "\033[1;31mError! A vertex with Z type 3 has type 4!\033[0m" << endl;
		vertex_q[i] = (-1) * vertex_z_type[i] + 2 * vertex_type[i];
	}
}

void calculateGroundState(int i)
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
			if (vertex_neighbor_id[i][j] != -1) {
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
				// if the particle is closer than the half of the pinning length, then:
				if (x_rotated * x_rotated + y_rotated * y_rotated <= pinning_lattice_ax * pinning_lattice_ax / 4.0)
				{
					vertex_type[i] ++;
					particle_is_near_to_vertex[particle_id] = 1;
				} else {
					particle_is_near_to_vertex[particle_id] = 0;
				}
			} else {
				particle_is_near_to_vertex[particle_id] = 0;
			}
		}
        if (vertex_type[i] == 2) calculateGroundState(i);
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
	tab_start = r2;

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

	color[30] = 6;

	for(i = 0; i < verlet.size(); i++)
	{
		if (verlet.at(i).at(0) == 30) color[verlet.at(i).at(1)] = 10;
		if (verlet.at(i).at(1) == 30) color[verlet.at(i).at(0)] = 10;
	}
}

void calculatePairwiseForces()
{
	int it, i, j, tab_index;
	double f;
	double diffX, diffY, distance2;

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

		distance2 = diffX * diffX + diffY * diffY;
		if (distance2 < tab_start) tab_index = 0;
		else tab_index = (int) floor((distance2 - tab_start) / tab_measure);
  
		if (tab_index >= N_tabulate) f = 0.0;
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
	int j;

	#pragma omp parallel for simd
	for (j = 0; j < N; j++)
	{
		int i;
		double diffX, diffY;
		double x_rotated, y_rotated;
		double fx_rotated, fy_rotated;
		i = pinning_ID[j];

		diffX = x[j] - x_pinning_site[i];
		diffY = y[j] - y_pinning_site[i];

		if (diffX < -Sx_2) diffX += Sx;
		if (diffX > Sx_2) diffX -= Sx;
		if (diffY < -Sy_2) diffY += Sy;
		if (diffY > Sy_2) diffY -= Sy;

		x_rotated = diffX * cos_fi[i] + diffY * sin_fi[i];
		y_rotated = -diffX * sin_fi[i] + diffY * cos_fi[i];

		if (x_rotated <= -half_length_pinning_site) 
		{
			if (position_in_pinning_site[j] != -1) {
				position_in_pinning_site[j] = -1;
				spin_flip[j] ++;
			}
			x_rotated = x_rotated + half_length_pinning_site;
			fx_rotated = -K_max_pinning_site * x_rotated;
			fy_rotated = -K_max_pinning_site * y_rotated;
		} else 
			if (x_rotated >= half_length_pinning_site) 
			{
				if (position_in_pinning_site[j] != 1) {
					position_in_pinning_site[j] = 1;
					spin_flip[j] ++;
				}
				x_rotated = x_rotated - half_length_pinning_site;
				fx_rotated = -K_max_pinning_site * x_rotated;
				fy_rotated = -K_max_pinning_site * y_rotated;
			} else {
				fx_rotated = K_middle_pinning_site * (half_length_pinning_site - fabs(x_rotated));
				if (x_rotated < 0) fx_rotated = -fx_rotated;
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

	for (i = 0; i < N; i++) 
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
		if (y[i] > Sy) y[i] -= Sy;

		//only check on unset flag
		if (rebuild_verlet_flag == 0)
			if (x_so_far[i] * x_so_far[i] + y_so_far[i] * y_so_far[i] >= r_travel_max * r_travel_max)
				rebuild_verlet_flag = 1;
	}

	zeros(fx, N);
	zeros(fy, N);
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
	f << (N_pinning_sites - decimate_number) * 3 << endl;

    // for(j = 0; j < vertex_neighbor_number; j++)
	// {
    //     i = vertex_neighbor_id[8][j];
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
		if (!is_pinning_site_removed[i])
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

void zerosStatistics()
{
	zeros(vertex_average, vertex_type_number);
	charge_3_near_3 = 0;
	charge_3_near_4 = 0;
	charge_4_near_3 = 0;
	charge_4_near_4 = 0;
	N_z3_charged = 0;
	N_z4_charged = 0;
	stat_steps = 0;
	zeros(spin_flip, N);
}

void calculateStatisticsPerStep()
{
	int i, j, add;
	int actual_z_type, neighbor_z_type;
	int q, id;
	#pragma omp parallel for simd
	for (i = 0; i < N_vertex; i++)
	{
		if (vertex_z_type[i] == 3)
			add = 7;
		else
			add = 0;

		vertex_average[add + vertex_type[i]] ++;
	}

	for (i = 0; i < N_vertex; i++)
	{
		actual_z_type = vertex_z_type[i];
		if (vertex_q[i])
			if (vertex_z_type[i] == 3)
				N_z3_charged ++;
			else
				N_z4_charged ++;

		for (j = 0; j < vertex_neighbor_number; j++)
		{
			if (vertex_neighbor_id[i][j] != -1)
			{
				id = vertex_neighbor_id[i][j];
				neighbor_z_type = vertex_z_type[id];
				q = vertex_q[id];

				if (actual_z_type == 3)
					if (neighbor_z_type == 3)
						charge_3_near_3 += q;
					else
						charge_3_near_4 += q;
				else
					if (neighbor_z_type == 3)
						charge_4_near_3 += q;
					else
						charge_4_near_4 += q;
			}
		}
	}
	stat_steps ++;
}

void writeStatistics()
{
	int i;
	int charge_type_z4[5];
	int charge_type_z3[5];
	int total_charge_on_z3 = 0;
	int total_charge_on_z4 = 0;
	int total_flips = 0;
	ofstream stat(statfile, ios_base::app);
	stat << current_time << " ";
	for (i = 0; i < 7; i++)
		stat << (double) vertex_average[i] / (double) stat_steps / (double) N_vertex_z_type_4 << " ";
	for (i = 7; i < vertex_type_number; i++)
		stat << (double) vertex_average[i] / (double) stat_steps / (double) N_vertex_z_type_3 << " ";

	charge_type_z4[0] = vertex_average[2] + vertex_average[5] + vertex_average[6]; // q = 0;
	charge_type_z4[1] = vertex_average[1];	// q = -2;
	charge_type_z4[2] = vertex_average[3];	// q = +2;
	charge_type_z4[3] = vertex_average[0];	// q = -4;
	charge_type_z4[4] = vertex_average[4];	// q = +4;

	charge_type_z3[0] = vertex_average[9] + vertex_average[12] + vertex_average[13]; // q = +1;
	charge_type_z3[1] = vertex_average[8];	// q = -1;
	charge_type_z3[2] = vertex_average[10];	// q = +3;
	charge_type_z3[3] = vertex_average[7];	// q = -3;
	charge_type_z3[4] = vertex_average[11];	// 0. can't happen

	for (i = 0; i < 5; i++)
		stat << (double) charge_type_z4[i] / (double) stat_steps / (double) N_vertex_z_type_4 << " ";
	for (i = 0; i < 5; i++)
		stat << (double) charge_type_z3[i] / (double) stat_steps / (double) N_vertex_z_type_3 << " ";

	total_charge_on_z4 = charge_type_z4[1] * (-2) + charge_type_z4[2] * 2 
					   + charge_type_z4[3] * (-4) + charge_type_z4[4] * 4;
	total_charge_on_z3 = charge_type_z3[0] - charge_type_z3[1]
					   + charge_type_z3[2] * (-3) + charge_type_z3[3] * 3;
	stat << (double) total_charge_on_z4 / (double) stat_steps / (double) N_vertex_z_type_4 << " ";
	stat << (double) total_charge_on_z3 / (double) stat_steps / (double) N_vertex_z_type_3 << " ";

	stat << (double) N_z4_charged / (double) stat_steps / (double) N_vertex_z_type_4 << " ";
	stat << (double) N_z3_charged / (double) stat_steps / (double) N_vertex_z_type_3 << " ";

	if (N_z3_charged)
	{
		stat << (double) charge_3_near_3 / (double) N_z3_charged / 3.0 << " ";
		stat << (double) charge_4_near_3 / (double) N_z3_charged / 3.0 << " ";
	} else {
		stat << "0 0 ";
	}

	if (N_z4_charged)
	{
		stat << (double) charge_3_near_4 / (double) N_z4_charged / 4.0 << " ";
		stat << (double) charge_4_near_4 / (double) N_z4_charged / 4.0 << " ";
	} else {
		stat << "0 0 ";
	}

	for (i = 0; i < N; i++)
		total_flips += spin_flip[i];
	stat << (double) total_flips / (double) stat_steps << endl;
	zerosStatistics();
	stat.close();
}

int main(int argc, char* argv[]) 
{
	cout << "Generic BD simulation" << endl;
	cout << "Simulating spontaneous lane formation" << endl;

	FILE* f;
	moviefile = new char[30];
	statfile = new char[30];
	if (argc == 2) 
		readDataFromFile(argv[1]);
	initData();
	// erasing stat file
	ofstream stat(statfile);
	stat << "";
	stat.close();

	f = fopen(moviefile, "wb");

	initSquarePinningSites();
	initSquareVerteces();
	zeros(vertex_type, N_vertex);
	zeros(vertex_q, N_vertex);
	zerosStatistics();
	decimating();
	initParticles();
	writeGfile();
	cout << "Sx: " << Sx << endl;
	cout << "Sy: " << Sy << endl;
	cout << "N: " << N << endl;

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
		}

		calculatePairwiseForces();
		calculateThermalForce();
		calculateModifiedPinningiteForces();
		moveParticles();

		calculateStatisticsPerStep();

		if (current_time % 100 == 0)
        {
            calculateVertexType();
			calculateTopologicalChargeOnVertex();
			writeCmovie(f, current_time);
            multiplier = (double) current_time / (double) total_time * max_multiplier;
            writeStatistics();
        }
	}
	fclose(f);
	freeData();
	stop_timing();
	return 0;
}