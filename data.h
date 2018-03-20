#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <regex.h>
#include <stdio.h>

/*(N, number of  particles and also the box size, proportionally! if  you
        increase the box to be 2x as big in x and y direction N has to be 4x as much to have
the same density*/
double sX, sY;          //system size x,y direction
double sX2, sY2;        //half of the system size
int N;                  //number of particles
double dt;              //length of a single time step
Particle *particles;    //list of particle
int t;                  //time
//Cutoff
double r0;              //the cutoff distance
double rv;              //the verlet cutoff distance
double rvminr02;         //(rv-r0)^2
//Verlet lista
int N_verlet_list;      //szomszedok szama
int *vlist1;            //i
int *vlist2;            //j
int flag_to_rebuild_verlet; //ha ujra kell epiteni a szomszedsagi listat
int N_verlet_actual;    //hany szomszed van
//Tabulated force
int N_tabulated;        //hany f/r-et szamitunk ki
double *tabulated_f_per_r;//f/r ertekeinek eltarolasa
double tabulate_start;    //hol kezdodik
double tabulate_step;     //leptek

//time measuring, how much the simulation ran
time_t time_begin;
time_t time_end;

//statistic file, file with coordinates of the  particles
FILE *statistics_file;
FILE *moviefile;
