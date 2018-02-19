#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "particle.h"

double sX, sY;          //system size x,y direction
double sX2, sY2;        //half of the system size
int N;                  //number of particles
double dt;              //length of a single time step
Particle *particles;    //list of particle
int t;                  //time

//time measuring
time_t time_begin;
time_t time_end;

void initParticles(int nrParticles, double systemSize, double timeStep) {
    N = nrParticles;
    sX = sY = systemSize;
    sX2 = sY2 = systemSize / 2;
    dt = timeStep;
    particles = (Particle *) malloc(N * sizeof(Particle));
    if (!particles) {
        perror("Allocation problem!");
        exit(EXIT_FAILURE);
    }
    int i = 0, overlap;
    double dx, dy, tmpX, tmpY;

    while (i < N) {
        //printf("kor: %d\n", i);
        particles[i].id = i;
        do {
            tmpX = sX * rand() / (RAND_MAX + 1.0);
            tmpY = sY * rand() / (RAND_MAX + 1.0);
            for (int j = 0; j < i; j++) {
                dx = tmpX - particles[j].coord.x;
                dy = tmpY - particles[j].coord.y;
                //PBC fold back
                if (dx >= sX2) dx -= sX;
                if (dx < -sX2) dx += sX;
                if (dy >= sY2) dy -= sY;
                if (dy < -sY2) dy += sY;
                overlap = (sqrt(dx * dx + dy * dy) < 0.2); //checking if already taken
            }
        } while (overlap == 1); //regenerate until, the  coordinate is unused

        particles[i].color = ((rand() / (RAND_MAX + 1.0)) > 0.5); //color
        particles[i].coord.x = tmpX;
        particles[i].coord.y = tmpY;
        //particles[i].x = tmpX; //particle coordinate
        //particles[i].y = tmpY;
        i++;
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
    }
}

void start(){
    printf("Let's do it!");
    srand(time(NULL));
    time(&time_begin);
    printf("Program started: %s\n",asctime(localtime(&time_begin)));
}

void end(){
    time(&time_end);
    printf("Program started: %s\n",asctime(localtime(&time_begin)));
    printf("Program ended: %s\n",asctime(localtime(&time_end)));
    printf("Program ran: %lf seconds\n",difftime(time_end,time_begin));
}

int main() {
    start();
    Particle *particle = (Particle *) malloc(sizeof(Particle));
    (*particle).color = 4;
    (*particle).coord.x = 4;
    (*particle).coord.y = 2*4;
    particle->fx = 3;
    printf("%.2f\n", sqrt(particle->fx));
    printf("%.2f\n", particle->coord.x);
    printf("%.2f\n", particle->coord.y);

    initParticles(20, 20.0, 0.002);
    printf("N=%d sX= %.2f Sy=%.2f sX= %.2f Sy=%.2f  dt=%.4f\n", N, sX, sY, sX2, sY2, dt);
    particles[2].fx = 4;
    printf("tttttt%.2f\n", particles[2].fx);

//    int i = 0;
//    for (i = 0; i < 2; i++) {
//        printf("Generalt szam %d\n", rand() % 20);
//    }
//    int i=0;
//    while(i<1000000){
//        printf("Generalt szam %d\n", i++);
//    }

    printf("Hello, World!\n");
   // free(particles);
    //free(particle);
//     for(i=0; i<N;i++){
//            free(particles[i]);
//     }
    end();
    return 0;
}