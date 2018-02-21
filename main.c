#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "thpool.h"
#include "particle.h"

/*(N, number of  particles and also the box size, proportionally! if  you
increase the box to be 2x as big in x and y direction N has to be 4x as much to have
the same density*/
double sX, sY;          //system size x,y direction
double sX2, sY2;        //half of the system size
int N;                  //number of particles
double dt;              //length of a single time step
Particle *particles;    //list of particle
int t;                  //time

//thread
int interval;

//time measuring, how much the simulation ran
time_t time_begin;
time_t time_end;

//statistic file, file with coordinates of the  particles
FILE *statistics_file;
FILE *moviefile;

void keepParticleInSystem(double *dx, double *dy) {
    if (*dx > sX2) *dx -= sX;
    if (*dx < -sX2) *dx += sX;
    if (*dy > sY2) *dy -= sY;
    if (*dy < -sY2) *dy += sY;
}

void stowParticle(Coordinate *right_coord, int i) {
    Coordinate tmp;
    double dx, dy;
    int overlap = 0;
    do {
        tmp.x = sX * rand() / (RAND_MAX + 1.0);
        tmp.y = sY * rand() / (RAND_MAX + 1.0);
        for (int j = 0; j < i; j++) {
            dx = tmp.x - particles[j].coord.x;
            dy = tmp.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy); //back to the "box"
            overlap = (sqrt(dx * dx + dy * dy) < 0.2); //checking if the position is already taken
        }
    } while (overlap == 1); //regenerate until, the  coordinate is unused

    *right_coord = tmp;
}

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
    int i = 0;
    Coordinate position;

    while (i < N) {
        particles[i].id = i;
        particles[i].color = ((rand() / (RAND_MAX + 1.0)) > 0.5); //color
        stowParticle(&position, i);  //particle coordinate
        particles[i].coord.x = position.x;
        particles[i].coord.y = position.y;
        i++; //miert csak a kovetkezonek kell allitani az erot?
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
    }
}

void calculateExternalForces(int *multipler) {
    int i;

    int thread_from = (*multipler * interval) - interval;
    int thread_to = *multipler * interval;

    for (i = thread_from; i < thread_to; i++) {
        if (particles[i].color) {
            particles[i].fx -= 0.5;
        } else {
            particles[i].fx += 0.5;
        }
    }
}

void calculatePairwiseForces(int *multipler) {
    int i, j;
    double dx, dy, dr, dr2;
    double f, fx, fy;

    int thread_from = (*multipler * interval) - interval;
    int thread_to = *multipler * interval;

    for (i = thread_from; i < thread_to - 1; i++) {
        for (j = i + 1; j < thread_to; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;
            dr = sqrt(dr2);

            (dr < 0.2) ? (f = 100.0, printf("Warning!!!dr%f\n", dr)) : (f = 1 / dr2 * exp(-0.25 * dr));

            //project it to the axes get the fx, fy components
            fx = f * dx / dr;
            fy = f * dy / dr;

            particles[i].fx += fx;
            particles[i].fy += fy;

            particles[j].fx -= fx;
            particles[j].fy -= fy;

        }
    }
}

void moveParticles() {
    double dx, dy;
    for (int i = 0; i < N; i++) {
        dx = particles[i].fx * dt;
        dy = particles[i].fy * dt;
        particles[i].coord.x += dx;
        particles[i].coord.y += dy;
        //out of box
        if (particles[i].coord.x < 0) particles[i].coord.x += sX;
        if (particles[i].coord.y < 0) particles[i].coord.y += sY;
        if (particles[i].coord.x > sX) particles[i].coord.x -= sX;
        if (particles[i].coord.y > sY) particles[i].coord.y -= sY;

        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
    }
}

//this is for plot
void write_cmovie() {
    int i;
    float floatholder;
    int intholder;

    intholder = N;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    intholder = t;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    for (i = 0; i < N; i++) {
        intholder = particles[i].color + 2; //miert kell +2
        fwrite(&intholder, sizeof(int), 1, moviefile);
        intholder = particles[i].id;//ID
        fwrite(&intholder, sizeof(int), 1, moviefile);
        floatholder = (float) particles[i].coord.x;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = (float) particles[i].coord.y;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, moviefile);
    }

}

void write_statistics() {
    double avg_vx = 0.0;
    for (int i = 0; i < N; i++) {
        avg_vx += particles[i].fx;
    }

    avg_vx = avg_vx / (double) N;

    fprintf(statistics_file, "%d %f\n", t, avg_vx);

}

void start() {
    printf("\tLet's do it!\n");
    //  srand(time(NULL));
    time(&time_begin);
}

void end() {
    time(&time_end);
    printf("Program started: %s\n", asctime(localtime(&time_begin)));
    printf("Program ended: %s\n", asctime(localtime(&time_end)));
    printf("Program ran: %f seconds\n", difftime(time_end, time_begin));
}

int main() {
    start();
    initParticles(400, 20.0, 0.002);
    const int rate = 16;
    int inters[rate];

    for (int x = 0; x < rate ; x++) {
        inters[x] = x + 1;
    }
    interval = N / rate;

    moviefile = fopen("result.mvi", "w");
    statistics_file = fopen("statistics.txt", "wt");

    threadpool thpool = thpool_init(4);

    for (t = 0; t < 100000; t++) {
        for (int x = 0; x < rate; x++) {
            thpool_add_work(thpool, (void *) calculatePairwiseForces, &inters[x]);
        }

        thpool_wait(thpool);

        for (int x = 0; x < rate; x++) {
            thpool_add_work(thpool, (void *) calculateExternalForces, &inters[x]);
        }

        thpool_wait(thpool);

        write_statistics();
        moveParticles();
        if (t % 100 == 0) {
            write_cmovie();
        }
        if (t % 500 == 0) {
            printf("time = %d\n", t);
            fflush(stdout);
        }
    }
    fclose(statistics_file);
    fclose(moviefile);
    end();
    free(particles);
    thpool_destroy(thpool);
    return 0;
}