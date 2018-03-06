#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <regex.h>
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
    double dx, dy, dr;
    int overlap, j;
    do {
        overlap = 0;
        tmp.x = sX * rand() / (RAND_MAX + 1.0);
        tmp.y = sY * rand() / (RAND_MAX + 1.0);
        for (j = 0; j < i; j++) {
            dx = tmp.x - particles[j].coord.x;
            dy = tmp.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy); //back to the "box"
            dr = sqrt(dx * dx + dy * dy);
            if (dr<0.2) { //checking if the position is already taken
                overlap = 1;
                break;
            }
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
        particles[i].q = particles[i].color ? -1.0:1.0;
        stowParticle(&position, i);  //particle coordinate
        particles[i].coord.x = position.x;
        particles[i].coord.y = position.y;
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
        i++;

    }

    //particles initialized
    FILE *teszt;
    teszt=fopen("reszecskek.txt","wt");
    for(int i=0;i<N;i++)
        fprintf(teszt,"%lf %lf\n",particles[i].coord.x,particles[i].coord.y);
    fclose(teszt);

}

void calculateExternalForces() {
    int i;
    for (i = 0; i < N; i++) {
            particles[i].fx += particles[i].q * 0.5;
    }
}

void calculatePairwiseForces() {
    int i, j;
    double dx, dy, dr, dr2;
    double f, fx, fy;

    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;
            dr = sqrt(dr2);

            (dr < 0.2) ? (f = 100.0, printf("Warning!!!dr%f particles %d(%lf %lf) %d(%lf %lf)\n", dr,i,particles[i].coord.x,particles[i].coord.y,j,particles[j].coord.x,particles[j].coord.y)) : (f = 1 / dr2 * exp(-0.25 * dr));

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
        intholder = particles[i].color + 2;
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
        if (particles[i].color==0) avg_vx += particles[i].fx;
        if (particles[i].color==1) avg_vx -= particles[i].fx;
    }

    avg_vx = avg_vx / (double) N;

    fprintf(statistics_file, "%d %f\n", t, avg_vx);

}

void start() {
    printf("\tLet's do it!\n");
    time(&time_begin);
}

void end() {
    time(&time_end);
    double timeDiff = 0.0;
    printf("Program started: %s\n", asctime(localtime(&time_begin)));
    printf("Program ended: %s\n", asctime(localtime(&time_end)));
    timeDiff = difftime(time_end, time_begin);
    printf("Program ran: %.2f seconds, %.2f minutes.\n", timeDiff, timeDiff / 60);
}

//return 1 if ok,
int properFilename(char *filename) {
    regex_t regexCompiled;
    char const *PATTERN = "[a-zA-Z]+([a-z0-9A-Z_])*";
    int result;
    size_t maxMatches = 1; //Is the number of matches allowed.
    size_t maxGroups = 1;
    regmatch_t pmatch[maxGroups]; //When maxMatches is non-zero, points to an array with at least maxMatches elements.

    if (regcomp(&regexCompiled, PATTERN, REG_EXTENDED)) {
        printf("Could not compile regular expression. regcomp() failed, returning nonzero\n");
        exit(1);
    }

    //If successful, the regexec() function returns zero to indicate that string matched PATTERN
    result = !regexec(&regexCompiled, filename, maxMatches, pmatch, 0) ? (!(pmatch[maxGroups - 1].rm_so) &&
                                                                          (pmatch[maxGroups - 1].rm_eo ==
                                                                           strlen(filename))) : 0;

    regfree(&regexCompiled);
    return result;
}

int main(int argc, char *argv[]) {
    char filename[100] = {0};
    if (argc > 1) {
        char *inputParameter = argv[1];
        strncpy(filename, properFilename(inputParameter) ? strcat(inputParameter, ".mvi") : "result.mvi",
                sizeof(filename));
    } else {
        strncpy(filename, "result.mvi", sizeof(filename));
    }

    start();
    initParticles(200, 20.0, 0.002);

    moviefile = fopen(filename, "wb");
    statistics_file = fopen("statistics.txt", "wt");
    if (!moviefile || !statistics_file) {
        fprintf(stderr, "Failed to open file.\n");
        return 1;
    }
    for (t = 0; t < 100000; t++) {
        calculatePairwiseForces();
        calculateExternalForces();
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
    printf("The result(for plot) can be found in file: %s\n", filename);
    free(particles);
    return 0;
}