#include <stdio.h>

struct Coordinate {
    double x, y;
};

typedef struct Coordinate Coordinate;

struct Particle {
    int id;         //id of a particle
    int color;      //to distinguish it the particles
    Coordinate coord;// double x, y;    //x,y coordinate
    double fx, fy;  //fx,fy forces acting on the particle
    double q;       //toltes
    double drx;
    double dry;
    int pinningSiteId; //melyik pinning sitehoz
    int mortonId;       //morton kodok
};

typedef struct Particle Particle;