#include <stdio.h>

struct Coordinate {
    double x, y;
};

typedef struct Coordinate Coordinate;

struct Particle {
    u_short id;             //id of a particle
    u_char color;           //to distinguish it the particles
    Coordinate coord;       // double x, y;    //x,y coordinate
    double fx, fy;          //fx,fy forces acting on the particle
    float q;                //toltes
    u_short pinningSiteId;  //melyik pinning sitehoz
    u_short mortonId;       //morton kodok
};

typedef struct Particle Particle;