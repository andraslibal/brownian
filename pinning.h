//
// Created by rebeka on 27.03.2018.
//
#include <stdio.h>
#include "particle.h"


struct Pinning{
    int id;             //id of a pinning
    Coordinate coord;   //x,y coordinate
    double r;           //kor sugara
    double f_max;       //
};

typedef struct Pinning Pinning;
