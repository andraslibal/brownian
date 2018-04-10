//
// Created by rebeka on 27.03.2018.
//
#include "particle.h"


struct Pinning{
    int id;             //id of a pinning
    Coordinate coord;   //x,y coordinate
    double r;           //kor sugara
    double fMax;       //
    int particlesId; //melyik reszecske tartozik hozza, max 1

    double lx,ly;       //koztes resz hossza,szelessege
    double sinfi,cosfi; //mennyire van megdolve
    double middleHeight;//godrok kozti magassag
};

typedef struct Pinning Pinning;
