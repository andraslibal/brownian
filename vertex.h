//
// Created by rebeka on 15.05.2018.
//
//#include "particle.h"

struct Vertex{
    int id;             //id of vertex
    Coordinate coord;   //helyzete, x,y koord
    int type;           //tipus
    int chessColor;     //a GS megkulonboztetesere, sakktabla szeruen helyezkednek el(fekete/feher)
    int color;          //szin tipus szerint
    int particles[4];   //hozza tartozo reszecskek
    int nearparticles[4];//eltaroljuk hogy a reszecskek kozel vannak vagy tavol
    int pinningsites[4];//hozza tartozo pinnintSiteok
    int clusterId;      //melyik clusterhez tartozik
};

typedef struct Vertex Vertex;