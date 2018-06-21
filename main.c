#include "data.h"
#include <errno.h>
#include <stdint.h>

uint64_t u, v, w;
double norm_num;
int norm_saved;

void setseed(uint64_t seedmodifiervalue);

double Rand(void);

double gasdev(void);

uint64_t int64();

/***************************************************
	int64_t
***************************************************/
inline uint64_t int64() {
    uint64_t x;
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    x = u ^ (u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x + v) ^ w;
}

/***************************************************
	setSeed
***************************************************/
inline void setseed(uint64_t seedmodifiervalue) {
    v = 4101842887655102017LL;
    w = 1;
    u = seedmodifiervalue ^ v;
    int64();
    v = u;
    int64();
    w = v;
    int64();
    norm_saved = 0;
}

/***************************************************
	uniform distribution
***************************************************/\


inline double Rand() {
    return 5.42101086242752217E-20 * int64();
}

/***************************************************
	gaussian distribution
***************************************************/
inline double gasdev() {
    double x1, x2, w;
    if (!norm_saved) {
        do {
            x1 = 2.0 * Rand() - 1.0;
            x2 = 2.0 * Rand() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0 || w == 0);
        w = sqrt((-2.0 * log(w)) / w);
        norm_num = x1 * w;
        norm_saved = 1;
        return (x2 * w);
    } else {
        norm_saved = 0;
        return norm_num;
    }
}


//elore kiszamolt f/r, hogy a reszecskek egymasra hataskor ne kelljen szamolni, csak kiszedni a listabol az erteket
void tabulateForces() {
    double x_min = 0.1, x_max = rv, x2, x;
    double f = 0;

    tabulate_start = x_min * x_min;
    tabulate_step = (x_max * x_max - x_min * x_min) / (N_tabulated - 1.0);

    for (int i = 0; i < N_tabulated; i++) {
        x2 = i * tabulate_step + tabulate_start;
        x = sqrt(x2);
        f = 1 / x2 * exp(-(1 / r0) * x);
        tabulated_f_per_r[i] = f / x;
    }

}

void initVertexStat() {
    avgNCluster = 0;
    avgSmallCluster = 0;
    avgBigCluster = 0;
    avgAvgCluster = 0.0;
    for (int j = 0; j < sumVertexTypes; j++)
        vertexTypes[j] = 0;
}


//alapbeallitasok, Particle,system,verlet,tabulated force, temperature,multiply, Vertexszam
void init(int nrParticles, double systemSize, double timeStep, double cutOff, double verletCutOff, double temp,
          double mult, double maxmult) {
    N_particles = nrParticles;
    particles = NULL;
    particles = (Particle *) realloc(particles, N_particles * sizeof(Particle));
    if (!particles) {
        perror("Allocation problem! Particles");
        exit(EXIT_FAILURE);
    }
    sX = pinningNX * (pinningDistX);
    sY = pinningNY * (pinningDistY);
    sX2 = sX / 2;
    sY2 = sY / 2;
    dt = timeStep;

    //Cutoff distance
    r0 = cutOff;
    rv = verletCutOff;
    rvminr02 = (rv - r0) * (rv - r0);

    //mozgatas-homerseklet
    temperature = temp;
    multiply = mult;
    maxMultiply = maxmult;
    printf("Temperature=%.2lf multiply=%.2lf\n", temperature, multiply);

    //Verlet list
    ////N_verlet= (pi*rv^2*N^2)/(2*sX*sY)
    N_verlet_list = (int) ((3.14 * rv * rv * N_particles * N_particles) / (2 * sX * sY));

    vlist1 = NULL;
    vlist2 = NULL;
    vlist1 = (int *) realloc(vlist1, N_verlet_list * sizeof(int));
    vlist2 = (int *) realloc(vlist2, N_verlet_list * sizeof(int));
    if (vlist1 == NULL || vlist2 == NULL) {
        perror("Allocation problem with Verlet list, exiting");
        exit(EXIT_FAILURE);
    }
    flag_to_rebuild_verlet = 0;

    //Tabulated force
    N_tabulated = 50000;
    tabulated_f_per_r = NULL;
    tabulated_f_per_r = (double *) realloc(tabulated_f_per_r, N_tabulated * sizeof(double));
    if (!tabulated_f_per_r) {
        perror("Allocation problem! Tabulated force(f/r) list.");
        exit(EXIT_FAILURE);
    }

    tabulateForces();

    ///Vertex
    N_vertex = pinningNX * pinningNY;
    vertex = (Vertex *) malloc(N_vertex * sizeof(Vertex));
    if (!vertex) {
        perror("Allocation problem! Vertex");
        exit(EXIT_FAILURE);
    }

    ///Cluster
    N_clusters = 0;
    clusters = (int *) malloc(N_clusters * sizeof(int));
    if (clusters == NULL) {
        perror("Allocation problem with Cluster list, exiting");
        exit(EXIT_FAILURE);
    }

    ///Statistics
    sumVertexTypes = 7;
    vertexTypes = (int *) malloc(sumVertexTypes * sizeof(int));
    if (clusters == NULL) {
        perror("Allocation problem with VertexTypes list, exiting");
        exit(EXIT_FAILURE);
    }
    initVertexStat();
    N_average = 0;


    printf("General Init done!\n N_particles=%d N_pinning=%d, N_vertex=%d, sX=%.2f sY=%.2f\n", N_particles, N_pinning,
           N_vertex, sX, sY);

}

void initPinning(int nrX, int nrY, double distX, double distY, double lx2, double ly2, double r, double fmax,
                 double middleHeight) {
    //Pinnings
    N_pinning = nrX * nrY * 2;
    pinnings = NULL;
    pinnings = (Pinning *) realloc(pinnings, N_pinning * sizeof(Pinning));
    if (!pinnings) {
        perror("Allocation problem! PinningSites");
        exit(EXIT_FAILURE);
    }

    pinningNX = nrX;
    pinningNY = nrY;
    pinningDistX = distX;
    pinningDistY = distY;
    pinningLx2 = lx2;
    pinningLy2 = ly2;
    pinningR = r;
    pinningFMax = fmax;
    pinningMiddleHeight = middleHeight;
    pinningK = pinningFMax / pinningR;
}


uint64_t generateMort(uint32_t x, uint32_t y) {
    uint64_t key = 0;
    int level = 0, left_bit = 0, right_bit = 0, value = 0;

    while ((x > 0) || (y > 0)) {
/*   split off the row (left_bit) and column (right_bit) bits and
     then combine them to form a bit-pair representing the
     value                                                  */

        left_bit = x % 2;
        right_bit = y % 2;
        value = right_bit + 2 * left_bit;

        key += value << (2 * level);

        x /= 2;
        y /= 2;
        level++;

    }
    return (key);

}

//ha kilepett a dobozbol, visszarakjuk
void keepParticleInSystem(double *dx, double *dy) {
    if (*dx > sX2) *dx -= sX;
    if (*dx < -sX2) *dx += sX;
    if (*dy > sY2) *dy -= sY;
    if (*dy < -sY2) *dy += sY;
}

//reszecskek helyenek a megkeresese, ne fedjek le egymast es maradjon a dobozban
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
            if (dr < 0.2) { //checking if the position is already taken
                overlap = 1;
                break; //ha mar egy kozel kerult, a tobbit nem szamoljuk hanem uj koordinatat generalunk neki
            }
        }
    } while (overlap == 1); //regenerate until, the  coordinate is unused

    *right_coord = tmp;
}


void initParticles() {
    int i = 0;
    Coordinate position;

    while (i < N_particles) {
        particles[i].id = (u_short) i;
        particles[i].color = (u_char) ((rand() / (RAND_MAX + 1.0)) > 0.5); //color
        particles[i].q = particles[i].color ? -1.0 : 1.0;
        stowParticle(&position, i);  //particle coordinate
        particles[i].coord.x = position.x;
        particles[i].coord.y = position.y;
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
        particles[i].pinningSiteId = 0;
        i++;
    }
}

void initSquareParticles() {
    for (int i = 0; i < N_pinning; i++) {
        particles[i].id = (u_short) i;
        particles[i].color = 1; //color
        particles[i].q = 1;
        particles[i].coord.x = pinnings[i].coord.x;
        particles[i].coord.y = pinnings[i].coord.y;
        particles[i].coord.x += pinnings[i].lx * pinnings[i].cosfi;
        particles[i].coord.y += pinnings[i].lx * pinnings[i].sinfi;
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
        particles[i].mortonId = (u_short) generateMort((uint32_t) particles[i].coord.x,
                                                       (uint32_t) particles[i].coord.y);
        particles[i].pinningSiteId = pinnings[i].id;
        pinnings[i].particlesId = particles[i].id;
    }
}

void initSquareParticlesRandom() {

    for (int i = 0; i < N_pinning; i++) {
        particles[i].id = (u_short) i;
        particles[i].color = 1; //color
        particles[i].q = 1;
        particles[i].coord.x = pinnings[i].coord.x;
        particles[i].coord.y = pinnings[i].coord.y;
        if (Rand() < 0.5) {
            particles[i].coord.x += (-1) * pinnings[i].lx * pinnings[i].cosfi;
            particles[i].coord.y += (-1) * pinnings[i].lx * pinnings[i].sinfi;
        } else {
            particles[i].coord.x += pinnings[i].lx * pinnings[i].cosfi;
            particles[i].coord.y += pinnings[i].lx * pinnings[i].sinfi;
        }

        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
        particles[i].mortonId = (u_short) generateMort((uint32_t) particles[i].coord.x,
                                                       (uint32_t) particles[i].coord.y);
        particles[i].pinningSiteId = pinnings[i].id;
        pinnings[i].particlesId = particles[i].id;
    }
}

void initSquareVertex() {
    u_short k = 0; // hanyadik vertexnel tartunk
    int pinningNX2 = 2 * pinningNX;
    int i0, i1, i2, i3;

    for (int i = 0; i < pinningNX; i++)
        for (int j = 0; j < pinningNY; j++) {

            vertex[k].id = k;
            vertex[k].coord.x = (i + 0.1) * pinningDistX;
            vertex[k].coord.y = (j + 0.1) * pinningDistY;
            vertex[k].chessColor = (i + j) % 2;
            vertex[k].type = 0;
            vertex[k].color = 2;
            vertex[k].clusterId = -1;

            //hozza tartozo reszecske es pinningSite indexenek meghatarozasa

            if (j == 0)
                i0 = (i + 1) * pinningNX2 - 1;
            else
                i0 = i * pinningNX2 + j * 2 - 1;
            if (i == 0)
                i1 = (pinningNX - 1) * pinningNX2 + j * 2;
            else
                i1 = i * pinningNX2 - pinningNX2 + j * 2;
            i2 = i * pinningNX2 + j * 2;
            i3 = i * pinningNX2 + j * 2 + 1;

            for (int l = 0; l < 4; l++) {
                vertex[k].nearparticles[l] = 0;
            }
            //hozza tartozo reszecskek
            vertex[k].particles[0] = pinnings[i0].particlesId;
            vertex[k].particles[1] = pinnings[i1].particlesId;
            vertex[k].particles[2] = pinnings[i2].particlesId;
            vertex[k].particles[3] = pinnings[i3].particlesId;

            //hozza tartozo pinningSiteok
            vertex[k].pinningsites[0] = pinnings[i0].id;
            vertex[k].pinningsites[1] = pinnings[i1].id;
            vertex[k].pinningsites[2] = pinnings[i2].id;
            vertex[k].pinningsites[3] = pinnings[i3].id;

            k++;

        }
}

void printVertex() {
    for (int i = 0; i < N_vertex; i++) {
        printf("id %d,szin %d,pin0 %d pin1 %d pin2 %d pin3 %d, part0 %d part1 %d part2 %d part3 %d mort0 %d mort1 %d mort2 %d mort3 %d\n",
               vertex[i].id, vertex[i].color, vertex[i].pinningsites[0], vertex[i].pinningsites[1],
               vertex[i].pinningsites[2], vertex[i].pinningsites[3], vertex[i].particles[0],
               vertex[i].particles[1], vertex[i].particles[2], vertex[i].particles[3],
               particles[vertex[i].particles[0]].mortonId, particles[vertex[i].particles[1]].mortonId,
               particles[vertex[i].particles[2]].mortonId, particles[vertex[i].particles[3]].mortonId);
    }
}

void printPinning() {
    for (int i = 0; i < N_pinning; i++) {
        printf("id %d, c: %.2f %.2f,r %.2f, f %.2f, id: %d, m:%.2f cos:%.2f, szin:%.2f,lx: %.2f ly: %.2f\n",
               pinnings[i].id, pinnings[i].coord.x,
               pinnings[i].coord.y, pinnings[i].r, pinnings[i].fMax, pinnings[i].particlesId, pinnings[i].middleHeight,
               pinnings[i].cosfi, pinnings[i].sinfi, pinnings[i].lx, pinnings[i].ly);
    }
}

void printParticle() {
    for (int i = 0; i < N_particles; i++) {
        printf("i %d id %d,%d, %.2f %.2f color %d pinningID %d pinPart %d\n", i, particles[i].id, particles[i].mortonId,
               particles[i].coord.x, particles[i].coord.y, particles[i].color, particles[i].pinningSiteId,
               pinnings[particles[i].pinningSiteId].particlesId);
    }
}

//pinning helyenek a megkeresese
void stowPinning(Coordinate *right_coord, int i) {
    Coordinate tmp;
    double dx, dy, dr;
    int overlap, j;
    do {
        overlap = 0;
        tmp.x = sX * rand() / (RAND_MAX + 1.0);
        tmp.y = sY * rand() / (RAND_MAX + 1.0);
        for (j = 0; j < i; j++) {
            dx = tmp.x - pinnings[j].coord.x;
            dy = tmp.y - pinnings[j].coord.y;
            keepParticleInSystem(&dx, &dy); //back to the "box"
            dr = sqrt(dx * dx + dy * dy);
            if (dr < 2.2) { //checking if the position is already taken
                overlap = 1;
                break; //ha mar egy kozel kerult, a tobbit nem szamoljuk hanem uj koordinatat generalunk neki
            }
        }
    } while (overlap == 1); //regenerate until, the  coordinate is unused

    *right_coord = tmp;
}

void initPinningSites() {
    u_short k = 0;

    for (int i = 0; i < pinningNX; i++)
        for (int j = 0; j < pinningNY; j++) {
            //horizontal
            pinnings[k].id = k;
            pinnings[k].coord.x = (i + 0.5 + 0.1) * pinningDistX;
            pinnings[k].coord.y = (j + 0.1) * pinningDistY;
            pinnings[k].lx = pinningLx2;
            pinnings[k].ly = pinningLy2;
            pinnings[k].sinfi = 0.0;
            pinnings[k].cosfi = 1.0;
            pinnings[k].middleHeight = pinningMiddleHeight;
            pinnings[k].fMax = pinningFMax;
            pinnings[k].particlesId = 0;
            pinnings[k].r = pinningR;

            //vertical
            k++;
            pinnings[k].id = k;
            pinnings[k].coord.x = (i + 0.1) * pinningDistX;
            pinnings[k].coord.y = (j + 0.5 + 0.1) * pinningDistY;
            pinnings[k].lx = pinningLx2;
            pinnings[k].ly = pinningLy2;
            pinnings[k].sinfi = 1.0;
            pinnings[k].cosfi = 0.0;
            pinnings[k].middleHeight = pinningMiddleHeight;
            pinnings[k].fMax = pinningFMax;
            pinnings[k].particlesId = 0;
            pinnings[k].r = pinningR;
            k++;
        }
}

//meghatarozza, hogy a ket GS(5/6 tipusu a vertex) kozul melyikhez tartozik
void calculateVertexTypeGS(int i) {

    if (vertex[i].nearparticles[0] == 1 && vertex[i].nearparticles[3] == 1)
        vertex[i].type = 5 + vertex[i].chessColor;
    if (vertex[i].nearparticles[1] == 1 && vertex[i].nearparticles[2] == 1)
        vertex[i].type = 6 - vertex[i].chessColor;
}

//vertex tipusok
void calculateVertexTypes() {
    double dx, dy, dr2;

    for (int i = 0; i < N_vertex; i++) {
        vertex[i].type = 0;
        for (int j = 0; j < 4; j++) {
            vertex[i].nearparticles[j] = 0;
            dx = vertex[i].coord.x - particles[vertex[i].particles[j]].coord.x;
            dy = vertex[i].coord.y - particles[vertex[i].particles[j]].coord.y;
            dx = dx * pinnings[vertex[i].pinningsites[j]].cosfi;
            dy = dy * pinnings[vertex[i].pinningsites[j]].sinfi;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;
            if (dr2 < ((pinningDistX * pinningDistX) / 4.0)) {
                vertex[i].type++;
                vertex[i].nearparticles[j] = 1;
                particles[vertex[i].particles[j]].color = 3;
            }
        }

        if (vertex[i].type == 2)
            calculateVertexTypeGS(i);
        vertex[i].color = vertex[i].type + 4;
    }
}

void calculateThermalForces() {
    double fx, fy;
    for (int i = 0; i < N_particles; i++) {
        fx = temperature * gasdev();
        fy = temperature * gasdev();
        particles[i].fx += fx;
        particles[i].fy += fy;
    }
}

void calculatePinningForces() {
    double dx, dy, fx, fy;
    double dxp, dyp, fxp = 0, fyp = 0;
    int j;

    for (int i = 0; i < N_particles; i++) {
        j = i;
        dx = particles[i].coord.x - pinnings[j].coord.x;
        dy = particles[i].coord.y - pinnings[j].coord.y;
        keepParticleInSystem(&dx, &dy);
        //forgatas vizszintesbe
        dxp = pinnings[j].cosfi * dx - pinnings[j].sinfi * dy;
        dyp = pinnings[j].sinfi * dx + pinnings[j].cosfi * dy;

        if (dxp <= (-pinningLx2)) {
            fxp = -pinningK * (dxp + pinningLx2);
            fyp = -pinningK * dyp;
//            elso eset dxp<=-lx2\n";
        } else if (dxp > (-pinningLx2) && (dxp < pinningLx2)) {
            fxp = (pinningLx2 - fabs(dxp)) * pinningMiddleHeight;
            if (dxp < 0) fxp = -fxp;
//            fxp=0.0;
            fyp = -pinningK * dyp;
//            "masodik eset -lx2<dxp<lx2\n";
        } else if (dxp >= pinningLx2) {
            fxp = -pinningK * (dxp - pinningLx2);
            fyp = -pinningK * dyp;
//            "harmadik eset dxp>=lx2\n";
        }

        fx = pinnings[j].cosfi * fxp + pinnings[j].sinfi * fyp;
        fy = -pinnings[j].sinfi * fxp + pinnings[j].cosfi * fyp;
        particles[i].fx += fx;
        particles[i].fy += fy;
    }

}

void calculateExternalForces() {
    for (int i = 0; i < N_particles; i++) {
        particles[i].fx += (double) t / 10000 * 5.0;
    }
}

void calculatePairwiseForces() {
    int i, j;
    double dx, dy, dr2;
    double fx, fy;
    int tabulatedForceIndex;

    for (i = 0; i < N_particles - 1; i++) {
        for (j = i + 1; j < N_particles; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;

            tabulatedForceIndex = (int) ((dr2 - tabulate_start) / tabulate_step);
            if (tabulatedForceIndex >= N_tabulated) {
                fx = 0.0;
                fy = 0.0;
            } else {
                fx = tabulated_f_per_r[tabulatedForceIndex] * dx; ///f*(dx/dr)
                fy = tabulated_f_per_r[tabulatedForceIndex] * dy;
            }

            particles[i].fx += fx;
            particles[i].fy += fy;

            particles[j].fx -= fx;
            particles[j].fy -= fy;

        }
    }
}

void calculatePairwiseForcesWithoutTabulatedForce() {
    int i, j;
    double dx, dy, dr, dr2;
    double f, fx, fy;

    for (i = 0; i < N_particles - 1; i++) {
        for (j = i + 1; j < N_particles; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;


            ////ha igy szamoljuk az erot, az koltseges ezert inkabb a tabulated force listabol vesszuk ki
            dr = sqrt(dr2);

            (dr < 0.2) ? (f = 100.0, printf("Warning!!!dr%f particles %d(%lf %lf) %d(%lf %lf)\n", dr, i,
                                            particles[i].coord.x, particles[i].coord.y, j, particles[j].coord.x,
                                            particles[j].coord.y)) : (f = 1 / dr2 * exp(-0.25 * dr));

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

void calculatePairwiseForcesWithVerlet() {
    int i, j;
    double dx, dy, dr2;
    double fx, fy;
    int tabulatedForceIndex;

    for (int k = 0; k < N_verlet_actual; k++) {
        i = vlist1[k];
        j = vlist2[k];

        dx = particles[i].coord.x - particles[j].coord.x;
        dy = particles[i].coord.y - particles[j].coord.y;

        keepParticleInSystem(&dx, &dy);

        dr2 = dx * dx + dy * dy;

        if (dr2 < 0.01) tabulatedForceIndex = 0;
        else tabulatedForceIndex = (int) ((dr2 - tabulate_start) / tabulate_step);

        if (tabulatedForceIndex >= N_tabulated) {
            fx = 0.0;
            fy = 0.0;
        } else {
            fx = multiply * tabulated_f_per_r[tabulatedForceIndex] * dx; ///f*(dx/dr)
            fy = multiply * tabulated_f_per_r[tabulatedForceIndex] * dy;
        }

        particles[i].fx += fx;
        particles[i].fy += fy;

        particles[j].fx -= fx;
        particles[j].fy -= fy;

    }
}

void calculateClusters() {
    for (int i = 0; i < N_vertex; i++) {
        vertex[i].clusterId = -1; //meg nem tartozik egyik clusterhez sem
    }
    N_clusters = 0;

    int currentVertex = 0;
    int nrVertex = 0; //adott clusterhez tartoz reszecskek szama
    int currentType = 0;
    int nrGoodToCheck = 0; //hany vertexszomszedot jo meg leellenorizni
    int iInside = 0;
    int neighborToCheck = 0;
    int *indicesToCheck;
    indicesToCheck = (int *) malloc(N_vertex * sizeof(int));
    if (indicesToCheck == NULL) {
        perror("Allocation problem with indicesToCheck list, exiting");
        exit(EXIT_FAILURE);
    }
    int currentCluster = 0;

    while (currentVertex < N_vertex) {
        nrVertex = 0;
        nrGoodToCheck = 0;
        if ((vertex[currentVertex].clusterId == -1) &&
            (vertex[currentVertex].type == 5 || vertex[currentVertex].type == 6)) {
            N_clusters++;
            if (N_clusters != currentCluster) {
                currentCluster = N_clusters;
            }
            nrVertex++;
            vertex[currentVertex].clusterId = N_clusters - 1; //az illeto reszecskenek beallitjuk a clusterIDjat
            currentType = vertex[currentVertex].type;
            indicesToCheck[nrGoodToCheck] = currentVertex;
            nrGoodToCheck++;
            iInside = 0;

            while (iInside < nrGoodToCheck) {
                ///neighbor up
                //ha a legfelso vertexrol van szo, meg kell nezni a doboz masik felen levot
                if (indicesToCheck[iInside] % pinningNY == (pinningNY - 1))
                    neighborToCheck = indicesToCheck[iInside] - (pinningNY - 1);
                else neighborToCheck = indicesToCheck[iInside] + 1;

                //ellenorizzuk, hogy az a szomszed ahhoz a clusterhez tartozik-e
                if (vertex[neighborToCheck].type == currentType && vertex[neighborToCheck].clusterId == -1) {
                    nrVertex++;
                    vertex[neighborToCheck].clusterId = N_clusters - 1;
                    indicesToCheck[nrGoodToCheck] = neighborToCheck;
                    nrGoodToCheck++;
                }

                ///neighbor to the right
                if ((indicesToCheck[iInside] / pinningNX) == (pinningNX - 1))
                    neighborToCheck = indicesToCheck[iInside] - (pinningNX - 1) * pinningNY;
                else neighborToCheck = indicesToCheck[iInside] + pinningNY;

                if (vertex[neighborToCheck].type == currentType && vertex[neighborToCheck].clusterId == -1) {
                    nrVertex++;
                    vertex[neighborToCheck].clusterId = N_clusters - 1;
                    indicesToCheck[nrGoodToCheck] = neighborToCheck;
                    nrGoodToCheck++;
                }

                ///neighbor to the left
                if ((indicesToCheck[iInside]) / pinningNX == 0)
                    neighborToCheck = indicesToCheck[iInside] + (pinningNX - 1) * pinningNY;
                else neighborToCheck = indicesToCheck[iInside] - pinningNY;

                if (vertex[neighborToCheck].type == currentType && vertex[neighborToCheck].clusterId == -1) {
                    nrVertex++;
                    vertex[neighborToCheck].clusterId = N_clusters - 1;
                    indicesToCheck[nrGoodToCheck] = neighborToCheck;
                    nrGoodToCheck++;
                }

                ///neighbor down
                if (indicesToCheck[iInside] % pinningNY == 0)
                    neighborToCheck = indicesToCheck[iInside] + (pinningNY - 1);
                else neighborToCheck = indicesToCheck[iInside] - 1;

                if (vertex[neighborToCheck].type == currentType && vertex[neighborToCheck].clusterId == -1) {
                    nrVertex++;
                    vertex[neighborToCheck].clusterId = N_clusters - 1;
                    indicesToCheck[nrGoodToCheck] = neighborToCheck;
                    nrGoodToCheck++;
                }

                iInside++;

            }

        }

        if (N_clusters >= currentCluster && nrVertex > 0) {
            clusters = (int *) realloc(clusters, N_clusters * sizeof(int));
            if (clusters == NULL) {
                perror("Allocation problem with Cluster list, exiting");
                exit(EXIT_FAILURE);
            }
            clusters[N_clusters - 1] = nrVertex;
        }
        currentVertex++;
    }

    free(indicesToCheck);
}

void calculateMultiply(int currentTime, int totalTime) {
    multiply = (double) currentTime / (double) totalTime * maxMultiply;
}

void buildVerletList() {
    N_verlet_actual = 0; //szomszedos reszecskek szama,verlet lista hossz
    double dx, dy, dr, dr2;


    for (int i = 0; i < N_particles; i++) {
        for (int j = i + 1; j < N_particles; j++) {

            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;

            keepParticleInSystem(&dx, &dy);

            dr2 = dx * dx + dy * dy;
            dr = sqrt(dr2);

            if (dr <= rv) {

                N_verlet_actual++;
                if (N_verlet_actual >= N_verlet_list) {

                    N_verlet_list += 10;
                    vlist1 = (int *) realloc(vlist1, N_verlet_list * sizeof(int));
                    vlist2 = (int *) realloc(vlist2, N_verlet_list * sizeof(int));
                    if ((vlist1 == NULL) || (vlist2 == NULL)) {
                        printf("Problem reallocating \n");
                        exit(-1);
                    }

                }

                vlist1[N_verlet_actual - 1] = i;
                vlist2[N_verlet_actual - 1] = j;

            }
        }

    }

    flag_to_rebuild_verlet = 0;
}

void moveParticles() {
    double dx, dy;

    for (int i = 0; i < N_particles; i++) {
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

    intholder = N_particles + N_vertex;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    intholder = t;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    for (i = 0; i < N_particles; i++) {
        intholder = particles[i].color;
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

    for (i = 0; i < N_vertex; i++) {
        intholder = vertex[i].color;
        fwrite(&intholder, sizeof(int), 1, moviefile);
        intholder = vertex[i].id;//ID
        fwrite(&intholder, sizeof(int), 1, moviefile);
        floatholder = (float) vertex[i].coord.x;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = (float) vertex[i].coord.y;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, moviefile);
    }

}

void write_gfile() {
    FILE *f;
    f = fopen("gfile", "wt");
    if (f == NULL) {
        printf("Gfile not created= %d\n", errno);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "set xrange %d %d\n", -1, (int) sX + 1);
    fprintf(f, "set yrange %d %d\n", -1, (int) sY + 1);
    fprintf(f, "loadcolormap colors.txt\n");
    fprintf(f, "loadcontour contour.txt\n");
    fprintf(f, "cmovie\n");
    fclose(f);
    fflush(stdout);
}

void write_contour_file() {
    FILE *f;
    f = fopen("contour.txt", "wt");
    if (f == NULL) {
        printf("Contour file not created, errno = %d\n", errno);
        exit(EXIT_FAILURE);
    }
    fprintf(f, "%d\n", N_pinning);
    for (int i = 0; i < N_pinning; i++) {
        fprintf(f, "%e\n", pinnings[i].coord.x);
        fprintf(f, "%e\n", pinnings[i].coord.y);
        fprintf(f, "%e\n", pinnings[i].r);
        fprintf(f, "%e\n", pinnings[i].r);
        fprintf(f, "%e\n", pinnings[i].r);
    }

    fclose(f);
    fflush(stdout);

}

void write_statistics() {
    double avg_vx = 0.0;

    for (int i = 0; i < N_particles; i++) {
        if (particles[i].color == 0) avg_vx += particles[i].fx;
        if (particles[i].color == 1) avg_vx -= particles[i].fx;
    }

    avg_vx = avg_vx / (double) N_particles;

    fprintf(statistics_file, "%d %f\n", t, avg_vx);
}

void calculateVertexStatistics() {
    int maxcluster = 0;
    int mincluster = N_vertex;
    double average = 0.0;

    for (int i = 0; i < N_vertex; i++) {
        vertexTypes[vertex[i].type]++;
    }

    avgNCluster += N_clusters;
    for (int i = 0; i < N_clusters; i++) {
        average += clusters[i];
        if (clusters[i] > maxcluster) {
            maxcluster = clusters[i];
        }
        if (mincluster > clusters[i]) {
            mincluster = clusters[i];
        }
    }

    avgAvgCluster += average / (double) N_clusters;
    avgBigCluster += maxcluster;
    avgSmallCluster += mincluster;

}

void write_cluster_stat() {
    ///time, cluster_type, cluster_type_avg, nrofcluster,bigcluster,smallcluster,avgCluster
    fprintf(statistics_file, "%d ", t);
    for (int i = 0; i < sumVertexTypes; i++)
        fprintf(statistics_file, "%.2lf  ", (double) (vertexTypes[i] / (double) N_average));


    fprintf(statistics_file, " %.2lf", (double) (avgNCluster / (double) N_average));
    fprintf(statistics_file, " %.2lf", avgBigCluster / (double) N_average);
    fprintf(statistics_file, " %.2lf", avgSmallCluster / (double) N_average);
    fprintf(statistics_file, " %.2lf\n", avgAvgCluster / (double) N_average);
}

void write_cluster_statistics() {
    int maxcluster = 0, maxclusterindex = 0;
    int mincluster = N_vertex, minclusterindex = 0;
    double average = 0.0;
    fprintf(statistics_file, "Number of clusters = %d \n", N_clusters);
    for (int i = 0; i < N_clusters; i++) {
        printf("%d cluster with %d vertex\n", i, clusters[i]);
        fprintf(statistics_file, "%d cluster with %d vertex\n", i, clusters[i]);
        average += clusters[i];
        if (clusters[i] > maxcluster) {
            maxcluster = clusters[i];
            maxclusterindex = i;
        }
        if (mincluster > clusters[i]) {
            mincluster = clusters[i];
            minclusterindex = i;
        }
    }
    if (N_clusters >= 2) {
        fprintf(statistics_file,
                "The biggest cluster %d with %d vertex\n minim cluster %d with %d vertex\n average cluster %lf\n",
                maxclusterindex, maxcluster, minclusterindex, mincluster, average / (double) N_vertex);
    } else {
        fprintf(statistics_file, "The biggest cluster %d with %d vertex\n average cluster size %lf\n", maxclusterindex,
                maxcluster,
                average / (double) N_clusters);
    }

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
    fprintf(statistics_file, "Nr.Partic %d, System size: %.2f,\nProgram ran: %.2f seconds, %.2f minutes.\n",
            N_particles, sX,
            timeDiff, timeDiff / 60);
}

void sortParticles() {
    Particle part;
    for (int i = 0; i < N_particles - 1; i++) {
        for (int j = i + 1; j < N_particles; j++) {
            if (particles[i].mortonId > particles[j].mortonId) {
                part = particles[j];
                particles[j] = particles[i];
                particles[i] = part;
            }
        }
    }
}

//return 1 if ok, file name(from commandline) check
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

void allInit() {

    ///init Pinningnek:Nx,Ny,distX,distY,lx2,ly2,r,fmax,middleHeight
    initPinning(20, 20, 2.0, 2.0, 0.6, 0.2, 0.2, 2.0, 0.15);

    ///General init: nr Particles, nr Verlet list, nr tabulate force, nr Vertex
    ///init Particle: N_particles, sX+sY, dt, r, rv, temperature, multiply, maxmultiply
    init(N_pinning, 40.0, 0.002, 4.0, 6.0, 3.0, 0.0, 1.0);

    initPinningSites();
    write_contour_file();

    //write_gfile();

    // initParticles();
    //initSquareParticles();
    initSquareParticlesRandom();
    //sortParticles();

    initSquareVertex();
}

void simulation(int time, int statisticTime) {
    for (t = 0; t < time; t++) {
        calculatePairwiseForcesWithVerlet();
        calculatePairwiseForces();
//        calculateExternalForces();
        calculateThermalForces();
        calculatePinningForces();

        moveParticles();
        if (flag_to_rebuild_verlet) buildVerletList();

        ///itt allitani a moovie file-nak
        if (t % 300 == 0) {
            calculateVertexTypes();
            calculateClusters();
            calculateVertexStatistics();
            write_cmovie();
        }
        if (t % statisticTime == 0) {
            write_cluster_stat();
            initVertexStat();
            calculateMultiply(t, time);
        }
        if (t % 500 == 0) {
            printf("time = %d\n", t);
            fflush(stdout);
        }
    }
}

int main(int argc, char *argv[]) {
    char filename[100] = {0};
    char statfile[100] = {0};
    char *be = (char *) malloc(100);
    uint64_t seed = 0;
    int lepes = 0;

    ///result, statistic,
    if (argc > 1) {
        strcpy(be, argv[1]);
        strncpy(filename, properFilename(be) ? strcat(be, ".mvi") : "result.mvi",
                sizeof(filename));
        strcpy(be, argv[2]);
        strncpy(statfile, properFilename(be) ? strcat(be, ".txt") : "statistics.txt",
                sizeof(statfile));
        seed = (uint64_t) argv[3];
        lepes = atoi(argv[4]);
    } else {
        strncpy(filename, "result.mvi", sizeof(filename));
        strncpy(statfile, "statistics.txt", sizeof(statfile));
    }


    setseed(seed);
    allInit();
    moviefile = fopen(filename, "wb");
    statistics_file = fopen(statfile, "wt");
    if (!moviefile || !statistics_file) {
        perror("Failed to open file.\n");
        return 1;
    }


    buildVerletList();
    calculateVertexTypes();

    start();
    
    if (!lepes) lepes = 100000;

    N_average = 100;
    simulation(lepes, N_average);
    calculateClusters();
    fclose(moviefile);
    end();
    fclose(statistics_file);
    printf("The result(for plot) can be found in file: %s\n", filename);
    free(particles);
    free(pinnings);

    return 0;
}