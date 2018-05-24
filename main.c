#include "data.h"

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
        //printf("kor %d %.2f %.2f %.2f %.8f\n", i, x, x2, f, tabulated_f_per_r[i]);
    }
}

//alapbeallitasok
void init(int nrParticles, double systemSize, int nrPinnings, double timeStep, double cutOff, double verletCutOff) {
    N = nrParticles;
    particles = (Particle *) malloc(N * sizeof(Particle));
    if (!particles) {
        perror("Allocation problem! Particles");
        exit(EXIT_FAILURE);
    }
    sX = sY = systemSize;
    sX2 = sY2 = systemSize / 2;
    dt = timeStep;

    //Pinnings
    N_pinning = nrPinnings;
    pinnings = (Pinning *) malloc(N_pinning * sizeof(Pinning));
    if (!pinnings) {
        perror("Allocation problem! PinningSites");
        exit(EXIT_FAILURE);
    }

    //Cutoff distance
    r0 = cutOff;
    rv = verletCutOff;
    rvminr02 = (rv - r0) * (rv - r0);

    //Verlet list
    ////N_verlet= (pi*rv^2*N^2)/(2*sX*sY)
    N_verlet_list = ((int) 3.14 * (int) rv * (int) rv * N * N) / (2 * (int) sX * (int) sY);
    vlist1 = (int *) malloc(N_verlet_list * sizeof(int));
    vlist2 = (int *) malloc(N_verlet_list * sizeof(int));
    if (!vlist1 || !vlist2) {
        perror("Allocation problem! Verlet list.");
        exit(EXIT_FAILURE);
    }
    flag_to_rebuild_verlet = 0;

    //Tabulated force
    N_tabulated = 50000;
    tabulated_f_per_r = (double *) malloc(N_tabulated * sizeof(double));
    if (!tabulated_f_per_r) {
        perror("Allocation problem! Tabulated force(f/r) list.");
        exit(EXIT_FAILURE);
    }
    tabulateForces();
    printf("Init done!\n");

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

    while (i < N) {
        particles[i].id = i;
        particles[i].color = ((rand() / (RAND_MAX + 1.0)) > 0.5); //color
        particles[i].q = particles[i].color ? -1.0 : 1.0;
        stowParticle(&position, i);  //particle coordinate
        particles[i].coord.x = position.x;
        particles[i].coord.y = position.y;
        particles[i].fx = 0.0; //force
        particles[i].fy = 0.0;
        particles[i].drx = 0.0;
        particles[i].dry = 0.0;
        i++;
    }
    printf("Init particles done!\n");
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
    int i = 0;
    Coordinate position;

    while (i < N_pinning) {

        pinnings[i].id = i;
        stowPinning(&position, i);  //pinning coordinate
        pinnings[i].coord.x = position.x;
        pinnings[i].coord.y = position.y;
        pinnings[i].r = 1.0; //radius
        pinnings[i].f_max = 2.0;
        i++;
    }
    printf("Init pinnings done!\n");
}

void calculatePinningForces() {
    double dx, dy, dr2, f, fx, fy;

    for (int i = 0; i < N_pinning; i++)
        for (int j = 0; j < N_pinning; j++) {
            dx = particles[i].coord.x - pinnings[j].coord.x;
            dy = particles[i].coord.y - pinnings[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;
            if (dr2 < (pinnings[j].r * pinnings[j].r)) {
                f = 1.0 / pinnings[j].r * pinnings[j].f_max;
                fx = -f * dx;
                fy = -f * dy;

                particles[i].fx += fx;
                particles[i].fy += fy;
            }
        }
}

void calculateExternalForces() {
    int i;
    for (i = 0; i < N; i++) {
        particles[i].fx += (double)t/10000* 5.0;
    }
}

void calculatePairwiseForces() {
    int i, j;
    double dx, dy, dr, dr2;
    double f, fx, fy;
    int tabulatedForceIndex;

    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;

            tabulatedForceIndex = (int) floor((dr2 - tabulate_start) / tabulate_step);
            if (tabulatedForceIndex >= N_tabulated) {
                fx = 0.0;
                fy = 0.0;
            } else {
                fx = tabulated_f_per_r[tabulatedForceIndex] * dx; ///f*(dx/dr)
                fy = tabulated_f_per_r[tabulatedForceIndex] * dy;
            }

            ////ha igy szamoljuk az erot, az koltseges ezert inkabb a tabulated force listabol vesszuk ki
//            dr = sqrt(dr2);
//
//            (dr < 0.2) ? (f = 100.0, printf("Warning!!!dr%f particles %d(%lf %lf) %d(%lf %lf)\n", dr, i,
//                                            particles[i].coord.x, particles[i].coord.y, j, particles[j].coord.x,
//                                            particles[j].coord.y)) : (f = 1 / dr2 * exp(-0.25 * dr));
//
//            //project it to the axes get the fx, fy components
//            fx = f * dx / dr;
//            fy = f * dy / dr;

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
        tabulatedForceIndex = (int) floor((dr2 - tabulate_start) / tabulate_step);
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

void buildVerletList() {
    //printf("Build verlet %d idopillanatban\n",t);
    N_verlet_actual = 0; //szomszedos reszecskek szama,verlet lista hossz
    double dx, dy, dr, dr2;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            dx = particles[i].coord.x - particles[j].coord.x;
            dy = particles[i].coord.y - particles[j].coord.y;
            keepParticleInSystem(&dx, &dy);
            dr2 = dx * dx + dy * dy;
            dr = sqrt(dr2);
            if (dr <= rv) {
                N_verlet_actual++;
                if (N_verlet_actual >= N_verlet_list) {
                    N_verlet_list++;
                    vlist1 = (int *) realloc(vlist1, N_verlet_list * sizeof(int));
                    vlist2 = (int *) realloc(vlist2, N_verlet_list * sizeof(int));
                }
                vlist1[N_verlet_actual - 1] = i;
                vlist2[N_verlet_actual - 1] = j;

            }
        }
        particles[i].drx = 0.0;
        particles[i].dry = 0.0;
    }

    for (int i; i < N; i++) {
        if (particles[i].q == 1.0) particles[i].color = 0;
        else particles[i].color = 1;
    }

    int i, j;
    for (int k = 0; k < N_verlet_actual; k++) {
        i = vlist1[k];
        j = vlist2[k];
        if (i == 30) particles[j].color = 4;
        if (j == 30) particles[i].color = 5;
    }
    flag_to_rebuild_verlet = 0;
}

void moveParticles() {
    double dx, dy;

    for (int i = 0; i < N; i++) {
        dx = particles[i].fx * dt;
        dy = particles[i].fy * dt;
        particles[i].coord.x += dx;
        particles[i].coord.y += dy;
        particles[i].drx += dx;
        particles[i].dry += dy;
        //out of box
        if (particles[i].coord.x < 0) particles[i].coord.x += sX;
        if (particles[i].coord.y < 0) particles[i].coord.y += sY;
        if (particles[i].coord.x > sX) particles[i].coord.x -= sX;
        if (particles[i].coord.y > sY) particles[i].coord.y -= sY;

        if ((particles[i].drx * particles[i].drx + particles[i].dry * particles[i].dry) >= rvminr02)
            flag_to_rebuild_verlet = 1;
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
        intholder = particles[i].color + 4;
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

void write_contour_file() {
    int i;
    FILE *f;

    f = fopen("contour.txt", "wt");

    fprintf(f, "%d\n", N_pinning);


    for (i = 0; i < N_pinning; i++) {
        fprintf(f, "%e\n", pinnings[i].coord.x);
        fprintf(f, "%e\n", pinnings[i].coord.y);
        fprintf(f, "%e\n", pinnings[i].r);
        fprintf(f, "%e\n", pinnings[i].r);
        fprintf(f, "%e\n", pinnings[i].r);
    }
    fclose(f);
}

void write_statistics() {
    double avg_vx = 0.0;
    for (int i = 0; i < N; i++) {
        if (particles[i].color == 0) avg_vx += particles[i].fx;
        if (particles[i].color == 1) avg_vx -= particles[i].fx;
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
    fprintf(statistics_file, "Nr.Partic %d, System size: %.2f,\nProgram ran: %.2f seconds, %.2f minutes.\n", N, sX,
            timeDiff, timeDiff / 60);
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
    ///inicializalas: N, sX+sY, N_pinning, dt, r, rv
    init(400, 40.0, 200, 0.002, 4.0, 6.0);

    initParticles();
    initPinningSites();
    write_contour_file();

    buildVerletList();
    //printf("%d aktualis %d\n", N_verlet_list, N_verlet_actual);

    moviefile = fopen(filename, "wb");
    statistics_file = fopen("statistics.txt", "wt");
    if (!moviefile || !statistics_file) {
        fprintf(stderr, "Failed to open file.\n");
        return 1;
    }

    for (t = 0; t < 10000; t++) {
        calculatePairwiseForcesWithVerlet();
//        calculatePairwiseForces();

        calculateExternalForces();
        calculatePinningForces();
        write_statistics();
        moveParticles();
        if (flag_to_rebuild_verlet) buildVerletList();

        if (t % 100 == 0) {
            write_cmovie();
        }
        if (t % 500 == 0) {
            printf("time = %d\n", t);
            fflush(stdout);
        }
    }

    fclose(moviefile);
    end();
    fclose(statistics_file);
    printf("The result(for plot) can be found in file: %s\n", filename);
    free(particles);
    free(pinnings);
    return 0;
}
