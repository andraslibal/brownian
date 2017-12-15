
/* Simulation Course
 
 Assignment 1
 Molecular Dynamics Simulation / Brownian Dynamics Simulation
 
 Simulate the motion of the particles using brownian dynamics
 Pedestrian crossing problem
 
 Assignment 2
 Melting of a crystal
 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

struct particle_struct
{
    double drx_so_far,dry_so_far;
    double x,y;                         //x,y coordinate of the particles
    double fx,fy;                       //fx,fy forces acting on the particle
    int color;                          //this is to distinguish the particles
    int ID;                             //ID of a particle
} *particles;

struct pinning_struct
{
    double x,y; //center of the pinning site
    double r;   //radius - this is how far it affects particles
    double f_max; //maximum force (at the boundary)
} *pinningsites;

double SX,SY;           //system size x,y direction
double SX2,SY2;         //half of the system size x,y direction
int N;                  //number of particles
int N_pins;             //number of pinningsites
double dt;              //length of a single time step
int t;                  //time - time steps so far
FILE *moviefile;        //file to store the coordinates of the particles
FILE *statistics_file;

//Verlet list variables
int *vlist1=NULL;
int *vlist2=NULL;
int N_vlist;

//this flag will tell me whether
//I need to rebuild the Verlet list
int flag_to_rebuild_Verlet;

/*
What the Verlet list stores:

 vlist1[5] vlist2[5] - these are the i,j numbers of 2 particles that interact
 vlist1[17] vlist2[17] - these are the i,j numbers of a different two particles that interact
 
 this is how I read the i,j from the Verlet list
 i = vlist1[7]
 j = vlist2[7]
 
 particles[vlist1[7]].x
 particles[vlist2[7]].x
 particles[vlist1[7]].y
 particles[vlist2[7]].y
 */

//variables for tabulating the force
double *tabulated_f_per_r;
int N_tabulated;
double tabulalt_start, tabulalt_lepes;

//these are for time keeping purposes
//to count how many seconds the simulation ran
time_t 	time_start;
time_t 	time_end;

void program_timing_begin()
{
    time(&time_start);
}

void program_timing_end()
{
    double time_difference;
    struct tm *timeinfo;
    
    time(&time_end);
    
    time_difference = difftime(time_end,time_start);
    
    timeinfo = localtime(&time_start);
    printf("Program started at: %s",asctime(timeinfo));
    
    timeinfo = localtime(&time_end);
    printf("Program ended at: %s",asctime(timeinfo));
    
    printf("Program running time = %lf seconds\n", time_difference);
    printf("%d %lf\n",N,time_difference);
}

void tabulate_forces()
{
    int i;
    double x_min,x_max;
    double x2,x;
    double f;
    
    x_min = 0.1;
    x_max = 6.0;
    
    double x_min2 = x_min * x_min;
    double x_max2 = x_max * x_max;
    
    N_tabulated = 50000;
    tabulated_f_per_r = (double *) malloc(N_tabulated*sizeof(double));
    for(i=0;i<N_tabulated;i++)
    {
        x2 = i*(x_max2-x_min2)/(N_tabulated-1.0) + x_min2;
        x = sqrt(x2);
        f = 1/x2 * exp(-0.25*x);
        tabulated_f_per_r[i] = f/x;
        //printf("%d %lf %lf %lf %lf\n",i,x,x2,f,tabulated_f_per_r[i]);
    }
    
    tabulalt_start = x_min2;
    tabulalt_lepes = (x_max2-x_min2)/(N_tabulated-1.0);
    printf("Tabulalt start = %lf, lepes = %lf\n",tabulalt_start,tabulalt_lepes);
    
    /*
     0 xmin^2                    f/r
     1 xmin^2 + 1 * lepes        f/r
     2 xmin^2 + 2 * lepes        f/r
     
     lepes = (x_max*x_max-x_min*x_min)/(N_tabulated-1.0)
     
     
     dr^2
     
     index = (int)floor( ( dr^2 - tabulalt_start) / tabulalt_lepes )
     
     
     */
}

void rebuild_verlet_list()
{
    int i,j;
    double dx,dy,dr2;
    
    //printf("rebuilding Verlet\n");fflush(stdout);
    
    N_vlist = 0;
    vlist1 = (int *) realloc(vlist1,N_vlist*sizeof(int));
    vlist2 = (int *) realloc(vlist2,N_vlist*sizeof(int));
    
    for(i=0;i<N;i++)
        for(j=i+1;j<N;j++)
        {
            dx = particles[i].x - particles[j].x;
            dy = particles[i].y - particles[j].y;
            
            //PBC check
            //(maybe the neighbor cell copy is closer)
            if (dx>SX2) dx -=SX;
            if (dx<-SX2) dx +=SX;
            if (dy>SY2) dy -=SY;
            if (dy<-SY2) dy +=SY;
            
            dr2 = dx*dx+dy*dy;
            
            //dr=sqrt(dr2)
            //if dr<6.0
            
            if (dr2<=36.0) //instead of 4*4 I will take 6*6
            {
                N_vlist++;
                vlist1 = (int *) realloc(vlist1,N_vlist*sizeof(int));
                vlist2 = (int *) realloc(vlist2,N_vlist*sizeof(int));
                vlist1[N_vlist-1] = i;
                vlist2[N_vlist-1] = j;
            }
        }
    /*
     for(i=0;i<N;i++)
     printf("%d %d \n",vlist1[i],vlist2[i]);
     */
    
    //once I rebuilt the Verlet list,
    //I can start counting the distances again
    for(i=0;i<N;i++)
    {
        particles[i].drx_so_far = 0.0;
        particles[i].dry_so_far = 0.0;
    }
    flag_to_rebuild_Verlet = 0;
    //printf("Verlet rebuilt at t=%d\n",t);
}

void initialize_particles()
{
    int i,j,overlap;
    double dx,dy,dr,dr2;
    double tempx,tempy;
    double multiplier;
    
    //multiplier = 512.0;
    
    SX = 60.0;//*sqrt(2.0);//2.0*sqrt(multiplier);//2.5*sqrt(multiplier);
    SY = 60.0;//*sqrt(2.0);//2.0*sqrt(multiplier);
    SX2 = SX/2.0;
    SY2 = SY/2.0;
    
    N = 2000;
    
    particles = (struct particle_struct *) malloc(N*sizeof(struct particle_struct));
    
    dt = 0.002;
    /*
    for(i=0;i<23;i++)
        for(j=0;j<23;j++)
            {
            tempx = ((double)i+(j%2)*0.5) * 2.0;
            tempy = (double)j * 2.0;
            
            //printf("%lf %lf\n",tempx,tempy);fflush(stdout);
            
            particles[ii].x = tempx;
            particles[ii].y = tempy;
    */
    for(i=0;i<N;i++)
    {
        particles[i].ID =i;
        
        do
        {
            overlap = 0;
            tempx = SX * rand()/(RAND_MAX+1.0);
            tempy = SY * rand()/(RAND_MAX+1.0);
            
            for(j=0;j<i;j++)
            {
                dx = tempx-particles[j].x;
                dy = tempy-particles[j].y;
                //PBC fold back
                if (dx>=SX2) dx-= SX;
                if (dx<-SX2) dx+= SX;
                if (dy>=SY2) dy-= SY;
                if (dy<-SY2) dy+= SY;
                
                dr2 = dx*dx + dy*dy;
                dr = sqrt(dr2);
                if (dr<0.2) 
                {
                	overlap=1;
                	break;
                }
            }
        }
        while (overlap==1);
        //N_probalkozas++;
        //a better version would check how many attempts were made
        //after too many attempts quit gracefully (and not get stuck)
        
        particles[i].x = tempx;
        particles[i].y = tempy;
        
        //*/
     
        //solve the problem: two particles should never be on top of each other!!!
        // 0.2 safe distance, I know the force at 0.2, it's not that big (around 100.0)
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        
        particles[i].drx_so_far = 0.0;
        particles[i].dry_so_far = 0.0;
        
        particles[i].color = 0;
        
        //if ( rand()/(RAND_MAX+1.0) < 0.5)   particles[i].color = 0;
        //else                                particles[i].color = 1;
        
        
        //printf("%d",particles[i].color);
        //rand()%100 Never ever use this when generating random numbers
        //has very bad properties
    }
    
}

void initialize_pinning_sites()
{
    int i,j,overlap;
    double dx,dy,dr,dr2;
    double tempx,tempy;

    N_pins = 400;
    
    pinningsites = (struct pinning_struct *) malloc(N*sizeof(struct pinning_struct));
    
    for(i=0;i<N_pins;i++)
    {
        
        do
        {
            overlap = 0;
            tempx = SX * rand()/(RAND_MAX+1.0);
            tempy = SY * rand()/(RAND_MAX+1.0);
            
            for(j=0;j<i;j++)
            {
                dx = tempx-pinningsites[j].x;
                dy = tempy-pinningsites[j].y;
                dr2 = dx*dx + dy*dy;
                dr = sqrt(dr2);
                if (dr<2.2) 
                {
                	overlap=1;
                	break;
                }
            }
        }
        while (overlap==1);
        //a better version would check how many attempts were made
        //after too many attempts quit gracefully (and not get stuck)
        
        pinningsites[i].x = tempx;
        pinningsites[i].y = tempy;
        
        pinningsites[i].f_max = 2.0;
        pinningsites[i].r = 1.0;
        
        }
    
}

void write_contour_file()
{
	int i;
	FILE *f;

	f=fopen("contour.txt","wt");

	fprintf(f,"%d\n",N_pins);


	for(i=0;i<N_pins;i++)
		{
		fprintf(f,"%e\n",pinningsites[i].x);
		fprintf(f,"%e\n",pinningsites[i].y);
		fprintf(f,"%e\n",pinningsites[i].r);
		fprintf(f,"%e\n",pinningsites[i].r);
		fprintf(f,"%e\n",pinningsites[i].r);
		}
	fclose(f);
}


void write_particles()
{
    FILE *f;
    int i;
    
    f = fopen("test.txt","wt");
    for(i=0;i<N;i++)
        fprintf(f,"%lf %lf\n",particles[i].x,particles[i].y);
    fclose(f);
}

void calculate_thermal_force()
{
    int i;
    
    for(i=0;i<N;i++)
    {
        //rand() gives an integer 0 ... RAND_MAX
        //rand()/(RAND_MAX+1.0) this is a double between [0,1)
        //this is a well behaving random number
        //gaussian distributed random number generator
        //Numerical Recipes in C - gasdev()
        //particles[i].fx += 3.0 * gasdev();
        
        particles[i].fx += 5.0 * (rand()/(RAND_MAX+1.0)-0.5);
        particles[i].fy += 5.0 * (rand()/(RAND_MAX+1.0)-0.5);
    }
}

void calculate_external_forces()
{
    int i;
    
    for(i=0;i<N;i++)
    {
        //if (particles[i].color==0)    particles[i].fx += 0.5;
        //if (particles[i].color==1)    particles[i].fx -= 0.5;
        
        particles[i].fx += 2.0*(double)t/100000.0;
    }
}

void calculate_pinning_force()
{
	int i,j;
	double dx,dy,dr2,dr,f,fx,fy;

	 //inefficient simple way fo calculating the pinning force
	 //best way would be with a Verlet lookup grid
	 for(i=0;i<N;i++)
		for(j=0;j<N_pins;j++)
	    {
	        dx = particles[i].x - pinningsites[j].x;
	        dy = particles[i].y - pinningsites[j].y;
	        
	        //PBC check
	        //maybe the neighbor cell copy of j is closer
	        if (dx>SX2) dx -=SX;
	        if (dx<-SX2) dx +=SX;
	        if (dy>SY2) dy -=SY;
	        if (dy<-SY2) dy +=SY;
	        
	        dr2 = dx*dx+dy*dy;
	        
	        //we are calculating the forces directly
	        //dr = sqrt(dr2);
	        
	        if (dr2<(pinningsites[j].r*pinningsites[j].r))
	            {
	            
	                /*
	                
	                f = pinningsites[j].f_max * 1.0/pinningsites[j].r;
	                fx = f * dx;
	                fy = f * dy;
	                
	                */
	                
	                f = pinningsites[j].f_max/pinningsites[j].r;
	                fx = -f * dx;
	                fy = -f * dy;
	                
	                particles[i].fx += fx;
	                particles[i].fy += fy;
	                
	            }
	    }
}

void calculate_pairwise_forces()
{
    int i,j;
    double dx,dy;
    double dr2,dr;
    double f,fx,fy;
    
    for(i=0;i<N-1;i++)
        for(j=i+1;j<N;j++)
        {
            dx = particles[i].x - particles[j].x;
            dy = particles[i].y - particles[j].y;
            
            //PBC check
            //maybe the neighbor cell copy of j is closer
            if (dx>SX2) dx -=SX;
            if (dx<-SX2) dx +=SX;
            if (dy>SY2) dy -=SY;
            if (dy<-SY2) dy +=SY;
            
            dr2 = dx*dx+dy*dy;
            
            //we are calculating the forces directly
            dr = sqrt(dr2);
            
            if (dr<0.2) {f = 100.0;printf("WARNING!!!\n");}
            //nice way to do this: give a warning or exit if this happens
            else
                //check if dr>4.0 I can cut off the force
            {
                f = 1/dr2 * exp(-0.25*dr);
            }
            
            //project it to the axes get the fx, fy components
            fx = f*dx/dr;
            fy = f*dy/dr;
            
            
            particles[i].fx += fx;
            particles[i].fy += fy;
            
            particles[j].fx -= fx;
            particles[j].fy -= fy;
        }
}

void calculate_pairwise_forces_with_verlet()
{
    int i,j,ii;
    double dx,dy;
    double dr2,dr;
    double f,fx,fy;
    int tab_index;
    
    for(ii=0;ii<N_vlist;ii++)
    {
        i = vlist1[ii];
        j = vlist2[ii];
        //printf("%d %d\n",i,j);
        dx = particles[i].x - particles[j].x;
        dy = particles[i].y - particles[j].y;
        
        //PBC check
        if (dx>SX2) dx -=SX;
        if (dx<-SX2) dx +=SX;
        if (dy>SY2) dy -=SY;
        if (dy<-SY2) dy +=SY;
        
        dr2 = dx*dx+dy*dy;
        
        //recall the tabulated value of the force
        
        
         tab_index = (int)floor((dr2-tabulalt_start)/(tabulalt_lepes));
         if ((tab_index>=N_tabulated))
         {
		     //printf("tab_index = %d\n",tab_index);
		     //tab_index = N_tabulated-1;
		     //exit(1);
		     fx = 0.0;
		     fy = 0.0;
         }
         else
         {
		     fx = tabulated_f_per_r[tab_index] * dx;  //f/dr is what I recalled
		     fy = tabulated_f_per_r[tab_index] * dy;
         }
        
        
        
        //direct calculation of the force
        /*
        dr = sqrt(dr2); //this is EVIL
        
        if (dr<0.1)
        {
            f = 97.53;
            printf("Warning! Particles %d and %d too close at time %d\n",i,j,t);
        }
        else
            //check if dr>4.0 I can cut off the force
        {
            f = 1/dr2 * exp(-0.25*dr);
        }
        
        fx = f*dx/dr;
        fy = f*dy/dr;
        //printf("%lf %lf hasonlitva %lf %lf\n",fx,fy,f*dx/dr,f*dy/dr);
        */
        
        particles[i].fx += fx;
        particles[i].fy += fy;
        particles[j].fx -= fx;
        particles[j].fy -= fy;
    }
}


void move_particles()
{
    int i;
    double deltax,deltay;
    
    for(i=0;i<N;i++)
    {
        //brownian dynamics
        //the particle is in a highly viscous environment
        deltax = particles[i].fx * dt;
        deltay = particles[i].fy * dt;
        
        particles[i].x += deltax;
        particles[i].y += deltay;
        
        particles[i].drx_so_far += deltax;
        particles[i].dry_so_far += deltay;
        
        
        if ((particles[i].drx_so_far*particles[i].drx_so_far +
             particles[i].dry_so_far*particles[i].dry_so_far) >= 4.0)
            flag_to_rebuild_Verlet = 1;
        
        //PBC check - check if they left the box
        //Box: 0,0 to SX, SY
        if (particles[i].x > SX) particles[i].x -=SX;
        if (particles[i].y > SY) particles[i].y -=SY;
        if (particles[i].x < 0) particles[i].x +=SX;
        if (particles[i].y < 0) particles[i].y +=SY;
        
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
    }
}

//this is for plot(linux plotter)
void write_cmovie()
{
    int i;
    float floatholder;
    int intholder;
    
    intholder = N;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    intholder = t;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    for (i=0;i<N;i++)
    {
        intholder = particles[i].color+2;
        fwrite(&intholder,sizeof(int),1,moviefile);
        intholder = i;//ID
        fwrite(&intholder,sizeof(int),1,moviefile);
        floatholder = (float)particles[i].x;
        fwrite(&floatholder,sizeof(float),1, moviefile);
        floatholder = (float)particles[i].y;
        fwrite(&floatholder,sizeof(float),1,moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder,sizeof(float),1,moviefile);
    }
    
}

void write_statistics()
{
	int i;
	double avg_vx;

	avg_vx = 0.0;
	for (i=0;i<N;i++)
		{
		    avg_vx += particles[i].fx;
		}

	avg_vx = avg_vx/(double)N;

	fprintf(statistics_file,"%d %lf\n",t,avg_vx);

}


int main()
{
    printf("Simulation 1 run\n");
    program_timing_begin();
    
    tabulate_forces();
    
    initialize_particles();
    initialize_pinning_sites();
    write_contour_file();
    rebuild_verlet_list();
    
    moviefile = fopen("results.mvi","w");
    //write_movie_header();
    statistics_file = fopen("stat.txt","wt");
    
    for(t=0;t<100000;t++)
    {
        calculate_pairwise_forces_with_verlet();

        
        //calculate_thermal_force();
        //calculate_pairwise_forces();
        
        
        calculate_external_forces();
        calculate_pinning_force();
        
        //right now I have all the information
        //time to calculate some statistics
        write_statistics();
        
        move_particles();
        
        
        if (flag_to_rebuild_Verlet)
            rebuild_verlet_list();
        
        
        if (t%100==0)
            write_cmovie();
        //write_movie_frame();
        
        if (t%500==0)
        {
            printf("time = %d\n",t);
            fflush(stdout);
        }
    }
    
    fclose(moviefile);
    fclose(statistics_file);
    program_timing_end();
    return 0;
}

