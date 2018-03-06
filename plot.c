/*Special version for showing tmovies.  3D movies*/
/*Revisions log:
 *4.1.05 For network project, adding a new type of network contour.
 *4.1.05 del-plot9.  Adding some new commands: xmagnify and ymagnify, which
         scale the input values, as well as xshift and yshift, for moving
	 the location of periodic boundaries from within the plotting
	 program on command.
 *12.11.03 Changing functionality of the onevortex command so that it
           will draw trails behind only one vortex.
 *12.10.03 Adding triangular pin plotting capability.
 *8.13.03 Delaunay plot works but can't image force chains the way I wanted.
          Thus I move the analysis outside of this program.  For convenience,
	  I introduce a new movie format which replaces "zmovie".  This
	  format is like kmovie except it has 2 ints and 3 floats.  The
	  extra int can (at last!) be used to properly store color
	  information.  I call it "cmovie" format.
 *8.13.03 Adjusting xwidth and ywidth for delaunay plot.
 *8.13.03 Adding command "delaunay" which toggles Delaunay triangulation
          plot (rather than the Voronoi polygon plot).  The triangles
	  are unfilled.
 *8.12.03 Cleaning code and using ANSI standard function call format.
          Adding Voronoi memory leak fix.
 *8.11.03 Attempting to assemble a version of del-plot that contains all
          of the scattered improvements made over time.  Removing the
	  y-offset that was still in del-plot7.  Adding code to handle
	  two (or more) particle species with overloaded tmovie format.
	  Adding special code so that, for a system containing only 3
	  particles, they will be plotted with different sizes.
	  Adding stripe plotting capability for kmovie format.
	  Adding grain plotting capability.  Rather than increasing 
	  size of the old_pos structure to handle the artificially large 
	  ID numbers of granular dimers and trimers (as in del-plotG2),
	  I include the following choice implemented in del-plotG3:
	  Historically, old_pos is indexed using particle ID number.  This
 	  was necessary for vortex avalanche type runs, where particles
 	  were being added to and removed from the data structure frequently.
 	  For the grain runs, particles will (so far) be only added and
 	  not removed, so it is safe to use absolute array position rather
 	  than particle ID number as the index on old_pos.  This way
 	  the huge ID numbers given to some trimer grains will not cause
 	  continuing overflows of old_pos.
	  Finally, adding granular chain plotting capability.
 *5.27.03 Contains color fix for linux boxes.  A few modifications
          relative to del-plot6CT.c.
 *3.1.01 Version for pins of varying radii.
 *3.2.00 Special setting: offsetting lattice by SY/2.
 *6.14.99 Special 2-layer setting
 *3.5.99 Setting to read michigan style 3D contours
 *2.21.99 Adding greater flexibility to voronoi plotting: can now plot
          choice of layer in 3D systems.
 *2.10.99 For deep systems.
 *2.5.99 Adding pinfile pin position plotting capability.  For now, I
         take just the top layer of pinning and plot as open circles
	 using the existing 2D routine.
 *12.1.98 Adding code for viewing a single *vortex* in all layers.
 *9.23.98 Adding zmovie capability
 *9.22.98 Altering 3D layer-radius calculation.  Increasing the radius
          only on every 2nd layer.  Used for systems with 30 layers.
 *5.18.98 Program spends over 70% of its time in my brute force matching
          routine.  Replaced this routine with something much simpler and
	  faster.
 *5.3.98  Attempting to color-code defect sites in voronoi mode.  Adding
          array sidenum[] to hold information on the number of sides.
 *5.2.98  Code will now plot smovies, kmovies, and tmovies.  Next step:
          the triangulations.  Random old stuff: the "clear" command from
	  command window doesn't do anything; only moving the whole
	  window with the mouse works.
 *5.1.98  Attempting to adapt code to animate the att Delaunay triangulation
          data.  The goal is to have a multi-mode program.  Normal mode
	  will be to plot smovies, kmovies, or tmovies, as specified.
	  In the new second mode, the program will generate the Delaunay
	  triangulation of the (one node) data *on the fly* and plot the
	  results.  To work well, this second mode requires high-resolution
	  movies.
 *10.13.97 Testing method to plot antivortices with a reversed
          color scheme.
 *10.12.97 Version to read in parallel tmovies.  Specify number of
          nodes with command "nodenum".  Default is 1.  Can handle
	  up to 16 nodes.
 *10.8.97 Continue command does work just fine.  Fixed clock read-out.
          Help screen working.  Last item: Zoom.
 *10.8.97 One more shot at unbreakable rewind.  Fixed this time
          with a simpler scheme.  Several other things are still
	  broken, including zoom and continue.  Still trying to
	  get everything fixed.  Adding a help screen.
 *10.7.97 Unbreakable rewind needs more work.  Disabling.
 *10.7.97 Setting up unbreakable rewind, using a structure that holds
          the sizes of the frames as they are read in when going
	  forward.  This rewind is limited by file size, so it
	  can be disabled when plotting extremely large files.  In
	  this case the breakable rewind is used instead.
	  This is controlled by the variable deactivate_rewind;
	  initially set to zero, the variable can never be
	  reset to zero once it has been triggered.  Initial
	  frame_num is now 1 instead of 0.
 *10.6.97 Considerable revisions.  Fixed the rewind, although my
          method will break when particles start being removed
	  from the system.  Fixed fast forward.  Also added a few
	  features I always thought would be nice.
 *10.5.97 Continuing to clean and restructure.  Presently working
          on get_string, out of get_command.  current_pen is never
	  used; removing.  Adding default system size.
 *10.4.97 Thoroughly linted, finally.  Note that the X11 libraries
          live in /usr/openwin/include/X11
	  Killing a large number of global variables.
 *10.4.97 Further revised from bin_plot3D.  Having another go at
          reorganizing program into something easier to follow.
	  TO FIX:
            Possible rewinding problem in plot_frame.
	    Check use of number of pins generally.
 *10.4.97 Single layer works fine.  Working with color multiple layer.
          The variable BOX_SIZE is used to set the number of pixels
	  wide each vortex is.  I want to alter this so that vortices
	  on lower layers are drawn larger, so that they can be seen
	  beneath their upper neighbors.  Working.  Note that the 
	  "canonical" value for BOX_SIZE is 3.  I decide on the
	  following change:  Previously, plot_object assumed everybody
	  was the same size, and had Box and Box2 as static variables.
	  Now I'm having them passed to the routine with properly
	  set values.  When the number of layers increases significantly
	  in the future, perhaps cijol can develop a "flag" method
	  whereby layer depth is indicated by a stick emerging from the
	  vortex at a radial angle proportional to depth in sample.
 *10.2.97 Adapted from bin-plot to show three-dimensional movies.  Each
          layer is displayed in a different color.  User can specify
	  single-layer mode, if desired.
 *6.10.96 On line 1135 or so, in the routine draw_contour, you will find
          code for altering the way the pinning site size is determined.
	  At this location, the size can be either fixed or scaled according
	  to the actual pinning radius.  If the latter, is chosen, the
	  program will read the pinning radius in from a file into the
	  variable PinSize.
 *4.23.96 load_special_contour completed.
 *4.23.96 All code related to drawing old-style bitmap contours removed.
          New code to draw background based on special_contour files
          in development.
 *4.23.96 In keeping with the radical changes, I rename this program.
          Formerly x-plot_smovieII.c, it is now bin-plot.c
 *4.23.96 Lots of organizational improvements.  plot_pic, plot_file combined.
          cul-de-sac call to check_events eliminated.
 *4.23.96 Working on routine draw_contour.  Basic interface
          set up.  Now I add load_special_contour to read in info.
 *4.23.96 Status:  toggle_traj works for current driven movies, but not 
          quite for aval configurations.  I leave it for the moment and
          address next concern:  drawing pinning sites.
 *4.23.96 Routines not currently in use removed.
 *4.23.96 Code linted and functioning for the most part.  Slight problem
          with toggle_traj will now be addressed.  Inspect use of
          variable traj_on.
 *4.22.96 Code reorganized.  Special lint setting.
 *4.22.96 check_input eliminated.
 *4.22.96 Further tinkering by CIJOL.  Goal:  get trajectories working.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <termios.h>
#include <math.h>
#include "X11/Xlib.h"
#include "X11/Xutil.h"
#include "X11/Xos.h"
#include "X11/cursorfont.h"

#define MAX(x,y) (x>y) ? x : y
#define MIN(x,y) (x<y) ? x : y
#define ABS(x)   (x<0) ? -x: x
#define PI 3.14159265

#define MAX_OBJECT 50000
/* For Voronoi */
#define MAXNUM 20000

#define SMOVIE 0
#define KMOVIE 1
#define TMOVIE 2
#define CMOVIE 3

/* CIJOL For overloaded tmovie format, LNUM also sets the change in size */
/* for the different particles, so set it to 4 for 4 different types of */
/* particles, for example. */
#define LNUM 8 /* Max. number of layers */

#define BOX_SIZE 3

/* The data necessary for including Jared's "(VERY) lame" icon */

#define icon_width 40
#define icon_height 40
static char icon_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x80, 0x3f, 0x00, 0xf0, 0x01, 0xc0, 0x7f, 0x00, 0xfc,
   0x07, 0xe0, 0xff, 0x00, 0xfe, 0x0f, 0xf0, 0xff, 0x01, 0xfe, 0x0f, 0xf0,
   0xff, 0x01, 0xff, 0x1f, 0xf0, 0xff, 0x01, 0xff, 0x1f, 0xf0, 0xff, 0x01,
   0xff, 0x1f, 0xf0, 0xff, 0x01, 0xff, 0x1f, 0xf0, 0xff, 0x01, 0xff, 0x1f,
   0xf0, 0xff, 0x01, 0xfe, 0x0f, 0xe0, 0xff, 0x00, 0xfe, 0x0f, 0xc0, 0x7f,
   0x00, 0xfc, 0x07, 0x80, 0x3f, 0x00, 0xf0, 0x01, 0x00, 0x0c, 0x00, 0x60,
   0x00, 0x00, 0x18, 0x00, 0x30, 0x00, 0x00, 0x30, 0x00, 0x18, 0x00, 0xc0,
   0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff,
   0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x03,
   0xc0, 0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x7f, 0xc0, 0xff,
   0xff, 0xff, 0xff, 0xc0, 0xff, 0xff, 0xff, 0xff, 0xc0, 0xff, 0xff, 0xff,
   0xff, 0xc0, 0xff, 0xff, 0xff, 0x7f, 0xc0, 0xff, 0xff, 0xff, 0x03, 0xc0,
   0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x03, 0xc0, 0xff, 0xff,
   0xff, 0x03, 0xc0, 0xff, 0xff, 0xff, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

/* Display and screen_num are used as arguments to nearly every Xlib routine, 
 * so it simplifies routine calls to declare them global.  If there were 
 * additional source files, these variables would be declared extern in
 * them. */

Display *display;
Window win;
int screen_num;
Screen *screen_ptr;
Pixmap pixmaps[256];
unsigned long c_map[256];
/* pixel values */
unsigned long foreground_pixel, background_pixel, border_pixel;
Pixmap pixmap,back_pixmap,traj_pixmap;
GC pix_gc[256];
GC pixmap_gc,back_gc,traj_gc;
GC gc,xor_gc;
XFontStruct *font_info;
Cursor cursor;
Pixmap icon_pixmap;

#define EVENT_MASK (ExposureMask | ButtonPressMask | ButtonReleaseMask | \
		    Button1MotionMask | StructureNotifyMask | KeyPressMask)

#define BORDER  10
#define BORDERY 10

unsigned int win_width,win_height;

float x_scale,y_scale;
float x_offset,y_offset;

struct smdata {
  int   p_num;
  float x;
  float y;
} smdata;

struct kmdata {
  int   p_num;
  float x;
  float y;
  float cum_disp;
  float cum_x;
} kmdata;

struct cmdata {
  int   color;
  int   p_num;
  float x;
  float y;
  float cum_disp;
} cmdata;

struct tmdata {
  int   layr;
  int   p_num;
  float x;
  float y;
} tmdata;

struct smdata smframes[MAX_OBJECT];
struct kmdata kmframes[MAX_OBJECT];
struct cmdata cmframes[MAX_OBJECT];
struct tmdata tmframes[LNUM][MAX_OBJECT];

int *rewind_dat[16];  /* Used for unbreakable rewind */
int disable_rewind;

struct pinning_site {
  float x;
  float y;
  float Ap;
  float depth;
  float radiusx;
  float radiusy;
  int id;
  struct pinning_site *next;
  /* For the triangles */
  float px[3],py[3];
};

struct pinning_site *pin_sites = NULL;

FILE *picfile[16];

/*============================== main =============================*/
main(int argc,char *argv[])
{
  int initialize(int samples,float *x_pts,float *y_pts,int *c_pts,int *p_pts,
		 char *display_name,int *do_color,int monochrome,int argc,
		 char *argv[],char execname[],int *do_file,float *uxmin,
		 float *uxmax,float *uymin,float *uymax,char last_file_name[]);
  void process_command(char comm_line[],char pl_file[],int *do_plot,
		       int *stop_now,int *ipos,int *continue_plotting,
		       int *num_pins,int *one_layer,int *one_vortex,
		       int *vortex_num,int *layer_num,int argc,char **argv,
		       int *mode256,int *max_color,int *samples,float *x_pts,
		       float *y_pts,int *c_pts,int *p_pts,int *do_color,
		       int *monochrome,char execname[],float *uxmin,
		       float *uxmax,float *uymin,float *uymax,
		       unsigned int *delay_counter,int *frame_num,
		       int *do_S_contour,int *do_T_contour,
		       float *PinSize,int *do_file,
		       int *do_clear,int *traj_on,float *orig_xmin,
		       float *orig_xmax,float *orig_ymin,float *orig_ymax,
		       int *num_nodes,int *do_voronoi,int *voronoi_layer,
		       int *movie_type,int *ntypes,int *do_stripe,
		       float *xmagnify,float *ymagnify,int *do_xshift,
		       int *do_yshift,float *syssizex,float *syssizey,
		       float *xshift,float *yshift,
		       int *do_net_contour,float *centerx,float *centery);
  int plot_frame(int *do_rewind,int *do_fast_forward,int *one_layer,
		 int layer_num,int one_vortex,int vortex_num,int do_color,
		 unsigned int delay_counter,int *frame_num,int do_S_contour,
		 int do_T_contour,
		 int *do_clear,int old_pos[][MAX_OBJECT][2],int traj_on,
		 int num_nodes,int movie_type,int do_voronoi,int voronoi_layer,
		 int ntypes,int do_stripe,int,float xmagnify,float ymagnify,
		 int do_xshift,int do_yshift,float syssizex,float syssizey,
		 float xshift,float yshift,int do_net_contour);
  int get_arg(char *line,int *ipos,char *arg);
  void setwindow(float xmin,float ymin,float xmax,float ymax,int num_pins,
		 int argc,char **argv,int max_color,int *do_color,
		 int monochrome,int do_S_contour,int do_T_contour,
		 float PinSize,int *do_clear,
		 int traj_on,float *orig_xmin,float *orig_xmax,
		 float *orig_ymin,float *orig_ymax,int do_net_contour,
		 float centerx,float centery,int do_xshift,int do_yshift,
		 float xshift,float yshift,float xmagnify,float ymagnify,
		 float syssizex,float syssizey,int redrawcontour);
  void check_events(int *key_press,int *do_rewind,int *do_fast_forward,
		    int num_pins,int argc,char **argv,int max_color,
		    int *do_color,int monochrome,float uxmin,float uxmax,
		    float uymin,float uymax,int frame_num,int do_S_contour,
		    int do_T_contour,
		    float PinSize,int *do_clear,int traj_on,int *paused,
		    float *orig_xmin,float *orig_xmax,float *orig_ymin,
		    float *orig_ymax,int do_net_contour,float centerx,
		    float centery,int do_xshift,int do_yshift,float xshift,
		    float yshift,float xmagnify,float ymagnify,
		    float syssizex,float syssizey);
  void pause_events(int num_pins,int argc,char **argv,int max_color,
		    int *do_color,int monochrome,float *uxmin,float *uxmax,
		    float *uymin,float *uymax,int do_S_contour,
		    int do_T_contour,float PinSize,
		    int *do_clear,int traj_on,int *paused,float *orig_xmin,
		    float *orig_xmax,float *orig_ymin,float *orig_ymax,
		    int do_net_contour,float centerx,float centery,
		    int do_xshift,int do_yshift,float xshift,float yshift,
		    float xmagnify,float ymagnify,float syssizex,
		    float syssizey);

  char plot_name[200];
  char execname[200];
  int do_plot;            /* Boolean: to plot or not to plot */
  int stop_now;           /* Boolean: if true, stop plotting and exit */
  int error;              /* Boolean: triggered by errors.*/
  int ipos;
  int key_press=0;        /* Boolean: user pressed key while plot running */
  int do_file;             /* indicates if executing commands from file */
  int do_rewind = 0;       /* indicates if we're rewinding */
  int do_fast_forward = 0; /* indicates if we're fast forwarding */
  int continue_plotting = 0;
  int num_pins = 0;        /* Number of pins in contour file, if used.*/
  int do_S_contour = 0;    /* indicates if we're doing a special contour*/
  int do_T_contour = 0;   /* indicates whether to draw triangular pins. */
  int do_net_contour = 0; /* indicates whether to draw network contour */
  int do_clear = 0;
  int traj_on = 0;
  int paused = 0;
  int one_layer = 0;       /* Plot only a single layer */
  int layer_num = 0;       /* Layer to use in single layer plot.*/
  int one_vortex = 0;      /* Plot only a single vortex */
  int vortex_num = 0;      /* Vortex to use in single vortex plot.*/
  unsigned int delay_counter = 0; /* variable for controlling frame rate */
  int frame_num = 1;
  int num_nodes = 1;
  /* This variable keeps track of the number of particle types, normally 1. */
  int ntypes = 1;
  /* CIJOL Adding functionality to do_voronoi: 1 indicates polygons, */
  /* 2 indicates Delaunay triangulation. */
  int do_voronoi = 0;  /* Toggle Voronoi triangulation mode on or off. */
  int voronoi_layer = 0; /* Which layer to use when computing voronoi const. */
  int movie_type = SMOVIE;  /* What type of movie to plot (3 supported) */

  char *display_name = NULL;
  int mode256 = 0;         /* indicating whether to use 256 or 128 colors */
  int max_color = 128;
  int samples = 5000; /* max number of points in the file ; user specified */
     /* pointer to memory block holding the x-pos of all points:*/
  float *x_pts = NULL; 
  float *y_pts = NULL; /* pointer for y */
  int   *c_pts = NULL; /* pointer to color info in file, if there is any */
  int   *p_pts = NULL; /* pointer to which point drawing */
  int do_color = 0;  /* Indicates if color available. If so, contains depth */
  int monochrome = 0;
  int old_pos[LNUM][MAX_OBJECT][2];
  int do_stripe = 0; /* toggles drawing direction of angle theta for stripes */
  float xmagnify=1.0; /* Factors by which x,y coordinates will be multiplied */
  float ymagnify=1.0;
  /* These commands handle a periodic boundary condition shift */
  int do_xshift=0;
  int do_yshift=0;
  float syssizex;
  float syssizey;
  float xshift;
  float yshift;
  float centerx,centery; /* Middle of circle for the network */
  
  char comm_line[255];
  int i, j, k, l, l1, matches;
  float uxmin,uxmax,uymin,uymax; /* x/y ranges specified by user */
  float orig_xmin,orig_xmax,orig_ymin,orig_ymax;
  float PinSize;
  int file_not_open = 1;
  FILE *execfile;
  int len;
  char filename[200];
  char last_file_name[200];
  char id_string[80];
  char temp[30];

  FILE* shell;
  char command[255];
  char command2[255];
  char strs[1000][255];
  char last_word[255];
  char* strptr;
  char c;
  int ok;
  struct termios old, new;

  error = initialize(samples,x_pts,y_pts,c_pts,p_pts,display_name,&do_color,
		     monochrome,argc,argv,execname,&do_file,&uxmin,&uxmax,
		     &uymin,&uymax,last_file_name);
  switch(error) {
  case 1:
    printf("Not enough memory to run!\n");
    exit(100);
    break;
  case 2:
    exit(200);
    break;
  }
  
  /*Process commands, first from file, then from the user:*/
  /*Clear command string */
  for (j=0; j<255; j++) comm_line[j] = '\0';

  do{

    /* If the program was passed the name of a file full of */
    /* commands to execute, process that file now: */
    if((do_file)&&(file_not_open)){
      /*Open file.  If file not found, move directly to command line*/
      if((execfile = fopen(execname,"r")) == NULL) {
	/* We were unable to open requested file */
	printf("Unknown filename!\n");
	do_file = 0;
      }
      else
	file_not_open = 0;
    }

    if(do_file){
      /*Clear command string */
      for (j=0; j<255; j++) comm_line[j] = '\0';
      /* Read single command from file, ignoring lines beginning with '#' */
      if (fgets(comm_line,80,execfile) != NULL) {
	if (comm_line[0] == '#') {
	  comm_line[0] = '\0';
	} 
	else{
	  len = strlen(comm_line);
	  *(comm_line + len-1) = '\0'; /* erase the \n charecter */
	}
      }
      else {
	/*We have hit the end of the file.  Move to command line.*/
	fclose(execfile);
	file_not_open = 1;
	do_file = 0;
      }
    }

    if(!do_file){
      /* Print prompt */
      printf("del-plot> ");
      fflush(stdout);
      /*Clear command string */
      for (j=0; j<100; j++) comm_line[j] = '\0';
      /* Read command from user.*/

    tcgetattr(0,&old);
  	new = old;
  	new.c_lflag&=~(ICANON | ECHO);
  	tcsetattr(0,TCSANOW,&new);

  	i = 0;
  	j = 0;
  	strcpy(last_word, "");

  	read(0, &c, 1);
  	if ((c != '\t') && (c != 127))
	{ 
		printf("%c", c);
		fflush(stdout);
	}

    while (c != '\n')
    {
      	if (c == '\t')
      	{
      		if (strlen(last_word) != 0)
      		{
      			last_word[j] = '\0';
      			strptr = last_word;
      			strptr += strlen(last_word) - 1;

      			while ((*strptr != '/') && (strptr != last_word))
      				strptr--;

      			if (strptr != last_word)
      			{
	      			strcpy(command, "ls ");
	      			*strptr = '\0';
	      			strcat(command, last_word);
	      			*strptr = '/';
	      			strcat(command, "| grep \"^");
	      			last_word[j] = '\0';
	      			l1 = strlen(strptr + 1);
	      			strcat(command, strptr + 1);
	      			strcat(command, "\"");
				}
				else
				{
					strcpy(command, "ls | grep \"^");
	      			last_word[j] = '\0';
	      			strcat(command, last_word);
	      			strcat(command, "\"");
	      			l1 = strlen(last_word);
				}

      			l = 0;
      			while (l < strlen(command))
      			{
      				if (command[l] == '.')
      				{	
      					memcpy(command2, command, l);
      					command2[l] = '\0';
      					strcat(command2, "\\.");
      					strptr = command;
      					strptr += l + 1;
      					strcat(command2, strptr);
      					strcpy(command, command2);
      					l++;
      				}

      				l++;
      			}
      			//printf("%s\n", command);
      			shell = popen(command, "r");

      			k = 0;
      			c = 'x';
      			
      			while (c != EOF)
      			{	
      				c = fgetc(shell);
      				l = 0;
      				strcpy(strs[k], "");
      				while ((c != '\n') && (c != EOF))
      				{
      					strs[k][l++] = c;
      					c = fgetc(shell);
      				}

      				strs[k][l] = '\0';
      				if (c != EOF)
      					k++; 
      			}
      			
      			l = strlen(last_word);

      			ok = 1;
      			
      			while ((ok) && (l1 < strlen(strs[0])))
      			{
      				c = strs[0][l1];
      				for (matches = 1; matches < k; matches++)
      					if (strs[matches][l1] != c)
      						ok = 0;

      				if (ok)
      				{
      					last_word[j++] = c;
      					comm_line[i++] = c;
      					printf("%c", c);
      					fflush(stdout);
      					l1++;
      				}
      			}
      		}
      	}
      	else
      		if ((c == '\b') || (c == 127))
      		{
      			if (i != 0)
      			{
      				i--;
      				comm_line[i] = '\0';
      				printf("\b \b"); fflush(stdout);

      				if (j != 0)
      				{
      					j--;
      					last_word[j] = '\0';
      				}
      				else
      				{
      					strptr = comm_line;
      					strptr += strlen(comm_line) - 2;
      					while ((*strptr != ' ') && (strptr != comm_line))
      						strptr--;

      					while (*strptr != '\0')
      					{
      						last_word[j++] = *strptr;
      						strptr++;
      					}
      				}
      			}
      		}
      		else
	      		if (c == ' ')
	      		{
	      			comm_line[i++] = c;
	      			strcpy(last_word, "");
	      			j = 0; 
	      		}
	      		else
	      		{
	      			comm_line[i++] = c;
	      			last_word[j++] = c;
	      		}

      	read(0, &c, 1);
      	if ((c != '\t') && (c != 127))
  	  	{ 
  	  		printf("%c", c);
  	  		fflush(stdout);
  	  	}
    }
    tcsetattr(0,TCSANOW,&old);
    }

    process_command(comm_line,plot_name,&do_plot,&stop_now,&ipos,
		&continue_plotting,&num_pins,&one_layer,&one_vortex,
		&vortex_num,&layer_num,argc,argv,
		&mode256,&max_color,&samples,x_pts,y_pts,c_pts,p_pts,
		&do_color,&monochrome,execname,&uxmin,&uxmax,&uymin,&uymax,
		&delay_counter,&frame_num,&do_S_contour,&do_T_contour,&PinSize,
		&do_file,&do_clear,&traj_on,&orig_xmin,&orig_xmax,&orig_ymin,
		&orig_ymax,&num_nodes,&do_voronoi,&voronoi_layer,&movie_type,
		    &ntypes,&do_stripe,&xmagnify,&ymagnify,&do_xshift,
		    &do_yshift,&syssizex,&syssizey,&xshift,&yshift,
		    &do_net_contour,&centerx,&centery);

    /*If the command is to plot, plot:*/
    if(do_plot){

      /*Are we setting up new file?*/
      if(do_plot == 1)  
	/* Read name of file to plot into filename*/
	get_arg(comm_line,&ipos,filename);

      do_rewind = do_fast_forward = 0;

      if (!continue_plotting) {
	/* We are beginning to plot the file.*/

	/* Did the user specify "replot" without using "plot" first? */
	/* If not, reopen previous file */
	if(do_plot == 2){
	  strcpy(filename,last_file_name);
	  if (filename[0] == '\0') {
	    printf("Need to plot a file first!\n");
	    continue;
	  }
	}
	else strcpy(last_file_name,filename);
  
	/* Open either one file or a series of files, depending on */
	/* how many nodes we are working with. */
	if(num_nodes == 1){
	  /* Single node or normal useage mode */
	  if ( (picfile[0] = fopen(filename,"r")) == NULL) {
	    printf("Unknown file\n");
	    continue;
	  }
	}
	else{
	  /* Multiple nodes.  User gave id_string in call. */
	  strcpy(id_string,filename);
	  for(j=0;j<num_nodes;j++){
	    /* Construct file name */
	    strcpy(filename,id_string);
	    sprintf(temp,"%d",j);
	    strcat(filename,"n");
	    strcat(filename,temp);

	    if((picfile[j] = fopen(filename,"r"))==NULL){
	      printf("File %s not found\n",filename);
	      continue;
	    }
	  }
	}
	
	setwindow(uxmin,uymin,uxmax,uymax,num_pins,argc,argv,
		  max_color,&do_color,monochrome,do_S_contour,
		  do_T_contour,
		  PinSize,&do_clear,traj_on,&orig_xmin,&orig_xmax,
		  &orig_ymin,&orig_ymax,do_net_contour,centerx,
		  centery,do_xshift,do_yshift,xshift,yshift,
		  xmagnify,ymagnify,syssizex,syssizey,0);
	
	do_clear = 1;
      }

      while ( plot_frame(&do_rewind,&do_fast_forward,&one_layer,
			 layer_num,one_vortex,vortex_num,
                         do_color,delay_counter,&frame_num,
			 do_S_contour,do_T_contour,
			 &do_clear,old_pos,traj_on,num_nodes,
			 movie_type,do_voronoi,voronoi_layer,ntypes,
			 do_stripe,max_color,xmagnify,ymagnify,
			 do_xshift,do_yshift,syssizex,syssizey,
			 xshift,yshift,do_net_contour) ) {
	check_events(&key_press,&do_rewind,&do_fast_forward,num_pins,
		     argc,argv,max_color,&do_color,monochrome,uxmin,
		     uxmax,uymin,uymax,frame_num,do_S_contour,
		     do_T_contour,
		     PinSize,&do_clear,traj_on,&paused,&orig_xmin,
		     &orig_xmax,&orig_ymin,&orig_ymax,do_net_contour,
		     centerx,centery,do_xshift,do_yshift,xshift,yshift,
		     xmagnify,ymagnify,syssizex,syssizey);
	if (key_press){
	  printf("You pressed a key!\n");
	  do_file = 0;
	  break;
	}
	if (paused) pause_events(num_pins,argc,argv,max_color,&do_color,
				 monochrome,&uxmin,&uxmax,&uymin,&uymax,
				 do_S_contour,do_T_contour,PinSize,&do_clear,
				 traj_on,&paused,&orig_xmin,&orig_xmax,
				 &orig_ymin,&orig_ymax,do_net_contour,
				 centerx,centery,do_xshift,do_yshift,
				 xshift,yshift,xmagnify,ymagnify,
				 syssizex,syssizey);
      }
    }

    /*If the command is to exit, exit:*/
  }while(!stop_now);

  /*Program exits*/
  free((char *)x_pts);
  free((char *)y_pts);
  if (c_pts != NULL) free((char *)c_pts);
}
/*==================== initialize ========================*/
/* Set up memory to run program, connect to X server, etc. */
/* Also initializes unbreakable rewind pointer */

initialize(int samples,float *x_pts,float *y_pts,int *c_pts,int *p_pts,
	   char *display_name,int *do_color,int monochrome,int argc,
	   char *argv[],char execname[],int *do_file,float *uxmin,
	   float *uxmax,float *uymin,float *uymax,char last_file_name[])
{
  void open_display1(int *do_color,int monochrome);

  int j;

  last_file_name = '\0';
  /* Set default values for max, min x,y */
  *uxmin = 0.0;
  *uymin = 0.0;
  *uxmax = 36.0;
  *uymax = 36.0;

  /* Unbreakable rewind pointer */
  for(j=0;j<16;j++){
    rewind_dat[j] = (int *) malloc((unsigned)(4000*sizeof(int)));
    if(rewind_dat[j] == NULL){
      printf("Out of memory for rewind_dat\n");
      exit(-1);
    }
  }
  disable_rewind = 0;

  /*Pointers to x,y,color, and current point.*/
  
  /* Assign the pointers.*/
  /* Pointers to x positions */
  if ( (x_pts = (float *) malloc((unsigned)(samples*sizeof(float)))) == NULL) {
    printf("Insufficent memory!\n");
    return(1);
  }
  
  /* Pointers to y positions */
  if ( (y_pts = (float *) malloc((unsigned)(samples*sizeof(float)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    return(1);
  }
  
  /* Pointers to colors */
  if ( (c_pts = (int *) malloc((unsigned)(samples*sizeof(int)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    free((char *)y_pts);
    return(1);
  }
    
  /* Pointers to the current point */
  if ( (p_pts = (int *) malloc((unsigned)(samples*sizeof(int)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    free((char *)y_pts);
    free((char *)c_pts);
    return(1);
  }

  /* Connect to X server */
  if ( (display=XOpenDisplay(display_name)) == NULL ) {
    fprintf(stderr, "basicwin: cannot connect to X server %s\\n",
	    XDisplayName(display_name));
    return(2);
  }

  /* Open the Display window */
  open_display1(do_color,monochrome);

  /* If there were any arguments passed to program, set variable do_file.*/
  if(argc == 2) {
    strcpy(execname,argv[1]);
    *do_file = 1;
  }

  return(0);
}

/*============================= process_command =======================*/
/* This routine takes control after main.  It receives a command line */
/* which it processes. */
void process_command(char comm_line[],char pl_file[],int *do_plot,
		     int *stop_now,int *ipos,int *continue_plotting,
		     int *num_pins,int *one_layer,int *one_vortex,
		     int *vortex_num,int *layer_num,int argc,char **argv,
		     int *mode256,int *max_color,int *samples,float *x_pts,
		     float *y_pts,int *c_pts,int *p_pts,int *do_color,
		     int *monochrome,char execname[],float *uxmin,
		     float *uxmax,float *uymin,float *uymax,
		     unsigned int *delay_counter,int *frame_num,
		     int *do_S_contour,int *do_T_contour,
		     float *PinSize,int *do_file,
		     int *do_clear,int *traj_on,float *orig_xmin,
		     float *orig_xmax,float *orig_ymin,float *orig_ymax,
		     int *num_nodes,int *do_voronoi,int *voronoi_layer,
		     int *movie_type,int *ntypes,int *do_stripe,
		     float *xmagnify,float *ymagnify,int *do_xshift,
		     int *do_yshift,float *syssizex,float *syssizey,
		     float *xshift,float *yshift,int *do_net_contour,
		     float *centerx,float *centery)
{
  int get_arg(char *line,int *ipos,char *arg);
  void setwindow(float xmin,float ymin,float xmax,float ymax,int num_pins,
		 int argc,char **argv,int max_color,int *do_color,
		 int monochrome,int do_S_contour,int do_T_contour,
		 float PinSize,int *do_clear,
		 int traj_on,float *orig_xmin,float *orig_xmax,
		 float *orig_ymin,float *orig_ymax,int do_net_contour,
		 float centerx,float centery,int do_xshift,int do_yshift,
		 float xshift,float yshift,float xmagnify,float ymagnify,
		 float syssizex,float syssizey,int redrawcontour);
  void set_color_mode(char *line,int pos,int *mode256,int *max_color);
  void set_sample(char *line,int pos,int *samples,float *x_pts,float *y_pts,
		  int *c_pts,int *p_pts);
  void set_xrange(char *line,int pos,float *uxmin,float *uxmax);
  void set_yrange(char *line,int pos,float *uymin,float *uymax);
  void set_delay(char *line,int pos,unsigned int *delay_counter);
  void load_color_map(char *line,int pos,int *max_color,int do_color);
  void load_special_contour(char *line,int pos,int *num_pins,int *do_S_contour,
			    int *do_T_contour,
			    float *PinSize);
  void load_triangle_contour(char *line,int pos,int *num_pins,
			     int *do_S_contour,int *do_T_contour,
			    float *PinSize);
  void load_network_contour(char *line,int pos,int *num_pins,
			int *do_S_contour,int *do_net_contour,
			float *PinSize,float *centerx,float *centery);
  void load_3D_contour(char *line,int pos,int *num_pins,int *do_S_contour,
		       int *do_T_contour,
		       float *PinSize);
  void print_help_screen();

  char comm[100],comm2[100];
  int command_called; /* this variable is set to true if
                       a command is called. To check for
                       unknown commands */

  /* Initialize some variables */
  *ipos = 0;
  command_called = 0;
  *do_plot = 0;
  *stop_now = 0;
  
  /*If the command line was not empty, interpret it.*/
  if (comm_line[0] != '\0') {
    get_arg(comm_line,ipos,comm);
  }
  else {
    comm[0] = '\0';
  }

  /* Now find out which command was issued. */
  /* Toggles: */

  /* one_layer = Number of layers to display */
  if(strcmp(comm,"onelayer") == 0){
    *one_layer = 1;
    if (get_arg(comm_line,ipos,comm2)) {
      *layer_num = atoi(comm2);
      command_called = 1;
    }
    else{
      printf("Enter number of layer: ");
      scanf("%d",layer_num);
      command_called = 1;
    }
  }
  if(strcmp(comm,"alllayer")==0){
    *one_layer = 0;
    command_called = 1;
  }
  
  /* one_vortex = Number of vortex to display */
  if(strcmp(comm,"onevortex")==0){
    *one_vortex = 1;
    if (get_arg(comm_line,ipos,comm2)){
      *vortex_num=atoi(comm2);
      command_called = 1;
    }
    else{
      printf("Enter number of vortex: ");
      scanf("%d",vortex_num);
      command_called=1;
    }
  }
  if(strcmp(comm,"allvortex")==0){
    *one_vortex=0;
    command_called=1;
  }


  /* Toggle trajectories */
  if (strcmp(comm,"toggle_traj") == 0) {
    command_called = 1;
    *traj_on = ((*traj_on != 0)?0:1);
    if(*traj_on)
      printf("Trajectory mode ON\n");
    else
      printf("Trajectory mode OFF\n");
    setwindow(*uxmin,*uymin,*uxmax,*uymax,*num_pins,argc,argv,
	      *max_color,do_color,*monochrome,*do_S_contour,*do_T_contour,
	      *PinSize,do_clear,*traj_on,orig_xmin,orig_xmax,
	      orig_ymin,orig_ymax,*do_net_contour,*centerx,*centery,
	      *do_xshift,*do_yshift,*xshift,*yshift,*xmagnify,*ymagnify,
	      *syssizex,*syssizey,0);
  }
  
  /* Toggle Voronoi mode */
  if (strcmp(comm,"toggle_voronoi") == 0) {
    command_called = 1;
    *do_voronoi = ((*do_voronoi != 0)?0:1);
    if(*do_voronoi){
      printf("Voronoi mode ON\n");
      if(*movie_type==TMOVIE){
	printf("Select layer for Voronoi construction: ");
	scanf("%d",voronoi_layer);
      }
      else
	*voronoi_layer=0;
    }
    else
      printf("Voronoi/Delaunay mode OFF\n");
  }

  /* Toggle Delaunay triangulation mode */
  if (strcmp(comm,"delaunay") == 0) {
    command_called = 1;
    *do_voronoi = ((*do_voronoi != 0)?0:2);
    if(*do_voronoi){
      printf("Delaunay mode ON\n");
      if(*movie_type==TMOVIE){
	printf("Select layer for Delaunay construction: ");
	scanf("%d",voronoi_layer);
      }
      else
	*voronoi_layer=0;
    }
    else
      printf("Voronoi/Delaunay mode OFF\n");
  }

  /* Shut off unbreakable rewind -- THIS IS NOT A TOGGLE */
  /* Once shut off, this can never be turned back on.  Period.  */
  if (strcmp(comm,"disable_rewind") == 0){
    disable_rewind = 1;
    command_called = 1;
    printf("Rewind disabled; cannot be reenabled\n");
  }

  /* Set number of nodes */
  if (strcmp(comm,"nodenum") == 0){
    if (get_arg(comm_line,ipos,comm2)) {
      *num_nodes = atoi(comm2);
      command_called = 1;
    }
    else{
      printf("Enter number of nodes: ");
      scanf("%d",num_nodes);
      command_called = 1;
    }
    if(*num_nodes>16){
      printf("Maximum number of nodes is 16\n");
      *num_nodes = 1;
    }
  }

  /* Set number of types of particles (default 1) */
  if (strcmp(comm,"ntypes") == 0){
    if (get_arg(comm_line,ipos,comm2)) {
      *ntypes = atoi(comm2);
      command_called = 1;
    }
    else{
      printf("Enter number of particle types: ");
      scanf("%d",ntypes);
      command_called = 1;
    }
    if(*ntypes>LNUM){
      printf("Maximum number of particle types is %d\n",LNUM);
      *ntypes = 1;
    }
  }
  
  /* Set monochrome */
  if (strcmp(comm,"monochrome") == 0) {
    command_called = 1;
    printf("Monochrome mode ON\n");
    *monochrome = 1;
    *do_color = 0;
  }

  /* Set stripe mode */
  if (strcmp(comm,"do_stripe") == 0) {
    command_called = 1;
    printf("Stripe mode ON\n");
    *do_stripe = 1;
  }

  /* Magnification commands, new in version 9 */
  /* Note that I don't know how to use the parser to read a float. */
  if (strcmp(comm,"xmagnify")==0){
    command_called = 1;
    printf("Enter magnification factor: ");
    scanf("%f",xmagnify);
    printf("X magnification factor set to %f\n",*xmagnify);
    /* This will only work with contours if we reset the window. */
    setwindow(*uxmin,*uymin,*uxmax,*uymax,*num_pins,argc,argv,
	      *max_color,do_color,*monochrome,*do_S_contour,*do_T_contour,
	      *PinSize,do_clear,*traj_on,orig_xmin,orig_xmax,
	      orig_ymin,orig_ymax,*do_net_contour,*centerx,*centery,
	      *do_xshift,*do_yshift,*xshift,*yshift,*xmagnify,*ymagnify,
	      *syssizex,*syssizey,1);
  }
  if (strcmp(comm,"ymagnify")==0){
    command_called = 1;
    printf("Enter magnification factor: ");
    scanf("%f",ymagnify);
    printf("Y magnification factor set to %f\n",*ymagnify);
    /* This will only work with contours if we reset the window. */
    setwindow(*uxmin,*uymin,*uxmax,*uymax,*num_pins,argc,argv,
	      *max_color,do_color,*monochrome,*do_S_contour,*do_T_contour,
	      *PinSize,do_clear,*traj_on,orig_xmin,orig_xmax,
	      orig_ymin,orig_ymax,*do_net_contour,*centerx,*centery,
	      *do_xshift,*do_yshift,*xshift,*yshift,*xmagnify,*ymagnify,
	      *syssizex,*syssizey,1);
  }

  /* Periodic boundary shift commands, new in version 9 */
  if(strcmp(comm,"xshift")==0){
    command_called=1;
    printf("Enter actual system size in x direction: ");
    scanf("%f",syssizex);
    printf("Enter distance to shift: ");
    scanf("%f",xshift);
    *do_xshift=1;
    /* This will only work with contours if we reset the window. */
    setwindow(*uxmin,*uymin,*uxmax,*uymax,*num_pins,argc,argv,
	      *max_color,do_color,*monochrome,*do_S_contour,*do_T_contour,
	      *PinSize,do_clear,*traj_on,orig_xmin,orig_xmax,
	      orig_ymin,orig_ymax,*do_net_contour,*centerx,*centery,
	      *do_xshift,*do_yshift,*xshift,*yshift,*xmagnify,*ymagnify,
	      *syssizex,*syssizey,1);
  }
  if(strcmp(comm,"yshift")==0){
    command_called=1;
    printf("Enter actual system size in y direction: ");
    scanf("%f",syssizey);
    printf("Enter distance to shift: ");
    scanf("%f",yshift);
    *do_yshift=1;
    /* This will only work with contours if we reset the window. */
    setwindow(*uxmin,*uymin,*uxmax,*uymax,*num_pins,argc,argv,
	      *max_color,do_color,*monochrome,*do_S_contour,*do_T_contour,
	      *PinSize,do_clear,*traj_on,orig_xmin,orig_xmax,
	      orig_ymin,orig_ymax,*do_net_contour,*centerx,*centery,
	      *do_xshift,*do_yshift,*xshift,*yshift,*xmagnify,*ymagnify,
	      *syssizex,*syssizey,1);
  }

  /* Set movie type */
  if (strcmp(comm,"smovie")==0){
    command_called = 1;
    printf("Entering smovie mode\n");
    *movie_type = SMOVIE;
    *num_nodes = 1;
  }
  if (strcmp(comm,"kmovie")==0){
    command_called = 1;
    printf("Entering kmovie mode\n");
    *movie_type = KMOVIE;
    *num_nodes = 1;
  }
  if(strcmp(comm,"tmovie")==0){
    command_called = 1;
    printf("Entering tmovie mode\n");
    *movie_type = TMOVIE;
  }
  if (strcmp(comm,"cmovie")==0){
    command_called = 1;
    printf("Entering cmovie mode\n");
    *movie_type = CMOVIE;
    *num_nodes = 1;
  }
    
  /* set = Set value of variable */
  if (strcmp(comm,"set") == 0) {
    
    if (get_arg(comm_line,ipos,comm2)) {
      if (strcmp(comm2,"colormode") == 0) 
	set_color_mode(comm_line,*ipos,mode256,max_color);
      if ((strcmp(comm2,"sample") == 0)||(strcmp(comm2,"samples")==0)) 
	set_sample(comm_line,*ipos,samples,x_pts,y_pts,c_pts,p_pts);
      if (strcmp(comm2,"xrange") == 0) 
	set_xrange(comm_line,*ipos,uxmin,uxmax);
      if (strcmp(comm2,"yrange") == 0) 
	set_yrange(comm_line,*ipos,uymin,uymax);
      if (strcmp(comm2,"delay") == 0) 
	set_delay(comm_line,*ipos,delay_counter);
      command_called = 1;
    }
    else printf("need to know what to set!\n");
  }

  /* Show help screen */
  if (strcmp(comm,"help")==0){
    print_help_screen();
    command_called = 1;
  }

  /* Executable commands: */

  /* Plot file */
  if (strcmp(comm,"plot") == 0){
    strcpy(pl_file,comm_line);
    *do_plot = 1;
    return;
  }

  /* Replot file */
  if (strcmp(comm,"replot") == 0) {
    *do_plot = 2;
    *frame_num = 1;
    return;
  }
    
  /* Load color map */
  if (strcmp(comm,"loadcolormap") == 0){
    load_color_map(comm_line,*ipos,max_color,*do_color);
    command_called = 1;
  }

  /* Load contour */
  if (strcmp(comm,"loadcontour") == 0){
    load_special_contour(comm_line,*ipos,num_pins,do_S_contour,do_T_contour,
			 PinSize);
    command_called = 1;
  }

  /* Load triangle contour */
  if (strcmp(comm,"loadtriangle") == 0){
    load_triangle_contour(comm_line,*ipos,num_pins,do_S_contour,do_T_contour,
			 PinSize);
    command_called = 1;
  }

  /* Load 3D contour (pinfile format) */
  if (strcmp(comm,"load3Dpins") == 0){
    load_3D_contour(comm_line,*ipos,num_pins,do_S_contour,do_T_contour,
		    PinSize);
    command_called = 1;
  }

  /* Load network contour */
  if (strcmp(comm,"loadcircnetwork") == 0){
    load_network_contour(comm_line,*ipos,num_pins,do_S_contour,do_net_contour,
		    PinSize,centerx,centery);
    if(*do_net_contour) *do_net_contour=2;
    command_called = 1;
  }
  if (strcmp(comm,"loadnetwork") == 0){
    load_network_contour(comm_line,*ipos,num_pins,do_S_contour,do_net_contour,
		    PinSize,centerx,centery);
    command_called = 1;
  }

  /* Read commands from file */
  if (strcmp(comm,"execute") == 0){
    /* Do not allow execution of a file from a file */
    if(*do_file){
      printf("Error: cannot execute file from file\n");
      return;
    }
    if (get_arg(comm_line,ipos,execname)) {
      *do_file = 1;
      return;
    }
    else printf("need a filename to execute!\n");
    command_called = 1;
  }

  /* Stop */
  if (strcmp(comm,"exit") == 0) command_called = 1;
  if (strcmp(comm,"quit") == 0){
    *stop_now = 1;
    return;
  }

  /* Clear */
  if (strcmp(comm,"") == 0) command_called = 1;
  if (strcmp(comm,"clear") == 0) {
    command_called = 1;
    *do_clear = 1;
  }

  /* Continue */
  if (strcmp(comm,"cont") == 0){
    *continue_plotting = 1;
    *do_plot = 2;
    return;
  }
  
  
  if (command_called == 0) printf("Unknown command\n");
}
/*============================= plot_frame ========================*/
/* This function reads in data a frame at a time. When the frame is */
/* done, it outputs it to the screen, setups up for the next frame, */
/* and returns a true, if there is a frame following this one.      */
/* CIJOL modification 10.4.97 Now calculates the desired pixel size */
/* for each vortex, and passes it to plot_object in the variables */
/* Box and Box2 (Box2 = 2*Box). */
/* 5.2.98 Attempting to add multiple modes (smovie,kmovie,tmovie) */
/* Called by: main$.*/
int plot_frame(int *do_rewind,int *do_fast_forward,int *one_layer,
	       int layer_num,int one_vortex,int vortex_num,int do_color,
	       unsigned int delay_counter,int *frame_num,int do_S_contour,
	       int do_T_contour,
	       int *do_clear,int old_pos[][MAX_OBJECT][2],int traj_on,
	       int num_nodes,int movie_type,int do_voronoi,int voronoi_layer,
	       int ntypes,int do_stripe,int max_color,float xmagnify,
	       float ymagnify,int do_xshift,int do_yshift,float syssizex,
	       float syssizey,float xshift,float yshift,int do_net_contour)
{
  void plot_object(float x,float y,int color,int p_num,int Box,int Box2,
		   int do_color,int do_clear,int old_pos[][MAX_OBJECT][2],
		   int traj_on,int layr,float theta,int do_stripe,
		   int one_vortex,int vortex_num);
  void rewind_movie(int *frame_num,int *do_rewind,int disable_rewind,
		    int num_nodes,int movie_type,int ntypes);
  void plot_Voronoi(int movie_type,int num_pars,int sidenum[],
		    int voronoi_layer,int,int);

  int sidenum[MAX_OBJECT];
  int i,j;
  static int num_pars = 0;
  long int steps;
  int md_time;
  int check_number;
  int color;
  int num_layers=1;
  int tot_num_layers;
  int ll,layr;
  int Box,Box2;
  int Boxindex[LNUM],Box2index[LNUM];
  static int timecount = 0;
  int index;
  int p_num;
  float x,y;
  int p_layr;
  float theta;
  int even;
  int par_x,par_y,dimer_x,dimer_y;
  int width,height;
  float xdimer,ydimer;
  int ii;

  /* Set up the indexes that will be used to determine pixel */
  /* sizes for plotting.  CIJOL this will probably need quite */
  /* a bit of adjustment, later. */
  /* Adding modification that can handle different particle types */
  /* in an overloated tmovie format. 8.11.03 */
  /* CIJOL Now adjusting.  Increasing the radius only every 2 layers. */
  /* CIJOL Further adjustments.  Putting maximum size limit. */
  if(movie_type == TMOVIE){
    for(ll=0;ll<LNUM;ll++){
      if((LNUM==2)||(ntypes==2)){
	Boxindex[ll] = 3 + ll*3;
      }
      else if((LNUM<10)||(ntypes<10))
	Boxindex[ll] = 2 + ll;
      else if((LNUM<25)||(ntypes<25))
	Boxindex[ll] = 2 + ll/2;
      else
	Boxindex[ll] = 2 + ll/4;
      if(Boxindex[ll]>13) Boxindex[ll]=13;
      Box2index[ll] = Boxindex[ll]*2;
    }
  }

  if (*do_rewind){
    rewind_movie(frame_num,do_rewind,disable_rewind,num_nodes,movie_type,
		 ntypes);
  }
    
  tot_num_layers = 0;
  num_layers = 1;
  for(j=0;j<num_nodes;j++){
    fread(&num_pars,sizeof(int),1,picfile[j]);
    fread(&md_time,sizeof(int),1,picfile[j]);
    if(movie_type==TMOVIE){
      fread(&num_layers,sizeof(int),1,picfile[j]);
      tot_num_layers += num_layers;
    }
    else
      tot_num_layers = 1;
  }
  
  if((*one_layer)&&(tot_num_layers<=layer_num)){
    printf("Single layer %d does not exist; plotting all %d layers\n",
	   *one_layer,tot_num_layers);
    *one_layer = 0;
  }

  for(j=0;j<num_nodes;j++){
    if((!disable_rewind)&&(!(*do_rewind))){
      rewind_dat[j][*frame_num] = num_pars;
      if((*frame_num)&&(!(*frame_num % 4000))){
	/* Allocate more memory to rewind_dat */
	rewind_dat[j] = 
	  (int *) realloc(rewind_dat[j],
			  (unsigned)((*frame_num+4000)*sizeof(int)));
	/* Check allocation */
	if(&(rewind_dat[j][*frame_num+1]) == NULL){
	  printf("Unable to allocate more memory for rewind\n");
	  disable_rewind = 1;
	}
      }
    }
  }

  if (*do_clear || (*frame_num == 1 )) {
    for(i=0; i<MAX_OBJECT; i++) {
      for(j=0;j<num_nodes;j++){
	for(ll=0;ll<num_layers;ll++){
	  index = ll + num_layers*j;
	  old_pos[index][i][0] = 2*win_width;
	  old_pos[index][i][1] = 2*win_height;
	}
      }
    }

    if (traj_on && (do_S_contour||do_T_contour||do_net_contour)) {
      XCopyArea(display,back_pixmap,traj_pixmap,traj_gc,0,0,win_width,
		win_height,0,0);     
    }
  }
  
  /* This little bit of code fast forwards 5 times faster. */
  if (*do_fast_forward) {
    for (i=0; i<5; i++) {
      /*Do this for each file */
      (*frame_num)++;
      for(j=0;j<num_nodes;j++){
	/* Go forward by layers; check the number of particles at each frame*/
	switch(movie_type){
	case SMOVIE:
	  steps = (long) num_pars*sizeof(smdata);
	  break;
	case KMOVIE:
	  steps = (long) num_pars*sizeof(kmdata);
	  break;
	case CMOVIE:
	  steps = (long) num_pars*sizeof(cmdata);
	  break;
	case TMOVIE:
	  steps = (long) num_layers*num_pars*sizeof(tmdata);
	  break;
	}
	fseek(picfile[j],steps,1);
	fread(&num_pars,sizeof(int),1,picfile[j]);
	fread(&md_time,sizeof(int),1,picfile[j]);
	if(movie_type==TMOVIE){
	  fread(&num_layers,sizeof(int),1,picfile[j]);
	  if(num_layers>LNUM) {
	    printf("Max number of layers exceeded. Alter code.\n");
	    exit(-1);
	  }
	}

	if(!disable_rewind){
	  rewind_dat[j][*frame_num] = num_pars;
	  if((*frame_num)&&(!(*frame_num % 4000))){
	    /* Allocate more memory to rewind_dat */
	    rewind_dat[j] = (int *) realloc(rewind_dat[j],
				 (unsigned)((*frame_num+4000)*sizeof(int)));
	    /* Check allocation */
	    if(&(rewind_dat[j][*frame_num+1]) == NULL){
	      printf("Unable to allocate more memory for rewind\n");
	      disable_rewind = 1;
	    }
	  }
	}
      }
    }
  }

  /* Read in the actual data on the vortices to be plotted: */
  for(j=0;j<num_nodes;j++){
    switch(movie_type){
    case SMOVIE:
      check_number = fread(smframes,sizeof(smdata),num_pars,
			   picfile[j]);
      /* Adding periodic boundary shifts here, BEFORE magnification */
      /* Adding the multiplication factor here*/
      for(ii=0;ii<num_pars;ii++){
	if(do_xshift){
	  smframes[ii].x+=xshift;
	  if(smframes[ii].x>syssizex) smframes[ii].x-=syssizex;
	  if(smframes[ii].x<0) smframes[ii].x+=syssizex;
	}
	if(do_yshift){
	  smframes[ii].y+=yshift;
	  if(smframes[ii].y>syssizey) smframes[ii].y-=syssizey;
	  if(smframes[ii].y<0) smframes[ii].y+=syssizey;
	}
	smframes[ii].x*=xmagnify;
	smframes[ii].y*=ymagnify;
      }
      break;
    case KMOVIE:
      check_number = fread(kmframes,sizeof(kmdata),num_pars,
			   picfile[j]);
      /* Adding periodic boundary shifts here, BEFORE magnification */
      /* Adding the multiplication factor here*/
      for(ii=0;ii<num_pars;ii++){
	if(do_xshift){
	  kmframes[ii].x+=xshift;
	  if(kmframes[ii].x>syssizex) kmframes[ii].x-=syssizex;
	  if(kmframes[ii].x<0) kmframes[ii].x+=syssizex;
	}
	if(do_yshift){
	  kmframes[ii].y+=yshift;
	  if(kmframes[ii].y>syssizey) kmframes[ii].y-=syssizey;
	  if(kmframes[ii].y<0) kmframes[ii].y+=syssizey;
	}
	kmframes[ii].x*=xmagnify;
	kmframes[ii].y*=ymagnify;
      }
      break;
    case CMOVIE:
      check_number = fread(cmframes,sizeof(cmdata),num_pars,
			   picfile[j]);
      /* Adding periodic boundary shifts here, BEFORE magnification */
      /* Adding the multiplication factor here*/
      for(ii=0;ii<num_pars;ii++){
	if(do_xshift){
	  cmframes[ii].x+=xshift;
	  if(cmframes[ii].x>syssizex) cmframes[ii].x-=syssizex;
	  if(cmframes[ii].x<0) cmframes[ii].x+=syssizex;
	}
	if(do_yshift){
	  cmframes[ii].y+=yshift;
	  if(cmframes[ii].y>syssizey) cmframes[ii].y-=syssizey;
	  if(cmframes[ii].y<0) cmframes[ii].y+=syssizey;
	}
	cmframes[ii].x*=xmagnify;
	cmframes[ii].y*=ymagnify;
      }
      break;
    case TMOVIE:
      for(ll=0;ll<num_layers;ll++){
	index = ll + num_layers*j;
	check_number = fread(tmframes[index],sizeof(tmdata),num_pars,
			     picfile[j]);
	/* Adding periodic boundary shifts here, BEFORE magnification */
	/* Adding the multiplication factor here*/
	for(ii=0;ii<num_pars;ii++){
	  if(do_xshift){
	    tmframes[index][ii].x+=xshift;
	    if(tmframes[index][ii].x>syssizex) tmframes[index][ii].x-=syssizex;
	    if(tmframes[index][ii].x<0) tmframes[index][ii].x+=syssizex;
	  }
	  if(do_yshift){
	    tmframes[index][ii].y+=yshift;
	    if(tmframes[index][ii].y>syssizey) tmframes[index][ii].y-=syssizey;
	    if(tmframes[index][ii].y<0) tmframes[index][ii].y+=syssizey;
	}
	  tmframes[index][ii].x*=xmagnify;
	  tmframes[index][ii].y*=ymagnify;
	}
      }
      break;
    }
    if (check_number != num_pars) {
      return 0;
    }
  }

  /* Pause for a specified amount of time */
  usleep(delay_counter);

  timecount++;
  //if(!(timecount%5))//changed this to be able to see which frame I'm at
  {
    printf("\r Frame number: %8d",*frame_num);
    fflush(stdout);
    timecount = 0;
  }

  /* Plot Voronoi construction if this mode is set. */
  if(do_voronoi){
    plot_Voronoi(movie_type,num_pars,sidenum,voronoi_layer,do_voronoi,
		 max_color);
  }

  /* CIJOL: Adding logic here to color the particles according to */
  /* frame number if in multiple layer mode.  If in one layer mode, */
  /* suppressing all vortices except those in desired layer. */
  /* ALSO Running through the layers backwards so that layer 0 is */
  /* displayed as the "top".  Easy to invert. */
  /* Finally, drawing lower layers with larger vortices (i.e., greater */
  /* pixel radius) so they can be seen from above.*/
  /* Plot ALL layers that have been read in */
  for(layr=(num_layers*num_nodes-1);layr>=0;layr--){
    if(movie_type==TMOVIE){
      /* CIJOL For multiple particle types, this is overridden below.*/
      Box = Boxindex[layr];
      Box2 = Box2index[layr];
    }
    else{
      /* For movies with only 3 particles, this is overridden below. */
      Box = 3;
      Box2 = 6;
    }
    for (i=0; i<num_pars; i++) {
      /* We need to add 1 to the layer number to get the color because */
      /* the first layer is zero, and the zeroth value of the colormap */
      /* is being used for the background. */
      /* ALSO Some vortices are actually antivortices, as indicated by */
      /* their negative id number.  Antivortices are to be plotted */
      /* with a reversed color scheme.*/
      switch(movie_type){
      case SMOVIE:
	color = 1;
	p_num = smframes[i].p_num;
	x = smframes[i].x;
	y = smframes[i].y;
	theta = 0.0;
	p_layr = layr;
	// Special: For movies containing only 3 particles
	if(3==num_pars){
	  Box = Boxindex[i];
	  Box2 = Box2index[i];
	}
	break;
      case KMOVIE:
	color = 1;
	p_num = kmframes[i].p_num;
	x = kmframes[i].x;
	y = kmframes[i].y;
	theta = kmframes[i].cum_disp;
	p_layr = layr;
	break;
      case CMOVIE:
	color = cmframes[i].color;
	p_num = cmframes[i].p_num;
	x = cmframes[i].x;
	y = cmframes[i].y;
	theta = cmframes[i].cum_disp;
	p_layr = layr;
	break;
      case TMOVIE:
	if(tmframes[layr][i].p_num<0){
	  color = (unsigned int)(num_layers*num_nodes-tmframes[layr][i].layr);
	  p_num = -tmframes[layr][i].p_num;
	}
	else{
	  color = (unsigned int)(tmframes[layr][i].layr) + 1;
	  p_num = tmframes[layr][i].p_num;
	}
	x = tmframes[layr][i].x;
	y = tmframes[layr][i].y;
	theta = 0.0;

	p_layr = tmframes[layr][i].layr;
	/* CIJOL Adjusting particle radius according to particle type. */
	if(ntypes>1){
	  Box = Boxindex[p_layr];
	  Box2 = Box2index[p_layr];
	  p_layr = layr;
	}
	break;
      }

      /* For grain plotting: */
      /* Logic takes advantage of my convention that the second particle */
      /* in each dimer pair has an artificially high ID number */
      even = p_num>=num_pars;

      /* CIJOL Disabling this */
#ifdef NOTUSING
      /* New: For Voronoi mode, the color is determined by the number */
      /* of sides the object has. */
      if(do_voronoi){
	color = sidenum[i];
	if (!color) 
	  color = 1;
      }
#endif

      /* Draw the frame either if all layers are being drawn or if */
      /* this particle is on the correct layer */
      if((movie_type != TMOVIE) || (!*one_layer) ||
	 (tmframes[layr][i].layr == layer_num)){
	/* Old functionality of one_vortex: *plot* only the single vortex. */
	/*if((!one_vortex)||(p_num==vortex_num)){*/
	{
	  /* Changed 8.12.03: Don't send id number, send absolute array */
	  /* position here. (pass i instead of p_num).  This avoids old_pos */
	  /* overflow problem in overloaded smovie format. */
	  plot_object(x,y,color,i,Box,Box2,do_color,*do_clear,
		    old_pos,traj_on,p_layr,theta,do_stripe,
		      one_vortex,vortex_num);
	  /* For grains: */
	  /* Also draw line between dimer pairs */
	  if(even){
	    /* draw a line between x,y and xdimer,ydimer */
	    /* Convert x,y location into plotting coordinates */
	    par_x = (int) ( (x + x_offset) * x_scale + BORDER );
	    par_y = (int) (((float)win_height) - (y + y_offset) * y_scale - 
			   BORDERY);
	    dimer_x = (int) ( (xdimer + x_offset) * x_scale + BORDER );
	    dimer_y = (int) (((float)win_height) - 
			     (ydimer + y_offset) * y_scale - BORDERY);
	    width  = ABS((dimer_x - par_x));
	    height = ABS((dimer_y - par_y));
    
	    if (((width<(win_width/3)) && (height<(win_height/3))) 
		&& ((width>0) || (height>0))) {
	      XSetForeground(display,pixmap_gc,
			     (BlackPixel(display,screen_num)));
	      XDrawLine(display,pixmap,pixmap_gc,
			(dimer_x),(dimer_y),(par_x),(par_y));
	    }
	    /* CIJOL Altering to draw line along length of chain: */
	    /* the following two lines accomplish this. */
	    xdimer = x;
	    ydimer = y;
	  }
	  else{
	    xdimer=x;
	    ydimer=y;
	  }
	  /* end of grain segment */
	}
      }
    }
  }

  if (num_pars) {
    XCopyArea(display,pixmap,win,gc,0,0,win_width,win_height,0,0);
    XFlush(display);
    (*frame_num)++;
  }

  *do_clear = 0; /* if do_clear was set, it doesn't need to be anymore. */

  if (traj_on) {
    XCopyArea(display,traj_pixmap,pixmap,gc,0,0,win_width,win_height,0,0);
  }
  else if (do_S_contour||do_T_contour||do_net_contour) {
    XCopyArea(display,back_pixmap,pixmap,gc,0,0,win_width,win_height,0,0);
  }

  else if ( (!traj_on) && (!do_S_contour) && (!do_T_contour) && 
	    (!do_net_contour)){
    if (do_color) {
      XSetForeground(display,pixmap_gc,c_map[0]);
    }
    else {
      XSetForeground(display,pixmap_gc,(WhitePixel(display,screen_num)));
    }
	
    XFillRectangle(display,pixmap,pixmap_gc,0,0,win_width,win_height);
    *do_clear = 0;
  }
    
  return num_pars;
}
/*========================= plot_object ==================*/
/* This function actually draws the objects on the screen. */
/* It is also responsible for drawing the trajectories.    */
/* Called by: plot_frame$.*/
void plot_object(float x,float y,int color,int p_num,int Box,int Box2,
		 int do_color,int do_clear,int old_pos[][MAX_OBJECT][2],
		 int traj_on,int layr,float theta,int do_stripe,
		 int one_vortex,int vortex_num)
{
  int par_x,par_y,old_x,old_y;
  int par_x2,par_y2;
  int width,height;
  unsigned int win_Box2;
  float x2,y2;
  float linelength;

  /* CIJOL Hardwiring linelength */
  linelength = 2.0;

  win_Box2 = Box2;

  /* If frame_num = 0, or we are clearing the screen, the previous */
  /* trajectory info is invalid. So we set the old position data   */
  /* beyond the screen boundaries. Thus no path is drawn.          */

  /* Convert x,y location into plotting coordinates */
  par_x = (int) ( (x + x_offset) * x_scale + BORDER );
  par_y = (int) (((float)win_height) - (y + y_offset) * y_scale - BORDERY);
  
  /* Clip particles that fall outside plotting region */
  
  if (par_x < (BORDER - 10))    return;
  if (par_y < (BORDERY - 10))   return;
  if (par_x >= ((int) win_width) - (BORDER - 10))    return;
  if (par_y >= ((int) win_height) - (BORDERY - 10))  return;

  /* Now we restore the background under the objects */
  /* If the window has been cleared, no need to do a */
  /* restore.                                        */

  if ( (traj_on) && (!do_clear) ) {
    old_x = old_pos[layr][p_num][0];
    old_y = old_pos[layr][p_num][1];
    XCopyArea(display,traj_pixmap,pixmap,gc,old_x-Box,old_y-Box,win_Box2,
	      win_Box2,old_x-Box,old_y-Box); 
  }
  
  /* Now draw the trajectory info. We do this before drawing the */
  /* particle because the trajectory is UNDER the particle.      */
  
  /* Our first task is to find where the particle was previously. */
  old_x = old_pos[layr][p_num][0];
  old_y = old_pos[layr][p_num][1];          
    
  /* What if a particle wraps around the screen, due to periodic  */
  /* boundary conditions? We DON'T want lines wrapping around the */
  /* screen. This code checks for position changes that are too   */
  /* large, and prevents the trajectory being drawn, if they are. */
  
  if (traj_on) {
    /* CIJOL New functionality for one_vortex: plot trails for only */
    /* that vortex. 12.11.03 */
    if((!one_vortex)||(p_num==vortex_num)){
      width  = ABS((old_x - par_x));
      height = ABS((old_y - par_y));
      
      if (((width<(win_width/3)) && (height<(win_height/3))) 
	  && ((width>0) || (height>0))) {
	XSetForeground(display,traj_gc,(BlackPixel(display,screen_num)));
	XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));
	XDrawLine(display,pixmap,pixmap_gc,(old_x),(old_y),(par_x),(par_y));
	XDrawLine(display,traj_pixmap,traj_gc,(old_x),(old_y),(par_x),(par_y));
      }
    }
  }
  
  /* Now draw the particles. For now, they're just spheres. */
  
  if (color > 125) color = 125;
  
  // CIJOL Messing with this 3.10.05
  //  if (do_color) {
  //    /*XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));*/
  //XSetForeground(display,pixmap_gc,c_map[1]);
  //XDrawArc(display,pixmap,pixmap_gc,par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
  //XFillArc(display,pixmap,pix_gc[color],par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
  //}
  if (do_color) {
    /*XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));*/
    //3.11.05
    XSetForeground(display,pix_gc[color],c_map[color]);
    XDrawArc(display,pixmap,pix_gc[color],par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
    XFillArc(display,pixmap,pix_gc[color],par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
  }
  else {
    XSetForeground(display,pixmap_gc,(WhitePixel(display,screen_num)));
    XDrawArc(display,pixmap,pixmap_gc,par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
    XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));
    XFillArc(display,pixmap,pixmap_gc,par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
  }
  
  /* Store the position of the particle, so we can draw trajectories. */
  
  old_pos[layr][p_num][0] = par_x;
  old_pos[layr][p_num][1] = par_y;

  /* Lastly, draw a line indicating the direction of theta for this */
  /* particle, if do_stripe mode is set. */
  if(do_stripe){
    x2 = x+linelength*(float)cos((double)theta);
    y2 = y+linelength*(float)sin((double)theta);
    par_x2 = (int) ((x2+x_offset) * x_scale + BORDER);
    par_y2 = (int) (((float)win_height) - (y2 + y_offset) * y_scale - BORDERY);

    XSetForeground(display,traj_gc,(BlackPixel(display,screen_num)));
    XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));
    XDrawLine(display,pixmap,pixmap_gc,(par_x2),(par_y2),(par_x),(par_y));
    XDrawLine(display,traj_pixmap,traj_gc,(par_x2),(par_y2),(par_x),(par_y));
  }
}

/*================================= get_arg =====================*/
/* this function is THE KEY i/o function! It takes a command */
/* line (typed by the user), a position in the buffer, and   */
/* a buffer to store the arguement it finds. It returns a    */
/* space seperated word from the buffer, and it updates the  */
/* current buffer position. If it can't find an arguement it */
/* returns a false value to the caller. It is designed to    */
/* take the " symbol as a command to include all charecters  */
/* in the arg buffer until the closing quote is found.       */
/* Called by: get_command$, do_modifier$,set_comm$,          */
/*            set_pen$, set_color$, set_color_mode$,         */
/*            set_sample$, set_xrange$, set_yrange$,         */
/*            set_delay$, load_color_map$     */
/*            file_color$, execute$*/
/*            plot_pic$*/
int get_arg(char *line,int *ipos,char *arg)
{
  char temp[100];
  int j,include_all;
  char ch;
  
  include_all = 0;
  j = 0;
  strcpy(temp,line);
  
  /*Read through command line starting from position i.  Find a word.*/
  do {
    ch = line[(*ipos)];
    if (ch == '"') {  
      if (include_all == 0) {
	include_all = 1;
      }
      else {
	include_all = 0;
      }
    }
    if ( ( ((ch != ' ') && (ch != '\0')) || (include_all) ) && (ch != '"')) {
      arg[j] = ch;
      j++;
    }
    (*ipos)++;
  } while ( ( (!((ch == ' ') && (j != 0))) || (include_all)) && (ch != '\0') );
  arg[j] = '\0';
  
  if (j == 0) 
    return(0); /* the arguement was NOT found */
  else
    return(1); /* an arguement was found */
}

/*======================== get_sample_mem ======================*/
/* this routine allocates the necessary space to read */
/* in the info for a file. It also deallocates mem    */
/* previously allocated. If memory routine worked it  */
/* returns a 0.                                       */
/* Called by: set_sample$, prep_colors_in_file$ */

int get_sample_mem(int samples,float *x_pts,float *y_pts,int *c_pts,
		   int *p_pts)
{

  /*Pointers to x,y,color, and current point.*/
  /* Free the already assigned pointers. */
  free((char *)x_pts);
  free((char *)y_pts);
  free((char *)c_pts);
  free((char *)p_pts);
  
  /* Now assign the pointers.*/
  if ( (x_pts = (float *) malloc((unsigned)(samples*sizeof(float)))) == NULL) {
    printf("Insufficent memory!\n");
    return(1);
  }
  
  if ( (y_pts = (float *) malloc((unsigned)(samples*sizeof(float)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    return(1);
  }
  
  if ( (c_pts = (int *) malloc((unsigned)(samples*sizeof(int)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    free((char *)y_pts);
    return(1);
  }
  
  if ( (p_pts = (int *) malloc((unsigned)(samples*sizeof(int)))) == NULL) {
    printf("Insufficent memory!\n");
    free((char *)x_pts);
    free((char *)y_pts);
    free((char *)c_pts);
    return(1);
  }

  return(0);
}

/*======================= open_display =====================*/
/* This routine opens a display window.*/
/* Called by: main$ */
void open_display1(int *do_color,int monochrome)
{

  Visual *visual;
  int depth;

  /* get screen_num size from display structure macro */
  screen_num = DefaultScreen(display);
  screen_ptr = DefaultScreenOfDisplay(display);
  
  border_pixel     = BlackPixel(display,screen_num);
  background_pixel = BlackPixel(display,screen_num);
  foreground_pixel = border_pixel;
  
  visual = DefaultVisual(display,screen_num);
  
  /* CIJOL TEST 12.18.01 Removing this line. */
  /* CIJOL This fixes the color! */
  /*if (visual->class != PseudoColor) {
    *do_color = 0;
  }
  else {*/
  {
    depth = DefaultDepth(display,screen_num);
    if (depth == 1) {
      *do_color = 0;
    }
    else {
      *do_color = depth;
    }
  }
  
  if (monochrome) *do_color = 0;
  /* DEBUG */
  printf("Depth: %d\n",*do_color);
}

/*========================== setwindow =======================*/
/* this procedure sets up a window if it hasn't been already */
/* if it has been opened, it resets the ranges if user has   */
/* specified new ones.                                       */
/* Called by: get_command$,plot_pic$,check_events$,pause_events$*/
void setwindow(float xmin,float ymin,float xmax,float ymax,int num_pins,
	       int argc,char **argv,int max_color,int *do_color,int monochrome,
	       int do_S_contour,int do_T_contour,
	       float PinSize,int *do_clear,int traj_on,
	       float *orig_xmin,float *orig_xmax,float *orig_ymin,
	       float *orig_ymax,int do_net_contour,float centerx,float centery,
	       int do_xshift,int do_yshift,float xshift,float yshift,
	       float xmagnify,float ymagnify,float syssizex,float syssizey,
	       int redrawcontour)
{
  void open_win(float xmin,float ymin,float xmax,float ymax,int argc,
		char *argv[],int max_color,int *do_color,int monochrome);
  void draw_contour(int num_pins,float PinSize,int do_xshift,int do_yshift,
		    float xshift,float yshift,float xmagnify,float ymagnify,
		    float syssizex,float syssizey);
  void draw_triangle_contour(int num_pins,int do_xshift,int do_yshift,
			     float xshift,float yshift,float xmagnify,
			     float ymagnify,float syssizex,float syssizey);
  void draw_network_contour(int num_pins,float,float,int do_xshift,int
			    do_yshift,float xshift,float yshift,float xmagnify,
			    float ymagnify,float syssizex,float syssizey);
  void draw_network_contourcirc(int num_pins,float,float,int do_xshift,int
			    do_yshift,float xshift,float yshift,float xmagnify,
			    float ymagnify,float syssizex,float syssizey);
  void get_GC(GC *gc, XFontStruct *font_info);
  void get_traj_GC(GC *gc,XFontStruct *font_info);

  float x_size,y_size;
  unsigned int depth;
  int size_change,range_change;
  static unsigned int old_width  = 0;
  static unsigned int old_height = 0;
  static float old_xmin,old_xmax,old_ymin,old_ymax;
  XWindowAttributes *win_atb;

  if (old_width == 0) {
    old_xmin = xmin;
    old_xmax = xmax;
    old_ymin = ymin;
    old_ymax = ymax;
    
    *orig_xmin = xmin;
    *orig_ymin = ymin;
    *orig_xmax = xmax;
    *orig_ymax = ymax;
  }

  open_win(xmin,ymin,xmax,ymax,argc,argv,max_color,do_color,monochrome);

  win_atb=(XWindowAttributes *)malloc((unsigned)(sizeof(XWindowAttributes)));
  XGetWindowAttributes(display,win,win_atb);
  win_width  = win_atb->width;
  win_height = win_atb->height;
  
  depth = DefaultDepth(display,screen_num);
  
  x_offset   = 0 - xmin;
  y_offset   = 0 - ymin;

  x_size     = xmax - xmin;
  y_size     = ymax - ymin;

  x_scale    = (win_width-2*BORDER)/x_size;
  y_scale    = (win_height-2*BORDERY)/y_size;

  size_change = !( (old_width == win_width) && (old_height == win_height) );
  range_change = !((old_xmin==xmin)&&(old_xmax==xmax)&&(old_ymin==ymin)&&(old_ymax == ymax));
  
  /* CIJOL This will be repeated not only for size changes but also */
  /* for other alterations. 4.3.05. */
  if ( size_change || do_xshift || do_yshift || redrawcontour)
    {
    pixmap = XCreatePixmap(display,win,win_width,win_height,depth);
    get_GC(&(pixmap_gc),font_info);  
    
    back_pixmap = XCreatePixmap(display,win,win_width,win_height,depth);
    get_GC(&(back_gc),font_info);  

    traj_pixmap = XCreatePixmap(display,win,win_width,win_height,depth);
    get_traj_GC(&(traj_gc),font_info);  

    if(do_T_contour)
      draw_triangle_contour(num_pins,do_xshift,do_yshift,xshift,yshift,
			    xmagnify,ymagnify,syssizex,syssizey);
    else if(2==do_net_contour)
      draw_network_contourcirc(num_pins,centerx,centery,do_xshift,do_yshift,
			   xshift,yshift,xmagnify,ymagnify,syssizex,
			   syssizey);
    else if(1==do_net_contour)
      draw_network_contour(num_pins,centerx,centery,do_xshift,do_yshift,
			   xshift,yshift,xmagnify,ymagnify,syssizex,
			   syssizey);
    else
      draw_contour(num_pins,PinSize,do_xshift,do_yshift,xshift,yshift,
		   xmagnify,ymagnify,syssizex,syssizey);

    old_width  = win_width;
    old_height = win_height;
  }

  if (range_change) {
    old_xmin = xmin;
    old_xmax = xmax;
    old_ymin = ymin;
    old_ymax = ymax;
  }

  if ( (!(size_change)) && range_change){ 
    if(do_T_contour)
      draw_triangle_contour(num_pins,do_xshift,do_yshift,xshift,yshift,
			    xmagnify,ymagnify,syssizex,syssizey);
    else if(2==do_net_contour)
      draw_network_contourcirc(num_pins,centerx,centery,do_xshift,do_yshift,
			   xshift,yshift,xmagnify,ymagnify,syssizex,
			   syssizey);
    else if(1==do_net_contour)
      draw_network_contour(num_pins,centerx,centery,do_xshift,do_yshift,
			   xshift,yshift,xmagnify,ymagnify,syssizex,
			   syssizey);
    else
      draw_contour(num_pins,PinSize,do_xshift,do_yshift,xshift,yshift,
		   xmagnify,ymagnify,syssizex,syssizey);
  }

  if (traj_on && (do_S_contour||do_T_contour||do_net_contour)) {
    
    XCopyArea(display,back_pixmap,traj_pixmap,traj_gc,0,0,
	      win_width,win_height,0,0);     

    XCopyArea(display,back_pixmap,pixmap,pixmap_gc,0,0,
	      win_width,win_height,0,0);     
  }

  if (traj_on && (!do_S_contour) && (!do_T_contour)&&(!do_net_contour)) {

    if (*do_color) {
      XSetForeground(display,traj_gc,c_map[0]);
    }
    else {
      XSetForeground(display,traj_gc,(WhitePixel(display,screen_num)));    
    }
    XFillRectangle(display,traj_pixmap,traj_gc,0,0,win_width,win_height);

    *do_clear = 1;
  }

  if (*do_color) {
    XSetForeground(display,gc,c_map[0]);
  }
  else {
    XSetForeground(display,gc,(WhitePixel(display,screen_num)));    
  }

  XFillRectangle(display,win,gc,0,0,win_width,win_height);

  *do_clear = 1;

  free((char *)(win_atb));
}/*end setwindow*/

/*========================== open_win ========================*/
/* Opens a window. */
/* Called by: setwindow$. */
void open_win(float xmin,float ymin,float xmax,float ymax,int argc,
	      char *argv[],int max_color,int *do_color,int monochrome)
{
  void load_font(XFontStruct **font_info);
  void get_GC(GC *gc, XFontStruct *font_info);
  void get_xorGC(GC *gc,XFontStruct *font_info);
  void make_boxes(int max_color,int do_color);

  int x,y;
  float x_size,y_size;
  unsigned int width, height;      /* window size and position */
  unsigned int borderwidth = 4;	 /* four pixels */
  unsigned int display_width, display_height;
  /* unsigned int icon_width, icon_height; */
  char window_name[128];
  char *icon_name = "x-plot";
  Pixmap icon_pixmap;
  XSizeHints size_hints;
  XEvent report;
  Visual *visual;

  static int opened = 0; /* Indicating if window already opened */

  strcpy(window_name,"Movie");
  
  /* get screen_num size from display structure macro */
  screen_num = DefaultScreen(display);
  screen_ptr = DefaultScreenOfDisplay(display);
  
  display_width  = DisplayWidth(display, screen_num) - 100;
  display_height = DisplayHeight(display, screen_num) - 100;
  
  /* place window */
  x = 0, y = 0;
  
  /* size window with enough room for text */

  x_size = xmax - xmin;
  y_size = ymax - ymin;

  if (x_size > y_size) {
    width  = display_width;
    height = (display_width*y_size)/x_size;
    
    if (height > display_height) {
      width  = (display_width*display_height)/height;
      height = display_height;
    }
  }
  else {
    height = display_height;
    width  = display_height*(x_size/y_size);
    
    if (width > display_width) {
      height = display_width*(display_width/width);
      width  = display_width;
    }
  }

  if ( (height == win_height) && (width == win_width) ) return;
  
  border_pixel     = BlackPixel(display,screen_num);
  background_pixel = BlackPixel(display,screen_num);
  foreground_pixel = border_pixel;
  
  visual = DefaultVisual(display,screen_num);
  
  /* CIJOL TEST 12.19.01 Removing this */
  /* CIJOL This fixes the color! */
  /*if (visual->class != PseudoColor) {
    *do_color = 0;
  }
  else {*/
  {
    int depth;
    depth = DefaultDepth(display,screen_num);
    if (depth == 1) {
      *do_color = 0;
    }
    else {
      *do_color = depth;
    }
  }
  
  if (monochrome) *do_color = 0;
  
  win = XCreateSimpleWindow(display,(RootWindow(display,screen_num)), x, y,
			    width, height, borderwidth, border_pixel,
			    background_pixel);
  
  icon_pixmap = XCreateBitmapFromData(display,win,icon_bits,icon_width,icon_height);
  
  /* Set resize hints */
  size_hints.flags = PPosition | PSize | PMinSize;
  size_hints.x = x;
  size_hints.y = y;
  size_hints.width = width;
  size_hints.height = height;
  size_hints.min_width = 2*BORDER;
  size_hints.min_height = 2*BORDER;
  
  /* set Properties for window manager (always before mapping) */
  XSetStandardProperties(display, win, window_name, icon_name, 
			 icon_pixmap, argv, argc, &size_hints);
  
  /* Select event types wanted */
  XSelectInput(display, win, EVENT_MASK);
  
  load_font(&font_info);
  
  /* create GC for text and drawing */
  get_GC(&gc, font_info);
  get_xorGC(&xor_gc, font_info);
  
  if (!opened) {
    XMapWindow(display, win);
    XNextEvent(display, &report); 
    opened = 1;
    make_boxes(max_color,*do_color); 
  }
}
/*======================== draw_contour ================*/
/* This routine takes an ASCII special_contour file, which lists*/
/* locations and sizes of parabolic pinning sites, and draws them*/
/* as open circles on the screen.  */
/* CIJOL Dropping use of variable PinSize, instead taking the size */
/* of the pinning site directly out of the stored structure. */
/* Called by: setwindow$*/
void draw_contour(int num_pins,float PinSize,int do_xshift,int do_yshift,
		  float xshift,float yshift,float xmagnify,float ymagnify,
		  float syssizex,float syssizey)
{
  float x,y;
  int par_x,par_y;
  int i;
  int pin_radiusx,pin_radiusy;
  float LPinSizex;
  float LPinSizey;

  XSetForeground(display,back_gc,(WhitePixel(display,screen_num)));    
  XFillRectangle(display,back_pixmap,back_gc,0,0,win_width,win_height); 

  XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));    
  for(i=0;i<num_pins;i++){
    x = pin_sites[i].x;
    y = pin_sites[i].y;
    /* Take into account possible position shifts */
    if(do_xshift){
      x+=xshift;
      if(x>syssizex) x-=syssizex;
      if(x<0) x+=syssizex;
    }
    if(do_yshift){
      y+=yshift;
      if(y>syssizey) y-=syssizey;
      if(y<0) y+=syssizey;
    }
    x*=xmagnify;
    y*=ymagnify;
    /*CIJOL NOTE I have not attempted to adjust radius to the rescaling. */
    LPinSizex = pin_sites[i].radiusx*2.0;
    LPinSizey = pin_sites[i].radiusy*2.0;

    par_x = (int) ( (x + x_offset) * x_scale + BORDER );
    par_y = (int) ((float)win_height - (y + y_offset) * y_scale - BORDERY);
    
    if (par_x >= (int) win_width - (BORDER  - 10))   continue;
    if (par_y >= (int) win_height - (BORDERY - 10))  continue;
    if (par_y < (BORDERY - 10))              continue;
    if (par_x < (BORDER - 10))               continue;
    
    /*In this call, to alter the size of the pinning site which is*/
    /*drawn, alter the integer variable pin_radius.*/
    /*To draw the pin site to scale, use these lines:*/
    {
      /* CIJOL For pins all the same size */
      /*pin_radius = (int)(PinSize * x_scale);*/
      /* For pins of varying size */
      pin_radiusx = (int)(LPinSizex * x_scale);
      pin_radiusy = (int)(LPinSizey * x_scale);
      if(pin_radiusx <= 0) pin_radiusx = 1;
      if(pin_radiusy <= 0) pin_radiusy = 1;

    }
    /*Otherwise, use a fixed radius:*/
    /*pin_radius = 8;*/
    /*XDrawArc(display,back_pixmap,back_gc,par_x-4,par_y-4,
      pin_radius,pin_radius,0,360*64);*/

    //this is for circular pins used in the sparse
    XDrawArc(display,back_pixmap,back_gc,par_x-pin_radiusx/2.0,par_y-pin_radiusy/2.0,pin_radiusx,pin_radiusy,0,360*64);
/*XDrawArc(display,back_pixmap,back_gc,par_x-pin_radius/2,par_y-pin_radius/2,
	     pin_radius,pin_radius,0,360*64);*/
  }
}
/*======================== draw_triangle_contour ================*/
/* Version of draw_contour for triangular pinning sites, which are */
/* drawn quite differently from the circles. */
/* Called by: setwindow$*/
void draw_triangle_contour(int num_pins,int do_xshift,int do_yshift,
			   float xshift,float yshift,float xmagnify,
			   float ymagnify,float syssizex,float syssizey)
{
  float x1,x2,x3,y1,y2,y3;
  int j;
  /* CIJOL: type XPoint consists of pairs of shorts, x y. */
  XPoint points[3];

  XSetForeground(display,back_gc,(WhitePixel(display,screen_num)));    
  XFillRectangle(display,back_pixmap,back_gc,0,0,win_width,win_height); 

  XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));    
  for(j=0;j<num_pins;j++){
    x1=pin_sites[j].px[0];
    y1=pin_sites[j].py[0];
    x2=pin_sites[j].px[1];
    y2=pin_sites[j].py[1];
    x3=pin_sites[j].px[2];
    y3=pin_sites[j].py[2];

    /* Take into account possible position shifts. */
    if(do_xshift){
      x1+=xshift;
      x2+=xshift;
      x3+=xshift;
      if(x1>syssizex) x1-=syssizex;
      if(x2>syssizex) x2-=syssizex;
      if(x3>syssizex) x3-=syssizex;
      if(x1<0) x1+=syssizex;
      if(x2<0) x2+=syssizex;
      if(x3<0) x3+=syssizex;
    }
    if(do_yshift){
      y1+=yshift;
      y2+=yshift;
      y3+=yshift;
      if(y1>syssizey) y1-=syssizey;
      if(y2>syssizey) y2-=syssizey;
      if(y3>syssizey) y3-=syssizey;
      if(y1<0) y1+=syssizey;
      if(y2<0) y2+=syssizey;
      if(y3<0) y3+=syssizey;
    }
    x1*=xmagnify;
    x2*=xmagnify;
    x3*=xmagnify;
    y1*=xmagnify;
    y2*=xmagnify;
    y3*=xmagnify;

    /* Convert x,y location into plotting coordinates */
    points[0].x = (short) ( (x1 + x_offset) * x_scale + BORDER );
    points[0].y = (short) (((float)win_height)-(y1+y_offset)*y_scale -BORDERY);
    points[1].x = (short) ( (x2 + x_offset) * x_scale + BORDER );
    points[1].y = (short) (((float)win_height)-(y2+y_offset)*y_scale -BORDERY);
    points[2].x = (short) ( (x3 + x_offset) * x_scale + BORDER );
    points[2].y = (short) (((float)win_height)-(y3+y_offset)*y_scale -BORDERY);
  
    /* Clip lines that fall outside plotting region */
    if (points[0].x < (BORDER - 10))    continue;
    if (points[0].y < (BORDERY - 10))   continue;
    if (points[0].x >= ((int) win_width) - (BORDER - 10))    continue;
    if (points[0].y >= ((int) win_height) - (BORDERY - 10))  continue;
    if (points[1].x < (BORDER - 10))    continue;
    if (points[1].y < (BORDERY - 10))   continue;
    if (points[1].x >= ((int) win_width) - (BORDER - 10))    continue;
    if (points[1].y >= ((int) win_height) - (BORDERY - 10))  continue;
    if (points[2].x < (BORDER - 10))    continue;
    if (points[2].y < (BORDERY - 10))   continue;
    if (points[2].x >= ((int) win_width) - (BORDER - 10))    continue;
    if (points[2].y >= ((int) win_height) - (BORDERY - 10))  continue;

    /* For black lines, use next line only: */
    XDrawLine(display,back_pixmap,back_gc,(points[0].x),(points[0].y),
	      (points[1].x),(points[1].y));
    XDrawLine(display,back_pixmap,back_gc,(points[0].x),(points[0].y),
	      (points[2].x),(points[2].y));
    XDrawLine(display,back_pixmap,back_gc,(points[2].x),(points[2].y),
	      (points[1].x),(points[1].y));
  }
}
/*======================== draw_network_contourcirc ================*/
/* Version of draw_contour for network sites */
/* Called by: setwindow$*/
/* I try drawing "V" shapes instead of circles. */
/* For circle geometry */
void draw_network_contourcirc(int num_pins,float centerx,float centery,
			  int do_xshift,int do_yshift,float xshift,
			  float yshift,float xmagnify,float ymagnify,
			  float syssizex,float syssizey)
{
  float x1,x2,y1,y2;
  float dx,dy,rx,ry,dr;
  float cx,cy;
  float radius;
  int j;
  int lBox,lBox2;
  float leftx,rightx,lefty,righty;
  /* CIJOL: type XPoint consists of pairs of shorts, x y. */
  XPoint points[3];

  XSetForeground(display,back_gc,(WhitePixel(display,screen_num)));    
  XFillRectangle(display,back_pixmap,back_gc,0,0,win_width,win_height); 

  XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));    
  for(j=0;j<num_pins;j++){
    /* These are the coordinates of the two connected nodes. */
    x1=pin_sites[j].px[0];
    y1=pin_sites[j].py[0];
    x2=pin_sites[j].px[1];
    y2=pin_sites[j].py[1];

    /* For drawing a circle centered halfway between, find circle center. */
    /* Note, all coordinates are ASSUMED positive, as usual. */
    dx=0.5*(x1-x2);
    dy=0.5*(y1-y2);
    rx=x2+dx;
    ry=y2+dy;

    /* I now try plotting a "V" shape instead.  For this I need to find */
    /* a point on the same radial line as the circle center but shifted */
    /* inward by some fixed amount. */
    dx=rx-centerx;
    dy=ry-centery;
    dr=(float)sqrt((double)(dx*dx+dy*dy));
    cx=rx-(dx/dr)*0.4;
    cy=ry-(dy/dr)*0.4;

    /* Take into account possible position shifts */
    if(do_xshift){
      x1+=xshift;
      cx+=xshift;
      x2+=xshift;
      if(x1>syssizex) x1-=syssizex;
      if(cx>syssizex) cx-=syssizex;
      if(x2>syssizex) x2-=syssizex;
      if(x1<0) x1+=syssizex;
      if(cx<0) cx+=syssizex;
      if(x2<0) x2+=syssizex;
    }
    if(do_yshift){
      y1+=yshift;
      cy+=yshift;
      y2+=yshift;
      if(y1>syssizey) y1-=syssizey;
      if(cy>syssizey) cy-=syssizey;
      if(y2>syssizex) y2-=syssizey;
      if(y1<0) y1+=syssizey;
      if(cy<0) cy+=syssizey;
      if(y2<0) y2+=syssizey;
    }
    x1*=xmagnify;
    cx*=xmagnify;
    x2*=xmagnify;
    y1*=ymagnify;
    cy*=ymagnify;
    y2*=ymagnify;
    /* Convert x,y location into plotting coordinates */
    points[0].x = (short) ( (x1 + x_offset) * x_scale + BORDER );
    points[0].y = (short) (((float)win_height)-(y1+y_offset)*y_scale -BORDERY);
    points[1].x = (short) ( (cx + x_offset) * x_scale + BORDER );
    points[1].y = (short) (((float)win_height)-(cy+y_offset)*y_scale -BORDERY);
    points[2].x = (short) ( (x2 + x_offset) * x_scale + BORDER );
    points[2].y = (short) (((float)win_height)-(y2+y_offset)*y_scale -BORDERY);
  
    /* Clip lines that fall outside plotting region */
    if (points[0].x < (BORDER - 10))    continue;
    if (points[0].y < (BORDERY - 10))   continue;
    
    XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));
    XDrawLine(display,back_pixmap,back_gc,points[0].x,points[0].y,
	      points[1].x,points[1].y);
    XDrawLine(display,back_pixmap,back_gc,points[1].x,points[1].y,
	      points[2].x,points[2].y);

    /* Unused from circles. */
    /* Now, figure out the local box size. */
    //radius=(float)sqrt(dx*dx+dy*dy);
    //lBox2=(int)(radius*2*x_scale);
    //if(lBox2<=0) lBox2=1;
    //lBox=lBox2/2;
    //XDrawArc(display,back_pixmap,back_gc,points.x-lBox,points.y-lBox,
    //     lBox2,lBox2,(int)startangle*64,(int)endangle*64); 
  }
}
/*======================== draw_network_contour ================*/
/* Version of draw_contour for network sites */
/* Called by: setwindow$*/
/* For linear geometry */
void draw_network_contour(int num_pins,float centerx,float centery,
			  int do_xshift,int do_yshift,float xshift,
			  float yshift,float xmagnify,float ymagnify,
			  float syssizex,float syssizey)
{
  float x1,x2,y1,y2;
  int j;
  /* CIJOL: type XPoint consists of pairs of shorts, x y. */
  XPoint points[2];

  XSetForeground(display,back_gc,(WhitePixel(display,screen_num)));    
  XFillRectangle(display,back_pixmap,back_gc,0,0,win_width,win_height); 

  XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));    
  for(j=0;j<num_pins;j++){
    /* These are the coordinates of the two connected nodes. */
    x1=pin_sites[j].px[0];
    y1=pin_sites[j].py[0];
    x2=pin_sites[j].px[1];
    y2=pin_sites[j].py[1];

    /* In this case, simply connect the dots. */

    /* Take into account possible position shifts */
    if(do_xshift){
      x1+=xshift;
      x2+=xshift;
      if(x1>syssizex) x1-=syssizex;
      if(x2>syssizex) x2-=syssizex;
      if(x1<0) x1+=syssizex;
      if(x2<0) x2+=syssizex;
    }
    if(do_yshift){
      y1+=yshift;
      y2+=yshift;
      if(y1>syssizey) y1-=syssizey;
      if(y2>syssizex) y2-=syssizey;
      if(y1<0) y1+=syssizey;
      if(y2<0) y2+=syssizey;
    }
    x1*=xmagnify;
    x2*=xmagnify;
    y1*=ymagnify;
    y2*=ymagnify;
    /* Convert x,y location into plotting coordinates */
    points[0].x = (short) ( (x1 + x_offset) * x_scale + BORDER );
    points[0].y = (short) (((float)win_height)-(y1+y_offset)*y_scale -BORDERY);
    points[1].x = (short) ( (x2 + x_offset) * x_scale + BORDER );
    points[1].y = (short) (((float)win_height)-(y2+y_offset)*y_scale -BORDERY);
  
    /* Clip lines that fall outside plotting region */
    if (points[0].x < (BORDER - 10))    continue;
    if (points[0].y < (BORDERY - 10))   continue;
    
    XSetForeground(display,back_gc,(BlackPixel(display,screen_num)));
    XDrawLine(display,back_pixmap,back_gc,points[0].x,points[0].y,
	      points[1].x,points[1].y);
  }
}

/*============================= load_font =======================*/
/* Called by: open_win$. */
void load_font(XFontStruct **font_info)
{
  char *fontname = "6x13";
  
  /* Access font */
  if ((*font_info = XLoadQueryFont(display,fontname)) == NULL)
    {
      (void) fprintf( stderr, "Basic: Cannot open 6x13 font\\n");
      exit( -1 );
    }
}
/*============================ get_GC ==========================*/
/* Called by setwindow$,open_win$, make_boxes$*/
void get_GC(GC *gc, XFontStruct *font_info)
{
  unsigned int line_width = 1;
  int dash_offset = 0;
  static char dash_list[] = {
    12, 24	};
  int list_length = 2;
  unsigned long valuemask = GCGraphicsExposures; 
  XGCValues values;
  int line_style = LineSolid;
  int cap_style = CapRound;
  int join_style = JoinRound;
  
  values.graphics_exposures = False;
  
  /* Create default Graphics Context */
  *gc = XCreateGC(display, win, valuemask, &values);
  
  /* specify font */
  XSetFont(display, *gc, font_info->fid);
  
  /* specify black foreground since default may be white on white */
  XSetForeground(display, *gc, foreground_pixel);
  
  /* set line attributes */
  XSetLineAttributes(display, *gc, line_width, line_style, cap_style, 
		     join_style);
  
  /* set dashes to be line_width in length */
  XSetDashes(display, *gc, dash_offset, dash_list, list_length);
}
/*============================ get_xorGC =========================*/
/* Called by: open_win$.*/
void get_xorGC(GC *gc,XFontStruct *font_info)
{
  unsigned int line_width = 1;
  int dash_offset = 0;
  static char dash_list[] = {
    12, 24	};
  int list_length = 2;
  unsigned int fore,back;
  XGCValues values;
  unsigned long valuemask = GCFunction | GCGraphicsExposures; 
  int line_style = LineSolid;
  int cap_style = CapRound;
  int join_style = JoinRound;
  
  fore = BlackPixel(display,screen_num);
  back = WhitePixel(display,screen_num);
  
  valuemask = valuemask | GCBackground | GCForeground;
  
  values.graphics_exposures = False;
  values.function           = GXinvert;
  values.foreground         = fore ^ back;
  values.background         = 0;
  
  /* Create default Graphics Context */
  *gc = XCreateGC(display, win, valuemask, &values);
  
  /* specify font */
  XSetFont(display, *gc, font_info->fid);
  
  /* set line attributes */
  XSetLineAttributes(display, *gc, line_width, line_style, cap_style, 
		     join_style);
  
  /* set dashes to be line_width in length */
  XSetDashes(display, *gc, dash_offset, dash_list, list_length);
}
/*======================= make_boxes() =======================*/
/* Called by: open_win$. */
void make_boxes(int max_color,int do_color)
{
  int i,x;
  unsigned int depth;
  
  depth = DefaultDepth(display,screen_num);
  
  for (i=0; i<max_color; i++) {
    pixmaps[i] = XCreatePixmap(display,win,15,15,depth);
    get_GC(&(pix_gc[i]),font_info);
    if (do_color) {
      XSetForeground(display,pix_gc[i],c_map[i]);
    }
    else {
      if (i == 0) {
	XSetForeground(display,pix_gc[i],(WhitePixel(display,screen_num)));
      }
      else {
	XSetForeground(display,pix_gc[i],(BlackPixel(display,screen_num)));
      }
    }
    XFillRectangle(display,pixmaps[i],pix_gc[0],0,0,15,15); 
    XFillArc(display,pixmaps[i],pix_gc[i],0,0,6,6,0,360*64);
  }

  for (i=0; i<128; i++) {
    x = i*10;
    XCopyArea(display,pixmaps[i],win,pix_gc[i],0,0,10,10,x,100);
  }
  XFlush(display);
}
/*===================== get_traj_GC =======================*/
/* Called by: setwindow$ */
void get_traj_GC(GC *gc,XFontStruct *font_info)
{
  unsigned long valuemask = GCGraphicsExposures; 
  XGCValues values;
  unsigned int line_width = 1;
  int line_style = LineSolid;
  int cap_style = CapRound;
  int join_style = JoinRound;
  int dash_offset = 0;
  static char dash_list[] = {
    12, 24	};
  int list_length = 2;
  
  values.graphics_exposures = False;
  
  /* Create default Graphics Context */
  *gc = XCreateGC(display, win, valuemask, &values);
  
  /* specify font */
  XSetFont(display, *gc, font_info->fid);
  
  /* specify black foreground since default may be white on white */
  XSetForeground(display, *gc, foreground_pixel);
  
  /* set line attributes */
  XSetLineAttributes(display, *gc, line_width, line_style, cap_style, 
		     join_style);
  
  /* set dashes to be line_width in length */
  XSetDashes(display, *gc, dash_offset, dash_list, list_length);
}
/*======================= BusyCursor ====================*/
/* switch between default cursor and BUSY one */
void BusyCursor()
{
  static int busy = 0;

  if (!(busy)) {
    cursor = XCreateFontCursor(display,XC_gumby);
    if ( cursor != (Cursor)None) {   
      XDefineCursor(display,win,cursor);
      busy = 1;
    }
  }
  else {
    if ( cursor != (Cursor)None) {
      XUndefineCursor(display,win);
      XFreeCursor(display,cursor);
      busy = 0;
    }
  }
}	 
/*========================= check_events() ======================*/
/* Called by: main$*/
void check_events(int *key_press,int *do_rewind,int *do_fast_forward,
		  int num_pins,int argc,char **argv,int max_color,
		  int *do_color,int monochrome,float uxmin,float uxmax,
		  float uymin,float uymax,int frame_num,int do_S_contour,
		  int do_T_contour,
		  float PinSize,int *do_clear,int traj_on,int *paused,
		  float *orig_xmin,float *orig_xmax,float *orig_ymin,
		  float *orig_ymax,int do_net_contour,
		  float centerx,float centery,int do_xshift,int do_yshift,
		  float xshift,float yshift,float xmagnify,float ymagnify,
		  float syssizex,float syssizey)
{
  int button;
  XEvent report;

  *key_press = 0;
  while (XCheckWindowEvent(display,win,(EVENT_MASK),&report)) {
    switch  (report.type) {
    case Expose:
      /* get all other Expose events on the queue */
      while (XCheckTypedEvent(display, Expose, &report));
      break;
    case ConfigureNotify:
      /* window has been resized, change width and
       * height to send to place_text and place_graphics
       * in next Expose */

      setwindow(uxmin,uymin,uxmax,uymax,num_pins,argc,argv,
		max_color,do_color,monochrome,do_S_contour,do_T_contour,
		PinSize,
		do_clear,traj_on,orig_xmin,orig_xmax,orig_ymin,orig_ymax,
		do_net_contour,centerx,centery,do_xshift,do_yshift,
		xshift,yshift,xmagnify,ymagnify,syssizex,syssizey,0); 
      break;
    case KeyPress:
      *key_press = 1;
      break;        
    case ButtonPress:
      button = report.xbutton.button;

      if(button ==1){
	if(*do_fast_forward){
	  *do_fast_forward = 0;
	  break;
	}
	else if(frame_num > 16)
	  *do_rewind = 1;
      }

      if (button == 2) {
        if (*do_rewind) {
          *do_rewind = 0;
          break;
        }
        if (*do_fast_forward) {
          *do_fast_forward = 0;
          break; 
        }  
	if (!(*paused)) {
	  *paused = 1;
	}
      }

      if (button == 3) {
	if(*do_rewind){
	  *do_rewind = 0;
	  break;
	}
	else
	  *do_fast_forward = 1;
      }  
    
      break;
    default:
      /* all events selected by StructureNotifyMask
       * except ConfigureNotify are thrown away here,
       * since nothing is done with them */
      break;
    } /* end switch */
  }
}
/*============================= pause_events =====================*/
/* Called by: main$.*/
void pause_events(int num_pins,int argc,char **argv,int max_color,
		  int *do_color,int monochrome,float *uxmin,float *uxmax,
		  float *uymin,float *uymax,int do_S_contour,
		  int do_T_contour,float PinSize,
		  int *do_clear,int traj_on,int *paused,float *orig_xmin,
		  float *orig_xmax,float *orig_ymin,float *orig_ymax,
		  int do_net_contour,float centerx,float centery,
		  int do_xshift,int do_yshift,float xshift,float yshift,
		  float xmagnify,float ymagnify,float syssizex,
		  float syssizey)
{
  void get_box(int start_x,int start_y,int *x1,int *x2,int *y1,int *y2);

  int x,y;
  int button;
  int do_resize = 0;
  int x1,x2,y1,y2;
  static int resized = 0;
  XEvent report;

  while(*paused) {

    XNextEvent(display,&report);
    
    switch  (report.type) {
    case Expose:
      /* get all other Expose events on the queue */
      while (XCheckTypedEvent(display, Expose, &report));
      break;
      
    case ConfigureNotify:
      /* window has been resized, change width and
       * height to send to place_text and place_graphics
       * in next Expose */
      
      setwindow(*uxmin,*uymin,*uxmax,*uymax,num_pins,argc,argv,
		max_color,do_color,monochrome,do_S_contour,do_T_contour,
		PinSize,
		do_clear,traj_on,orig_xmin,orig_xmax,orig_ymin,orig_ymax,
		do_net_contour,centerx,centery,do_xshift,do_yshift,
		xshift,yshift,xmagnify,ymagnify,syssizex,syssizey,0);
      break;
      
    case KeyPress:
      break;
      
    case ButtonPress:
      x = (int)(report.xbutton.x);
      y = (int)(report.xbutton.y);
      button = report.xbutton.button;

      /* Zoom feature can be used when movie is paused */
      /* CIJOL DISABLED 10.8.97 since it doesn't work */
      if (button == 1) {
	/*do_resize = 1;
	get_box(x,y,&x1,&x2,&y1,&y2);*/
	printf("Zoom feature not supported\n");
      }
      
      /* Unpause movie and set zoom or unzoom if requested */
      if (button == 2) {
	*paused = 0;
	if (do_resize) {
	  *uxmin = ((float)(x1 - BORDER))/x_scale - x_offset;
	  *uxmax = ((float)(x2 - BORDER))/x_scale - x_offset;
	  *uymax = ((float)(win_height - y1 - BORDER))/y_scale - y_offset;
	  *uymin = ((float)(win_height - y2 - BORDER))/y_scale - y_offset;
	  setwindow(*uxmin,*uymin,*uxmax,*uymax,num_pins,
		    argc,argv,max_color,do_color,monochrome,do_S_contour,
		    do_T_contour,
		    PinSize,do_clear,traj_on,orig_xmin,orig_xmax,
		    orig_ymin,orig_ymax,do_net_contour,centerx,centery,
		    do_xshift,do_yshift,xshift,yshift,xmagnify,ymagnify,
		    syssizex,syssizey,0);
	  resized = 1;
	}
      }
      
      /* Unzoom (back to original size). */
      /* CIJOL DISABLED 10.8.97 since it doesn't work */
      if (button == 3) {
	/*if ( (!(do_resize)) && (resized) ) {
	  *uxmin = *orig_xmin;
	  *uymin = *orig_ymin;
	  *uxmax = *orig_xmax;
	  *uymax = *orig_ymax;

	  setwindow(*uxmin,*uymin,*uxmax,*uymax,num_pins,
		    argc,argv,max_color,do_color,monochrome,do_S_contour,
		    do_T_contour,
		    PinSize,do_clear,traj_on,orig_xmin,orig_xmax,
		    orig_ymin,orig_ymax,do_net_contour,centerx,centery,
		    do_xshift,do_yshift,xshift,yshift,xmagnify,ymagnify,
		    syssizex,syssizey,0);
	  resized = 0;
	}
	do_resize = 0;*/
	printf("Unzoom feature not supported\n");
      }
      
      break;
      
    case ButtonRelease:
      x = report.xbutton.x;
      y = report.xbutton.y;
      button = report.xbutton.button;
      
      if (button == 1) {
	printf("Button one was RELEASED!\n");
      }
      
      
      break;
      
    default:
      /* all events selected by StructureNotifyMask
       * except ConfigureNotify are thrown away here,
       * since nothing is done with them */
      break;
    } /* end switch */
    
  }

}

/*============================ set_color_mode ====================*/
/* Called by: set_comm$.*/
void set_color_mode(char *line,int pos,int *mode256,int *max_color)
{
  char tempstring[10];
  int num;
  
  if (get_arg(line,&pos,tempstring)) {
    num = atof(tempstring);
    if ( (num == 256) || (num == 128) ) {
      if ((num == 256) && (*mode256 == 0)) {
	*mode256 = 1;
	*max_color = 256;
      }
      if ((num == 128) && (*mode256 == 1)) {
	*mode256 = 0;
	*max_color = 128;
      }
    }
    else printf("Illegal color mode value (128 or 256 only)\n");
  }
  else printf("need to know which color mode to set to\n");
}
/*==================== set_sample ============================*/
/* this function gets the size the user expects the file to be */
/* This is done so the required memory can be grabbed all at   */
/* once, in order to prevent fragmentation.                    */
/* That is to say, this function can be used to increase the */
/* total allowed number of particles plotted at any given time. */
/*Called by set_comm$.*/
void set_sample(char *line,int pos,int *samples,float *x_pts,float *y_pts,
		int *c_pts,int *p_pts)
{
  char samplestring[10];
  int old_samples,error;
  
  /* Keep track of the former maximum value */
  old_samples = *samples;

  /* Get desired new value */
  if (get_arg(line,&pos,samplestring)) {
    *samples = atoi(samplestring);
    /* Attempt to allocate enough memory */
    error = get_sample_mem(*samples,x_pts,y_pts,c_pts,p_pts);
    if (error) {
      printf("Couldn't allocate enough memory\n");
      /*If we fail, go back to former amount of memory */
      *samples = old_samples;
      error = get_sample_mem(*samples,x_pts,y_pts,c_pts,p_pts);
      /* If this fails, serious problem. */
      if (error) {
	printf("Memory shortage critical!\nExiting program!\n");
	exit(100);
      }
    }
  }
  else
    printf("Need to know # of samples you wish\n");
  
}
/*======================= set_xrange ===================*/
/* sets the range of plotting and replotting functions */
/* Called by: set_comm$.*/
void set_xrange(char *line,int pos,float *uxmin,float *uxmax)
{
  float a,b;
  char tempstring[10];
    
  if (get_arg(line,&pos,tempstring)) 
    a = atof(tempstring);
  else 
    return;

  if (get_arg(line,&pos,tempstring)) 
    b = atof(tempstring);
  else {
    printf("Need BOTH the max and min values!\n");
    return;
  }
  if (a<b) {
    *uxmin = a;
    *uxmax = b;
  }
  else printf("ERROR: xmin is BIGGER then xmax!\n");
}
/*======================== set_yrange ===========================*/
/* sets yrange */
/* Called by: set_comm$.*/
void set_yrange(char *line,int pos,float *uymin,float *uymax)
{
  float a,b;
  char tempstring[10];
  
  if (get_arg(line,&pos,tempstring)) 
    a = atof(tempstring);
  else 
    return;

  if (get_arg(line,&pos,tempstring)) 
    b = atof(tempstring);
  else {
    printf("Need BOTH the max and min values!\n");
    return;
  }
  if (a<b) {
    *uymin = a;
    *uymax = b;
  }
  else printf("ERROR: xmin is BIGGER then xmax!\n");
}
/*=========================== set_delay ====================*/
/* Called by: set_comm$.*/
void set_delay(char *line,int pos,unsigned int *delay_counter)
{ 
  char tempstring[10];
  
  if (get_arg(line,&pos,tempstring))
    *delay_counter = 1000 * atof(tempstring);
  else
    *delay_counter = 200*1000;
}
/*=========================== load_color_map ===================*/
/* this routine reads in a file containing data in the   */
/* format R G B\n... NOTE: there is no checking yet      */
/* to prevent the pen # getting to big.                 */
/* Called by: get_command$.*/
void load_color_map(char *line,int pos,int *max_color,int do_color)
{
  char filename[200];
  float red,green,blue;
  FILE *mapfile;
  int i,error;
  XColor in_out;
  
  if (!(do_color)) return;
  i = 0;

  if (get_arg(line,&pos,filename)) {
    if ( (mapfile = fopen(filename,"r")) == NULL) {
      printf("Unknown file!\n");
      return;
    }
    while (fscanf(mapfile,"%f %f %f",&red,&green,&blue) != EOF) {
      
      in_out.green = green * 65535;
      in_out.red   = red * 65535;
      in_out.blue  = blue * 65535;
      error = XAllocColor(display,DefaultColormap(display,screen_num),&in_out);
      if (!(error)) break;
      c_map[i] = in_out.pixel;
      i++;
    }
  }
  else
    printf("Need a file name!\n");

  *max_color = i-1;
}
/*=============================== load_special_contour ===========*/
/* Reads in a special_contour file (ASCII). */
/* Called by get_command$ */
void load_special_contour(char *line,int pos,int *num_pins,int *do_S_contour,
			  int *do_T_contour,
			  float *PinSize)
{
  char filename[200];
  FILE *contour;
  int j;

  if(get_arg(line,&pos,filename)){
    if((contour=fopen(filename,"r"))==NULL){
      printf("Contour file %s not found\n",filename);
      *do_S_contour = 0;
      return;
    }

    *do_S_contour = 1;
    fscanf(contour,"%d",num_pins);

    pin_sites = (struct pinning_site *)malloc(sizeof(struct pinning_site)*(*num_pins+10));
    if (pin_sites == NULL){
      printf("Could not allocate pinning site memory\n");
      exit(10);
    }
    for(j=0;j<*num_pins;j++){
      fscanf(contour,"%f",&(pin_sites[j].x));
      fscanf(contour,"%f",&(pin_sites[j].y));
      fscanf(contour,"%f",&(pin_sites[j].Ap));
      //fscanf(contour,"%f",&(pin_sites[j].depth));
      fscanf(contour,"%f",&(pin_sites[j].radiusx));
      fscanf(contour,"%f",&(pin_sites[j].radiusy));
      *PinSize = pin_sites[j].radiusx*2.0;
    }	     
    if(*PinSize < 0.01){
	printf("at pin %d radius is %f \n",j,pin_sites[j].radiusx);
      printf("Please enter pinning diameter: ");
      scanf("%f",PinSize);
    }
  }
  else
    printf("Need a file name\n");
}

/*=============================== load_triangle_contour ===========*/
/* Reads in a contour file containing locations of the vertices for */
/* triangular pinning sites. */
/* Called by get_command$ */
void load_triangle_contour(char *line,int pos,int *num_pins,int *do_S_contour,
			  int *do_T_contour,
			  float *PinSize)
{
  char filename[200];
  FILE *contour;
  int j;

  if(get_arg(line,&pos,filename)){
    if((contour=fopen(filename,"r"))==NULL){
      printf("Triangle contour file %s not found\n",filename);
      *do_T_contour = 0;
      return;
    }

    *do_T_contour = 1;
    fscanf(contour,"%d",num_pins);
    printf("Read %d pins\n",*num_pins);

    pin_sites = (struct pinning_site *)malloc(sizeof(struct pinning_site)*(*num_pins+10));
    if (pin_sites == NULL){
      printf("Could not allocate pinning site memory\n");
      exit(10);
    }
    for(j=0;j<*num_pins;j++){
      /* Input is in the form of triples of points */
      fscanf(contour,"%f %f %f %f %f %f",&(pin_sites[j].px[0]),
	     &(pin_sites[j].py[0]),&(pin_sites[j].px[1]),
	     &(pin_sites[j].py[1]),&(pin_sites[j].px[2]),
	     &(pin_sites[j].py[2]));
    }
  }
  else
    printf("Need a file name\n");
}


/*=============================== load_network_contour ===========*/
/* Reads in a contour file containing locations of the vertices for */
/* long-range network connections. */
/* Called by get_command$ */
void load_network_contour(char *line,int pos,int *num_pins,int *do_S_contour,
			  int *do_net_contour,float *PinSize,
			  float *centerx,float *centery)
{
  char filename[200];
  FILE *contour;
  int j;

  if(get_arg(line,&pos,filename)){
    if((contour=fopen(filename,"r"))==NULL){
      printf("Network contour file %s not found\n",filename);
      *do_net_contour = 0;
      return;
    }

    *do_net_contour = 1;
    fscanf(contour,"%d %f %f",num_pins,centerx,centery);
    printf("Read %d links\n",*num_pins);

    pin_sites = (struct pinning_site *)malloc(sizeof(struct pinning_site)*(*num_pins+10));
    if (pin_sites == NULL){
      printf("Could not allocate pinning site memory\n");
      exit(10);
    }
    for(j=0;j<*num_pins;j++){
      /* Input is in the form of pairs of points */
      /* which are separated so they can be plotted by xmgrace */
      fscanf(contour,"%f %f",&(pin_sites[j].px[0]),&(pin_sites[j].py[0]));
      fscanf(contour,"%f %f",&(pin_sites[j].px[1]),&(pin_sites[j].py[1]));
    }
  }
  else
    printf("Need a file name\n");
}

/*=============================== load_3D_contour ===========*/
/* Reads in a 3D pinfile (ASCII). */
/* Called by get_command$ */
void load_3D_contour(char *line,int pos,int *num_pins,int *do_S_contour,
		     int *do_T_contour,
		     float *PinSize)
{
  char filename[200];
  FILE *contour;
  int j;
  int layer,num_layers;
  float another_number;

  if(get_arg(line,&pos,filename)){
    if((contour=fopen(filename,"r"))==NULL){
      printf("Pin file %s not found\n",filename);
      *do_S_contour = 0;
      return;
    }

    *do_S_contour = 1;

    /* For michigan style only */
    fscanf(contour,"%d",&num_layers);

    /* Both styles */
    fscanf(contour,"%d",num_pins);

    printf("Reading in only the first layer\n");
    pin_sites = (struct pinning_site *)malloc(sizeof(struct pinning_site)*(*num_pins+10));
    if (pin_sites == NULL){
      printf("Could not allocate pinning site memory\n");
      exit(10);
    }
    for(j=0;j<(*num_pins);j++){
      /* Michigan style */
      fscanf(contour,"%f",&(pin_sites[j].x));
      fscanf(contour,"%f",&(pin_sites[j].y));
      fscanf(contour,"%f",&(pin_sites[j].depth));
      fscanf(contour,"%f",&another_number);
      fscanf(contour,"%f",&another_number);
      /* van O. style*/
     /*fscanf(contour,"%d %f %f %f",&layer,&(pin_sites[j].x),&(pin_sites[j].y),
	&(pin_sites[j].depth));*/
    }	     
    printf("Please enter pinning radius: ");
    scanf("%f",PinSize);
  }
  else
    printf("Need a file name\n");
}
/*========================= get_box ====================*/
/* Called by pause_events$.*/
void get_box(int start_x,int start_y,int *x1,int *x2,int *y1,int *y2)
{
  int x,y;
  int button;
  int xmin,xmax,ymin,ymax,length,width;
  int old_xmin,old_ymin,old_length,old_width;
  int done = 0;
  XEvent report;

  old_length = old_width = 0;
  old_ymin = old_xmin = 0;
  while(!(done)) {
    XNextEvent(display,&report);

    switch  (report.type) {
    case ButtonRelease:
      x = report.xbutton.x;
      y = report.xbutton.y;
      button = report.xbutton.button;
      
      if (button == 1) {
	done = 1;
	xmin = MIN(start_x,x);
	xmax = MAX(start_x,x);
	ymin = MIN(start_y,y);
	ymax = MAX(start_y,y);
	length = xmax - xmin;
	width  = ymax - ymin;
	*x1 = xmin;
	*x2 = xmax;
	*y1 = ymin;
	*y2 = ymax;

	XSetForeground(display,gc,(WhitePixel(display,screen_num)));    
	XDrawRectangle(display,win,gc,xmin,ymin,length,width);
	XFlush(display);
      }
      break;

    case MotionNotify:
      x = report.xmotion.x;
      y = report.xmotion.y;
      xmin = MIN(start_x,x);
      xmax = MAX(start_x,x);
      ymin = MIN(start_y,y);
      ymax = MAX(start_y,y);
      length = xmax - xmin;
      width  = ymax - ymin;

      if ( (old_length) || (old_width) ) {
	XDrawRectangle(display,win,xor_gc,old_xmin,old_ymin,old_length,old_width);
      }

      old_xmin   = xmin;
      old_ymin   = ymin;
      old_length = length;
      old_width  = width;
      
      XDrawRectangle(display,win,xor_gc,xmin,ymin,length,width);
      break;		     
		     
    default:
      /* all events selected by StructureNotifyMask
     * except ConfigureNotify are thrown away here,
     * since nothing is done with them */
      break;
    } /* end switch */
  }
}

/*======================= print_help_screen ===================*/
/* Prints a help screen containing summary of available commands */
void print_help_screen()
{
  char trash[10];
  printf("Summary of commands:\n");
  printf(" set xrange <min> <max>\tSets x values for plotting region\n");
  printf(" set yrange <min> <max>\tSets y values for plotting region\n");
  printf(" set sample <#>\tSets max. number of particles. Default 5000.\n");
  printf(" set delay <#>\tDelay between drawing screens.  Default 0.\n");
  printf(" set colormode <128 | 256>\tLeave set at default value 128\n");
  printf(" loadcolormap <colormap>\tReads in colormap file\n");
  printf(" monochrome\tDisables use of color; for monochrome monitors\n");
  printf(" smovie\tMode for plotting smovies (default)\n");
  printf(" kmovie\tGo into kmovie (river format) plotting mode\n");
  printf(" cmovie\tGo into cmovie (colorized format) plotting mode\n");
  printf(" tmovie\tGo into tmovie (3D format) plotting mode\n");
  printf(" do_stripe\tStripe mode: draws lines indicating spin directions\n");
  printf(" loadcontour <contourfile>\tReads contour, enables contour plot\n");
  printf(" loadtriangle <contourfile>\tReads triangle contour\n");
  printf(" load3Dpins <contourfile>\tReads 3D contour, top layer only\n");
  printf(" loadnetwork <contourfile>\tReads network contour\n");
  printf(" toggle_traj\tToggles drawing of trajectories\n");
  printf(" toggle_voronoi\tToggles Voronoi polygon mode\n");
  printf(" delaunay\tToggles Delaunay triangulation mode\n");
  printf(" onelayer <#>\tPlot only specified layer number\n");
  printf(" alllayer\tPlot all layers (default)\n");
  printf(" onevortex <#>\tPlot only specified vortex\n");
  printf(" allvortex\tPlot all vortices (default)\n");
  printf(" xmagnify, ymagnify\tMultiply all x or y coordinates by a factor\n");
  printf(" xshift, yshift\t Shift system through periodic boundary\n");
  printf(" ntypes <#>\tSet number of particle types for overloaded tmovie\n");
  printf(" numnodes <#>\tSet number of nodes\n");
  printf(" execute <file>\tRead commands from file\n");
  printf(" plot <file>\tPlot 3D tm* file\n");
  printf(" replot\tPlot same file again\n");
  printf(" quit\tEnd program\n");
  printf(" clear\tClear screen; moving window with mouse has same effect\n");
  printf(" disable_rewind\tDisables unbreakable rewind (for large files)\n");
  printf(" cont\tContinue plotting movie from this point\n");
  printf(" help\tShow this screen\n");
  printf("KEYBOARD: Press key in movie window to get command prompt\n");
  printf("MOUSE COMMANDS: While movie is running\n");
  printf(" Left mouse button: Rewind, or stop fast forward\n");
  printf(" Middle mouse button: Pause, or stop rewind or fast forward\n");
  printf(" Right mouse button: Fast forward, or stop rewind\n");
  printf("MOUSE COMMANDS: While movie is paused\n");
  printf(" Left mouse button: Select zoom region (NOT SUPPORTED)\n");
  printf(" Middle mouse button: Un-pause\n");
  printf(" Right mouse button: Auto-size (un-zoom) (NOT SUPPORTED)\n");
}

/*====================== rewind_movie ===========================*/
/* There are two rewind methods.  The first is "unbreakable," */
/* and uses the data structure rewind_dat.  It is, however, */
/* limited to files that are not tremendously long.  The second */
/* can be used on files of any length, but will not work correctly */
/* in a sample from which vortices are removed.  The first method */
/* is the default; it can be disabled by hand or by lack of memory, */
/* in which case the second method is used. */
void rewind_movie(int *frame_num,int *do_rewind,int disable_rewind,
		  int num_nodes,int movie_type,int ntypes)
{
  int i,j;
  long int steps;
  struct smdata smtempdat;
  struct kmdata kmtempdat;
  struct cmdata cmtempdat;
  struct tmdata tmtempdat;
  int num_layers;
  int num_pars;

  if (*frame_num < 15) {
    /* We aren't far enough into the file to rewind.*/
    *do_rewind = 0;
    return;
  }
  /* The first rewind method (unbreakable) */
  if(!disable_rewind){
    /* Back up one position and find the number of layers. */
    /* Rewind all files this way */
    for(j=0;j<num_nodes;j++){
      switch(movie_type){
      case SMOVIE:
	steps = sizeof(smdata);
	fseek(picfile[j],-steps,1);
	fread(&smtempdat,sizeof(smdata),1,picfile[j]);
	num_layers = 1;
	/* Go backwards five frames.*/
	for(i=0;i<5;i++){
	  steps = 
	    (long)((rewind_dat[j][(*frame_num-1-i)])*sizeof(smdata) 
		   + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	}
	break;
      case KMOVIE:
	steps = sizeof(kmdata);
	fseek(picfile[j],-steps,1);
	fread(&kmtempdat,sizeof(kmdata),1,picfile[j]);
	num_layers = 1;
	/* Go backwards five frames.*/
	for(i=0;i<5;i++){
	  steps = 
	    (long)((rewind_dat[j][(*frame_num-1-i)])*sizeof(kmdata) 
		   + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	}
	break;
      case CMOVIE:
	steps = sizeof(cmdata);
	fseek(picfile[j],-steps,1);
	fread(&cmtempdat,sizeof(cmdata),1,picfile[j]);
	num_layers = 1;
	/* Go backwards five frames.*/
	for(i=0;i<5;i++){
	  steps = 
	    (long)((rewind_dat[j][(*frame_num-1-i)])*sizeof(cmdata) 
		   + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	}
	break;
      case TMOVIE:
	steps = sizeof(tmdata);
	fseek(picfile[j],-steps,1);
	fread(&tmtempdat,sizeof(tmdata),1,picfile[j]);
	num_layers = tmtempdat.layr + 1;
	/* Go backwards five frames.*/
	for(i=0;i<5;i++){
	  steps = 
            (long)(num_layers*(rewind_dat[j][(*frame_num-1-i)])*sizeof(tmdata) 
		   + 3*sizeof(int));
	  fseek(picfile[j],-steps,1);
	}
	break;
      }
    }
    *frame_num -= 5;
  }
  
  /* The second rewind method (breakable) */
  else{
    /* Go backwards five frames.*/
    /* Rewind all files this way */
    for(j=0;j<num_nodes;j++){
      for (i=0; i<5; i++) {
	
	/* Back up one position and find out how many particles there are, */
	/* based on the ID number.  (Far from foolproof method) */
	/* This will also give us the number of layers there are. */
	switch(movie_type){
	case SMOVIE:
	  steps = sizeof(smdata);
	  fseek(picfile[j],-steps,1);
	  fread(&smtempdat,sizeof(smdata),1,picfile[j]);
	  if(smtempdat.p_num<0)
	    num_pars = -smtempdat.p_num + 1;
	  else
	    num_pars = smtempdat.p_num + 1;
	  num_layers = 1;
	  steps = (long)((num_pars)*sizeof(smdata) + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	  break;
	case KMOVIE:
	  steps = sizeof(kmdata);
	  fseek(picfile[j],-steps,1);
	  fread(&kmtempdat,sizeof(kmdata),1,picfile[j]);
	  if(kmtempdat.p_num<0)
	    num_pars = -kmtempdat.p_num + 1;
	  else
	    num_pars = kmtempdat.p_num + 1;
	  num_layers = 1;
	  steps = (long)((num_pars)*sizeof(kmdata) + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	  break;
	case CMOVIE:
	  steps = sizeof(cmdata);
	  fseek(picfile[j],-steps,1);
	  fread(&cmtempdat,sizeof(cmdata),1,picfile[j]);
	  if(cmtempdat.p_num<0)
	    num_pars = -cmtempdat.p_num + 1;
	  else
	    num_pars = cmtempdat.p_num + 1;
	  num_layers = 1;
	  steps = (long)((num_pars)*sizeof(cmdata) + 2*sizeof(int));
	  fseek(picfile[j],-steps,1);
	  break;
	case TMOVIE:
	  steps = sizeof(tmdata);
	  fseek(picfile[j],-steps,1);
	  fread(&tmtempdat,sizeof(tmdata),1,picfile[j]);
	  if(tmtempdat.p_num<0)
	    num_pars = -tmtempdat.p_num + 1;
	  else
	    num_pars = tmtempdat.p_num + 1;
	  num_layers = tmtempdat.layr + 1;
	  /* Back over the layers of data and the three integers.*/
	  steps = (long)(num_layers*(num_pars)*sizeof(tmdata)+3*sizeof(int));
	  fseek(picfile[j],-steps,1);
	  break;
	}
      }
    }
    (*frame_num) -= 5;
  }
}

/*===================== Voronoi routines ========================*/
/* Includes memory leak fix 8.12.03 */
/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */
/*--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==--==*/
/*  Defined constants and global variables for Voronoi routines */
#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

struct  Freenode        {
struct  Freenode        *nextfree;
};
struct  Freelist        {
struct  Freenode        *head;
int                     nodesize;
};

float xmin, xmax, ymin, ymax, deltax, deltay;
float pxmin, pxmax, pymin, pymax;

struct Point    {
float x,y;
};
 
/* structure used both for sites and for vertices */
struct Site     {
struct  Point   coord;
int             sitenbr;
int             refcnt;
};
 
/* CIJOL new: making into an array 8.12.03 */
struct  Site    sites[MAX_OBJECT];
int             nsites;
int             siteidx;
int             sqrt_nsites;
int             nvertices;

/* CIJOL new: prepare an array of free Sites 8.12.03 */
struct Site sfl_avail[MAXNUM];
int next_free_sfl;

struct  Freelist sfl;
struct  Site    *bottomsite;
 
struct Edge     {
float           a,b,c;
struct  Site    *ep[2];
struct  Site    *reg[2];
int             edgenbr;
};
#define le 0
#define re 1
int nedges;
struct  Freelist efl;
/* CIJOL new: prepare an array of free Edges 8.12.03 */
struct Edge efl_avail[MAXNUM];
int next_free_efl;
 
int has_endpoint(),right_of();
struct Site *intersect();
float dist();
struct Point PQ_min();
struct Halfedge *PQextractmin();
struct Edge *bisect();
 
struct Halfedge {
struct Halfedge *ELleft, *ELright;
struct Edge     *ELedge;
int             ELrefcnt;
char            ELpm;
struct  Site    *vertex;
float           ystar;
struct  Halfedge *PQnext;
};

/* CIJOL new: prepare an array of free Halfedges 8.12.03 */
struct Halfedge hfl_avail[MAXNUM];
int next_free_hfl;
 
struct   Freelist       hfl;
struct  Halfedge *ELleftend, *ELrightend;
int     ELhashsize;
/* This structure needs to be freed when done: ELhash */
/* CIJOL turning into array 8.12.03 */
/*struct  Halfedge **ELhash;*/
struct Halfedge *ELhash[MAXNUM];
struct  Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
struct  Site *leftreg(), *rightreg();
 
int PQhashsize;
/* This structure needs to be freed when done: PQhash */
/* CIJOL turning into array 8.12.03 */
/*struct  Halfedge *PQhash;*/
struct Halfedge PQhash[MAXNUM];
struct  Halfedge *PQfind();
int PQcount;
int PQmin;
int PQempty();

/*============================= plot_Voronoi ======================*/
/* Plots Voronoi construction data (animation) */
void plot_Voronoi(int movie_type,int num_pars,int sidenum[],int voronoi_layer,
		  int flag,int max_color)
{
  void calculate_voronoi(int nVin,float xdatin[],float ydatin[],float linea[],
			 float lineb[],float linec[],int *numlines,float rx[],
			 float ry[],int *site_counter,float vertx[],
			 float verty[],int *num_vert,int cjoleft[],
			 int cjoright[],int *num_ep,int tripA[],int tripB[],
			 int tripC[],int *num_trip);
  void plot_poly_line(float x1,float y1,float x2,float y2);
  void plot_polygon(int index,int sidenum,int polylist[][20],float vertx[],
		    float verty[],float xcenter,float ycenter);
  void plot_triangle(float x1,float y1,float x2,float y2,float x3,float y3,
		     int sideA,int sideB,int sideC,float dmin,float dmax,
		     float dstore[][3],int,int);
  /*void plot_vertex();*/

  float xdat[MAX_OBJECT],ydat[MAX_OBJECT];
  float rx[MAX_OBJECT],ry[MAX_OBJECT];
  float vertx[MAXNUM],verty[MAXNUM];
  int cjoleft[MAXNUM],cjoright[MAXNUM];
  int tripA[MAXNUM],tripB[MAXNUM],tripC[MAXNUM];
  float linea[MAXNUM],lineb[MAXNUM],linec[MAXNUM];
  int numlines;
  int site_counter;
  int num_vert;
  int num_ep;
  int num_trip;
  int i,j;
  float x[3],y[3];
  float xwidth,ywidth;
  int polynum[20000],polylist[20000][20]; /* Used to reconstruct polygons.*/
  int index;
  float dx,dy,dist;
  /* Can't bear to calculate the distances twice...*/
  float dstore[MAXNUM][3];
  float dmin,dmax;
  float avgdist;
  int count;

  /* Zero vertex structures */
  for(i=0;i<20000;i++){
    vertx[i] = 0;
    verty[i] = 0;
    polynum[i] = 0;
  }

  site_counter = 0;
  num_vert = 0;
  num_ep = 0;
  num_trip = 0;

  /* Pull required data out of frame: x,y positions only. */
  for(i=0;i<num_pars;i++){
    switch(movie_type){
    case SMOVIE:
      xdat[i]=smframes[i].x;
      ydat[i]=smframes[i].y;
      break;
    case KMOVIE:
      xdat[i]=kmframes[i].x;
      ydat[i]=kmframes[i].y;
      break;
    case CMOVIE:
      xdat[i]=cmframes[i].x;
      ydat[i]=cmframes[i].y;
      break;
    case TMOVIE:
      xdat[i]=tmframes[voronoi_layer][i].x;
      ydat[i]=tmframes[voronoi_layer][i].y;
      break;
    }
  }

  /* Now pass this info to the att voronoi program. */
  calculate_voronoi(num_pars,xdat,ydat,linea,lineb,linec,&numlines,
		    rx,ry,&site_counter,vertx,verty,&num_vert,
		    cjoleft,cjoright,&num_ep,
		    tripA,tripB,tripC,&num_trip);

  /* Count the number of sides each vortex has */
  for(i=0;i<MAX_OBJECT;i++)
    sidenum[i] = 0;

  for(i=0;i<num_trip;i++){
    sidenum[tripA[i]]++;
    sidenum[tripB[i]]++;
    sidenum[tripC[i]]++;
  }

  /* At this point we have all the required data: the locations of the */
  /* sides of all the Voronoi triangles.  We next need to plot it to screen*/
  /* Here, I plot the Voronoi polyhedra.  It would also be possible to */
  /* plot the Delaunay triangulation. */
  /* Use the value of flag to determine what to plot.  flag is 1 for */
  /* Voronoi polyhedra and 2 for Delaunay triangles. */
  switch(flag){
  case 1:
    /* We reconstruct the polygons based on data received from the program. */
    for(i=0;i<num_trip;i++){
      index = tripA[i];
      polylist[index][polynum[index]] = i;
      polynum[index]++;
      index = tripB[i];
      polylist[index][polynum[index]] = i;
      polynum[index]++;
      index = tripC[i];
      polylist[index][polynum[index]] = i;
      polynum[index]++;
    }

    /* Now plot these polygons directly, filled or unfilled according to  */
    /* the number of sides. */
    xwidth = (xmax - xmin)/2;
    ywidth = (ymax - ymin)/2;
    for(i=0;i<num_trip;i++){
      plot_polygon(i,polynum[i],polylist,vertx,verty,rx[i],ry[i]);
    }

    /* Also draw lines around the outsides of the polygons. */
    for(i=0;i<num_ep;i++){
      x[0] = vertx[cjoleft[i]];
      x[1] = vertx[cjoright[i]];
      y[0] = verty[cjoleft[i]];
      y[1] = verty[cjoright[i]];
      /* Clip lines that are longer than 1/2 the system size in x or y */
      if(fabs(x[0]-x[1])>xwidth) continue;
      if(fabs(y[0]-y[1])>ywidth) continue;
      plot_poly_line(x[0],y[0],x[1],y[1]);
    }
    break;
  case 2:
    /* Code to be used to plot the Delaunay triangulation */

    /* New: Coloring the lines according to length in order to */
    /* resolve the stresses in the system.  First I need to */
    /* identify the longest and shortest lengths. */
    /* I also store the distances since I can't bring myself to */
    /* calculate them twice. */
    xwidth = (xmax - xmin)/2;
    ywidth = (ymax - ymin)/2;
    dmin=10000.0;
    dmax=0.0;
    avgdist=0;
    count=0;
    for(i=0;i<num_trip;i++){
      x[0] = rx[tripA[i]];
      x[1] = rx[tripB[i]];
      x[2] = rx[tripC[i]];
      y[0] = ry[tripA[i]];
      y[1] = ry[tripB[i]];
      y[2] = ry[tripC[i]];
      dx=x[0]-x[1];
      dy=y[0]-y[1];
      if((fabs(x[0]-x[1])>xwidth)||(fabs(y[0]-y[1])>ywidth)){
	/* Skip this line, it's an artifact. */
	dstore[i][0]=-1;
      }
      else{
	dist=(float)sqrt((double)(dx*dx+dy*dy));
	dstore[i][0]=dist;
	if(dist<dmin) dmin=dist;
	if(dist>dmax) dmax=dist;
	avgdist += dist;
	count++;
      }
      dx=x[1]-x[2];
      dy=y[1]-y[2];
      if((fabs(x[1]-x[2])>xwidth)||(fabs(y[1]-y[2])>ywidth)){
	/* Skip this line, it's an artifact. */
	dstore[i][1]=-1;
      }
      else{
	dist=(float)sqrt((double)(dx*dx+dy*dy));
	dstore[i][1]=dist;
	if(dist<dmin) dmin=dist;
	if(dist>dmax) dmax=dist;
	avgdist += dist;
	count++;
      }
      dx=x[2]-x[0];
      dy=y[2]-y[0];
      if((fabs(x[2]-x[0])>xwidth)||(fabs(y[2]-y[0])>ywidth)){
	/* Skip this line, it's an artifact. */
	dstore[i][2]=-1;
      }
      else{
	dist=(float)sqrt((double)(dx*dx+dy*dy));
	dstore[i][2]=dist;
	if(dist<dmin) dmin=dist;
	if(dist>dmax) dmax=dist;
	avgdist += dist;
	count++;
      }
    }
    avgdist /= (float)count;
    for(i=0;i<num_trip;i++){
      x[0] = rx[tripA[i]];
      x[1] = rx[tripB[i]];
      x[2] = rx[tripC[i]];
      y[0] = ry[tripA[i]];
      y[1] = ry[tripB[i]];
      y[2] = ry[tripC[i]];
      if((dstore[i][0]<0)||(dstore[i][1]<0)||(dstore[i][2]<0)) continue;
      /* Using 1.5*avgdist as an upper distance limit to avoid problems */
      /* with edge triangles. */
      if((dstore[i][0]>1.5*avgdist)||(dstore[i][1]>1.5*avgdist)
	 ||(dstore[i][2]>1.5*avgdist)) continue;
      plot_triangle(x[0],y[0],x[1],y[1],x[2],y[2],sidenum[tripA[i]],
		    sidenum[tripB[i]],sidenum[tripC[i]],dmin,avgdist*1.5,
		    dstore,i,max_color);
    }
    break;
  }
}

/*======================= calculate_voronoi ============*/
/* CIJOL making this a routine (was main()): */
/* nVin:  Number of particles */
/* xdatin[], ydatin[]: X,Y positions of particles */
/* linea[], lineb[], linec[]:  Equations of lines */
/* numlines: Number of lines */
/* vertx[], verty[]: Locations of vertices of Voronoi polygons */
/* num_vert: Number of vertices */
/* cjoleft[], cjoright[]: Used to index left and right edges of each line */
/*                        segment in the Voronoi construction */
/* num_ep: Number of these line segments */
/* tripA[], tripB[], tripC[]: The vortex positions forming a triangle around */
/*                            each of the Voronoi polygon vertex positions. */
/* num_trip: The number of these triangles, probably equals num_vert? */

void calculate_voronoi(int nVin,float xdatin[],float ydatin[],float linea[],
		       float lineb[],float linec[],int *numlines,float rx[],
		       float ry[],int *site_counter,float vertx[],
		       float verty[],int *num_vert,int cjoleft[],
		       int cjoright[],int *num_ep,int tripA[],int tripB[],
		       int tripC[],int *num_trip)
{
  void read_program_sites(int *num_part,int nvin,float xdatin[],
			  float ydatin[]);
  void voronoi(struct Site *(*nextsite)(),float rx[],float ry[],
	       int *site_counter,float vertx[],float verty[],int *num_vert,
	       float linea[],float lineb[],float linec[],int *numlines,
	       int cjoleft[],int cjoright[],int *num_ep,int tripA[],
	       int tripB[],int tripC[],int *num_trip);
  void geominit();
  struct Site *nextone();
  struct Site *(*next)();
  void read_program_sites();
 
  /*CIJOL added:*/
  int num_part;
 
  *numlines = 0;

  freeinit(&sfl, sizeof *sites);
 
  /* CIJOL: Now receiving data directly from calling program.*/
  read_program_sites(&num_part,nVin,xdatin,ydatin);
  next = nextone;
 
  /* We have now read in all of the data on particle positions. */
  siteidx = 0;
  geominit();
 
  voronoi(next,rx,ry,site_counter,vertx,verty,num_vert,
          linea,lineb,linec,numlines,cjoleft,cjoright,num_ep,
          tripA,tripB,tripC,num_trip);
 
  return;
}
 
/*=============== scomp ===================*/
/* sort sites on y, then x, coord */
int scomp(struct Point *s1,struct Point *s2)
{
  if(s1 -> y < s2 -> y) return(-1);
  if(s1 -> y > s2 -> y) return(1);
  if(s1 -> x < s2 -> x) return(-1);
  if(s1 -> x > s2 -> x) return(1);
  return(0);
}

/*============= nextone =============*/
/* return a single in-storage site */
struct Site *nextone()
{
  struct Site *s;
  if(siteidx < nsites){
    s = &sites[siteidx];
    siteidx += 1;
    return(s);
  }
  else
    return( (struct Site *)NULL);
}
/*================== read_program_sites ===================*/
/* read all sites, sort, and compute xmin, xmax, ymin, ymax */
/* Assumes that it has already been passed the needed data.*/
/* It just funnels this data into the correct structures.*/
void read_program_sites(int *num_part,int nvin,float xdatin[],float ydatin[])
{
  int i;
 
  nsites=0;
  /* CIJOL altering sites to be a fixed-size array, 8.12.03 */
  /*sites = (struct Site *) myalloc(4000*sizeof *sites);*/
  
  *num_part = nvin;
  nsites = nvin;
  for(i=0;i<*num_part;i++){
    sites[i].coord.x = xdatin[i];
    sites[i].coord.y = ydatin[i];
    sites[i].sitenbr = i;
    sites[i].refcnt = 0;
  };

  qsort(sites, nsites, sizeof *sites, scomp);
  xmin=sites[0].coord.x;
  xmax=sites[0].coord.x;
  for(i=1; i<nsites; i+=1){
    if(sites[i].coord.x < xmin) xmin = sites[i].coord.x;
    if(sites[i].coord.x > xmax) xmax = sites[i].coord.x;
  }
  ymin = sites[0].coord.y;
  ymax = sites[nsites-1].coord.y;
}

/*======================= voronoi ===============*/
/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
   deltax, deltay (can all be estimates).
   Performance suffers if they are wrong; better to make nsites,
   deltax, and deltay too big than too small.  (?) */
 
void voronoi(struct Site *(*nextsite)(),float rx[],float ry[],
	     int *site_counter,float vertx[],float verty[],int *num_vert,
	     float linea[],float lineb[],float linec[],int *numlines,
	     int cjoleft[],int cjoright[],int *num_ep,int tripA[],
	     int tripB[],int tripC[],int *num_trip)
{
  int out_ep(struct Edge *e,int cjoleft[],int cjoright[],int *num_ep);
  int out_site(struct Site *s,float rx[],float ry[],int *site_counter);
  int out_triple(struct Site *s1,struct Site *s2,struct Site *s3,int tripA[],
		 int tripB[],int tripC[],int *num_trip);
  int ELinsert(struct Halfedge *lb,struct Halfedge *new);
  int ELdelete(struct Halfedge *he);
  int endpoint(struct Edge *e,int lr,struct Site *s,int cjoleft[],
	       int cjoright[],int *num_ep);
  int PQinsert(struct Halfedge *he,struct Site *v,float offset);

  struct Site *newsite, *bot, *top, *temp, *p;
  struct Site *v;
  struct Point newintstar;
  int pm;
  struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
  struct Edge *e;

  /* CIJOL new 8.12.03 */
  next_free_hfl=0;
  next_free_efl=0;
  next_free_sfl=0;
 
  PQinitialize();
  bottomsite = (*nextsite)();
 
  /* Write out the coordinates of this "bottom" site:*/
  out_site(bottomsite,rx,ry,site_counter);
  
  ELinitialize();
  
  newsite = (*nextsite)();
  while(1)
    {
      if(!PQempty()) newintstar = PQ_min();
      
      if (newsite != (struct Site *)NULL
	  && (PQempty()
              || newsite -> coord.y < newintstar.y
              || (newsite->coord.y == newintstar.y
                  && newsite->coord.x < newintstar.x)))
        {/* new site is smallest */
	  
          /* Write out the coordinates of this "new" site: */
          out_site(newsite,rx,ry,site_counter);
 
          lbnd = ELleftbnd(&(newsite->coord));
	  rbnd = ELright(lbnd);
          bot = rightreg(lbnd);
          e = bisect(bot, newsite,linea,lineb,linec,numlines);
          bisector = HEcreate(e, le);
          ELinsert(lbnd, bisector);
          if ((p = intersect(lbnd, bisector)) != (struct Site *) NULL)
            {   PQdelete(lbnd);
                PQinsert(lbnd, p, dist(p,newsite));
              };
          lbnd = bisector;
          bisector = HEcreate(e, re);
          ELinsert(lbnd, bisector);
          if ((p = intersect(bisector, rbnd)) != (struct Site *) NULL)
            {   PQinsert(bisector, p, dist(p,newsite));
              };
          newsite = (*nextsite)();
        }
      else if (!PQempty())      /* intersection is smallest */
        {
          lbnd = PQextractmin();
	  llbnd = ELleft(lbnd);
          rbnd = ELright(lbnd);
          rrbnd = ELright(rbnd);
          bot = leftreg(lbnd);
          top = rightreg(rbnd);
 
          /* Write out this triple of numbers: bot, top, and rightreg.*/
          out_triple(bot, top, rightreg(lbnd),tripA,tripB,tripC,num_trip);
 
          v = lbnd->vertex;

          /* Create a vertex; also write it out.*/
          makevertex(v,vertx,verty,num_vert);
 
          /* Write out the two endpoints of this segment? */
          endpoint(lbnd->ELedge,lbnd->ELpm,v,cjoleft,cjoright,num_ep);
          endpoint(rbnd->ELedge,rbnd->ELpm,v,cjoleft,cjoright,num_ep);
 
          ELdelete(lbnd);
          PQdelete(rbnd);
          ELdelete(rbnd);
          pm = le;
          if (bot->coord.y > top->coord.y)
            {   temp = bot; bot = top; top = temp; pm = re;}
 
          /* This will write out bisector location */
          e = bisect(bot, top,linea,lineb,linec,numlines);
 
          bisector = HEcreate(e, pm);
          ELinsert(llbnd, bisector);
 
          /*This will write out endpoints*/
          endpoint(e, re-pm, v,cjoleft,cjoright,num_ep);
 
          deref(v);
          if((p = intersect(llbnd, bisector)) != (struct Site *) NULL)
            {   PQdelete(llbnd);
                PQinsert(llbnd, p, dist(p,bot));
              };
          if ((p = intersect(bisector, rrbnd)) != (struct Site *) NULL)
            {   PQinsert(bisector, p, dist(p,bot));
              };
        }
      else break;
    };
 
  for(lbnd=ELright(ELleftend); lbnd != ELrightend; lbnd=ELright(lbnd))
    {   e = lbnd -> ELedge;
        /* Write out endpoints */
        out_ep(e,cjoleft,cjoright,num_ep);
      };
}
/*================= out_bisector ==========*/
/* CIJOL Adding a pass of the lines data back to calling program.*/
out_bisector(struct Edge *e,float linea[],float lineb[],float linec[],
	     int *numlines)
{
  int num;
 
  num = *numlines;
 
  /*fprintf(out,"l %f %f %f", e->a, e->b, e->c);*/
  /*Put this data into structure.*/
  linea[num] = e->a;
  lineb[num] = e->b;
  linec[num] = e->c;
  num++;
  if(num>10000){
    printf("Not enough room for lines\n");
    exit(-1);
  }
  *numlines = num;
 
}
 
/*================ out_ep =============*/
/* Writes out an end point.*/
/* Each edge is numbered, and the other two numbers give the */
/* vertices at either end of this edge.*/
out_ep(struct Edge *e,int cjoleft[],int cjoright[],int *num_ep)
{
  int clip_line(struct Edge *e);
  int tmp;
 
  clip_line(e);
 
  /*fprintf(out,"e %d", e->edgenbr);
  fprintf(out," %d ", e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1);
  fprintf(out,"%d\n", e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1);*/
 
  tmp = *num_ep;
  cjoleft[e->edgenbr]=
   (e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1);
  cjoright[e->edgenbr]=
    (e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1);
  tmp++;
  *num_ep = tmp;
}
 
/*======================== out_vertex ===================*/
/*This routine outputs the x and y coordinates of a vertex*/
/* Passes back to program vertices indexed by vertex number.*/
out_vertex(struct Site *v,float vertx[],float verty[],int *num_vert)
{
  int tmp;
  /*fprintf (out,"v %f %f\n", v->coord.x, v->coord.y);*/
  vertx[v->sitenbr] = v->coord.x;
  verty[v->sitenbr] = v->coord.y;
  tmp = *num_vert;
  tmp++;
  *num_vert = tmp;
 
}
 
/*======================= out_site ==============================*/
/*This routine outputs the x and y coordinates of an input particle*/
/* CIJOL changing a bit... should now pass back coords. indexed*/
/* by site number.*/
out_site(struct Site *s,float rx[],float ry[],int *site_counter)
{
  int tmp;
 
 /*fprintf(out,"s %f %f\n", s->coord.x, s->coord.y);*/
  rx[s->sitenbr] = s->coord.x;
  ry[s->sitenbr] = s->coord.y;
  tmp = *site_counter;
  tmp++;
  *site_counter = tmp;
 
}
 
/*============== out_triple ===========*/
/* For delaunay triangulation */
out_triple(struct Site *s1,struct Site *s2,struct Site *s3,int tripA[],
	   int tripB[],int tripC[],int *num_trip)
{
  int tmp;
  /*fprintf(out,"%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);*/
  tmp = *num_trip;
  tripA[tmp] = s1->sitenbr;
  tripB[tmp] = s2->sitenbr;
  tripC[tmp] = s3->sitenbr;
  tmp++;
  *num_trip = tmp;
 
}
 
/*============== clip_line =============*/
int clip_line(struct Edge *e)
{
  struct Site *s1, *s2;
  float x1,x2,y1,y2;
 
  if(e -> a == 1.0 && e ->b >= 0.0){
    s1 = e -> ep[1];
    s2 = e -> ep[0];
  }
  else{
    s1 = e -> ep[0];
    s2 = e -> ep[1];
  };
 
  if(e -> a == 1.0){
    y1 = pymin;
    if (s1!=(struct Site *)NULL && s1->coord.y > pymin)
      y1 = s1->coord.y;
    if(y1>pymax) return(0);
    x1 = e -> c - e -> b * y1;
    y2 = pymax;
    if (s2!=(struct Site *)NULL && s2->coord.y < pymax)
      y2 = s2->coord.y;
    if(y2<pymin) return(0);
    x2 = e -> c - e -> b * y2;
    if ((x1> pxmax & x2>pxmax) | (x1<pxmin&x2<pxmin)) return(0);
    if(x1> pxmax){
      x1 = pxmax; y1 = (e -> c - x1)/e -> b;};
    if(x1<pxmin){
      x1 = pxmin; y1 = (e -> c - x1)/e -> b;};
    if(x2>pxmax){
      x2 = pxmax; y2 = (e -> c - x2)/e -> b;};
    if(x2<pxmin){
      x2 = pxmin; y2 = (e -> c - x2)/e -> b;};
  }
  else{
    x1 = pxmin;
    if (s1!=(struct Site *)NULL && s1->coord.x > pxmin)
      x1 = s1->coord.x;
    if(x1>pxmax) return(0);
    y1 = e -> c - e -> a * x1;
    x2 = pxmax;
    if (s2!=(struct Site *)NULL && s2->coord.x < pxmax)
      x2 = s2->coord.x;
    if(x2<pxmin) return(0);
    y2 = e -> c - e -> a * x2;
    if ((y1> pymax & y2>pymax) | (y1<pymin&y2<pymin)) return(0);
    if(y1> pymax){
      y1 = pymax; x1 = (e -> c - y1)/e -> a;
    };
    if(y1<pymin){
      y1 = pymin; x1 = (e -> c - y1)/e -> a;
    };
    if(y2>pymax){
      y2 = pymax; x2 = (e -> c - y2)/e -> a;
    };
    if(y2<pymin){
      y2 = pymin; x2 = (e -> c - y2)/e -> a;
    };
  };
 
}
 
/*============== ELinitialize ================*/
ELinitialize()
{
  int i;
 
  freeinit(&hfl, sizeof **ELhash);
  ELhashsize = 2 * sqrt_nsites;
  /* CIJOL Using a fixed size array 8.12.03 */
  /*ELhash = (struct Halfedge **) myalloc ( sizeof *ELhash * ELhashsize);*/
  for(i=0; i<ELhashsize; i +=1) ELhash[i] = (struct Halfedge *)NULL;
  ELleftend = HEcreate( (struct Edge *)NULL, 0);
  ELrightend = HEcreate( (struct Edge *)NULL, 0);
  ELleftend -> ELleft = (struct Halfedge *)NULL;
  ELleftend -> ELright = ELrightend;
  ELrightend -> ELleft = ELleftend;
  ELrightend -> ELright = (struct Halfedge *)NULL;
  ELhash[0] = ELleftend;
  ELhash[ELhashsize-1] = ELrightend;
}
 
/*============ HEcreate ==============*/
struct Halfedge *HEcreate(struct Edge *e,int pm)
{
  struct Halfedge *answer;
  /* CIJOL Use precreated array 8.12.03 */
  /*answer = (struct Halfedge *) getfree(&hfl);*/
  answer=&(hfl_avail[next_free_hfl]);
  next_free_hfl++;

  answer -> ELedge = e;
  answer -> ELpm = pm;
  answer -> PQnext = (struct Halfedge *) NULL;
  answer -> vertex = (struct Site *) NULL;
  answer -> ELrefcnt = 0;
  return(answer);
}
 
/*========== ELinsert ==============*/
ELinsert(struct Halfedge *lb,struct Halfedge *new)
{
  new -> ELleft = lb;
  new -> ELright = lb -> ELright;
  (lb -> ELright) -> ELleft = new;
  lb -> ELright = new;
}
 
/*============ ELgethash ===============*/
/* Get entry from hash table, pruning any deleted nodes */
struct Halfedge *ELgethash(int b)
{
  struct Halfedge *he;
 
  if(b<0 || b>=ELhashsize) return((struct Halfedge *) NULL);
  he = ELhash[b];
  if (he == (struct Halfedge *) NULL ||
      he -> ELedge != (struct Edge *) DELETED ) return (he);
 
  /* Hash table points to deleted half edge.  Patch as necessary. */
  ELhash[b] = (struct Halfedge *) NULL;
  if ((he -> ELrefcnt -= 1) == 0) makefree(he, &hfl);
  return ((struct Halfedge *) NULL);
}
 
/*============== ELleftbnd =================*/
struct Halfedge *ELleftbnd(struct Point *p)
{
  int i, bucket;
  struct Halfedge *he;
 
  /* Use hash table to get close to desired halfedge */
  bucket = (p->x - xmin)/deltax * ELhashsize;
  if(bucket<0) bucket =0;
  if(bucket>=ELhashsize) bucket = ELhashsize - 1;
  he = ELgethash(bucket);
  if(he == (struct Halfedge *) NULL){
    for(i=1; 1 ; i += 1){
      if ((he=ELgethash(bucket-i)) != (struct Halfedge *) NULL) break;
      if ((he=ELgethash(bucket+i)) != (struct Halfedge *) NULL) break;
    };
  };
  /* Now search linear list of halfedges for the corect one */
  if (he==ELleftend  || (he != ELrightend && right_of(he,p))){
    do{
      he = he -> ELright;
    } while (he!=ELrightend && right_of(he,p));
    he = he -> ELleft;
  }
  else
    do {
      he = he -> ELleft;
    } while (he!=ELleftend && !right_of(he,p));
 
  /* Update hash table and reference counts */
  if(bucket > 0 && bucket <ELhashsize-1){
    if(ELhash[bucket] != (struct Halfedge *) NULL)
      ELhash[bucket] -> ELrefcnt -= 1;
    ELhash[bucket] = he;
    ELhash[bucket] -> ELrefcnt += 1;
  };
  return (he);
}
 
/*=============== ELdelete =================*/
/* This delete routine can't reclaim node, since pointers from hash
   table may be present.   */
ELdelete(struct Halfedge *he)
{
  (he -> ELleft) -> ELright = he -> ELright;
  (he -> ELright) -> ELleft = he -> ELleft;
  he -> ELedge = (struct Edge *)DELETED;
}
 
/*============= ELright =============*/
struct Halfedge *ELright(struct Halfedge *he)
{
  return (he -> ELright);
}

/*=========== ELleft =================*/
struct Halfedge *ELleft(struct Halfedge *he)
{
  return (he -> ELleft);
}
 
/*=============== leftreg ===============*/
struct Site *leftreg(struct Halfedge *he)
{
  if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
  return( he -> ELpm == le ?
         he -> ELedge -> reg[le] : he -> ELedge -> reg[re]);
}
 
/*=============== rightreg ===================*/
struct Site *rightreg(struct Halfedge *he)
{
  if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
  return( he -> ELpm == le ?
         he -> ELedge -> reg[re] : he -> ELedge -> reg[le]);
}

/*===================== geominit ================*/
void geominit()
{
  struct Edge e;
  float sn;
 
  freeinit(&efl, sizeof e);
  nvertices = 0;
  nedges = 0;
  sn = nsites+4;
  sqrt_nsites = sqrt(sn);
  deltay = ymax - ymin;
  deltax = xmax - xmin;
}
 
/*===================== bisect ==================*/
struct Edge *bisect(struct Site *s1,struct Site *s2,float linea[],
		    float lineb[],float linec[],int *numlines)
{
  float dx,dy,adx,ady;
  struct Edge *newedge;
 
  /* CIJOL Use precreated array 8.12.03 */
  /*newedge = (struct Edge *) getfree(&efl);*/
  newedge = &(efl_avail[next_free_efl]);
  next_free_efl++;
 
  newedge -> reg[0] = s1;
  newedge -> reg[1] = s2;
  ref(s1);
  ref(s2);
  newedge -> ep[0] = (struct Site *) NULL;
  newedge -> ep[1] = (struct Site *) NULL;
 
  dx = s2->coord.x - s1->coord.x;
  dy = s2->coord.y - s1->coord.y;
  adx = dx>0 ? dx : -dx;
  ady = dy>0 ? dy : -dy;
  newedge -> c = s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5;
  if (adx>ady){
    newedge -> a = 1.0; newedge -> b = dy/dx; newedge -> c /= dx;
  }
  else{
    newedge -> b = 1.0; newedge -> a = dx/dy; newedge -> c /= dy;
  }
 
  newedge -> edgenbr = nedges;
 
  /* Write out location of the bisecting line, part of Voronoi construct.*/
  out_bisector(newedge,linea,lineb,linec,numlines);
 
  nedges += 1;
  return(newedge);
}
 
/*======================== intersect ==============*/
struct Site *intersect(struct Halfedge *el1,struct Halfedge *el2)
{
  struct        Edge *e1,*e2, *e;
  struct  Halfedge *el;
  float d, xint, yint;
  int right_of_site;
  struct Site *v;
 
  e1 = el1 -> ELedge;
  e2 = el2 -> ELedge;
  if(e1 == (struct Edge*)NULL || e2 == (struct Edge*)NULL)
    return ((struct Site *) NULL);
  if (e1->reg[1] == e2->reg[1]) return ((struct Site *) NULL);
 
  d = e1->a * e2->b - e1->b * e2->a;
  if (-1.0e-10<d && d<1.0e-10) return ((struct Site *) NULL);
 
  xint = (e1->c*e2->b - e2->c*e1->b)/d;
  yint = (e2->c*e1->a - e1->c*e2->a)/d;
 
  if( (e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
     (e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
      e1->reg[1]->coord.x < e2->reg[1]->coord.x) ){
    el = el1; e = e1;
  }
  else{
    el = el2; e = e2;
  }
  right_of_site = xint >= e -> reg[1] -> coord.x;
  if ((right_of_site && el -> ELpm == le) ||
      (!right_of_site && el -> ELpm == re)) return ((struct Site *) NULL);
 
  /* CIJOL Use preconstructed array 8.12.03 */
  /*v = (struct Site *) getfree(&sfl);*/
  v = &(sfl_avail[next_free_sfl]);
  next_free_sfl++;

  v -> refcnt = 0;
  v -> coord.x = xint;
  v -> coord.y = yint;
  return(v);
}
 
/*====================== right_of ================*/
/* returns 1 if p is to right of halfedge e */
int right_of(struct Halfedge *el,struct Point *p)
{
  struct Edge *e;
  struct Site *topsite;
  int right_of_site, above, fast;
  float dxp, dyp, dxs, t1, t2, t3, yl;
 
  e = el -> ELedge;
  topsite = e -> reg[1];
  right_of_site = p -> x > topsite -> coord.x;
  if(right_of_site && el -> ELpm == le) return(1);
  if(!right_of_site && el -> ELpm == re) return (0);
 
  if (e->a == 1.0){
    dyp = p->y - topsite->coord.y;
    dxp = p->x - topsite->coord.x;
    fast = 0;
    if ((!right_of_site &e->b<0.0) | (right_of_site&e->b>=0.0) ){
     above = dyp>= e->b*dxp;
      fast = above;
   }
    else{
      above = p->x + p->y*e->b > e-> c;
      if(e->b<0.0) above = !above;
      if (!above) fast = 1;
    };
    if (!fast){
      dxs = topsite->coord.x - (e->reg[0])->coord.x;
      above = e->b * (dxp*dxp - dyp*dyp) <
        dxs*dyp*(1.0+2.0*dxp/dxs + e->b*e->b);
      if(e->b<0.0) above = !above;
    };
  }
  else{  /*e->b==1.0 */
    yl = e->c - e->a*p->x;
    t1 = p->y - yl;
    t2 = p->x - topsite->coord.x;
    t3 = yl - topsite->coord.y;
    above = t1*t1 > t2*t2 + t3*t3;
  };
  return (el->ELpm==le ? above : !above);
}
 
/*=================== endpoint ==================*/
endpoint(struct Edge *e,int lr,struct Site *s,int cjoleft[],int cjoright[],
	 int *num_ep)
{
  e -> ep[lr] = s;
  ref(s);
  if(e -> ep[re-lr]== (struct Site *) NULL) return(0);
  out_ep(e,cjoleft,cjoright,num_ep);
  deref(e->reg[le]);
  deref(e->reg[re]);
  makefree(e, &efl);
}
 
/*==================== dist =======================*/
float dist(struct Site *s,struct Site *t)
{
  float dx,dy;
  dx = s->coord.x - t->coord.x;
  dy = s->coord.y - t->coord.y;
  return(sqrt(dx*dx + dy*dy));
}
 
/*===================== makevertex ===================*/
int makevertex(struct Site *v,float vertx[],float verty[],int *num_vert)
{
  v -> sitenbr = nvertices;
  nvertices += 1;
  out_vertex(v,vertx,verty,num_vert);
}
 
/*================== deref ===============*/
deref(struct Site *v)
{
  v -> refcnt -= 1;
  if (v -> refcnt == 0 ) makefree(v, &sfl);
}
 
/*=================== ref ==============*/
ref(struct Site *v)
{
  v -> refcnt += 1;
}
 
/*================ PQinsert ===============*/
PQinsert(struct Halfedge *he,struct Site *v,float offset)
{
  struct Halfedge *last, *next;
 
  he -> vertex = v;
  ref(v);
  he -> ystar = v -> coord.y + offset;
  last = &PQhash[PQbucket(he)];
  while ((next = last -> PQnext) != (struct Halfedge *) NULL &&
         (he -> ystar   > next -> ystar  ||
          (he -> ystar == next -> ystar && v -> coord.x > next->vertex->coord.x))){
    last = next;
  }
  he -> PQnext = last -> PQnext;
  last -> PQnext = he;
  PQcount += 1;
}
 
/*============== PQdelete ===========*/
PQdelete(struct Halfedge *he)
{
  struct Halfedge *last;
 
  if(he ->  vertex != (struct Site *) NULL){
    last = &PQhash[PQbucket(he)];
    while (last -> PQnext != he) last = last -> PQnext;
    last -> PQnext = he -> PQnext;
    PQcount -= 1;
    deref(he -> vertex);
    he -> vertex = (struct Site *) NULL;
  };
}
 
/*=============== PQbucket =============*/
int PQbucket(struct Halfedge *he)
{
  int bucket;
 
  bucket = (he->ystar - ymin)/deltay * PQhashsize;
  if (bucket<0) bucket = 0;
  if (bucket>=PQhashsize) bucket = PQhashsize-1 ;
  if (bucket < PQmin) PQmin = bucket;
  return(bucket);
}
 
 
/*============== PQempty ==============*/
int PQempty()
{
  return(PQcount==0);
}
 
/*============== PQ_min ===============*/
struct Point PQ_min()
{
  struct Point answer;
 
  while(PQhash[PQmin].PQnext == (struct Halfedge *)NULL) {PQmin += 1;};
  answer.x = PQhash[PQmin].PQnext -> vertex -> coord.x;
  answer.y = PQhash[PQmin].PQnext -> ystar;
  return (answer);
}
 
/*================ PQextractmin ============*/
struct Halfedge *PQextractmin()
{
  struct Halfedge *curr;
 
  curr = PQhash[PQmin].PQnext;
  PQhash[PQmin].PQnext = curr -> PQnext;
  PQcount -= 1;
  return(curr);
}

/*===================== PQinitialize ===========*/
PQinitialize()
{
  int i;
 
  PQcount = 0;
  PQmin = 0;
  PQhashsize = 4 * sqrt_nsites;
  /* CIJOL using fixed array 8.12.03 */
  /*PQhash = (struct Halfedge *) myalloc(PQhashsize * sizeof *PQhash);*/
  for(i=0; i<PQhashsize; i+=1) PQhash[i].PQnext = (struct Halfedge *)NULL;
}
 
/*================ freeinit ===========*/
freeinit(struct Freelist *fl,int size)
{
  fl -> head = (struct Freenode *) NULL;
  fl -> nodesize = size;
}
 
/*=========== makefree ============*/
makefree(struct Freenode *curr,struct Freelist *fl)
{
  curr -> nextfree = fl -> head;
  fl -> head = curr;
}
/*========================= plot_poly_line ==================*/
/* This function draws lines for Voronoi construction on screen. */
/* Called by: plot_voronoi$.*/
void plot_poly_line(float x1,float y1,float x2,float y2)
{
  int par_x1,par_y1,par_x2,par_y2;

  /* Convert x,y location into plotting coordinates */
  par_x1 = (int) ( (x1 + x_offset) * x_scale + BORDER );
  par_y1 = (int) (((float)win_height) - (y1 + y_offset) * y_scale - BORDERY);
  par_x2 = (int) ( (x2 + x_offset) * x_scale + BORDER );
  par_y2 = (int) (((float)win_height) - (y2 + y_offset) * y_scale - BORDERY);
  
  /* Clip lines that fall outside plotting region */
  if (par_x1 < (BORDER - 10))    return;
  if (par_y1 < (BORDERY - 10))   return;
  if (par_x1 >= ((int) win_width) - (BORDER - 10))    return;
  if (par_y1 >= ((int) win_height) - (BORDERY - 10))  return;
  if (par_x2 < (BORDER - 10))    return;
  if (par_y2 < (BORDERY - 10))   return;
  if (par_x2 >= ((int) win_width) - (BORDER - 10))    return;
  if (par_y2 >= ((int) win_height) - (BORDERY - 10))  return;

  XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));
  XDrawLine(display,pixmap,pixmap_gc,(par_x1),(par_y1),(par_x2),(par_y2));
}

/*====================== plot_polygon ==========================*/
/* Plots Voronoi polygons as solid 2D objects, rather than with the older */
/* line-by-line disconnected method. */
void plot_polygon(int index,int sidenum,int polylist[][20],float vertx[],
		  float verty[],float xcenter,float ycenter)
{
  float radial(float xin,float yin,float xcenter,float ycenter);

  int i,j;
  float x1,y1,a,b,c;
  /* CIJOL: type XPoint consists of pairs of shorts, x y. */
  float xlist[20],ylist[20];
  XPoint points[20];
  float angle[20];

  /* Generate list of vertices */
  for(i=0;i<sidenum;i++){
    xlist[i] = vertx[polylist[index][i]];
    ylist[i] = verty[polylist[index][i]];
  }

  /* To plot the polygons properly, the vertices need to be in */
  /* counterclockwise order.  Convert x and y coordinates of each vertex */
  /* to a positive angle. */
  for(i=0;i<sidenum;i++){
    angle[i] = radial(xlist[i],ylist[i],xcenter,ycenter);
  }
  /* Sort angles into ascending order using standard algorithm for */
  /* small N from Numerical Recipes */
  for(j=1;j<sidenum;j++){
    a = angle[j];
    b = xlist[j];
    c = ylist[j];
    i = j-1;
    while(i>-1&&angle[i]>a){
      angle[i+1]=angle[i];
      xlist[i+1]=xlist[i];
      ylist[i+1]=ylist[i];
      i--;
    }
    angle[i+1]=a;
    xlist[i+1]=b;
    ylist[i+1]=c;
  }

  /* Convert x,y location into plotting coordinates */
  for(i=0;i<sidenum;i++){
    points[i].x =(short)((xlist[i] + x_offset) * x_scale + BORDER );
    points[i].y =(short)(((float)win_height)-
			 (ylist[i]+y_offset)*y_scale - BORDERY);
    
    /* Clip lines that fall outside plotting region */
    if (points[i].x < (BORDER - 10))    return;
    if (points[i].y < (BORDERY - 10))   return;
    if (points[i].x >= ((int) win_width) - (BORDER - 10))    return;
    if (points[i].y >= ((int) win_height) - (BORDERY - 10))  return;
  }

  switch(sidenum){
  case 6:
    XSetForeground(display,pixmap_gc,(WhitePixel(display,screen_num)));
    break;
  case 5:
    XSetForeground(display,pixmap_gc,c_map[4]);
    break;
  case 7:
    XSetForeground(display,pixmap_gc,c_map[5]);
    break;
  default:
    XSetForeground(display,pixmap_gc,c_map[6]);
    break;
  }
  XSetFillStyle(display,pixmap_gc,FillSolid);
  XFillPolygon(display,pixmap,pixmap_gc,points,sidenum,Nonconvex,
	       CoordModeOrigin); 
}
float radial(float xin,float yin,float xcenter,float ycenter)
{
  float theta;
  float x,y;

  /* Use relative coordinates centered on vortex in middle of polygon. */
  x = xin-xcenter;
  y = yin-ycenter;

  /* Check for x-values close to zero and set to +- pi/2 */
  if(fabs(x)<0.0001){
    if(y>0)
      return(PI/2.0);
    else
      return(3.0*PI/2.0);
  }

  /* For non-zero x-values: */
  theta = (float)atan((double)(y/x));

  /*Arctan returns values between +- pi/2; move values to correct */
  /* quadrant and make all angle values positive. */
  if(x<0)
    theta += PI;
  else if (y<0)
    theta += 2.0*PI;
  
  return(theta);
}
  

/*========================= plot_triangle ==================*/
/* This function draws Delaunay triangles on screen, and colors them */
/* according to number of neighbors. */
/* Called by: plot_voronoi$.*/
void plot_triangle(float x1,float y1,float x2,float y2,float x3,float y3,
		   int sideA,int sideB,int sideC,float dmin,float dmax,
		   float dstore[][3],int index,int max_color)
{
  /* CIJOL: type XPoint consists of pairs of shorts, x y. */
  XPoint points[3];
  int color[3];
  float diffd=dmax-dmin;

  /* New: Determine the color of each of the three lines based on */
  /* their relative lengths compared to dmin,dmax. */
  /* Depending on the colormap you're using you may want a reversed value. */
  /*color[0]=(1.0-(dstore[index][0]-dmin)/diffd)*max_color;*/
  /*color[1]=(1.0-(dstore[index][1]-dmin)/diffd)*max_color;*/
  /*color[2]=(1.0-(dstore[index][2]-dmin)/diffd)*max_color;*/
  color[0]=((dstore[index][0]-dmin)/diffd)*(max_color-1)+1;
  color[1]=((dstore[index][1]-dmin)/diffd)*(max_color-1)+1;
  color[2]=((dstore[index][2]-dmin)/diffd)*(max_color-1)+1;

  /* Convert x,y location into plotting coordinates */
  points[0].x = (short) ( (x1 + x_offset) * x_scale + BORDER );
  points[0].y = (short) (((float)win_height)-(y1+y_offset)*y_scale - BORDERY);
  points[1].x = (short) ( (x2 + x_offset) * x_scale + BORDER );
  points[1].y = (short) (((float)win_height)-(y2+y_offset)*y_scale - BORDERY);
  points[2].x = (short) ( (x3 + x_offset) * x_scale + BORDER );
  points[2].y = (short) (((float)win_height)-(y3+y_offset)*y_scale - BORDERY);
  
  /* Clip lines that fall outside plotting region */
  if (points[0].x < (BORDER - 10))    return;
  if (points[0].y < (BORDERY - 10))   return;
  if (points[0].x >= ((int) win_width) - (BORDER - 10))    return;
  if (points[0].y >= ((int) win_height) - (BORDERY - 10))  return;
  if (points[1].x < (BORDER - 10))    return;
  if (points[1].y < (BORDERY - 10))   return;
  if (points[1].x >= ((int) win_width) - (BORDER - 10))    return;
  if (points[1].y >= ((int) win_height) - (BORDERY - 10))  return;
  if (points[2].x < (BORDER - 10))    return;
  if (points[2].y < (BORDERY - 10))   return;
  if (points[2].x >= ((int) win_width) - (BORDER - 10))    return;
  if (points[2].y >= ((int) win_height) - (BORDERY - 10))  return;

  /* For black lines, use next line only, and remove other XSetForegrounds: */
  /*XSetForeground(display,pixmap_gc,(BlackPixel(display,screen_num)));*/
  XSetForeground(display,pixmap_gc,c_map[color[0]]);
  XDrawLine(display,pixmap,pixmap_gc,(points[0].x),(points[0].y),
	    (points[1].x),(points[1].y));
  XSetForeground(display,pixmap_gc,c_map[color[2]]);
  XDrawLine(display,pixmap,pixmap_gc,(points[0].x),(points[0].y),
	    (points[2].x),(points[2].y));
  XSetForeground(display,pixmap_gc,c_map[color[1]]);
  XDrawLine(display,pixmap,pixmap_gc,(points[2].x),(points[2].y),
	    (points[1].x),(points[1].y));

  /* CIJOL Disabling the triangle coloring. */
#ifdef NOTUSING
  /* Now we need to figure out how many neighbors the vortices at the */
  /* vertices of this triangle have.  The idea is that we want to fill the */
  /* triangle with a contrasticng color IF one or more of the vortices is */
  /* a defect site.  The numbering of the triangles is conveyed through the */
  /* tripleA,B,C arrays.  A single element of each array is passed to this */
  /* routine, which can then directly use the sidenum structure to find the */
  /* number of sides associated with each vertex of the triangle.*/
  if((sideA!=6)||(sideB!=6)||(sideC!=6)){
    XSetForeground(display,pixmap_gc,c_map[4]);
  }
  else{
    XSetForeground(display,pixmap_gc,(WhitePixel(display,screen_num)));
  }

  XSetFillStyle(display,pixmap_gc,FillSolid);
  XFillPolygon(display,pixmap,pixmap_gc,points,3,Nonconvex,CoordModeOrigin); 
#endif
}

void plot_vertex(float x,float y,int sideA,int sideB,int sideC)
{
  int par_x,par_y;
  int Box = 3;
  int Box2 = 6;
  int color;

  color = 4;

  /* Convert x,y location into plotting coordinates */
  par_x = (int) ((x + x_offset) * x_scale) + BORDER;
  par_y = (int) (((float)win_height)-(y+y_offset)*y_scale - BORDERY);
  
  /* Clip vertices that fall outside plotting region */
  if (par_x < (BORDER - 10))    return;
  if (par_y < (BORDERY - 10))   return;
  if (par_x >= ((int) win_width) - (BORDER - 10))    return;
  if (par_y >= ((int) win_height) - (BORDERY - 10))  return;

  if((sideA!=6)||(sideB!=6)||(sideC!=6)){
    XSetForeground(display,pixmap_gc,c_map[4]);
    XDrawArc(display,pixmap,pixmap_gc,par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
    XFillArc(display,pixmap,pix_gc[color],par_x-Box,par_y-Box,Box2,Box2,0,360*64); 
  }
  /*else{
    XSetForeground(display,pixmap_gc,(WhitePixel(display,screen_num)));
  }*/
}

