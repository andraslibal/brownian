# brownian
This is a generic Brownian Dynamics code for colloidal and granular systems

First version: no optimization
2D system of point like-particles
Each particle has x,y,fx,fy, q, color
Initialization: random initialization
Periodic Boundary condition in both x,y directions
1/r^2 exp(-r/r0) repulsion between particles
q*fx0 driving force in the x direction (simulating 50-50% going in opposite directions)
Movie output: creating a movie file with all x,y positions written out in each frame
Statistics output: output the average speed of movement avg(q*fx) to a statistics file
avergaed over all particles and averaged over the time steps the average considers (ex.1000)  
