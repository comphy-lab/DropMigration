/**
 * Simulation of a two-phase droplet system with surfactant effects
 * 
 * This code simulates a droplet with surfactant-modified surface tension
 * using the CLSVOF (Coupled Level Set and Volume of Fluid) method.
 * The system is non-dimensionalized with the following parameters:
 * - Oh (Ohnesorge number): Ratio of viscous to inertial and surface tension forces
 * - Ma (Marangoni number): This is more like the Peclet number
 * - GammaR: Base surface tension coefficient
 * - AcNum: Activity number for surfactant effects
 */

#define MIN_LEVEL 0
#define MAX_LEVEL 7

#define VelErr 1e-3
#define FErr 1e-3
#define cErr 1e-3
#define KErr 1e-3

#define tsnap 1e-1

#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase-clsvof.h"
#include "integral.h"
#include "activity.h"
// #include "curvature.h"

/**
 * Global variables and boundary conditions
 * cL: Surfactant concentration field
 * sigmaf: Surface tension coefficient field
 */
scalar cL[],  *stracers = {cL};
#define c0 0.0

cL[top] = neumann(0.);
cL[right] = neumann(0.);
cL[left] = neumann(0.);
cL[bottom] = neumann(0.);
f[bottom] = dirichlet(0.0);

u.t[top] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

scalar * list = NULL;
int ny, nx; 
double Deltay, Deltax;
double dtmax, tmax; 

scalar sigmaf[];

/**
 * Non-dimensional parameters
 */
#define Oh 1e0
#define GammaR 1e0
#define Ma 1e1
#define AcNum 1e0

/**
 * Main function: Sets up and runs the simulation
 */
int main(){
  stokes = true;
  // dtmax = 1e-2;
  L0 = 8.0;
  origin (-0.5*L0, -0.5*L0);
  // origin (0.,0.);
  N = 1 << MAX_LEVEL;
  init_grid (N);

  d.sigmaf = sigmaf;

  rho1 = 1e0/sq(Oh); rho2 = 1e0/sq(Oh);
  tmax = 50.;
  mu1 = 1.0; mu2 = 1.0;

  cL.inverse = true;
  cL.A = AcNum;
  cL.D = 1e0/Ma;

  char comm[160];
  sprintf (comm, "rm -rf intermediate");
  system(comm);
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  run();  

}

/**
 * Initialization event
 * Sets up initial conditions for:
 * - Distance function (d)
 * - Velocity fields (u.x, u.y)
 * - Surfactant concentration (cL)
 * - Surface tension coefficient (sigmaf)
 */
event init (i = 0) {
  // fraction (fphi, sq(rin) - sq(x) - sq(y));
  foreach() {
    d[] = 1. - sqrt (sq(x) + sq(y));
    u.x[] = 0.0;
    u.y[] = 0.0;
    cL[] = c0; //sq(x) + sq(y) > sq(1.) ? (sq(x) + sq(y) > sq(rout) ? c0 : c0 - cL.A*log(sqrt(sq(x) + sq(y))/rout)) : 0.;
    sigmaf[] = GammaR + cL[];
  }
}

/**
 * Properties update event
 * Updates the surface tension coefficient based on surfactant concentration
 */
event properties(i++){
  foreach(){
    sigmaf[] = GammaR + cL[];
  }
}

scalar KAPPA[];
event adapt(i++){  

  foreach()
    KAPPA[] = distance_curvature(point, d);
  
  adapt_wavelet({f, u.x, u.y, cL, KAPPA}, (double[]){FErr, VelErr, VelErr, cErr, KErr}, MAX_LEVEL, MIN_LEVEL);
}

/**
 * Output event
 * Saves simulation snapshots at regular intervals
 */
event outputs (t = 0.; t += tsnap; t <= tmax) { 
  char dumpFile[160];
  sprintf (dumpFile, "intermediate/snapshot-%5.4f", t);
  dump (file = dumpFile);
}
scalar cTest[];

/**
 * Logging event
 * Computes and logs kinetic energy of the system
 * Also performs assertions to check simulation stability
 */
event logWriting (i++) {
  
  double ke = 0.;
  foreach(reduction(+:ke)){
    ke += 0.5*rho(f[])*(sq(u.x[])+sq(u.y[]))*sq(Delta);
  }

  static FILE * fp;
  if (pid() == 0){
    if (i == 0){
      fprintf (ferr, "i t ke\n");
      fp = fopen ("log.dat", "w");
      fprintf (fp, "i t ke\n");
      fclose (fp);
    }
    fprintf (ferr, "%d %g %5.5e\n", i, t, ke);
    fp = fopen ("log.dat", "a");
    fprintf (fp, "%d %g %5.5e\n", i, t, ke);
    fclose (fp);
  }

  if (i > 10){
  assert(ke > 1e-6);
  assert(ke < 1e3);
  }
}
