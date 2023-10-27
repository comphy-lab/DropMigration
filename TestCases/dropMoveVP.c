// same as dropMove.c but proper non-dimensionalization

#define MIN_LEVEL 0
#define MAX_LEVEL 8

#define VelErr 1e-3
#define FErr 1e-3
#define cErr 1e-3
#define KErr 1e-3

#define tsnap 1e-1

#include "navier-stokes/centered.h"
#define FILTERED
#include "../src/two-phase-clsvof-VP.h"
#include "integral.h"
#include "../src/activity.h"
// #include "curvature.h"

scalar cL[],  *stracers = {cL};
#define c0 0.0

cL[top] = dirichlet(c0);
cL[right] = dirichlet(c0);
cL[left] = dirichlet(c0);
f[bottom] = dirichlet(0.0);

scalar * list = NULL;
int ny, nx; 
double Deltay, Deltax;
double dtmax, tmax; 

scalar sigmaf[];

#define Oh 1e0
#define GammaR 1e0
#define Ma 1e1
#define AcNum 1e0
#define J 1e-1

int main(){
  stokes = true;
  // dtmax = 1e-2;
  L0 = 16.0;
  origin (-0.5*L0, -0.5*L0);
  // origin (0.,0.);
  N = 1 << MAX_LEVEL;
  init_grid (N);

  d.sigmaf = sigmaf;

  rho1 = 1e0/sq(Oh); rho2 = 1e0/sq(Oh);
  tmax = 25.;
  mu1 = 1.0; mu2 = 1.0;
  mumax = (1e4)*mu2;
  tauy = J;


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

event outputs (t = 0.; t += tsnap; t <= tmax) { 
  char dumpFile[160];
  sprintf (dumpFile, "intermediate/snapshot-%5.4f", t);
  dump (file = dumpFile);
}
scalar cTest[];

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
  // assert(ke > 1e-10);
  assert(ke < 1e3);
  }
}