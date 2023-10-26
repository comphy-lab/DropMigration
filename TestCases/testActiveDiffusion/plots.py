# -*- coding: utf-8 -*-
# #define MIN_LEVEL 0
# #define MAX_LEVEL 10
# 
# #define VelErr 1e-3
# #define FErr 1e-3
# #define cErr 1e-3
# #define KErr 1e-3
# 
# #define tsnap 1e0
# 
# #include "navier-stokes/centered.h"
# #define FILTERED
# #include "two-phase-clsvof.h"
# // #include "integral.h"
# #include "curvature.h"
# #include "../src/activity.h"
# 
# scalar cL[],  *stracers = {cL};
# #define c0 0.0
# 
# cL[top] = dirichlet(c0);
# cL[right] = dirichlet(c0);
# 
# scalar * list = NULL;
# int ny, nx; 
# double Deltay, Deltax;
# double dtmax, tmax; 
# 
# int main(){
#   // dtmax = 1e-2;
#   L0 = 8.0;
#   // origin (-0.5*L0, -0.5*L0);
#   origin (0.,0.);
#   N = 1 << MAX_LEVEL;
#   init_grid (N);
# 
#   // d.sigma = 0.;
#   rho1 = 1e0; rho2 = 1e0;
#   tmax = 100.;
#   // mu1 = 1.0; mu2 = 1.0;
# 
#   cL.inverse = true;
#   cL.A = 1e0;
#   cL.D = 1e0;
# 
#   char comm[160];
#   sprintf (comm, "rm -rf intermediate");
#   system(comm);
#   sprintf (comm, "mkdir -p intermediate");
#   system(comm);
# 
#   // Output features
#   ny = 1000;
#   Deltay = (double)((L0 - 0.)/(ny));
#   // fprintf(ferr, "%g\n", Deltay);
#   nx = (int)((L0 - 0.)/Deltay);
#   // fprintf(ferr, "%d\n", nx);
#   Deltax = (double)((L0 - 0.)/(nx));
#   // fprintf(ferr, "%g\n", Deltax);
#   list = list_add (list, f);
#   list = list_add (list, cL);
#   
#   run();  
# 
# }
# 
# 
# event init (i = 0) {
#   // fraction (fphi, sq(rin) - sq(x) - sq(y));
#   double rout = 0.999*L0;
#   foreach() {
#     d[] = 1. - sqrt (sq(x) + sq(y));
#     u.x[] = 0.0;
#     u.y[] = 0.0;
#     cL[] = sq(x) + sq(y) > sq(1.) ? (sq(x) + sq(y) > sq(rout) ? c0 : c0 - cL.A*log(sqrt(sq(x) + sq(y))/rout)) : 0.;
#   }
# }
# 
# event properties(i++){
#   double rout = 0.999*L0;
#   foreach(){
#     cL[] = sq(x) + sq(y) > sq(rout) ? c0 : cL[];
#   }
# }
# 
# scalar KAPPA[];
# event adapt(i++){  
#   curvature(f, KAPPA);
#   adapt_wavelet({f, u.x, u.y, cL, KAPPA}, (double[]){FErr, VelErr, VelErr, cErr, KErr}, MAX_LEVEL, MIN_LEVEL);
# }
# 
# event outputs (t = 0.; t += tsnap; t <= tmax) { 
# 
#   char dumpFile[160];
#   sprintf (dumpFile, "intermediate/snapshot-%5.4f", t);
#   dump (file = dumpFile);
# 
#   char filenameOut[160];
#   sprintf (filenameOut, "intermediate/fieldData-%5.4f.dat", t);
# 
#   FILE * fp = fopen (filenameOut, "w");
#   fprintf (fp, "y f cL\n");
# 
#   for (int j = 0; j < ny; j++){
#     double y = Deltay*(j) + 0.;
#     fprintf (fp, "%g", y);
# 
#     for (scalar s in list){
#       fprintf (fp, " %g", interpolate(s, 0.0, y));
#     }
#     fprintf (fp, "\n");
#   }
# 
#   fflush (fp);
#   fclose (fp);
#   fprintf (ferr, "Data written to %s\n", filenameOut);
#   // return 1;
# }
# 
# scalar cTest[];
# 
# event logWriting (i++) {
#   
#   double MaxErr = 0.;
#   if (i == 0){
#   foreach()
#     cTest[] = cL[];
#   } else {
#     foreach()
#       cTest[] = cL[] - cTest[];
#     
#     MaxErr = normf(cTest).avg;
#     assert(MaxErr > 1e-3);
# 
#     foreach()
#       cTest[] = cL[];
#   }
# 
# 
# 
#   static FILE * fp;
#   if (i == 0) {
#     fprintf (ferr, "i dt t H\n");
#     fp = fopen ("log.dat", "w");
#     fprintf (fp, "i dt t H\n");
#     fprintf (fp, "%d %g %g %6.5e\n", i, dt, t, MaxErr);
#     fclose(fp);
#   } else {
#     fp = fopen ("log.dat", "a");
#     fprintf (fp, "%d %g %g %6.5e\n", i, dt, t, MaxErr);
#     fclose(fp);
#   }
#   fprintf (ferr, "%d %g %g %6.5e\n", i, dt, t, MaxErr);
# 
#   // dump (file = "dump");
#   // if (ke < 1e-6 && i > 100){
#   //   fprintf(ferr, "kinetic energy too small now! Stopping!\n");
#   //   fp = fopen ("log", "a");
#   //   fprintf(fp, "kinetic energy too small now! Stopping!\n");
#   //   fclose(fp);
#   //   return 1;
#   // }
# }
# 
