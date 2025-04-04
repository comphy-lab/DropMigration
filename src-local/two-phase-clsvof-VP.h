/**
# Two-phase interfacial flows with coupled VOF and levelset

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The signed distance field `d` is advected as a tracer and is
relaxed toward the VOF-defined interface.

This coupling ensures a mass conservation at least as good as that of
the [pure VOF solver](two-phase.h).

This solver can be combined with the [integral formulation of surface
tension](integral.h).

The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"
#include "tracer.h"

scalar d[], f[], * interfaces = {f}, * tracers = {d};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double mumax = 0., tauy = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(muTemp, muDrop, f)  (clamp(f,0.,1.)*(muDrop - muTemp) + muTemp)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

#include "fractions.h"

event properties (i++)
{
  /**
This is part where we have made changes.
$$D_{11} = \frac{\partial u}{\partial x}$$
$$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{22} = \frac{\partial v}{\partial y}$$
The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
$$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
the equivalent viscosity is
$$\mu_{eq}= \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{N-1} + \frac{\tau_y}{\sqrt{2} D_2 }$$
**Note:** $\|D\| = D_2/\sqrt{2}$.<br/>
For Bingham Fluid, N = 1:
$$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$.

The fluid flows always, it is not a solid, but a very viscous fluid.
$$ \mu = min\left(\mu_{eq}, \mu_{\text{max}}\right) $$
*/
  double muTemp = mu2;
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      double D11 = (u.x[] - u.x[-1,0]);
      double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
      double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
      double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
      if (D2 > 0. && tauy > 0.) {
        double temp = tauy/(sqrt(2.0)*D2) + mu2;
        muTemp = min(temp, mumax);
      } else {
        if (tauy > 0.){
          muTemp = mumax;
        } else {
          muTemp = mu2;
        }
      }
      /**
      Note that only the heavier fluid is Viscoplastic.
      */
      muv.x[] = fm.x[]*mu(muTemp, mu1, ff);
    }
  }
  
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}


/**
The initial volume fraction is computed from the initial distance
field, which must be initialised by the user. */

event init (i = 0)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  fractions (phi, f);
}

/**
The distance function is reinitialised at each timestep. */

#include "redistance.h"

event properties (i++)
{

  /**
  In interfacial cells, the signed distance is obtained directly from
  the VOF reconstruction of the interface. This distance is combined
  with the existing distance using a small weight, thus ensuring
  exponential time relaxation of the signed distance toward its VOF
  value. */

  double weight = 0.1;
  foreach()
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      normalize (&n);
      double alpha = plane_alpha (f[], n);
      d[] = (1. - weight)*d[] + weight*Delta*alpha;
    }

  /**
  The redistancing operation itself is quite expensive. */
  
  redistance (d, imax = 3);
}

/**
## See also

* [Two-phase interfacial flows with VOF](two-phase.h)
* [Two-phase interfacial flows with levelset](two-phase-levelset.h)
*/
