#include<iostream>
#include<cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

// GSL functions for solving ordinary differential equations
#include <gsl/gsl_odeiv2.h>

using namespace std;

/* 
 Define the array of right-hand-side functions y[i] to be integrated. 
 For single ode:
    y'(x) = y²(x)
 y(x) --> y[0]
 dy[0]/dx=f[0]  =  y[0]*y[0];
*/
int rhs (double, const double y[], double f[], void*) {
	f[0] = y[0]*y[0];
	return GSL_SUCCESS; /* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

/*
 * Exact (analytical) solution with boundary condition y(0)=y0
 */
double yexact(double x, double y0) {
    return 1.0 / ( (1.0/y0)  - x );
}

int main (void)
{
    double y0 = -2.0;  /* initial condition, y(x=0) */
    
    const double h = 1.0e-3;       /* starting step size for ode solver */
    const int n_fixed_steps = 1;   /* steps performed by applying driver once */
    const int n_loops = 1000;
    
    const double eps_abs = 1.e-3;  /* absolute error requested  */
    const double eps_rel = 1.e-6;  /* relative error requested  */

    const int dimension = 1;   /* number of differential equations */
    double y[dimension];	        /* current solution vector */
    y[0] = y0;	/* initial value of y(x) */
    
    gsl_odeiv2_system ode_system;
    ode_system.function = rhs;       /* the right-hand-side of equation */
    ode_system.dimension = dimension;  /* number of diffeq's */
    ode_system.params = NULL;   /* parameters to pass to dfunc */
    ode_system.jacobian = NULL; /* the Jacobian df[i]/dy[j] ; needed only by some methods */

    gsl_odeiv2_driver *d =
        gsl_odeiv2_driver_alloc_y_new (&ode_system, gsl_odeiv2_step_rk4, h, eps_abs, eps_rel);
    
    double x = 0.0;

    printf ("# x y_num y_exact err_abs err_rel\n");
    printf ("%7.6g %.6e %.6e %.6e %.6e\n", x, y0, y0, 0.0, 0.0 );

    int status;
    for (int i = 0; i < n_loops; i++)
    {
        status = gsl_odeiv2_driver_apply_fixed_step (d, &x, h, n_fixed_steps, y);
        if (status != GSL_SUCCESS)
        {
            printf ("### error: driver returned %d\n", status);
            break;
        }
        
        double y_num = y[0];
        double y_exact = yexact(x,y0);
        double err_abs = y_num - y_exact;
        double err_rel = err_abs / y_exact;

        printf ("%7.6g %.6e %.6e %.6e %.6e\n", x, y_num, y_exact, err_abs, err_rel );
    }

  gsl_odeiv2_driver_free (d);
  return status;
}
