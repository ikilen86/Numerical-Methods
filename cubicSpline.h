#ifndef __CUBIC_SPLINE_H_INCLUDED__
#define __CUBIC_SPLINE_H_INCLUDED__

/* A class for spline interpolation of 1d functions on a uniform grid.
 * This is designed for use when multiple evlautions is executed on a single curve.
 * The endpoint derivatives are estimated from the endpoint values using 5 order or lower finite difference formulas
 * Example: The following code interpolates a function defined by vectors x[Nx], y[Nx]
 * CubicSpline *test = new CubicSpline(Nx, x);
 * test->prepare(y);
 * for(int i = 0; i < numX; i++)
 * 	yp[i] = test->evaluate(xp[i])
 **/

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <complex>

using namespace std;

class CubicSpline
{
	public:
		//! An empty constructor
		/*! Initialize the spline object with null pointers. No arrays are allocated.
		 **/
		CubicSpline()
		{
			cubic_spline_diag_a = NULL;
			cubic_spline_diag_b = NULL;
			cubic_spline_diag_c = NULL;
			cubic_spline_diag_d = NULL;
			cubic_spline_z = NULL;
			cubic_spline_coeff_p = NULL;
			grid_x = NULL;
			grid_dx = NULL;
			num_x = 0;
		}
		
		//! A constructor
		/*! Initialize the spline object
		 * \param Nx Number of x-grid points
		 * \param x x-grid to evaluate over 
		 **/
		CubicSpline(int Nx, double *x)
		{
			cubic_spline_diag_a = NULL;
			cubic_spline_diag_b = NULL;
			cubic_spline_diag_c = NULL;
			cubic_spline_diag_d = NULL;
			cubic_spline_z = NULL;
			cubic_spline_coeff_p = NULL;
			grid_x = NULL;
			grid_dx = NULL;

			num_x = Nx;
			initialize_arrays();
			
			for(int i = 0; i < num_x; i++)
			{
				grid_x[i] = x[i];
			}
			
			for(int i =0; i< num_x-1;i++)
			{
				grid_dx[i] = grid_x[i+1]-grid_x[i];
			}

			// Tmp array
			double dX_cs[num_x-2];
			for(int i =0; i< num_x-2;i++)
			{
				dX_cs[i] = grid_dx[i+1]+grid_dx[i];
			}
			
			// Setup diagonal matrix
			cubic_spline_diag_b[0]       = 2.0*grid_dx[0];
			cubic_spline_diag_b[num_x-1] = 2.0*grid_dx[num_x-2];
			for(int i = 1; i < num_x-1; i++)
			{
				cubic_spline_diag_b[i] = 2.0*dX_cs[i-1];
			}
			memcpy(cubic_spline_diag_a, grid_dx, (num_x-1)*sizeof(double));
			memcpy(cubic_spline_diag_c, grid_dx, (num_x-1)*sizeof(double));
			
			// Prepare diag_c for Tri. Diag. Solver
			cubic_spline_diag_c[0] = cubic_spline_diag_c[0]/cubic_spline_diag_b[0];
			for(int i = 1; i < num_x-1; i++)
			{
				cubic_spline_diag_c[i] = cubic_spline_diag_c[i]/(cubic_spline_diag_b[i] - cubic_spline_diag_a[i-1]*cubic_spline_diag_c[i-1]);
			}
			
			for(int i = 0; i < num_x; i++)
			{
				cubic_spline_diag_d[i] = 0.0;
				cubic_spline_z[i] = 0.0;
			}
			
			// Prepare storage for coeff.
			for(int i = 0; i < num_x-1; i++)
			{
				for(int j = 0; j < 4; j++)
				{
					cubic_spline_coeff_p[i][j] = 0.0;
				}
			}
		}
		
		//! Destructor
		~CubicSpline()
		{
			if (cubic_spline_diag_a != NULL)
			{
				delete [] cubic_spline_diag_a;
				delete [] cubic_spline_diag_b;
				delete [] cubic_spline_diag_c;
				delete [] cubic_spline_diag_d;
				delete [] cubic_spline_z;
				delete [] grid_x;
				delete [] grid_dx;
				for(int i = 0; i < num_x-1; i++)
				{
					delete [] cubic_spline_coeff_p[i];
				}
				delete [] cubic_spline_coeff_p;
			}
		}
		
		//! Copy constructor
		CubicSpline(const CubicSpline &obj)
		{
			num_x = obj.num_x;
			if (num_x > 0)
			{
				initialize_arrays();
				
				for(int i = 0; i < num_x-1; i++)
				{
					cubic_spline_diag_a[i]  = obj.cubic_spline_diag_a[i];
					cubic_spline_diag_c[i]  = obj.cubic_spline_diag_c[i];
					grid_dx[i] 		= obj.grid_dx[i];
				}
				
				for(int i = 0; i < num_x; i++)
				{
					cubic_spline_diag_b[i]	= obj.cubic_spline_diag_b[i];
					cubic_spline_diag_d[i]	= obj.cubic_spline_diag_d[i];
					cubic_spline_z[i]	= obj.cubic_spline_z[i];
					grid_x[i]		= obj.grid_x[i];
				}
				
				for(int i = 0; i < num_x-1; i++)
				{
					for(int j = 0; j < 4; j++)
					{
						cubic_spline_coeff_p[i][j] = obj.cubic_spline_coeff_p[i][j];
					}
				}
			}
		}
		
		
		//! Generate spline interpolation object for discreete values
		/*! Estimate derivatives at endpoints, fill/solve tri-diagonal matrix, and prepare interpolation coeff.
		 *  Call evaluate(x) to get the interpolated values at a given point.
		 * \param y function values to be interpolated. Length must be equal to the values set during initialization.
		 * \sa evaluate
		 */
		void prepare(double *y)
		{
			// Estimate derivatives at endpoints
			double dx = grid_dx[0];
			double df_0, df_N;
			if ( num_x == 2)
			{
				df_0 = (y[0] 	    - y[1]	)/dx;
				df_N = -(y[num_x-1] - y[num_x-2])/dx;
			} else if (num_x == 3)
			{
				// 2nd Order
				df_0 =  (-3.0*y[0]       + 4.0*y[1]       - y[2]	)/(2.0*dx);	
				df_N = -(-3.0*y[num_x-1] + 4.0*y[num_x-2] - y[num_x-3]  )/(2.0*dx);
			} else {
				// 4th Order
				df_0 =  (-25.0*y[0]       + 48.0*y[1]       - 36.0*y[2]       + 16.0*y[3]       - 3.0*y[4]      )/(12.0*dx); 
				df_N = -(-25.0*y[num_x-1] + 48.0*y[num_x-2] - 36.0*y[num_x-3] + 16.0*y[num_x-4] - 3.0*y[num_x-5])/(12.0*dx);
			}

			// prepare array
			cubic_spline_diag_d[0]       = 6.0*(y[1]-y[0])/grid_dx[0] - 6.0*df_0;
			cubic_spline_diag_d[num_x-1] = 6.0*df_N - 6.0*(y[num_x-1]-y[num_x-2])/grid_dx[num_x-2];
			for(int i = 1; i < num_x-1; i++)
			{
				cubic_spline_diag_d[i] = 6.0*(y[i+1]-y[i])/grid_dx[i] - 6.0*(y[i]-y[i-1])/grid_dx[i-1];
			}
			
			// solver: fill cubic_spline_z
			tri_diag_solver();
			
			// Store coeff.
			for(int i = 0; i < num_x-1; i++)
			{
				cubic_spline_coeff_p[i][0] = cubic_spline_z[i]/(6.0*grid_dx[i]);
				cubic_spline_coeff_p[i][1] = cubic_spline_z[i+1]/(6.0*grid_dx[i]);
				cubic_spline_coeff_p[i][2] = y[i+1]/grid_dx[i] - cubic_spline_z[i+1]*grid_dx[i]/6.0;
				cubic_spline_coeff_p[i][3] = y[i]/grid_dx[i] - cubic_spline_z[i]*grid_dx[i]/6.0;
			}
		}
		
		//! Evalutate spline interpolation object
		/*! Should be called after evaluate(). 
		 * \param kp x-value where interpolation object is evaluated.
		 * \return Spline polynomial evaluated at given point
		 * \sa prepare
		 **/
		double evaluate(double kp)
		{
			// Find correct spline polynomial
			int i0;
			if (kp < grid_x[num_x-1]) // select correct interval
			{
				// uniform grid
				double kp_indx = (kp-grid_x[0])/grid_dx[0];
				i0 = floor(kp_indx);
				
			} else if (kp < grid_x[0]) {
				i0 = 0; // First polynomial
			} else {
				i0 = num_x-2; // Last polynomial
			}
			
			// Evaluate polynomial 
			double kp_x0_1 = (kp  -grid_x[i0]);
			double x1_kp_1 = (grid_x[i0+1]-kp);
			double kp_x0_3 = kp_x0_1*kp_x0_1*kp_x0_1;
			double x1_kp_3 = x1_kp_1*x1_kp_1*x1_kp_1;
			
			double y0 =  cubic_spline_coeff_p[i0][0]*x1_kp_3 + cubic_spline_coeff_p[i0][1]*kp_x0_3 + cubic_spline_coeff_p[i0][2]*kp_x0_1 + cubic_spline_coeff_p[i0][3]*x1_kp_1;
			
			return y0;
		}
		
		//! Evalutate derivative of spline interpolation object
		/*! Should be called after evaluate(). 
		 * \param kp x-value where interpolation object is evaluated.
		 * \return Derivative of spline polynomial evaluated at given point
		 * \sa prepare
		 **/
		double evaluate_df(double kp)
		{
			// Find correct spline polynomial
			int i0;
			if (kp < grid_x[num_x-1]) // select correct interval
			{
				// uniform grid
				double kp_indx = (kp-grid_x[0])/grid_dx[0];
				i0 = floor(kp_indx);
				
			} else if (kp < grid_x[0]) {
				i0 = 0; // First polynomial
			} else {
				i0 = num_x-2; // Last polynomial
			}
			
			// Evaluate derivative of polynomial 
			double kp_x0_1 = (kp  -grid_x[i0]);
			double x1_kp_1 = (grid_x[i0+1]-kp);
			double kp_x0_3 = 3.0*kp_x0_1*kp_x0_1;
			double x1_kp_3 = -3.0*x1_kp_1*x1_kp_1;
			
			double y0 =  cubic_spline_coeff_p[i0][0]*x1_kp_3 + cubic_spline_coeff_p[i0][1]*kp_x0_3 + cubic_spline_coeff_p[i0][2] - cubic_spline_coeff_p[i0][3];
			
			return y0;
		}
		
	private:
		//! Allocate space for interpolation object
		void initialize_arrays()
		{
			if (cubic_spline_diag_a == NULL)
			{
				cubic_spline_diag_a = new double[num_x-1];
				cubic_spline_diag_b = new double[num_x];
				cubic_spline_diag_c = new double[num_x-1];
				cubic_spline_diag_d = new double[num_x];
				cubic_spline_z = new double[num_x];
				grid_x = new double[num_x];
				grid_dx = new double[num_x-1];
				cubic_spline_coeff_p = new double*[num_x-1];
				for(int i = 0; i < num_x-1; i++)
				{
					cubic_spline_coeff_p[i] = new double[4];
					for(int j = 0; j < 4; j++)
					{
						cubic_spline_coeff_p[i][j] = 0.0;
					}
				}
			}
		}
	
		//! Thomas algorithm for simplified Gaussian elimination and back subst.
		void tri_diag_solver()
		{
			cubic_spline_diag_d[0] = cubic_spline_diag_d[0]/cubic_spline_diag_b[0];
			for(int i = 1; i < num_x; i++)
			{
				cubic_spline_diag_d[i] = (cubic_spline_diag_d[i] - cubic_spline_diag_a[i-1]*cubic_spline_diag_d[i-1])/(cubic_spline_diag_b[i]-cubic_spline_diag_a[i-1]*cubic_spline_diag_c[i-1]);
			}
			
			// Back subst.
			cubic_spline_z[num_x-1] = cubic_spline_diag_d[num_x-1];
			for(int i = num_x-2; i >=0; i--)
			{
				cubic_spline_z[i] = cubic_spline_diag_d[i] - cubic_spline_diag_c[i]*cubic_spline_z[i+1];
			}
		}
			
		int num_x;
		double *grid_x;
		double *grid_dx;
		double *cubic_spline_diag_a;
		double *cubic_spline_diag_b;
		double *cubic_spline_diag_c;
		double *cubic_spline_diag_d;
		double *cubic_spline_z;
		double **cubic_spline_coeff_p;
};


#endif
