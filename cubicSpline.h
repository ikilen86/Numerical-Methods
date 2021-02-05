
#ifndef __CUBIC_SPLINE_H_INCLUDED__
#define __CUBIC_SPLINE_H_INCLUDED__

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
			
			double dx0_tmp = grid_dx[0];

/*
			// GUPPE always has uniform grid
			for(int i =0; i< num_x-1;i++)
			{
				if (abs(grid_dx[i]-dx0_tmp) > 1.0e-12)
				{
					// Mostly because of speed, simply uncomment the lines
					//i0 = misc_get_k_index(kp, grid_x, num_x);
					// in the evaluate functions to enable support for non-uniform grids
					cout << "CubicSpline(num_x, x): Currently not evalulating on non-uniform grid_x" << endl;
					for(int k = 0; k < Nx; k++)
					{
						cout << "x = " << grid_x[k] << endl;
					}
					
					cout << "quitting.." << endl;
					exit(-1);
				}
			}
*/
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
		
		
		//! Generate spline interpolation object
		/*! Estimate derivatives at endpoints, fill/solve tri-diagonal matrix, and prepare interpolation coeff.
		 * \param y function values to be interpolated. Length must be equal to the values set during initialization.
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
		
		double evaluate(double kp)
		{
			int i0;
			if (kp < grid_x[num_x-1]) // select correct interval
			{
				double kp_indx = (kp-grid_x[0])/grid_dx[0];
				i0 = floor(kp_indx);
				//i0 = misc_get_k_index(kp, grid_x, num_x);
				
			} else {
				i0 = num_x-2; // Last interval
			}

			//i0 = misc_get_k_index(kp, grid_x, num_x);
			
			double kp_x0_1 = (kp  -grid_x[i0]);
			double x1_kp_1 = (grid_x[i0+1]-kp);
			double kp_x0_3 = kp_x0_1*kp_x0_1*kp_x0_1;
			double x1_kp_3 = x1_kp_1*x1_kp_1*x1_kp_1;
			
			double y0 =  cubic_spline_coeff_p[i0][0]*x1_kp_3 + cubic_spline_coeff_p[i0][1]*kp_x0_3 + cubic_spline_coeff_p[i0][2]*kp_x0_1 + cubic_spline_coeff_p[i0][3]*x1_kp_1;
			
			return y0;
		}
		
		double evaluate_df(double kp)
		{
			int i0;
			if (kp < grid_x[num_x-1]) // select correct interval
			{
				double kp_indx = (kp-grid_x[0])/grid_dx[0];
				i0 = floor(kp_indx);
				
			} else {
				i0 = num_x-2; // Last interval
			}
//			i0 = misc_get_k_index(kp, grid_x, num_x);
			
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
	
		/* INPUT: A - array of pointers to rows of a square matrix having dimension N
		 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
		 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
		 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
		 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
		 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
		 */
		int LUPDecompose(double **A, int N, double Tol, int *P) 
		{

			int i, j, k, imax; 
			double maxA, *ptr, absA;

			for (i = 0; i <= N; i++)
				P[i] = i; //Unit permutation matrix, P[N] initialized with N

			for (i = 0; i < N; i++) {
				maxA = 0.0;
				imax = i;

				for (k = i; k < N; k++)
					if ((absA = fabs(A[k][i])) > maxA) { 
						maxA = absA;
						imax = k;
					}

				if (maxA < Tol) return 0; //failure, matrix is degenerate

				if (imax != i) {
					//pivoting P
					j = P[i];
					P[i] = P[imax];
					P[imax] = j;

					//pivoting rows of A
					ptr = A[i];
					A[i] = A[imax];
					A[imax] = ptr;

					//counting pivots starting from N (for determinant)
					P[N]++;
				}

				for (j = i + 1; j < N; j++) {
					A[j][i] /= A[i][i];

					for (k = i + 1; k < N; k++)
						A[j][k] -= A[j][i] * A[i][k];
				}
			}

			return 1;  //decomposition done 
		}


		/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
		 * OUTPUT: x - solution vector of A*x=b
		 * Example:
			A = Mat(N,N);
			int permutations[N+1] = {0, 0, 0, 0};
			int sig = LUPDecompose(A, N, 1e-8, permutations);
			if (sig!=1)
			{
				cout << "CubicSpline::prepare()LU:Decomposition faliure.." << endl;
				exit(-1);
			}
			A_inv = Mat(N,N);
			LUPInvert(A, permutations, N, A_inv);
		 */
		void LUPSolve(double **A, int *P, double *b, int N, double *x) 
		{

			for (int i = 0; i < N; i++) {
				x[i] = b[P[i]];

				for (int k = 0; k < i; k++)
					x[i] -= A[i][k] * x[k];
			}

			for (int i = N - 1; i >= 0; i--) {
				for (int k = i + 1; k < N; k++)
					x[i] -= A[i][k] * x[k];

				x[i] = x[i] / A[i][i];
			}
		}

		/* INPUT: A,P filled in LUPDecompose; N - dimension
		 * OUTPUT: IA is the inverse of the initial matrix
		 */
		void LUPInvert(double **A, int *P, int N, double **IA)
		{
			for (int j = 0; j < N; j++) {
				for (int i = 0; i < N; i++) {
					if (P[i] == j) 
						IA[i][j] = 1.0;
					else
						IA[i][j] = 0.0;
					for (int k = 0; k < i; k++)
						IA[i][j] -= A[i][k] * IA[k][j];
				}
				
				for (int i = N - 1; i >= 0; i--) {
					for (int k = i + 1; k < N; k++)
						IA[i][j] -= A[i][k] * IA[k][j];
					IA[i][j] = IA[i][j] / A[i][i];
				}
			}
		}

		/* INPUT: A,P filled in LUPDecompose; N - dimension. 
		 * OUTPUT: Function returns the determinant of the initial matrix
		 */
		double LUPDeterminant(double **A, int *P, int N) {

			double det = A[0][0];

			for (int i = 1; i < N; i++)
				det *= A[i][i];

			if ((P[N] - N) % 2 == 0)
				return det; 
			else
				return -det;
		}
	
		/* Return index 'i' of kp in array x such that x[i] <= kp < x[i+1]
		 * Assumes array is ordered such that K[0] < K[1] < ... < K[numEl-1]
		 * Does NOT check if kp is contained inside x[0] <= kp < x[numEl-1]
		 */
		int misc_get_k_index(double kp, double *x, int numEl)
		{
				int ia = 0;
				int ib = numEl-1;
				int ic = 0;
				int counter = 0;
				int counter_max = 10000;
				while (counter<counter_max)
				{
						ic = floor((ia+ib)/2.0);
						if ((x[ic] <= kp)&&(kp < x[ic+1]))
						{
								return ic;
						}

						if (x[ic] < kp)
						{
								ia = ic;
						} else {
								ib = ic;
						}

						counter++;
				}

				cout << "misc_get_k_index() Could not find element inside array!" << endl;
				cout << "kp = " << kp << endl;
				cout << "x[0] = " << x[0] << ", x[" << numEl-1 << "] = " << x[numEl-1] << endl;
				exit(-1);
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
