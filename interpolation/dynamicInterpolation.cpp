#ifndef __DYNAMICINTERPOLATION_CPP_INCLUDED__
#define __DYNAMICINTERPOLATION_CPP_INCLUDED__

/* Interpolate a function T = f(x,y,z) where size(T) = Nk
 * Uses Trilinear interpolation, and is thus exact for functions
 * of the form: 1+x+y+z+xy+xz+yz+xyz
 * 
 * When file_save() (or file_load()) is called
 * the program saves (or loads) the current interpolation grid to file (or from file).
 * The filename is the variable 'newName' given to the program at init.
 * All points are compared up to ERROR_TOL and no duplicates are written to file
 * 
 * Tip: When increasing the resolution use the same grids as before 
 * but with dx/n where n is a whole number, then the old points are loaded from file.
 * 
 * Warning: 
 * 1. When evaluating a point that contains the upper domain boundary (x1,y1 or z1) the program crashes 
 *  	This can be fixed with cases, or one keeps away from the boundary.
 * 2. The Program crashes if a point outside of the domain is used
 *  	However, it isnt much of a problem to expand the domain dynamically
 * */

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "fileIO.cpp"


// Enable to debug
//#define DYANMIC_INTERPOLATION_DIAGNOSTIC

using namespace std;


class DynamicInterpolation
{
	public:
		//! An empty constructor
		DynamicInterpolation()
		{
			//F = NULL;
			FXYZ = NULL;
			fxyz = NULL;

			x_target_ind = -1;
			y_target_ind = -1;
			z_target_ind = -1;
			x_target = NULL;
			y_target = NULL;
			z_target = NULL;
		}
		
		//! A constructor
		/*! Construct interpolation space
		 * \param newName Filename used for storing interpolation space variables
		 * \param Nk Number of output variables
		 * \param x0 Input parameter 1 min value
		 * \param x1 Input parameter 1 max value
		 * \param dx Input parameter 1 interpolation resolution
		 * \param y0 Input parameter 2 min value
		 * \param y1 Input parameter 2 max value
		 * \param dy Input parameter 2 interpolation resolution
		 * \param z0 Input parameter 3 min value
		 * \param z1 Input parameter 3 max value
		 * \param dz Input parameter 3 interpolation resolution
		 **/
		DynamicInterpolation(const std::string & newName, int Nk, double x0, double x1, double dx, double y0, double y1, double dy, double z0, double z1, double dz)
		{
			ERROR_TOL = 1.0e-12;
			
			if (dx <= 0)
			{
				cout << "DynamicInterpolation:: Need dx > 0" << endl;
				cout << "dx = " << dx << endl;
				exit(-1);
			}
			
			if (dy <= 0)
			{
				cout << "DynamicInterpolation:: Need dy > 0" << endl;
				cout << "dy = " << dy << endl;
				exit(-1);
			}
			
			if (dz <= 0)
			{
				cout << "DynamicInterpolation:: Need dz > 0" << endl;
				cout << "dz = " << dz << endl;
				exit(-1);
			}
			
			if (x0 >= x1)
			{
				cout << "DynamicInterpolation:: Need x0 < x1" << endl;
				cout << "x0 = " << x0 << endl;
				cout << "x1 = " << x1 << endl;
				cout << "dx = " << dx << endl;
				exit(-1);
			}
			
			if (y0 >= y1)
			{
				cout << "DynamicInterpolation:: Need y0 < y1" << endl;
				cout << "y0 = " << y0 << endl;
				cout << "y1 = " << y1 << endl;
				cout << "dy = " << dy << endl;
				exit(-1);
			}
			
			if (z0 >= z1)
			{
				cout << "DynamicInterpolation:: Need z0 < z1" << endl;
				cout << "z0 = " << z0 << endl;
				cout << "z1 = " << z1 << endl;
				cout << "dz = " << dz << endl;
				exit(-1);
			}
			
			file_name = newName;
			
			//F = f;
			
			X0 = x0;
			X1 = x1;
			DX = dx;
			
			Y0 = y0;
			Y1 = y1;
			DY = dy;
			
			Z0 = z0;
			Z1 = z1;
			DZ = dz;
			
			num_x = floor((X1-X0)/DX)+1;
			num_y = floor((Y1-Y0)/DY)+1;
			num_z = floor((Z1-Z0)/DZ)+1;
			
			num_output = Nk;
			FXYZ = new double***[num_x];
			for(int i = 0; i < num_x; i++)
			{
				FXYZ[i] = new double**[num_y];
				for(int j = 0; j < num_y; j++)
				{
					FXYZ[i][j] = new double*[num_z];
					for(int k = 0; k < num_z; k++)
					{
						FXYZ[i][j][k] = NULL;
					}
				}
			}
			
			// The edges of the interpolation cube and their index
			x_target = new double[2];
			y_target = new double[2];
			z_target = new double[2];
			x_target_ind = -1;
			y_target_ind = -1;
			z_target_ind = -1;
			
			// TMP storage for interpolation
			fxyz = new double***[2];
			for(int i = 0; i < 2; i++)
			{
				fxyz[i] = new double**[2];
				for(int j = 0; j < 2; j++)
				{
					fxyz[i][j] = new double*[2];
					for(int k = 0; k < 2; k++)
					{
						fxyz[i][j][k] = NULL;
					}
				}
			}
			
			// Try to find as many old data points as possible
			cout << "DynamicInterpolation:: Trying to load as many values as possible from file: " << file_name << endl;
			int num = file_load();
			cout << "DynamicInterpolation:: Found and loaded " << num << " points from file" << endl;
		}
		
		//! Destructor
		~DynamicInterpolation()
		{
			if (x_target != NULL)
			{
				delete [] x_target;
				delete [] y_target;
				delete [] z_target;
				
				for(int i = 0; i < 2; i++)
				{
					for(int j = 0; j < 2; j++)
					{
						delete [] fxyz[i][j];
						
					}
					delete [] fxyz[i];
				}
				delete [] fxyz;
				
				for(int i = 0; i < num_x; i++)
				{
					for(int j = 0; j < num_y; j++)
					{
						for(int k = 0; k < num_z; k++)
						{
							if (FXYZ[i][j][k] != NULL)
							{
								delete [] FXYZ[i][j][k];
							}
						}
						delete [] FXYZ[i][j];
					}
					delete [] FXYZ[i];
				}
				delete [] FXYZ;
			}
			
		}
	
		//! Copy constructor
		DynamicInterpolation(const DynamicInterpolation &obj)
		{
			file_name = obj.file_name;
			X0 = obj.X0;
			X1 = obj.X1;
			DX = obj.DX;
			num_x = obj.num_x;
		
			Y0 = obj.Y0;
			Y1 = obj.Y1;
			DY = obj.DY;
			num_y = obj.num_y;

			Z0 = obj.Z0;
			Z1 = obj.Z1;
			DZ = obj.DZ;
			num_z = obj.num_z;
			
			x_target_ind = obj.x_target_ind;
			y_target_ind = obj.y_target_ind;
			z_target_ind = obj.z_target_ind;
			x_target = NULL;
			y_target = NULL;
			z_target = NULL;
			
			num_output = obj.num_output;
			ERROR_TOL = obj.ERROR_TOL;

			FXYZ = NULL;
			fxyz = NULL;

			if (obj.FXYZ != NULL)
			{
				
				FXYZ = new double***[num_x];
				for(int i = 0; i < num_x; i++)
				{
					FXYZ[i] = new double**[num_y];
					for(int j = 0; j < num_y; j++)
					{
						FXYZ[i][j] = new double*[num_z];
						for(int k = 0; k < num_z; k++)
						{
							FXYZ[i][j][k] = NULL;
						}
					}
				}
				
				// The edges of the interpolation cube and their index
				x_target = new double[2];
				y_target = new double[2];
				z_target = new double[2];
				x_target_ind = -1;
				y_target_ind = -1;
				z_target_ind = -1;
				
				// TMP storage for interpolation
				fxyz = new double***[2];
				for(int i = 0; i < 2; i++)
				{
					fxyz[i] = new double**[2];
					for(int j = 0; j < 2; j++)
					{
						fxyz[i][j] = new double*[2];
						for(int k = 0; k < 2; k++)
						{
							fxyz[i][j][k] = NULL;
						}
					}
				}
				
			}

		}
		
		//! Initialize interpolation object
		/*! Given that WE want to calculate at the point (x,y,z), the interpolation requests that we also evaluate the points in new_points. At most num_num_points need to be calculated.
		 * \param x the first coordinate we wish to evaluate
		 * \param y the second coordinate we wish to evalute
		 * \param z the third coordinate we wish to evaluate
		 * \param new_points A [8x6] list of new points to evaluate, where only the first 'num_new_points' points along first dimension are required. The second dimension contains [x,y,z,x_ind,y_ind,z_ind] where the last 3 are needed for update_interpolation_grid()
		 * \param num_new_points The number of new points to evaluate
		 * \sa update_interpolation_grid
		 */
		void prepare_interpolation(double x, double y, double z, double **new_points, int *num_new_points)
		{
			if ((x < X0)||(x >= X1))
			{
				cout << "DynamicInterpolation::prepare_interpolation() x is out of bounds, increase domain size" << endl;
				cout << "x  = " << x << endl;
				cout << "x0 = " << X0 << endl;
				cout << "x1 = " << X1 << endl;
				exit(-1);
			}
			
			if ((y < Y0)||(y >= Y1))
			{
				cout << "DynamicInterpolation::prepare_interpolation() y is out of bounds, increase domain size" << endl;
				cout << "y  = " << y << endl;
				cout << "y0 = " << Y0 << endl;
				cout << "y1 = " << Y1 << endl;
				exit(-1);
			}
			
			if ((z < Z0)||(z >= Z1))
			{
				cout << "DynamicInterpolation::prepare_interpolation() z is out of bounds, increase domain size" << endl;
				cout << "z  = " << z << endl;
				cout << "z0 = " << Z0 << endl;
				cout << "z1 = " << Z1 << endl;
				exit(-1);
			}
			
			// Find correct index
			x_target_ind 	= floor((x-X0)/DX);
			x_target[0] 	= X0 + x_target_ind*DX;
			x_target[1] 	= X0 + (x_target_ind+1.0)*DX;
			
			
			y_target_ind 	= floor((y-Y0)/DY);
			y_target[0] 	= Y0 + y_target_ind*DY;
			y_target[1] 	= Y0 + (y_target_ind+1.0)*DY;
			
			z_target_ind 	= floor((z-Z0)/DZ);
			z_target[0] 	= Z0 + z_target_ind*DZ;
			z_target[1] 	= Z0 + (z_target_ind+1.0)*DZ;
			
			
			// 1. Check if x CAN be interpolated from previous data
			*num_new_points = 0;
			for(int i =  0; i < 2; i++)
			{
				for(int j = 0; j < 2; j++)
				{
					for(int k = 0; k < 2; k++)
					{
						if (FXYZ[x_target_ind+i][y_target_ind+j][z_target_ind+k] == NULL)
						{
						#ifdef DYANMIC_INTERPOLATION_DIAGNOSTIC
							cout << "DynamicInterpolation::prepare_interpolation() Request new point at (x"<< i <<", y"<< j <<", z"<< k <<") = (" << x_target[i] << ", " << y_target[j] << ", " << z_target[k] << ") ... " << endl;
						#endif
							// Generate new point
							new_points[*num_new_points][0] = x_target[i];
							new_points[*num_new_points][1] = y_target[j];
							new_points[*num_new_points][2] = z_target[k];
							new_points[*num_new_points][3] = x_target_ind+i;
							new_points[*num_new_points][4] = y_target_ind+j;
							new_points[*num_new_points][5] = z_target_ind+k;
							*num_new_points += 1;
						#ifdef DYANMIC_INTERPOLATION_DIAGNOSTIC
						} else {
							cout << "DynamicInterpolation::prepare_interpolation() Found old point (x"<< i <<", y"<< j <<", z"<< k <<") = (" << x_target[i] << ", " << y_target[j] << ", " << z_target[k] << ") " << endl;
						#endif
						}
					}
				}
			}
		}
		
		//! Add data to the interpolation grid at x,y,z
		/*! Update stored interpolation object with new data
		 * \param x_ind First coordinate index of a single gridpoint found from prepare_interpolation()
		 * \param y_ind Second coordinate index of a single gridpoint found from prepare_interpolation()
		 * \param z_ind Third coordinate index of a single gridpoint found from prepare_interpolation()
		 * \param new_data A vector of length num_output of data that should be stored.
		 * \sa prepare_interpolation
		 * */
		void update_interpolation_grid(int x_ind, int y_ind, int z_ind, double *new_data)
		{
			if (FXYZ[x_ind][y_ind][z_ind] == NULL)
			{
			#ifdef DYANMIC_INTERPOLATION_DIAGNOSTIC
				cout << "DynamicInterpolation::update_interpolation_grid() Adding new point (" << x << ", " << y << ", " << z << ")" << endl;
			#endif
				FXYZ[x_ind][y_ind][z_ind] = new double[num_output];
				for(int i = 0; i < num_output; i++)
				{
					FXYZ[x_ind][y_ind][z_ind][i] = new_data[i];
				}
				
			} else {
				cout << "DynamicInterpolation::update_interpolation_grid() Point already calculated, cannot update" << endl;
				cout << "x = " << x_ind << endl;
				cout << "y = " << y_ind << endl;
				cout << "z = " << z_ind << endl;
				exit(-1);
			}
			
		}
		
		//! Evaluate the function at the point x,y,z
		/*! Requires that one calls:
		 * 1. prepare_interpolation(): In order to figure out if new points are needed in the interpolation grid
		 * 2. User calculates all new points outside of this program
		 * 3. update_interpolation_grid(): To update the function
		 * then call this function with:
		 * \param x The first coordinate that we want to calculate
		 * \param y The second coordinate that we want to calculate
		 * \param z The third coordinate that we want to calculate
		 * \return output Output vector of length num_output
		 * */
		void evalF(double *output, double x, double y, double z)
		{
			// Find correct index
			x_target_ind 	= floor((x-X0)/DX);
			x_target[0] 	= X0 + x_target_ind*DX;
			x_target[1] 	= X0 + (x_target_ind+1.0)*DX;
			
			
			y_target_ind 	= floor((y-Y0)/DY);
			y_target[0] 	= Y0 + y_target_ind*DY;
			y_target[1] 	= Y0 + (y_target_ind+1.0)*DY;
			
			z_target_ind 	= floor((z-Z0)/DZ);
			z_target[0] 	= Z0 + z_target_ind*DZ;
			z_target[1] 	= Z0 + (z_target_ind+1.0)*DZ;
			
			for(int i =  0; i < 2; i++)
			{
				for(int j = 0; j < 2; j++)
				{
					for(int k = 0; k < 2; k++)
					{	
						// Store pointer to target data array
						fxyz[i][j][k] = FXYZ[x_target_ind+i][y_target_ind+j][z_target_ind+k];
					}
				}
			}
			
			
			// 3. Interpolate and return
			// 3D - Trilinear interpolation
			double dy0 = y - y_target[0];
			double dy1 = y_target[1]-y;
			double dx0 = x - x_target[0];
			double dx1 = x_target[1]-x;
			
			double DXDY = DX*DY; // For speed
			double dz_target = (z_target[1]-z_target[0]); // For speed
			for(int i = 0; i < num_output; i++)
			{
				double interp_z0 = dy1*(fxyz[0][0][0][i]*dx1 + fxyz[1][0][0][i]*dx0) + dy0*(fxyz[0][1][0][i]*dx1 + fxyz[1][1][0][i]*dx0);
				double interp_z1 = dy1*(fxyz[0][0][1][i]*dx1 + fxyz[1][0][1][i]*dx0) + dy0*(fxyz[0][1][1][i]*dx1 + fxyz[1][1][1][i]*dx0);
				
				double r1 = (z-z_target[0])/dz_target;
				double interp = interp_z0 + (interp_z1 - interp_z0)*r1;
				output[i] = interp/(DXDY);
			}
		}
		
		//! Load from file, any data points that correspond to the current grid
		//! File has data stored as: "x y z data_array" on each line
		int file_load()
		{
			int num_points_found = 0;
			
			// Test if file exists, otherwise do nothing
			if (fileExists(file_name))
			{
				// Read lines one by one and testgrid
				FILE *fid = fopen(file_name.c_str(),"r");
				if (fid == NULL)
				{
					cout << "DynamicInterpolation::file_load() Cannot open file" << endl;
					exit(-1);
				}
				
				double xf, yf, zf;
				double *tmp = new double[num_output];
				while (fscanf(fid, "%lg %lg %lg", &xf, &yf, &zf) == 3)
				{
					// Ensure that the gridpoint fits into the current grid
					bool in_grid = true;
					if ((xf < X0)||(xf >= X1))
					{
						in_grid = false;
					}
					
					if ((yf < Y0)||(yf >= Y1))
					{
						in_grid = false;
					}
					
					if ((zf < Z0)||(zf >= Z1))
					{
						in_grid = false;
					}
					
					if (in_grid)
					{	
					
						// Is point in my current grid?
						double ind_x = (xf-X0)/DX;
						double ind_y = (yf-Y0)/DY;
						double ind_z = (zf-Z0)/DZ;
						
						if (((fabs(round(ind_x)-ind_x)<ERROR_TOL)&&(fabs(round(ind_y)-ind_y))<ERROR_TOL)&&(fabs(round(ind_z)-ind_z)<ERROR_TOL))
						{
							// Read rest of line
							for(int i = 0; i < num_output; i++)
							{
								int id = fscanf(fid," %lg",&(tmp[i]));
							}
							
							num_points_found++;
							
							//cout << "DynamicInterpolation::file_load(): "<< num_points_found << " found (" << xf << ", " << yf << ", " << zf << ") = (" << round(ind_x) << ", " << round(ind_y) << ", " << round(ind_z) << ")" << endl;
							
							// Put into correct location
							update_interpolation_grid(round(ind_x), round(ind_y), round(ind_z), tmp);
							
							
						} else {
							cout << "DynamicInterpolation::file_load() Gridpoint doesnt match. skipping: " << ind_x << ", " << ind_y << ", " << ind_z << endl;
							cout << "x: " << ind_x << ", floor(x) = " << round(ind_x) << ", diff = " << fabs(ind_x-round(ind_x)) << endl; 
							cout << "y: " << ind_y << ", floor(y) = " << round(ind_y) << ", diff = " << fabs(ind_y-round(ind_y)) << endl; 
							cout << "z: " << ind_z << ", floor(z) = " << round(ind_z) << ", diff = " << fabs(ind_z-round(ind_z)) << endl; 

							// Skip rest of line
							int id = fscanf(fid, "%*[^\n]\n", NULL);
						}
					}
				}
				
				delete [] tmp;
				
			}
			
			return num_points_found;
		}
		
		//! Save current grid to file, do not overwrite old data
		void file_save()
		{
			// Append new data to back of old file, or create if does not exist
			FILE *fid = fopen(file_name.c_str(),"a+");
			if (fid != NULL)
			{
				// Check if each point exists in file already
				for(int i = 0; i < num_x; i++)
				{
					for(int j = 0; j < num_y; j++)
					{
						for(int k = 0; k < num_z; k++)
						{
							// Check if point exists in grid
							if (FXYZ[i][j][k] != NULL)
							{
								double x = X0 + i*DX;
								double y = Y0 + j*DY;
								double z = Z0 + k*DZ;
								
								// Check if point is in the file already
								if (!file_does_point_exist(x,y,z))
								{
									fprintf(fid,"% .16e % .16e % .16e",x,y,z);
									for(int m = 0; m < num_output; m++)
									{
										fprintf(fid," % .16e",FXYZ[i][j][k][m]);
									}
									fprintf(fid,"\n");
								}
							}
						}
					}
				}
				
			} else {
				cout << "DynamicInterpolation::file_save() Could not open file for saving, out of storage??" << endl;
			}
			fclose(fid);
			
		}
		
	private:
		
		//! Return true if point (x,y,z) exists in the file
		//! Otherwise return false
		bool file_does_point_exist(double x,double y,double z)
		{
			
			FILE *fid = fopen(file_name.c_str(),"r");
			if (fid == NULL)
			{
				cout << "DynamicInterpolation::file_does_not_exist() Cannot open file for checking" << endl;
				exit(-1);
			}
			
			double xf, yf, zf;
			while (fscanf(fid, "%lg %lg %lg", &xf, &yf, &zf) == 3)
			{
				double dist = sqrt((xf-x)*(xf-x) + (yf-y)*(yf-y) + (zf-z)*(zf-z));
				if (dist<ERROR_TOL)
				{
					fclose(fid);
					return true;
				}
				
				// Read until newline, put whatever is there nowhere
				int id = fscanf(fid, "%*[^\n]\n", NULL);
				
			}
			
			fclose(fid);
			return false;
		}

		std::string file_name;
	
		double X0;
		double X1;
		double DX;
		int num_x;
		
		double Y0;
		double Y1;
		double DY;
		int num_y;
		
		double Z0;
		double Z1;
		double DZ;
		int num_z;
		
		int x_target_ind;
		int y_target_ind;
		int z_target_ind;
		double *x_target;
		double *y_target;
		double *z_target;
		
		int num_output;
		double ****FXYZ;
		double ****fxyz;
		
		double ERROR_TOL; // File IO error tolerance
};

#endif
