#ifndef __solver3_h__
#define __solver3_h__

#include <iostream>
#include <string>
#include <cstdlib> // exit()
#include <cstring> //memset

using namespace std;


// Flags
//#define DBG_ISAK_REPORT_ODE_SOLVER_STEPS // Print ode solver step sizes and # function evaluations to file

/*
 *
 *
 */


/* Error estimation safety coefficient for methods that use step  
 * doubling for error estimates. Error estimates are multiplied by  
 * this constant to ensure that the error of a step is not  
 * underestimated.  
 *  
 * The default safety value of 8.0 ensures 90% of samples lie within  
 * the error (assuming a Gaussian distribution with prior  
 * p(sigma) = 1 / sigma). Value of 1.0 conforms to equation  
 * by Ascher and Petzold (reference: Ascher, U.M., Petzold, L.R.,  
 * Computer methods for ordinary differential and  
 * differential-algebraic equations, SIAM, Philadelphia, 1998).  */ 
const double solver3_ERR_SAFETY = 8.0;


class ODE_Solver {
	public:

		ODE_Solver() {
			method_name = "";
			dimension = 0;
			absolute_error_tol = 0;
			relative_error_tol = 0;
			adaptive = false;
			num_func_eval = 0.0;
			function = NULL;
			jacobian = NULL;
			selected_step_method = NULL;
			selected_step_method_LTE_order = 0.0;
			selected_step_method_max_step_increase = 0.0;
			selected_step_method_storage_size = 0;
			maximal_step_size = HUGE_VAL;
			minimal_step_size = -HUGE_VAL;

			num_steps = 0.0;
			num_steps_min = 0.0;
			avg_step_size = 0.0;
			file_step_output = NULL;
			file_step_output_num = 0;

			prev_err_ratio = -1.0;
			prev_step_size = HUGE_VAL;
			fsal_time = HUGE_VAL;
			fsal_p = NULL;
			

			solver_state.ytmp = NULL;
			solver_state.y0 = NULL;
			solver_state.y = NULL;
			solver_state.t0 = HUGE_VAL;
			solver_state.t = HUGE_VAL;
			solver_state.y_int = NULL;
			solver_state.yerr = NULL;
			solver_state.y0_save = NULL;
			solver_state.interp_t0 = HUGE_VAL;
			solver_state.interp_t1 = HUGE_VAL;
		}

		~ODE_Solver() {

			if (file_step_output != NULL)
			{
				fprintf(file_step_output,"Total number of steps = %.0f\n",num_steps);
				fprintf(file_step_output,"Avg. step size steps  = %.3e\n",avg_step_size);
				fprintf(file_step_output,"number of min-steps   = %.0f\n",num_steps_min);
				fprintf(file_step_output,"number of func. eval  = %.0f\n",num_func_eval);
				fclose(file_step_output);
				file_step_output = NULL;
			}

			if (solver_state.yerr != NULL)
			{
				delete [] solver_state.yerr;
			}

			if (solver_state.ytmp != NULL)
			{
				delete [] solver_state.ytmp;
				delete [] solver_state.y;
				delete [] solver_state.y0;
			}

			if (solver_state.y0_save)
			{
				delete [] solver_state.y0_save;
			}
	

			if (solver_state.y_int != NULL)
			{
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					delete [] solver_state.y_int[i];
				}
				delete [] solver_state.y_int;
			}
		}

		void Reset() {
			num_func_eval = 0.0;

			num_steps = 0.0;
			num_steps_min = 0.0;
			avg_step_size = 0.0;

			prev_err_ratio = -1.0;
			prev_step_size = HUGE_VAL;
			fsal_time = HUGE_VAL;
			fsal_p = NULL;

			solver_state.t = HUGE_VAL;
			solver_state.t0 = HUGE_VAL;
			solver_state.interp_t0 = HUGE_VAL;
			solver_state.interp_t1 = HUGE_VAL;

			#ifdef DBG_ISAK_REPORT_ODE_SOLVER_STEPS
			if (file_step_output != NULL)
			{
				file_output_close();
			}
			std::ostringstream fname("");
			fname << "ode_solver_"<< method_name << "_step_info_" << file_step_output_num << ".dat";
			file_step_output = fopen(fname.str().c_str(),"w+");
			file_step_output_num += 1;
			file_step_output_num = (file_step_output_num % 2); // [0,1]
			#endif


			if (solver_state.ytmp != NULL)
			{
				memset(solver_state.ytmp,0.0,dimension*sizeof(double));
			}

			if (solver_state.y0 != NULL)
			{
				memset(solver_state.y0,0.0,dimension*sizeof(double));
			}

			if (solver_state.y != NULL)
			{
				memset(solver_state.y,0.0,dimension*sizeof(double));
			}

			if (solver_state.yerr != NULL)
			{
				memset(solver_state.yerr,0.0,dimension*sizeof(double));
			}

			if (solver_state.y0_save)
			{
				memset(solver_state.y0_save,0.0,dimension*sizeof(double));
			}

			
			if (solver_state.y_int != NULL)
			{
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					memset(solver_state.y_int[i],0.0,dimension*sizeof(double));
				}
			}
		}

		void Init(const char *method,  double tol_abs, double tol_rel, double min_step, double max_step, void *pars, int n, int (*func)(double , const double *, double *, void *),  int (*jac)(double , const double *, double *, double *, void *), bool adapt) 
		{

			method_name = method;
			dimension = n;
			function = func;
			jacobian = jac;
			if (jac != NULL)
			{
				std::cerr << "ERROR: ode solver must have jac=NULL in this version. Quitting!" << endl;
				exit(-1);
			}
	
			adaptive = adapt;
			maximal_step_size = max_step;
			minimal_step_size = min_step;

			absolute_error_tol = tol_abs;
			relative_error_tol = tol_rel;

			if (adaptive)
			{
				for(int i = 0; i < dimension; i++)
				{
					if ((absolute_error_tol < 0.0) || ( relative_error_tol < 0.0))
					{
						std::cerr << "ERROR: Init() Adaptive error tolerances must be positive. Quitting!" << endl;
						exit(-1);
					}

					if ((absolute_error_tol == 0.0) && ( relative_error_tol == 0.0))
					{
						std::cerr << "ERROR: Init() Adaptive error tolerances cannot both be 0.0. Quitting!" << endl;
						exit(-1);

					}
				}
			}

			solver_state.params = pars;

			solver_state.y 		= new double[dimension];
			solver_state.y0 	= new double[dimension];
			solver_state.ytmp 	= new double[dimension];
			solver_state.yerr 	= new double[dimension];
			solver_state.y0_save 	= new double[dimension]; // Save space

			#ifdef DBG_ISAK_REPORT_ODE_SOLVER_STEPS
			std::ostringstream fname("");
			fname << "ode_solver_"<< method_name << "_step_info_" << file_step_output_num << ".dat";
			file_step_output = fopen(fname.str().c_str(),"w+");
			file_step_output_num += 1;
			file_step_output_num = (file_step_output_num % 2); // [0,1]
			#endif
			

			if (method_name == "rkf12") {
				
				selected_step_method_LTE_order = 1.0;
				selected_step_method_max_step_increase = 5.0;
				selected_step_method_storage_size = 3;
				selected_step_method = &ODE_Solver::rkf12_step;
				
				solver_state.y_int 	= new double*[selected_step_method_storage_size];
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					solver_state.y_int[i] 	= new double[dimension];
				}
			} else if (method_name == "rk23"){

				selected_step_method_LTE_order = 3.0;
				selected_step_method_max_step_increase = 5.0;
				selected_step_method_storage_size = 4;
				selected_step_method = &ODE_Solver::rk23_step;

				solver_state.y_int 	= new double*[selected_step_method_storage_size];
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					solver_state.y_int[i] 	= new double[dimension];
				}
				
			
			} else if ( method_name == "rk4"   ) {
				selected_step_method_LTE_order = 4.0; // used with step doubling for error estimate
				selected_step_method_max_step_increase = 5.0;
				selected_step_method_storage_size = 3;
				if (adaptive)
				{
					selected_step_method = &ODE_Solver::rk4_step_err;
				} else {
					selected_step_method = &ODE_Solver::rk4_step;
				}

				solver_state.y_int 	= new double*[selected_step_method_storage_size];
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					solver_state.y_int[i] 	= new double[dimension];
				}


			} else if( method_name == "rkf45" ) {

				selected_step_method_LTE_order = 5.0;
				selected_step_method_max_step_increase = 5.0;
				selected_step_method_storage_size = 6;
				selected_step_method = &ODE_Solver::rkf45_step;

				solver_state.y_int 	= new double*[selected_step_method_storage_size];
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					solver_state.y_int[i] 	= new double[dimension];
				}

			} else if (method_name == "rkdp45")
			{
				selected_step_method_LTE_order = 5.0;
				selected_step_method_max_step_increase = 5.0;
				selected_step_method_storage_size = 7;
				selected_step_method = &ODE_Solver::rkdp45_step;

				solver_state.y_int 	= new double*[selected_step_method_storage_size];
				for(int i = 0; i < selected_step_method_storage_size; i++)
				{
					solver_state.y_int[i] 	= new double[dimension];
				}
				
			} else {
				std::cerr << "ERROR: ode solver name unknown. Quitting!" << endl;
				exit(-1);
			}

			// zero all arrays
			Reset();
		}

		/* Evolve position (t,y(t)) to (t1, y(t1)) with adaptive steps. Initial value and result is stored in y_out
		 * 
		 * Input:
		 *   t_out     -> Initial time, will be updated to t1 once method has finished
		 *   t1        -> Target time, t1 > t, algorithm will run and produce output at t=t1
		 *   step_size -> Requested step size, will be adjusted by function call.
		 *   y_out     -> initial position and output storage at t = t_out
		 */
		void Step_Adaptive(double *t_out, double t1, double *step_size, double *y_out)
		{
			// Some convenient warnings
			if (t1 == *t_out)
			{
				return;
			} else if (( t1 < *t_out) || (*step_size <= 0.0))
			{
				std::cerr << "ERROR: Step_Adaptive cannot integrate backwards (t_final < t_start) or with non-positive timesteps (dt <= 0.0). Quitting." << endl;
				std::cerr << "t_out = " << *t_out << ", t1 = " << t1 << ", step_size = " << *step_size << endl;
				exit(-1);
			}


			// Cannot progess without these
			if (!adaptive)
			{
				std::cerr << "ERROR: Step_Adaptive called, but not initialized as adaptive solver. Qutting!" << endl;
				exit(-1);
			}

			// First iteration copies initialization
			if (num_steps == 0.0)
			{
				prev_err_ratio = -1.0;
				prev_step_size = *step_size;
				solver_state.t  = *t_out;
				solver_state.t0 = *t_out;
				memcpy(solver_state.y, y_out, dimension*sizeof(double));
				memcpy(solver_state.y0_save, y_out, dimension*sizeof(double));
			}
		
			double *y = solver_state.y;
			double *t = &(solver_state.t);

			int status;
			double initial_step_size = *step_size;
			double h0 = *step_size;
			double max_err_ratio;
			while ((num_steps == 0.0)||((t1 - solver_state.t0)/(*t - solver_state.t0) > 1.0))
			{
				// Save initial state in y0
				solver_state.t0 = *t;
				memcpy(solver_state.y0_save, y, dimension*sizeof(double));

				// Compute step with error estimates
				status = (this->*(this->selected_step_method))(*t, *step_size, y);

				// Find new step
				double new_step_size;
				int step_status = adaptive_estimate_error(*step_size, &new_step_size, y, &max_err_ratio);

				// Find appropriate first step size to fit error
				h0 = *step_size;
				int tr_cnt = 0;
				while (step_status == -1)
				{
					h0 = new_step_size;

					// Recover initial state
					memcpy(y, solver_state.y0_save, dimension*sizeof(double));
				
					// Compute step with error estimates
					status = (this->*(this->selected_step_method))(*t, h0, y);

					// Find new step
					step_status = adaptive_estimate_error(h0, &new_step_size, y, &max_err_ratio);

			//		cout << "im stuch here, why? h = " << h0 << ", h_new = " << new_step_size << ", diff = " << new_step_size - h0 << ", step_status = " << step_status << endl;

					if (h0/initial_step_size < 1e-12)
					{
						std::cerr << "ERROR: Step_adaptive is attempting to find a good stepsize to reach your error tolerances, but has reduced the step size to 1e-12*h_start. This could imply a singularity, but progress has stopped. Quitting." << endl;
						exit(-1);
					} else if (h0/initial_step_size < 1e-6)
					{
						double curr_ratio = h0/initial_step_size;
						std::cerr << "WARNING: Step_adaptive has reduced step size to "<< curr_ratio << "*h_start. This could imply a singularity." << endl;
					}

					if (tr_cnt == 100)
					{
						std::cerr << "WARNING: Step_adaptive has used 100 attempts to improve the error of a single step, possible problem with too low error tolerances and/or method has too low order." << endl;

					} else if (tr_cnt == 1000)
					{
						std::cerr << "WARNING: Step_adaptive has used 1000 attempts to improve the error of a single step. Are your error tolerances too low and/or your method too low order!? No more warnings..." << endl;
					}
					tr_cnt += 1;
				}

				// Set next step
				*t = *t + h0;
				*step_size = new_step_size;
	
				// Keep track of previous step size 
				prev_step_size = h0;
				prev_err_ratio = max_err_ratio;

				// Keep track of some details
				num_steps++;
				avg_step_size += (new_step_size-avg_step_size)/num_steps;

				#ifdef DBG_ISAK_REPORT_ODE_SOLVER_STEPS
				fprintf(file_step_output,"%.16e %.16e %.0f\n",*t-h0,h0,num_func_eval);
				#endif
			}

			// Interpolate solution at t_out=t1
			*t_out = t1;
			interpolate_hermite(solver_state.t0, solver_state.t, solver_state.y0_save, solver_state.y, *t_out, y_out);
		}


		/* Evolve position (t,y(t)) to (t1, y(t1)) with fixed steps. Initial value and result is stored in y_out
		 * 
		 * Input:
		 *   t_out     -> Initial time, will be updated to t1 once method has finished
		 *   t1        -> Target time, t1 > t, algorithm will step until t=t1 with requested step size
		 *   step_size -> Requested step size, will not be changed by this method. Can be longer than t1-*t
		 *   y_out     -> initial position and output storage
		 */
		void Step_Fixed(double *t_out, double t1, double *step_size, double *y_out) 
		{
			
			if (t1 == *t_out)
			{
				return;
			} else if (( t1 < *t_out) || (*step_size <= 0.0))
			{
				std::cerr << "ERROR: Step_Fixed cannot integrate backwards (t_final < t_start) or with non-positive timesteps (dt <= 0.0). Quitting." << endl;
				exit(-1);
			}

			// First iteration copies initialization
			if (num_steps == 0.0)
			{
				solver_state.t  = *t_out;
				solver_state.t0 = *t_out;
				memcpy(solver_state.y, y_out, dimension*sizeof(double));
				memcpy(solver_state.y0_save, y_out, dimension*sizeof(double));
			}
		
			double *y = solver_state.y;
			double *t = &(solver_state.t);
	
			int status;
			while ((num_steps == 0.0)||((t1 - solver_state.t0)/(*step_size) > 1.0))
			{
				if (t1 <= *t + *step_size)
				{
					// Save initial state in y0_save
					solver_state.t0 = *t;
					memcpy(solver_state.y0_save, y, dimension*sizeof(double));
				}
				
				status = (this->*(this->selected_step_method))(*t, *step_size, y);
				*t += *step_size;

				// Keep track of some details
				num_steps++;
				avg_step_size += (*step_size-avg_step_size)/num_steps;

				#ifdef DBG_ISAK_REPORT_ODE_SOLVER_STEPS
				fprintf(file_step_output,"%.16e %.16e %.0f\n",*t-*step_size,*step_size,num_func_eval);
				#endif
			}

			// Interpolate solution at t_out=t1
			*t_out = t1;
			interpolate_hermite(solver_state.t0, solver_state.t, solver_state.y0_save, solver_state.y, *t_out, y_out);

		}

		// Return number of function evaluations taken by algorithm
		double number_function_evaluations()
		{
			return num_func_eval;
		}

		// Return average step size used for calculations
		double average_step_size()
		{
			return avg_step_size;
		}

		// Force output file finish writing and close
		void file_output_close()
		{
			#ifdef DBG_ISAK_REPORT_ODE_SOLVER_STEPS
			fprintf(file_step_output,"Total number of steps = %.0f\n",num_steps);
			fprintf(file_step_output,"Avg. step size steps  = %.3e\n",avg_step_size);
			fprintf(file_step_output,"number of min-steps   = %.0f\n",num_steps_min);
			fprintf(file_step_output,"number of func. eval  = %.0f\n",num_func_eval);
			fclose(file_step_output);
			file_step_output = NULL;
			#endif
			
		}


	private:
		std::string method_name;
		int dimension;
		double absolute_error_tol;
		double relative_error_tol;
		double num_func_eval; // count number of function evaluations
		double maximal_step_size; // user defined limit
		double minimal_step_size; // user defined limit
	
		bool adaptive;
		double prev_step_size;
		double prev_err_ratio;
		double fsal_time; // First same as last record time
		double *fsal_p; // Pointer to fsal storage
			
		int (*function)(double , const double *, double *, void *);
		int (*jacobian)(double , const double *, double *, double *, void *);
		
		int (ODE_Solver::*selected_step_method)(double, double, double *);
		double selected_step_method_LTE_order;
		double selected_step_method_max_step_increase;
		int selected_step_method_storage_size;

		// Counters to keep track of average and number of adaptive steps
		double num_steps;
		double num_steps_min; // counter for number of steps taken on minimal step_size
		double avg_step_size;
		FILE *file_step_output; // diagnostic output for stepsizes
		int file_step_output_num;


		struct state {
			double *y0;	// intial y value of integration
			double *y;	// Current y-value of integration
			double t0;	// time for y0 value
			double t;	// time for y value 
			double *ytmp; 	// temporary storage for y values
			double **y_int;	// storage for intermediate steps in integration
			void *params; 	// pointer to potential parameter values used in evaluation
			double *yerr;	// Estimate for error in adaptive methods
			double *y0_save;// intial y value of adaptive integration, for when you need to roll-back changes
			double interp_t0; // Stores t-values for interpolation method
			double interp_t1; // Stores t-values for interpolation method
		} solver_state;

		int rkf12_step(double t, double step_size, double *y)
		{
			// Runge-Kutta-Fehlberg coefficients. Zero elements left out
			const double ah[] = { step_size*1.0 / 2.0, step_size*1.0};

			const double c1 = step_size*1.0 / 256.0;
			const double c2 = step_size*255.0 / 256.0;

			// These are the differences of second and first order coefficientsfor error estimation, 0 elements removed
			const double ec[] = {step_size*1.0 / 512.0,-step_size*1.0 / 512.0};

			// get state names
			int status;
			void *params = solver_state.params;
			double *ytmp = solver_state.ytmp;
			double *yerr = solver_state.yerr;
			double *k1 = solver_state.y_int[0];
			double *k2 = solver_state.y_int[1];
			double *k3 = solver_state.y_int[2];
			
			// Save y0 in intial state vector
			//memcpy(solver_state.y0_save, y, dimension*sizeof(double));

			// k1 step
			if ( fsal_time - t == 0.0)
			{
				if (fsal_p != solver_state.y_int[0])
				{
					k1 = solver_state.y_int[2];
					k3 = solver_state.y_int[0];
				}
				
			} else {
				status = function(t, y, k1, params); num_func_eval++;
			}
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + ah[0]*k1[i];
			}

			// k2 step
			status = function(t + ah[0], ytmp, k2, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += c1*k1[i] + c2*k2[i];
			}

			// k3 step
			status = function(t + ah[1], y, k3, params); num_func_eval++;
			fsal_time = t + step_size; // save correct time
			fsal_p = k3; // save pointer

			// Error estimate
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				yerr[i] = ec[0]*k1[i] + ec[1]*k3[i];
			}
			return 0;
		}

		int rk23_step(double t, double step_size, double *y)
		{
			// Runge-Kutta Bogacki–Shampine coefficients. Zero elements left out
			const double ah[] = { step_size*1.0 / 2.0, step_size*3.0/4.0, step_size*1.0};
			const double b3[] = { step_size*3.0 / 4.0 };

			const double c1 = step_size*2.0 / 9.0;
			const double c2 = step_size*1.0 / 3.0;
			const double c3 = step_size*4.0 / 9.0;

			// These are the differences of second and first order coefficientsfor error estimation, 0 elements removed
			const double ec[] = {-step_size*5.0 / 72.0, step_size*1.0 / 12.0, step_size*1.0 / 9.0, -step_size*1.0/8.0};

			// get state names
			int status;
			void *params = solver_state.params;
			double *ytmp = solver_state.ytmp;
			double *yerr = solver_state.yerr;
			double *k1 = solver_state.y_int[0];
			double *k2 = solver_state.y_int[1];
			double *k3 = solver_state.y_int[2];
			double *k4 = solver_state.y_int[3];
			
			// Save y0 in intial state vector
			//memcpy(solver_state.y0_save, y, dimension*sizeof(double));

			// k1 step
			if ( fsal_time - t == 0.0)
			{
				if (fsal_p != solver_state.y_int[0])
				{
					k1 = solver_state.y_int[3];
					k4 = solver_state.y_int[0];
				}
				
			} else {
				status = function(t, y, k1, params); num_func_eval++;
			}
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + ah[0]*k1[i];
			}

			// k2 step
			status = function(t + ah[0], ytmp, k2, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b3[0]*k2[i]);
			}

			// k3 step
			status = function(t + ah[1], ytmp, k3, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += c1*k1[i] + c2*k2[i] + c3*k3[i];
			}

			status = function(t + ah[2], y, k4, params); num_func_eval++;
			fsal_time = t + step_size; // save correct time
			fsal_p = k4; // save pointer

			// Error estimate
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				yerr[i] = ec[0]*k1[i] + ec[1]*k2[i] + ec[2]*k3[i] + ec[3]*k4[i];
			}
			return 0;
		}

		// Single step of RK4 algorithm
		int rk4_step(double t, double step_size, double *y)
		{
			int status;
	
			// Save y0 in intial state vector
			memcpy(solver_state.y0, y, dimension*sizeof(double));

			// get state names
			const double *y0 = solver_state.y0;
			double *ytmp = solver_state.ytmp;
			double *k = solver_state.y_int[0];
			void *params = solver_state.params;


			// k1 step
			status = function(t, y0, k, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += (step_size/6.0)*k[i];
				ytmp[i] = y0[i] + 0.5*step_size*k[i];
			}
			

			// k2 step
			status = function(t+0.5*step_size, ytmp, k, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += (step_size/3.0)*k[i];
				ytmp[i] = y0[i] + 0.5*step_size*k[i];
			}

			// k3 step
			status = function(t+0.5*step_size, ytmp, k, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += (step_size/3.0)*k[i];
				ytmp[i] = y0[i] + 0.5*step_size*k[i];
			}

			// k4 step
			status = function(t+step_size, ytmp, k, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += (step_size/6.0)*k[i];
			}

			return 0;
		}

		// RK4 step with substeps in order to estimate error
		// ODEINT algorithm from "Numerical Recipes in Pascal" by W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling, Cambridge U. Press, Cambridge 1989 (with a few minor changes)
		int rk4_step_err(double t, double step_size, double *y)
		{
			int status;
	
			// Save y0 in intial state vector
			//memcpy(solver_state.y0_save, y, dimension*sizeof(double));

			// get state names
			double *k = solver_state.y_int[0];
			double *k1 = solver_state.y_int[1];
			double *y0_onestep = solver_state.y_int[2];
			void *params = solver_state.params;
			double *yerr = solver_state.yerr;


			/* Error estimating based on step doubling */

			// take one big step, saved to y0_onestep
			memcpy(y0_onestep, y, dimension*sizeof(double));
			status = rk4_step(t, step_size, y0_onestep);

			// Take two small steps
			status = rk4_step(t                , step_size/2.0, y);
			status = rk4_step(t + step_size/2.0, step_size/2.0, y);

			// error estimating
			for(int i = 0; i < dimension; i++)
			{
				yerr[i] = solver3_ERR_SAFETY * (y[i] - y0_onestep[i])/15.0;
			}

			return 0;
		}

		int rkf45_step(double t, double step_size, double *y)
		{
			// Runge-Kutta-Fehlberg coefficients. Zero elements left out
			const double ah[] = { step_size*1.0 / 4.0, step_size*3.0 / 8.0, step_size*12.0 / 13.0, step_size*1.0, step_size*1.0 / 2.0 };
			const double b3[] = { step_size*3.0 / 32.0, step_size*9.0 / 32.0 };
			const double b4[] = { step_size*1932.0 / 2197.0, -step_size*7200.0 / 2197.0, step_size*7296.0 / 2197.0 };
			const double b5[] = { step_size*8341.0 / 4104.0, -step_size*32832.0 / 4104.0, step_size*29440.0 / 4104.0, -step_size*845.0 / 4104.0 };
			const double b6[] = { -step_size*6080.0 / 20520.0, step_size*41040.0 / 20520.0, -step_size*28352.0 / 20520.0,step_size*9295.0 / 20520.0, -step_size*5643.0 / 20520.0};

			const double c1 = step_size*902880.0 / 7618050.0;
			const double c3 = step_size*3953664.0 / 7618050.0;
			const double c4 = step_size*3855735.0 / 7618050.0;
			const double c5 = -step_size*1371249.0 / 7618050.0;
			const double c6 = step_size*277020.0 / 7618050.0;

			// These are the differences of fifth and fourth order coefficientsfor error estimation, 0 elements removed
			const double ec[] = {step_size*1.0 / 360.0,-step_size*128.0 / 4275.0,-step_size*2197.0 / 75240.0,step_size*1.0 / 50.0,step_size*2.0 / 55.0};

			// get state names
			int status;
			void *params = solver_state.params;
			double *ytmp = solver_state.ytmp;
			double *yerr = solver_state.yerr;
			double *k1 = solver_state.y_int[0];
			double *k2 = solver_state.y_int[1];
			double *k3 = solver_state.y_int[2];
			double *k4 = solver_state.y_int[3];
			double *k5 = solver_state.y_int[4];
			double *k6 = solver_state.y_int[5];
			
			// Save y0 in intial state vector
			//memcpy(solver_state.y0_save, y, dimension*sizeof(double));

			// k1 step
			status = function(t, y, k1, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + ah[0]*k1[i];
			}

			// k2 step
			status = function(t + ah[0], ytmp, k2, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b3[0]*k1[i] + b3[1]*k2[i]);
			}

			// k3 step
			status = function(t + ah[1], ytmp, k3, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b4[0]*k1[i] + b4[1]*k2[i] + b4[2]*k3[i]);
			}

			// k4 step
			status = function(t + ah[2], ytmp, k4, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b5[0]*k1[i] + b5[1]*k2[i] + b5[2]*k3[i] + b5[3]*k4[i]);
			}

			// k5 step
			status = function(t + ah[3], ytmp, k5, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b6[0]*k1[i] + b6[1]*k2[i] + b6[2]*k3[i] + b6[3]*k4[i] + b6[4]*k5[i]);
			}

			// k6 step
			status = function(t + ah[4], ytmp, k6, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				double d_i = c1*k1[i] + c3*k3[i] + c4*k4[i] + c5*k5[i] + c6*k6[i];
				y[i] += d_i;
			}

			// Error estimate
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				yerr[i] = (ec[0]*k1[i] + ec[1]*k3[i] + ec[2]*k4[i] + ec[3]*k5[i] + ec[4]*k6[i]);
			}
			return 0;
		}

		int rkdp45_step(double t, double step_size, double *y)
		{
	//		GSL parameters
			// Runge-Kutta-Fehlberg coefficients. Zero elements left out
			const double ah[] ={ step_size*1.0 / 5.0, step_size*3.0 / 10.0, step_size*4.0 / 5.0, step_size*8.0 / 9.0, step_size*1.0 };
			const double b3[] = { step_size*3.0 / 40.0, step_size*9.0 / 40.0 };
			const double b4[] = { step_size*44.0 / 45.0, -step_size*56.0 / 15.0, step_size*32.0 / 9.0 };
			const double b5[] = { step_size*19372.0 / 6561.0, -step_size*25360.0 / 2187.0, step_size*64448.0 / 6561.0, -step_size*212.0 / 729.0 };
			const double b6[] = { step_size*9017.0 / 3168.0, -step_size*355.0 / 33.0, step_size*46732.0 / 5247.0, step_size*49.0 / 176.0, -step_size*5103.0 / 18656.0};
			const double b7[] = { step_size*35.0 / 384.0 , step_size*500.0 / 1113.0, step_size*125.0 / 192.0, -step_size*2187.0 / 6784.0, step_size*11.0/84.0};

			// 5th order method, remember to include in calculation below!
			const double c1 = step_size*35.0 / 384.0;
			const double c3 = step_size*500.0 / 1113.0;
			const double c4 = step_size*125.0 / 192.0;
			const double c5 = -step_size*2187.0 / 6784.0;
			const double c6 = step_size*11.0/84.0;

			// These are the differences of fifth and fourth order coefficientsfor error estimation
			const double ec[] = {step_size*71.0 / 57600.0,-step_size*71.0 / 16695.0, step_size*71.0 / 1920.0,-step_size*17253.0 / 339200.0, step_size*22.0 / 525.0, -step_size*1.0/40.0};
			
			// get state names
			int status;
			void *params = solver_state.params;
			double *ytmp = solver_state.ytmp;
			double *yerr = solver_state.yerr;
			double *k1 = solver_state.y_int[0];
			double *k2 = solver_state.y_int[1];
			double *k3 = solver_state.y_int[2];
			double *k4 = solver_state.y_int[3];
			double *k5 = solver_state.y_int[4];
			double *k6 = solver_state.y_int[5];
			double *k7 = solver_state.y_int[6];
			
			// Save y0 in intial state vector
			//memcpy(solver_state.y0_save, y, dimension*sizeof(double));


			// k1 step
			if ( fsal_time - t == 0.0)
			{
				if (fsal_p != solver_state.y_int[0])
				{
					k1 = solver_state.y_int[6];
					k7 = solver_state.y_int[0];
				}
				
			} else {
				status = function(t, y, k1, params); num_func_eval++;
			}
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + ah[0]*k1[i];
			}

			// k2 step
			status = function(t + ah[0], ytmp, k2, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b3[0]*k1[i] + b3[1]*k2[i]);
			}

			// k3 step
			status = function(t + ah[1], ytmp, k3, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b4[0]*k1[i] + b4[1]*k2[i] + b4[2]*k3[i]);
			}


			// k4 step
			status = function(t + ah[2], ytmp, k4, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b5[0]*k1[i] + b5[1]*k2[i] + b5[2]*k3[i] + b5[3]*k4[i]);
			}

			// k5 step
			status = function(t + ah[3], ytmp, k5, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				ytmp[i] = y[i] + (b6[0]*k1[i] + b6[1]*k2[i] + b6[2]*k3[i] + b6[3]*k4[i] + b6[4]*k5[i]);
			}

			// k6 step and sum
			status = function(t + step_size, ytmp, k6, params); num_func_eval++;
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				y[i] += (b7[0]*k1[i] + b7[1]*k3[i] + b7[2]*k4[i] + b7[3]*k5[i] + b7[4]*k6[i]);
			}
			

			// k7 step for error estimation
			status = function(t + step_size, y, k7, params); num_func_eval++;
			fsal_time = t + step_size; // save correct time
			fsal_p = k7; // save pointer


			// Error estimate
			#pragma ivdep
			for(int i = 0; i < dimension; i++)
			{
				yerr[i] = (ec[0]*k1[i] + ec[1]*k3[i] + ec[2]*k4[i] + ec[3]*k5[i] + ec[4]*k6[i] + ec[5]*k7[i]);
			}
			return 0;
		}


		int adaptive_estimate_error(double step_size, double *new_step_size, double *y, double *max_err_ratio)
		{
			double q = selected_step_method_LTE_order; // local truncation error order of formula
			double max_h = selected_step_method_max_step_increase; // maximal increase factor

			*max_err_ratio = 0.0;
			for(int i = 0; i < dimension; i++)
			{
				// compute desired error level
				const double alpha = 1.0; // default error scaling
				double adaptive_err = absolute_error_tol + relative_error_tol*alpha*fabs(y[i]);

				double err_ratio = fabs(solver_state.yerr[i])/(adaptive_err);

				// Find maximal error ratio in array
				if (err_ratio > *max_err_ratio)
				{
					*max_err_ratio = err_ratio;
				}
			}

			// Initialize at current step
			*new_step_size = step_size;

			//==== Controller ====
			double fact = 1.0;
			fact = pow(*max_err_ratio,-1.0/q);
/*
			if (prev_err_ratio < 0.0)
			{
				// Use elementary controller [3]
				fact = pow(*max_err_ratio,-1.0/q);

			} else {
				// Use PI.4.2 controller [G.Soderlind. Autometatic control and adaptive time-stepping Numerical Algorithms 31:281-310 2002]
	//			fact = pow(*max_err_ratio,-3.0/(5.0*))*pow(prev_err_ratio,1.0/(5.0*q));
			
				// H211b digital filter (with b=4) [G. Soderlind. Digital filters in adaptive time-stepping ACM Tans. Math. Software 29:1-26 2003]
				fact = pow(*max_err_ratio,-1.0/(4.0*q))*pow(prev_err_ratio,-1.0/(4.0*q))*pow(step_size/prev_step_size, -1.0/4.0);
	//			cout  << "fact = " << fact << ", *max_err_ratio = " << *max_err_ratio << ", prev_err_ratio = " << prev_err_ratio << ", prev_step_size" << prev_step_size << endl;
			}
*/			
			// Limits
			if (fact < 0.1)
			{
				fact = 0.1;
			} else if (fact > max_h)
			{
				fact = max_h;
			}
			*new_step_size *= fact;


			// Maintain maximal step size
			if (*new_step_size > maximal_step_size)
			{
				*new_step_size = maximal_step_size;
			} else if (*new_step_size < minimal_step_size) 
			{
				*new_step_size = minimal_step_size;
				num_steps_min++;
				return 1; // Dont try to make this better
			}


			// Increase or reduce step size? 
			if (*max_err_ratio > 1.1) // If observed error is > 10% of error level then reduce step
			{
				double safe = 0.9; // safety factor for step size reduction (0,1)
				*new_step_size *= safe; // aditional reduction

				return (-1);

			} else if (*max_err_ratio < 0.5) // If observed error is < 50% of error level then increase step
			{

				return 1;
			}

			// Keep step size the same
			*new_step_size = step_size;
			return 0;
		}


		/* Interpolate data from points (t0,y0) and (t1,y1) using Hermite splines.
		 * The derivatives (t0,f0) and (t1,f1) are collected from the ODE solvers
		 * and might require additional function evaluations (depending on ODE solver used).
		 * t_out in  [t0,t1] is the target time
		 * y_out is where to store 
		 **/
		void interpolate_hermite(double t0, double t1, double *y0, double *y1, double t_out, double *y_out)
		{
			// Interpolate solution at t_out=t1
			// Hermite Spline interpolant
			double dt = t1-t0;
			double theta = (t_out - t0)/dt; // in [0,1]

			double a0 = (1.0-theta)*(1.0-theta)*(1.0 +2.0*theta);
			double a1 = theta*theta*(3.0-2.0*theta);
			double b0 = theta*(1.0-theta)*(1.0-theta)*dt;
			double b1 = theta*theta*(theta-1.0)*dt;

			double *f0, *f1;
			if (method_name == "rkf12") {
				f0 = solver_state.y_int[0];
				f1 = solver_state.y_int[2];
				if (fsal_p != solver_state.y_int[2])
				{
					f0 = solver_state.y_int[2];
					f1 = solver_state.y_int[0];
				}
			} else if(method_name == "rk23"){
				f0 = solver_state.y_int[0];
				f1 = solver_state.y_int[3];
				if (fsal_p != solver_state.y_int[3])
				{
					f0 = solver_state.y_int[3];
					f1 = solver_state.y_int[0];
				}
				
			} else if( method_name == "rk4") 
			{
				if (solver_state.interp_t0 != t0)
				{
					void *params = solver_state.params;
					int status = function(t0, y0, solver_state.y_int[0], params); num_func_eval++;
					solver_state.interp_t0 = t0;
				}
				if (solver_state.interp_t1 != t1)
				{
					void *params = solver_state.params;
					int status = function(t1, y1, solver_state.y_int[1], params); num_func_eval++;
					solver_state.interp_t1 = t1;
				}

				f0 = solver_state.y_int[0];
				f1 = solver_state.y_int[1];

			} else if (method_name == "rkf45")
			{
				if (solver_state.interp_t1 != t1)
				{
					void *params = solver_state.params;
					int status = function(t1, y1, solver_state.y_int[5], params); num_func_eval++;
					solver_state.interp_t1 = t1;
				}

				f0 = solver_state.y_int[0];
				f1 = solver_state.y_int[5];

			} else if (method_name == "rkdp45")
			{
				
				f0 = solver_state.y_int[0];
				f1 = solver_state.y_int[6];
				if (fsal_p != solver_state.y_int[6])
				{
					f0 = solver_state.y_int[6];
					f1 = solver_state.y_int[0];
				}
			}

			// update output y(t_out)
			for(int i = 0; i < dimension; i++)
			{
				y_out[i] = a0*y0[i] + a1*y1[i] + b0*f0[i] + b1*f1[i];
			}

		}
};




#endif // __solver3_h__
































