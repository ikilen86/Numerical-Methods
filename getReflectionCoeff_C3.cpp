
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <complex>

using namespace std;
std::complex<double> I(0,1);
const double c0 = 2.99792458E+08;

#define MATLAB_RUN
//#define USE_OPENMP

#ifdef MATLAB_RUN
	#include "matrix.h"
	#include "mex.h"
#endif

#ifdef USE_OPENMP
	#include <omp.h>
#endif

/** Find Reflection and tranmission coefficient of structure
	r,t 	-> Where we store reflection and transmission coefficient (complex)
	n	-> List of indices of layers (light incoming from the left to the right 0,1,2,3...)
	h	-> List of layer thickness
	NumLayers	-> Number of layers in structure (not including air/first layer and substrate)
	na		-> Infinite medium to the left, where incoming pulse is located
	nb		-> Infinite medium to the right
	M_end		-> If parts of transfer matrix is precalculated, then this is on the right, (otherwise leave as NULL)
*/


void getRT(std::complex<double> *r, std::complex<double> *t,double *n_re, double *n_im, double *h, int NumLayers, double freq, double na, double nb, std::complex<double> *G0, std::complex<double> *G_end)
{
	//======================
	// Find transfer matrix
	//======================
	double tmp = freq/c0;

	std::complex<double> ni, nim;
	std::complex<double> rho;
	std::complex<double> G;
	// Add on any precalculated layers
	if (G_end != NULL)
	{
		// Initialize M0 matrix
		G = G_end[0];
	} else {
		
		ni 	= n_re[NumLayers-1] + I*n_im[NumLayers-1];
		rho = (ni-nb)/(nb+ni);
		G = rho;
	}

	std::complex<double> psi, e_psi;
	for(int i = NumLayers-1; i >= 1 ;i--)
	{
		ni 		= n_re[i] 	+ I*n_im[i];
		nim 	= n_re[i-1] + I*n_im[i-1];
		rho 	= (nim - ni)/(nim + ni);
		psi 	= I*tmp*ni*h[i];
		e_psi 	= exp(-2.0*psi);
		
		G = (rho + G*e_psi)/(1.0+rho* G*e_psi);
	}
	if (G0 != NULL)
	{
		*G0 = G;
	}
	
	
	// Do not include the final transition into G0
	ni 	= n_re[0] 	+ I*n_im[0];
	rho 	= (na - ni)/(na + ni);
	psi 	= I*tmp*ni*h[0];
	e_psi 	= exp(-2.0*psi);
	
	*r = (rho+G*e_psi)/(1.0+rho*G*e_psi);
	*t = 0;
}

#ifdef MATLAB_RUN
/**
 * Usage: 	[r,t,G0] = (na,nb,freq,n,h,G_init,n_qw_ind,n_qw)
 *			[r,t,G0] = (na,nb,freq,n,h,[],n_qw_ind,n_qw)
 *			[r,t,G0] = (na,nb,freq,n,h,G_init)
 * 			[r,t,G0] = (na,nb,freq,n,h)
 * 
 *  Desc: 	Light coming from material 'a' on the left propagating into
 *  		material 'b' on the right. Each element of the lists 'n' and 'h' 
 * 			represent the refractive index and length of a material
 * 			between material 'a' and 'b'.
 * 
 * Input:
 * na,nb 	-> Refractive indices of material 'a' and 'b'
 * freq 	-> Frequency of light, can be a list. Units: [1/s]
 *  
 * n		-> List of refractive indices of material layers
 * h		-> List of lengths of material layers. Units: [m]
 * G		-> A precomputed reflection coefficient, for each freq.
 * 				When this is used, the final layer index and width has to be
 * 				included in and the final element in n and h.
 * 				ex: First compute G using n = [n1,n2,n3,n4];
 * 					Then want to use G while adding on a single layer on the left
 * 					The new n vector is n=[n_new, n1]; and similar for h
 * n_qw_ind -> List of indices of which layers of n to apply n_qw to.
 * n_qw		-> List of frequency dependent refractive indices. NumQW x NumFreq matrix
 * 
 * 
 * Output:
 * r,t		-> Amplitude Reflection and Transmission coefficients (complex).
 * 				Do R = abs(r)^2 to get intensity reflection.
 * 				Do T = (nb/na)*abs(t)^2 to get intensity transmission.
 * G0		-> The computed reflection coefficient for each frequency
 * */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!(((nrhs == 5)||(nrhs ==6))||(nrhs ==8)))
	{
		mexPrintf("nrhs = %d\n",nrhs);
		mexErrMsgTxt("getReflectionCoefficient(): Need 5,6 or 8 input arguments.");
	}
	
	if (!(((nlhs == 1)||(nlhs == 2))||(nlhs == 3)))
	{
		mexErrMsgTxt("getReflectionCoefficient(): Need 1, 2 or 3 output arguments");
	}
	
	
	int NumLayers;
	double na 			= *(double*)mxGetData(prhs[0]);
	double nb 			= *(double*)mxGetData(prhs[1]);
	
	int NumFreqs;
	int freq_rows,freq_cols;
	//===============
	// input: freq
	//===============
	int number_of_dims = mxGetNumberOfDimensions(prhs[2]);
	if (number_of_dims > 2)
	{
		mexErrMsgTxt("getReflectionCoefficient():freq has to be a 1D list");
	}
	const int  *dim_array = mxGetDimensions(prhs[2]);
	if (dim_array[0]==1)
	{
		freq_rows = 1;
		freq_cols = dim_array[1];
		NumFreqs = dim_array[1];
	} else if (dim_array[1]==1)
	{
		freq_rows = dim_array[0];
		freq_cols = 1;
		NumFreqs = dim_array[0];
	} else {
		mexErrMsgTxt("getReflectionCoefficient():freq has to be a 1D list");
	}
	
	double freq[NumFreqs];
	if (mxGetPr(prhs[2]) != NULL)
	{
		for(int i =0; i < NumFreqs; i++) 
		{
			freq[i] = mxGetPr(prhs[2])[i];
		}
	} else {
		mexErrMsgTxt("getReflectionCoefficient():freq must be real");
	}
	
	//============
	// input: n
	//============
	number_of_dims = mxGetNumberOfDimensions(prhs[3]);
	if (number_of_dims > 2)
	{
		mexErrMsgTxt("getReflectionCoefficient():n has to be a 1D list");
	}
	dim_array = mxGetDimensions(prhs[3]);
	if ((dim_array[0]==0)&&(dim_array[1]==0))
	{
		NumLayers = 0;
		mexErrMsgTxt("getReflectionCoefficient():n has to contain at least 1 layer");
	} else if (dim_array[0]==1)
	{
		NumLayers = dim_array[1];
	} else if (dim_array[1]==1)
	{
		NumLayers = dim_array[0];
	} else {
		mexErrMsgTxt("getReflectionCoefficient():n has to be a 1D list");
	}
	
	double n_re[NumLayers];
	
	if (mxGetPr(prhs[3]) != NULL)
	{
		for(int i =0; i < NumLayers; i++) 
		{
			n_re[i] = mxGetPr(prhs[3])[i];
		}
	} else {
		for(int i =0; i < NumLayers; i++) 
		{
			n_re[i] = 0;
		}
	}
	
	double n_im[NumLayers];
	if (mxGetPi(prhs[3]) != NULL)
	{
		for(int i =0; i < NumLayers; i++) 
		{
			n_im[i] = mxGetPi(prhs[3])[i];
		}
	} else {
		for(int i =0; i < NumLayers; i++) 
		{
			n_im[i] = 0.0;
		}
	}
	
	//==========
	// input: h
	//==========
	int h_number_of_dims = mxGetNumberOfDimensions(prhs[4]);
	if (h_number_of_dims > 2)
	{
		mexErrMsgTxt("getReflectionCoefficient():h has to be a 1D list");
	}
	dim_array = mxGetDimensions(prhs[4]);
	if ((dim_array[0]==0)&&(dim_array[1]==0))
	{
		if (0 != NumLayers)
		{
			mexErrMsgTxt("getReflectionCoefficient():n and h has to be of same length");
		}
		//NumLayers = 0;
	} else if (dim_array[0]==1)
	{
		//NumLayers = dim_array[1];
		if (dim_array[1] != NumLayers)
		{
			mexErrMsgTxt("getReflectionCoefficient():n and h has to be of same length");
		}
	} else if (dim_array[1]==1)
	{
		if (dim_array[0] != NumLayers)
		{
			mexErrMsgTxt("getReflectionCoefficient():n and h has to be of same length");
		}
		//NumLayers = dim_array[0];
	} else {
		mexErrMsgTxt("getReflectionCoefficient():h has to be a 1D list");
	}
	
	double h[NumLayers];
	for(int i =0; i < NumLayers; i++) 
	{
		h[i] = mxGetPr(prhs[4])[i];
	}

	//==========================
	// input: G_final
	//==========================
	std::complex<double> *G_final = NULL;
	if (nrhs >= 6)
	{
		// Check input G_final matrix
		int G_number_of_dims = mxGetNumberOfDimensions(prhs[5]);
		if (G_number_of_dims != 2)
		{
			mexErrMsgTxt("getReflectionCoefficient(): G_final has to be a 1D list");
		}
		dim_array = mxGetDimensions(prhs[5]);
		
		if (!((dim_array[0]==0)&&(dim_array[1]==0)))
		{
			if (dim_array[0]==1)
			{
				//NumLayers = dim_array[1];
				if (dim_array[1] != NumFreqs)
				{
					mexErrMsgTxt("getReflectionCoefficient():G_final and h has to be of same length as freq");
				}
			} else if (dim_array[1]==1)
			{
				if (dim_array[0] != NumFreqs)
				{
					mexErrMsgTxt("getReflectionCoefficient():G_final and h has to be of same length as freq");
				}
				//NumLayers = dim_array[0];
			} else {
				mexErrMsgTxt("getReflectionCoefficient():G_final and h has to be a 1D list");
			}
			
			G_final = new std::complex<double>[NumFreqs];
		
			for(int i =0; i < NumFreqs; i++)
			{
				G_final[i] 	= mxGetPr(prhs[5])[i];
			}
			
			if (mxGetPi(prhs[5]) != NULL)
			{
				for(int i =0; i < NumFreqs; i++)
				{
					G_final[i] 	+= I*mxGetPi(prhs[5])[i];
				}
			}
			
		} else {
			//mexErrMsgTxt("getReflectionCoefficient():G_final empty??");
		}
	}
	
	//=================
	// input: n_qw_ind
	//=================
	int NumQWs = 0;
	int *n_qw_ind = NULL;
	if (nrhs >= 8)
	{
		number_of_dims = mxGetNumberOfDimensions(prhs[6]);
		if (h_number_of_dims > 2)
		{
			mexErrMsgTxt("getReflectionCoefficient():n_qw_ind has to be a 1D list");
		}
		dim_array = mxGetDimensions(prhs[6]);
		if ((dim_array[0]==0)&&(dim_array[1]==0))
		{
			NumQWs = 0;
		} else if (dim_array[0]==1)
		{
			NumQWs = dim_array[1];
		} else if (dim_array[1]==1)
		{
			NumQWs = dim_array[0];
		} else {
			mexErrMsgTxt("getReflectionCoefficient():n_qw_ind has to be a 1D list");
		}
		
		if (NumQWs > 0)
		{
			n_qw_ind = new int[NumQWs];
			for(int i =0; i < NumQWs; i++) 
			{
				n_qw_ind[i] = mxGetPr(prhs[6])[i];
			}
		}
	}
	
	//==============
	// input: n_qw
	//==============
	std::complex<double> **n_qw = NULL;
	if (nrhs >= 8)
	{
		// Check input G_final matrix
		number_of_dims = mxGetNumberOfDimensions(prhs[7]);
		if (number_of_dims != 2)
		{
			mexErrMsgTxt("getReflectionCoefficient(): n_qw has to be a 2D array");
		}
		dim_array = mxGetDimensions(prhs[7]);
		
		if (!((dim_array[0]==0)&&(dim_array[1]==0)))
		{
			n_qw = new std::complex<double> *[NumFreqs];
			if ((dim_array[0]!=NumQWs)||(dim_array[1] != NumFreqs))
			{
				mexErrMsgTxt("getReflectionCoefficient(): n_qw has to have size: Num_QW x Num_Freq");
			}
		
			for(int j =0; j < NumFreqs; j++)
			{
				n_qw[j] = new std::complex<double> [NumQWs];
				for(int i =0; i < NumQWs; i++)
				{
					n_qw[j][i] = mxGetPr(prhs[7])[j*NumQWs + i];// + I*mxGetPi(prhs[7])[j + i*NumFreqs];
				}
				if (mxGetPi(prhs[7]) != NULL)
				{
					for(int i =0; i < NumQWs; i++)
					{
						n_qw[j][i] += I*mxGetPi(prhs[7])[j*NumQWs + i];
					}
				}
			}
		} else {
			if (NumQWs != 0)
			{
				mexErrMsgTxt("getReflectionCoefficient(): n_qw doesnt have any QWs, but n_qw_ind does...");
			}
		}
	}
	
	//===============================
	// output: G0 for each frequency
	//===============================
	std::complex<double> *G0 = NULL;
	if (nlhs >= 3)
	{
		G0 = new std::complex<double>[NumFreqs];
	}
	
	//==============
	// Calculations
	//==============
	#ifdef USE_OPENMP
		omp_set_num_threads(4);
	#endif

	std::complex<double> r[NumFreqs];
	std::complex<double> t[NumFreqs];

	if (nrhs < 6)
	{
		if (nlhs >= 3)
		{
			// No QW indices included
			#ifdef USE_OPENMP
			#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb,G0)
			#endif
			for(int i =0 ;i < NumFreqs; i++)
			{
				getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,&(G0[i]),NULL);
			}
			
		} else {
			
			#ifdef USE_OPENMP
			#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb)
			#endif
			for(int i =0 ;i < NumFreqs; i++)
			{
				getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,NULL,NULL);
			}
		}
		
	} 
	else if (nrhs == 6)
	{
		if (G_final != NULL)
		{
			if (nlhs >= 3)
			{
				// No QW indices included
				#ifdef USE_OPENMP
				#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb,G0,G_final)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,&(G0[i]),&(G_final[i]));
				}
				
			} else {
				
				#ifdef USE_OPENMP
				#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb,G_final)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,NULL,&(G_final[i]));
				}
			}
		} else {
			
			if (nlhs >= 3)
			{
				// No QW indices included
				#ifdef USE_OPENMP
				#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb,G0,G_final)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,&(G0[i]),NULL);
				}
				
			} else {
				
				#ifdef USE_OPENMP
				#pragma omp parallel for shared(r,t,n_re,n_im,h,NumLayers,freq,na,nb,G_final)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,NULL,NULL);
				}
			}
		}
		
	} 
	else {
		if (G_final != NULL)
		{
			if (nlhs >= 3)
			{
				#ifdef USE_OPENMP
				#pragma omp parallel for private(n_re,n_im) shared(r,t,h,NumLayers,freq,na,nb,G0,G_final,n_qw_ind,n_qw)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					for(int j = 0; j < NumQWs; j++)
					{
						n_re[n_qw_ind[j]] = real(n_qw[i][j]);
						n_im[n_qw_ind[j]] = imag(n_qw[i][j]);
					}
					
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,&(G0[i]),&(G_final[i]));
				}
				
			} else {
				
				#ifdef USE_OPENMP
				#pragma omp parallel for private(n_re,n_im) shared(r,t,h,NumLayers,freq,na,nb,G_final,n_qw_ind,n_qw)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					for(int j = 0; j < NumQWs; j++)
					{
						n_re[n_qw_ind[j]] = real(n_qw[i][j]);
						n_im[n_qw_ind[j]] = imag(n_qw[i][j]);
					}
					
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,NULL,&(G_final[i]));
				}
			}
		} else {
			
			if (nlhs >= 3)
			{
				#ifdef USE_OPENMP
				#pragma omp parallel for private(n_re,n_im) shared(r,t,h,NumLayers,freq,na,nb,G0,G_final,n_qw_ind,n_qw)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					for(int j = 0; j < NumQWs; j++)
					{
						n_re[n_qw_ind[j]] = real(n_qw[i][j]);
						n_im[n_qw_ind[j]] = imag(n_qw[i][j]);
					}
					
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,&(G0[i]),NULL);
				}
				
			} else {
				
				#ifdef USE_OPENMP
				#pragma omp parallel for private(n_re,n_im) shared(r,t,h,NumLayers,freq,na,nb,G_final,n_qw_ind,n_qw)
				#endif
				for(int i =0 ;i < NumFreqs; i++)
				{
					for(int j = 0; j < NumQWs; j++)
					{
						n_re[n_qw_ind[j]] = real(n_qw[i][j]);
						n_im[n_qw_ind[j]] = imag(n_qw[i][j]);
					}
					
					getRT(&(r[i]),&(t[i]),n_re,n_im,h,NumLayers,freq[i],na,nb,NULL,NULL);
				}
			}
		}
	}
	

	// Return values
	plhs[0] = mxCreateDoubleMatrix(freq_rows,freq_cols,mxCOMPLEX); // Allocate r
	double *r_re = mxGetPr(plhs[0]);
	double *r_im = mxGetPi(plhs[0]);
	for(int i =0 ;i < NumFreqs; i++)
	{
		 r_re[i] = real(r[i]);
		 r_im[i] = imag(r[i]);
	}
	
	if (nlhs >= 2)
	{
		plhs[1] = mxCreateDoubleMatrix(freq_rows,freq_cols,mxCOMPLEX); // Allocate t
		double *t_re = mxGetPr(plhs[1]); 
		double *t_im = mxGetPi(plhs[1]); 
		for(int i =0 ;i < NumFreqs; i++)
		{
			t_re[i] = real(t[i]);
			t_im[i] = imag(t[i]);
		}
	}
	
	if (nlhs >= 3)
	{
        plhs[2] = mxCreateDoubleMatrix(freq_rows,freq_cols,mxCOMPLEX); // Allocate G0
        
        double *R_re = mxGetPr(plhs[2]); 
		double *R_im = mxGetPi(plhs[2]); 
		for(int i =0 ;i < NumFreqs; i++)
		{
			R_re[i] = real(G0[i]);
			R_im[i] = imag(G0[i]);
		}
	}
	
	//==========
	// clean up
	//==========
	if (nlhs >= 3)
	{
		delete [] G0;
	}
	
	
	if (nrhs >= 6)
	{
		if (G_final != NULL)
		{
			delete [] G_final;
		}
	}
	
	
	if (nrhs >= 8)
	{
		if (n_qw_ind != NULL)
		{
			delete [] n_qw_ind;
		}
		
		if (n_qw != NULL)
		{
			for(int i =0 ;i < NumFreqs; i++)
			{
				delete [] n_qw[i];
			}
			delete [] n_qw;
		}
	}
}

#else
int main()
{
	cout << "Hello world" << endl;
	std::complex<double> r;
	std::complex<double> t;

	int N_layers = 10;
	double lambda = 1000.0e-9;
	double na = 1;
	double nb = 1;


	double *n_re = new double[N_layers];
	double *n_im = new double[N_layers];
	double *h = new double[N_layers];

	for(int i = 0; i < N_layers; i++)
	{
		n_re[i] = 1.0 + (rand() / ((double)RAND_MAX));
		n_im[i] = 0;
		h[i] = 100.0e-9;
	}

	getRT(&r,&t,n_re,n_im,h,N_layers,lambda,na,nb,NULL);
	cout << "r = " << r << endl;
	cout << "t = " << t << endl;
	cout << "|r|^2 + |t|^2 = " << r*conj(r) + t*conj(t) << endl;

	cout << "done.." << endl;

	return 0;
}
#endif


