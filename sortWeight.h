/*
	Sort weight indices using quicksort.
	Isak Kilen
*/

#ifndef __SORTWEIGHT_H_INCLUDED__
#define __SORTWEIGHT_H_INCLUDED__


#include <iostream>

using namespace std;


class sortWeight
{
	public:
		//! Am empty constructor
		sortWeight()
		{
			a_1 = NULL;
			index = NULL;
			array_length = 0;
		}
		
		//! A constructor
		/*! \param A_1 Weights that will be sorted
		 *  \param N Number of elements in A_1
		 **/
		sortWeight(double *A_1, int N)
		{
			a_1 = A_1;
			index = NULL;
			array_length = N;
			
			// Initialize index array
			index = new int[array_length];
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}
		
		//! Re-initialize sorting, can be used to switch the input weights or length
		/*! \param A_1 Weights that will be sorted
		 *  \param N Number of elements in A_1
		 **/
		void initSort(double *A_1,int N)
		{
			array_length = N;
			a_1 = A_1;
			
			// Initialize index array
			if(index==NULL)
			{
				index = new int[array_length];
			} else {
				
				delete [] index;
				index = new int[array_length];
			}
			
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}
		
		//! Compute index I for sorted weights. Does not alter original weights.
		//! \return pointer sorted index array I. such that A(I) = B is sorted
		//! \sa initSort 
		int *runSort()
		{
			quick_sort(index, array_length);
			return index;
		}
		
		//! Apply sorting to a (int) array using an external temp array.
		/*! Sort array A, after running runSort to compute index set.
		 * \param array target array to sort
		 * \param dummy an array to assist with sorting.
		 * \se runSort
		 **/
		void applyIndexSort_int(int *array, int *dummy)
		{
			for(int i = 0; i < array_length; i++)
			{
				dummy[i] = array[index[i]];
			}
			
			for(int i = 0; i < array_length; i++)
			{
				array[i] = dummy[i];
			}
		}
		
		//! Apply sorting to a (double) array using an external temp array.
		/*! Sort array A, after running runSort to compute index set.
		 * \param array target array to sort
		 * \param dummy an array to assist with sorting.
		 * \se runSort
		 **/
		void applyIndexSort_double(double *array, double *dummy)
		{
			for(int i = 0; i < array_length; i++)
			{
				dummy[i] = array[index[i]];
			}
			
			for(int i = 0; i < array_length; i++)
			{
				array[i] = dummy[i];
			}
		}
		
		
	
	private:
	
		//! User defined comparison function for sorting: is f(a) < f(b)
		//! Can be switched to finding 'f(a) > f(b)' or more complex expressions
		bool my3Sort(int i, int j) 
		{
			return (abs(a_1[i]) < abs(a_1[j]));
		}
		
		
		//! Basic quick sort partition function.
		int partition(int list[], int p, int r)
		{
			int index, exchange_temp;
			int pivot = list[r];
			index = p - 1;
			for(int i = p; i < r; i++)
			{
				if (my3Sort(list[i],pivot))
				{
					index++;
					exchange_temp = list[i];
					list[i] = list[index];
					list[index] = exchange_temp;
				}
			}
			exchange_temp = list[r];
			list[r] = list[index+1];
			list[index+1] = exchange_temp;
			return index+1;
		}

		//! Quicksort helper function
		void quicksort_aux(int list[], int p, int r)
		{
			int q;
			if(p<r)
			{
				q = partition(list, p, r);
				quicksort_aux(list, p, q-1);
				quicksort_aux(list, q+1, r);
			}
		}

		//! Quick sort main function
		void quick_sort(int list[], int size)
		{
			quicksort_aux(list,0, size-1);
		}
		//======================
	
		double *a_1;
		int *index;
		int array_length;
		
};

#endif
