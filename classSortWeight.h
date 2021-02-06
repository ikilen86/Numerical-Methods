/*
	Sort by Weight into array
	* 
	Isak Kilen @ 2016
*/

#ifndef __SORTWEIGHT_H_INCLUDED__
#define __SORTWEIGHT_H_INCLUDED__


#include <iostream>

using namespace std;


class sortWeight
{
	public:
		sortWeight()
		{
			a_1 = NULL;
			index = NULL;
			array_length = 0;
		}
		
		sortWeight(double *A_1, int N)
		{
			a_1 = A_1;
			index = NULL;
			array_length = N;
			

			index = new int[array_length];
			
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}
		
		void initSort(double *A_1,int N)
		{
			array_length = N;
			
			a_1 = A_1;
			
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
		
		int *runSort()
		{
			quick_sort(index, array_length);
			return index;
		}
		
		// Apply sorting to array[] via a temp array
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
	
		// is |F[i]| < |F[j]| ?
		bool my3Sort(int i, int j) 
		{
			return (abs(a_1[i]) > abs(a_1[j]));
		}
		
		
		//========================
		// Basic quick sort algo.
		int partition(int list[], int p, int r)
		{
			int index, exchange_temp;
			int pivot = list[r];
			index = p - 1;
			for(int i = p; i < r; i++)
			{
				//if(list[i] <= pivot)
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
