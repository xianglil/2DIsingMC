//
//  Array3d.hpp
//  3dArray
//
//  Created by Xianglin Liu on 12/11/18.
//  Copyright © 2018 Xianglin Liu. All rights reserved.
//

//	the multi-dimensional arrays defined here adopt the column major ordering (as in Fortran)
//	to maximize the efficiency of the LAPACK routine.
//	Note the index still starts with 0 instead of 1.

#ifndef Array3d_hpp
#define Array3d_hpp

#include <cassert>
#include <iostream>

using std::cout; 
using std::endl;
using std::string;

namespace Array3d
{
	template <class T>
	class Array2d
	{
	public:
		//		Default constructor: Array2d<double> mat; mat=mat_in;
		Array2d(): col_size(0), row_size(0), tot_size(0), data(nullptr){;}
		
		//		Constructor for direct initialization: Array2d<double> mat(2,2);
		Array2d( size_t col_in, size_t row_in ):
		col_size(col_in), row_size(row_in), tot_size(col_in*row_in)
		{
			assert( col_in>=0 && row_in>=0 ); // make sure the lengths >=0
			try
			{
				data = new T[tot_size];
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			
		}
		
		//		Copy constructor: Array2d<double> mat(mat_in) or Array2d<double> mat=mat_in;
		Array2d( const Array2d<T> &a):
		col_size(a.col_size), row_size(a.row_size), tot_size(a.tot_size)
		{
			try
			{
				data = new T[tot_size];
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = a.data[i];
			}
		}
		
		~Array2d()
		{
			delete [] data;
			data = nullptr;
		}
		
		//		Overload operator () to access matrix elements. Note no range-check is performed.
		T& operator()( size_t col, size_t row )
		{
			return data[ col+row*col_size ];
		}
		
		//		Overload operator [] to access matrix as a 1d array.
		T& operator[](size_t i)
		{
			return data[i];
		}
		
		//		Overload copy assignment operator.
		Array2d& operator=( const Array2d<T>& a )
		{
			if (this == &a) return *this;
			if (data)
			{
				delete [] data;
			}
			col_size = a.col_size;
			row_size = a.row_size;
			tot_size = a.tot_size;
			try
			{
				data = new T[tot_size];
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = a.data[i];
			}
			return *this;
		}
		
		//		Overload operator = to assign value to the whole matrix.
		Array2d& operator=( const T& val )
		{
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = val;
			}
			return *this;
		}
		
		//		Methods to access the information about the matrix.
		size_t colsize() const { return col_size; }
		size_t rowsize() const { return row_size; }
		size_t totsize() const { return tot_size; }
		
		//	private:
		size_t col_size;
		size_t row_size;
		size_t tot_size;
		T* data;
		const char* classname() { return "Array2d"; }
		
	};

	template <class T>
	void print_test_Array2d( Array2d<T>& mat)
	{
		cout << "***********************************************" << endl;
		cout << " Test of the Array2d class" << endl;
		cout << "mat.colsize()=" << mat.colsize() << endl;
		cout << "mat.rowsize()=" << mat.rowsize() << endl;
		cout << "mat.totsize()=" << mat.totsize() << endl;
		
		for ( int m=0; m<mat.totsize(); m++ )
		{
			printf( "initialized mat(%i)=%.15f \n", m, mat[m] );
			mat[m] = double(m);
		}
		for ( int j=0; j<mat.rowsize(); j++)
			for ( int i=0; i<mat.colsize(); i++)
			{
				printf( "assigned mat(%i, %i)=%.15f \n", i, j, mat(i,j) );
			}
		cout << "Success run !\n";
		cout << "***********************************************" << endl;
	}


	template <class T>
	class Array3d
	{
	public:
		//		Default constructor: Array3d<double> mat; mat=mat_in;
		Array3d(): col_size(0), row_size(0), dep_size(0), tot_size(0), data(nullptr){;}
		
		//		Constructor for direct initialization: Array3d<double> mat(2,2,3);
		Array3d( size_t col_in, size_t row_in, size_t dep_in ):
		col_size(col_in), row_size(row_in), dep_size(dep_in), tot_size(col_in*row_in*dep_in)
		{
			assert( col_in>=0 && row_in>=0 && dep_in>=0 ); // make sure the lengths >=0
			try
			{
				data = new T[tot_size] ;
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			
		}
		
		//		Copy constructor: Array3d<double> mat(mat_in) or Array3d<double> mat=mat_in;
		Array3d( const Array3d<T> &a):
		col_size(a.col_size), row_size(a.row_size), dep_size(a.dep_size), tot_size(a.tot_size)
		{
			try
			{
				data = new T[tot_size];
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = a.data[i];
			}
		}
		
		~Array3d()
		{
			delete [] data;
			data = nullptr;
		}
		
		//		Overload operator () to access matrix elements. Note no range-check is performed.
		T& operator()( size_t col, size_t row, size_t dep )
		{
			return data[ col + row * col_size + dep * col_size * row_size ];
		}
		
		//		Overload operator [] to access matrix as a 1d array.
		T& operator[](size_t i)
		{
			return data[i];
		}
		
		//		Overload copy assignment operator.
		Array3d& operator=( const Array3d<T>& a )
		{
			if (this == &a) return *this;
			if (data)
			{
				delete [] data;
			}
			col_size = a.col_size;
			row_size = a.row_size;
			dep_size = a.dep_size;
			tot_size = a.tot_size;
			try
			{
				data = new T[tot_size];
			}
			catch (std::bad_alloc&)
			{
				cout << "Error allocating memory !" << endl;
				cout << "At class " << classname() << endl;
				std::exit(EXIT_FAILURE);
			}
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = a.data[i];
			}
			return *this;
		}
		
		//		Overload operator = to assign value to the whole matrix.
		Array3d& operator=( const T& val )
		{
			for (size_t i=0; i<tot_size; i++)
			{
				data[i] = val;
			}
			return *this;
		}
		
		//		Methods to access the information about the matrix.
		size_t colsize() const { return col_size; }
		size_t rowsize() const { return row_size; }
		size_t depsize() const { return dep_size; }
		size_t totsize() const { return tot_size; }
		
		//	private:
		size_t col_size;
		size_t row_size;
		size_t dep_size;
		size_t tot_size;
		T* data;
		const char* classname() { return "Array3d"; }
		
	};

	template <class T>
	void print_test_Array3d( Array3d<T>& mat)
	{
		cout << "***********************************************" << endl;
		cout << " Test of the Array3d class" << endl;
		cout << "mat.colsize()=" << mat.colsize() << endl;
		cout << "mat.rowsize()=" << mat.rowsize() << endl;
		cout << "mat.depsize()=" << mat.depsize() << endl;
		cout << "mat.totsize()=" << mat.totsize() << endl;
		
		for ( int m=0; m<mat.totsize(); m++ )
		{
			printf( "initialized mat(%i)=%.15f \n", m, mat[m] );
			mat[m] = double(m);
		}
		for ( int k=0; k<mat.depsize(); k++)
			for ( int j=0; j<mat.rowsize(); j++)
				for ( int i=0; i<mat.colsize(); i++)
				{
					printf( "assigned mat(%i, %i, %i)=%.15f \n", i, j, k, mat(i,j,k) );
				}
		cout << "Success run !\n";
		cout << "***********************************************" << endl;
	}
}

#endif
