#ifndef MODEL_HH
#define MODEL_HH

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

#define INTENS_NORM 2000.

struct model_t {
	double **in, **out, **weight ;
	double intens_norm ;
	int size ;
	
	void print_intens() {
		fstream fp ;
		int x, y ;
		
		fp.open("finish.dat", ios::out) ;
		
		for (x = 0 ; x < size ; ++x) {
			for (y = 0 ; y < size ; ++y)
				fp << out[x][y] << " ";
				
			fp << endl ;
		}
		
		fp.close() ;
	}
	
	double diff_intens() {
		int x, y ;
		double d, change = 0. ;
		
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y) {
			d = in[x][y] - out[x][y] ;
			change += d*d ;
			in[x][y] = out[x][y] ;
		}
		
		return sqrt(change / (size * size)) ;
	}
	
	int init() {
		fstream fp ;
		int num, i, j, qmax ;
		
		fp.open("det.dat", ios::in) ;
		if (!fp.is_open()) {
			cerr << "cannot open det.dat\n" ;
			return 0 ;
		}
		fp >> num >> qmax ;
		size = 2 * qmax + 2 ;
		fp.close() ;
		
		intens_norm = INTENS_NORM ; 	// Chosen to make the image easily visible on a scale of 0-1
		
		in = new double*[size] ;
		out = new double*[size] ;
		weight = new double*[size] ;
		for (i = 0 ; i < size ; ++i) {
			in[i] = new double[size] ;
			out[i] = new double[size] ;
			weight[i] = new double[size] ;
			for (j = 0 ; j < size ; ++j)
				in[i][j] = ((double) rand()) / RAND_MAX ;
		}
		
		return 1 ;
	}
	
	~model_t() {
		int i ;
		
		for (i = 0 ; i < size ; ++i) {
			delete[] in[i] ;
			delete[] out[i] ;
			delete[] weight[i] ;
		}
		
		delete[] in ;
		delete[] out ;
		delete[] weight ;
	}
} ;

#endif // MODEL_HH
