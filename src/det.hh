#ifndef DET_HH
#define DET_HH

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

struct det_t {
	double **pix ;
	int num_pix ;
	
	int init() {
		fstream fp ;
		int t, i, size ;
		
		fp.open("det.dat", ios::in) ;
		if (!fp.is_open()) {
			cerr << "cannot open det.dat\n" ;
			return 0 ;
		}
		fp >> num_pix >> size ;
		pix = new double*[num_pix] ;
		for (t = 0 ; t < num_pix ; ++t) {
			pix[t] = new double[2] ;
			for (i = 0 ; i < 2 ; ++i)
				fp >> pix[t][i] ;
		}
		fp.close() ;
		
		return 1 ;
	}
	
	~det_t() {
	    int i ;
	    
	    for (i = 0 ; i < num_pix ; ++i)
	        delete[] pix[i] ;
	    delete[] pix ;
	}
} ;

#endif // DET_HH
