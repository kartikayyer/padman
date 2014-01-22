#ifndef VIEW_HH
#define VIEW_HH

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

#define NUM_ROT 100

struct view_t {
	double **in, **out ;
	int num_pix, num_rot ;
	
	int init() {
		fstream fp ;
		int r ;
		
		fp.open("det.dat", ios::in) ;
		if (!fp.is_open()) {
			cerr << "cannot open det.dat\n" ;
			return 0 ;
		}
		fp >> num_pix ;
		fp.close() ;
		
		num_rot = NUM_ROT ;
		in = new double*[num_rot] ;
		out = new double*[num_rot] ;
		
		for (r = 0 ; r < num_rot ; ++r) {
			in[r] = new double[num_pix] ;
			out[r] = new double[num_pix] ;
		}
		
		return 1 ;
	}
	
	~view_t() {
        int r ;
        
        for (r = 0 ; r < num_rot ; ++r) {
            delete[] in[r] ;
            delete[] out[r] ;
        }
        delete[] in ;
        delete[] out ;
	}
} ;

#endif // VIEW_HH
