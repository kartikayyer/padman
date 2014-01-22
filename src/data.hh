#ifndef DATA_HH
#define DATA_HH

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

struct frame_t {
	int ones, multi ;
	int *place_ones, *place_multi, *count ;
} ;

struct data_t {
	frame_t *frame ;
	double mean_count ;
	int num_data ;
	
	int init() {
		fstream fp ;
		int d, t ;
		mean_count = 0. ;
		
		fp.open("photons.dat", ios::in) ;
		if (!fp.is_open()) {
			cerr << "cannot open photons.dat\n" ;
			return 0 ;
		}
		
		fp >> num_data ;
		
		frame = new frame_t[num_data] ;
		
		for (d = 0 ; d < num_data ; ++d) {
			fp >> frame[d].ones ;
			mean_count += (double) frame[d].ones ;
			
			frame[d].place_ones = new int[frame[d].ones] ;
			
			for (t = 0 ; t < frame[d].ones ; ++t)
				fp >> frame[d].place_ones[t] ;
			
			fp >> frame[d].multi ;
			
			frame[d].place_multi = new int[frame[d].multi] ;
			frame[d].count = new int[frame[d].multi] ;
			
			for (t = 0 ; t < frame[d].multi ; ++t) {
				fp >> frame[d].place_multi[t] >> frame[d].count[t] ;
				mean_count += (double) frame[d].count[t] ;
			}
		}
		
		mean_count /= (double) num_data ;
		
		return 1 ;
	}
    
    ~data_t() {
        int d ;
        
        for (d = 0 ; d < num_data ; ++d) {
            delete[] frame[d].place_ones ;
            delete[] frame[d].place_multi ;
            delete[] frame[d].count ;
        }
        delete[] frame ;
	}
} ;

#endif // DATA.HH
