#ifndef STRUCTS_HH
#define STRUCTS_HH

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

#define NUM_ROT 100
#define INTENS_NORM 2000.

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
} ;

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
} ;

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
} ;

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
} ;

void expand(det_t*, model_t*, view_t*) ;
void maximize(data_t*, view_t*) ;
void compress(det_t*, model_t*, view_t*) ;

#endif // STRUCTS_HH
