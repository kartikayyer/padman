/* Reconstruction of mask pattern from randomly oriented patterns with low photon count

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace std ;

#define EULER 0.57721566490153286
#define P_MIN .001
#define PI 3.14159265

double **pix ;
struct frame_t {
	int ones, multi ;
	int *place_ones, *place_multi, *count ;
} *frame ;
double **intens1, **intens2, **inter_weight, *s, *u, *p ;
double **w_in, **w_out ;
int q_max, size, num_pix, num_rot, num_data ;
double mean_count, rescale, intens_norm ;

int setup() ;
double emc() ;
void expand(double**, double**, double**) ;
void maximize() ;
void compress(double**, double**, double**, double**) ;
void print_intens() ;
void free_mem() ;
void make_rot(int, double[][2]) ;


int main(int argc, char* argv[]) {
	int iter, i ;
	double rms_change ;
	fstream fp ;
	
	if ( argc == 2 )
		iter = atoi(argv[1]) ;
	else {
		printf("expected one argument: iter\n") ;
		return 0 ;
	}
		
	if (!setup())
		return 0 ;
	
	fp.open("EMC.log", ios::out) ;
	fp << "num_rot = " << num_rot << "\tnum_data = " << num_data << "\tnum_pix = " << num_pix << endl ;
	fp << "mean_count = " << mean_count << endl << endl ; 
	fp << "iter\trms_change\titer_time\n" ;
	fp.close() ;
	
	for (i = 1 ; i <= iter ; ++i) {
		rms_change = emc() ;
		
		fp.open("EMC.log", ios::out|ios::app) ;
		fp << i << "\t" << rms_change << "\t" ;
		fp.close() ;
	}
	
	print_intens() ;
	
	free_mem() ;
	
	return 0 ;
}
	
	
double emc() {
	double info, change = 0., d ;
	int x, y ;
	time_t t1, t2 ;
	fstream fp ;
	
	time(&t1) ;
	
	expand(pix, intens1, w_in) ;
	maximize() ;
	compress(pix, w_out, intens2, inter_weight) ;
	
	time(&t2) ;
		
	fp.open("EMC.log", ios::out|ios::app) ;
	fp << difftime(t2, t1) << endl ;
	fp.close() ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y) {
		d = intens1[x][y] - intens2[x][y] ;
		change += d*d ;
		intens1[x][y] = intens2[x][y] ;
	}
	
	return (num_pix / intens_norm) * sqrt( change / (size*size) ) ;
}


void maximize() {
	int d, r, t ;
	double max_exp, p_sum ;
	double rot_intens, total_intens = 0. ;
	double weight ;
	double w ;
	frame_t *frame_d ;
	double *w_in_r, *w_out_r ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		rot_intens = 0. ;
		for (t = 0 ; t < num_pix ; ++t)
			rot_intens += w_in[r][t] ;
			
		total_intens += rot_intens / num_rot ;
	}
		
	rescale = mean_count / total_intens ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		s[r] = 0. ;
		u[r] = 0. ;
		
		for (t = 0 ; t < num_pix ; ++t) {
			w = w_in[r][t] * rescale ;
			u[r] -= w ;
			w_in[r][t] = log(w) ;
			
			w_out[r][t] = 0. ;
		}
	}
	
    	for (d = 0 ; d < num_data ; ++d) {
		frame_d = &frame[d] ;
		
		max_exp = -100. * num_pix ;
		for (r = 0 ; r < num_rot ; ++r) {
			p[r] = u[r] ;
			
			w_in_r = w_in[r] ;
			
			for (t = 0 ; t < frame[d].ones ; ++t)
				p[r] += w_in_r[frame_d->place_ones[t]] ;
				
			for (t = 0 ; t < frame[d].multi ; ++t)
				p[r] += frame_d->count[t] * w_in_r[frame_d->place_multi[t]] ;
				
			if (p[r] > max_exp)
				max_exp = p[r] ;
		}
		
		p_sum = 0. ;
		for (r = 0 ; r < num_rot ; ++r) {
			p[r] = exp(p[r] - max_exp) ;
			p_sum += p[r] ;
		}
		
        	for (r = 0 ; r < num_rot ; ++r) {
			p[r] /= p_sum ;
			
			if (p[r] < P_MIN)
				continue ;
				
			s[r] += p[r] ;
				
			w_out_r = w_out[r] ;
			
			for (t = 0 ; t < frame[d].ones ; ++t)
				w_out_r[frame_d->place_ones[t]] += p[r] ;
				
			for (t = 0 ; t < frame[d].multi ; ++t)
				w_out_r[frame_d->place_multi[t]] += p[r] * frame_d->count[t] ;
		}
	}
	
    	for (r = 0 ; r < num_rot ; ++r) {
		if (s[r] == 0.)
			continue ;
			
		for (t = 0 ; t < num_pix ; ++t)
			w_out[r][t] /= s[r] ;
	}
}


void expand(double **pixel, double **intens, double **view) {
	int r, t, i, j ;
	double rot_pix[2] ;
	double rot[2][2] ;
	int x, y ;
	double tx, ty, fx, fy, cx, cy ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		make_rot(r, rot) ;
		
		for (t = 0 ; t < num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*pixel[t][j] ;
			}
				
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			view[r][t] = cx*cy*intens[x][y] + cx*fy*intens[x][y+1] + fx*cy*intens[x+1][y] + fx*fy*intens[x+1][y+1] ;
		}
	}
}


void compress(double **pixel, double **view, double **intens, double **weight) {
	int r, t, i, j ;
	double rot_pix[2] ;
	double rot[2][2] ;
	int x, y ;
	double tx, ty, fx, fy, cx, cy ;
	double w, f ;
	double norm ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y) {
		weight[x][y] = 0. ;
		intens[x][y] = 0. ;
	}
	
	for (r = 0 ; r < num_rot ; ++r) {
//		if (s[r] == 0.)
//			continue ;
			
		make_rot(r, rot) ;
		
		for (t = 0 ; t < num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*pixel[t][j] ;
			}
			
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			w = view[r][t] ;
			
			f = cx*cy ;
			weight[x][y] += f ;
			intens[x][y] += f * w ;
			
			f = cx*fy ;
			weight[x][y+1] += f ;
			intens[x][y+1] += f * w ;
			
			f = fx*cy ;
			weight[x+1][y] += f ;
			intens[x+1][y] += f * w ;
			
			f = fx*fy ;
			weight[x+1][y+1] += f ;
			intens[x+1][y+1] += f * w ;
		}
	}
		
	norm = intens_norm / mean_count ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
		if (weight[x][y] != 0.)
			intens[x][y] *= norm / weight[x][y] ;
}
	
	
void print_intens() {
	fstream fp ;
	int x, y ;
	
	fp.open("finish.dat", ios::out) ;
	
	for (x = 0 ; x < size ; ++x) {
		for (y = 0 ; y < size ; ++y)
			fp << intens2[x][y] ;
			
		fp << endl ;
	}
		
	fp.close() ;
}
	
	
int setup() {
	fstream fp ;
	int i, j, k, d, r, t ;
	int rand_start = 0 ;
	char buffer[50] ;
	
	num_rot = 1000 ;
	intens_norm = 2000. ; 		// Chosen to make the image easily visible on a scale of 0-1
	
	fp.open("det.dat", ios::in) ;
	if (!fp) {
		printf("cannot open det.dat\n") ;
		return 0 ;
	}
		
	fp >> num_pix >> size ;
	
	q_max = (size - 1) / 2 ;
	
	pix = new double*[num_pix] ;
	
	for (t = 0 ; t < num_pix ; ++t) {
		pix[t] = new double[2] ;
		for (i = 0 ; i < 2 ; ++i) 
			fp >> pix[t][i] ;
	}
	
	fp.close() ;
	
	size += 1 ;	
	intens1 = new double*[size] ;
	intens2 = new double*[size] ;
	inter_weight = new double*[size] ;
	for (i = 0 ; i < size ; ++i) {
		intens1[i] = new double[size] ;
		intens2[i] = new double[size] ;
		inter_weight[i] = new double[size] ;
		for (j = 0 ; j < size ; ++j)
			intens1[i][j] = ((double) rand()) / RAND_MAX ;
	}
	
	fp.open("photons.dat", ios::in) ;
	if (!fp) {
		printf("cannot open photons.dat\n") ;
		return 0 ;
	}
	
	fp >> num_data ;
	
	frame = new frame_t[num_data] ;
	
	for (d = 0 ; d < num_data ; ++d) {
		fp >> frame[d].ones ;
		mean_count = frame[d].ones ;
		
		frame[d].place_ones = new int[frame[d].ones] ;
		
		for (t = 0 ; t < frame[d].ones ; ++t)
			fp >> frame[d].place_ones[t] ;
		
		fp >> frame[d].multi ;
		
		frame[d].place_multi = new int[frame[d].multi] ;
		frame[d].count = new int[frame[d].multi] ;
	
		for (t = 0 ; t < frame[d].multi ; ++t) {
			fp >> frame[d].place_multi[t] >> frame[d].count[t] ;
			mean_count += frame[d].count[t] ;
		}
	}
	
	mean_count /= num_data ;

	fp.close() ;
	
	w_in = new double*[num_rot] ;
	w_out = new double*[num_rot] ;
	s = new double[num_rot] ;
	u = new double[num_rot] ;
	p = new double[num_rot] ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		w_in[r] = new double[num_pix] ;
		w_out[r] = new double[num_pix] ;
	}
	
	return 1 ;
}


void free_mem() {
	int i, j, d, r ;
	
	delete[] pix ;
	
	for (i = 0 ; i < size ; ++i) {
		delete[] intens1[i] ;
		delete[] intens2[i] ;
		delete[] inter_weight[i] ;
	}
	delete[] intens1 ;
	delete[] intens2 ;
	delete[] inter_weight ;
	
	for (d = 0 ; d < num_data ; ++d) {
		delete[] frame[d].place_ones ;
		delete[] frame[d].place_multi ;
		delete[] frame[d].count ;
	}
	delete[] frame ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		delete[] w_in[r] ;
		delete[] w_out[r] ;
	}
	delete[] w_in ;
	delete[] w_out ;
	
	delete[] s ;
	delete[] u ;
	delete[] p ;
}


void make_rot(int r, double r_matrix[][2]) {
	double a ;
	
	a = r * 2 * PI / num_rot ;

	r_matrix[0][0] = cos(a) ;
	r_matrix[0][1] = sin(a) ;
	r_matrix[1][0] = -sin(a);
	r_matrix[1][1] = cos(a) ;
}
