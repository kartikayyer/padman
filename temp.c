/* Reconstruction of mask pattern from randomly oriented patterns with low photon count

Compile:
gcc recon.c -o recon -lm -O3

Needs:
photons.dat, det.dat

Generates:
finish.dat

Usage:
./recon <iter>
where <iter> is the number of iteration you want to perform

Log file:
EMC.log

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define EULER 0.57721566490153286
#define P_MIN .001
#define PI 3.14159265

double rot[2][2] ;
double (*pix)[2] ;
int *ones, *multi, **place_ones, **place_multi, **count ;
double **intens1, **intens2, **inter_weight, *s, *u, *p ;
float **w_in, **w_out ;
int q_max, size, num_pix, num_rot, num_data ;
double mean_count, rescale, intens_norm, frac, w_b ;

int setup() ;
double emc() ;
void expand() ;   // E
void maximize() ; // M
void compress() ; // C
void print_intens() ;
void free_mem() ;
void make_rot( int ) ;


int main(int argc, char* argv[]) {
	int iter, i ;
	double rms_change ;
	FILE *fp ;
	
	if ( argc == 2 ) iter = atoi(argv[1]) ;
	else {
		printf("expected one argument: iter\n") ;
		return 0 ;
	}
		
	if (!setup()) return 0 ;
	fprintf(stderr, "Completed setup\n") ;
	
	fp = fopen("EMC.log", "w") ;
	fprintf(fp, "num_rot = %d    num_data = %d    num_pix = %d\n", num_rot, num_data, num_pix) ;
	fprintf(fp, "total_mean_count = %f  signal_count = %f\n\n", mean_count * (1 + frac), mean_count) ; 
	fclose(fp) ;
	
	for (i = 1 ; i <= iter ; ++i) {
		rms_change = emc() ;
		
		fp = fopen("EMC.log", "a") ;
		fprintf(fp, "iter = %d    rms_change = %f\n\n", i, rms_change) ;
		fclose(fp) ;
		
		print_intens(i) ;
	}
	
	free_mem() ;
	
	return 0 ;
}
	
	
double emc() {
	double info, change = 0., d ;
	int x, y ;
	time_t t1, t2 ;
	FILE *fp ;
	
	time(&t1) ;
	
	expand() ;
	maximize() ;
	compress() ;
	
	time(&t2) ;
		
	fp = fopen("EMC.log", "a") ;
	fprintf(fp, "	iter_time = %.0f sec\n", difftime(t2, t1)) ;
	fclose(fp) ;
	
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
	double max_exp, p_sum, weight, w, h ;
	double rot_intens, total_intens = 0. ;
	int *place_ones_d, *place_multi_d, *count_d ;
	float *w_in_r, *w_out_r ;
	FILE *fp ;
	
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
			w_in[r][t] = log(w + w_b) ;
			
			w_out[r][t] = 0. ;
		}
	}
	
    	for (d = 0 ; d < num_data ; ++d) {
		place_ones_d = place_ones[d] ;
		place_multi_d = place_multi[d] ;
		count_d = count[d] ;
		
		max_exp = -100. * num_pix ;
		for (r = 0 ; r < num_rot ; ++r) {
			p[r] = u[r] ;
			
			w_in_r = w_in[r] ;
			
			for (t = 0 ; t < ones[d] ; ++t)
				p[r] += w_in_r[place_ones_d[t]] ;
				
			for (t = 0 ; t < multi[d] ; ++t)
				p[r] += count_d[t] * w_in_r[place_multi_d[t]] ;
				
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
			if (p[r] < P_MIN) continue ;
				
			s[r] += p[r] ;
				
			w_out_r = w_out[r] ;
			w_in_r = w_in[r] ;
			
			for (t = 0 ; t < ones[d] ; ++t) {
				h = exp(w_in_r[place_ones_d[t]]) ;
				w_out_r[place_ones_d[t]] += p[r] * h / (h + w_b) ;
			}
			
			for (t = 0 ; t < multi[d] ; ++t) {
				h = exp(w_in_r[place_multi_d[t]]) ;
				w_out_r[place_multi_d[t]] += p[r] * count_d[t] * h / (h + w_b) ;
			}
		}
	}
    	
	for (r = 0 ; r < num_rot ; ++r) {
		if (s[r] == 0.) continue ;
			
		for (t = 0 ; t < num_pix ; ++t)
			w_out[r][t] /= s[r] ;
	}
}


void expand() {
	int r, t, i, j ;
	double rot_pix[3] ;
	int x, y ;
	double tx, ty, fx, fy, cx, cy ;
	
	FILE *fp ;
	char fname[100] ;
	double **array ;
	array = malloc(size * sizeof(double*)) ;
	for (i = 0 ; i < size ; ++i) 
		array[i] = calloc(size, sizeof(double)) ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		make_rot(r) ;
		
		for (t = 0 ; t < num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*pix[t][j] ;
			}
				
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			w_in[r][t] = cx*cy*intens1[x][y] + cx*fy*intens1[x][y+1] + fx*cy*intens1[x+1][y] + fx*fy*intens1[x+1][y+1] ;
		}
		
		if (r % 5 == 0) {
			for (t = 0 ; t < num_pix ; ++t)
				array[(int)pix[t][0] + q_max][(int)pix[t][1] + q_max] = w_in[r][t] ;
			sprintf(fname, "expand%.3d.dat", r) ;
			fp = fopen(fname, "w") ;
			for (i = 0 ; i < size ; ++i) {
				for (j = 0 ; j < size ; ++j)
					fprintf(fp, "%.10e ", array[i][j]) ;
				fprintf(fp, "\n") ;
			}
		}
	}
	
	for (i = 0 ; i < size ; ++i)
		free(array[i]) ;
	free(array) ;
}


void compress() {
	int r, t, i, j ;
	double rot_pix[2] ;
	int x, y ;
	double tx, ty, fx, fy, cx, cy ;
	double w, f ;
	double norm ;
	
	FILE *fp ;
	char fname[100] ;
	double **array ;
	array = malloc(size * sizeof(double*)) ;
	for (i = 0 ; i < size ; ++i) 
		array[i] = calloc(size, sizeof(double)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y) {
		inter_weight[x][y] = 0. ;
		intens2[x][y] = 0. ;
	}
	
	for (r = 0 ; r < num_rot ; ++r) {
		if (s[r] == 0.) continue ;
			
		make_rot(r) ;
		
		for (t = 0 ; t < num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*pix[t][j] ;
			}
			
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			w = w_out[r][t] ;
			
			f = cx*cy ;
			inter_weight[x][y] += f ;
			intens2[x][y] += f * w ;
			
			f = cx*fy ;
			inter_weight[x][y+1] += f ;
			intens2[x][y+1] += f * w ;
			
			f = fx*cy ;
			inter_weight[x+1][y] += f ;
			intens2[x+1][y] += f * w ;
			
			f = fx*fy ;
			inter_weight[x+1][y+1] += f ;
			intens2[x+1][y+1] += f * w ;
			
		}
		
		if (r % 5 == 0) {
			for (t = 0 ; t < num_pix ; ++t)
				array[(int)pix[t][0] + q_max][(int)pix[t][1] + q_max] = w_out[r][t] ;
			sprintf(fname, "compress%.3d.dat", r) ;
			fp = fopen(fname, "w") ;
			for (i = 0 ; i < size ; ++i) {
				for (j = 0 ; j < size ; ++j)
					fprintf(fp, "%.10e ", array[i][j]) ;
				fprintf(fp, "\n") ;
			}
		}
	}
		
	norm = intens_norm / mean_count ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
		if (inter_weight[x][y] != 0.) intens2[x][y] *= norm / inter_weight[x][y] ;
	
	for (i = 0 ; i < size ; ++i)
		free(array[i]) ;
	free(array) ;
}
	
	
void print_intens(int iter) {
	FILE *fp ;
	int x, y ;
	char fname[100] ;
	
	sprintf(fname, "finish%.3d.dat", iter) ;
	fp = fopen(fname, "w") ;
	
	for (x = 0 ; x < size ; ++x) {
		for (y = 0 ; y < size ; ++y)
			fprintf(fp, "%21.15e ", intens2[x][y]) ;
			
		fprintf(fp, "\n") ;
	}
	
	fclose(fp) ;
}
	
	
int setup() {
	FILE *fp ;
	int i, j, k, d, r, t ;
	int rand_start = 0 ;
	char buffer[50] ;
	double photon_count, total_count = 0. ;
	srand(time(NULL)) ;
		
	fp = fopen("det.dat", "r") ;
	if (!fp) {
		printf("cannot open det.dat\n") ;
		return 0 ;
	}
	
	fscanf(fp, "%d %d", &num_pix, &q_max) ;
	
	size = 2*q_max + 1 ;
	
	pix = malloc(num_pix * sizeof(*pix)) ;
	
	for (t = 0 ; t < num_pix ; ++t)
	for (i = 0 ; i < 2 ; ++i)
		fscanf(fp, "%lf", &pix[t][i]) ;
		
	fclose(fp) ;
	
	fp = fopen("start_model.dat", "r") ;
	size++ ;	
	intens1 = malloc(size * sizeof(double*)) ;
	intens2 = malloc(size * sizeof(double*)) ;
	inter_weight = malloc(size * sizeof(double*)) ;
	for (i = 0 ; i < size ; ++i) {
		intens1[i] = malloc(size * sizeof(double)) ;
		intens2[i] = malloc(size * sizeof(double)) ;
		inter_weight[i] = malloc(size * sizeof(double)) ;
		for (j = 0 ; j < size ; ++j)
			fscanf(fp, "%lf", &intens1[i][j]) ;
	}
	fclose(fp) ;
	
	num_rot = 45 ;
	
	fp = fopen("photons.dat", "r") ;
	if (!fp) {
		printf("cannot open photons.dat\n") ;
		return 0 ;
	}
	
	fscanf(fp, "%d  %lf", &num_data, &frac) ;
	intens_norm = 2000. ; 		// Chosen to make the image easily visible on a scale of 0-1
	
	ones = malloc(num_data * sizeof(int)) ;
	place_ones = malloc(num_data * sizeof(int*)) ;
	multi = malloc(num_data * sizeof(int)) ;
	place_multi = malloc(num_data * sizeof(int*)) ;
	count = malloc(num_data * sizeof(int*)) ;
	
	for (d = 0 ; d < num_data ; ++d) {
		fscanf(fp, "%d", &ones[d]) ;
		
		place_ones[d] = malloc(ones[d] * sizeof(int)) ;
		
		for (t = 0 ; t < ones[d] ; ++t)
			fscanf(fp, "%d", &place_ones[d][t]) ;
		
		
		fscanf(fp, "%d", &multi[d]) ;
		
		place_multi[d] = malloc(multi[d] * sizeof(int)) ;
		count[d] = malloc(multi[d] * sizeof(int)) ;
		
		for (t = 0 ; t < multi[d] ; ++t)
			fscanf(fp, "%d %d", &place_multi[d][t], &count[d][t]) ;
	}
	
	fclose(fp) ;
	
	w_in = malloc(num_rot * sizeof(double*)) ;
	w_out = malloc(num_rot * sizeof(double*)) ;
	s = malloc(num_rot * sizeof(double)) ;
	u = malloc(num_rot * sizeof(double)) ;
	p = malloc(num_rot * sizeof(double)) ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		w_in[r] = malloc(num_pix * sizeof(double)) ;
		w_out[r] = malloc(num_pix * sizeof(double)) ;
	}
	
	for (d = 0 ; d < num_data ; ++d) {
		photon_count = ones[d] ;
		for (t = 0 ; t < multi[d] ; ++t)
			photon_count += count[d][t] ;
		
		total_count += photon_count ;
	}
		
	mean_count = total_count / num_data / (1 + frac) ;
	w_b = frac * mean_count / num_pix ;
	
	return 1 ;
}


void free_mem() {
	int i, j, d, r ;
	
	free(pix) ;
	
	for (i = 0 ; i < size ; ++i) {
		free(intens1[i]) ;
		free(intens2[i]) ;
		free(inter_weight[i]) ;
	}
	free(intens1) ;
	free(intens2) ;
	free(inter_weight) ;
	
	for (d = 0 ; d < num_data ; ++d) {
		free(place_ones[d]) ;
		free(place_multi[d]) ;
		free(count[d]) ;
	}
	free(ones) ;
	free(multi) ;
	free(place_ones) ;
	free(place_multi) ;
	free(count) ;
	
	for (r = 0 ; r < num_rot ; ++r) {
		free(w_in[r]) ;
		free(w_out[r]) ;
	}
	free(w_in) ;
	free(w_out) ;
	
	free(s) ;
	free(u) ;
	free(p) ;
}


void make_rot( int r ) {
	double a ;
	
	a = r * 2 * PI / num_rot ;

	rot[0][0] = cos(a) ;
	rot[0][1] = sin(a) ;
	rot[1][0] = -sin(a);
	rot[1][1] = cos(a) ;
}
