/* Data generation for simulation of 2-d reconstruction of mask pattern

This program generates diffraction patterns by randomly rotating and sampling
the data according to the Poisson distribution.

Compile:
gcc -o datagen.c -o datagen -lm -O3

Needs:
mask.dat, det.dat

Generates:
photons.dat

Usage:
./datagen <num> <mean_count> <frac>
where
	<num> represents number of patterns to be generated
	<mean_count> represents mean number of signal photons per hit
	<frac> represents number of background photons as a fraction of signal photons
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.14159265

double rot[2][2] ;
double (*pix)[2] ;
int ones, multi, *place_ones, *place_multi, *count ;
double **intens ;
int q_max, size, m_pix ;
double mean_count, frac, intens_water ;

int setup() ;
void free_mem() ;
void make_rot() ;
int poisson( double ) ;


int main(int argc, char* argv[]) {
	int num, d, t ;
	int i, j, k ;
	double rot_pix[3] ;
	int x, y, z ;
	double tx, ty, tz, fx, fy, fz, cx, cy, cz ;
	int photons ;
	double intens_val, intens_ave = 0. ;
	int m_ave = 1000 ;
	FILE *fp ;
	
	if ( argc == 4 ) {
		num = atoi(argv[1]) ;
		mean_count = atof(argv[2]) ;
		frac = atof(argv[3]) ;
	}
	else {
		printf("expected three arguments: num, mean_count, frac\n") ;
		return 0 ;
	}
	
	if (!setup())
		return 0 ;
	
	intens_water = mean_count * frac / m_pix ;
	
	srand( time(0) ) ;
	
	fp = fopen("photons.dat", "w") ;
	
	for (d = 0 ; d < num + m_ave ; ++d) {
		make_rot() ;
		
		ones = 0 ;
		multi = 0 ;
		
		if (d == m_ave) {
			intens_ave /= m_ave ;
			
			fprintf(fp, "%d  %f\n\n", num, frac) ;
			
			for (i = 0 ; i < size ; ++i)
			for (j = 0 ; j < size ; ++j)
				intens[i][j] *= mean_count / intens_ave ;
		}
		
		for (t = 0 ; t < m_pix ; ++t) {
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
			
			intens_val = fx*fy*intens[x][y] + fx*cy*intens[x][y+1] + cx*fy*intens[x+1][y] + cx*cy*intens[x+1][y+1] ;
			
			if (d < m_ave) {
				intens_ave += intens_val ;
				continue ;
			}
			
			photons = poisson( intens_val + intens_water) ;
			
			if (!photons)
				continue ;
			
			if (photons == 1) {
				place_ones[ones] = t ;
				++ones ;
			}
			else {
				place_multi[multi] = t ;
				count[multi] = photons ;
				++multi ;
			}
		}
		
		if (d < m_ave)
			continue ;
		
		fprintf(fp, "%d\n", ones) ;
		for (t = 0 ; t < ones ; ++t)
			fprintf(fp, "%d ", place_ones[t]) ;
		fprintf(fp, "\n") ;
		
		fprintf(fp, "%d\n", multi) ;
		for (t = 0 ; t < multi ; ++t)
			fprintf(fp, "%d %d  ", place_multi[t], count[t]) ;
		fprintf(fp, "\n\n") ;
	}
	
	fclose(fp) ;
	
	free_mem() ;
	
	return 0 ;
}


int setup() {
	FILE *fp ;
	int i, j, k, t ;
	
	fp = fopen("det.dat", "r") ;
	if (!fp) {
		printf("cannot open det.dat\n") ;
		return 0 ;
	}
	
	fscanf(fp, "%d ", &size) ;
	
	m_pix = size * size ;
	q_max = (size - 1) / 2 ;
		
	place_ones = malloc(m_pix * sizeof(int*)) ;
	place_multi = malloc(m_pix * sizeof(int*)) ;
	count = malloc(m_pix * sizeof(int*)) ;
	
	pix = malloc(m_pix * sizeof(*pix)) ;
	
	for (t = 0 ; t < m_pix ; ++t)
	for (i = 0 ; i < 2 ; ++i)
		fscanf(fp, "%lf", &pix[t][i]) ;
	
	fclose(fp) ;
	
	fp = fopen("mask.dat", "r") ;
	if (!fp) {
		printf("cannot open mask.dat\n") ;
		return 0 ;
	}
	
	intens = malloc(size * sizeof(double*)) ;
	for (i = 0 ; i < size ; ++i) {
		intens[i] = malloc(size * sizeof(double)) ;
		for (j = 0 ; j < size ; ++j)
			fscanf(fp, "%lf", &intens[i][j]) ;
	}
	
	fclose(fp) ;
	
	return 1 ;
}
	
	
void free_mem() {
	int i ;
	
	free(place_ones) ;
	free(place_multi) ;
	free(count) ;
	
	free (pix) ;
	
	for (i = 0 ; i < size ; ++i)
		free(intens[i]) ;
	free(intens) ;
}
	
	
void make_rot() {
	double a ;
	
	a = ((double) rand()) / RAND_MAX * 2 * PI ;
	
	rot[0][0] = cos(a) ;
	rot[0][1] = sin(a) ;
	rot[1][0] = -sin(a) ;
	rot[1][1] = cos(a) ;
}


int poisson( double m ) {
	int i = 0 ;
	double p, q, r ;
	
	r = exp(-m) ;
	p = r ;
	q = ((double) rand()) / RAND_MAX ;
	
	while (p < q) {
		++i ;
		r *= m / i ;
		p += r ;
	}
		
	return i ;
}
