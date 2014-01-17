#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

int **mask ; 
double centers[10][2] ;

int main() {
	int i, j, k ;
	int size = 101 ;
	double radius, r ;
	radius = (double) sqrt(3) * size / 20 ;
	double point[2] ;
	FILE *fp ;
	
	srand( time(NULL) ) ;
	
	mask = malloc(size * sizeof(int*)) ;
	for (i = 0 ; i < size ; ++i)
		mask[i] = calloc(size, sizeof(int)) ;
	
	k = 0 ;
	while (1) {
		if (k > 9) break ;
		
		point[0] = (((double) rand() / RAND_MAX) * 2. - 1.) * size / 2;
		point[1] = (((double) rand() / RAND_MAX) * 2. - 1.) * size / 2;
		
		if (point[1]*point[1] + point[0]*point[0] > size * size / 4)
			continue ;
		
		centers[k][0] = point[0] ;
		centers[k][1] = point[1] ;
		
		k++ ;
	}
	
	for (k = 0 ; k < 10 ; ++k) {
		r = sqrt(centers[k][0]*centers[k][0] + centers[k][1]*centers[k][1]) ;
		centers[k][0] *= (1 - 4 * radius / size) + radius / r ;
		centers[k][1] *= (1 - 4 * radius / size) + radius / r ;
	}
	
	for (i = 0 ; i < size ; ++i)
	for (j = 0 ; j < size ; ++j) 
	for (k = 0 ; k < 10 ; ++k) 
	if ( sqrt((i - size/2 - centers[k][0])*(i - size/2 - centers[k][0]) + (j - size/2 - centers[k][1])*(j - size/2 - centers[k][1])) < radius)
		mask[i][j] = 1 ;
		
	fp = fopen("mask.dat", "w") ;
	for (i = 0 ; i < size ; ++i) {
		for (j = 0 ; j < size ; ++j) 
			fprintf(fp, "%d ", mask[i][j]) ;
		fprintf(fp, "\n") ;
	}
	fclose(fp) ;
	
	for (i = 0 ; i < size ; ++i)
		free(mask[i]) ;
	free(mask) ;
	
	return 0 ;
}
