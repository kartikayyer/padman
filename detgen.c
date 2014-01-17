#include <stdio.h>

int main() {
	int i, j, size, q_max ;
	FILE *fp ;
	
	printf("Specify detector size as an odd number: ") ;
	scanf("%d", &size) ;
	q_max = (size - 1) / 2 ;
		
	fp = fopen("det.dat", "w") ;
	fprintf(fp, "%d\n", size) ;
	for (i = 0 ; i < size ; ++i) 
	for (j = 0 ; j < size ; ++j)
		if ( (i - q_max)*(i - q_max) + (j - q_max)*(j - q_max) <= q_max*q_max )
			fprintf(fp, "%d %d \n", i - q_max, j - q_max) ;
	fclose(fp) ;
	
	return 0 ;
}
