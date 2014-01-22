#include <stdio.h>

int main() {
	int i, j, size, q_max ;
	FILE *fp ;
	
	printf("Specify detector size as an odd number: ") ;
	scanf("%d", &size) ;
	q_max = (size - 1) / 2 ;
		
	fp = fopen("det.dat", "w") ;
	fprintf(fp, "%d\n", q_max) ;
	for (i = -q_max ; i <= q_max ; ++i) 
	for (j = -q_max ; j <= q_max ; ++j)
		if (i*i + j*j <= q_max*q_max)
			fprintf(fp, "%d %d\n", i, j) ;
	fclose(fp) ;
	
	return 0 ;
}
