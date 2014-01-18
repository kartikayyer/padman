/* Reconstruction of mask pattern from randomly oriented patterns with low photon count

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "structs.hpp"
using namespace std ;

#define EULER 0.57721566490153286

double diff_intens(model_t*) ;
void print_intens(model_t*) ;

int main(int argc, char* argv[]) {
	int iter, i ;
	double rms_change ;
	time_t t1, t2 ;
	fstream fp ;
	view_t View ;
	data_t Data ;
	model_t Model ;
	det_t Det ;
	
	if ( argc == 2 )
		iter = atoi(argv[1]) ;
	else {
		printf("expected one argument: iter\n") ;
		return 0 ;
	}
	
	if (!setup(&Data, &Det, &Model, &View))
		return 0 ;
	
	fp.open("EMC.log", ios::out) ;
	fp << "num_rot = " << View.num_rot << "\tnum_data = " << Data.num_data << "\tnum_pix = " << View.num_pix << endl ;
	fp << "mean_count = " << Data.mean_count << endl << endl ; 
	fp << "iter\trms_change\titer_time\n" ;
	fp.close() ;
	
	for (i = 1 ; i <= iter ; ++i) {
		time(&t1) ;
		
		expand(&Det, &Model, &View) ;
		maximize(&Data, &View) ;
		compress(&Det, &Model, &View) ;
		
		rms_change = diff_intens(&Model) * View.num_pix / Model.intens_norm * Data.mean_count ;
		print_intens(&Model) ;
		
		time(&t2) ;
		
		fp.open("EMC.log", ios::out|ios::app) ;
		fp << i << "\t" << rms_change << "\t" << difftime(t2, t1) << endl ;
		fp.close() ;
	}
	
	free_mem(&Data, &Det, &Model, &View) ;
	
	return 0 ;
}


void print_intens(model_t *model) {
	fstream fp ;
	int x, y ;
	
	fp.open("finish.dat", ios::out) ;
	
	for (x = 0 ; x < model->size ; ++x) {
		for (y = 0 ; y < model->size ; ++y)
			fp << model->out[x][y] ;
			
		fp << endl ;
	}
		
	fp.close() ;
}


double diff_intens(model_t *model) {
	int x, y ;
	double d, change = 0. ;
	
	for (x = 0 ; x < model->size ; ++x)
	for (y = 0 ; y < model->size ; ++y) {
		d = model->in[x][y] - model->out[x][y] ;
		change += d*d ;
		model->in[x][y] = model->out[x][y] ;
	}
	
	return sqrt(change / (model->size * model->size)) ;
}
