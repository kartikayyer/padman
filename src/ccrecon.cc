/* Reconstruction of mask pattern from randomly oriented patterns with low photon count

*/

#include <iomanip>
#include <time.h>
#include <math.h>
#include "structs.hh"
using namespace std ;

void free_mem(data_t*, det_t*, model_t*, view_t*) ;

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
		cerr << "expected one argument: iter\n" ;
		return 0 ;
	}
	
	if (!View.init())
		return 1 ;
	if (!Data.init())
		return 1 ;
	if (!Model.init())
		return 1 ;
	if (!Det.init())
		return 1 ;
	
	Model.intens_norm /= Data.mean_count ;
	
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
		
		rms_change = Model.diff_intens() * View.num_pix / Model.intens_norm * Data.mean_count ;
		Model.print_intens() ;
		
		time(&t2) ;
		
		fp.open("EMC.log", ios::out|ios::app) ;
		fp << i << "\t" << rms_change << "\t" << difftime(t2, t1) << endl ;
		fp.close() ;
	}
	
	free_mem(&Data, &Det, &Model, &View) ;
	
	return 0 ;
}


void free_mem(data_t *data, det_t *det, model_t *model, view_t *view) {
	int i, d, r ;
	
	delete[] det->pix ;
	
	for (i = 0 ; i < model->size ; ++i) {
		delete[] model->in[i] ;
		delete[] model->out[i] ;
		delete[] model->weight[i] ;
	}
	delete[] model->in ;
	delete[] model->out ;
	delete[] model->weight ;
	
	for (d = 0 ; d < data->num_data ; ++d) {
		delete[] data->frame[d].place_ones ;
		delete[] data->frame[d].place_multi ;
		delete[] data->frame[d].count ;
	}
	delete[] data->frame ;
	
	for (r = 0 ; r < view->num_rot ; ++r) {
		delete[] view->in[r] ;
		delete[] view->out[r] ;
	}
	delete[] view->in ;
	delete[] view->out ;
}
