/* Reconstruction of mask pattern from randomly oriented patterns with low photon count

*/

#include <iomanip>
#include <time.h>
#include <math.h>
#include "data.hh"
#include "det.hh"
#include "model.hh"
#include "view.hh"
using namespace std ;

void expand(det_t*, model_t*, view_t*) ;
void maximize(data_t*, view_t*) ;
void compress(det_t*, model_t*, view_t*) ;

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
	fp << "size = " << Model.size << endl ;
	fp << "mean_count = " << Data.mean_count << "\tintens_norm = " << Model.intens_norm << endl << endl ; 
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
		fp << i << "\t" << rms_change << "  \t" << difftime(t2, t1) << " s" << endl ;
		fp.close() ;
	}
	
	return 0 ;
}
