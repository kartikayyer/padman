#ifndef STRUCTS_H
#define STRUCTS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std ;

struct frame_t {
	int ones, multi ;
	int *place_ones, *place_multi, *count ;
} ;

struct data_t {
	frame_t *frame ;
	double mean_count ;
	int num_data ;
} Data ;

struct det_t {
	double **pix ;
	int num_pix ;
} Det ;

struct model_t {
	double **in, **out, **weight ;
	double intens_norm ;
	int size ;
} Model ;

struct view_t {
	double **in, **out ;
	int num_pix, num_rot ;
} View ;

void expand(det_t*, model_t*, view_t*) ;
void maximize(data_t*, view_t*) ;
void compress(det_t*, model_t*, view_t*) ;

int setup(data_t *data, det_t *det, model_t *model, view_t *view) {
	fstream fp ;
	int i, j, k, d, r, t ;
	int rand_start = 0 ;
	
	view->num_rot = 1000 ;
	model->intens_norm = 2000. ; 		// Chosen to make the image easily visible on a scale of 0-1
	
	fp.open("det.dat", ios::in) ;
	if (!fp.is_open()) {
		cerr << "cannot open det.dat\n" ;
		return 0 ;
	}
	
	fp >> view->num_pix >> model->size ;
	model->size += 1 ;
	det->num_pix = view->num_pix ;
	
	det->pix = new double*[det->num_pix] ;
	
	for (t = 0 ; t < det->num_pix ; ++t) {
		det->pix[t] = new double[2] ;
		for (i = 0 ; i < 2 ; ++i) 
			fp >> det->pix[t][i] ;
	}
	
	fp.close() ;
	
	model->in = new double*[model->size] ;
	model->out = new double*[model->size] ;
	model->weight = new double*[model->size] ;
	for (i = 0 ; i < model->size ; ++i) {
		model->in[i] = new double[model->size] ;
		model->out[i] = new double[model->size] ;
		model->weight[i] = new double[model->size] ;
		for (j = 0 ; j < model->size ; ++j)
			model->in[i][j] = ((double) rand()) / RAND_MAX ;
	}
	
	fp.open("photons.dat", ios::in) ;
	if (!fp.is_open()) {
		cerr << "cannot open photons.dat\n" ;
		return 0 ;
	}
	
	fp >> data->num_data ;
	
	data->frame = new frame_t[data->num_data] ;
	
	for (d = 0 ; d < data->num_data ; ++d) {
		fp >> data->frame[d].ones ;
		data->mean_count = data->frame[d].ones ;
		
		data->frame[d].place_ones = new int[data->frame[d].ones] ;
		
		for (t = 0 ; t < data->frame[d].ones ; ++t)
			fp >> data->frame[d].place_ones[t] ;
		
		fp >> data->frame[d].multi ;
		
		data->frame[d].place_multi = new int[data->frame[d].multi] ;
		data->frame[d].count = new int[data->frame[d].multi] ;
		
		for (t = 0 ; t < data->frame[d].multi ; ++t) {
			fp >> data->frame[d].place_multi[t] >> data->frame[d].count[t] ;
			data->mean_count += data->frame[d].count[t] ;
		}
	}
	
	data->mean_count /= data->num_data ;
	model->intens_norm /= data->mean_count ;
	
	fp.close() ;
	
	view->in = new double*[view->num_rot] ;
	view->out = new double*[view->num_rot] ;
	
	for (r = 0 ; r < view->num_rot ; ++r) {
		view->in[r] = new double[view->num_pix] ;
		view->out[r] = new double[view->num_pix] ;
	}
	
	return 1 ;
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

#endif // STRUCTS_H
