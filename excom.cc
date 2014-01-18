#ifndef EXCOM_HPP
#define EXCOM_HPP

#include <math.h>
#include "structs.hpp"

#define PI 3.14159265


void make_rot(double a, double r_matrix[][2]) {
	r_matrix[0][0] = cos(a) ;
	r_matrix[0][1] = sin(a) ;
	r_matrix[1][0] = -sin(a);
	r_matrix[1][1] = cos(a) ;
}


void expand(det_t *det, model_t *model, view_t *view) {
	int r, t, i, j ;
	double rot_pix[2] ;
	double rot[2][2] ;
	int x, y, q_max = (model->size - 3) / 2 ;
	double tx, ty, fx, fy, cx, cy ;
	
	for (r = 0 ; r < view->num_rot ; ++r) {
		make_rot(r * 2 * PI / view->num_rot, rot) ;
		
		for (t = 0 ; t < view->num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*det->pix[t][j] ;
			}
			
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			view->in[r][t] = cx*cy*model->in[x][y] + cx*fy*model->in[x][y+1] + fx*cy*model->in[x+1][y] + fx*fy*model->in[x+1][y+1] ;
		}
	}
}


void compress(det_t *det, model_t *model, view_t *view) {
	int r, t, i, j ;
	double rot_pix[2] ;
	double rot[2][2] ;
	int x, y, q_max = (model->size - 3) / 2 ;
	double tx, ty, fx, fy, cx, cy ;
	double w, f ;
	
	for (x = 0 ; x < model->size ; ++x)
	for (y = 0 ; y < model->size ; ++y) {
		model->weight[x][y] = 0. ;
		model->out[x][y] = 0. ;
	}
	
	for (r = 0 ; r < view->num_rot ; ++r) {
//		if (s[r] == 0.)
//			continue ;
			
		make_rot(r * 2 * PI / view->num_rot, rot) ;
		
		for (t = 0 ; t < view->num_pix ; ++t) {
			for (i = 0 ; i < 2 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 2 ; ++j)
					rot_pix[i] += rot[i][j]*det->pix[t][j] ;
			}
			
			tx = rot_pix[0] + q_max ;
			ty = rot_pix[1] + q_max ;
			
			x = tx ;
			y = ty ;
			
			fx = tx - x ;
			fy = ty - y ;
			
			cx = 1. - fx ;
			cy = 1. - fy ;
			
			w = view->out[r][t] ;
			
			f = cx*cy ;
			model->weight[x][y] += f ;
			model->out[x][y] += f * w ;
			
			f = cx*fy ;
			model->weight[x][y+1] += f ;
			model->out[x][y+1] += f * w ;
			
			f = fx*cy ;
			model->weight[x+1][y] += f ;
			model->out[x+1][y] += f * w ;
			
			f = fx*fy ;
			model->weight[x+1][y+1] += f ;
			model->out[x+1][y+1] += f * w ;
		}
	}
	
	for (x = 0 ; x < model->size ; ++x)
	for (y = 0 ; y < model->size ; ++y)
		if (model->weight[x][y] != 0.)
			model->out[x][y] *= model->intens_norm / model->weight[x][y] ;
}

#endif // EXCOM_HPP
