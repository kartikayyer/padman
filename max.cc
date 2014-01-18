#ifndef MAX_HPP
#define MAX_HPP

#include <math.h>
#include "structs.hpp"

#define P_MIN 0.001


void maximize(data_t *data, view_t *view) {
	int d, r, t ;
	double max_exp, p_sum, rot_intens, total_intens = 0., rescale, rescaled_view ;
	frame_t *frame_d ;
	double *w_in_r, *w_out_r ;
	double *s, *u, *p ;
	
	s = new double[view->num_rot] ;
	u = new double[view->num_rot] ;
	p = new double[view->num_rot] ;
	
	for (r = 0 ; r < view->num_rot ; ++r) {
		rot_intens = 0. ;
		for (t = 0 ; t < view->num_pix ; ++t)
			rot_intens += view->in[r][t] ;
			
		total_intens += rot_intens / view->num_rot ;
	}
	
	rescale = data->mean_count / total_intens ;
	
	for (r = 0 ; r < view->num_rot ; ++r) {
		s[r] = 0. ;
		u[r] = 0. ;
		
		for (t = 0 ; t < view->num_pix ; ++t) {
			rescaled_view = view->in[r][t] * rescale ;
			u[r] -= rescaled_view ;
			view->in[r][t] = log(rescaled_view) ;
			
			view->out[r][t] = 0. ;
		}
	}
	
    	for (d = 0 ; d < data->num_data ; ++d) {
		frame_d = &(data->frame[d]) ;
		
		max_exp = -100. * view->num_pix ;
		for (r = 0 ; r < view->num_rot ; ++r) {
			p[r] = u[r] ;
			
			w_in_r = view->in[r] ;
			
			for (t = 0 ; t < frame_d->ones ; ++t)
				p[r] += w_in_r[frame_d->place_ones[t]] ;
			
			for (t = 0 ; t < frame_d->multi ; ++t)
				p[r] += frame_d->count[t] * w_in_r[frame_d->place_multi[t]] ;
			
			if (p[r] > max_exp)
				max_exp = p[r] ;
		}
		
		p_sum = 0. ;
		for (r = 0 ; r < view->num_rot ; ++r) {
			p[r] = exp(p[r] - max_exp) ;
			p_sum += p[r] ;
		}
		
        	for (r = 0 ; r < view->num_rot ; ++r) {
			p[r] /= p_sum ;
			
			if (p[r] < P_MIN)
				continue ;
				
			s[r] += p[r] ;
				
			w_out_r = view->out[r] ;
			
			for (t = 0 ; t < frame_d->ones ; ++t)
				w_out_r[frame_d->place_ones[t]] += p[r] ;
			
			for (t = 0 ; t < frame_d->multi ; ++t)
				w_out_r[frame_d->place_multi[t]] += p[r] * frame_d->count[t] ;
		}
	}
	
    	for (r = 0 ; r < view->num_rot ; ++r) {
		if (s[r] == 0.)
			continue ;
			
		for (t = 0 ; t < view->num_pix ; ++t)
			view->out[r][t] /= s[r] ;
	}
	
	delete[] s ;
	delete[] u ;
	delete[] p ;
}

#endif // MAX_H
