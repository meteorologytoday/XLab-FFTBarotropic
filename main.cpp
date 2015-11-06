#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include <fftw3.h>
#include "configuration.h"

#define M_2PI (2.0*M_PI)

float dx, dy, Lx, Ly;


float *vort, *u, *v, *dvortdx, *dvortdy, *dvortdt;
fftwf_plan p_fwd_vort,    p_bwd_vort,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_fwd_dvortdt;

fftwf_complex *vort_c0, *vort_c, *tmp_c, *psi_c, *rk1_c, *rk2_c, *rk3_c, *rk4_c;

float *gradx_coe, *grady_coe, *laplacian_coe;

int total_steps;
float dt;

int main(){

	// initiate variables
	vort    = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	u       = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v       = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdx = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdy = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdt = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	// complex numbers
	vort_c0   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	vort_c    = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	tmp_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	psi_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk1_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk2_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk3_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk4_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);

	// initializing plan
	p_fwd_vort = fftwf_plan_dft_r2c_2d(XPTS, YPTS, vort, vort_c, FFTW_FORWARD);
	p_fwd_dvortdt = fftwf_plan_dft_r2c_2d(XPTS, YPTS, dvortdt, tmp_c, FFTW_FORWARD);

	p_bwd_vort       = fftwf_plan_dft_c2r_2d(XPTS, YPTS, vort_c, vort, FFTW_BACKWARD);
	p_bwd_dvortdx    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdx, FFTW_BACKWARD);
	p_bwd_dvortdy    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdy, FFTW_BACKWARD);
	p_bwd_u          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, u, FFTW_BACKWARD);
	p_bwd_v          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, v, FFTW_BACKWARD);

	// read input
	Lx = 600000.0;
	Ly = 600000.0;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

	float centerx = Lx / 2.0, centery = Ly / 2.0, epsilon = 0.7, lambda = 2.0, zeta0 = 5e-3, r_i = 30000.0, r_o = 60000.0;
	float r, r_i_alpha, r_o_alpha, r_prime;

	auto radius = [centerx, centery](float x, float y) -> float {
		return sqrtf(pow(x-centerx,2) + pow(y-centery,2));
	};
	auto alpha = [centerx, centery, epsilon, radius](float x, float y) -> float {
		float c = (y-centery) / radius(x,y);
		return sqrtf((1.0 - pow(epsilon,2)) / (1.0 - pow(epsilon*c,2)));
	};

	int x, y;
	for(int i=0; i<XPTS; ++i) {
		x = i * dx;
		for(int j=0; j<YPTS; ++j) {
			y = j * dy;
			r = radius(x,y);
			r_i_alpha = r_i * alpha(x,y);
			r_o_alpha = r_o * alpha(x,y);

			if(r <= r_i_alpha) {
				vort[IDX(i,j)] = zeta0;
			} else if (r <= r_o_alpha) {
				r_prime = (r - r_i_alpha)/(r_o_alpha - r_i_alpha);
				vort[IDX(i,j)] = zeta0 * (1.0 - exp( - lambda / r_prime * exp(1.0 / (r_prime - 1))));
			} else {
				vort[IDX(i,j)] = 0;
			}
		}
	}

	// setup gradx_coe, grady_coe, laplacian_coe
	gradx_coe = (float*) fftwf_malloc(sizeof(float) * HALF_XPTS);
	grady_coe = (float*) fftwf_malloc(sizeof(float) * HALF_YPTS);
	laplacian_coe = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);

	for(int i=0; i<HALF_XPTS; ++i) {
		gradx_coe[i] = M_2PI * ((float) i) / Lx;
	}

	for(int j=0; j<HALF_YPTS; ++j) {
			grady_coe[j] = M_2PI * ((float) j) / Ly;
	}

	for(int i=0; i<HALF_XPTS; ++i) {
		for(int j=0; j<HALF_YPTS; ++j) {
			laplacian_coe[IDX(i,j)] = - (pow(gradx_coe[i],2) + pow(grady_coe[j],2));
		}
	}
	laplacian_coe[IDX(0,0)] = 1.0; // this coe is special

	// step 01
	fftwf_execute(p_fwd_vort);
	memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

	// step 02 Everythings ready!!
	for(int step = 0; step < total_steps; ++step) {

		// step 03 take dvortdx, save as tmp_c
		for(int j=0; j<HALF_YPTS; ++j) {
			for(int i=0; i<HALF_XPTS; ++i) {
				tmp_c[IDX(i,j)][0] = - vort_c[IDX(i,j)][1] * gradx_coe[i];
				tmp_c[IDX(i,j)][1] = vort_c[IDX(i,j)][0] * gradx_coe[i];
			}
		}

		// step 04
		fftwf_execute(p_bwd_dvortdx);

		// step 05 take dvortdy, save as tmp_c
		for(int j=0; j<HALF_YPTS; ++j) {
			for(int i=0; i<HALF_XPTS; ++i) {
				tmp_c[IDX(i,j)][0] = - vort_c[IDX(i,j)][1] * grady_coe[j];
				tmp_c[IDX(i,j)][1] = vort_c[IDX(i,j)][0] * grady_coe[j];
			}
		}

		// step 06
		fftwf_execute(p_bwd_dvortdy);

		// step 07 get psi_c
		for(int j=0; j<HALF_YPTS; ++j) {
			for(int i=0; i<HALF_XPTS; ++i) {
				psi_c[IDX(i,j)][0] = - vort_c[IDX(i,j)][0] / laplacian_coe[IDX(i,j)];
				psi_c[IDX(i,j)][1] = - vort_c[IDX(i,j)][1] / laplacian_coe[IDX(i,j)];
			}
		}

		//fftwf_destroy_plan(p);
		//fftwf_free(in); fftwf_free(out);
	}







	return 0;
}
