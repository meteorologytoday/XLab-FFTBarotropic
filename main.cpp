#define NDEBUG

#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fftw3.h>
#include "configuration.h"
#include "fftwf_operation.h"

#include <assert.h>

float dx, dy, Lx, Ly;


float *vort, *u, *v, *dvortdx, *dvortdy, *dvortdt, *workspace;
fftwf_plan p_fwd_vort,    p_bwd_vort,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_fwd_dvortdt, p_bwd_psi;

fftwf_complex *vort_c0, *vort_c, *lvort_c, *tmp_c, *psi_c, *rk1_c, *rk2_c, *rk3_c, *copy_for_c2r;

int total_steps = 100 * (30 / 5);
float dt = 3.0f;

fftwf_operation<XPTS,YPTS> fop(LX, LY);

char filename[256];

void writeField(const char * filename, float *data) {

	FILE * file = fopen(filename, "wb");
	int flg;
	if((flg = fwrite(data, sizeof(float), GRIDS, file)) < 0) {
		printf("ERROR! %d\n", flg);
	}
	fclose(file);

	printf("Output %s\n", filename);
}

void fftwf_backward_normalize(float *data) {
	for(int i=0; i<GRIDS; ++i) {
		data[i] /= GRIDS;
	}
}

float sumSqr(fftwf_complex *c) {
	float strength = 0;
	for(int i=0; i<HALF_GRIDS;++i){
		strength += pow(c[i][0],2) + pow(c[i][1],2);
	}
	return strength;
}

void print_spectrum(fftwf_complex *in) {
	for(int i=0; i<XPTS;++i){
		for(int j=0; j<HALF_YPTS;++j){
			printf("(%+6.2f, %+6.2f) ", in[HIDX(i,j)][0], in[HIDX(i,j)][1]);
		}
		printf("\n");
	}
}

int main(){
	printf("Start project.\n");
	std::cout << "C++ feature -- 1st time" << std::endl;
	// initiate variables
	vort      = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	u         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdx   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdy   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdt   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	workspace = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	// complex numbers
	vort_c0   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	vort_c    = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	lvort_c   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	tmp_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	psi_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk1_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk2_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk3_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	copy_for_c2r = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);

	// initializing plan
	p_fwd_vort       = fftwf_plan_dft_r2c_2d(XPTS, YPTS, vort, vort_c, FFTW_ESTIMATE);
	p_fwd_dvortdt    = fftwf_plan_dft_r2c_2d(XPTS, YPTS, dvortdt, tmp_c, FFTW_ESTIMATE);

	p_bwd_vort       = fftwf_plan_dft_c2r_2d(XPTS, YPTS, vort_c, vort, FFTW_ESTIMATE);
	p_bwd_dvortdx    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdx, FFTW_ESTIMATE);
	p_bwd_dvortdy    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdy, FFTW_ESTIMATE);
	p_bwd_u          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, u, FFTW_ESTIMATE);
	p_bwd_v          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, v, FFTW_ESTIMATE);

	p_bwd_psi        = fftwf_plan_dft_c2r_2d(XPTS, YPTS, psi_c, workspace, FFTW_ESTIMATE);

	// read input
	Lx = LX;
	Ly = LY;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

	float centerx = Lx / 2.0, centery = Ly / 2.0, epsilon = 0.7, lambda = 2.0, zeta0 = .005f, r_i = 30000.0, r_o = 60000.0;
	float r, r_i_alpha, r_o_alpha, r_prime;

	auto radius = [centerx, centery](float x, float y) -> float {
		return sqrtf(pow(x-centerx,2) + pow(y-centery,2));
	};
	auto alpha = [centerx, centery, epsilon, radius](float x, float y) -> float {
		float c;
		float r = radius(x,y);
		if(r == 0.0f) {
			c = 0;
		} else {
			c = (y-centery) / radius(x,y);
		}
		return sqrtf((1.0 - pow(epsilon,2)) / (1.0 - pow(epsilon*c,2)));
	};

	float x, y;
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
	sprintf(filename, "initial.bin");
	writeField(filename, vort);

	auto getDvortdt = [&](bool debug, int step){
		// step ?? take lvort_c
		fop.laplacian(vort_c, lvort_c);

		// step 03 take dvortdx, save as tmp_c
		fop.gradx(vort_c, tmp_c);

		// step 04
		fftwf_execute(p_bwd_dvortdx); fftwf_backward_normalize(dvortdx);
		if(debug) {
			sprintf(filename, "dvortdx_step_%d.bin", step);
			writeField(filename, dvortdx);
		}

		// step 05 take dvortdy, save as tmp_c
		fop.grady(vort_c, tmp_c);

		// step 06
		fftwf_execute(p_bwd_dvortdy); fftwf_backward_normalize(dvortdy);

		if(debug) {
			sprintf(filename, "dvortdy_step_%d.bin", step);
			writeField(filename, dvortdy);
		}

		// step 07 get psi_c
		fop.invertLaplacian(vort_c, psi_c);

		if(debug) {
			 // backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, psi_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_psi); fftwf_backward_normalize(workspace);
			sprintf(filename, "psi_step_%d.bin", step);
			writeField(filename, workspace);

			 // restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(psi_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
		}

		// step 08 get u_c
		fop.grady(psi_c, tmp_c);
		// step 09
		fftwf_execute(p_bwd_u); fftwf_backward_normalize(u);
		for(int i=0; i<GRIDS;++i) { u[i] = -u[i]; }

		if(debug) {
			sprintf(filename, "u_step_%d.bin", step);
			writeField(filename, u);
		}

		// step 10 get v_c
		fop.gradx(psi_c, tmp_c);
		// step 11
		fftwf_execute(p_bwd_v); fftwf_backward_normalize(v);

		if(debug) {
			sprintf(filename, "v_step_%d.bin", step);
			writeField(filename, v);
		}

		// step 12 get dvortdt
		for(int i=0; i<GRIDS;++i) {
			dvortdt[i] = - u[i] * dvortdx[i] - v[i] * dvortdy[i];
		}

		if(debug) {
			sprintf(filename, "dvortdt_step_%d.bin", step);
			writeField(filename, dvortdt);
		}

		// step 13 get dvortdt_c and save in tmp_c
		fftwf_execute(p_fwd_dvortdt);

		// step ?? add laplacian term
		for(int i=0; i<HALF_GRIDS;++i) {
			tmp_c[i][0] += lvort_c[i][0] * NU;
			tmp_c[i][1] += lvort_c[i][1] * NU;
		}
	};

	auto evolve = [&](fftwf_complex *rk, float dt) {
		for(int i=0; i<HALF_GRIDS; ++i){
			vort_c[i][0] = vort_c0[i][0] + rk[i][0] * dt;
			vort_c[i][1] = vort_c0[i][1] + rk[i][1] * dt;
		}
	};

	printf("Initialization complete.\n");

	// step 01
	fftwf_execute(p_fwd_vort);

	// step 02 Everything is ready!!
	for(int step = 0; step < total_steps; ++step) {
		printf("# Step %d\n", step+1);

		memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

		for(int k = 0 ; k < 4; ++k) {

			getDvortdt((step % 100 == 0) && k==0, step);

			// step 14+15 dealiasing dvortdt_c and save to rk?_c)
			// DEPENDS ON RK?
			switch(k) {
				case 0:
					fop.dealiase(tmp_c, rk1_c);	evolve(rk1_c, dt / 2.0f);
					break;
				case 1:
					fop.dealiase(tmp_c, rk2_c);	evolve(rk2_c, dt / 2.0f);
					break;
				case 2:
					fop.dealiase(tmp_c, rk3_c);	evolve(rk3_c, dt);
					break;
				case 3:
					// Actually the variable rk4_c is not needed because this is the last call
					// we can simply replace rk4_c by tmp_c.
					fop.dealiase(tmp_c, tmp_c);

					// step 24: get new vort_c
					for(int i=0; i<HALF_GRIDS; ++i){
						vort_c[i][0] = vort_c0[i][0] + (rk1_c[i][0] + 2.0f * rk2_c[i][0] + 2.0f * rk3_c[i][0] + tmp_c[i][0]) * dt / 6.0f;
						vort_c[i][1] = vort_c0[i][1] + (rk1_c[i][1] + 2.0f * rk2_c[i][1] + 2.0f * rk3_c[i][1] + tmp_c[i][1]) * dt / 6.0f;
					}

					break;

			}
		}

		if((step+1) % 100 == 0) { // every 5 min

			 // backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "vort_%d.bin", (step+1) / 20);
			writeField(filename, vort);

			 // restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
		}


		//fftwf_destroy_plan(p);
		//fftwf_free(in); fftwf_free(out);
	}
	printf("Program ends. Congrats!\n");
	return 0;
}
