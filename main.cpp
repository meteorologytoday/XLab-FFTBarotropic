#include <cstdio>
#include <iostream>
#include <fftw3.h>
#include "configuration.h"

float* psi, vort, u, v;
fftw_plan *p_fwd, *p_bwd;

float *in;
fftwf_complex *out;
fftwf_plan p;

int total_steps;
float dt;

int main(){

	// read input

	// initiate variables
	psi = new float[GRIDS]();
	u   = new float[GRIDS]();

	p = fftwf_plan_dft_r2c_2d(XPTS, YPTS, in, out, FFTW_FORWARD);

    in = (float*) fftwf_malloc(sizeof(float) * GRIDS);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * GRIDS);

    fftwf_execute(p); /* repeat as needed */

    fftwf_destroy_plan(p);
    fftwf_free(in); fftwf_free(out);

	for(int step = 0; step < total_steps; ++step) {
		// invert psi by fourier transform

		//

	}








	return 0;
}
