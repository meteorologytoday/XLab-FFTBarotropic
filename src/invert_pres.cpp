/*
 * Author: Hsu, Tien-Yiao
 *
 * Description: 
 *
 * This program uses the psi output from main program to
 * invert pressure anomaly.
 *
 */

#define NDEBUG

#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fftw3.h>
#include <errno.h>
#include <assert.h>

#include "configuration.hpp"
#include "fftwfop.cpp"
#include "fieldio.hpp"


float *pres, *psi, *dpsidx2, *dpsidy2, *dpsidxdy, *gaus_curv;

fftwf_complex *tmp_c, *psi_c, *dpsidx2_c, *dpsidy2_c, *dpsidxdy_c, *lap_pres_c;
fftwf_plan p_fwd_psi,      p_fwd_gaus_curv2lap_pres_c,
		   p_bwd_dpsidx2,  p_bwd_dpsidy2,
		   p_bwd_dpsidxdy, p_bwd_tmp2pres;

fftwf_operation<XPTS,YPTS> fop(LX, LY);

void fftwf_backward_normalize(float *data) {
	for(int i=0; i < GRIDS; ++i) {
		data[i] /= GRIDS;
	}
}

void trim(char * str) {
	char * p = str + strlen(str) - 1;
	while(p >= str) {
		if(strcmp(p,"\n") == 0) {
			*p = '\0';
		}
		--p;
	}
}

int main(int argc, char* args[]) {

	// initiate variables
	pres      = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	psi       = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dpsidx2   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dpsidy2   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dpsidxdy  = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	gaus_curv = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	// complex numbers
	tmp_c      = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	psi_c      = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	dpsidx2_c  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	dpsidy2_c  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	dpsidxdy_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	lap_pres_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);

	// initializing plan
	p_fwd_psi        = fftwf_plan_dft_r2c_2d(XPTS, YPTS, psi, psi_c, FFTW_ESTIMATE);
	p_fwd_gaus_curv2lap_pres_c = fftwf_plan_dft_r2c_2d(XPTS, YPTS, gaus_curv, lap_pres_c, FFTW_ESTIMATE);

	p_bwd_dpsidx2    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, dpsidx2_c, dpsidx2, FFTW_ESTIMATE);
	p_bwd_dpsidy2    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, dpsidy2_c, dpsidy2, FFTW_ESTIMATE);
	p_bwd_dpsidxdy   = fftwf_plan_dft_c2r_2d(XPTS, YPTS, dpsidxdy_c, dpsidxdy, FFTW_ESTIMATE);

	p_bwd_tmp2pres   = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, pres, FFTW_ESTIMATE);

	// read input
	char filename[1024], from_file[1024], to_file[1024];
	char * sep_beg;
	char sep[] = "=>";
	while(fgets(filename, 1024, stdin) != NULL) {
		trim(filename);	
		if( (sep_beg  = strstr (filename, sep)) != NULL ) {
			int l;

			l = sep_beg - filename;
			memcpy(from_file, filename, l);
			from_file[l] = '\0';

			l = strlen(filename) - strlen(from_file) - strlen(sep);
			memcpy(to_file, sep_beg + strlen(sep), l);
			to_file[l] = '\0';

		} else {
			printf("Error reading input: %s. Continue next line...\n", filename);
			continue;
		}

		printf("Ready to read\n");
		readField(from_file, psi, GRIDS);
		printf("Already read\n");


		fftwf_execute(p_fwd_psi);

		// ### gaussian curvature ###

		fop.gradx(psi_c, tmp_c);
		fop.gradx(tmp_c, dpsidx2_c);

		fop.grady(psi_c, tmp_c);
		fop.grady(tmp_c, dpsidy2_c);

		fop.gradx(tmp_c, dpsidxdy_c);

		
		fop.dealiase(dpsidx2_c,  dpsidx2_c);
		fop.dealiase(dpsidy2_c,  dpsidy2_c);
		fop.dealiase(dpsidxdy_c, dpsidxdy_c);

		// back to physical space
		fftwf_execute(p_bwd_dpsidx2);  fftwf_backward_normalize(dpsidx2);
		fftwf_execute(p_bwd_dpsidy2);  fftwf_backward_normalize(dpsidy2);
		fftwf_execute(p_bwd_dpsidxdy); fftwf_backward_normalize(dpsidxdy);



		for(int i=0; i<GRIDS;++i) { gaus_curv[i] = dpsidx2[i] * dpsidy2[i] - pow(dpsidxdy[i], 2.0f); }
		// ### Source term of lap_prec ### 
		fftwf_execute(p_fwd_gaus_curv2lap_pres_c);


		fop.laplacian(psi_c, tmp_c);

		for(int i=0; i<HALF_GRIDS;++i) {
			 lap_pres_c[i][0] = rho * ( f * tmp_c[i][0] + 2.0 * lap_pres_c[i][0] );
			 lap_pres_c[i][1] = rho * ( f * tmp_c[i][1] + 2.0 * lap_pres_c[i][1] );
		}
		
		fop.invertLaplacian(lap_pres_c, tmp_c);

		fftwf_execute(p_bwd_tmp2pres); fftwf_backward_normalize(pres);
		writeField(to_file, pres, GRIDS);
	}

	printf("Program ends. Congrats!\n");
	return 0;
}
