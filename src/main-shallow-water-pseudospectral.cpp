#define NDEBUG

#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <errno.h>
#include <assert.h>
#include <unistd.h> // getopt

#include "configuration.hpp"
#include "fieldio.hpp"
#include "vector_opeartion.hpp"


using namespace std;
using namespace VORT_SRC_READER;

float dx, dy, Lx, Ly;
float *vort, *u, *v, *dvortdx, *dvortdy, *dvortdt, *workspace, *vort_src;

string vort_src_filename = "";
int has_vort_src = 0;
enum RECIPE_TYPE recipe_type = EMPTY;

fftwf_plan p_fwd_vort,    p_bwd_vort,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_fwd_dvortdt, p_bwd_psi;

fftwf_complex *vort_c0, *vort_c, *lvort_c, *tmp_c, *psi_c, *rk1_c, *rk2_c, *rk3_c, *copy_for_c2r;

fftwf_operation<XPTS,YPTS> fop(LX, LY);

VectorOperation<GRIDS> vop;

char filename[256];

void fftwf_backward_normalize(float *data) {
	for(int i=0; i < GRIDS; ++i) {
		data[i] /= GRIDS;
	}
}

float sumSqr(fftwf_complex *c) {
	float strength = 0;
	for(int i=0; i < HALF_GRIDS;++i){
		strength += pow(c[i][0],2) + pow(c[i][1],2);
	}
	return strength;
}

void print_spectrum(fftwf_complex *in) {
	for(int i=0; i < XPTS;++i){
		for(int j=0; j < HALF_YPTS;++j){
			printf("(%+6.2f, %+6.2f) ", in[HIDX(i,j)][0], in[HIDX(i,j)][1]);
		}
		printf("\n");
	}
}

void print_error(char * str) {
	printf("Error: %s\n", str);
}


int main(int argc, char* args[]) {
	char opt;

	while ((opt = getopt(argc, args, "I:O:i:s:f:")) != EOF) {
		switch(opt) {
			case 'I':
				input = optarg;
				break;
			case 'O':
				output = optarg;
				break;
			case 'i':
				init_file = optarg;
				break;
			case 's':
				vort_src_filename = optarg;
				recipe_type = SCRIPT;
				break;
			case 'f':
				vort_src_filename = optarg;
				recipe_type = FIFO;
				break;
		}
	}


	printf("##### Model setting #####\n");
	printf("Initial file          : %s \n", init_file.c_str());
	printf("Input folder          : %s \n", input.c_str());
	printf("Output folder         : %s \n", output.c_str());
	printf("Length X              : %.3f [m]\n", LX);
	printf("Length Y              : %.3f [m]\n", LY);
	printf("Spatial Resolution dx : %.3f [m]\n", dx);
	printf("Spatial Resolution dy : %.3f [m]\n", dy);
	printf("Time Resolution dt    : %.3f [s]\n", dt);
	printf("#########################\n\n\n");
	
	printf("Start project.\n");
	
	// open log
	FILE *log_fd = fopen("log", "w");
	if(log_fd == NULL) {
		perror("Open log file");
	}

	// initiate variables
	
	pv        = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	div       = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	h         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	f         = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	u_rot     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v_rot     = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	u_div     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v_div     = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	
	u         = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	v         = (float*) fftwf_malloc(sizeof(float) * GRIDS);


	// initializing plan
	p_fwd_vort       = fftwf_plan_dft_r2c_2d(XPTS, YPTS, vort, vort_c, FFTW_ESTIMATE);
	p_fwd_dvortdt    = fftwf_plan_dft_r2c_2d(XPTS, YPTS, dvortdt, tmp_c, FFTW_ESTIMATE);

	p_bwd_vort       = fftwf_plan_dft_c2r_2d(XPTS, YPTS, vort_c, vort, FFTW_ESTIMATE);
	p_bwd_dvortdx    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdx, FFTW_ESTIMATE);
	p_bwd_dvortdy    = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, dvortdy, FFTW_ESTIMATE);
	p_bwd_u          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, u, FFTW_ESTIMATE);
	p_bwd_v          = fftwf_plan_dft_c2r_2d(XPTS, YPTS, tmp_c, v, FFTW_ESTIMATE);

	p_bwd_psi        = fftwf_plan_dft_c2r_2d(XPTS, YPTS, psi_c, workspace, FFTW_ESTIMATE);


	dvortdx   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdy   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	dvortdt   = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	workspace = (float*) fftwf_malloc(sizeof(float) * GRIDS);
	vort_src  = (float*) fftwf_malloc(sizeof(float) * GRIDS);

	// complex numbers
	vort_c0   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	vort_c    = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	lvort_c   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);  // laplacian vorticity complex
	tmp_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	psi_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk1_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk2_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	rk3_c     = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);
	copy_for_c2r = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * HALF_GRIDS);

	// read input
	Lx = LX;
	Ly = LY;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

	sprintf(filename, "%s/%s", input.c_str(), init_file.c_str());
	readField(filename, vort, GRIDS);

	auto getDvortdt = [&](bool debug, int step){
		
		// 1. Take laplacian (diffusion terms)
		fop.laplacian(vort_c, lvort_c); 
		fop.laplacian(divg_c, ldivg_c); 
		fop.laplacian(h_c,       lh_c);

		// 2. Invert divergent and rotational flow
		fop.invertLaplacian(vort_c, psi_c);
		fop.invertLaplacian(divg_c, chi_c);
		
		// - rotation flow
		// -- calculate u_rot
		fop.grady(psi_c, tmp_c);
		fftwf_execute(p_bwd_u_vort); fftwf_backward_normalize(u_vort);
		for(int i=0; i<GRIDS;++i) { u_vort[i] = -u_vort[i]; }

		// -- calculate v_rot
		fop.gradx(psi_c, tmp_c);
		fftwf_execute(p_bwd_v_vort); fftwf_backward_normalize(v_vort);

		// - divergent flow
		// -- calculate u_div
		fop.gradx(chi_c, tmp_c);
		fftwf_execute(p_bwd_u_divg); fftwf_backward_normalize(u_divg);

		// -- calculate v_div
		fop.grady(chi_c, tmp_c);
		fftwf_execute(p_bwd_v_divg); fftwf_backward_normalize(v_divg);

		// - add together
		vop.add(u, u_vort, u_divg);
		vop.add(v, v_vort, v_divg);

		// 3. Calculate nonlinear multiplication in physical space

		// - invert vort, h
		fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
		fftwf_execute(p_bwd_h   ); fftwf_backward_normalize(h   );

		// - vort to absvort
		vop.add(absvort, vort, bg_vort);

		// - calculate multiplication
		vop.mul(absvort_u, absvort, u);
		vop.mul(absvort_v, absvort, v);

		vop.mul(h_u,  h,  u);
		vop.mul(h_v,  h,  v);

		vop.mul(v2,   v,     v);
		vop.mul(u2,   u,     u);
		vop.add( K , u2,    v2);
		vop.add( E ,  K,  geop);
		
		// - forward transformation
		fftwf_execute(p_fwd_absvort_u);
		fftwf_execute(p_fwd_absvort_v);
		fftwf_execute(p_fwd_h_u);
		fftwf_execute(p_fwd_h_v);
		fftwf_execute(p_fwd_E);

		// - get spectral derivative

		// -- [vort]
		vop.set(dvortdt_c, 0f);

		// --- diffusion, Rayleigh friction
		vop.mul(tmp_c, lvort_c, nu);   vop.isub(dvortdt_c, tmp_c);
		vop.mul(tmp_c,  vort_c, mu);   vop.iadd(dvortdt_c, tmp_c);

		// --- rest
		fop.gradx(absvort_u_c, tmp_c); vop.isub(dvortdt_c, tmp_c);
		fop.grady(absvort_v_c, tmp_c); vop.isub(dvortdt_c, tmp_c);

		// -- [divg]
		vop.set(ddivgdt_c, 0f);

		// --- diffusion, Rayleigh friction
		vop.mul(tmp_c, ldivg_c, nu);   vop.isub(ddivgdt_c, tmp_c);
		vop.mul(tmp_c,  divg_c, mu);   vop.iadd(ddivgdt_c, tmp_c);

		// --- rest
		fop.gradx(absvort_v_c, tmp_c); vop.iadd(ddivgdt_c, tmp_c);
		fop.grady(absvort_u_c, tmp_c); vop.isub(ddivgdt_c, tmp_c);
		fop.laplacian(E_c, tmp_c);     vop.isub(ddivgdt_c, tmp_c);

		// -- [h]
		vop.set(dhdt_c, 0f);

		// --- Source
		// vop.isub(dh_c, Q_c);

		// --- diffusion
		vop.mul(tmp_c, lh_c, nu);   vop.isub(dhdt_c, tmp_c);

		// --- rest
		fop.gradx(h_u_c, tmp_c);       vop.isub(dhdt_c, tmp_c);
		fop.grady(h_v_c, tmp_c);       vop.isub(dhdt_c, tmp_c);

		/*
		#ifdef OUTPUT_GRAD_VORT
		if(debug) {
			sprintf(filename, "%s/dvortdx_step_%d.bin", output.c_str(), step);
			writeField(filename, dvortdx, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);
		}
		#endif
		*/

	};

	auto RK4_evolve_single = [&](fftwf_complex * updated, fftw_complex * updated0, fftwf_complex * rk4_term, float coe) {
		for(int i=0; i<HALF_GRIDS; ++i) {
			updated[i][0] = updated0[i][0] + rk4_term[k][i][0] * coe;
			updated[i][1] = updated0[i][1] + rk4_term[k][i][1] * coe;
		}
	}

	auto RK4_last_step = [&](fftwf_complex * updated, fftw_complex * updated0, fftwf_complex** rk4_term) {
		for(int i=0; i<HALF_GRIDS; ++i) {
			updated[i][0] = updated0[i][0] + (rk4_term[0][i][0] + 2.0f * rk4_term[1][i][0] + 2.0f * rk4_term[2][i][0] + rk4_term[3][i][0]) / 6.0f;
			updated[i][1] = updated0[i][1] + (rk4_term[0][i][1] + 2.0f * rk4_term[1][i][1] + 2.0f * rk4_term[2][i][1] + rk4_term[3][i][1]) / 6.0f;
		}
	}



	float RK4_step_coe[3] = {0.5f, 0.5f, 1.0f};
	auto RK4_run = [&]() {
		memcpy(vort_c0, vort_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(divg_c0, divg_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup
		memcpy(h_c0,       h_c, sizeof(fftwf_complex) * HALF_GRIDS); // backup

		for(int k = 0; k < 4 ; ++k) {
			getDvortdt();
			vop.mul(rk4_vort_c[k], dvortdt_c, dt) 
			vop.mul(rk4_divg_c[k], ddivgdt_c, dt) 
			vop.mul(rk4_h_c[k],    dhdt_c,    dt)
 
			if(k==3) { continue; }

			evolve_single(vort_c, vort_c0, rk4_vort_c[k], RK4_step_coe[k]);
			evolve_single(divg_c, divg_c0, rk4_divg_c[k], RK4_step_coe[k]);
			evolve_single(h_c,       h_c0, rk4_h_c[k],    RK4_step_coe[k]);
		}

		RK4_last_step(vort_c, vort_c0, rk4_vort_c);
		RK4_last_step(divg_c, divg_c0, rk4_divg_c);
		RK4_last_step(h_c,       h_c0, rk4_h_c);

	};

	/*
	auto output_field_from_complex = [&](fftw_complex * data, char * filename) {
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

		
	};
	*/

	printf("Program initialization complete.\n");

	// Preparation: 
	fftwf_execute(p_fwd_vort);
	fftwf_execute(p_fwd_divg);
	fftwf_execute(p_fwd_h);
	fftwf_execute(p_fwd_Q);

	int record_flag = 0;
	for(int step = 0; step < total_steps; ++step) {

		printf("# Step %d, time = %.2f", step, step * dt);
		if( (record_flag = ((step % record_step) == 0)) ) { printf(", record now!");}
		printf("\n");

		if(record_flag) {

			// Output vort
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, vort_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_vort); fftwf_backward_normalize(vort);
			sprintf(filename, "%s/vort_step_%d.bin", output.c_str(), step);
			writeField(filename, vort, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(vort_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

			
			// Output divg
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, divg_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_divg); fftwf_backward_normalize(divg);
			sprintf(filename, "%s/divg_step_%d.bin", output.c_str(), step);
			writeField(filename, divg, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(divg_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);
	
			// Output h
			// backup vort_c because c2r must destroy input (NO!!!!!)
			memcpy(copy_for_c2r, h_c, sizeof(fftwf_complex) * HALF_GRIDS);

			fftwf_execute(p_bwd_h); fftwf_backward_normalize(h);
			sprintf(filename, "%s/h_step_%d.bin", output.c_str(), step);
			writeField(filename, h, GRIDS);
			fprintf(log_fd, "%s\n", filename); fflush(log_fd);

			// restore vort_c because c2r must destroy input (NO!!!!!)
			memcpy(h_c, copy_for_c2r, sizeof(fftwf_complex) * HALF_GRIDS);

		}


		RK4_run();

		//fftwf_destroy_plan(p);
		//fftwf_free(in); fftwf_free(out);
	}

	fclose(log_fd);
	printf("Program ends. Congrats!\n");
	return 0;
}
