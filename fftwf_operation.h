#ifndef FFTWF_OPERATION_H
#define FFTWF_OPERATION_H
#include <fftw3.h>

#include <cmath>

const float TWOPI = (acos(-1.0f) * 2.0f);

template<int XPTS, int YPTS> class fftwf_operation {
private:
	float *gradx_coe, *grady_coe, *laplacian_coe, *laplacian_coe_inverse, *dealiasing_mask;
	int dealiase_xwavenumber, dealiase_ywavenumber;
	const int HALF_XPTS = (int)(XPTS/2) + 1,
			  HALF_YPTS = (int)(YPTS/2) + 1,
			  HALF_GRIDS = XPTS*HALF_YPTS;
public:
	fftwf_operation(float Lx, float Ly);
	~fftwf_operation();

	void gradx(fftwf_complex *in, fftwf_complex *out);
	void grady(fftwf_complex *in, fftwf_complex *out);
	void laplacian(fftwf_complex *in, fftwf_complex *out);
	void invertLaplacian(fftwf_complex *in, fftwf_complex *out);
	void dealiase(fftwf_complex *in, fftwf_complex *out);

	inline int reflectedXWavenumberIndex(int i) {return XPTS - i;};
	inline int HIDX  (int i, int j) {return HALF_YPTS*i + j;};
	inline int R_HIDX(int i, int j) {return HIDX(this->reflectedXWavenumberIndex(i),j);};
};



template<int XPTS, int YPTS> fftwf_operation<XPTS,YPTS>::fftwf_operation(float Lx, float Ly) {
	this->gradx_coe = (float*) fftwf_malloc(sizeof(float) * XPTS);
	this->grady_coe = (float*) fftwf_malloc(sizeof(float) * HALF_YPTS);
	this->laplacian_coe_inverse = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->laplacian_coe = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->dealiasing_mask = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->dealiase_xwavenumber = (int) ceil(((float)XPTS)/3.0);
	this->dealiase_ywavenumber = (int) ceil(((float)YPTS)/3.0);

	for(int i=0; i < this->HALF_XPTS; ++i) {
		this->gradx_coe[i] = TWOPI * ((float) i) / Lx;
		this->gradx_coe[this->reflectedXWavenumberIndex(i)] = - this->gradx_coe[i];
	}

	for(int j=0; j < this->HALF_YPTS; ++j) {
		this->grady_coe[j] = TWOPI * ((float) j) / Ly;
	}

	for(int i=0; i<this->HALF_XPTS; ++i) {
		for(int j=0; j<this->HALF_YPTS; ++j) {
			this->laplacian_coe_inverse[this->HIDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
			this->laplacian_coe_inverse[this->R_HIDX(i,j)] = this->laplacian_coe_inverse[this->HIDX(i,j)];

			this->laplacian_coe[this->HIDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
			this->laplacian_coe[this->R_HIDX(i,j)] = this->laplacian_coe[this->HIDX(i,j)];
		}
	}
	this->laplacian_coe_inverse[this->HIDX(0,0)] = 1.0; // this coe is special

	float generalized_wavenumber_square = pow(this->dealiase_xwavenumber,2) + pow(this->dealiase_ywavenumber,2);

	for(int i = 0; i < HALF_XPTS ; ++i) {
		for(int j = 0; j < HALF_YPTS; ++j) {
			this->dealiasing_mask[this->HIDX(i,j)] = (pow(i,2) + pow(j,2) >= generalized_wavenumber_square) ? 0.0f : 1.0f;
			this->dealiasing_mask[this->R_HIDX(i,j)] = this->dealiasing_mask[this->HIDX(i,j)];
		}
	}


/*	for(int i = 0; i < XPTS ; ++i) {
		for(int j = 0; j < HALF_YPTS; ++j) {
			printf("%d ", (int)(this->dealiasing_mask[this->HIDX(i,j)]));
		}
		printf("\n");
	}
	*/
}

template<int XPTS, int YPTS> fftwf_operation<XPTS,YPTS>::~fftwf_operation() {
	fftwf_free(this->gradx_coe);
	fftwf_free(this->grady_coe);
	fftwf_free(this->laplacian_coe);
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::gradx(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<this->HALF_YPTS; ++j) {
		for(int i=0; i<this->HALF_XPTS; ++i) {
			out[this->HIDX(i,j)][0] = - in[this->HIDX(i,j)][1] * this->gradx_coe[i];
			out[this->HIDX(i,j)][1] =   in[this->HIDX(i,j)][0] * this->gradx_coe[i];

			out[this->R_HIDX(i,j)][0] = - in[this->R_HIDX(i,j)][1] * this->gradx_coe[this->reflectedXWavenumberIndex(i)];
			out[this->R_HIDX(i,j)][1] =   in[this->R_HIDX(i,j)][0] * this->gradx_coe[this->reflectedXWavenumberIndex(i)];
		}
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::grady(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<this->HALF_YPTS; ++j) {
		for(int i=0; i<this->HALF_XPTS; ++i) {
			out[this->HIDX(i,j)][0]   = - in[this->HIDX(i,j)][1] * this->grady_coe[j];
			out[this->HIDX(i,j)][1]   =   in[this->HIDX(i,j)][0] * this->grady_coe[j];

			out[this->R_HIDX(i,j)][0] = - in[this->R_HIDX(i,j)][1] * this->grady_coe[j];
			out[this->R_HIDX(i,j)][1] =   in[this->R_HIDX(i,j)][0] * this->grady_coe[j];
		}
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::laplacian(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = - in[i][0] * this->laplacian_coe[i];
		out[i][1] = - in[i][1] * this->laplacian_coe[i];
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::invertLaplacian(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = - in[i][0] / this->laplacian_coe_inverse[i];
		out[i][1] = - in[i][1] / this->laplacian_coe_inverse[i];
	}
}

template<int XPTS, int YPTS> void fftwf_operation<XPTS,YPTS>::dealiase(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<this->HALF_GRIDS; ++i) {
		out[i][0] = in[i][0] * this->dealiasing_mask[i];
		out[i][1] = in[i][1] * this->dealiasing_mask[i];
	}
}
#endif
