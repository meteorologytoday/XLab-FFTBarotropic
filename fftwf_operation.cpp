#include "fftwf_operation.h"

#ifndef FFTWF_OPERATION_CPP
#define FFTWF_OPERATION

fftwf_operation::fftwf_operation() {
	this->gradx_coe = (float*) fftwf_malloc(sizeof(float) * HALF_XPTS);
	this->grady_coe = (float*) fftwf_malloc(sizeof(float) * HALF_YPTS);
	this->laplacian_coe = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->dealiasing_mask = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);

	for(int i=0; i<HALF_XPTS; ++i) {
		this->gradx_coe[i] = M_2PI * ((float) i) / Lx;
	}

	for(int j=0; j<HALF_YPTS; ++j) {
		this->grady_coe[j] = M_2PI * ((float) j) / Ly;
	}

	for(int i=0; i<HALF_XPTS; ++i) {
		for(int j=0; j<HALF_YPTS; ++j) {
			this->laplacian_coe[IDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
		}
	}
	this->laplacian_coe[IDX(0,0)] = 1.0; // this coe is special

	//for(int )
}

fftwf_operation::~fftwf_operation() {
	fftw_free(this->gradx_coe);
	fftw_free(this->grady_coe);
	fftw_free(this->laplacian_coe);
}

void fftwf_operation::gradx(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<HALF_YPTS; ++j) {
		for(int i=0; i<HALF_XPTS; ++i) {
			out[HIDX(i,j)][0] = - in[HIDX(i,j)][1] * gradx_coe[i];
			out[HIDX(i,j)][1] =   in[HIDX(i,j)][0] * gradx_coe[i];
		}
	}
}

void fftwf_operation::grady(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<HALF_YPTS; ++j) {
		for(int i=0; i<HALF_XPTS; ++i) {
			out[HIDX(i,j)][0] = - in[HIDX(i,j)][1] * grady_coe[j];
			out[HIDX(i,j)][1] =   in[HIDX(i,j)][0] * grady_coe[j];
		}
	}
}

void fftwf_operation::invertLaplacian(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<HALF_YPTS; ++j) {
		for(int i=0; i<HALF_XPTS; ++i) {
			out[HIDX(i,j)][0] = - in[HIDX(i,j)][0] / laplacian_coe[HIDX(i,j)];
			out[HIDX(i,j)][1] = - in[HIDX(i,j)][1] / laplacian_coe[HIDX(i,j)];
		}
	}
}

void fftwf_operation::dealiase(fftwf_complex *inout){
	for(int i=0; i<HALF_GRIDS; ++i) {
		inout[i] *= dealiasing_mask[i];
	}
}

#endif
