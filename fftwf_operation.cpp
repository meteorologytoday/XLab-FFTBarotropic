#include "fftwf_operation.h"

#ifndef FFTWF_OPERATION_CPP
#define FFTWF_OPERATION

fftwf_operation::fftwf_operation() {
	this->gradx_coe = (float*) fftwf_malloc(sizeof(float) * HALF_XPTS);
	this->grady_coe = (float*) fftwf_malloc(sizeof(float) * HALF_YPTS);
	this->laplacian_coe_inverse = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->laplacian_coe = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->dealiasing_mask = (float*) fftwf_malloc(sizeof(float) * HALF_GRIDS);
	this->dealiase_xwavenumber = (int) ceil(((float)XPTS)/3.0);
	this->dealiase_ywavenumber = (int) ceil(((float)YPTS)/3.0);

	for(int i=0; i<HALF_XPTS; ++i) {
		this->gradx_coe[i] = TWOPI * ((float) i) / LX;
	}

	for(int j=0; j<HALF_YPTS; ++j) {
		this->grady_coe[j] = TWOPI * ((float) j) / LY;
	}

	for(int i=0; i<HALF_XPTS; ++i) {
		for(int j=0; j<HALF_YPTS; ++j) {
			this->laplacian_coe_inverse[IDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
			this->laplacian_coe[IDX(i,j)] = - (pow(this->gradx_coe[i],2) + pow(this->grady_coe[j],2));
		}
	}
	this->laplacian_coe_inverse[IDX(0,0)] = 1.0; // this coe is special

	float generalized_wavenumber_square = pow(this->dealiase_xwavenumber,2) + pow(this->dealiase_ywavenumber,2);
	for(int i = this->dealiase_xwavenumber; i < HALF_XPTS ; ++i) {
		for(int j = this->dealiase_ywavenumber; j < HALF_YPTS; ++j) {
			this->dealiasing_mask[HIDX(i,j)] = (pow(i,2) + pow(j,2) >= generalized_wavenumber_square) ? 0.0f : 1.0f;
		}
	}
}

fftwf_operation::~fftwf_operation() {
	fftwf_free(this->gradx_coe);
	fftwf_free(this->grady_coe);
	fftwf_free(this->laplacian_coe);
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

void fftwf_operation::laplacian(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<HALF_YPTS; ++j) {
		for(int i=0; i<HALF_XPTS; ++i) {
			out[HIDX(i,j)][0] = - in[HIDX(i,j)][0] * laplacian_coe[HIDX(i,j)];
			out[HIDX(i,j)][1] = - in[HIDX(i,j)][1] * laplacian_coe[HIDX(i,j)];
		}
	}
}

void fftwf_operation::invertLaplacian(fftwf_complex *in, fftwf_complex *out){
	for(int j=0; j<HALF_YPTS; ++j) {
		for(int i=0; i<HALF_XPTS; ++i) {
			out[HIDX(i,j)][0] = - in[HIDX(i,j)][0] / laplacian_coe_inverse[HIDX(i,j)];
			out[HIDX(i,j)][1] = - in[HIDX(i,j)][1] / laplacian_coe_inverse[HIDX(i,j)];
		}
	}
}

void fftwf_operation::dealiase(fftwf_complex *in, fftwf_complex *out){
	for(int i=0; i<HALF_GRIDS; ++i) {
		out[i][0] = in[i][0] * dealiasing_mask[i];
		out[i][1] = in[i][1] * dealiasing_mask[i];
	}
}

#endif
