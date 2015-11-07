#ifndef FFTWF_OPERATION_H
#define FFTWF_OPERATION_H
#include "configuration.h"
#include <fftw3.h>

#include <cmath>

const float TWOPI = (acos(0.0f) * 2.0f);

class fftwf_operation {
private:
	float *gradx_coe, *grady_coe, *laplacian_coe, *laplacian_coe_inverse, *dealiasing_mask;
	int dealiase_xwavenumber, dealiase_ywavenumber;
public:
	fftwf_operation();
	~fftwf_operation();

	void gradx(fftwf_complex *in, fftwf_complex *out);
	void grady(fftwf_complex *in, fftwf_complex *out);
	void laplacian(fftwf_complex *in, fftwf_complex *out);
	void invertLaplacian(fftwf_complex *in, fftwf_complex *out);
	void dealiase(fftwf_complex *in, fftwf_complex *out);
};

#endif
