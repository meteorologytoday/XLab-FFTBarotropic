#ifndef FFTWF_OPERATION_H
#define FFTWF_OPERATION_H
#include "configuration.h"
#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <cmath>

#define M_2PI (2.0*M_PI)

class fftwf_operation {
private:
	float *gradx_coe, *grady_coe, *laplacian_coe, *dealiasing_mask;
public:
	fftwf_operation();
	~fftwf_operation();

	void gradx(fftwf_complex *in, fftwf_complex *out);
	void grady(fftwf_complex *in, fftwf_complex *out);
	void invertLaplacian(fftwf_complex *in, fftwf_complex *out);
	void dealiase(fftwf_complex *inout);
};

#endif
