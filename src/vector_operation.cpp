#include "vector_operation.hpp"

#ifndef VECTOR_OPERATION_CPP
#define VECTOR_OPERATION_CPP
template<int PTS> void VectorOperation<PTS>::add(float *out, float *in1, float *in2) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in1[i] + in2[i];
	}
}

template<int PTS> void VectorOperation<PTS>::sub(float *out, float *in1, float *in2) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

template<int PTS> void VectorOperation<PTS>::mul(float *out, float *in1, float *in2) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in1[i] * in2[i];
	}
}

template<int PTS> void VectorOperation<PTS>::div(float *out, float *in1, float *in2) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in1[i] / in2[i];
	}
}


// With scalar

template<int PTS> void VectorOperation<PTS>::add(float *out, float *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in[i] + k;
	}
}
template<int PTS> void VectorOperation<PTS>::sub(float *out, float *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in[i] - k;
	}
}
template<int PTS> void VectorOperation<PTS>::mul(float *out, float *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in[i] * k;
	}
}
template<int PTS> void VectorOperation<PTS>::div(float *out, float *in, float k) {
	for(int i=0; i < GRIDS; ++i) {
		out[i] = in[i] / k;
	}
}


#endif
