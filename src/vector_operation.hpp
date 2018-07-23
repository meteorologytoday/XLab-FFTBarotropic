#ifndef VECTOR_OPERATION_HPP
#define VECTOR_OPERATION_HPP

template<int PTS> class VectorOperation {
private:
public:
	//VectorOperation();
	//~VectorOperation();
	void add(float *out, float *in1, float *in2);
	void sub(float *out, float *in1, float *in2);
	void mul(float *out, float *in1, float *in2);
	void div(float *out, float *in1, float *in2);
	
	void add(float *out, float *in, float k);
	void sub(float *out, float *in, float k);
	void mul(float *out, float *in, float k);
	void div(float *out, float *in, float k);

	inline void iadd(float *out, float *in) { this->add(out, out, in) };
	inline void isub(float *out, float *in) { this->sub(out, out, in) };
	inline void imul(float *out, float *in) { this->mul(out, out, in) };
	inline void idiv(float *out, float *in) { this->div(out, out, in) };
	
	inline void iadd(float *out, float k) { this->add(out, out, k) };
	inline void isub(float *out, float k) { this->sub(out, out, k) };
	inline void imul(float *out, float k) { this->mul(out, out, k) };
	inline void idiv(float *out, float k) { this->div(out, out, k) };
	
};


#endif
