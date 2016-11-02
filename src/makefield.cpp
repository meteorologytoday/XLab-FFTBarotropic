#include <cstdio>

#define _USE_MATH_DEFINES
#include <cmath>
#include "configuration.h"
#include <errno.h>

int main() {
	
	float centerx = Lx / 2.0, centery = Ly / 2.0, epsilon = 0.7, lambda = 2.0, zeta0 = .005f, r_i = 30000.0, r_o = 60000.0;
	float r, r_i_alpha, r_o_alpha, r_prime;

	auto radius = [centerx, centery](float x, float y) -> float {
		return sqrtf(pow(x-centerx,2) + pow(y-centery,2));
	};
	auto alpha = [centerx, centery, epsilon, radius](float x, float y) -> float {
		float c;
		float r = radius(x,y);
		if(r == 0.0f) {
			c = 0;
		} else {
			c = (y-centery) / radius(x,y);
		}
		return sqrtf((1.0 - pow(epsilon,2)) / (1.0 - pow(epsilon*c,2)));
	};

	float x, y;
	for(int i=0; i<XPTS; ++i) {
		x = i * dx;
		for(int j=0; j<YPTS; ++j) {
			y = j * dy;
			r = radius(x,y);
			r_i_alpha = r_i * alpha(x,y);
			r_o_alpha = r_o * alpha(x,y);

			if(r <= r_i_alpha) {
				vort[IDX(i,j)] = zeta0;
			} else if (r <= r_o_alpha) {
				r_prime = (r - r_i_alpha)/(r_o_alpha - r_i_alpha);
				vort[IDX(i,j)] = zeta0 * (1.0 - exp( - lambda / r_prime * exp(1.0 / (r_prime - 1))));
			} else {
				vort[IDX(i,j)] = 0;
			}
		}
	}
	sprintf(filename, "initial.bin");
	writeField(filename, vort);

	return 0;
}

