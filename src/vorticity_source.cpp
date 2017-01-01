#include "configuration.hpp"

#ifndef VORTICITY_SOURCE
#define VORTICITY_SOURCE

/* A vorticity source input can be
 * 1. Script file with format specified in [#1]
 * 2. A file descriptor which reads [GRIDS] float point in binary each step
 *
 * #1 A script file's format: 
 *    [time] [binary filename]<newline>
 *
 *
 */

struct vort_receipe {
	size_t step,
	struct vort_script * next_recipe
}

void initVortSrcScript() {
	
}

void readVortSrc(float * data, size_t step) {
	if(step == 0) {
		initVortSrcScript();
	}
}

#endif
