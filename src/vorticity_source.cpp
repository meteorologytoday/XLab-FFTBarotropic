#include "configuration.hpp"

#ifndef VORTICITY_SOURCE
#define VORTICITY_SOURCE

namespace VORT_SRC_READER {

enum VORTICITY_RECIPE_TYPE { SCRIPT, FIFO, EMPTY };

/* A vorticity source input can be
 * 1. Script file with format specified in [#1]
 * 2. A file descriptor which reads [GRIDS] float point in binary each step
 *
 * #1 A script file's format: 
 *    [time] [binary filename]<newline>
 *
 *
 */


struct Recipe {
	float time,
	char filename[256]
}

class VortSrcRecipe {
	private:
	RECIPE_TYPE recipe_type;
	string filename;
	Recipe * recipe;
	float * vort_src;
	FILE * fd;

	public:
	VortSrcRecipe(RECIPE_TYPE recipe_type, string filename, float * vort_src) {
		if(recipe_type != SCRIPT && recipe_type != FIFO && recipe_type != EMPTY) {
			printf("ERROR: vorticity source recipe type is not valid.\n");
		}
		this->recipe_type = recipe_type;
		this->filename = filename;
		this->vort_src = vort_src;

		if(recipe_type == SCRIPT) {
			this->readScript();
		}

		if((fd = fopen(filename.c_str(), "rb")) == NULL) {
			printf("ERROR: cannot open file [%s].\n", filename.c_str());
		}

		
	}

	void readScript() {
		// open file
	}

	void read(float time) {
		
	}
};

#endif
