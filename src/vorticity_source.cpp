#include <cstdio>
#include <cstdlib>
#include <string>
#include "configuration.hpp"

#ifndef VORTICITY_SOURCE
#define VORTICITY_SOURCE

namespace VORT_SRC_READER {

enum RECIPE_TYPE { SCRIPT, FIFO, EMPTY };

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
	float time;
	char filename[256];
};

void stripComment(char * str, size_t len) {
	for(size_t i = 0; i < len; ++i) {
		if(str[i] == '#') { // comment found
			str[i] = '\0';
			break;
		}
	}
}

int readLine(FILE * file, char * str, size_t len) {
	if(fgets(str, len, file) != NULL) {
		stripComment(str, len);
		return 0;
	} else {
		return -1;
	}
}


template<size_t GRIDS> class VortSrcRecipeReader {
	private:
	enum RECIPE_TYPE recipe_type;
	std::string filename;
	struct Recipe * recipe_head;
	float * vort_src;
	FILE * fifo;

	public:
	VortSrcRecipeReader() {
	}

	~VortSrcRecipeReader() {
		fclose(fifo);
	}

	int read(float time) {
		int flag = 0;
		switch(this->recipe_type) {
			case SCRIPT:
				flag = this->readScript();
				break;
			case FIFO:
				flag = this->readFIFO();
				break;
			case EMPTY:
				flag = 0;
				break;
		}
		return flag;
	}

	void init(RECIPE_TYPE recipe_type, std::string filename, float * vort_src) {
		if(recipe_type != SCRIPT && recipe_type != FIFO && recipe_type != EMPTY) {
			printf("ERROR: vorticity source recipe type is not valid.\n");
		}
		this->recipe_type = recipe_type;
		this->filename = filename;
		this->vort_src = vort_src;

		printf("INITIALIZATION\n"); fflush(stdout);

		if(recipe_type == SCRIPT) {
			this->readScript();
		} else if(recipe_type == FIFO) {
			if((this->fifo = fopen(this->filename.c_str(), "rb")) == NULL) {
				printf("ERROR: cannot open file [%s].\n", this->filename.c_str());
			}
		}
		
	}

	private:

	int readScript() {
		// open file
		// constructing
		FILE * fd;
		float time;
		char buf[1024];
		if((fd = fopen(this->filename.c_str(), "r")) == NULL) {
			printf("ERROR: cannot open file [%s].\n", this->filename.c_str());
		}
		return 0;
	}

	int readFIFO() {
		char new_flag;
		
		if(fread(&new_flag, sizeof(char), 1, fifo) != 1) {
			fprintf(stderr, "Don't change vorticity source input.\n");
			return 1; // error while reading new_flag
		}

		if(((unsigned int) new_flag) == 1) { // new input
			fprintf(stderr, "Reading new vorticity source input...");
			if(fread(this->vort_src, sizeof(float), GRIDS, fifo) != GRIDS) {
				// error while reading vort_src
				return 2;
			}
			fprintf(stderr, " done.\n");
		}

		return 0;
	}

};

}
#endif
