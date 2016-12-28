#include <cmath>
#include <iostream>
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

const float rho = 1.0f;
const float f = 1e-5;

const float L = 600000.0f;

const float LX = L;
const float LY = L;
const float NU = 50.0f;
const int NPTS = 768;

const int XPTS = NPTS;
const int YPTS = NPTS;

const int GRIDS = XPTS * YPTS;
const int HALF_XPTS = ((int)(XPTS/2))+1;
const int HALF_YPTS = ((int)(YPTS/2))+1;
const int HALF_GRIDS = XPTS * HALF_YPTS;

inline int IDX(int i, int j) { return (YPTS*i + j); }
inline int HIDX(int i, int j) { return (HALF_YPTS*i + j); }

const int total_steps = 10; //100.0 * 60.0 / 3.0;
const float dt = 3.0f;

std::string input = "input";
std::string output = "output";
std::string init_file = "initial_vorticity.bin";

#endif
