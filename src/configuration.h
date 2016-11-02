#include <cmath>
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

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

#endif
