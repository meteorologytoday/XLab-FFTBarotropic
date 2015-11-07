#define LX (600000.0f)
#define LY (600000.0f)
#define NU (50.0f)

#define NPTS (768)


// In this model XPTS == YPTS by default because of dealiasing issue.
#define XPTS (NPTS)
#define YPTS (NPTS)

#define GRIDS (XPTS*YPTS)

#define HALF_XPTS ((int)((XPTS/2)+1))
#define HALF_YPTS ((int)((YPTS/2)+1))
#define HALF_GRIDS (HALF_XPTS*HALF_YPTS)

#define IDX(i,j) (XPTS*j+i)
#define HIDX(i,j) (HALF_XPTS*j+i)

#define SWAP(a,b) {a^=b,b^=a,a^=b;}
