#define XPTS 768
#define YPTS 768

#define GRIDS (XPTS*YPTS)

#define HALF_XPTS ((int)((XPTS/2)+1))
#define HALF_YPTS ((int)((YPTS/2)+1))

#define HALF_GRIDS (HALF_XPTS*HALF_YPTS)

#define IDX(i,j) (XPTS*j+i)

#define SWAP(a,b) {a^=b,b^=a,a^=b; }
