#define LX (600000f)
#define LY (600000f)

#define XPTS (768)
#define YPTS (768)

#define GRIDS (XPTS*YPTS)

#define HALF_XPTS ((int)((XPTS/2)+1))
#define HALF_YPTS ((int)((YPTS/2)+1))
#define HALF_GRIDS (HALF_XPTS*HALF_YPTS)
#define DEALIASE_XWAVENUMBER ((int)(XPTS/3))
#define DEALIASE_YWAVENUMBER ((int)(YPTS/3))

#define IDX(i,j) (XPTS*j+i)
#define HIDX(i,j) (HALF_XPTS*j+i)

#define SWAP(a,b) {a^=b,b^=a,a^=b;}
