float dx, dy, Lx, Ly;

float *vort, *divg, *h, *bg_vort, \
      *u2, *u, *u_divg, *u_vort, \
      *v2, *v, *v_divg, *v_vort, \
      *K, *E, *geop, \
      *absvort, *absvort_u, *absvort_v, \
      *geop_u, *geop_v, \
       *dvortdx, *dvortdy, *dvortdt, *workspace, *vort_src;

fftwf_complex *vort_c0, *vort_c, *lvort_c, \
              *divg_c0, *divg_c, *ldivg_c, \
              *geop_c0,    *geop_c,    *lgeop_c,    \
              *absvort_u_c, *absvort_v_c,
              *geop_u_c, *geop_v_c,
              *E_c,
              *chi_c, *psi_c, \
              *tmp_c, **rk4_vort_c[4], **rk4_divg_c[4], **rk4_geop_c[4],
              *copy_for_c2r,
              *dvortdt_c, *ddivgdt_c, *dgeopdt_c;

fftwf_plan p_fwd_vort,    p_bwd_vort,
           p_fwd_divg,    p_bwd_divg,
           p_fwd_geop,    p_bwd_geop,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_bwd_u_vort,  p_bwd_v_vort,
		   p_bwd_u_divg,  p_bwd_v_divg,
		   p_fwd_dvortdt, 
           p_bwd_chi,     p_bwd_psi,
           p_fwd_absvort_u,
           p_fwd_absvort_v,
           p_fwd_geop_u, p_fwd_geop_v,
           p_fwd_E;



fftwf_operation<XPTS,YPTS> fop(LX, LY);

VectorOperation<GRIDS> vop;

char filename[256];


