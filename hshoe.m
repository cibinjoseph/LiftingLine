function vel = hshoe(P,A,B,C,D,gam)
vel=vortxl(P,A,B,gam);
vel=vel+vortxl(P,B,C,gam);
vel=vel+vortxl(P,C,D,gam);
return;
