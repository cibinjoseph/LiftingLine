function vel = hshoe(P, A, B, C, D, gam)
    vel = vortxl(P, A, B, gam) + vortxl(P, B, C, gam) + vortxl(P, C, D, gam);
    return;
