function vel = vortxl(P,A,B,gam)
r0 = B-A;
r1 = P-A;
r2 = P-B;

vel = gam*cross(r1,r2)/(4*pi*(norm(cross(r1,r2)))^2);
vel = vel*dot(r0,r1/norm(r1)-r2/norm(r2));
return;
