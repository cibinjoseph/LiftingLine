function vel = vortxl(P,A,B,gam)
% r0 = B-A;
r1 = P-A;
r2 = P-B;

% vel = gam*cross(r1,r2)/(4*pi*(norm(cross(r1,r2)))^2);
% vel = vel*dot(r0,r1/norm(r1)-r2/norm(r2));

r1_abs = norm(r1);
r2_abs = norm(r2);

vel = (r1_abs+r2_abs)*cross(r1,r2)/(r1_abs*r2_abs*(r1_abs*r2_abs+dot(r1,r2)));
vel = gam/(4*pi)*vel;
return;
