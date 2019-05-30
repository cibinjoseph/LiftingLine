% Prandtl lifting line method (fourier coefficients method)
clc; clear; clf;

ns = 24;
nFourier = ns;

Vinf = 10;
cl_slope = 2*pi;
theta = 5*(pi/180);
span = 8;
chord_max = 1;
semispan = span/2;
rho = 1.2;

thetaSpan = linspace(0.001*pi, 0.999*pi, ns);
chord = chord_max*sin(thetaSpan);
y = -semispan*cos(thetaSpan);
dy = diff([-semispan 0.5*(y(2:end)+y(1:end-1)) semispan]);
wing_area = sum(chord.*dy);
chord_mean = wing_area/span;
AR = span*span/wing_area;

% Solving for fourier coefficients in each span station
a = zeros(ns,nFourier);
for ispan = 1:ns
  for n = 1:nFourier
    a(ispan,n) = sin(n*thetaSpan(ispan))*(sin(thetaSpan(ispan)) + n*cl_slope*chord(ispan)/(8*semispan));
  end
end

RHS = zeros(ns,1);
for ispan = 1:ns
  RHS(ispan) = ((cl_slope*chord(ispan))/(8*semispan))*sin(thetaSpan(ispan))*theta;
end

x = a\RHS;

% Circulation
gam = zeros(ns,1);
for ispan = 1:ns
  for n = 1:length(x)
    sumterm = 4*semispan*Vinf*x(n)*sin(n*thetaSpan(ispan));
    gam(ispan) = gam(ispan) + sumterm;
  end
end

% Lift per unit span
lift = rho*Vinf*gam;
lift_nondim = lift/(0.5*rho*Vinf*Vinf*chord_mean*cl_slope*theta);

% Induced velocity
wi = zeros(ns,1);
for ispan = 1:ns
  for n = 1:length(x)
    sumterm = Vinf*(n*x(n)*sin(n*thetaSpan(ispan)))/sin(thetaSpan(ispan));
    wi(ispan) = wi(ispan) + sumterm;
  end
end

alfaind = wi/Vinf;
cl_local = cl_slope*(theta - alfaind);

% Theoretical
lift_nondim_theoretical = 4*AR/(pi*AR+cl_slope)*sqrt(1-(y/semispan).^2);

plot(y,lift_nondim,'bo')
hold on;
plot(y,lift_nondim_theoretical,'r-')
grid on;
legend('numerical - Fourier','Prandtl')
title(['Aspect ratio   ' num2str(AR)])
xlabel('spanwise stations')
ylabel('Non-dim lift')
