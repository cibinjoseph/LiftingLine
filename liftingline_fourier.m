% Prandtl lifting line method (fourier coefficients method)
clc; clear; clf;

ns = 24;
nFourier = ns;

Vinf = 10;
CLslope = 2*pi;
theta = 5*(pi/180);
span = 8;
chordMax = 1;
semiSpan = span/2;
rho = 1.2;

thetaSpan = linspace(0.001*pi, 0.999*pi, ns);
chord = chordMax*sin(thetaSpan);
y = -semiSpan*cos(thetaSpan);
dy = diff([-semiSpan 0.5*(y(2:end)+y(1:end-1)) semiSpan]);
wingArea = sum(chord.*dy);
chordMean = wingArea/span;
AR = span*span/wingArea;

% Solving for fourier coefficients in each span station
a = zeros(ns,nFourier);
for ispan = 1:ns
  for n = 1:nFourier
    a(ispan,n) = sin(n*thetaSpan(ispan))*(sin(thetaSpan(ispan)) + n*CLslope*chord(ispan)/(8*semiSpan));
  end
end

RHS = zeros(ns,1);
for ispan = 1:ns
  RHS(ispan) = ((CLslope*chord(ispan))/(8*semiSpan))*sin(thetaSpan(ispan))*theta;
end

x = a\RHS;

% Circulation
gam = zeros(ns,1);
for ispan = 1:ns
  for n = 1:nFourier
    sumterm = 4*semiSpan*Vinf*x(n)*sin(n*thetaSpan(ispan));
    gam(ispan) = gam(ispan) + sumterm;
  end
end

% Lift per unit span
lift = rho*Vinf*gam;
lift_nondim = lift/(0.5*rho*Vinf*Vinf*chordMean*CLslope*theta);

% Induced velocity
Vi = zeros(ns,1);
for ispan = 1:ns
  for n = 1:nFourier
    sumterm = Vinf*(n*x(n)*sin(n*thetaSpan(ispan)))/sin(thetaSpan(ispan));
    Vi(ispan) = Vi(ispan) + sumterm;
  end
end

alfaind = Vi/Vinf;
cl_local = CLslope*(theta - alfaind);

% Theoretical
liftNonDimTheoretical = getLiftNonDimTheoretical(AR,CLslope,semiSpan,y);

plot(y,lift_nondim,'bo')
hold on;
plot(y,liftNonDimTheoretical,'r-')
grid on;
legend('numerical - Fourier','Prandtl')
title(['Aspect ratio   ' num2str(AR)])
xlabel('spanwise stations')
ylabel('Non-dim lift')
