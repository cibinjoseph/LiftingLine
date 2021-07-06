% Lifting Line theory (horseshoe vortices method)
clc; clear; clf;

% Parameters
span = 6;
chordMax = 1.0;
ns = 40;
theta = 5 * (pi / 180);
cA = [100; 0; 0];
cD = [100; 0; 0];
Vinf = [10; 0; 0];
rho = 1.2;
CLslope = 2 * pi;

% Create geometry
semiSpan = span * 0.5;
thetaSpan = linspace(0.001, 0.999 * pi, ns);
chord = chordMax * sin(thetaSpan);
y = - semiSpan * cos(thetaSpan);
yvec = [- semiSpan 0.5 * (y(1:end - 1) + y(2:end)) semiSpan];
dy = diff([- semiSpan 0.5 * (y(2:end) + y(1:end - 1)) semiSpan]);
wingArea = sum(chord .* dy);
chordMean = wingArea / span;
AR = span * span / wingArea;

cB = zeros(3, ns); cB(1, :) = 0.25 * chord;
cC = cB;
cB(2, :) = yvec(1:end - 1);
cC(2, :) = yvec(2:end);

cp = zeros(3, ns);
cp(1, :) = 0.75 * chord;
cp(2, :) = y;

indVelcp = zeros(3, ns);
indVelcp(1, :) = 0.25 * chord;
indVelcp(2, :) = y;

ncap = [sin(theta); 0; cos(theta)];

% % Plot wing and cp
% plot(cB(1,:),cB(2,:),'ro')
% hold on;
% plot(cC(1,:),cC(2,:),'b*')
% plot(cp(1,:),cp(2,:),'bo')

% AIC
AIC = zeros(ns, ns);
for is = 1:ns % panel
  cA(2) = cB(2, is);
  cD(2) = cC(2, is);
  for icp = 1:ns % cp
    AIC(is, icp) = dot(hshoe(cp(:, icp), cA, cB(:, is), cC(:, is), cD, 1), ncap);
  end
end
RHS = ones(ns, 1) * dot(- Vinf, ncap);
gam = AIC \ RHS;

% Calculate induced velocity
indVel = zeros(1, ns);
for is = 1:ns
  for ipanel = 1:ns
    cA(2) = cB(2, ipanel);
    cD(2) = cC(2, ipanel);
    indVel(is) = indVel(is) + dot(vortxl(indVelcp(:, is), cA, cB(:, ipanel), gam(ipanel)), [0; 0; 1]);
    indVel(is) = indVel(is) + dot(vortxl(indVelcp(:, is), cC(:, ipanel), cD, gam(ipanel)), [0; 0; 1]);
    if (is ~= ipanel)
      indVel(is) = indVel(is) + dot(vortxl(indVelcp(:, is), cB(:, ipanel), cC(:, ipanel), gam(ipanel)), [0; 0; 1]);
    end
  end
end

% Lift from gam
lift_gam = rho * norm(Vinf) * gam;
lift_nondim_gam = lift_gam / (0.5 * rho * norm(Vinf) ^ 2 * chordMean * CLslope * theta);
cl_gam = lift_gam / (0.5 * rho * (norm(Vinf)) ^ 2 * span * chordMean);

% ==================================================================

% Calculate alpha
alpha = theta - atan(abs(indVel) / norm(Vinf));

% Lift from alpha
cl_alpha = CLslope * alpha;
lift_alpha = cl_alpha .* (0.5 * rho * (norm(Vinf)) ^ 2 * dy * chordMean);

% Theoretical
liftNonDimTheoretical = getLiftNonDimTheoretical(AR, CLslope, semiSpan, y);

plotOpt = 3;
switch plotOpt
  case 1
    plot(y, lift_nondim_gam, 'bo');
    hold on;
    plot(y, liftNonDimTheoretical, 'r-');
    grid on
    legend('numerical - Hshoe', 'Prandtl')
    title(['Aspect ratio   ' num2str(AR)])
    xlabel('spanwise stations')
    ylabel('Non-dim lift')

  case 2
    plot(y, gam, 'bo');
    grid on
    title(['Aspect ratio   ' num2str(AR)])
    xlabel('spanwise stations')
    ylabel('Gam')

  case 3
    plot(y, abs(indVel))
    grid on
    title(['Aspect ratio   ' num2str(AR)])
    xlabel('spanwise stations')
    ylabel('Vi')
end
