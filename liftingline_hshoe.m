% Lifting Line theory (horseshoe vortices method)
clc; clear; clf;

ll_fourier
hold on

% Parameters
span = 6;
chord0 = 1;
ns = 100;
theta = 5*(pi/180);
cA = [100 ; 0 ; 0];
cD = [100 ; 0 ; 0];
Vinf = [10 ; 0 ; 0];
rho = 1.2;

% Create geometry
thetaSpan = linspace(0.001, 0.999*pi, ns+1);
chord = chord0*sin(thetaSpan); chord = 0.5*(chord(1:end-1)+chord(2:end));
yvec = -0.5*span*cos(thetaSpan); 
yElem = 0.5*(yvec(1:end-1)+yvec(2:end));
dy = yvec(2:end)-yvec(1:end-1);
chord_mean = sum(abs(yElem).*chord)/span;

cB = zeros(3,ns); cB(1,:) = 0.25*chord;
cC = cB;
cB(2,:) = yvec(1:end-1);
cC(2,:) = yvec(2:end);

cp = zeros(3,ns);
cp(1,:) = 0.75*chord;
cp(2,:) = yElem;

indVelcp = zeros(3,ns);
indVelcp(1,:) = 0.25*chord;
indVelcp(2,:) = yElem;

ncap = [sin(theta); 0; cos(theta)];

% % Plot wing and cp
% plot(cB(1,:),cB(2,:),'ro')
% hold on;
% plot(cC(1,:),cC(2,:),'b*')
% plot(cp(1,:),cp(2,:),'bo')

% AIC
AIC = zeros(ns,ns);
for is = 1:ns   % panel
  cA(2) = cB(2,is);
  cD(2) = cC(2,is);
  for icp = 1:ns   % cp
    AIC(is,icp) = dot(hshoe(cp(:,icp),cA,cB(:,is),cC(:,is),cD,1),ncap);
  end
end
RHS = ones(ns,1)*dot(-Vinf,ncap);
gam = AIC\RHS;

% Lift from gam
lift_gam = rho*norm(Vinf)*gam.*dy;
cl_gam = lift_gam/(0.5*rho*(norm(Vinf))^2*span*chord_mean);

% ==================================================================

% Calculate induced velocity
indVel = zeros(1,ns);
for is = 1:ns
  for ipanel = 1:ns
    cA(2) = cB(2,ipanel);
    cD(2) = cC(2,ipanel);
    indVel(is) = indVel(is) + dot(vortxl(indVelcp(:,is),cA,cB(:,ipanel),gam(ipanel)),[0;0;1]);
    indVel(is) = indVel(is) + dot(vortxl(indVelcp(:,is),cC(:,ipanel),cD,gam(ipanel)),[0;0;1]);
  end
end

% Calculate alpha
alpha = theta-atan(abs(indVel)/norm(Vinf));

% Lift from alpha
cl_alpha = 2*pi*alpha;
lift_alpha = cl_alpha.*(0.5*rho*(norm(Vinf))^2*dy*chord_mean);

% plot(y,lift_gam,'bo-')
% hold on;
% plot(y,lift_alpha,'ro-')
% legend('gam','alpha')
% grid on

plot(yElem,abs(indVel),'ro-')
grid on
